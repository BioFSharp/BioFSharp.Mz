(**
---
title: SEQUEST-like scoring
category: Peptide identification
categoryindex: 2
index: 3
---
*)

(*** hide ***)

(*** condition: prepare ***)
#r "nuget: FSharpAux, 2.1.0"
#r "nuget: FSharpAux.IO, 2.1.0"
#r "nuget: FSharp.Stats, 0.6.0"
#r "nuget: BioFSharp, 2.0.0"
#r "../src/BioFSharp.Mz/bin/Release/netstandard2.0/BioFSharp.Mz.dll"

(*** condition: ipynb ***)
#if IPYNB
#r "nuget: BioFSharp.Mz, {{fsdocs-package-version}}"
#endif // IPYNB

(**
[![Binder]({{root}}img/badge-binder.svg)](https://mybinder.org/v2/gh/BioFSharp/BioFSharp.Mz/gh-pages?filepath={{fsdocs-source-basename}}.ipynb)&emsp;
[![Script]({{root}}img/badge-script.svg)]({{root}}{{fsdocs-source-basename}}.fsx)&emsp;
[![Notebook]({{root}}img/badge-notebook.svg)]({{root}}{{fsdocs-source-basename}}.ipynb)

# SEQUEST-like scoring

Identifying the peptide behind an MS2 spectrum runs through five steps: the
measured spectrum is preprocessed ([signal detection]({{root}}01_02_signal_detection.html)),
candidate sequences are selected by precursor mass
([peptide search databases]({{root}}02_02_search_databases.html)), a theoretical
spectrum is predicted for every candidate, the predictions are matched against the
measurement, and the closeness of fit is scored. This page walks steps three to five
for one spectrum, using the `SequestLike` module.

The scoring question sounds simple: which candidate's predicted peaks line up
best with the measured ones? Representing both sides as intensity vectors
over 1 Da bins makes a dot product the obvious measure of agreement. On its
own, that measure has a flaw. A measured spectrum with many peaks all over
the m/z axis gets a decent dot product with almost any candidate, simply
because most predictions land near some peak. SEQUEST's cross-correlation
score (xcorr) fixes this by asking how much better the prediction correlates
with the measured spectrum at zero offset than on average across a range of
shifted offsets. That mean off-position correlation is the background, and
subtracting it leaves only the agreement that is specific to the candidate's
peak pattern.

The `SequestLike` implementation follows this idea in a form that lets one
spectrum be scored against many candidates cheaply. The measured spectrum is
binned to 1 Da, and its intensities are square-root scaled and normalized to
the maximum within each of 10 windows. The average of the spectrum's shifted
copies over a delay range of 75 is then subtracted once. Algebraically that
folds the whole background subtraction into the measured side, so scoring a
candidate afterwards reduces to a dot product between the candidate's
predicted intensity vector and this preprocessed measured vector. The module is deliberately named SEQUEST-like.
Its background window differs from the symmetric average in the published
SEQUEST algorithm, so scores follow the same idea as the original without
matching its numbers exactly, and they are comparable within this
implementation.

## Loading the measured spectrum

The running example of these pages is the MS2 scan of the doubly charged
peptide ANLGMEVMHER in `ms2Example.mgf`. Its precursor m/z of 643.803548 and
charge 2 were recovered on the
[charge state determination]({{root}}01_03_charge_state_determination.html)
page. `PeakArray.zip` turns the two raw arrays of the reader into the
`PeakArray` the scorer consumes, as introduced on the
[peaks and peak arrays]({{root}}01_01_peaks_and_peak_arrays.html) page.
*)

open BioFSharp
open BioFSharp.FileFormats.MGF
open BioFSharp.IO
open BioFSharp.Mz

let ms2 =
    MGF.read (__SOURCE_DIRECTORY__ + "/data/ms2Example.mgf")
    |> List.head

let precursorMZ =
    match MGFEntry.tryGetPrecursorMZ ms2 with
    | Some mz -> mz
    | None -> failwith "no precursor m/z in the MS2 header"

let spectrum : PeakArray<Peak> = PeakArray.zip ms2.Mass ms2.Intensity

printfn "peaks: %i covering m/z %.2f to %.2f"
    spectrum.Length spectrum.[0].Mz spectrum.[spectrum.Length - 1].Mz

(**
```text
peaks: 971 covering m/z 100.67 to 1337.63
```

## Assembling the candidate list

In a production pipeline the candidates come out of a
[search database]({{root}}02_02_search_databases.html) query for all peptides
within a narrow mass window around the measured precursor mass. To keep this
page self-contained we build such a result list by hand with
`SearchDB.createLookUpResult`, which produces the same `LookUpResult` records
the database returns. The neutral mass of each candidate is computed in code
as its residue masses summed up plus one water for the termini.

The list holds the true peptide and three pretenders. Two are permutations
of ANLGMEVMHER that keep the C-terminal arginine, so they weigh exactly the
same. The third replaces the asparagine with two glycines, a substitution
that is also exactly isobaric because the asparagine residue and two glycine
residues share the elemental composition. A mass window query cannot
distinguish any of these by mass. Only the fragment pattern can.
*)

let mono : IBioItem -> float = BioItem.monoisoMass

let neutralMass (peptide: AminoAcids.AminoAcid list) =
    (peptide |> List.sumBy mono) + mono ModificationInfo.Table.H2O

// Stand-ins for what a search database mass window query would return;
// built by hand to keep the page self-contained.
let candidate modSeqId pepSeqId (sequence: string) =
    let bioSequence = BioList.ofAminoAcidString sequence
    let mass = neutralMass bioSequence
    let roundedMass = int64 (System.Math.Round(mass * 1000000.))
    SearchDB.createLookUpResult modSeqId pepSeqId mass roundedMass sequence bioSequence 0

let candidates =
    [ candidate 1 1 "ANLGMEVMHER"   // the peptide the scan was recorded from
      candidate 2 2 "MANGLEVMHER"   // permutation, same mass
      candidate 3 3 "EVMANLGMHER"   // permutation, same mass
      candidate 4 4 "AGGLGMEVMHER"  // Asn replaced by Gly-Gly, also isobaric
    ]

candidates
|> List.iter (fun c -> printfn "%-13s %10.4f Da" c.StringSequence c.Mass)

(**
```text
ANLGMEVMHER    1285.5907 Da
MANGLEVMHER    1285.5907 Da
EVMANLGMHER    1285.5907 Da
AGGLGMEVMHER   1285.5907 Da
```

All four candidates sit at 1285.59 Da, within about two millidaltons of
the neutral precursor mass determined from the MS1 scan.

## Predicting a theoretical spectrum per candidate

Each candidate is paired with its predicted fragment masses.
`Fragmentation.Series.fragmentMasses` produces them exactly as on the
[in silico fragmentation]({{root}}02_01_fragmentation.html) page: the b and y
series of the candidate sequence as `TargetMasses`, and the same series of
the reversed sequence as `DecoyMasses`. We use monoisotopic masses
throughout.

`SequestLike.getTheoSpecs` then converts every pair into a
`TheoreticalSpectrum`, holding one binned intensity vector for the target and
one for the decoy. The intensity of each predicted peak comes from a simple
model (`predictIntensitySimpleModel`): main series ions get the full
predicted intensity, loss peaks and minor series a fraction of it,
everything divided by the charge the ion is predicted at. Every fragment is
laid down once per charge from 1 up to the precursor charge, and the vector
is binned to 1 Da like the measured side, so the two can be compared bin by
bin.

The scan limits define the m/z range of that binning. We use 100 to 1300,
which covers the recorded peaks starting at m/z 100.67 as well as the
heaviest interesting fragment, the singly protonated full-length y ion at
m/z 1286.6. The upper border only cuts away a single stray peak at m/z
1337.6, and no peak of the scan rounds exactly onto a border, where the
binning would drop it.
*)

let scanlimits = 100., 1300.

let fragmentPairs =
    candidates
    |> List.map (fun c ->
        let fragments =
            Fragmentation.Series.fragmentMasses
                Fragmentation.Series.bOfBioList
                Fragmentation.Series.yOfBioList
                mono
                c.BioSequence
        c, fragments)

let theoSpecs = SequestLike.getTheoSpecs scanlimits 2 fragmentPairs

open FSharp.Stats

let countOccupied v =
    v |> Vector.toArray |> Array.filter (fun x -> x > 0.) |> Array.length

printfn "bins per vector: %i" (Vector.length theoSpecs.Head.TheoSpec)

theoSpecs
|> List.iter (fun ts ->
    printfn "%-13s target bins occupied: %i  decoy bins occupied: %i"
        ts.LookUpResult.StringSequence
        (countOccupied ts.TheoSpec)
        (countOccupied ts.DecoyTheoSpec))

(**
```text
bins per vector: 1200
AGGLGMEVMHER  target bins occupied: 90  decoy bins occupied: 90
EVMANLGMHER   target bins occupied: 104  decoy bins occupied: 105
MANGLEVMHER   target bins occupied: 97  decoy bins occupied: 96
ANLGMEVMHER   target bins occupied: 97  decoy bins occupied: 98
```

Every candidate now owns two binned prediction vectors, each occupying
around a hundred of the 1200 bins. `getTheoSpecs` builds its result list by
prepending, so the order is reversed relative to the input, which does not
matter because the scorer ranks by score anyway.

## Preprocessing the measured spectrum

The scorer performs the background subtraction on the measured side once, through
`spectrumToIntensityArrayMinusAutoCorrelation`. `calcSequestScore` calls it
internally, so this block is purely illustrative.
*)

let preprocessed =
    SequestLike.spectrumToIntensityArrayMinusAutoCorrelation scanlimits spectrum

let values = preprocessed |> Vector.toArray

printfn "vector length: %i" values.Length
printfn "positive entries: %i" (values |> Array.filter (fun x -> x > 0.) |> Array.length)
printfn "negative entries: %i" (values |> Array.filter (fun x -> x < 0.) |> Array.length)

(**
```text
vector length: 1200
positive entries: 407
negative entries: 793
```

The negative entries are the signature of the background subtraction. A bin
whose intensity is lower than the local self-correlation average now
penalizes a candidate that predicts a peak there, while a bin that stands
out above the background rewards it. The dot product with a prediction
vector therefore measures alignment beyond what shifted versions of the
spectrum would produce by chance.

## Scoring the candidates

`SequestLike.calcSequestScore` takes the scan limits, the measured spectrum,
the scan time, the precursor charge, the isolation window target m/z (the
precursor m/z from the header), the theoretical spectra and a free-text
spectrum identifier. It scores every target and every decoy vector against
the preprocessed measured vector and returns one `SearchEngineResult` per
scored spectrum, ranked by descending score.
*)

let results =
    SequestLike.calcSequestScore
        scanlimits spectrum 20.93 2 precursorMZ theoSpecs "ms2Example"

let printRanked (rs: SearchEngineResult.SearchEngineResult<float> list) =
    printfn "%-13s %-6s %8s %12s %8s" "sequence" "target" "score" "dBestToRest" "dNext"
    rs
    |> List.iter (fun r ->
        printfn "%-13s %-6b %8.4f %12.4f %8.4f"
            r.StringSequence r.IsTarget r.Score r.NormDeltaBestToRest r.NormDeltaNext)

printRanked results

(**
```text
sequence      target    score  dBestToRest    dNext
ANLGMEVMHER   true    13.0383       0.0000   0.0132
AGGLGMEVMHER  true    12.8665       0.0132   0.3777
MANGLEVMHER   true     7.9418       0.3909   0.2979
EVMANLGMHER   true     4.0572       0.6888   0.1270
ANLGMEVMHER   false    2.4008       0.8159   0.0282
AGGLGMEVMHER  false    2.0329       0.8441   0.1274
MANGLEVMHER   false    0.3723       0.9714   0.0244
EVMANLGMHER   false    0.0548       0.9958   0.0000
```

The target of ANLGMEVMHER comes out on top at 13.04. The two permutations
score far below it, at 7.94 and 4.06. Both end in the same MHER stretch and therefore
predict the true peptide's low y ions, which is why they still beat every
decoy, while their remaining predictions land in the wrong bins.
MANGLEVMHER even shares the y ladder up to y6 with the true sequence, and
that larger overlap is its lead over the other permutation.

The Gly-Gly candidate almost ties the true peptide at 12.87: two glycines weigh
exactly as much as one asparagine, so from b3 onward its b ladder reproduces the
true one bin for bin, and the shared C-terminal nine residues make y1 to y9
identical as well. Only a handful of bins differ, among them the extra b2 of the
glycine pair, and those few bins are the entire margin. Under any mass-based
fragment comparison a sequence whose fragments are isobaric with the true ones is
close to indistinguishable, and the small `dNext` of the top hit records that
ambiguity.

## Reading a SearchEngineResult

Printing the best hit in full shows everything the record carries.
*)

printfn "%A" results.Head

(**
```text
{ SearchEngine = SEQUESTLike
  SpectrumID = "ms2Example"
  ModSequenceID = 1
  PepSequenceID = 1
  GlobalMod = 0
  IsTarget = true
  ScanTime = 20.93
  StringSequence = "ANLGMEVMHER"
  PrecursorCharge = 2
  PrecursorMZ = 643.8035478
  MeasuredMass = 1285.592543
  TheoMass = 1285.590726
  PeptideLength = 11
  Score = 13.03831064
  NormDeltaBestToRest = 0.0
  NormDeltaNext = 0.01317481341 }
```

`SearchEngine` names the scorer that produced the record, so results of
different engines can share one result type. `SpectrumID` is the identifier
string passed into the call, and `ScanTime` travels along the same way, so a
PSM can be traced back to its scan. `ModSequenceID` and `PepSequenceID` are
copied from the `LookUpResult` and tie the hit back to the
[search database]({{root}}02_02_search_databases.html) rows it came from,
and the `GlobalMod` labeling flag rides along with them.
`IsTarget` distinguishes the candidate's own fragment prediction from that
of its reversed decoy sequence. Both records of such a pair share all
identifiers, since they stem from the same database entry.

`StringSequence`, `PrecursorCharge`, `PrecursorMZ` and `PeptideLength`
describe the match itself. `MeasuredMass` is the neutral mass computed from
the isolation window target m/z and the charge, while `TheoMass` is the
candidate's database mass, so the difference between the two is the
precursor mass error of the match. `Score` is the xcorr-style dot product.

The two delta fields put each score into the context of the whole ranking.
`calcNormDeltaBestToRest` fills `NormDeltaBestToRest` with (best score minus this
score) divided by the best score, so the best hit gets 0 and weaker hits approach 1:
MANGLEVMHER's 0.39 means it lost 39 percent of the top score. `calcNormDeltaNext`
fills `NormDeltaNext` with the gap to the next-ranked PSM, normalized by the best
score, with the last PSM getting 0: the top hit's 0.0132 is the Gly-Gly near-tie
discussed above. Both functions expect their input ranked by descending score, which
`calcSequestScore` arranges internally, and a best score of zero or below makes them
return sentinel values instead, 1 and 0 respectively for every PSM. How far a hit
stands out from its competitors is an input to
[false discovery rate control]({{root}}04_01_fdr_control.html).

## What the decoys are for

The ranked table contains eight PSMs for four candidates because every candidate
was also scored as its reversed decoy, the pairing introduced on the
[in silico fragmentation]({{root}}02_01_fragmentation.html) page. All four decoy
scores sit between 0.05 and 2.40, the level a wrong sequence of the right mass
reaches on this spectrum. A single xcorr value has no absolute meaning, and only the
comparison against a population of known wrong matches tells whether a 13 is
convincing. Collected over a whole run, the decoy scores estimate the score
distribution of chance matches, which is what FDR control is built on.

## Scoring many spectra

For a whole run with thousands of spectra, `calcSequestScoreParallel`
distributes the per-candidate scoring with `Async.Parallel`. It takes the
same arguments and produces the same ranking.
*)

let resultsParallel =
    SequestLike.calcSequestScoreParallel
        scanlimits spectrum 20.93 2 precursorMZ theoSpecs "ms2Example"

let bestSequential = results.Head
let bestParallel = resultsParallel.Head

printfn "sequential best: %-13s target=%b score=%.4f"
    bestSequential.StringSequence bestSequential.IsTarget bestSequential.Score
printfn "parallel best:   %-13s target=%b score=%.4f"
    bestParallel.StringSequence bestParallel.IsTarget bestParallel.Score

(**
```text
sequential best: ANLGMEVMHER   target=true score=13.0383
parallel best:   ANLGMEVMHER   target=true score=13.0383
```

## Where the scores go next

The library offers further scorers following a different scoring philosophy,
covered in
[Andromeda-like and X!Tandem-like scoring]({{root}}02_04_andromeda_like_scoring.html),
and the target and decoy scores collected across a run feed
[false discovery rate control]({{root}}04_01_fdr_control.html).
*)
