(**
---
title: Andromeda-like and X!Tandem-like scoring
category: Peptide identification
categoryindex: 2
index: 4
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

# Andromeda-like and X!Tandem-like scoring

The [SEQUEST-like scoring]({{root}}02_03_sequest_like_scoring.html) page
turned the agreement between a predicted and a measured spectrum into a dot
product of intensity vectors. The scorers on this page, `AndromedaLike` and
`XScoring`, ask a probability question.
Given that a candidate predicts n fragment peaks inside the measured range
and k of them coincide with a measured peak, how surprising is that under
random matching? The answer is a cumulative binomial probability, the chance
of seeing at least k hits among n tries when each try succeeds with a fixed
background probability. The score reports that probability as a -10 log10
value, so a higher score means the match is less likely to be chance. This
is the scoring idea of Andromeda, the search engine behind MaxQuant, and the
same matching machinery also yields an X!Tandem-style hyperscore.

Intensity enters through which measured peaks are allowed to match at all,
an idea called peak depth. `ratedSpectrum` rates every measured peak by
counting how many more intense peaks lie within a 100 Da window centered on
it. At depth q, only peaks that rank among the q most intense within their
own window are offered to the matcher, and q also sets the background
probability of a random hit to q/100 (capped at 0.5 in the implementation),
one window of width 100 Da holding q acceptable peaks. `countMatches` tallies n and k for every depth of a
user-given range in a single pass over the predictions. `scoreFuncImpl` then
turns each (n, k, q) triple into a score, and the best score across the
depths is kept. The trade-off behind the range: a small q admits only the
strongest peaks, so a random match is improbable, but genuine fragments of
modest intensity are missed. A large q admits more real fragments along with
more noise. Evaluating a range of depths lets every candidate be scored at
the depth that suits the spectrum best.

## Loading the measured spectrum

The measured side is the running example of these pages, the MS2 scan of the
doubly charged peptide ANLGMEVMHER from `ms2Example.mgf` with its precursor
at m/z 643.803548. The scorer locates matching peaks with a binary search,
so the measured spectrum must be sorted ascending by m/z, which this scan
is.
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

## Rebuilding the candidate list

The candidates are the same four isobaric sequences the
[SEQUEST-like page]({{root}}02_03_sequest_like_scoring.html) built and discussed in
detail. All four weigh 1285.5907 Da, so a mass window query cannot tell them apart.
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

printfn "%s" (candidates |> List.map (fun c -> c.StringSequence) |> String.concat ", ")

(**
```text
ANLGMEVMHER, MANGLEVMHER, EVMANLGMHER, AGGLGMEVMHER
```

## Predicting tagged fragment spectra

The fragment masses per candidate come from
`Fragmentation.Series.fragmentMasses` with the b and y series, exactly as on
the SEQUEST-like page. `AndromedaLike.getTheoSpecs` then converts every pair
into a `TheoreticalSpectrum` whose target and decoy sides are arrays of
`PeakFamily<TaggedPeak>`, predicted m/z values that keep their ion series
tags. The prediction model of this scorer family: a singly charged precursor
gets singly charged predicted fragments only, and a precursor of charge two
or more gets every fragment as a singly and a doubly charged predicted peak
(higher fragment charges are not emitted). The neutral-loss dependents ride
on the singly charged copies. The peaks carry no predicted
intensities, matching needs only the m/z and the tag. Restricting the
comparison to the scan limits happens in the matching step, which skips
predicted peaks outside the range.
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

let theoSpecs = AndromedaLike.getTheoSpecs scanlimits 2 fragmentPairs

theoSpecs
|> List.iter (fun ts ->
    let withLosses =
        ts.TheoSpec
        |> Array.filter (fun f -> not f.DependentPeaks.IsEmpty)
        |> Array.length
    printfn "%-13s peak families: %i  with loss dependents: %i"
        ts.LookUpResult.StringSequence ts.TheoSpec.Length withLosses)

(**
```text
AGGLGMEVMHER  peak families: 48  with loss dependents: 18
EVMANLGMHER   peak families: 44  with loss dependents: 22
MANGLEVMHER   peak families: 44  with loss dependents: 20
ANLGMEVMHER   peak families: 44  with loss dependents: 21
```

The 22 b and y fragments of an 11 residue candidate become 44 peak families,
one per charge state, and the extra residue of the Gly-Gly variant adds two
more fragments and four more families. The families that carry loss
dependents are the singly charged copies of fragments containing loss-prone
residues.

## Scoring at increasing peak depths

`AndromedaLike.calcAndromedaScore` takes the depth range, the scan limits,
the matching tolerance in ppm of the fragment m/z, the measured spectrum,
the scan time, the precursor charge and isolation window target m/z, the
theoretical spectra and a spectrum identifier. We evaluate depths 4 to 10 at
a matching tolerance of 100 ppm. The result is the familiar list of
`SearchEngineResult` records, one per target and decoy spectrum, ranked by
descending score with the two delta fields filled in, as described on the
[SEQUEST-like page]({{root}}02_03_sequest_like_scoring.html).
*)

let qMinAndMax = 4, 10

let andromedaResults =
    AndromedaLike.calcAndromedaScore
        qMinAndMax scanlimits 100. spectrum 20.93 2 precursorMZ theoSpecs "ms2Example"

let printRanked (rs: SearchEngineResult.SearchEngineResult<float> list) =
    printfn "%-13s %-6s %8s %12s %8s" "sequence" "target" "score" "dBestToRest" "dNext"
    rs
    |> List.iter (fun r ->
        printfn "%-13s %-6b %8.4f %12.4f %8.4f"
            r.StringSequence r.IsTarget r.Score r.NormDeltaBestToRest r.NormDeltaNext)

printRanked andromedaResults

(**
```text
sequence      target    score  dBestToRest    dNext
ANLGMEVMHER   true    86.4073       0.0000   0.2577
AGGLGMEVMHER  true    64.1383       0.2577   0.5796
MANGLEVMHER   true    14.0586       0.8373   0.0871
EVMANLGMHER   true     6.5320       0.9244   0.0756
ANLGMEVMHER   false    0.0000       1.0000   0.0000
MANGLEVMHER   false    0.0000       1.0000   0.0000
EVMANLGMHER   false    0.0000       1.0000   0.0000
AGGLGMEVMHER  false    0.0000       1.0000   0.0000
```

The target of ANLGMEVMHER wins at 86.4, and the Gly-Gly variant is again
the runner-up ahead of both permutations, the same shape the SEQUEST-like
ranking had. The probability view separates the top pair more clearly,
though. Where the dot product saw a near-tie with a `dNext` of 0.013, the
binomial score puts a quarter of the top score between them. The mechanism
is visible in the model: a predicted peak that finds no partner enlarges n
without enlarging k and lowers the score, so the extra predictions of the
12 residue variant, and the few fragments it does not share with the true
ladder, cost it directly. All four decoys land at exactly zero: their raw
scores fall below the fixed correction terms explained next, and the clamp
cuts them off.

## Reading an Andromeda-like score

The score is derived from -10 log10 of the cumulative binomial probability.
On that scale a value of 60 corresponds to a probability of 10^-6 under the
random-match model. The implementation adds a mass-dependent correction
computed from the precursor m/z plus constant modification and cleavage
correction terms taken from the original Andromeda release, then subtracts
100. Negative results are clamped to zero. Absolute values are therefore
calibrated for ranking within this implementation.

## Getting a hyperscore alongside

`XScoring.calcAndromedaAndXTandemScore` runs the same rating and matching
machinery once per candidate and returns a pair of result lists, an
Andromeda-like ranking whose records carry `SearchEngine = AndromedaLike`
and an X!Tandem-like ranking carrying `SearchEngine = XTandemLike`, each
independently score-sorted with its own delta fields. It accepts the same
arguments and the same theoretical spectra as `calcAndromedaScore`.

The X!Tandem score is the hyperscore: the summed intensity of all matched
measured peaks, taken as a logarithm, plus the log factorials of the number
of matched b ions and matched y ions, and a third log factorial for a separate
neutral-loss count. The factorials reward candidates that
match many ions of both series, since ten matched y ions weigh far more than
twice five. The counting is b/y-oriented, a matched peak contributes to the
count of the series flag it is tagged with, and since the loss peaks of this
library's generator carry their parent series flag, the separate neutral-loss
count stays empty here. The hyperscore is computed at
the deepest peak depth of the range, where the most measured peaks are
admitted to matching.
*)

let andromedaResults2, hyperscoreResults =
    XScoring.calcAndromedaAndXTandemScore
        qMinAndMax scanlimits 100. spectrum 20.93 2 precursorMZ theoSpecs "ms2Example"

printRanked hyperscoreResults

(**
```text
sequence      target    score  dBestToRest    dNext
ANLGMEVMHER   true    49.3979       0.0000   0.0994
AGGLGMEVMHER  true    44.4862       0.0994   0.5145
MANGLEVMHER   true    19.0709       0.6139   0.1745
EVMANLGMHER   true    10.4505       0.7884   0.1250
EVMANLGMHER   false    4.2767       0.9134   0.0082
ANLGMEVMHER   false    3.8712       0.9216   0.0000
MANGLEVMHER   false    3.8712       0.9216   0.0000
AGGLGMEVMHER  false    3.8712       0.9216   0.0000
```

Both rankings agree on the order of the four targets. The decoys behave differently
under the hyperscore: a decoy still matches a few peaks by chance, and the logarithm
of their summed intensity enters the score even when the factorial terms contribute
nothing, so the decoys settle at a nonzero floor around 4 where the corrected
Andromeda score clamped them to zero.

The Andromeda-like ranking returned alongside agrees with what
`AndromedaLike.calcAndromedaScore` produced on the same input, both modules
implement the same score.
*)

printfn "%-13s %-6s %14s %9s" "sequence" "target" "AndromedaLike" "XScoring"
List.zip andromedaResults andromedaResults2
|> List.iter (fun (a, x) ->
    printfn "%-13s %-6b %14.4f %9.4f" a.StringSequence a.IsTarget a.Score x.Score)

(**
```text
sequence      target  AndromedaLike  XScoring
ANLGMEVMHER   true          86.4073   86.4073
AGGLGMEVMHER  true          64.1383   64.1383
MANGLEVMHER   true          14.0586   14.0586
EVMANLGMHER   true           6.5320    6.5320
ANLGMEVMHER   false          0.0000    0.0000
MANGLEVMHER   false          0.0000    0.0000
EVMANLGMHER   false          0.0000    0.0000
AGGLGMEVMHER  false          0.0000    0.0000
```

The score columns are identical for all eight records, so results from the
combined call can stand in for a separate `AndromedaLike` run.

## Scoring whole runs

For batch searches over many spectra the library also offers
`SearchEngineGeneric.OrderedCache.generateTheoSpectra`. It wires the
database lookup and the spectrum predictors together: given the ion series
calculator, the mass function, a peptide lookup function and three caches
(one for lookup results, one for Andromeda-style and one for SEQUEST-style
theoretical spectra), it answers a precursor mass window with the
theoretical spectra for both engine families, keeping previously generated
spectra in memory across consecutive windows and clearing the caches when a
memory ceiling is exceeded. The per-engine `getTheoSpecs` calls shown on
this page and the SEQUEST-like page are the direct route for scoring
individual spectra.

## Where the scores go next

The target and decoy scores and their deltas collected across a run feed
[false discovery rate control]({{root}}04_01_fdr_control.html), and the
identified peptides move on to
[quantification]({{root}}03_01_quantification.html).
*)
