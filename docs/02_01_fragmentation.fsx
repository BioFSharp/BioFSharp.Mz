(**
---
title: In silico fragmentation
category: Peptide identification
categoryindex: 2
index: 1
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

# In silico fragmentation

An MS2 spectrum is the fragment spectrum of one selected peptide ion. Inside the
collision cell the peptide collides with gas molecules until its backbone breaks,
and the instrument records the m/z values of the resulting pieces. Reading the
sequence directly out of those peaks is hard. Search engines therefore go the
other way: take a candidate sequence and predict which fragment masses it would
produce, then check how many of the predictions show up in the measured
spectrum. The candidate whose predictions explain the spectrum best wins. This
page covers the prediction side, implemented in the `Fragmentation` module.

Under collision-induced dissociation (CID) and electron transfer dissociation (ETD)
the backbone cleaves between two residues, leaving an N-terminal piece and a
C-terminal piece. The Roepstorff-Fohlman-Biemann nomenclature names the
fragments by the exact bond that broke: N-terminal fragments are called a, b or
c ions, C-terminal fragments x, y or z ions. The index counts residues from the
fragment's own terminus, so b2 covers the first two residues and y1 is the last
residue alone. CID, the most common fragmentation method and the plausible origin
of our example spectrum, produces mostly b and y ions, which is why those two
series carry most identifications.

On top of the backbone fragments come satellite peaks. A fragment containing
arginine, lysine, asparagine or glutamine can additionally lose ammonia (17 Da),
and a fragment containing serine, threonine, glutamate or aspartate can lose
water (18 Da). These losses appear as smaller peaks slightly below the main
fragment peak. The generator knows both rules, the residue sets live in the
module as `Fragmentation.aminoLossSet` and `Fragmentation.waterLossSet`.

`ms2Example.mgf` holds a measured MS2 scan of the doubly charged peptide ANLGMEVMHER, and the
[charge state determination]({{root}}01_03_charge_state_determination.html) page
recovered its neutral mass. Here we predict the fragments of ANLGMEVMHER and
show that they actually sit in that spectrum.

## Building the peptide

The fragmentation functions work on an amino acid list from BioFSharp.
`BioList.ofAminoAcidString` parses the one letter sequence. Every `Series`
function additionally takes a mass function of type `IBioItem -> float`, which
is how the caller decides whether the ladder is computed from monoisotopic or
average masses. We pass `BioItem.monoisoMass`, since fragment matching at 0.05
Da tolerance needs monoisotopic values.
*)

open BioFSharp
open BioFSharp.Mz

let peptide = BioList.ofAminoAcidString "ANLGMEVMHER"

let mono : IBioItem -> float = BioItem.monoisoMass

printfn "residues: %i" (List.length peptide)
printfn "monoisotopic mass of Ala: %f" (mono AminoAcids.Ala)

(**
```text
residues: 11
monoisotopic mass of Ala: 71.037114
```

## Generating the b series

`Fragmentation.Series.bOfBioList` walks the sequence from the N-terminus and
accumulates residue masses. Each element of the result is a
`PeakFamily<TaggedMass>`, the shape introduced on the
[peaks and peak arrays]({{root}}01_01_peaks_and_peak_arrays.html) page: a main
peak tagged with its ion series flag plus a list of dependent loss peaks. The
masses are neutral fragment masses, no proton has been added yet. A small
printer makes the ladder readable.
*)

let printLadder (prefix: string) (families: PeakFamily<TaggedMass.TaggedMass> list) =
    families
    |> List.iteri (fun i family ->
        let losses =
            family.DependentPeaks
            |> List.map (fun d -> sprintf "%A at %.5f" d.Iontype d.Mass)
            |> String.concat "   "
        printfn "%s%-2i %A %10.5f   %s" prefix (i + 1) family.MainPeak.Iontype family.MainPeak.Mass losses)

let bSeries = Fragmentation.Series.bOfBioList mono peptide

printfn "b families: %i" bSeries.Length
printLadder "b" bSeries

(**
```text
b families: 11
b1  B   71.03711
b2  B  185.08004   B, lossNH3 at 168.05349
b3  B  298.16411   B, lossNH3 at 281.13756
b4  B  355.18557   B, lossNH3 at 338.15902
b5  B  486.22605   B, lossNH3 at 469.19950
b6  B  615.26865   B, lossNH3 at 598.24210   B, lossH2O at 597.25808
b7  B  714.33706   B, lossNH3 at 697.31051   B, lossH2O at 696.32650
b8  B  845.37755   B, lossNH3 at 828.35100   B, lossH2O at 827.36698
b9  B  982.43646   B, lossNH3 at 965.40991   B, lossH2O at 964.42589
b10 B 1111.47905   B, lossNH3 at 1094.45250   B, lossH2O at 1093.46849
b11 B 1267.58016   B, lossNH3 at 1250.55361   B, lossH2O at 1249.56960
```

The numbers check out by hand. b2 covers alanine and asparagine, and their
residue masses sum to 71.03711 + 114.04293 = 185.08004, exactly the printed
neutral mass. Protonation adds 1.00728, giving the singly charged ion at m/z
186.08732 that the peaks page constructed manually.

The generator emits one family per residue, so the 11 residue peptide yields 11 b
families, and the last one spans the full sequence. The loss dependents
are cumulative: once the ladder has passed a loss prone residue, every longer
fragment carries that loss peak. The asparagine at position 2 puts an ammonia
loss on b2 and on every later b ion, and the glutamate at position 6 adds a
water loss from b6 on. b1 has no dependents because alanine triggers neither
rule.

## Generating the y series

`Fragmentation.Series.yOfBioList` builds the C-terminal ladder the same way,
except that the accumulator starts with the mass of one water molecule, because
a y fragment keeps the C-terminal hydroxyl and picks up a hydrogen at the break.
The function returns the ladder longest first, the head of the list is the full
length y11. For display we reverse it so that y1 comes first, and the code shows
that reordering explicitly.
*)

let ySeries = Fragmentation.Series.yOfBioList mono peptide

printfn "y families: %i" ySeries.Length
printfn "mass of the first returned family: %.5f" ySeries.Head.MainPeak.Mass

let yAscending = ySeries |> List.rev

printLadder "y" yAscending

(**
```text
y families: 11
mass of the first returned family: 1285.59073
y1  Y  174.11168   Y, lossNH3 at 157.08513
y2  Y  303.15427   Y, lossNH3 at 286.12772   Y, lossH2O at 285.14370
y3  Y  440.21318   Y, lossNH3 at 423.18663   Y, lossH2O at 422.20262
y4  Y  571.25367   Y, lossNH3 at 554.22712   Y, lossH2O at 553.24310
y5  Y  670.32208   Y, lossNH3 at 653.29553   Y, lossH2O at 652.31151
y6  Y  799.36467   Y, lossNH3 at 782.33812   Y, lossH2O at 781.35411
y7  Y  930.40516   Y, lossNH3 at 913.37861   Y, lossH2O at 912.39459
y8  Y  987.42662   Y, lossNH3 at 970.40007   Y, lossH2O at 969.41606
y9  Y 1100.51069   Y, lossNH3 at 1083.48414   Y, lossH2O at 1082.50012
y10 Y 1214.55361   Y, lossNH3 at 1197.52706   Y, lossH2O at 1196.54305
y11 Y 1285.59073   Y, lossNH3 at 1268.56418   Y, lossH2O at 1267.58016
```

The cumulative loss rule reads differently from this side because the ladder
grows from the C-terminus. The C-terminal arginine puts an ammonia loss on y1
and with it on every y ion. The glutamate one position further in adds a water
loss from y2 on. y11 is the intact peptide plus water, and its 1285.59073 Da
match the 1285.59 Da neutral precursor mass determined on the charge state
page.

The module offers the remaining series of the nomenclature through the same
pattern, `aOfBioList`, `cOfBioList`, `xOfBioList` and `zOfBioList`, plus
combined variants such as `abOfBioList` or `yzOfBioList` that return several
series in one list.

## Matching the predictions against the measured spectrum

The point of all this is that predicted fragments of the correct peptide should
be findable in the measured spectrum. We load the MS2 scan and look up every
predicted b and y ion in its m/z array as a singly charged species.
`Mass.toMZ` from BioFSharp does the charge conversion, it adds one proton mass
of 1.00728 per charge and divides by the charge. A predicted ion counts as
matched when a measured peak lies within 0.05 Da of it, and we keep the closest
such peak.
*)

open BioFSharp.FileFormats.MGF
open BioFSharp.IO

let ms2 =
    MGF.read (__SOURCE_DIRECTORY__ + "/data/ms2Example.mgf")
    |> List.head

let measuredMz = ms2.Mass

let predictions =
    (bSeries |> List.mapi (fun i f -> sprintf "b%i" (i + 1), f.MainPeak.Mass))
    @ (yAscending |> List.mapi (fun i f -> sprintf "y%i" (i + 1), f.MainPeak.Mass))

let matched =
    predictions
    |> List.choose (fun (label, neutralMass) ->
        let predictedMz = Mass.toMZ neutralMass 1.
        let close = measuredMz |> Array.filter (fun mz -> abs (mz - predictedMz) <= 0.05)
        if Array.isEmpty close then None
        else Some (label, predictedMz, close |> Array.minBy (fun mz -> abs (mz - predictedMz))))

printfn "matched %i of %i predicted singly charged ions" matched.Length predictions.Length
printfn "ion  predicted   measured"
matched
|> List.iter (fun (label, predicted, measured) ->
    printfn "%-4s %9.4f  %9.4f" label predicted measured)

(**
```text
matched 18 of 22 predicted singly charged ions
ion  predicted   measured
b2    186.0873   186.0879
b3    299.1714   299.1703
b4    356.1928   356.1917
b5    487.2333   487.2292
b6    616.2759   616.2792
b7    715.3443   715.3339
b9    983.4437   983.4372
b10  1112.4863  1112.5331
y1    175.1190   175.1183
y2    304.1615   304.1598
y3    441.2205   441.2191
y4    572.2609   572.2601
y5    671.3294   671.3284
y6    800.3719   800.3729
y7    931.4124   931.4112
y8    988.4339   988.4343
y9   1101.5180  1101.5180
y10  1215.5609  1215.5282
```

18 of the 22 predictions find a measured peak, covering the y ladder from y1 to
y10 and most of the b ladder. The misses are plausible as well. b1 at m/z 72.04
lies below the first recorded peak of the scan at m/z 100.67, and the other
three predictions simply have no peak within the tolerance. For a wrong
candidate sequence only a few predictions would hit a peak by chance. Turning
this observation into a score is the job of the scoring pages.

## Generating a target and decoy pair

Search engines need to know what a random match looks like, otherwise there is
no way to tell a convincing score from a lucky one. The standard trick is to
score every spectrum against the real candidate (the target) and additionally
against a deliberately wrong candidate of the same length and composition (the
decoy), obtained by reversing the sequence. The decoy scores calibrate the
score distribution of random matches, which later drives
[false discovery rate control]({{root}}04_01_fdr_control.html).

`Fragmentation.Series.fragmentMasses` produces both at once. It takes the
N-terminal series function, the C-terminal series function, the mass function
and the peptide, and returns a `FragmentMasses` record with `TargetMasses` for
the given sequence and `DecoyMasses` computed from the reversed sequence.
*)

let fragments =
    Fragmentation.Series.fragmentMasses
        Fragmentation.Series.bOfBioList
        Fragmentation.Series.yOfBioList
        mono
        peptide

printfn "target families: %i" fragments.TargetMasses.Length
printfn "decoy families:  %i" fragments.DecoyMasses.Length

fragments.DecoyMasses
|> List.truncate 3
|> List.iter (fun f -> printfn "decoy %A ion at %.5f" f.MainPeak.Iontype f.MainPeak.Mass)

(**
```text
target families: 22
decoy families:  22
decoy B ion at 156.10111
decoy B ion at 285.14370
decoy B ion at 422.20262
```

Both lists concatenate the b and the y families, 22 entries each for our
peptide. The first decoy b ion weighs 156.10111 Da, the residue mass of
arginine, because the reversed sequence REHMVEMGLNA starts with the arginine
that ANLGMEVMHER ends with.

## Converting neutral masses to charged m/z ladders

The families so far carry neutral masses, while a spectrum shows charged ions.
In a real MS2 scan a fragment can appear at more than one charge when the
precursor charge allows it, our 2+ precursor can hand both protons to one
fragment, so the doubly charged ladder is worth predicting as well.
`Fragmentation.ladderElement` takes a peak family list and a charge list and
returns one family per input family and charge. Each resulting
`LadderedTaggedMass` carries the ion flag, the `MassOverCharge` value, the
`Number` of the ion within its series and the `Charge`. Internally the function
groups the input by ion type and sorts each group by mass, which is how the
number is assigned, so the mixed b and y list can go in as is.
*)

let laddered = Fragmentation.ladderElement (bSeries @ ySeries) [1.; 2.]

printfn "laddered families: %i" laddered.Length

laddered
|> List.filter (fun f -> f.MainPeak.Number <= 2)
|> List.iter (fun f ->
    let m = f.MainPeak
    printfn "%A ion %i at charge %.0f: m/z %9.4f" m.Iontype m.Number m.Charge m.MassOverCharge)

(**
```text
laddered families: 44
B ion 1 at charge 1: m/z   72.0444
B ion 1 at charge 2: m/z   36.5258
B ion 2 at charge 1: m/z  186.0873
B ion 2 at charge 2: m/z   93.5473
Y ion 1 at charge 1: m/z  175.1190
Y ion 1 at charge 2: m/z   88.0631
Y ion 2 at charge 1: m/z  304.1615
Y ion 2 at charge 2: m/z  152.5844
```

The 22 input families become 44 laddered ones, one per charge. The 1+ values repeat
what the matching section used, b2 at m/z 186.0873 for example, and dependent loss
peaks are laddered alongside their main peak with the same number and charge.

## Where the fragments go next

The `TheoreticalSpectra` module pairs each candidate peptide coming out of a database
search with its predicted target and decoy spectrum, and a scorer such as
[SEQUEST-like scoring]({{root}}02_03_sequest_like_scoring.html) then quantifies how
well the measured peaks agree with each prediction.
*)
