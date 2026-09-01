(**
---
title: Peaks and peak arrays
category: Spectrum processing
categoryindex: 1
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

# Peaks and peak arrays

A mass spectrometer separates ions by their mass-to-charge ratio (m/z) and records how
many of them arrive at the detector. A peak is the recorded signal of one ion species.
In raw profile data that signal is a bump spanning many measured points around the
ion's true m/z. Centroiding reduces each such bump to a single representative pair of
an m/z value and an intensity, a step covered on the
[signal detection]({{root}}01_02_signal_detection.html) page. The `Peak` type in
BioFSharp.Mz models exactly this centroided form, and it is the smallest unit
everything else is built from. Fragment matching, database search, spectrum scoring and
similarity comparisons all consume or produce peaks.

This page loads a measured MS2 spectrum from an mgf file, represents it as a `PeakArray`,
introduces tagged peaks and peak families (the vocabulary of fragment ion matching), and
finally converts the spectrum into the binned and sparse forms that spectrum scoring
works on.

## Reading a spectrum from an mgf file

Mascot generic format (mgf) is a plain text format for MS2 data. Every entry holds a
peak list plus metadata such as the precursor m/z (the m/z of the intact peptide ion
that was selected for fragmentation) and its charge state. BioFSharp ships a reader in
`BioFSharp.IO.MGF`, and the entry type with its metadata accessors lives in
`BioFSharp.FileFormats.MGF`. The example file contains a single MS2 scan of the doubly
charged peptide ANLGMEVMHER.
*)

open BioFSharp.FileFormats.MGF
open BioFSharp.IO
open BioFSharp.Mz

let ms2 =
    MGF.read (__SOURCE_DIRECTORY__ + "/data/ms2Example.mgf")
    |> List.head

printfn "%s" (MGFEntry.tryGetTitle ms2 |> Option.defaultValue "no title")

match MGFEntry.tryGetPrecursorMZ ms2, MGFEntry.tryGetPrecursorCharges ms2 with
| Some mz, Some charges -> printfn "precursor m/z %f at charge %A" mz charges
| _ -> printfn "no precursor information"

(**
```text
MS/MS scan at 20.93388333 min with Intensity: 501.0 and Sequence ANLGMEVMHER
precursor m/z 643.803548 at charge [2]
```

The reader hands back the m/z and intensity values as two separate arrays. `PeakArray.zip`
pairs them up into one `PeakArray`, an array of `Peak` values. `Peak` is a small struct
carrying `Mz` and `Intensity`, and it implements the `IPeak` interface that all peak-like
types in the library share.
*)

let spectrum : PeakArray<Peak> = PeakArray.zip ms2.Mass ms2.Intensity

printfn "peak count: %i" spectrum.Length

spectrum
|> Array.truncate 3
|> Array.iter (fun p -> printfn "m/z %f  intensity %f" p.Mz p.Intensity)

(**
```text
peak count: 971
m/z 100.666675  intensity 10.000000
m/z 103.057515  intensity 12.000000
m/z 103.068913  intensity 12.000000
```

Some of the 971 entries are adjacent samples of the same ion signal, so the list is denser than a fully centroided spectrum, a distinction the
[signal detection]({{root}}01_02_signal_detection.html) page develops. The low
intensity values are typical for a single MS2 scan of a low abundance precursor.

## Creating and transforming peaks

Single peaks are created with `Peaks.createPeak`. Because `PeakArray` is a plain array
under the hood, the familiar `Array` functions work on it, and `PeakArray.map` builds a
new peak array from an existing one. A common use is normalizing all intensities to the
most intense peak of the spectrum (the base peak), which makes spectra of different
absolute intensity comparable.
*)

let singlePeak = Peaks.createPeak 175.119 24.
printfn "m/z %f  intensity %f" singlePeak.Mz singlePeak.Intensity

let basePeakIntensity =
    spectrum |> Array.map (fun p -> p.Intensity) |> Array.max

let normalized =
    spectrum
    |> PeakArray.map (fun p -> Peaks.createPeak p.Mz (p.Intensity / basePeakIntensity))

normalized
|> Array.truncate 3
|> Array.iter (fun p -> printfn "m/z %f  relative intensity %f" p.Mz p.Intensity)

(**
```text
m/z 175.119000  intensity 24.000000
m/z 100.666675  relative intensity 0.058824
m/z 103.057515  relative intensity 0.070588
m/z 103.068913  relative intensity 0.070588
```

`PeakList` is the list-backed sibling of `PeakArray` with the same `zip` and `map`
functions. When you need the raw m/z and intensity sequences back, for example
to feed a plotting library, `PeakList.unzip` splits a peak list into its two component
lists.
*)

let mzList, intensityList =
    spectrum |> List.ofArray |> PeakList.unzip

printfn "m/z values:       %A" (mzList |> List.truncate 3)
printfn "intensity values: %A" (intensityList |> List.truncate 3)

(**
```text
m/z values:       [100.666675; 103.057515; 103.068913]
intensity values: [10.0; 12.0; 12.0]
```

## Tagging peaks with ion types

When a peptide is fragmented in the mass spectrometer, its backbone breaks at
characteristic positions. A fragment that keeps the N-terminus is called a b ion, one
that keeps the C-terminus a y ion (there are more series, a, c, x and z). A peak on its
own does not know which fragment it belongs to. For fragment matching the library
attaches that information as a tag: `Ions.IonTypeFlag` is a flags enum naming the ion
series, `TaggedMass` pairs a flag with a mass, and `TaggedPeak.TaggedPeak` pairs a flag
with an m/z and an intensity.

The tags become concrete with two fragments of ANLGMEVMHER, computed by hand from the
monoisotopic residue masses. The b2 ion covers the first two residues, the y1 ion is the
C-terminal arginine.
*)

// b2 of ANLGMEVMHER: A + N + proton = 71.03711 + 114.04293 + 1.00728
let b2 = TaggedMass.createTaggedMass Ions.IonTypeFlag.B 186.08732

// y1 of ANLGMEVMHER: R + H2O + proton = 156.10111 + 18.01056 + 1.00728
let y1 = TaggedMass.createTaggedMass Ions.IonTypeFlag.Y 175.11895

printfn "%A ion at m/z %f" b2.Iontype b2.Mass
printfn "%A ion at m/z %f" y1.Iontype y1.Mass

(**
```text
B ion at m/z 186.087320
Y ion at m/z 175.118950
```

Fragments often show up together with satellite peaks, for example the same ion minus a
water molecule (18.01056 Da lighter). `TaggedMass.createTaggedH2OLoss` builds such a
loss tag by combining the ion series flag with the `lossH2O` flag, and
`Peaks.createPeakFamily` groups a main peak with its dependent satellite peaks into a
`PeakFamily`. The [in silico fragmentation]({{root}}02_01_fragmentation.html) page
produces exactly this shape, `PeakFamily<TaggedMass>` values for every fragment of a
peptide, so the type is introduced here on a small example.
*)

let y1WaterLoss =
    TaggedMass.createTaggedH2OLoss Ions.IonTypeFlag.Y (175.11895 - 18.01056)

let y1Family = Peaks.createPeakFamily y1 [ y1WaterLoss ]

printfn "main peak:      %A at m/z %f" y1Family.MainPeak.Iontype y1Family.MainPeak.Mass
y1Family.DependentPeaks
|> List.iter (fun t -> printfn "dependent peak: %A at m/z %f" t.Iontype t.Mass)

(**
```text
main peak:      Y at m/z 175.118950
dependent peak: Y, lossH2O at m/z 157.108390
```

The dependent peak carries both flags at once, which is what a flags enum is for. A
matching function can later ask for the series with `Ions.hasFlag`.

Because ANLGMEVMHER is the peptide that was actually fragmented in our example scan, its
y1 ion should be present in the measured spectrum. It is.
*)

spectrum
|> Array.filter (fun p -> abs (p.Mz - y1.Mass) < 0.05)
|> Array.iter (fun p -> printfn "measured peak near y1: m/z %f  intensity %f" p.Mz p.Intensity)

(**
```text
measured peak near y1: m/z 175.118252  intensity 12.000000
measured peak near y1: m/z 175.120109  intensity 12.000000
measured peak near y1: m/z 175.121966  intensity 24.000000
```

A `TaggedPeak.TaggedPeak` adds an intensity to the tag and implements `IPeak` like every
other peak type, so tagged peaks fit into a `PeakArray` as well.
*)

let y1Peak = TaggedPeak.createTaggedPeak Ions.IonTypeFlag.Y 175.11895 24.
printfn "%A peak at m/z %f with intensity %f" y1Peak.Iontype y1Peak.Mz y1Peak.Intensity

(**
```text
Y peak at m/z 175.118950 with intensity 24.000000
```

## Binning a spectrum into unit dalton bins

SEQUEST-style scoring compares a measured spectrum against a predicted one as vectors.
For that the continuous m/z axis is discretized into bins of 1 Da width, and every peak
is assigned to its nearest bin. `PeakArray.peaksToNearestUnitDaltonBinVector` performs
this conversion between a lower and an upper mass border and returns an `FSharp.Stats`
vector. When several peaks fall into the same bin, the bin keeps the highest intensity.
The produced vector has length `maxMassBoarder - minMassBoarder`, here 900 bins covering
m/z 100 to 1000.
*)

open FSharp.Stats

let binned = PeakArray.peaksToNearestUnitDaltonBinVector spectrum 100.0 1000.0

printfn "vector length: %i" (Vector.length binned)

let occupied =
    binned
    |> Vector.toArray
    |> Array.indexed
    |> Array.filter (fun (_, intensity) -> intensity > 0.)

printfn "occupied bins: %i" occupied.Length

occupied
|> Array.truncate 5
|> Array.iter (fun (i, intensity) -> printfn "bin %i (m/z %i): %f" i (i + 100) intensity)

(**
```text
vector length: 900
occupied bins: 378
bin 1 (m/z 101): 10.000000
bin 3 (m/z 103): 12.000000
bin 4 (m/z 104): 12.000000
bin 10 (m/z 110): 24.000000
bin 12 (m/z 112): 10.000000
```

971 peaks collapse into 378 occupied bins because peaks closer together than 1 Da share
a bin. The bin index is simply the rounded m/z minus the lower mass border, so bin 0
collects the peaks around m/z 100.

## Comparing spectra with sparse peak arrays

Most of the 900 bins above are zero. `SparsePeakArray` stores only the occupied bins in
a dictionary from bin index to intensity, which saves memory and makes comparing two
spectra a walk over one dictionary with lookups in the other.
`SparsePeakArray.peaksToNearestBinVector` builds one from a peak array. The first two
arguments control the binning: with a bin width of 1.0 and an offset of 0.5 every peak
lands in its nearest unit bin. In the sparse form, peaks sharing a bin are summed. The
record also carries the two conversion functions `MzToBinIdx` and `BinIdxToMz` so you
can move between m/z values and bin indices.
*)

let sparse =
    spectrum
    |> SparsePeakArray.peaksToNearestBinVector 1.0 0.5 100.0 1000.0

printfn "occupied bins: %i" sparse.Data.Count
printfn "bin index of the y1 m/z: %i" (sparse.MzToBinIdx 175.11895)

(**
```text
occupied bins: 378
bin index of the y1 m/z: 175
```

`SparsePeakArray.dot` computes the dot product of two sparse spectra: for every bin both
spectra occupy, it multiplies the intensities and sums the products. The dot product
(usually after normalization) is the standard similarity measure between binned spectra,
and a spectrum compared with itself gives the maximal value.

The example file contains only one scan, so we simulate a spectrum of an unrelated
peptide by shifting every measured m/z by 7.3 Da. The shifted copy has the same
intensity distribution, only in different bins.
*)

// same peaks moved by 7.3 Da: a stand-in for a spectrum of an unrelated peptide
let shifted =
    spectrum
    |> PeakArray.map (fun p -> Peaks.createPeak (p.Mz + 7.3) p.Intensity)

let sparseShifted =
    shifted
    |> SparsePeakArray.peaksToNearestBinVector 1.0 0.5 100.0 1000.0

printfn "self  dot product: %f" (SparsePeakArray.dot sparse sparse)
printfn "cross dot product: %f" (SparsePeakArray.dot sparse sparseShifted)

(**
```text
self  dot product: 3692147.000000
cross dot product: 163836.000000
```

The spectrum agrees with itself far better than with the shifted copy, exactly what a
similarity measure should say. The remaining cross product comes from bins where a
shifted peak happens to land on another real peak.

The next step is generating predicted fragment peaks to match against, the topic of the
[in silico fragmentation]({{root}}02_01_fragmentation.html) page.
*)
