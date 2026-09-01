(**
---
title: Charge state determination
category: Spectrum processing
categoryindex: 1
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

# Charge state determination

A mass spectrometer reports m/z, the mass of an ion divided by the number of charges
it carries. The database search at the end of the pipeline works with neutral peptide
masses, so at some point the charge z has to be inferred. The stakes are concrete: a precursor observed at m/z 643.8 weighs
about 1285.6 Da if it carries two charges and about 1928.4 Da if it carries three.
Search the wrong mass range and the peptide is never found.

The information needed to infer z is already in the MS1 scan. Elements come in
naturally occurring isotopes, and with roughly one in a hundred carbon atoms being a
heavy ¹³C, a population of identical peptide molecules spreads out into an isotope
cluster: the monoisotopic peak, a peak one neutron heavier, one two neutrons heavier,
and so on. For a peptide with neutral mass M at charge z, the n-th isotope peak sits
at m/z = (M + z·H + n·1.00335) / z. Adjacent cluster peaks are therefore spaced
1.00335 / z apart on the m/z axis. A cluster whose peaks sit 1.0 apart is singly
charged, 0.5 apart means 2+, 0.33 apart means 3+, 0.25 apart means 4+. The
`ChargeState` module reads the
charge out of exactly this spacing, using 1 / z as an approximation of the
theoretical distance.

`ms2Example.mgf` holds a fragment
spectrum of the peptide ANLGMEVMHER, recorded at 20.93 min. `ms1Example.mgf` is the
MS1 survey scan recorded right before it, at 20.91 min in the same run, where the
fragmented precursor was measured intact, so its isotope cluster is sitting in that
scan's data. We determine the charge from the MS1 scan and check the resulting mass
against what the MS2 file states in its header.
*)

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

printfn "%s" (MGFEntry.tryGetTitle ms2 |> Option.defaultValue "no title")
printfn "precursor m/z: %f" precursorMZ

(**
```text
MS/MS scan at 20.93388333 min with Intensity: 501.0 and Sequence ANLGMEVMHER
precursor m/z: 643.803548
```

The header also carries a `PEPMASS` of 1285.591258 Da. That value is our ground
truth. If the charge determination works, we arrive at this mass using only the MS1
peaks and the precursor m/z.

## Centroiding the MS1 survey scan

Charge determination runs on centroided peaks, one peak per ion signal. We prepare
the MS1 scan exactly as on the
[signal detection]({{root}}01_02_signal_detection.html) page, with the same padding
and wavelet parameters, and refer there for what each of them does.
*)

let ms1 =
    MGF.read (__SOURCE_DIRECTORY__ + "/data/ms1Example.mgf")
    |> List.head

let paddingParams =
    SignalDetection.Padding.createPaddingParameters 0. None 0.05 150 95.

let waveletParams =
    SignalDetection.Wavelet.createWaveletParameters 3 0. 0.1 95. 1. false false

let centroidMz, centroidIntensity =
    SignalDetection.Padding.paddDataBy paddingParams ms1.Mass ms1.Intensity
    ||> SignalDetection.Wavelet.toCentroidWithRicker2D waveletParams

printfn "centroids: %i" centroidMz.Length

(**
```text
centroids: 1454
```

Before handing these arrays to the module, we look at the centroids in a small
window starting at the precursor m/z, together with their spacings.
*)

let clusterWindow =
    Array.zip centroidMz centroidIntensity
    |> Array.filter (fun (mz, _) -> mz >= 643.5 && mz <= 645.5)

clusterWindow
|> Array.iter (fun (mz, i) -> printfn "m/z %f  intensity %5.1f" mz i)

clusterWindow
|> Array.map fst
|> Array.pairwise
|> Array.iter (fun (a, b) -> printfn "spacing: %f" (b - a))

(**
```text
m/z 643.802071  intensity 501.0
m/z 644.311347  intensity 352.0
m/z 644.810133  intensity 150.0
m/z 645.305547  intensity 124.0
spacing: 0.509276
spacing: 0.498786
spacing: 0.495414
```

A falling intensity ladder spaced at about half a dalton, which by the formula above
is the signature of a doubly charged ion. The first peak's intensity of 501.0 is the
precursor intensity the MS2 title records, so we are looking at the right signal.
Reading the charge off by eye works here. The module does the same thing
systematically and with safeguards for messier neighborhoods.

## Generating charge state candidates

`ChargeState.putativePrecursorChargeStatesBy` needs a `ChargeDetermParams` record,
built with `createChargeDetermParams`. Its parameters, in order: the expected minimal
and maximal charge span the states that are considered at all. Tryptic peptides
mostly ionize at 2+ to 4+, one proton at the N-terminus and one at the C-terminal
lysine or arginine, sometimes more on internal basic residues, so 2 and 4 are
reasonable bounds. `Width` is the m/z window, in daltons, scanned to the right of the
precursor peak for cluster members. 1.1 covers one full isotope spacing at charge 1
and correspondingly more peaks at higher charges. The next two values decide which
peaks in that window count as cluster members: a peak must reach `MinIntensity`
(here 15 percent) of the start peak's intensity and `DeltaMinIntensity` (here 30
percent) of the peak before it, which cuts the cluster off once it has decayed into
noise. `NrOfRndSpectra` is the number of random clusters simulated for the hypothesis
test in the last section. The module author's design uses 10000.
*)

let chargeParams =
    ChargeState.createChargeDetermParams 2 4 1.1 0.15 0.3 10000

(**
The function locates the centroid closest to the given precursor m/z, collects the
window peaks that pass the two intensity criteria, and then considers every subset of
them that contains the start peak, because chemical noise can sit between real
isotope peaks and a fixed left-to-right reading would stumble over it. Each subset
gets the charge from the allowed range whose theoretical spacing best matches the
subset's mean peak spacing, a deviation measure (`MZChargeDev`) quantifying how well
the spacings fit that charge, and a score, the weighted sum of that deviation and a
length quotient, the fraction of window peaks the subset uses (weights from a linear
discriminant analysis, lower is better). The neutral mass `PutMass` is computed from the precursor m/z and the
assigned charge with `Mass.ofMZ` from BioFSharp. The result keeps the best subset per
distinct charge, sorted by score. One safeguard to know about: a window holding 15 or
more peaks is considered too noisy for spacing analysis, and the function then falls
back to enumerating every allowed charge instead.

The two string arguments are free-text identifiers for the MS1 and MS2 scan. They are
carried into the result record so a downstream consumer can trace which spectra an
assignment belongs to.
*)

let printCandidate (c: ChargeState.AssignedCharge) =
    printfn "charge %i+  score %8.4f  mzChargeDev %8.4f  putMass %9.4f  peaks %i"
        c.PrecCharge c.Score c.MZChargeDev c.PutMass c.SubSetLength

let putativeCharges =
    ChargeState.putativePrecursorChargeStatesBy
        chargeParams centroidMz centroidIntensity
        "MS1 at 20.91 min" "MS2 at 20.93 min" precursorMZ

putativeCharges |> List.iter printCandidate

(**
```text
charge 2+  score  -0.2906  mzChargeDev   0.0094  putMass 1285.5925  peaks 3
```

One candidate survives the subset competition: charge 2, backed by a cluster of three
peaks (the fourth ladder peak from above lies beyond the 1.1 Da window). The putative
mass of 1285.5925 Da lands within about a millidalton of the 1285.591258 Da stated in
the MS2 header, the small residual coming down to which mass constant is used for the
added protons. The mass the database search needs has been recovered from the MS1
peak pattern alone.

## Selecting among candidates

`ChargeState.removeSubSetsOfBestHit` takes the first list element
as the best hit, which the score-sorted output of `putativePrecursorChargeStatesBy`
provides, and drops every other candidate whose peak set is a subset of the best
hit's peaks. Such candidates explain no signal the best hit does not already explain.
Candidates built from different peaks survive, since they may describe a second,
overlapping cluster.
*)

let assignedCharges =
    ChargeState.removeSubSetsOfBestHit putativeCharges

assignedCharges |> List.iter printCandidate

(**
```text
charge 2+  score  -0.2906  mzChargeDev   0.0094  putMass 1285.5925  peaks 3
```

With a single candidate there is nothing to remove. The filter becomes visible on a
second cluster.

## Cross-checking on the cluster at m/z 435.19

The [signal detection]({{root}}01_02_signal_detection.html) page ended on the most
intense signal of this scan, four centroids around m/z 435.19 spaced about 0.334 Da
apart, and read a charge of 3 off that spacing by eye. The same call reproduces this
on a second, independent case.
*)

let putativeCharges435 =
    ChargeState.putativePrecursorChargeStatesBy
        chargeParams centroidMz centroidIntensity
        "MS1 at 20.91 min" "-" 435.189052

putativeCharges435 |> List.iter printCandidate

(**
```text
charge 3+  score  -0.2950  mzChargeDev   0.0050  putMass 1302.5453  peaks 3
charge 2+  score   0.0372  mzChargeDev   0.2372  putMass  868.3636  peaks 2
```

The best candidate is the expected 3+, built from three cluster peaks. The 2+
candidate is an artifact of the subset construction: the start peak paired with the
second isotope peak alone spans 0.668 Da, read as one poor charge-2 gap, which both
its deviation and its score reflect. Its two peaks are contained in the best hit's
three, so the subset filter disposes of it.
*)

ChargeState.removeSubSetsOfBestHit putativeCharges435
|> List.iter printCandidate

(**
```text
charge 3+  score  -0.2950  mzChargeDev   0.0050  putMass 1302.5453  peaks 3
```

## Testing the cluster against random spacings

`MZChargeDev` is a unitless deviation, and on its own it is hard to say how small is
small enough. The module frames the question as a hypothesis test. The null
hypothesis says the cluster is a random arrangement of legal charge spacings. When
random arrangements practically never fit the assigned charge as well as the observed
cluster does, the null is rejected and the cluster counts as a genuine isotope
pattern. First, `peakPosStdDevBy` estimates the spacing scatter of the current
measurement from the candidates' deviations between real and theoretical peak
spacings.
*)

let peakPosStdDev =
    ChargeState.peakPosStdDevBy putativeCharges

printfn "peak position standard deviation: %.4f" peakPosStdDev

(**
```text
peak position standard deviation: 0.0074
```

`initMzDevOfRndSpec` then builds a generator of random clusters: each simulated
cluster has as many peaks as the observed one, and every gap is computed as one over
a jittered charge (the implementation draws a charge from 1 up to one below the
configured maximum and adds normal jitter with the standard deviation just estimated
to the charge before taking the reciprocal). `empiricalPValueOfSim` generates `NrOfRndSpectra` such clusters
(memoized per peak count and charge, so repeated calls are cheap) and returns the
fraction of them whose deviation from the assigned charge's spacing is at most as
large as the observed one. That fraction is the empirical p-value. The design
described in the module author's thesis tests each candidate that survives the subset
filter this way, against a distribution of 10000 random clusters, and rejects the
null hypothesis at a significance criterion of 0.05.
*)

let rnd = System.Random()

let mzDevOfRandomCluster =
    ChargeState.initMzDevOfRndSpec rnd chargeParams peakPosStdDev

let bestHit = assignedCharges.Head

let pValue =
    ChargeState.empiricalPValueOfSim
        mzDevOfRandomCluster
        (bestHit.SubSetLength, float bestHit.PrecCharge)
        bestHit.MZChargeDev

printfn "empirical p-value for charge %i+: %.2f" bestHit.PrecCharge pValue

(**
<!-- stochastic: value varies between runs -->
```text
empirical p-value for charge 2+: 0.25
```

Part of the simulation draws its jitter from a random source the caller cannot seed,
so this number varies slightly between runs. It lies far above the 0.05 criterion, so
the null hypothesis stands and a pipeline applying the thesis criterion would drop
this candidate on the position metric alone. That is the expected verdict for a
cluster this short: every random gap is itself a legal charge spacing, and with only
two gaps to check, roughly a quarter of the random clusters draw the matching charge
twice and fit at least as well as the real measurement. Production pipelines
therefore also lean on the score ranking and assume typical charges such as 2+ and 3+
as a fallback when no candidate is convincing. The computed value can be stored in
the `PositionMetricPValue` field for exactly this kind of filtering.

## Where the charge goes next

Charge and mass together are what the next stage consumes, whether they come from a
ranked candidate like this one or from the 2+ and 3+ fallback: the mass selects the
candidate peptides to score against the fragment spectrum, the subject of
[peptide search databases]({{root}}02_02_search_databases.html).
*)
