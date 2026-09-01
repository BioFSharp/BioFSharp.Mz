(**
---
title: Peak quantification
category: Quantification
categoryindex: 3
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

# Peak quantification

The previous pages end with an identified peptide: the fragment spectrum recorded at
20.93 min matches ANLGMEVMHER, and its precursor sits at m/z 643.80 carrying two
charges (see [charge state determination]({{root}}01_03_charge_state_determination.html)).
Identification answers which peptide was measured. Quantification asks how much of it
was there, and the answer comes from a different slice through the data.

A peptide elutes from the LC column over a stretch of retention time, together with
many other peptides that happen to leave the column at the same moment, and every
MS1 survey scan taken during that stretch records its isotope cluster once more.
All signals of one peptide ion species therefore form a three-dimensional feature
along m/z, retention time and intensity. Extracting the intensity at one m/z across
consecutive MS1 scans flattens that feature into the extracted ion chromatogram
(XIC): a trace that climbs while the peptide elutes and falls back to the noise
level once it has passed. The area
under that chromatographic peak is proportional to the amount of peptide, and area
integration is a more stable and reproducible quantity than the bare apex height.

For the ANLGMEVMHER precursor one would extract the intensities around m/z 643.80
from every MS1 scan of the run. The example data in this repository holds a single
MS1 scan, a single time point on such a trace. This page therefore builds a
realistic XIC trace synthetically, which has an advantage for a walkthrough: the
generating parameters are known exactly, so every result can be checked against
ground truth.

## Building an XIC trace

XICs approximately follow a Gaussian shape, with one recurring deviation: peak
tailing, a positive skew commonly caused by column overloading or unwanted
retention mechanisms that keep part of the analyte on the column a little longer.
The exponentially modified Gaussian (EMG) describes such tailing peaks well, so the
main synthetic peak is generated from exactly that function. A second, smaller feature elutes 0.85 min later, and
seeded uniform noise imitates the positive chemical background an instrument
records between peaks.
*)

open System
open FSharp.Stats
open BioFSharp.Mz
open BioFSharp.Mz.Quantification

// Exponentially modified Gaussian: a Gaussian of height amplitude, center meanX
// and width sigma, convolved with an exponential decay of relaxation time tau.
// The tau term produces the right-sided tail typical of chromatographic peaks.
let emg amplitude meanX sigma tau x =
    (amplitude * sigma / tau) * sqrt (Math.PI / 2.)
    * exp (0.5 * (sigma / tau) ** 2. - (x - meanX) / tau)
    * SpecialFunctions.Errorfunction.Erfc ((sigma / tau - (x - meanX) / sigma) / sqrt 2.)

let gaussian amplitude meanX sigma x =
    amplitude * exp (-((x - meanX) ** 2.) / (2. * sigma ** 2.))

// Fixed seed: reruns of this page produce identical numbers.
let rnd = Random(1337)

// One MS1 scan every 0.7 s (0.0117 min) for three minutes around the 20.9 min
// elution of the running example. The main peak reaches the tens of thousands
// of counts, a realistic MS1 XIC intensity for a well ionizing peptide, and its
// sigma of 0.055 min (3.3 s) is a typical chromatographic peak width. The
// noise floor of up to 250 counts stays strictly positive, as measured
// intensities do.
let retentionTimes = [| 19.6 .. 0.0117 .. 22.6 |]

let intensities =
    retentionTimes
    |> Array.map (fun rt ->
        emg 42000. 20.9 0.055 0.09 rt
        + gaussian 9000. 21.75 0.05 rt
        + rnd.NextDouble() * 250.)

let apexIndex =
    intensities |> Array.findIndex ((=) (Array.max intensities))

printfn "%i points from %.2f to %.2f min" retentionTimes.Length retentionTimes.[0] (Array.last retentionTimes)
printfn "highest intensity: %.0f at %.3f min" intensities.[apexIndex] retentionTimes.[apexIndex]

(**
```text
257 points from 19.60 to 22.60 min
highest intensity: 27510 at 20.946 min
```

The apex lands a little right of the EMG center at 20.9 min because the exponential
tail shifts the maximum, exactly as it does on a real column. The generating EMG has
an analytic area of amplitude times sigma times the square root of two pi, about
5790 intensity·min, a number to keep in mind for the results below.

## Detecting peaks on the trace

Quantification starts by deciding where on the trace a peak begins and ends.
FSharp.Stats brings a second-derivative peak picker for this,
`Signal.PeakDetection.SecondDerivative.getPeaks`. It smooths the trace with a
Savitzky-Golay filter, takes the negative second derivative (which turns every
peak into a pronounced maximum), estimates the noise as the mean absolute
difference between the raw and the smoothed trace, and keeps candidate peaks that
rise above a signal-to-noise cutoff times that estimate. Its arguments are the
cutoff followed by the polynomial order and the window width in points of the
Savitzky-Golay filter.
*)

let peaks =
    Signal.PeakDetection.SecondDerivative.getPeaks 5. 2 13 retentionTimes intensities

peaks
|> Array.iter (fun p ->
    printfn "apex %.3f min, intensity %.0f, from %.3f to %.3f min, %i points"
        p.Apex.XVal p.Apex.YVal p.LeftEnd.XVal p.RightEnd.XVal p.XData.Length)

(**
```text
apex 20.946 min, intensity 27510, from 20.712 to 21.495 min, 68 points
apex 21.753 min, intensity 9114, from 21.542 to 21.917 min, 33 points
```

Both synthetic features are found and nothing else is. The cutoff carries that
result: with 0.5 instead of 5 the same call reports 15 peaks, most of them noise
ripples of around 200 counts. Each detected peak is an `IdentifiedPeak` holding the
apex, the left and right peak ends, the lift-off points where curvature indicates
the peak rising out of the baseline, and the slice of the measured data between the
ends. Note how the first peak's slice reaches from 20.71 to 21.50 min, well beyond
a symmetric window around the apex, because the border search followed the tail.
This record is the input type the quantification below consumes.

## Detecting peaks with the chromatographic wavelet

This library also ships its own chromatographic peak detector,
`PeakDetection.Wavelet.identifyPeaks`, which correlates the trace with Mexican Hat
(Ricker) wavelets of increasing scale. The Mexican Hat is the negative normalized
second derivative of a Gaussian, so it matches the shape a chromatographic peak is
expected to have, and shape matching lets weak peaks stand out against noise of
similar amplitude. We hand the detector the trace on a seconds axis, where the
0.7 s sampling and peak widths of a few seconds give the wavelet a convenient
scale range. Two parameters deserve a note. `Borderpadding` must be generous
relative to the wavelet width, and 200 works well on traces of this size.
`MaxPeakLength` caps the widest wavelet scale at one sixth of its value, and 30 s
sits comfortably above six times the 3.3 s sigma of the synthetic peak.
*)

let traceInSeconds =
    Array.map2 (fun rt intensity -> rt * 60., intensity) retentionTimes intensities

let waveletParameters : PeakDetection.Wavelet.Parameters =
    { Borderpadding = 200
      BorderPadMethod = Signal.Padding.BorderPaddingMethod.Zero
      InternalPaddingMethod = Signal.Padding.InternalPaddingMethod.LinearInterpolation
      HugeGapPaddingMethod = Signal.Padding.HugeGapPaddingMethod.Zero
      HugeGapPaddingDistance = 100.
      MaxPeakLength = 30.
      NoiseQuantile = 0.5
      MinSNR = 1. }

let peakGroups =
    PeakDetection.Wavelet.identifyPeaks waveletParameters traceInSeconds

let coveredByGroup rtInMinutes =
    peakGroups
    |> List.exists (fun g -> g.Start <= rtInMinutes * 60. && rtInMinutes * 60. <= g.End)

printfn "peak groups found: %i" peakGroups.Length
printfn "a group covers the main apex: %b" (coveredByGroup peaks.[0].Apex.XVal)
printfn "a group covers the second apex: %b" (coveredByGroup peaks.[1].Apex.XVal)

(**
```text
peak groups found: 29
a group covers the main apex: true
a group covers the second apex: true
```

The detector finds both true apexes, and it detects generously: on this two-feature
trace it reports 29 peak groups, treating shoulders and noise structures next to the
real peaks as separate groups. A consumer therefore selects the group covering the
retention time of interest and disregards the rest. A group's `Data` holds the
padded trace region the wavelet worked on, and the fitted `Stdev` reflects the
wavelet scale the group responded to within the configured range. For the rest of
this page we stay with the second-derivative peaks, which carry exact slice
boundaries and feed directly into the quantification step.

## Fitting and integrating the main peak

The `HULQ` module turns an `IdentifiedPeak` into a quantity. `HULQ.getPeakBy`
selects from the detected peaks the one whose boundaries contain a target retention
time, falling back to the nearest apex when no peak contains it. The natural target
is the retention time of the identification, 20.93 min for our MS2 scan.

`HULQ.quantifyPeak` then runs the model fitting workflow. It estimates starting
parameters from the peak slice (by weighted moments, with a Caruana log-parabola
fit stepping in when the slice is truncated on one side), fits a Gaussian and, when
the moments allow, an EMG via Levenberg-Marquardt, scores both fits by their
standard error of prediction, keeps the better one, and integrates the fitted model
analytically for the area.
*)

let describeModel (q: HULQ.QuantifiedPeak) =
    match q.Model with
    | Some (HULQ.Gaussian _) -> "Gaussian"
    | Some (HULQ.EMG _) -> "EMG"
    | None -> "none, trapezoid fallback"

let printQuantification (q: HULQ.QuantifiedPeak) =
    printfn "model selected              : %s" (describeModel q)
    match q.Model, q.EstimatedParams with
    | Some (HULQ.Gaussian _), [| amp; meanX; sigma |] ->
        printfn "amplitude                   : %.1f" amp
        printfn "mean retention time         : %.3f min" meanX
        printfn "sigma                       : %.4f min" sigma
    | Some (HULQ.EMG _), [| amp; meanX; sigma; tau |] ->
        printfn "amplitude                   : %.1f" amp
        printfn "mean retention time         : %.3f min" meanX
        printfn "sigma                       : %.4f min" sigma
        printfn "tau                         : %.4f min" tau
    | _ -> ()
    printfn "standard error of prediction: %.1f" q.StandardErrorOfPrediction
    printfn "area under the fitted model : %.1f" q.Area
    printfn "measured apex intensity     : %.1f" q.MeasuredApexIntensity

let mainPeak = HULQ.getPeakBy peaks 20.93

let mainQuantification = HULQ.quantifyPeak mainPeak

printQuantification mainQuantification

(**
```text
model selected              : EMG
amplitude                   : 35762.5
mean retention time         : 20.906 min
sigma                       : 0.0629 min
tau                         : 0.0759 min
standard error of prediction: 601.2
area under the fitted model : 5635.9
measured apex intensity     : 27509.6
```

The tailed peak is recognized as an EMG. The fitted mean retention time of
20.906 min sits on the generating center of 20.9, and sigma and tau land near the
generating 0.055 and 0.09 min. The amplitude is the height parameter of the
pre-convolution Gaussian, so it neither matches the measured apex of 27510 counts
nor needs to. The `MeasuredApexIntensity` field preserves that raw apex alongside
the model.

How good is the model area? A direct numerical integration of the same slice gives
an assumption-free reference.
*)

let mainTrapezoid =
    Integration.trapezEstAreaOf mainPeak.XData mainPeak.YData

printfn "model area    : %.1f" mainQuantification.Area
printfn "trapezoid area: %.1f" mainTrapezoid

(**
```text
model area    : 5635.9
trapezoid area: 5874.7
```

Both land close to the analytic area of 5790 intensity·min the trace was generated
with. The trapezoid integral runs a little above it because it also sums the noise
floor riding on the slice, the model area a little below. The model area is the
value a pipeline stores: it is anchored to a fitted peak shape, which keeps it
comparable when slice boundaries or noise differ between runs.

## How the model is chosen

Whether an EMG is attempted at all is decided by the moment estimates.
`ParameterEstimation.estTau` derives the EMG tail parameter from the moments as
sigma times the cube root of half the skewness. For a peak with zero or negative
skew that fractional power of a non-positive number is NaN, the EMG candidate
becomes unavailable, and the Gaussian is used. That is the intended selection path:
a peak that is not right-skewed is a Gaussian case, and no tail parameter should be
invented for it. The moments of our two peaks show both sides of this rule.
*)

let secondPeak = HULQ.getPeakBy peaks 21.75

let printMoments name peak =
    match ParameterEstimation.estimateMoments peak with
    | Some m ->
        printfn "%s: mean %.3f min, sigma %.4f min, skew %.3f, tau %.4f" name m.MeanX m.Std m.Skew m.Tau
    | None -> printfn "%s: no moment estimate" name

printMoments "main peak  " mainPeak
printMoments "second peak" secondPeak

(**
```text
main peak  : mean 20.991 min, sigma 0.1068 min, skew 1.109, tau 0.0878
second peak: mean 21.748 min, sigma 0.0555 min, skew -0.212, tau NaN
```

The tailed main peak carries a clear positive skew and a defined tau, so both
models compete and the output above showed the EMG winning on standard error. The
second peak was generated symmetric, its slight negative skew is noise, and its tau
is NaN. For it the EMG leg never runs.

## Quantifying the second feature

The workflow repeats per peptide feature: select the peak at the feature's
retention time, quantify it.
*)

let secondQuantification = HULQ.quantifyPeak secondPeak

printQuantification secondQuantification

(**
```text
model selected              : Gaussian
amplitude                   : 9032.7
mean retention time         : 21.749 min
sigma                       : 0.0516 min
standard error of prediction: 121.0
area under the fitted model : 1167.3
measured apex intensity     : 9113.5
```

As the moments predicted, the symmetric feature is fitted as a Gaussian, and for a
Gaussian the amplitude is the apex height, so 9033 sits right at the generating
9000. Mean and sigma recover the generating 21.75 min and 0.05 min. Should neither
model converge on some degenerate slice, `quantifyPeak` still returns a usable
result, the trapezoid area of the slice with `Model = None`.

## From one peak to a pipeline

The error measure behind the model selection is
`NonLinearRegression'.standardErrorOfPrediction` from this library's FSharp.Stats
extension, the root of the summed squared residuals normalized by the residual
degrees of freedom, and the model with the lower value wins. Around this core a
real pipeline adds the comparative layer: for isotopic labeling such as 15N, the
light and the heavy XIC of the same peptide are extracted and their fitted areas
are ratioed within one run, while label-free workflows compare the areas of
matching features across runs.

Before peak areas are aggregated into protein quantities, the identifications they
hang on are filtered to a controlled error rate, which is the subject of
[FDR control]({{root}}04_01_fdr_control.html). And the m/z every XIC follows comes
from the precursor's charge assignment on the
[charge state determination]({{root}}01_03_charge_state_determination.html) page.
*)
