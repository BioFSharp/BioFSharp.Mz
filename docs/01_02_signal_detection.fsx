(**
---
title: Signal detection and centroiding
category: Spectrum processing
categoryindex: 1
index: 2
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

# Signal detection and centroiding

A mass spectrometer in profile mode writes down the detector signal as a
continuous-looking trace: many m/z and intensity pairs per true ion signal. Every ion
species shows up as a roughly Gaussian bump that is a few dozen points wide, with its
apex at the ion's m/z. Charge state determination, database search, spectrum scoring and
quantification all expect one peak per ion species. Centroiding reduces every bump to a
single representative peak, an m/z position plus an apex or summed intensity.

`SignalDetection` implements a wavelet-based peak picker modeled after French et al.
(doi: [10.1021/pr500886y](https://doi.org/10.1021/pr500886y)). The idea: a Mexican hat
wavelet (a Ricker wavelet, which looks like a Gaussian bump with dips on both sides) is
slid along the m/z axis at several widths and correlated with the intensity trace.
Wherever the trace contains a bump of matching width, the correlation is high. True
peaks score well across scales while single-point noise spikes do not, which is what
makes the approach less trigger-happy than simple local maxima. This page covers the
wavelet approach.

## Loading a profile MS1 spectrum

The example file holds one MS1 survey scan (a full scan over the precursor mass range,
recorded before the instrument picks ions for fragmentation) stored in mgf format. We
load it with the reader from BioFSharp, as on the
[peaks and peak arrays]({{root}}01_01_peaks_and_peak_arrays.html) page, and check
empirically that it really is profile data by looking at the number of points and the
spacing between neighboring points.
*)

open BioFSharp.FileFormats.MGF
open BioFSharp.IO
open BioFSharp.Mz

let ms1 =
    MGF.read (__SOURCE_DIRECTORY__ + "/data/ms1Example.mgf")
    |> List.head

let profileMz = ms1.Mass
let profileIntensity = ms1.Intensity

printfn "%s" (MGFEntry.tryGetTitle ms1 |> Option.defaultValue "no title")
printfn "points in scan: %i" profileMz.Length
printfn "m/z range:      %f to %f" profileMz.[0] (Array.last profileMz)

let spacings =
    profileMz |> Array.pairwise |> Array.map (fun (a, b) -> b - a)

spacings
|> Array.skip 2
|> Array.truncate 3
|> Array.iter (printfn "spacing between adjacent points: %f")

printfn "median spacing: %f" ((Array.sort spacings).[spacings.Length / 2])
printfn "largest gap:    %f" (Array.max spacings)

(**
```text
MS scan at 20.91351667 min
points in scan: 78170
m/z range:      400.001711 to 1249.609952
spacing between adjacent points: 0.005613
spacing between adjacent points: 0.005614
spacing between adjacent points: 0.002807
median spacing: 0.003504
largest gap:    1.481318
```

78170 points with a median spacing of about 0.0035 Da is profile data. A centroided MS1
scan of the same range would have a few hundred to a few thousand peaks, one per
detected ion signal. The largest gap tells the second half of the story: the
instrument leaves out stretches where it saw nothing, so the sampling is dense inside
peak regions and sparse between them.

## Padding sparse regions

The wavelet correlation assumes reasonably dense, gap-free sampling, because the width
of the wavelet is derived from the local point spacing. A 1.5 Da hole in the data would
inflate that estimate and distort the correlation around the hole. `Padding.paddDataBy`
therefore fills sparse regions with artificial points at a fixed base intensity before
the peak picker runs.

Its parameters, built with `Padding.createPaddingParameters`: `paddingYValue` is the
intensity written into every inserted point, typically 0. or the minimum intensity of
the data. `maximumPaddingPoints` caps how many points may be inserted into one gap,
and `None` fills gaps completely. `mzTolerance` sets the m/z distance at which the
inserted points are placed, here 0.05 Da, finer than any real peak width.
`windowSize` and `spacingPerc` set the points per window and the percentile of the
spacings inside it used to estimate the typical local point spacing. A spacing is
padded as a gap when it stands out against this local estimate and exceeds
`mzTolerance`.
*)

let paddingParams =
    SignalDetection.Padding.createPaddingParameters 0. None 0.05 150 95.

let paddedMz, paddedIntensity =
    SignalDetection.Padding.paddDataBy paddingParams profileMz profileIntensity

printfn "points before padding: %i" profileMz.Length
printfn "points after padding:  %i" paddedMz.Length

(**
```text
points before padding: 78170
points after padding:  83672
```

The added points all carry intensity 0. and only exist so the wavelet sees a continuous
baseline between the real peak regions.

## Centroiding with the wavelet peak picker

`Wavelet.toCentroidWithRicker2D` takes a `WaveletParameters` record and the padded
arrays and returns the centroided m/z and intensity arrays. The parameters, built with
`Wavelet.createWaveletParameters`:

`numberOfScales` is the number of wavelet widths tried at every position. The widths
span a fixed range from one to seven times the local point spacing, so larger values
sample that range more finely at more cost. `yThreshold` is an absolute intensity
floor below which a bump's apex is dropped, and with 0. the filtering is left to the
signal-to-noise criterion. `mzTolerance` is the grouping tolerance in m/z: candidate
peaks closer together than this compete, and only the best-scoring one survives.
`sNRS_Percentile` sets which percentile of the local correlation scores is taken as
the noise estimate, and `minSNR` is the minimum ratio between a candidate's
correlation score and that estimate. `refineMZ` reports the m/z as the
intensity-weighted mean over the bump (`true`) or the position of the apex point
(`false`), and `sumIntensities` reports the intensity as the sum over the bump
(`true`) or the apex intensity (`false`).

The values below are what a production proteomics pipeline uses for MS1 scans: 3
scales, mzTolerance 0.1, the 95th percentile for the noise estimate and a minimum SNR
of 1. For the denser features of MS2 spectra the same pipeline raises the number of
scales to 10. We report apex positions and apex intensities first.
*)

let waveletParams =
    SignalDetection.Wavelet.createWaveletParameters 3 0. 0.1 95. 1. false false

let centroidMz, centroidIntensity =
    SignalDetection.Wavelet.toCentroidWithRicker2D waveletParams paddedMz paddedIntensity

printfn "profile points: %i" profileMz.Length
printfn "centroids:      %i" centroidMz.Length

(**
```text
profile points: 78170
centroids:      1454
```

78170 profile points become 1454 centroids, one per detected ion signal, a reduction of
roughly a factor of fifty on this scan.

## Watching a bump collapse to a centroid

Numbers over the whole scan hide what happens locally, so we zoom into a 1.5 Da window
around the most intense signal of the scan. First the profile points near its apex.
*)

let inWindow lo hi (mz, _) = mz >= lo && mz <= hi

let profileWindow =
    Array.zip profileMz profileIntensity
    |> Array.filter (inWindow 434.9 436.4)

printfn "profile points between m/z 434.9 and 436.4: %i" profileWindow.Length

profileWindow
|> Array.filter (inWindow 435.18 435.21)
|> Array.iter (fun (mz, i) -> printfn "m/z %f  intensity %8.1f" mz i)

(**
```text
profile points between m/z 434.9 and 436.4: 398
m/z 435.180270  intensity   1791.0
m/z 435.183197  intensity   4332.0
m/z 435.186125  intensity   7143.0
m/z 435.189052  intensity  10025.0
m/z 435.191980  intensity   9749.0
m/z 435.194907  intensity   7993.0
m/z 435.197835  intensity   6344.0
m/z 435.200762  intensity   3923.0
m/z 435.203690  intensity   2902.0
m/z 435.206617  intensity   1836.0
m/z 435.209545  intensity   1337.0
```

The intensity rises point by point to an apex at m/z 435.189 and falls off again, the
Gaussian bump from the introduction in the flesh. Now the centroids the wavelet produced
for the same window.
*)

Array.zip centroidMz centroidIntensity
|> Array.filter (inWindow 434.9 436.4)
|> Array.iter (fun (mz, i) -> printfn "m/z %f  intensity %8.1f" mz i)

(**
```text
m/z 435.189052  intensity  10025.0
m/z 435.525777  intensity   5519.0
m/z 435.856773  intensity   2806.0
m/z 436.190825  intensity   1053.0
```

398 profile points collapse to four centroids. The first sits exactly on the apex we
just looked at. The four peaks are spaced about 0.334 Da apart with intensities falling
off from left to right, which is the signature of an isotope cluster. Each following
peak is the same ion species with one more heavy isotope, and since the mass
difference between carbon-13 and carbon-12 is about 1.0034 Da, a spacing of
1.0034 / 3 means the ion carries three charges.

Rerunning with `refineMZ` and `sumIntensities` turned on changes what each centroid
reports. The m/z becomes the intensity-weighted mean over the bump and the intensity
becomes the sum over all its profile points, which is the better quantity when peak
areas are compared later.
*)

let summingParams =
    SignalDetection.Wavelet.createWaveletParameters 3 0. 0.1 95. 1. true true

let refinedMz, refinedIntensity =
    SignalDetection.Wavelet.toCentroidWithRicker2D summingParams paddedMz paddedIntensity

Array.zip refinedMz refinedIntensity
|> Array.filter (inWindow 434.9 436.4)
|> Array.iter (fun (mz, i) -> printfn "m/z %f  intensity %8.1f" mz i)

(**
```text
m/z 435.192899  intensity  59306.0
m/z 435.526224  intensity  31209.0
m/z 435.859807  intensity  17631.0
m/z 436.193927  intensity   6746.0
```

Same four signals, now with summed intensities and slightly shifted m/z positions
because the weighted mean also feels the asymmetric tails of each bump.

## Where the centroids go next

The centroided arrays can be zipped back into a `PeakArray` (see
[peaks and peak arrays]({{root}}01_01_peaks_and_peak_arrays.html)). The isotope spacing
we read off by eye above is what
[charge state determination]({{root}}01_03_charge_state_determination.html) exploits
systematically to assign each precursor its charge, needed to turn an m/z value into a
peptide mass for the database search.
*)
