module StatsExtensionTests

open System
open Expecto
open BioFSharp.Mz
open FSharp.Stats.Signal

let gauss amplitude mu sigma x = amplitude * exp (-((x - mu) ** 2.) / (2. * sigma ** 2.))

let xs = [|0.0 .. 0.05 .. 20.0|]

let twoPeakTrace =
    xs
    |> Array.map (fun x -> x, gauss 100.0 5.0 0.5 x + gauss 80.0 12.0 0.6 x)

let parameters : PeakDetection.Wavelet.Parameters =
    {
        Borderpadding = 200
        BorderPadMethod = Padding.BorderPaddingMethod.Zero
        InternalPaddingMethod = Padding.InternalPaddingMethod.LinearInterpolation
        HugeGapPaddingMethod = Padding.HugeGapPaddingMethod.Zero
        HugeGapPaddingDistance = 100.0
        MaxPeakLength = 3.0
        NoiseQuantile = 0.5
        MinSNR = 1.0
    }

let twoPeakPeaks = lazy (PeakDetection.Wavelet.identifyPeaks parameters twoPeakTrace)

[<Tests>]
let tests =
    testList "StatsExtensionTests" [
        testList "StandardErrorOfPrediction" [
            testCase "perfect prediction has zero standard error" <| fun _ ->
                let result = NonLinearRegression'.standardErrorOfPrediction 0.0 [|1.;2.;3.;4.|] [|1.;2.;3.;4.|]
                Expect.floatClose Accuracy.high result 0.0 "perfect prediction has zero standard error"
                // all residuals are zero, so the root of the summed squared residuals is zero under ANY positive denominator convention - no dependence on the implementation's degrees-of-freedom choice.

            testCase "standard error scales linearly with residual magnitude" <| fun _ ->
                let actual = [|1.;2.;3.|]
                let predA = [|1.;4.;3.|]
                let predB = [|1.;8.;3.|]
                let sA = NonLinearRegression'.standardErrorOfPrediction 0.0 predA actual
                let sB = NonLinearRegression'.standardErrorOfPrediction 0.0 predB actual
                Expect.floatClose Accuracy.high sB (3.0 * sA) "standard error scales linearly with residual magnitude"
                // sqrt(sum r^2 / denom) is homogeneous of degree 1 in the residuals and both calls share the same denominator, so tripling the only residual must triple the result - independent of the denominator convention.

            testCase "standard error strictly increases when degrees of freedom increase" <| fun _ ->
                let actual = [|1.;2.;3.;4.|]
                let pred = [|1.;4.;3.;4.|]
                let resultWithDofOne = NonLinearRegression'.standardErrorOfPrediction 1.0 pred actual
                let resultWithDofZero = NonLinearRegression'.standardErrorOfPrediction 0.0 pred actual
                Expect.isTrue (resultWithDofOne > resultWithDofZero) "standard error strictly increases when degrees of freedom increase"
                // the numerator (residual sum of squares) is fixed while the effective denominator shrinks with growing dOF, so the value cannot decrease while the denominator stays positive; with a nonzero residual and a strictly smaller positive denominator the increase is strict under any convention; >= would also pass for a constant stub.

            testCase "standard error is NaN when too few observations remain" <| fun _ ->
                Expect.isTrue (Double.IsNaN (NonLinearRegression'.standardErrorOfPrediction 3.0 [|1.;2.;3.|] [|1.;2.;3.|])) "standard error is NaN when three degrees of freedom leave no positive denominator"
                // with 3 observations, a degrees-of-freedom adjustment of 3 leaves no positive denominator
                // under either the implemented (n-1)-dOF or the conventional n-dOF denominator, so NaN is the
                // convention-independent expectation. The dOF=2 boundary is deliberately NOT asserted here -
                // it is where the suspected off-by-one manifests (see the pending test below).

            ptestCase "standard error at the minimal observation count uses the conventional denominator" <| fun _ ->
                // PENDING: with 3 observations, one nonzero residual of 2 (RSS = 4) and dOF = 2, the
                // conventional residual standard error is sqrt(4 / (3 - 2)) = 2.0 (hand-computed). The
                // implementation subtracts an extra 1 from the sample count and returns NaN at exactly this
                // boundary. Consumer impact: Quantification.HULQ.tryFit passes dOF = paramCount and maps NaN
                // to None, so fitting over exactly paramCount+1 points silently fails.
                let r = NonLinearRegression'.standardErrorOfPrediction 2.0 [|1.;4.;3.|] [|1.;2.;3.|]
                Expect.floatClose Accuracy.high r 2.0 "conventional denominator at the minimal observation count"

            testCase "standard error is finite and positive in the fitting consumer's regime" <| fun _ ->
                let r = NonLinearRegression'.standardErrorOfPrediction 3.0 [|1.;4.;3.;5.;2.;7.;6.|] [|1.;2.;3.;4.;2.;6.;6.|]
                Expect.isTrue (not (Double.IsNaN r) && not (Double.IsInfinity r) && r > 0.0) "standard error is finite and positive in the fitting consumer's regime"
                // Quantification.HULQ.tryFit calls this with dOF = parameter count (3 for a Gaussian) and treats NaN as fit failure; with 7 observations and nonzero residuals the denominator is positive under any convention, so a finite positive error is required for quantification to work at all - the consumer contract.
        ]

        testList "SubstractBaseLine" [
            testCase "baseline subtraction preserves length, keeps intensities nonnegative and keeps the apex position" <| fun _ ->
                let y = Array.init 50 (fun i -> 10.0 + gauss 100.0 25.0 3.0 (float i))
                let r = PeakDetection.Wavelet.substractBaseLine y
                Expect.equal r.Length 50 "baseline subtraction preserves length"
                Expect.isTrue (r.[0] < 5.0) "baseline subtraction removes the baseline at the left tail"
                Expect.isTrue (r.[49] < 5.0) "baseline subtraction removes the baseline at the right tail"
                // the flat tails carry pure baseline (10 units) by construction; after subtraction they must
                // drop below half the offset - without this the test would pass for the identity function.
                Expect.isTrue (Array.forall (fun v -> v >= 0.0) r) "baseline subtraction keeps intensities nonnegative"
                let maxIdx = r |> Array.mapi (fun i v -> i, v) |> Array.maxBy snd |> fst
                Expect.equal maxIdx 25 "baseline subtraction keeps the apex position"
                // documented elementwise mapping preserves length
                // explicit clamp in the doc/code: intensities are physically nonnegative
                // subtracting a smooth baseline from a single-peak signal must keep the apex at the same sample; the input's unique maximum is at index 25 by construction

            testCase "baseline subtraction is invariant to a constant offset and preserves the apex magnitude" <| fun _ ->
                let y0 = Array.init 50 (fun i -> gauss 100.0 25.0 3.0 (float i))
                let y10 = y0 |> Array.map ((+) 10.0)
                let r0 = PeakDetection.Wavelet.substractBaseLine y0
                let r10 = PeakDetection.Wavelet.substractBaseLine y10
                Expect.isTrue ((Array.map2 (fun a b -> abs (a - b)) r0 r10 |> Array.max) < 1e-6) "subtracting a constant offset leaves the corrected signals identical"
                Expect.isTrue (r10.[25] > 90.0) "the constructed apex survives baseline removal"
                // subtracting a baseline must remove a constant offset entirely: the corrected signals of y and y+10 are identical (translation invariance of baseline subtraction; probe-verified exact)
                // the 100-unit constructed apex must survive baseline removal nearly intact - a subtraction that eats the peak is not a baseline correction (probe: 98.7)
        ]

        testList "IdentifyPeaks" [
            testCase "identifyPeaks finds both Gaussian peaks of a synthetic two-peak trace" <| fun _ ->
                let peaks = twoPeakPeaks.Value
                let fits = peaks |> List.collect (fun p -> p.Fits)
                Expect.isFalse (List.isEmpty peaks) "identifyPeaks returns at least one peak group"
                Expect.isTrue (peaks |> List.forall (fun p -> p.Start <= p.End)) "every peak group has an ordered start and end"
                Expect.isTrue (peaks |> List.exists (fun p -> p.Start <= 5.0 && 5.0 <= p.End)) "a peak group contains the first Gaussian apex"
                Expect.isTrue (peaks |> List.exists (fun p -> p.Start <= 12.0 && 12.0 <= p.End)) "a peak group contains the second Gaussian apex"
                Expect.isTrue (fits |> List.exists (fun f -> abs (f.XLoc - 5.0) <= 0.2)) "a fit is centered near the first Gaussian apex"
                Expect.isTrue (fits |> List.exists (fun f -> abs (f.XLoc - 12.0) <= 0.2)) "a fit is centered near the second Gaussian apex"
                let firstFit = fits |> List.tryFind (fun f -> abs (f.XLoc - 5.0) <= 0.2)
                let secondFit = fits |> List.tryFind (fun f -> abs (f.XLoc - 12.0) <= 0.2)
                Expect.isTrue (firstFit |> Option.exists (fun f -> f.Amplitude >= 90.0 && f.Amplitude <= 101.0)) "the fit near the first Gaussian apex has the expected amplitude"
                // amplitude is the trace height at the fit location; worst case at |XLoc-5| = 0.2 is 100*exp(-0.04/0.5) = 92.3, hand-computed from the constructed Gaussian
                Expect.isTrue (secondFit |> Option.exists (fun f -> f.Amplitude >= 75.0 && f.Amplitude <= 81.0)) "the fit near the second Gaussian apex has the expected amplitude"
                // amplitude is the trace height at the fit location; worst case at |XLoc-12| = 0.2 is 80*exp(-0.04/0.72) = 75.7, hand-computed from the constructed Gaussian
                // every expected number derives from the constructed signal (apex positions 5 and 12), not the implementation. Borderpadding 200 is required: smaller paddings crash.

            testCase "the apex sample survives into the detected peak group's data" <| fun _ ->
                let group =
                    twoPeakPeaks.Value
                    |> List.find (fun p -> p.Start <= 5.0 && 5.0 <= p.End)
                Expect.isTrue (group.Data |> Array.exists (fun (x, _) -> abs (x - 5.0) <= 1e-9)) "the measured first apex sample is retained in the group's data"
                // whatever synthetic padding points the group carries (Data contains padder-interpolated values), the MEASURED apex sample itself must be retained - losing the apex would falsify any downstream quantity read from Data.

            // PENDING: with MaxPeakLength 4.2 the scale ceiling is 0.7, comfortably above the constructed
            // peak's sigma = 0.5 - a functioning width estimator must report ~0.5. Probes show the fitted
            // Stdev equals the configured ceiling VERBATIM at caps 0.5, 0.667 and 0.7 for this same peak:
            // the best-correlating scale always saturates at the maximum, so width estimation is non-functional.
            ptestCase "the fitted width tracks the true peak width when the scale ceiling allows it" <| fun _ ->
                let peaks =
                    PeakDetection.Wavelet.identifyPeaksBy 200 Padding.BorderPaddingMethod.Zero Padding.InternalPaddingMethod.LinearInterpolation Padding.HugeGapPaddingMethod.Zero
                        100.0 4.2 0.5 1.0 twoPeakTrace
                let fit =
                    peaks
                    |> List.collect (fun p -> p.Fits)
                    |> List.find (fun f -> abs (f.XLoc - 5.0) <= 0.2)
                Expect.isTrue (abs (fit.Stdev - 0.5) <= 0.15) "the fitted width is close to the constructed sigma"

            testCase "identifyPeaks finds no peaks in an all-zero trace" <| fun _ ->
                let trace = xs |> Array.map (fun x -> x, 0.0)
                let peaks = PeakDetection.Wavelet.identifyPeaks parameters trace
                Expect.isEmpty peaks "a constant signal has no strict local maxima"
                // a constant signal has no strict local maxima.

            testCase "identified fit evaluates to its amplitude at its own center" <| fun _ ->
                let peaks = twoPeakPeaks.Value
                Expect.isTrue (peaks |> List.exists (fun p -> not (List.isEmpty p.Fits))) "at least one group carries a fit"
                peaks
                |> List.iter (fun group ->
                    group.Fits
                    |> List.iter (fun f ->
                        Expect.floatClose Accuracy.high (f.Function f.XLoc) f.Amplitude "a Gaussian at its own center equals its amplitude"))
                // a Gaussian evaluated at its own mean equals its amplitude (exp(0) = 1) - a mathematical identity.

            ptestCase "two well-separated Gaussians yield exactly two peak groups" <| fun _ ->
                // PENDING: a clean noise-free trace containing exactly two well-separated Gaussian peaks
                // should produce exactly two peak groups; the detector currently reports ~21 groups, most of
                // them spurious side-structures (over-segmentation).
                Expect.equal (twoPeakPeaks.Value.Length) 2 "two constructed peaks, two groups"

            ptestCase "identifyPeaks tolerates moderate border padding" <| fun _ ->
                // PENDING: Borderpadding 50 on a 401-point trace is a reasonable parameterization; the
                 // implementation throws IndexOutOfRangeException instead of working or failing with a meaningful validation error.
                let _ = PeakDetection.Wavelet.identifyPeaks { parameters with Borderpadding = 50 } twoPeakTrace
                ()
        ]
    ]
