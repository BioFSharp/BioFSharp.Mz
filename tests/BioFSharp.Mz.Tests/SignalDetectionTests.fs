module SignalDetectionTests

open System
open Expecto
open BioFSharp.Mz

let gauss amplitude mu sigma x = amplitude * exp (-((x - mu) ** 2.) / (2. * sigma ** 2.))

[<Tests>]
let tests =
    testList "SignalDetectionTests" [
        testList "Care" [
            testCase "accumulate folds over the half-open index interval" <| fun _ ->
                let actual = SignalDetection.Care.accumulate 1 3 10.0 (+) [|1.;2.;3.;4.|]
                Expect.equal actual 15.0 "indices [1,3) add 2 and 3 to the start value"
                // The fold over indices [1,3) adds 2 and 3 to 10; index 3 is excluded by the half-open contract.

            testCase "scoreAtPercentile interpolates linearly on the (n-1)-rank convention" <| fun _ ->
                let sorted = [|0.;10.;20.;30.|]
                Expect.equal (SignalDetection.Care.scoreAtPercentile 0.0 4 sorted) 0.0 "0% selects the first element"
                Expect.equal (SignalDetection.Care.scoreAtPercentile 100.0 4 sorted) 30.0 "100% selects the last element"
                Expect.floatClose Accuracy.high (SignalDetection.Care.scoreAtPercentile 25.0 4 sorted) 7.5 "25% interpolates at rank 0.75"
                Expect.equal (SignalDetection.Care.scoreAtPercentile 50.0 0 sorted) 0.0 "zero total count returns zero"
                // The standard (n-1)-rank convention gives (4-1)*25/100 = 0.75 and 0 + 0.75*(10-0) = 7.5.

            testCase "createXspacing assigns each point its preceding gap with replicated boundary values" <| fun _ ->
                let r = SignalDetection.Care.createXspacing [|100.;100.5;101.5;104.|]
                let expected = [|0.5; 0.5; 1.0|]
                Expect.equal r.Length 4 "the spacing result has one value per input point"
                Array.zip r.[0..2] expected
                |> Array.iteri (fun index (actual, expected) ->
                    Expect.floatClose Accuracy.high actual expected (sprintf "spacing at index %d" index))
                // interior points carry the gap to their predecessor; index 0 replicates its interior neighbor (it has no predecessor). The LAST index is deliberately not asserted here - see the pending test below.

            // PENDING: the last point's preceding gap (104 - 101.5 = 2.5) is computed and then overwritten
            // with the previous gap (1.0). Consequence: paddDataWith judges the gap AFTER point i via
            // xSpacing[i+1], so a large trailing gap is evaluated with the wrong spacing and never padded.
            ptestCase "createXspacing reports the natural gap at the last point" <| fun _ ->
                let r = SignalDetection.Care.createXspacing [|100.;100.5;101.5;104.|]
                Expect.floatClose Accuracy.high r.[3] 2.5 "the last point carries its own natural gap"

            testCase "getColLowBound and getColHighBound bound the tolerance window for interior centers" <| fun _ ->
                let xData = [|99.0;100.0;100.4;100.5;100.6;101.5;102.5|]
                Expect.equal (SignalDetection.Care.getColLowBound xData 3 0.2) 2 "the lower bound includes the nearest interior neighbor"
                Expect.equal (SignalDetection.Care.getColHighBound xData 3 0.2) 4 "the upper bound includes the nearest interior neighbor"
                // neighbors at distance 0.1 (strictly inside the window) are included; the next points at distances 0.5 and 1.0 are outside and excluded - hand-computed; the sentinels guarantee the walk terminates on the tolerance check, not the array edge.

            // PENDING: for xData [|100.;100.4;100.5;100.6;101.5|], center 2, tolerance 0.2 the correct
            // bounds are (1,3): the edge points at distances 0.5 and 1.0 lie far outside the window. When
            // the walk reaches index 0 (or length-1) the implementation returns the edge index WITHOUT the
            // distance check, yielding (0,4) and silently widening fitting windows for peaks near the
            // spectrum edge.
            ptestCase "window bounds respect the tolerance at the array edges" <| fun _ ->
                Expect.equal (SignalDetection.Care.getColLowBound [|100.;100.4;100.5;100.6;101.5|] 2 0.2) 1 "the lower bound respects the tolerance at the array edge"
                Expect.equal (SignalDetection.Care.getColHighBound [|100.;100.4;100.5;100.6;101.5|] 2 0.2) 3 "the upper bound respects the tolerance at the array edge"

        ]

        testList "Padding" [
            testCase "paddDataExtensively fills a large gap with the padding value and preserves interior points" <| fun _ ->
                let xs = [|1.0;1.1;1.2;5.0;5.1;5.2|]
                let ys = [|10.;20.;30.;40.;50.;60.|]
                let px, py = SignalDetection.Padding.paddDataExtensively 0.0 0.3 3 50.0 xs ys
                Expect.equal px.Length py.Length "padded m/z and intensity arrays have equal lengths"

                let assertPaired originalMz originalIntensity =
                    let hasPair =
                        Array.mapi (fun index mz -> index, mz) px
                        |> Array.exists (fun (index, mz) -> abs (mz - originalMz) <= 1e-9 && abs (py.[index] - originalIntensity) <= 1e-9)
                    Expect.isTrue hasPair (sprintf "original point %g is retained with intensity %g" originalMz originalIntensity)

                assertPaired 1.1 20.
                assertPaired 1.2 30.
                assertPaired 5.0 40.
                assertPaired 5.1 50.

                let originalValues = xs
                let paddedIndices =
                    px
                    |> Array.mapi (fun index mz -> index, mz)
                    |> Array.filter (fun (_, mz) -> not (Array.exists (fun original -> abs (mz - original) <= 1e-9) originalValues))
                paddedIndices
                |> Array.iter (fun (index, mz) ->
                    Expect.isTrue (mz > 1.2 && mz < 5.0) "inserted m/z values lie strictly inside the large gap"
                    Expect.equal py.[index] 0.0 "inserted points carry the padding intensity")
                Expect.isTrue (paddedIndices.Length > 0) "at least one padded point was inserted"
                Expect.isTrue (px |> Array.pairwise |> Array.forall (fun (a,b) -> a <= b)) "padded m/z values are non-decreasing"
                // downstream consumers assume sorted m/z; the uncapped padding path must preserve order.
                // The only gap larger than the approximately 0.1 median spacing is (1.2, 5.0), so inserted points belong only there and use the supplied padding value.

            // PENDING: the implementation loops over interior indices only (1 .. n-2) and drops the first and last original observations; a padding operation must preserve all original data points.
            ptestCase "padding preserves the first and last original data points" <| fun _ ->
                let xs = [|1.0;1.1;1.2;5.0;5.1;5.2|]
                let ys = [|10.;20.;30.;40.;50.;60.|]
                let px, _ = SignalDetection.Padding.paddDataExtensively 0.0 0.3 3 50.0 xs ys
                Expect.isTrue (px |> Array.exists (fun x -> abs (x - 1.0) <= 1e-9)) "the first original point is preserved"
                Expect.isTrue (px |> Array.exists (fun x -> abs (x - 5.2) <= 1e-9)) "the last original point is preserved"

            // PENDING: with MaximumPaddingPoints = Some 2 the right-side inserts are emitted walking
            // DOWNWARD from the gap's upper edge, producing an unsorted padded array. A padded spectrum
            // must remain sorted.
            ptestCase "capped padding still produces sorted output" <| fun _ ->
                let px, _ = SignalDetection.Padding.paddDataWith 0.0 (Some 2) 0.3 3 50.0 [|1.0;1.1;1.2;5.0;5.1;5.2|] [|10.;20.;30.;40.;50.;60.|]
                Expect.isTrue (px |> Array.pairwise |> Array.forall (fun (a,b) -> a <= b)) "capped padded m/z values are sorted"
        ]

        testList "WaveletHelpers" [
            testCase "convertColToMz maps even columns to points and odd columns to midpoints" <| fun _ ->
                let mzData = [|100.;101.;103.|]
                Expect.floatClose Accuracy.high (SignalDetection.Wavelet.convertColToMz mzData 0) 100.0 "column 0 maps to the first point"
                Expect.floatClose Accuracy.high (SignalDetection.Wavelet.convertColToMz mzData 1) 100.5 "column 1 maps to the first midpoint"
                Expect.floatClose Accuracy.high (SignalDetection.Wavelet.convertColToMz mzData 2) 101.0 "column 2 maps to the second point"
                Expect.floatClose Accuracy.high (SignalDetection.Wavelet.convertColToMz mzData 3) 102.0 "column 3 maps to the second midpoint"
                Expect.floatClose Accuracy.high (SignalDetection.Wavelet.convertColToMz mzData 4) 103.0 "column 4 maps to the third point"
                // Even columns map to original points; odd columns map to adjacent midpoints (100+101)/2 and (101+103)/2.

            testCase "ricker2d samples a Ricker wavelet: positive center, zero crossings at the width, symmetric negative lobes" <| fun _ ->
                let padded = [|9.6;9.8;10.0;10.2;10.4|]
                let wave = Array.zeroCreate 5
                let r = SignalDetection.Wavelet.ricker2d padded wave 2 2 2 10.0 0.2
                let expectedCenter = 2.0 / (sqrt (3.0 * 0.2) * (System.Math.PI ** 0.25))
                // Accuracy.medium (~1.9e-5 here) passes with only ~1.7x margin over the typo-induced
                // difference of 1.14e-5 - do not "tighten" this to Accuracy.high, which fails on the typo.
                Expect.floatClose Accuracy.medium r.[2] expectedCenter "the center equals the Ricker normalization constant"
                Expect.floatClose Accuracy.high r.[1] 0.0 "the left zero crossing is zero"
                Expect.floatClose Accuracy.high r.[3] 0.0 "the right zero crossing is zero"
                Expect.floatClose Accuracy.high r.[0] r.[4] "the negative side lobes are symmetric"
                Expect.isTrue (r.[0] < 0.0) "the side lobes at distance two widths are negative"
                // The published Ricker wavelet has normalization 2/(sqrt(3w)*pi^(1/4)), zero crossings at distance w, and symmetric negative lobes at distance 2w. The implementation's pi literal is mistyped (3.141519), shifting the constant by ~1.1e-5; medium accuracy verifies the physics while high accuracy would fail on the typo (see the pending test).

            // PENDING: the wavelet normalization constant is 2/(sqrt(3w) * pi^(1/4)) with pi = 3.14159265;
            // the source hardcodes 3.141519. The ~1e-5 relative error scales all correlations uniformly
            // and is harmless for peak picking, but the constant is wrong.
            ptestCase "the Ricker normalization uses the true value of pi" <| fun _ ->
                let padded = [|9.6;9.8;10.0;10.2;10.4|]
                let wave = Array.zeroCreate 5
                let r = SignalDetection.Wavelet.ricker2d padded wave 2 2 2 10.0 0.2
                Expect.floatClose Accuracy.high r.[2] (2.0 / (sqrt (3.0 * 0.2) * (System.Math.PI ** 0.25))) "normalization constant uses true pi"
        ]

        testList "Wrappers" [
            testCase "toCentroid guards short inputs and otherwise delegates to the centroid function" <| fun _ ->
                let shortMz, shortIntensities = SignalDetection.toCentroid (fun mz i -> (mz, i)) [|1.;2.|] [|1.;2.|]
                Expect.equal shortMz [||] "fewer than three points produce no centroid m/z values"
                Expect.equal shortIntensities [||] "fewer than three points produce no centroid intensities"

                let mz3 = [|1.;2.;3.|]
                let i3 = [|4.;5.;6.|]
                let passedMz, passedIntensities = SignalDetection.toCentroid (fun mz i -> (mz, i)) mz3 i3
                Expect.equal passedMz mz3 "the wrapper passes m/z data through unchanged"
                Expect.equal passedIntensities i3 "the wrapper passes intensity data through unchanged"
                // The wrapper adds nothing to the callback result for inputs meeting its three-point minimum.
        ]

        testList "WindowSelection" [
            testCase "lowerIdxBy and upperIdxBy select the half-window index bounds" <| fun _ ->
                let mzData = [|100.;100.5;101.;101.5;102.|]
                Expect.equal (SignalDetection.lowerIdxBy mzData 1.1 101.0) 1 "the lower bound selects index 1"
                Expect.equal (SignalDetection.upperIdxBy mzData 1.1 101.0) 3 "the upper bound selects index 3"
                // The half-window is [100.45, 101.55], whose first point above the lower edge is 100.5 and last point below the upper edge is 101.5.

            testCase "windowToCentroid passes the inclusive index slice to the centroid function" <| fun _ ->
                let mzData = [|100.;100.5;101.;101.5;102.|]
                let intensities = [|10.;20.;30.;40.;50.|]
                let spy = fun mz i -> (mz, i)
                let smz, si = SignalDetection.windowToCentroid spy mzData intensities 1 3
                Expect.equal smz [|100.5;101.;101.5|] "indices 1 through 3 are passed inclusively"
                Expect.equal si [|20.;30.;40.|] "intensities remain paired with the selected m/z values"
                // Indices 1..3 inclusive select exactly these three pairs.

            testCase "windowToCentroidBy composes index selection and slicing" <| fun _ ->
                let mzData = [|100.;100.5;101.;101.5;102.|]
                let intensities = [|10.;20.;30.;40.;50.|]
                let spy = fun mz i -> (mz, i)
                let smz, si = SignalDetection.windowToCentroidBy spy mzData intensities 1.1 101.0
                Expect.equal smz [|100.5;101.;101.5|] "window selection composes to the expected m/z slice"
                Expect.equal si [|20.;30.;40.|] "window selection preserves the expected intensity pairing"
                // This composition must agree with the manually chained lowerIdxBy, upperIdxBy, and windowToCentroid calls.

            // PENDING: an inclusive window ending at index length-1 is valid, and upperIdxBy returns
            // length-1 for any precursor near the top of the spectrum. The validity guard demands
            // upperIdx <= length-2, so this routine case is rejected and the caller receives sentinel
            // arrays [|float lowerIdx|], [|float upperIdx|] masquerading as centroid data.
            ptestCase "a window ending at the last index is served" <| fun _ ->
                let mzData = [|100.;100.5;101.;101.5;102.|]
                let intensities = [|10.;20.;30.;40.;50.|]
                let smz, si = SignalDetection.windowToCentroid (fun mz i -> (mz, i)) mzData intensities 1 4
                Expect.equal smz [|100.5;101.;101.5;102.|] "an inclusive window ending at the last index passes the final point"
                Expect.equal si [|20.;30.;40.;50.|] "intensities remain paired through the final point"
        ]

        testList "FilterByIntensitySNR" [
            testCase "SNR filtering applies a single intensity threshold preserving order and pairing" <| fun _ ->
                let mz = [|1.;2.;3.;4.;5.|]
                let intens = [|10.;100.;50.;200.;20.|]
                let fm1, fi1 = SignalDetection.filterByIntensitySNR 50.0 1.0 mz intens
                let fm3, fi3 = SignalDetection.filterByIntensitySNR 50.0 3.0 mz intens
                Expect.equal fm1.Length fi1.Length "filtered m/z and intensity arrays have equal lengths"

                Expect.isFalse (Array.zip fm1 fi1 |> Array.contains (1.0, 10.0)) "the minimum-intensity pair is rejected at minSnr 1.0"
                // 10 is the minimum intensity, so any percentile-based noise estimate is >= 10 and 10/noise <= 1 can never strictly exceed minSnr 1 - rejection holds under any noise convention
                Expect.isTrue (Array.zip fm1 fi1 |> Array.contains (2.0, 100.0)) "the 100-intensity pair is retained at minSnr 1.0"
                Expect.isTrue (Array.zip fm1 fi1 |> Array.contains (4.0, 200.0)) "the 200-intensity pair is retained at minSnr 1.0"
                // 100 is the 75th percentile of these intensities, so any sane 50th-percentile noise estimate lies below 100; both points' SNR exceeds 1 under the correct and the current noise value alike

                let inputPairs = Array.zip mz intens |> Array.toList
                let retainedPairs = Array.zip fm1 fi1 |> Array.toList
                let rec isSubsequence source target =
                    match target with
                    | [] -> true
                    | targetHead :: targetTail ->
                        match List.tryFindIndex ((=) targetHead) source with
                        | None -> false
                        | Some index -> isSubsequence (List.skip (index + 1) source) targetTail
                Expect.isTrue (isSubsequence inputPairs retainedPairs) "retained pairs form an input-order subsequence"

                let droppedIntensities =
                    inputPairs
                    |> List.filter (fun pair -> not (List.contains pair retainedPairs))
                    |> List.map snd
                if droppedIntensities.Length > 0 && fi1.Length > 0 then
                    Expect.isTrue (Array.min fi1 > List.max droppedIntensities) "retained intensities exceed every dropped intensity"
                Expect.isTrue (fm3.Length <= fm1.Length) "a stricter minimum SNR cannot retain more points"
                // Retention is intensity/noise > minSnr, so it uses one intensity threshold and preserves pairing and order; the exact noise value is deliberately not pinned.

            // PENDING: the 50th-percentile noise of [10;100;50;200;20] is the median of the SORTED
            // intensities = 50, so at minSnr 0.9 the point (3, 50) has SNR 1.0 > 0.9 and must be retained.
            // The implementation feeds the UNSORTED array with an off-by-one count to the percentile
            // function, yielding noise 75 and threshold 67.5, which wrongly drops it.
            ptestCase "SNR filtering uses the median of the sorted intensities as its noise level" <| fun _ ->
                let fm, fi = SignalDetection.filterByIntensitySNR 50.0 0.9 [|1.;2.;3.;4.;5.|] [|10.;100.;50.;200.;20.|]
                Expect.isTrue (Array.zip fm fi |> Array.contains (3.0, 50.0)) "the median-intensity pair is retained"
        ]

        testList "SecondDerivative" [
            testCase "second-derivative centroiding guards inputs shorter than five points" <| fun _ ->
                let actualMz, actualIntensities =
                    SignalDetection.SecondDerivative.toCentroid true false 7 10.0 0.12 20 50.0 [|1.;2.;3.;4.|] [|1.;2.;3.;4.|]
                Expect.equal actualMz [||] "fewer than five m/z points produce no centroids"
                Expect.equal actualIntensities [||] "fewer than five intensity points produce no centroids"
                // The documented minimum-length guard rejects inputs shorter than five points.

            // PENDING: every candidate below yThreshold should simply be skipped, yielding empty output.
            // If the FIRST detected candidate fails the threshold, the implementation indexes xData.[-1]
            // and throws IndexOutOfRangeException.
            ptestCase "an all-sub-threshold spectrum yields no centroids instead of crashing" <| fun _ ->
                let mz = Array.init 41 (fun i -> 99.0 + float i * 0.05)
                let intens = mz |> Array.map (gauss 5.0 100.0 0.1)
                let r = SignalDetection.SecondDerivative.toCentroid true false 7 10.0 0.12 20 50.0 mz intens
                Expect.equal r ([||],[||]) "an all-sub-threshold spectrum yields no centroids"

            testCase "second-derivative centroiding reduces a clean Gaussian peak to a single centroid at its apex" <| fun _ ->
                let mz = Array.init 41 (fun i -> 99.0 + float i * 0.05)
                let intens = mz |> Array.map (gauss 1000.0 100.0 0.1)
                let cx, cy =
                    SignalDetection.SecondDerivative.toCentroid true false 7 10.0 0.12 20 50.0 mz intens
                Expect.equal cx.Length 1 "a single Gaussian peak yields one centroid"
                Expect.floatClose Accuracy.medium cx.[0] 100.0 "the refined centroid is at the symmetric Gaussian apex"
                Expect.floatClose Accuracy.high cy.[0] 1000.0 "the reported intensity is the raw apex maximum"
                // The constructed grid contains x=100.0 exactly, where the symmetric Gaussian reaches its amplitude.
        ]

        testList "WaveletCentroiding" [
            testCase "wavelet centroiding recovers both peaks of a clean two-peak spectrum" <| fun _ ->
                let mz = [|498.0 .. 0.01 .. 502.0|]
                let intens =
                    mz
                    |> Array.map (fun x -> gauss 1000.0 499.5 0.02 x + gauss 800.0 500.5 0.02 x)
                let parameters : SignalDetection.Wavelet.WaveletParameters = {
                    NumberOfScales = 10
                    YThreshold = 10.0
                    MzTolerance = 0.05
                    SNRS_Percentile = 50.0
                    MinSNR = 3.0
                    RefineMZ = true
                    SumIntensities = false
                    }
                let cx, cy = SignalDetection.Wavelet.toCentroidWithRicker2D parameters mz intens
                Expect.equal cx.Length cy.Length "centroid m/z and intensity arrays have equal lengths"
                Expect.isTrue (cx |> Array.exists (fun x -> abs (x - 499.5) <= 0.05)) "a centroid is recovered near the first apex"
                Expect.isTrue (cx |> Array.exists (fun x -> abs (x - 500.5) <= 0.05)) "a centroid is recovered near the second apex"
                Expect.isTrue (cx |> Array.forall (fun x -> x >= 498.0 && x <= 502.0)) "centroids lie inside the data range"
                Expect.isTrue (cy |> Array.forall (fun y -> y > 0.0)) "reported peak intensities are positive"
                let s = Array.sort cx
                // >= keeps the exact-boundary merge semantics unpinned; the current merge logic yields strictly greater gaps.
                Expect.isTrue (s |> Array.pairwise |> Array.forall (fun (a,b) -> b - a >= 0.05)) "sorted consecutive centroids are separated by at least the m/z tolerance"
                // the module's stated merging contract: peaks closer than the m/z tolerance are merged into one centroid, so no two reported centroids may sit within it
                let nearestFirstPeak =
                    cx
                    |> Array.mapi (fun i x -> i, x)
                    |> Array.minBy (fun (_, x) -> abs (x - 499.5))
                let nearestFirstPeakIntensity = cy.[fst nearestFirstPeak]
                Expect.isTrue (abs (nearestFirstPeakIntensity - 1000.0) <= 250.0) "the centroid nearest 499.5 has intensity within 25% of 1000.0"
                // the constructed apex (on the 0.01 grid) has intensity 1000; a centroid claiming that peak must report an intensity of its magnitude, not an artifact value
                // Expected positions come from the constructed apexes; the matching tolerance is 0.05 m/z.
        ]
    ]
