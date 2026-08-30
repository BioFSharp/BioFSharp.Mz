module SparsePeakArrayTests

open Expecto
open BioFSharp.Mz

[<Tests>]
let tests =
    testList "SparsePeakArrayTests" [
        testCase "initMzToBinIdx and initBinIdxToMz are consistent for exact bin centers" <| fun _ ->
            let width = 0.5
            let offset = 0.4
            let binIdx = SparsePeakArray.initMzToBinIdx width offset 100.3
            let mz = SparsePeakArray.initBinIdxToMz width offset binIdx
            // Exact bin-center conversion must preserve the expected index and m/z.
            Expect.equal binIdx 201 "100.3 maps to bin 201"
            Expect.floatClose Accuracy.high mz 100.3 "bin 201 maps back to m/z 100.3"

        testCase "peaksToNearestBinVector sums intensities of peaks falling in the same bin" <| fun _ ->
            let peaks = [| Peak(10.2, 5.0); Peak(10.4, 3.0); Peak(12.7, 2.0) |]
            let result = SparsePeakArray.peaksToNearestBinVector 1.0 0.0 0.0 1000.0 peaks
            let bin10 = SparsePeakArray.initMzToBinIdx 1.0 0.0 10.2
            let bin12 = SparsePeakArray.initMzToBinIdx 1.0 0.0 12.7
            // Peaks sharing a bin are aggregated by summing their intensities.
            Expect.equal result.Data.Count 2 "only two occupied bins are present"
            Expect.isTrue (result.Data.ContainsKey bin10) "bin 10 is present"
            Expect.isTrue (result.Data.ContainsKey bin12) "bin 12 is present"
            Expect.floatClose Accuracy.high result.Data.[bin10] 8.0 "bin 10 contains the summed intensity"
            Expect.floatClose Accuracy.high result.Data.[bin12] 2.0 "bin 12 contains its peak intensity"

        testCase "peaksToNearestBinVector excludes peaks at or beyond the mass borders" <| fun _ ->
            let peaks = [| Peak(10.0, 1.0); Peak(20.0, 2.0); Peak(15.0, 4.0) |]
            let result = SparsePeakArray.peaksToNearestBinVector 1.0 0.0 10.0 20.0 peaks
            let bin15 = SparsePeakArray.initMzToBinIdx 1.0 0.0 15.0
            // The lower and upper mass borders are strict exclusion boundaries.
            Expect.equal result.Data.Count 1 "only the interior peak's bin is present"
            Expect.isTrue (result.Data.ContainsKey bin15) "the interior peak's bin is present"
            Expect.floatClose Accuracy.high result.Data.[bin15] 4.0 "the interior peak is retained"

        testCase "dot multiplies matching bins and ignores disjoint ones" <| fun _ ->
            let x =
                SparsePeakArray.peaksToNearestBinVector 1.0 0.0 0.0 1000.0 [| Peak(10.2, 2.0); Peak(12.7, 3.0) |]
            let y =
                SparsePeakArray.peaksToNearestBinVector 1.0 0.0 0.0 1000.0 [| Peak(10.3, 4.0); Peak(30.0, 7.0) |]
            let disjointX =
                SparsePeakArray.peaksToNearestBinVector 1.0 0.0 0.0 1000.0 [| Peak(100.0, 2.0) |]
            let disjointY =
                SparsePeakArray.peaksToNearestBinVector 1.0 0.0 0.0 1000.0 [| Peak(200.0, 4.0) |]
            // Only bin 10 is shared, so the dot product is 2.0 * 4.0.
            Expect.floatClose Accuracy.high (SparsePeakArray.dot x y) 8.0 "matching bins are multiplied and summed"
            Expect.floatClose Accuracy.high (SparsePeakArray.dot disjointX disjointY) 0.0 "disjoint bins contribute nothing"
    ]
