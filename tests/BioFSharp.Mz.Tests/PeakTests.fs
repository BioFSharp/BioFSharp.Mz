module PeakTests

open Expecto
open BioFSharp.Mz

[<Tests>]
let tests =
    testList "PeakTests" [
        testList "Ions" [
            testCase "hasFlag detects a contained single ion flag" <| fun _ ->
                let combined = Ions.IonTypeFlag.B ||| Ions.IonTypeFlag.Y
                // IonTypeFlag values are distinct powers of two (B=8, Y=64, X=32); a contained flag shares a bit (AND <> 0), a non-contained flag shares none.
                Expect.isTrue (Ions.hasFlag combined Ions.IonTypeFlag.B) "B is contained in B ||| Y"
                Expect.isFalse (Ions.hasFlag combined Ions.IonTypeFlag.X) "X is not contained in B ||| Y"

            testCase "hasFlag uses any-overlap semantics for composite queries" <| fun _ ->
                // Domain contract: a composite query expresses set membership - "does this ion carry ANY of the
                // queried types (e.g. is it a b- or y-ion)?" All-bits semantics on a composite query would ask
                // "is this ion simultaneously b and y", which is domain nonsense. This contract is deliberate and
                // forward-looking: no production call site currently issues composite queries (they all pass a
                // single flag), so nothing downstream relies on the divergence from .NET Enum.HasFlag yet.
                Expect.isTrue (Ions.hasFlag (Ions.IonTypeFlag.B ||| Ions.IonTypeFlag.Y) (Ions.IonTypeFlag.X ||| Ions.IonTypeFlag.Y)) "B ||| Y overlaps X ||| Y through Y"
                Expect.isFalse (Ions.hasFlag (Ions.IonTypeFlag.B ||| Ions.IonTypeFlag.Y) (Ions.IonTypeFlag.X ||| Ions.IonTypeFlag.Z)) "B ||| Y does not overlap X ||| Z"

            testCase "createIonTypeList enumerates exactly the component flags of a combined value" <| fun _ ->
                let combined = Ions.IonTypeFlag.B ||| Ions.IonTypeFlag.lossH2O
                let actual = Ions.createIonTypeList combined |> Set.ofSeq
                // flags are independent bits, so decomposing a combined value must list exactly its components.
                Expect.equal actual (set [Ions.IonTypeFlag.B; Ions.IonTypeFlag.lossH2O]) "combined B and lossH2O decomposes into exactly those flags"
        ]

        testList "TaggedMassAndPeak" [
            testCase "createTaggedH2OLoss adds the water-loss flag and preserves the mass" <| fun _ ->
                let tm = TaggedMass.createTaggedH2OLoss Ions.IonTypeFlag.B 445.5
                // B=8 and lossH2O=256 are disjoint bits, so the enum addition 8+256 sets exactly both flags.
                Expect.floatClose Accuracy.high tm.Mass 445.5 "mass is preserved"
                Expect.isTrue (Ions.hasFlag tm.Iontype Ions.IonTypeFlag.B) "B flag is preserved"
                Expect.isTrue (Ions.hasFlag tm.Iontype Ions.IonTypeFlag.lossH2O) "lossH2O flag is added"
                Expect.isFalse (Ions.hasFlag tm.Iontype Ions.IonTypeFlag.lossNH3) "lossNH3 flag is not added"

            testCase "createTaggedNH3Loss adds the ammonia-loss flag and preserves the mass" <| fun _ ->
                let tm = TaggedMass.createTaggedNH3Loss Ions.IonTypeFlag.Y 512.3
                // Y=64 and lossNH3=512 are disjoint bits, so 64+512 sets exactly both flags.
                Expect.floatClose Accuracy.high tm.Mass 512.3 "mass is preserved"
                Expect.isTrue (Ions.hasFlag tm.Iontype Ions.IonTypeFlag.Y) "Y flag is preserved"
                Expect.isTrue (Ions.hasFlag tm.Iontype Ions.IonTypeFlag.lossNH3) "lossNH3 flag is added"
                Expect.isFalse (Ions.hasFlag tm.Iontype Ions.IonTypeFlag.lossH2O) "lossH2O flag is not added"

            testCase "createTaggedPeakOf places the peak at the tagged mass with predicted intensity" <| fun _ ->
                let tm = TaggedMass.createTaggedMass Ions.IonTypeFlag.Y 300.0
                let predictor = fun flag -> if Ions.hasFlag flag Ions.IonTypeFlag.Y then 42.0 else 0.0
                let tp = TaggedPeak.createTaggedPeakOf tm predictor
                // documented contract - the theoretical peak sits at the tagged mass and its intensity is the prediction for its ion type; the predictor here answers 42.0 exactly for Y-flagged input.
                Expect.floatClose Accuracy.high tp.Mz 300.0 "m/z is the tagged mass"
                Expect.floatClose Accuracy.high tp.Intensity 42.0 "intensity is the predicted intensity"
                Expect.equal tp.Iontype tm.Iontype "ion type is preserved"
        ]

        testList "PeakArray" [
            testCase "zip pairs mz and intensity arrays elementwise" <| fun _ ->
                let arr = PeakArray.zip [|100.;200.;300.|] [|1.;2.;3.|]
                // documented elementwise pairing; values are copied unchanged.
                Expect.equal arr.Length 3 "zipped array has three peaks"
                Expect.isTrue (arr.[0].Mz = 100. && arr.[0].Intensity = 1.) "first pair is copied unchanged"
                Expect.isTrue (arr.[2].Mz = 300. && arr.[2].Intensity = 3.) "last pair is copied unchanged"

            testCase "zipMzInt agrees with zip" <| fun _ ->
                let pairs = [|(100.,1.);(200.,2.)|]
                let fromPairs = PeakArray.zipMzInt pairs
                let fromSeparateArrays = PeakArray.zip (Array.map fst pairs) (Array.map snd pairs)
                // both construct one Peak per (mz,intensity) pair - cross-consistency between the two constructors.
                Expect.equal fromSeparateArrays.Length fromPairs.Length "constructors produce arrays of the same length"
                Expect.isTrue (Array.forall2 (fun (p: Peak) (q: Peak) -> p.Mz = q.Mz && p.Intensity = q.Intensity) fromPairs fromSeparateArrays) "constructors produce identical peaks elementwise"

            testCase "zip throws on unequal input lengths" <| fun _ ->
                // no truncation/padding is documented; Array.map2 requires equal lengths.
                Expect.throws (fun () -> PeakArray.zip [|1.;2.|] [|1.|] |> ignore) "unequal input lengths throw"

            testCase "binToUpperIntergerMass bins to the next upper integer mass and keeps the maximum on collision" <| fun _ ->
                let pk = PeakArray.zipMzInt [|(100.1,3.);(100.2,5.);(100.9,4.);(102.0,7.)|]
                let bins = PeakArray.binToUpperIntergerMass pk 100 104
                // bin index = ceil(mz) - minMassBoarder per the documented rule.
                Expect.floatClose Accuracy.high bins.[1] 5.0 "stronger colliding peak is retained" // three peaks collide in bin 101 (ceil 100.1 = ceil 100.2 = ceil 100.9 = 101 -> index 1) with the maximum in the MIDDLE of the input order: first-write-wins would give 3, last-write-wins 4, sum 12 - only true max-retention yields 5 (hand-computed)
                Expect.floatClose Accuracy.high bins.[2] 7.0 "102.0 is placed in bin 102" // ceil 102.0 = 102 -> index 2
                Expect.floatClose Accuracy.high bins.[0] 0.0 "no peak ceils to 100" // no peak ceils to 100

            testCase "peaksToNearestUnitDaltonBin bins to the nearest dalton and keeps the maximum on collision" <| fun _ ->
                let pk = PeakArray.zipMzInt [|(100.4,2.);(101.6,1.);(101.7,4.);(102.4,2.)|]
                let bins = PeakArray.peaksToNearestUnitDaltonBin pk 100 104
                // nearest-integer rounding rule from the doc, hand-computed; inputs avoid .5 midpoints so rounding mode is irrelevant.
                Expect.floatClose Accuracy.high bins.[0] 2.0 "100.4 rounds to 100" // round 100.4 = 100 -> index 0
                Expect.floatClose Accuracy.high bins.[2] 4.0 "stronger colliding peak is retained" // three peaks collide in bin 102 (round 101.6 = round 101.7 = round 102.4 = 102) with the maximum in the middle: first-wins gives 1, last-wins 2, sum 7 - only max-retention yields 4 (hand-computed; no .5 midpoints, so rounding mode is irrelevant)
                Expect.floatClose Accuracy.high bins.[1] 0.0 "nothing rounds to 101" // nothing rounds to 101

            testCase "peaks outside the mass borders are silently dropped" <| fun _ ->
                let pk = PeakArray.zipMzInt [|(50.0,9.);(500.0,9.)|]
                let a = PeakArray.binToUpperIntergerMass pk 100 104
                let b = PeakArray.peaksToNearestUnitDaltonBin pk 100 104
                // peaks far below/above the [minMassBoarder, maxMassBoarder] window must not contribute to any bin (documented filtering); the guard protects the production path where full measured spectra are binned against scan limits. Contribution is asserted via all-zero bins, deliberately NOT via the array length, whose off-by-one is under dispute (see the pending border test).
                Expect.isTrue (Array.forall (fun v -> v = 0.0) a) "peaks outside the mass borders do not contribute to upper-integer bins"
                Expect.isTrue (Array.forall (fun v -> v = 0.0) b) "peaks outside the mass borders do not contribute to nearest-dalton bins"

            // PENDING: the vector variant's doc says peaks with mz > maxMassBoarder are filtered out, i.e.
            // equality with the border is kept. The current allocation (maxMassBoarder - minMassBoarder bins)
            // and guard `index < maxIndex-1` drop a peak whose ceiling equals the upper border, and make the
            // border bin unrepresentable. Argued-correct behavior: bins cover [min, max] inclusively
            // (length max-min+1) and the border peak lands in the last bin.
            ptestCase "a peak at the upper mass border is retained" <| fun _ ->
                let bins = PeakArray.binToUpperIntergerMass (PeakArray.zipMzInt [|(104.0,6.)|]) 100 104
                Expect.equal bins.Length 5 "upper border is represented in the final bin"
                Expect.floatClose Accuracy.high bins.[4] 6.0 "a peak at the upper mass border is retained"

            testCase "peaksToNearestUnitDaltonBinVector agrees with the array version for integral borders" <| fun _ ->
                let pk = PeakArray.zipMzInt [|(100.4,2.);(101.6,1.);(101.7,4.);(102.4,2.)|]
                let v = PeakArray.peaksToNearestUnitDaltonBinVector pk 100.0 104.0
                let a = PeakArray.peaksToNearestUnitDaltonBin pk 100 104
                // both implement the same documented rule; the vector variant truncates float borders to int first, so integral borders must agree exactly.
                Expect.equal v.Length a.Length "array and vector have the same length"
                Expect.isTrue (Array.forall (fun i -> a.[i] = v.[i]) [|0 .. a.Length - 1|]) "array and vector values agree elementwise"

            // PENDING: PeakArray.unzip iterates `for i = 0 to n` where n is the input length (one past the end) and copies whole Peak values instead of .Mz/.Intensity, so it always throws IndexOutOfRangeException and its return type is not float arrays. Documented behavior is returning the mz and intensity arrays. Values are asserted through a generic float cast so the test compiles against the current (wrong) generic return type, fails on any non-float element type, and passes only for the true round trip.
            ptestCase "unzip returns the mz and intensity arrays of a zipped array" <| fun _ ->
                let arr = PeakArray.zip [|100.;200.|] [|1.;2.|]
                let a, b = PeakArray.unzip arr
                let mzs = a |> Seq.cast<float> |> List.ofSeq
                let ints = b |> Seq.cast<float> |> List.ofSeq
                Expect.equal mzs [100.;200.] "unzip returns the original mz values"
                Expect.equal ints [1.;2.] "unzip returns the original intensity values"
        ]

        testList "PeakList" [
            testCase "PeakList unzip inverts zip" <| fun _ ->
                let mzs, ints = PeakList.unzip (PeakList.zip [100.;200.;300.] [10.;20.;30.])
                // zip and unzip are documented inverses; values are copied unchanged so exact float equality holds.
                Expect.isTrue (mzs = [100.;200.;300.] && ints = [10.;20.;30.]) "unzip returns the original m/z and intensity lists"

            testCase "PeakList.zipMzInt agrees with PeakList.zip" <| fun _ ->
                let pairs = [(100.,1.);(200.,2.)]
                let fromPairs = PeakList.zipMzInt pairs
                let fromSeparateLists = PeakList.zip (List.map fst pairs) (List.map snd pairs)
                // cross-consistency; both construct one Peak per pair.
                Expect.equal fromSeparateLists.Length fromPairs.Length "constructors produce lists of the same length"
                Expect.isTrue (List.forall2 (fun (p: Peak) (q: Peak) -> p.Mz = q.Mz && p.Intensity = q.Intensity) fromPairs fromSeparateLists) "constructors produce identical peaks elementwise"

            testCase "PeakList.zip throws on unequal input lengths" <| fun _ ->
                // List.map2 requires equal lengths; no truncation documented.
                Expect.throws (fun () -> PeakList.zip [1.;2.] [1.] |> ignore) "unequal input lengths throw"
        ]
    ]
