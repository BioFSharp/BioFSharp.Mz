module ChargeStateTests

open System
open Expecto
open BioFSharp.Mz

let stdParams = ChargeState.createChargeDetermParams 1 4 1.1 0.05 0.05 20

let expectFloatListClose accuracy (actual: float list) (expected: float list) message =
    Expect.equal actual.Length expected.Length (sprintf "%s length" message)
    List.iter2
        (fun actualValue expectedValue -> Expect.floatClose accuracy actualValue expectedValue message)
        actual
        expected

let expectFloatArrayClose accuracy (actual: float []) (expected: float []) message =
    Expect.equal actual.Length expected.Length (sprintf "%s length" message)
    Array.iter2
        (fun actualValue expectedValue -> Expect.floatClose accuracy actualValue expectedValue message)
        actual
        expected

let expectPeakPairsClose tolerance (actual: (float * float) list) (expected: (float * float) list) message =
    let actualSorted = actual |> List.sortBy fst
    let expectedSorted = expected |> List.sortBy fst
    Expect.equal actualSorted.Length expectedSorted.Length (sprintf "%s length" message)
    List.iter2
        (fun (actualMz, actualIntensity) (expectedMz, expectedIntensity) ->
            Expect.isTrue
                (abs (actualMz - expectedMz) <= tolerance)
                (sprintf "%s m/z; expected %g, got %g" message expectedMz actualMz)
            Expect.isTrue
                (abs (actualIntensity - expectedIntensity) <= tolerance)
                (sprintf "%s intensity; expected %g, got %g" message expectedIntensity actualIntensity))
        actualSorted
        expectedSorted

let makeAssignedCharge score (peaks: Peak list) =
    ChargeState.createAssignedCharge
        "p"
        "f"
        500.
        2
        999.
        0.
        score
        []
        peaks.Length
        1.
        (Set.ofList peaks)
        None

[<Tests>]
let tests =
    testList "ChargeStateTests" [
        testList "ChargeSelection" [
            testCase "getChargeBy maps observed isotope spacing to the charge state" <| fun _ ->
                Expect.equal (ChargeState.getChargeBy stdParams 1.0) 1 "1.0 Da spacing identifies charge 1"
                Expect.equal (ChargeState.getChargeBy stdParams 0.5) 2 "0.5 Da spacing identifies charge 2"
                Expect.equal (ChargeState.getChargeBy stdParams 0.25) 4 "0.25 Da spacing identifies charge 4"
                // Adjacent isotopologue peaks are spaced ~1 Da divided by the charge (published isotope-envelope relation), so mean spacings 1.0, 0.5, 0.25 identify charges 1, 2, 4.

            testCase "getChargeBy selects the nearest allowed charge for imperfect spacing" <| fun _ ->
                Expect.equal (ChargeState.getChargeBy stdParams 0.34) 3 "0.34 Da spacing is nearest to charge 3"
                // real centroid spacings are never exact: |0.34 - 1/3| = 0.0067 versus 0.09 for charge 4 and 0.16 for charge 2 - nearest-spacing selection is the required behavior, hand-computed.

            testCase "mzChargeDeviationBy is zero for exact spacings, permutation-invariant, and grows with deviation" <| fun _ ->
                Expect.floatClose
                    Accuracy.high
                    (ChargeState.mzChargeDeviationBy [0.5; 0.5] 0.5)
                    0.0
                    "exact agreement must yield zero under any deviation measure"
                Expect.floatClose
                    Accuracy.high
                    (ChargeState.mzChargeDeviationBy [0.6; 0.4] 0.5)
                    (ChargeState.mzChargeDeviationBy [0.4; 0.6] 0.5)
                    "a deviation measure over a set of spacings cannot depend on their order"
                Expect.isTrue
                    (ChargeState.mzChargeDeviationBy [0.6; 0.4] 0.5 > ChargeState.mzChargeDeviationBy [0.55; 0.45] 0.5)
                    "larger deviations from the theoretical spacing must score larger"
                // Measure axioms (identity, symmetry, monotonicity); no exact formula value is pinned.

            testCase "getScore rewards explained peaks and penalizes deviation" <| fun _ ->
                Expect.isTrue
                    (ChargeState.getScore 3 5 0.2 > ChargeState.getScore 3 5 0.1)
                    "higher deviation is worse"
                Expect.isTrue
                    (ChargeState.getScore 4 5 0.1 < ChargeState.getScore 2 5 0.1)
                    "a subset explaining more of the source peaks scores better, i.e. lower"
                // Documented scoring intent - lower is better; direction only, no constants pinned.
        ]

        testList "ClusterExtraction" [
            testCase "getRelPeakPosInWindowBy returns start-normalized in-window peaks" <| fun _ ->
                let mz = [|100.; 100.3; 100.6; 102.|]
                let intensity = [|1000.; 800.; 600.; 500.|]
                let startInt, cluster = ChargeState.getRelPeakPosInWindowBy mz intensity 1.0 0.1 0.1 0
                Expect.floatClose Accuracy.high startInt 1000.0 "the start intensity is returned"
                let actualPairs = cluster.Peaks |> List.map (fun peak -> peak.Mz, peak.Intensity)
                expectPeakPairsClose
                    1e-9
                    actualPairs
                    [(0.3, 0.8); (0.6, 0.6)]
                    "relative in-window peaks"
                Expect.equal cluster.SourceSetLength 3 "the source length counts the start peak"
                Expect.equal cluster.SubSetLength 3 "the subset length counts the start peak"
                // Peaks within the 1.0-wide window right of the start peak (100.3 and 100.6; 102 is outside) are represented relative to the start peak's m/z and normalized by its intensity: (100.3-100, 800/1000) and (100.6-100, 600/1000), hand-computed. Both pass the intensity thresholds by construction - the checks compare RAW intensities: 800 > 0.1*1000 (vs start) and 800 > 0.1*1000 (vs prior); 600 > 0.1*1000 and 600 > 0.1*800. Lengths count the start peak plus accepted peaks.
                // The implementation returns the accepted peaks in descending order; comparison above is set/order agnostic as specified.

            testCase "intensity thresholds reject noise peaks during cluster extraction" <| fun _ ->
                let _, rejected = ChargeState.getRelPeakPosInWindowBy [|500.0; 500.4|] [|1000.0; 40.0|] 1.0 0.05 0.01 0
                Expect.isEmpty rejected.Peaks "a peak below the start-relative threshold is rejected"
                Expect.equal rejected.SourceSetLength 1 "the rejected peak is omitted from the source count"
                Expect.equal rejected.SubSetLength 1 "the rejected peak is omitted from the subset count"

                let _, acceptedThenRejected =
                    ChargeState.getRelPeakPosInWindowBy
                        [|500.0; 500.4; 500.8|]
                        [|1000.0; 800.0; 300.0|]
                        1.0
                        0.05
                        0.5
                        0
                Expect.equal acceptedThenRejected.Peaks.Length 1 "only the first candidate peak passes both thresholds"
                let acceptedPeak = List.head acceptedThenRejected.Peaks
                Expect.floatClose Accuracy.high acceptedPeak.Mz 0.4 "the accepted peak has start-relative m/z"
                Expect.floatClose Accuracy.high acceptedPeak.Intensity 0.8 "the accepted peak has start-relative intensity"
                Expect.equal acceptedThenRejected.SourceSetLength 2 "the accepted count includes the start peak"
                Expect.equal acceptedThenRejected.SubSetLength 2 "the subset count includes the start peak"
                // the two thresholds guard against noise creating false charge hypotheses; each case isolates one threshold.

            testCase "powerSetOf enumerates all subsets anchored at the start peak" <| fun _ ->
                let cluster = ChargeState.createPutativeIsotopeCluster [Peak(0.1, 0.5); Peak(0.2, 0.3)] 3 3
                let subsets = ChargeState.powerSetOf cluster
                Expect.equal subsets.Length 4 "2^2, powerset combinatorics"
                subsets
                |> List.iter (fun subset ->
                    Expect.isTrue
                        (subset.Peaks |> List.exists (fun peak -> peak.Mz = 0.0 && peak.Intensity = 1.0))
                        "every subset contains the start-peak anchor"
                    Expect.equal subset.SourceSetLength 3 "source length is preserved"
                    Expect.equal subset.SubSetLength subset.Peaks.Length "subset length matches its peak list")
                let canonicalSubsets =
                    subsets
                    |> List.map (fun subset -> subset.Peaks |> List.map (fun peak -> peak.Mz) |> List.sort)
                    |> Set.ofList
                Expect.equal
                    canonicalSubsets
                    (set [ [0.0]; [0.0; 0.1]; [0.0; 0.2]; [0.0; 0.1; 0.2] ])
                    "all distinct start-anchored membership combinations are present"
                // 2^n subsets of n non-anchor peaks is the definition of a powerset; the anchor represents the start peak, which every isotope-cluster hypothesis must contain.
                // count alone passes for duplicated subsets; distinct membership combinations are what the isotope-hypothesis search needs.

            testCase "mzDistancesOf returns adjacent gaps of a descending peak list" <| fun _ ->
                let actual = ChargeState.mzDistancesOf [Peak(3.0, 1.); Peak(1.0, 1.); Peak(0.0, 1.)]
                expectFloatListClose Accuracy.high actual [1.0; 2.0] "adjacent m/z gaps"
                // The output [1.0; 2.0] is in REVERSE adjacency order (head is the gap between the last two peaks, 1.0-0.0); the previous "ascending-position order" wording was wrong.
        ]

        testList "EmpiricalStatistics" [
            testCase "empiricalRightPValueOf counts the fraction of simulated values at or below the score" <| fun _ ->
                let dist = [|0.1; 0.2; 0.2; 0.5|]
                Expect.equal (ChargeState.empiricalRightPValueOf dist 0.05) 0.0 "0 of 4 are at or below 0.05"
                Expect.equal (ChargeState.empiricalRightPValueOf dist 0.2) 0.75 "3 of 4 are at or below 0.2"
                Expect.equal (ChargeState.empiricalRightPValueOf dist 0.5) 1.0 "all samples are at or below 0.5"
                Expect.equal (ChargeState.empiricalRightPValueOf dist 0.6) 1.0 "all samples are at or below 0.6"
                // Empirical CDF over 4 sorted samples: 0 of 4 are <= 0.05, 3 of 4 are <= 0.2, all are <= 0.5 - hand-counted.

            testCase "poissonProb is the Poisson probability mass function" <| fun _ ->
                Expect.floatClose Accuracy.high (ChargeState.poissonProb 2.0 0.0) (exp -2.0) "P(0) = e^-2"
                Expect.floatClose Accuracy.high (ChargeState.poissonProb 2.0 1.0) (2.0 * exp -2.0) "P(1) = 2e^-2"
                Expect.floatClose Accuracy.high (ChargeState.poissonProb 2.0 2.0) (2.0 * exp -2.0) "P(2) = 2e^-2"
                // Published Poisson PMF lambda^k e^-lambda / k!: for lambda=2, P(0)=e^-2, P(1)=2e^-2, P(2)=4e^-2/2=2e^-2.

            testCase "poissonEstofMassTrunc returns the normalized truncated Poisson distribution" <| fun _ ->
                let r = ChargeState.poissonEstofMassTrunc (fun _ -> 1.0) 3 999.0
                expectFloatArrayClose Accuracy.high r [|0.4; 0.4; 0.2|] "normalized truncated Poisson distribution"
                Expect.floatClose Accuracy.high (Array.sum r) 1.0 "truncated distribution sums to one"
                // Poisson(1) PMF over k = 0,1,2 is e^-1 * (1, 1, 1/2); normalizing (1, 1, 0.5)/2.5 = (0.4, 0.4, 0.2) - hand-computed from the published PMF; truncated normalization must sum to one.

            testCase "empiricalPValueOfSim delegates to the simulator and short-circuits trivial subsets" <| fun _ ->
                // mock-based contract test of the p-value entry point; the simulator is an injected dependency.
                let mockSim = fun (_: int * float) -> [|0.1;0.2;0.3;0.4|]
                Expect.floatClose
                    Accuracy.high
                    (ChargeState.empiricalPValueOfSim mockSim (3, 2.0) 0.25)
                    0.5
                    "2 of 4 simulated deviations are at or below 0.25"
                let anyScore = 0.25
                Expect.floatClose
                    Accuracy.high
                    (ChargeState.empiricalPValueOfSim mockSim (1, 2.0) anyScore)
                    1.0
                    "a subset of one peak carries no spacing information, so the p-value is 1.0"
                let mutable calls = 0
                let countingMock (_: int * float) =
                    calls <- calls + 1
                    [|0.1;0.2;0.3;0.4|]
                Expect.floatClose
                    Accuracy.high
                    (ChargeState.empiricalPValueOfSim countingMock (1, 2.0) anyScore)
                    1.0
                    "a single-peak subset short-circuits to a full p-value"
                Expect.equal calls 0 "the simulator is not invoked for a one-peak subset"

            testCase "KL divergence is zero for identical distributions, positive otherwise, None on length mismatch" <| fun _ ->
                match ChargeState.kullbackLeiblerDivergenceOf [|0.5; 0.5|] [|0.5; 0.5|] with
                | Some value -> Expect.floatClose Accuracy.high value 0.0 "identical distributions have zero divergence"
                | None -> failtest "identical distributions return Some divergence"
                match ChargeState.kullbackLeiblerDivergenceOf [|0.5; 0.5|] [|0.9; 0.1|] with
                | Some value ->
                    let expected = 0.9 * log (0.9 / 0.5) + 0.1 * log (0.1 / 0.5)
                    Expect.isTrue
                        (abs (value - expected) <= 1e-9)
                        (sprintf "the directed divergence has its defining value; expected %g, got %g" expected value)
                    Expect.isTrue (value > 0.0) "different distributions have positive divergence"
                | None -> failtest "same-length distributions return Some divergence"
                Expect.equal
                    (ChargeState.kullbackLeiblerDivergenceOf [|0.5; 0.5|] [|1.0|])
                    None
                    "mismatched supports return None"
                // Gibbs' inequality - D(p||q) >= 0 with equality iff p = q; mismatched supports are undefined (documented None).
                // the defining D(p||q) with q as first argument - the exact value pins direction and operand order, which zero/positive checks cannot.

            testCase "peakPosStdDevBy scales like a dispersion estimate" <| fun _ ->
                let mkAC residuals =
                    ChargeState.createAssignedCharge
                        "p"
                        "f"
                        500.0
                        2
                        999.0
                        0.0
                        1.0
                        residuals
                        (List.length residuals)
                        1.0
                        Set.empty
                        None
                let s1 = ChargeState.peakPosStdDevBy [mkAC [-0.1; 0.0; 0.1]]
                let s1' = ChargeState.peakPosStdDevBy [mkAC [0.1; -0.1; 0.0]]
                let s2 = ChargeState.peakPosStdDevBy [mkAC [-0.2; 0.0; 0.2]]
                Expect.isTrue (s1 > 0.0 && Double.IsFinite s1) "non-constant residuals have a finite positive spread"
                Expect.floatClose Accuracy.high s1' s1 "the spread is permutation-invariant"
                Expect.floatClose Accuracy.high s2 (2.0 * s1) "the spread is homogeneous under scaling"
                // any standard-deviation convention (sample or population) is positive for non-constant data, permutation-invariant, and homogeneous of degree 1 - convention-independent properties.
        ]

        testCase "heavier peptides get broader predicted isotope envelopes" <| fun _ ->
            let eLight = ChargeState.poissonEstofMassTrunc ChargeState.n14MassToLambda 10 500.0
            let eHeavy = ChargeState.poissonEstofMassTrunc ChargeState.n14MassToLambda 10 3000.0
            let eLightN15 = ChargeState.poissonEstofMassTrunc ChargeState.n15MassToLambda 10 500.0
            let eHeavyN15 = ChargeState.poissonEstofMassTrunc ChargeState.n15MassToLambda 10 3000.0
            [eLight; eHeavy; eLightN15; eHeavyN15]
            |> List.iter (fun envelope -> Expect.floatClose Accuracy.high (Array.sum envelope) 1.0 "the envelope is normalized")
            let expectedIndex envelope =
                envelope
                |> Array.mapi (fun i probability -> float i * probability)
                |> Array.sum
            Expect.isTrue
                (expectedIndex eHeavy > expectedIndex eLight)
                "the heavier n14 envelope has a greater expected isotopologue index"
            Expect.isTrue
                (expectedIndex eHeavyN15 > expectedIndex eLightN15)
                "the heavier n15 envelope has a greater expected isotopologue index"
            // more atoms -> more chances of a heavy isotope: the expected isotopologue index must grow with peptide mass - a domain invariant that leaves the empirical lambda coefficients unpinned.

        testList "CandidateHandling" [
            testCase "normalizePeaksByIntensitySum returns intensity fractions summing to one" <| fun _ ->
                let peaks = Set.ofList [Peak(100.0, 2.0); Peak(101.0, 3.0); Peak(102.0, 5.0)]
                let r = ChargeState.normalizePeaksByIntensitySum peaks
                Expect.floatClose Accuracy.high (Array.sum r) 1.0 "normalized intensities sum to one"
                expectFloatArrayClose Accuracy.high (Array.sort r) [|0.2; 0.3; 0.5|] "sorted normalized intensities"
                // Normalization by the total 10 preserves the 2:3:5 ratios - hand-computed; sorted comparison avoids pinning Set enumeration order.

            testCase "removeSubSetsOfBestHit drops subsets of the leading candidate on score-sorted input" <| fun _ ->
                let best = makeAssignedCharge 1.0 [Peak(500.0, 1.); Peak(500.5, 0.8); Peak(501.0, 0.5)]
                let subA = makeAssignedCharge 2.0 [Peak(500.0, 1.); Peak(500.5, 0.8)]
                let other = makeAssignedCharge 3.0 [Peak(505.0, 0.3)]
                let r = ChargeState.removeSubSetsOfBestHit [best; subA; other]
                Expect.isTrue (r |> List.exists (fun candidate -> candidate.Score = 1.0)) "the best candidate survives"
                Expect.isTrue (r |> List.exists (fun candidate -> candidate.Score = 3.0)) "the disjoint candidate survives"
                Expect.isFalse (r |> List.exists (fun candidate -> candidate.Score = 2.0)) "the subset candidate is removed"
                Expect.equal r.Length 2 "only the best and disjoint candidates remain"
                // the function treats the FIRST candidate as the best hit (it never consults Score), so this test supplies the score-sorted order that putativePrecursorChargeStatesBy emits; order-independent behavior is asserted by the pending test below.

            ptestCase "subset removal is anchored at the globally best candidate regardless of input order" <| fun _ ->
                // PENDING: the best hit is defined by the LOWEST score, not by list position. With the best
                // candidate not first, the implementation anchors at whatever comes first: subsets of the true
                // best survive. Argued-correct: reordering the input must not change which candidates are
                // removed.
                let best = makeAssignedCharge 1.0 [Peak(500.0, 1.); Peak(500.5, 0.8); Peak(501.0, 0.5)]
                let subA = makeAssignedCharge 2.0 [Peak(500.0, 1.); Peak(500.5, 0.8)]
                let other = makeAssignedCharge 3.0 [Peak(505.0, 0.3)]
                let r = ChargeState.removeSubSetsOfBestHit [subA; best; other]
                Expect.isTrue (r |> List.exists (fun candidate -> candidate.Score = 1.0)) "the global best survives"
                Expect.isTrue (r |> List.exists (fun candidate -> candidate.Score = 3.0)) "the disjoint candidate survives"
                Expect.isFalse (r |> List.exists (fun candidate -> candidate.Score = 2.0)) "the subset of the global best is removed"
                Expect.equal r.Length 2 "only the global best and the disjoint candidate remain"
        ]

        testList "RandomSimulation" [
            ptestCase "simulated isotope clusters respect the maximum cumulative distance" <| fun _ ->
                // PENDING: every peak of a simulated cluster must lie within maxDistance of the anchor - the
                // simulation models an isotope envelope of bounded width. At loopCounter = maxCount neither
                // overflow guard matches and an over-distance peak is accepted (probe: maxDistance 0.5
                // admitted a peak at Mz 2.0). Near-zero stdDev makes the outcome deterministic even though the
                // FSharp.Stats normal sampler is unseedable (only System.Random is injectable).
                let rnd = System.Random(42)
                let r = ChargeState.rndMzIntensityEntityCollectionBy rnd 0 1e-9 2 0.5 3
                match r with
                | Some cluster ->
                    Expect.isTrue
                        (cluster.Peaks |> List.forall (fun p -> p.Mz <= 0.5))
                        "every simulated peak must remain within maxDistance of the anchor"
                | None -> ()

            ptestCase "the charge simulation can produce the maximum expected charge" <| fun _ ->
                // PENDING: the simulated null distribution must cover every allowed charge including the
                // maximum; Random.Next's exclusive upper bound means ExpectedMaximumCharge is never drawn
                // (probe: 800 draws at maxCharge 4 produced charges 1-3 only), biasing empirical p-values.
                // With stdDev 1e-9 a drawn charge z gives spacings within 1e-6 of 1/z, so observing a spacing
                // near 0.25 is deterministic evidence of charge 4.
                let rnd = System.Random(7)
                let clusters =
                    [1..400]
                    |> List.choose (fun _ -> ChargeState.rndMzIntensityEntityCollectionBy rnd 100 1e-9 4 10.0 3)
                let hasMaximumChargeSpacing =
                    clusters
                    |> List.exists (fun cluster ->
                        cluster.Peaks
                        |> List.map (fun peak -> peak.Mz)
                        |> List.pairwise
                        |> List.map (fun (mzA, mzB) -> abs (mzA - mzB))
                        |> List.exists (fun spacing -> abs (spacing - 0.25) <= 0.01))
                Expect.isTrue
                    hasMaximumChargeSpacing
                    "at least one simulated cluster must contain a spacing near 0.25"
        ]

        testList "EndToEnd" [
            testCase "putativePrecursorChargeStatesBy identifies a charge-2 isotope cluster" <| fun _ ->
                let mz = [|500.0; 500.5; 501.0; 505.0|]
                let intensity = [|1000.; 800.; 500.; 100.|]
                let candidates = ChargeState.putativePrecursorChargeStatesBy stdParams mz intensity "precScan" "prodScan" 500.0
                Expect.isFalse (List.isEmpty candidates) "candidate list is non-empty"
                candidates
                |> List.iter (fun candidate ->
                    Expect.equal candidate.PrecursorMZ 500.0 "candidate precursor m/z is preserved"
                    Expect.equal candidate.PrecursorSpecID "precScan" "candidate precursor spectrum ID is preserved")
                let best = candidates |> List.minBy (fun candidate -> candidate.Score)
                Expect.equal best.PrecCharge 2 "the lowest-score candidate has charge 2"
                Expect.floatClose Accuracy.high best.PutMass 997.98545 "the charge-2 neutral mass"
                Expect.floatClose Accuracy.high best.MZChargeDev 0.0 "exact 0.5 Da spacing has zero charge deviation"
                let charges = candidates |> List.map (fun candidate -> candidate.PrecCharge)
                Expect.equal (Set.ofList charges).Count charges.Length "returned precursor charges are distinct"
                let chargeTwo = candidates |> List.find (fun candidate -> candidate.PrecCharge = 2)
                Expect.equal chargeTwo.SubSetLength 3 "the winning charge-2 candidate retains the full envelope"
                // The spectrum is constructed as a textbook 2+ isotope envelope (0.5 Da spacing); charge determination must prefer charge 2, and the derived neutral mass follows the fundamental m/z relation.
                // the function documents one best candidate per charge; the full three-peak envelope explains the most peaks and must be the retained charge-2 hypothesis.

            // PENDING: an uninterpretable dense window (>= 15 source peaks) makes the function fall
            // back to hedging across ALL allowed charges - one sentinel candidate per charge in
            // [ExpectedMinimalCharge .. ExpectedMaximumCharge]. Observed: charges [2;3;4] only; the
            // minimum charge 1 is omitted from the fallback.
            ptestCase "the dense-window fallback covers every allowed charge" <| fun _ ->
                let mz = Array.init 20 (fun i -> 499.5 + 0.05 * float i)
                let intensity = Array.create 20 100.0
                let candidates = ChargeState.putativePrecursorChargeStatesBy stdParams mz intensity "p" "f" 500.0
                let charges = candidates |> List.map (fun c -> c.PrecCharge) |> List.sort
                Expect.equal charges [1; 2; 3; 4] "one fallback candidate per allowed charge, including the minimum"
        ]
    ]
