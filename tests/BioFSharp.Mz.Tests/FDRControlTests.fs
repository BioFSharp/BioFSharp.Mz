module FDRControlTests

open System
open Expecto
open FSharp.Stats
open BioFSharp.Mz

let private expectWithin tolerance actual expected message =
    Expect.isTrue
        (abs (actual - expected) <= tolerance)
        (sprintf "%s; expected %g, got %g" message expected actual)

let private expectFloatClose actual expected message =
    Expect.floatClose Accuracy.high actual expected message

[<Tests>]
let tests =
    testList "FDRControlTests" [
        testList "MAYU" [
            testCase "estimatePi0HG is the hypergeometric expectation of false positives" <| fun _ ->
                expectWithin
                    1e-9
                    (FDRControl.MAYU.estimatePi0HG 4.0 2.0 2.0)
                    1.5
                    "estimatePi0HG matches the hypergeometric expectation"
                // hand-computed from the MAYU/ProteinFDREstimator-specific candidate construction (fp-dependent population split with renormalized weights) - the arithmetic is pure combinatorics, but the model is MAYU's, not the textbook hypergeometric distribution.

            testCase "estimatePi0HG is bounded by the decoy count and zero decoys give zero" <| fun _ ->
                let estimate = FDRControl.MAYU.estimatePi0HG 100.0 90.0 10.0
                Expect.isTrue
                    (estimate >= 0.0 && estimate <= 10.0)
                    "the expected number of false positives cannot exceed the number of observed decoy hits (cf) nor be negative - probability bounds"
                Expect.equal
                    (FDRControl.MAYU.estimatePi0HG 100.0 90.0 0.0)
                    0.0
                    "no decoy hits -> zero expected false positives; the implementation reaches this via its non-finite guard (all-zero probabilities -> 0/0 -> NaN -> 0), not a documented zero-draw branch"

            testCase "the large-count Stirling branch stays finite and bounded" <| fun _ ->
                let r = FDRControl.MAYU.estimatePi0HG 10000.0 9000.0 10.0
                Expect.isTrue
                    (not (Double.IsNaN r) && not (Double.IsInfinity r))
                    (sprintf "the large-count expectation is finite; observed %g" r)
                Expect.isTrue (0.0 < r && r < 10.0) (sprintf "the expectation is strictly inside its endpoints; observed %g" r)
                // above the exact-factorial cutoff the implementation switches to Stirling's approximation; for a valid nondegenerate candidate distribution every count 0..10 has positive weight, so the expectation lies strictly inside the endpoints - production databases exceed 1000 entries routinely.
        ]

        testList "Binning" [
            testCase "binningFunction estimates PEP and q of 1 for bins holding equal target and decoy counts" <| fun _ ->
                let inputs =
                    [
                        FDRControl.createQValueInput 0.1 false
                        FDRControl.createQValueInput 0.9 true
                        FDRControl.createQValueInput 1.1 true
                        FDRControl.createQValueInput 1.9 false
                    ]
                    |> Array.ofList
                let scores, peps, qVals =
                    FDRControl.binningFunction
                        1.0
                        1.0
                        (fun (x: FDRControl.QValueInput) -> x.Score)
                        (fun (x: FDRControl.QValueInput) -> x.IsDecoy)
                        inputs
                let actualScores = scores |> Seq.toArray
                let actualPeps = peps |> Seq.toArray
                let actualQVals = qVals |> Seq.toArray
                Expect.equal actualScores.Length 2 "the two bandwidth-1 bins produce two score entries"
                Array.iter2
                    (fun actual expected -> expectFloatClose actual expected "the representative bin score")
                    actualScores
                    [|0.5; 1.5|]
                Array.iter2
                    (fun actual expected -> expectFloatClose actual expected "the bin PEP")
                    actualPeps
                    [|1.0; 1.0|]
                Array.iter2
                    (fun actual expected -> expectFloatClose actual expected "the cumulative q-value")
                    actualQVals
                    [|1.0; 1.0|]
                // with pi0 = 1, a score bin containing exactly as many decoys as targets carries no discriminating signal, so its local false-discovery estimate (PEP) is 1; the same holds cumulatively (q) when every bin looks like that - the target-decoy FDR reading of the crafted data, not an implementation echo.

            ptestCase "a decoy-free dataset yields zero estimated error rates" <| fun _ ->
                // PENDING: with no decoys there is no evidence of noise, so the local and cumulative FDR
                // estimates of every bin should be 0. The implementation divides by the zero decoy total
                // (1/0 = infinity scaling) and returns NaN everywhere.
                let inputs =
                    [|
                        FDRControl.createQValueInput 0.2 false
                        FDRControl.createQValueInput 0.7 false
                        FDRControl.createQValueInput 1.2 false
                        FDRControl.createQValueInput 1.7 false
                    |]
                let _, peps, qs =
                    FDRControl.binningFunction
                        1.0
                        1.0
                        (fun (x: FDRControl.QValueInput) -> x.Score)
                        (fun (x: FDRControl.QValueInput) -> x.IsDecoy)
                        inputs
                Array.iter
                    (fun pep -> Expect.equal pep 0.0 "every local error estimate is zero")
                    (peps |> Seq.toArray)
                Array.iter
                    (fun q -> Expect.equal q 0.0 "every cumulative error estimate is zero")
                    (qs |> Seq.toArray)

            testCase "binningFunction emits raw unclamped local error estimates as pre-fit input" <| fun _ ->
                // binningFunction is a pre-fit helper: its per-bin estimates are the input points for
                // the logistic fit in calculateQValueLogReg, not final probabilities, so they are raw
                // unclamped ratios. A bin's local estimate is pi0 * 2 * scaledDecoys/binTotal; with
                // global decoy fraction 1/2 the decoy scale is 1/(2*1/2) = 1, so the low bin holding
                // 2 decoys of 2 entries gives 2*2/2 = 2.0 - hand-computed from the estimator, above 1
                // by design - and the decoy-free high bin gives 0.
                let inputs =
                    [
                        FDRControl.createQValueInput 0.5 true
                        FDRControl.createQValueInput 0.9 true
                        FDRControl.createQValueInput 1.5 false
                        FDRControl.createQValueInput 1.9 false
                    ]
                    |> Array.ofList
                let _, peps, _ =
                    FDRControl.binningFunction
                        1.0
                        1.0
                        (fun (x: FDRControl.QValueInput) -> x.Score)
                        (fun (x: FDRControl.QValueInput) -> x.IsDecoy)
                        inputs
                Array.iter2
                    (fun actual expected -> expectFloatClose actual expected "the raw pre-fit local error estimate")
                    (peps |> Seq.toArray)
                    [|2.0; 0.0|]

            testCase "binningFunction emits raw suffix q estimates without monotone correction" <| fun _ ->
                // The q column is the raw cumulative (suffix) decoy/total ratio per bin - pre-fit
                // input for the descending logistic fitted in calculateQValueLogReg, which is
                // monotone by model shape; no monotone or [0,1] correction is applied here. For
                // target/decoy/target in three unit bins the global decoy fraction is 1/3 (decoy
                // scale 1/(2*1/3) = 1.5), giving suffix ratios 2*1.5/3 = 1.0, 2*1.5/2 = 1.5 and
                // 0/1 = 0.0 - hand-computed from the estimator, non-monotone and above 1 by design.
                let inputs =
                    [|
                        FDRControl.createQValueInput 0.5 false
                        FDRControl.createQValueInput 1.5 true
                        FDRControl.createQValueInput 2.5 false
                    |]
                let scores, _, qs =
                    FDRControl.binningFunction
                        1.0
                        1.0
                        (fun (x: FDRControl.QValueInput) -> x.Score)
                        (fun (x: FDRControl.QValueInput) -> x.IsDecoy)
                        inputs
                let qValuesByScore =
                    Array.zip (scores |> Seq.toArray) (qs |> Seq.toArray)
                    |> Array.sortBy fst
                    |> Array.map snd
                Array.iter2
                    (fun actual expected -> expectFloatClose actual expected "the raw pre-fit suffix q estimate")
                    qValuesByScore
                    [|1.0; 1.5; 0.0|]

            testCase "binningFunction returns bins sorted by ascending average score" <| fun _ ->
                let orderedInputs =
                    [
                        FDRControl.createQValueInput 0.1 false
                        FDRControl.createQValueInput 0.9 true
                        FDRControl.createQValueInput 1.1 true
                        FDRControl.createQValueInput 1.9 false
                    ]
                    |> Array.ofList
                let shuffledInputs =
                    [
                        FDRControl.createQValueInput 1.9 false
                        FDRControl.createQValueInput 0.1 false
                        FDRControl.createQValueInput 1.1 true
                        FDRControl.createQValueInput 0.9 true
                    ]
                    |> Array.ofList
                let orderedScores, orderedPeps, orderedQVals =
                    FDRControl.binningFunction
                        1.0
                        1.0
                        (fun (x: FDRControl.QValueInput) -> x.Score)
                        (fun (x: FDRControl.QValueInput) -> x.IsDecoy)
                        orderedInputs
                let scores, peps, qVals =
                    FDRControl.binningFunction
                        1.0
                        1.0
                        (fun (x: FDRControl.QValueInput) -> x.Score)
                        (fun (x: FDRControl.QValueInput) -> x.IsDecoy)
                        shuffledInputs
                Expect.equal scores.Length 2 "the two bandwidth-1 bins produce two score entries"
                Array.iter2
                    (fun actual expected -> expectFloatClose actual expected "the ascending representative bin score")
                    (scores |> Seq.toArray)
                    [|0.5; 1.5|]
                Array.iter2
                    (fun actual expected -> expectFloatClose actual expected "the shuffled bin PEP matches the ordered run")
                    (peps |> Seq.toArray)
                    (orderedPeps |> Seq.toArray)
                Array.iter2
                    (fun actual expected -> expectFloatClose actual expected "the shuffled cumulative q-value matches the ordered run")
                    (qVals |> Seq.toArray)
                    (orderedQVals |> Seq.toArray)
                // bin membership and representative scores are set-properties of the data; input order must not matter.

            testCase "binning accumulates from the queried score upward and scales linearly with pi0" <| fun _ ->
                let inputs =
                    [|
                        FDRControl.createQValueInput 0.1 false
                        FDRControl.createQValueInput 0.2 true
                        FDRControl.createQValueInput 1.1 false
                        FDRControl.createQValueInput 1.2 false
                    |]
                let _, pep1, q1 =
                    FDRControl.binningFunction
                        1.0
                        0.5
                        (fun (x: FDRControl.QValueInput) -> x.Score)
                        (fun (x: FDRControl.QValueInput) -> x.IsDecoy)
                        inputs
                let _, pep2, q2 =
                    FDRControl.binningFunction
                        1.0
                        0.25
                        (fun (x: FDRControl.QValueInput) -> x.Score)
                        (fun (x: FDRControl.QValueInput) -> x.IsDecoy)
                        inputs
                let pep1 = pep1 |> Seq.toArray
                let q1 = q1 |> Seq.toArray
                let pep2 = pep2 |> Seq.toArray
                let q2 = q2 |> Seq.toArray
                Array.iter2 (fun actual expected -> expectFloatClose actual expected "the pi0=0.5 local error") pep1 [|1.0; 0.0|]
                Array.iter2 (fun actual expected -> expectFloatClose actual expected "the pi0=0.5 cumulative error") q1 [|0.5; 0.0|]
                Array.iter2 (fun actual expected -> expectFloatClose actual expected "the local error scales linearly with pi0") pep2 (pep1 |> Array.map (fun value -> value / 2.0))
                Array.iter2 (fun actual expected -> expectFloatClose actual expected "the cumulative error scales linearly with pi0") q2 (q1 |> Array.map (fun value -> value / 2.0))
                // hand-derived from the standard target-decoy estimator with decoy-to-total scaling (decoy fraction 1/4 -> scale 2) and the pi0 prior: the decoy-free high bin must carry ZERO local and cumulative error (pinning suffix-direction accumulation - a prefix accumulator would put the error on the high bin), and both estimates must be LINEAR in pi0 (the halving is convention-independent).
        ]

        testList "Storey" [
            testCase "calculateQValueStorey interpolates hand-counted cumulative decoy/target ratios" <| fun _ ->
                let data =
                    [|
                        FDRControl.createQValueInput 3.0 false
                        FDRControl.createQValueInput 2.0 true
                        FDRControl.createQValueInput 1.0 false
                    |]
                let f =
                    FDRControl.calculateQValueStorey
                        data
                        (fun x -> x.IsDecoy)
                        (fun x -> x.Score)
                        (fun x -> x.Score)
                expectWithin 1e-9 (f 3.0) 0.0 "the q-value at score 3.0"
                expectWithin 1e-9 (f 1.0) 0.5 "the q-value at score 1.0"
                expectWithin 1e-9 (f 2.0) 0.5 "the monotone-corrected q-value at score 2.0"
                expectWithin 1e-9 (f 2.5) 0.25 "linear interpolation between the score-2.0 and score-3.0 knots"
                // 0.25 pins the CONTINUOUS linear-spline interpolation the function advertises by returning a float -> float scorer (a step function would be the other standard choice); the surrounding bound f 2.0 >= f 2.5 >= f 3.0 is the convention-independent part.
                // all knot values hand-counted from the standard cumulative decoy/target estimator plus the non-increasing q-value correction.

            testCase "the decoy and target score selectors are routed by the decoy flag" <| fun _ ->
                let data =
                    [|
                        (3.0, -99.0, false)
                        (-99.0, 2.0, true)
                        (1.0, -99.0, false)
                    |]
                let f =
                    FDRControl.calculateQValueStorey
                        data
                        (fun (_, _, isDecoy) -> isDecoy)
                        (fun (_, decoyScore, _) -> decoyScore)
                        (fun (targetScore, _, _) -> targetScore)
                expectWithin 1e-9 (f 3.0) 0.0 "the target selector supplies the score-3 knot"
                expectWithin 1e-9 (f 2.0) 0.5 "the decoy selector supplies the score-2 knot"
                expectWithin 1e-9 (f 1.0) 0.5 "the target selector supplies the score-1 knot"
                expectWithin 1e-9 (f 2.5) 0.25 "selector-specific score knots interpolate linearly"
                // -99 sentinels in the wrong-side fields make selector misrouting (ignoring isDecoy, or using one accessor for both) produce wildly different knots; matching the hand-counted values proves per-class routing.

            testCase "tied target and decoy scores form one threshold" <| fun _ ->
                let dataA =
                    [|
                        FDRControl.createQValueInput 3.0 false
                        FDRControl.createQValueInput 3.0 true
                        FDRControl.createQValueInput 2.0 false
                    |]
                let dataB = Array.rev dataA
                let getQValues (data: FDRControl.QValueInput[]) =
                    FDRControl.calculateQValueStorey
                        data
                        (fun x -> x.IsDecoy)
                        (fun x -> x.Score)
                        (fun x -> x.Score)
                let fA = getQValues dataA
                let fB = getQValues dataB
                expectWithin 1e-9 (fA 3.0) 0.5 "the tied score-3 threshold has the capped q-value"
                expectWithin 1e-9 (fA 2.0) 0.5 "the score-2 threshold has the hand-counted q-value"
                expectWithin 1e-9 (fB 3.0) 0.5 "reversing tied observations preserves the score-3 q-value"
                expectWithin 1e-9 (fB 2.0) 0.5 "reversing tied observations preserves the score-2 q-value"
                // at threshold 3 the cumulative counts are 1 decoy / 1 target -> raw 1.0, monotone-capped by threshold 2's 1/2 -> 0.5 at both knots; input order of tied observations must not matter - hand-counted.

            ptestCase "Storey q-values are proportions even when decoys outscore all targets" <| fun _ ->
                // PENDING: a q-value estimates a proportion of false discoveries and cannot exceed 1. When
                // decoys outscore every target the implementation substitutes denominator 1 and returns the
                // raw decoy count (probe: 2.0); the monotone pass caps only from lower-score neighbors and
                // does not repair it.
                let data =
                    [|
                        FDRControl.createQValueInput 5.0 true
                        FDRControl.createQValueInput 4.0 true
                        FDRControl.createQValueInput 3.0 false
                    |]
                let f =
                    FDRControl.calculateQValueStorey
                        data
                        (fun x -> x.IsDecoy)
                        (fun x -> x.Score)
                        (fun x -> x.Score)
                Expect.isTrue (f 4.0 <= 1.0) "Storey q-values are proportions"

            testCase "Storey q-values are non-increasing in score and vanish when all decoys score below all targets" <| fun _ ->
                let targets =
                    [|10.0 .. 1.0 .. 19.0|]
                    |> Array.map (fun score -> FDRControl.createQValueInput score false)
                let decoys =
                    [|0.0 .. 1.0 .. 4.0|]
                    |> Array.map (fun score -> FDRControl.createQValueInput score true)
                let f =
                    FDRControl.calculateQValueStorey
                        (Array.append targets decoys)
                        (fun x -> x.IsDecoy)
                        (fun x -> x.Score)
                        (fun x -> x.Score)
                Expect.isTrue
                    (f 10.0 >= f 15.0 && f 15.0 >= f 19.0)
                    "a stricter score threshold cannot have a larger estimated FDR (monotonicity of the corrected estimator)"
                // cumulative decoy/target counting hand-derived; the descent happens where q actually varies.
                expectWithin 1e-9 (f 0.0) 0.5 "the q-value at score 0.0"
                expectWithin 1e-9 (f 2.0) 0.3 "the q-value at score 2.0"
                expectWithin 1e-9 (f 4.0) 0.1 "the q-value at score 4.0"
                Expect.isTrue
                    (f 0.0 > f 2.0 && f 2.0 > f 4.0)
                    "the q-value descends through the decoy region"
                expectWithin 1e-9 (f 19.0) 0.0 "the q-value at score 19.0"
                expectWithin 1e-9 (f 10.0) 0.0 "the q-value at score 10.0"
                // domain reading of a perfectly separated target/decoy score distribution.
        ]

        testList "PEP" [
            testCase "createTargetDecoyHis bins scores with half-bandwidth-centered labels and per-bin decoy counts" <| fun _ ->
                let data = [| (0.2, false); (0.7, true); (0.4, false); (1.3, false) |]
                let bins =
                    FDRControl.createTargetDecoyHis
                        1.0
                        snd
                        (fun (score, _) -> score)
                        (fun (score, _) -> score)
                        data
                    |> Array.sortBy (fun (bin, _, _, _) -> bin)
                Expect.equal bins.Length 2 "the data is split into two score bins"
                let bin0, count0, decoyCount0, median0 = bins.[0]
                expectFloatClose bin0 0.5 "the non-negative bin label is centered by half the bandwidth"
                Expect.equal count0 3 "the first bin contains three entries"
                Expect.equal decoyCount0 1 "the first bin contains one decoy"
                expectFloatClose median0 0.4 "the first bin median score"
                let bin1, count1, decoyCount1, median1 = bins.[1]
                expectFloatClose bin1 1.5 "the second non-negative bin label is centered by half the bandwidth"
                Expect.equal count1 1 "the second bin contains one entry"
                Expect.equal decoyCount1 0 "the second bin contains no decoys"
                expectFloatClose median1 1.3 "the second bin median score"
                let negativeBin, _, _, _ =
                    FDRControl.createTargetDecoyHis
                        1.0
                        snd
                        (fun (score, _) -> score)
                        (fun (score, _) -> score)
                        [| (-0.3, true) |]
                    |> Array.exactlyOne
                expectFloatClose negativeBin -0.5 "the negative bin label uses the negative half-bandwidth branch"

            testCase "calculatePEPValues returns score-sorted decoy/total ratios" <| fun _ ->
                let dataFreq = [| (2.0, 4.0, 1.0); (1.0, 2.0, 1.0) |]
                let actual =
                    FDRControl.calculatePEPValues
                        (fun (_, total, _) -> total)
                        (fun (_, _, decoy) -> decoy)
                        (fun (score, _, _) -> score)
                        dataFreq
                Expect.equal actual [ (1.0, 0.5); (2.0, 0.25) ] "PEP values are sorted ascending by score"

            testCase "logitTransformPepValues drops endpoint pep values and log10-logit-transforms the rest" <| fun _ ->
                let scores, pepValues =
                    FDRControl.logitTransformPepValues
                        [| 1.0; 2.0; 3.0; 4.0 |]
                        [| 0.0; 0.5; 0.9; 1.0 |]
                Expect.equal scores [| 2.0; 3.0 |] "endpoint PEP values are removed"
                expectFloatClose pepValues.[0] 0.0 "the PEP 0.5 logit is zero"
                expectFloatClose pepValues.[1] (log10 9.0) "the PEP 0.9 logit is log10(9)"

            testCase "initCalculateLin yields a monotonically usable PEP mapping on a separable target/decoy set" <| fun _ ->
                let targets = [| for score in 1.0 .. 0.25 .. 8.0 -> (score, false) |]
                let decoys =
                    Array.append
                        [| for score in -8.0 .. 0.25 .. -1.0 -> (score, true) |]
                        [| (1.1, true); (1.2, true); (2.1, true); (3.1, true) |]
                let data = Array.append targets decoys
                let msgs = ResizeArray<string>()
                let f =
                    FDRControl.initCalculateLin
                        msgs.Add
                        0.5
                        snd
                        (fun (score, _) -> score)
                        (fun (score, _) -> score)
                        data
                let atOne = f 1.0
                let atEight = f 8.0
                let atFive = f 5.0
                Expect.isTrue (not (Double.IsNaN atOne) && not (Double.IsInfinity atOne)) "the PEP at score 1.0 is finite"
                Expect.isTrue (not (Double.IsNaN atEight) && not (Double.IsInfinity atEight)) "the PEP at score 8.0 is finite"
                Expect.isTrue (atEight <= atOne) "the PEP does not increase for the better score"
                Expect.isTrue (0.0 <= atFive && atFive <= 1.0) "the PEP at score 5.0 is a probability"
                Expect.equal msgs.Count 3 "the initializer emits the two setup traces and chosen bandwidth trace"

            testCase "getLogisticRegressionFunction fits a descending logistic mapping" <| fun _ ->
                let f =
                    FDRControl.getLogisticRegressionFunction
                        (vector [| 1.0; 2.0; 3.0; 4.0; 5.0; 6.0 |])
                        (vector [| 1.0; 1.0; 1.0; 0.0; 0.0; 0.0 |])
                        0.0001
                let atMiddle = f 3.5
                Expect.isTrue (0.0 <= atMiddle && atMiddle <= 1.0) "the logistic prediction is a probability"
                Expect.isTrue (f 1.0 > f 6.0) "the logistic prediction is higher on the y=1 side"
        ]

        testList "LogisticRegression" [
            testCase "calculateQValueLogReg yields a finite, broadly descending q-value function on well-separated data" <| fun _ ->
                let targets =
                    [|0 .. 29|]
                    |> Array.map (fun i -> FDRControl.createQValueInput (10.0 + 0.5 * float i) false)
                let decoys =
                    [|0 .. 29|]
                    |> Array.map (fun i -> FDRControl.createQValueInput (0.0 + 0.2 * float i) true)
                let f =
                    FDRControl.calculateQValueLogReg
                        1.0
                        1.0
                        (Array.append targets decoys)
                        (fun x -> x.IsDecoy)
                        (fun x -> x.Score)
                        (fun x -> x.Score)
                let qAt1 = f 1.0
                let qAt8 = f 8.0
                let qAt20 = f 20.0
                Expect.isTrue
                    (not (Double.IsNaN qAt1) && not (Double.IsInfinity qAt1))
                    "the q-value at score 1.0 is finite"
                Expect.isTrue
                    (not (Double.IsNaN qAt8) && not (Double.IsInfinity qAt8))
                    "the q-value at score 8.0 is finite"
                Expect.isTrue
                    (not (Double.IsNaN qAt20) && not (Double.IsInfinity qAt20))
                    "the q-value at score 20.0 is finite"
                Expect.isTrue
                    (qAt1 >= qAt20)
                    "the fitted logistic descends from the decoy-dominated region to the target-dominated region"
                let sampledValues = [|1.0; 4.0; 8.0; 12.0; 16.0; 20.0|] |> Array.map (fun score -> score, f score)
                sampledValues
                |> Array.iter (fun (score, value) ->
                    Expect.isTrue
                        (value >= 0.0 && value <= 1.0)
                        (sprintf "the fitted q-value at score %g is a probability; observed %g" score value))
                Expect.isTrue (f 1.0 > f 20.0) "the fitted logistic is strictly descending between low and high scores"
                // only coarse, convergence-robust properties are asserted; exact fitted values depend on external Levenberg-Marquardt behavior and would be golden-value laundering.
        ]
    ]
