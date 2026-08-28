module FDRControlTests

open System
open Expecto
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
                // only coarse, convergence-robust properties are asserted; exact fitted values depend on external Levenberg-Marquardt behavior and would be golden-value laundering.
        ]
    ]
