module SearchEngineResultTests

open Expecto
open BioFSharp.Mz

let mkResultWith sequence score =
    SearchEngineResult.createSearchEngineResult SearchEngineResult.SearchEngine.SEQUESTLike
        "spec1" 1 1 0 true 10.0 sequence 2 500.0 998.0 998.0 7 score 9.9 9.9

let mkResult score =
    SearchEngineResult.createSearchEngineResult SearchEngineResult.SearchEngine.SEQUESTLike
        "spec1" 1 1 0 true 10.0 "PEPTIDE" 2 500.0 998.0 998.0 7 score 0.0 0.0

let expectFloatListClose accuracy (actual: float list) (expected: float list) message =
    Expect.equal actual.Length expected.Length (sprintf "%s length" message)
    List.iter2
        (fun actualValue expectedValue -> Expect.floatClose accuracy actualValue expectedValue message)
        actual
        expected

[<Tests>]
let tests =
    testList "SearchEngineResultTests" [
        testCase "calcNormDeltaBestToRest computes the SEQUEST deltaCn against the top hit" <| fun _ ->
            let results = [mkResult 100.0; mkResult 80.0; mkResult 50.0]
            let r =
                SearchEngineResult.calcNormDeltaBestToRest results
                |> List.map (fun x -> x.NormDeltaBestToRest)
            expectFloatListClose Accuracy.high r [0.0; 0.2; 0.5] "deltaCn against the top hit"
            // Published SEQUEST deltaCn definition: (XCorr_top - XCorr_i) / XCorr_top ->
            // (100-100)/100 = 0, (100-80)/100 = 0.2, (100-50)/100 = 0.5 (hand-computed).
            // The top hit's delta is 0 by definition.

        testCase "calcNormDeltaNext computes adjacent score gaps normalized by the best score" <| fun _ ->
            let results = [mkResult 100.0; mkResult 80.0; mkResult 50.0]
            let r =
                SearchEngineResult.calcNormDeltaNext results
                |> List.map (fun x -> x.NormDeltaNext)
            expectFloatListClose Accuracy.high r [0.2; 0.3; 0.0] "adjacent score gaps"
            // The adjacent-gap formula is documented in the source: (score_i - score_{i+1}) /
            // bestScore -> (100-80)/100 and (80-50)/100 (hand-computed). The trailing 0 for the
            // last PSM is NOT documented; the last PSM has no successor.

        testCase "adjacent deltas telescope to the best-to-rest deltas" <| fun _ ->
            let results = [mkResult 90.0; mkResult 60.0; mkResult 45.0; mkResult 30.0]
            let bestToRest =
                SearchEngineResult.calcNormDeltaBestToRest results
                |> List.map (fun x -> x.NormDeltaBestToRest)
            let next =
                SearchEngineResult.calcNormDeltaNext results
                |> List.map (fun x -> x.NormDeltaNext)

            Expect.equal bestToRest.Length next.Length "both delta lists preserve the result count"
            bestToRest
            |> List.iteri (fun i actual ->
                let expected = next |> List.take i |> List.sum
                Expect.floatClose
                    Accuracy.high
                    actual
                    expected
                    (sprintf "best-to-rest delta at index %d is the sum of preceding adjacent deltas" i))
            // The gap to the top hit is the telescoping sum of adjacent gaps, all normalized by
            // the same best score. This cross-checks the two functions against each other
            // without any implementation-derived number.

        testCase "both normalizations return the documented sentinel results when the best score is zero" <| fun _ ->
            let zeros = [mkResult 0.0; mkResult 0.0]

            Expect.equal
                (SearchEngineResult.calcNormDeltaBestToRest zeros |> List.map (fun x -> x.NormDeltaBestToRest))
                [1.0; 1.0]
                "zero best score gives a best-to-rest sentinel of 1.0"
            Expect.equal
                (SearchEngineResult.calcNormDeltaNext zeros |> List.map (fun x -> x.NormDeltaNext))
                [0.0; 0.0]
                "zero best score gives an adjacent-delta sentinel of 0.0"
            // Documented sentinel behavior for a zero top score (the formulas would divide by
            // zero); the doc comments on both functions state these values for the
            // "best Score equals 0" case.

        testCase "both normalizations map the empty list to the empty list" <| fun _ ->
            Expect.equal
                (SearchEngineResult.calcNormDeltaBestToRest [])
                ([] : SearchEngineResult.SearchEngineResult<float> list)
                "no PSMs produce no best-to-rest deltas"
            Expect.equal
                (SearchEngineResult.calcNormDeltaNext [])
                ([] : SearchEngineResult.SearchEngineResult<float> list)
                "no PSMs produce no adjacent deltas"
            // No PSMs, no deltas - boundary contract.

        testCase "a single PSM has zero deltas" <| fun _ ->
            let single = [mkResult 42.0]
            expectFloatListClose
                Accuracy.high
                (SearchEngineResult.calcNormDeltaBestToRest single |> List.map (fun x -> x.NormDeltaBestToRest))
                [0.0]
                "the sole PSM is the top hit"
            expectFloatListClose
                Accuracy.high
                (SearchEngineResult.calcNormDeltaNext single |> List.map (fun x -> x.NormDeltaNext))
                [0.0]
                "the sole PSM has no successor"

        testCase "normalization preserves every non-delta field" <| fun _ ->
            let results =
                [ mkResultWith "PEPA" 100.0
                  mkResultWith "PEPB" 80.0
                  mkResultWith "PEPC" 50.0 ]
            let r : SearchEngineResult.SearchEngineResult<float> list =
                SearchEngineResult.calcNormDeltaBestToRest results

            Expect.equal
                (r |> List.map (fun x -> x.StringSequence, x.Score))
                [("PEPA", 100.0); ("PEPB", 80.0); ("PEPC", 50.0)]
                "each PSM keeps its own sequence paired with its own score, in order"
            r
            |> List.iter (fun x ->
                Expect.equal x.SpectrumID "spec1" "spectrum ID passes through unchanged"
                Expect.equal x.IsTarget true "target flag passes through unchanged"
                Expect.equal x.PrecursorCharge 2 "precursor charge passes through unchanged")
            // The functions' contract is to fill exactly one delta field; PSM identity and scores
            // must pass through untouched - otherwise downstream FDR/inference would silently work
            // on corrupted PSMs. Distinct per-element sequences catch cross-element field mixing
            // that identical fixtures would hide.

        testCase "normalization changes exactly one field of every record" <| fun _ ->
            let inputs =
                [ SearchEngineResult.createSearchEngineResult SearchEngineResult.SearchEngine.SEQUESTLike "s1" 101 201 0 true 10.0 "PA" 2 401.0 801.0 801.1 2 100.0 9.9 9.9
                  SearchEngineResult.createSearchEngineResult SearchEngineResult.SearchEngine.AndromedaLike "s2" 102 202 1 false 20.0 "PB" 3 502.0 902.0 902.1 3 80.0 9.9 9.9
                  SearchEngineResult.createSearchEngineResult SearchEngineResult.SearchEngine.XTandemLike "s3" 103 203 2 true 30.0 "PC" 4 603.0 1003.0 1003.1 4 50.0 9.9 9.9 ]
            let bestToRest = SearchEngineResult.calcNormDeltaBestToRest inputs
            let next = SearchEngineResult.calcNormDeltaNext inputs
            let expectedBestToRest =
                [ { inputs.[0] with NormDeltaBestToRest = 0.0 }
                  { inputs.[1] with NormDeltaBestToRest = 0.2 }
                  { inputs.[2] with NormDeltaBestToRest = 0.5 } ]
            let expectedNext =
                [ { inputs.[0] with NormDeltaNext = 0.2 }
                  { inputs.[1] with NormDeltaNext = 0.3 }
                  { inputs.[2] with NormDeltaNext = 0.0 } ]
            Expect.equal bestToRest expectedBestToRest "best-to-rest normalization changes only its designated field"
            Expect.equal next expectedNext "adjacent normalization changes only its designated field"
            // whole-record equality against the input-with-one-field-changed pins EVERY field at once - IDs, masses, scan time, engine, the sibling delta - without enumerating them; corruption of any metadata fails.

        testCase "a zero best score triggers the sentinels even over negative tails" <| fun _ ->
            let results = [mkResult 0.0; mkResult -2.0; mkResult -5.0]
            let bestToRest =
                SearchEngineResult.calcNormDeltaBestToRest results
                |> List.map (fun result -> result.NormDeltaBestToRest)
            let next =
                SearchEngineResult.calcNormDeltaNext results
                |> List.map (fun result -> result.NormDeltaNext)
            expectFloatListClose Accuracy.high bestToRest [1.0; 1.0; 1.0] "zero best score gives best-to-rest sentinels"
            expectFloatListClose Accuracy.high next [0.0; 0.0; 0.0] "zero best score gives adjacent-delta sentinels"
            // the documented sentinel depends only on the BEST score being zero; mixed tails distinguish the documented policy from an all-zero special case or per-element sentinels.

        testCase "negative lower scores produce deltas above one, unclamped" <| fun _ ->
            let results = [mkResult 10.0; mkResult -5.0; mkResult -25.0]
            let bestToRest =
                SearchEngineResult.calcNormDeltaBestToRest results
                |> List.map (fun result -> result.NormDeltaBestToRest)
            let next =
                SearchEngineResult.calcNormDeltaNext results
                |> List.map (fun result -> result.NormDeltaNext)
            expectFloatListClose Accuracy.high bestToRest [0.0; 1.5; 3.5] "negative scores yield unclamped best-to-rest deltas"
            expectFloatListClose Accuracy.high next [1.5; 2.0; 0.0] "negative scores yield unclamped adjacent deltas"
            // hand-computed from the documented formulas: (10-s)/10 and adjacent gaps /10. A [0,1] clamp would silently distort these ranking features; this pins their absence for legitimate inputs.

        testCase "the production composition fills both delta fields without clobbering either" <| fun _ ->
            let results =
                [ mkResultWith "PEPA" 100.0
                  mkResultWith "PEPB" 80.0
                  mkResultWith "PEPC" 50.0 ]
            let r =
                results
                |> SearchEngineResult.calcNormDeltaBestToRest
                |> SearchEngineResult.calcNormDeltaNext

            expectFloatListClose
                Accuracy.high
                (r |> List.map (fun x -> x.NormDeltaBestToRest))
                [0.0; 0.2; 0.5]
                "best-to-rest deltas survive the second stage"
            expectFloatListClose
                Accuracy.high
                (r |> List.map (fun x -> x.NormDeltaNext))
                [0.2; 0.3; 0.0]
                "adjacent deltas are computed by the second stage"
            Expect.equal
                (r |> List.map (fun x -> x.StringSequence, x.Score))
                [("PEPA", 100.0); ("PEPB", 80.0); ("PEPC", 50.0)]
                "each PSM keeps its own sequence paired with its own score"
            // Every production caller (SequestLike, AndromedaLike, XScoring) pipes bestToRest into
            // next; the load-bearing contract is that the second normalization does not clobber
            // the first field and that record rebuilding keeps each PSM's identity with its score.
    ]
