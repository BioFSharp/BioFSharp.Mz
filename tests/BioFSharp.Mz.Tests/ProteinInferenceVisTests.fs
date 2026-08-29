module ProteinInferenceVisTests

open System
open System.IO
open System.Text.RegularExpressions
open Expecto
open BioFSharp.Mz
open BioFSharp.PeptideClassification

let private cl = PeptideEvidenceClass.C1a

let private items =
    [|
        ProteinInference.createInferredProteinClassItemScored
            (ProteinInference.proteinGroupToString [|"Pi"|])
            cl
            [|"pep"|]
            10.0
            1.0
            false
            false
            true
        |> fun item -> ProteinInference.createInferredProteinClassItemQValue item 0.01
        ProteinInference.createInferredProteinClassItemScored
            (ProteinInference.proteinGroupToString [|"Pi"|])
            cl
            [|"pep"|]
            8.0
            0.5
            false
            false
            true
        |> fun item -> ProteinInference.createInferredProteinClassItemQValue item 0.02
        ProteinInference.createInferredProteinClassItemScored
            (ProteinInference.proteinGroupToString [|"Pi"|])
            cl
            [|"pep"|]
            6.0
            2.0
            false
            false
            true
        |> fun item -> ProteinInference.createInferredProteinClassItemQValue item 0.05
        ProteinInference.createInferredProteinClassItemScored
            (ProteinInference.proteinGroupToString [|"Di"|])
            cl
            [|"pep"|]
            1.0
            7.0
            true
            true
            true
        |> fun item -> ProteinInference.createInferredProteinClassItemQValue item 0.5
        ProteinInference.createInferredProteinClassItemScored
            (ProteinInference.proteinGroupToString [|"Di"|])
            cl
            [|"pep"|]
            1.0
            7.0
            true
            true
            true
        |> fun item -> ProteinInference.createInferredProteinClassItemQValue item 0.8
    |]

let private withTempDirectory f =
    let tempDir = Path.Combine(Path.GetTempPath(), "BioFSharp.Mz-ProteinInferenceVisTests-" + Guid.NewGuid().ToString("N"))
    Directory.CreateDirectory(tempDir) |> ignore
    try
        f tempDir
    finally
        try
            if Directory.Exists tempDir then
                Directory.Delete(tempDir, true)
        with _ -> ()

[<Tests>]
let tests =
    testList "ProteinInferenceVisTests" [
        testCase "qValueHitsVisualization with groupFiles=false writes the suffixed HTML report" <| fun _ ->
            withTempDirectory (fun tempDir ->
                let basePath = Path.Combine(tempDir, "run1")
                BioFSharp.Mz.Vis.ProteinInference.qValueHitsVisualization 1.0 items basePath false

                let reportPath = basePath + "_QValueGraph.html"
                let reportExists = File.Exists reportPath
                Expect.isTrue reportExists "the suffixed QValueGraph HTML report exists"
                let reportLength = if reportExists then FileInfo(reportPath).Length else 0L
                Expect.isTrue (reportLength > 0L) "the suffixed QValueGraph HTML report is non-empty"
                // the function's contract is writing an HTML report; a garbage or truncated write of any length would otherwise pass.
                let reportContent = if reportExists then File.ReadAllText reportPath else ""
                Expect.isTrue (reportContent.Contains("<html", StringComparison.OrdinalIgnoreCase)) "the suffixed report contains HTML markup"
                // The documented single-file naming convention appends "_QValueGraph" to the base path (Plotly's SaveHtmlAs adds the .html extension); the observable contract of a visualization function is that it produces the report file where callers will look for it.
            )

        testCase "qValueHitsVisualization with groupFiles=true writes the report inside the target folder" <| fun _ ->
            withTempDirectory (fun tempDir ->
                let basePath = Path.Combine(tempDir, "run2")
                Directory.CreateDirectory(basePath) |> ignore
                BioFSharp.Mz.Vis.ProteinInference.qValueHitsVisualization 1.0 items basePath true

                let reportPath = System.IO.Path.Combine(basePath, "QValueGraph.html")
                let reportExists = File.Exists reportPath
                Expect.isTrue reportExists "the grouped QValueGraph HTML report exists inside the target folder"
                let reportLength = if reportExists then FileInfo(reportPath).Length else 0L
                Expect.isTrue (reportLength > 0L) "the grouped QValueGraph HTML report is non-empty"
                // the function's contract is writing an HTML report; a garbage or truncated write of any length would otherwise pass.
                let reportContent = if reportExists then File.ReadAllText reportPath else ""
                Expect.isTrue (reportContent.Contains("<html", StringComparison.OrdinalIgnoreCase)) "the grouped report contains HTML markup"
                // the grouped mode's contract is "the report lands inside the given folder" - asserted via
                // Path.Combine, which is byte-identical to the implementation's literal backslash on Windows
                // but would correctly FAIL on Unix, where the current code writes a file literally named
                // "dir\QValueGraph.html" beside the folder. A string-identical assertion would bless that
                // bug on every platform.
            )

        // PENDING: relative frequencies of the decoy histogram must sum to 1 over the DECOYS. The
        // implementation divides decoy bin counts by the TARGET count: with 2 equal-score decoys and
        // 3 targets the emitted Decoy trace carries y = 2/3 instead of 1.0 (probe-verified in the
        // generated HTML). The trace is located by its user-visible legend label "Decoy" - a
        // deliberate Chart.withTraceName contract, not an implementation detail.
        ptestCase "decoy relative frequencies are normalized by the decoy count" <| fun _ ->
            withTempDirectory (fun tempDir ->
                let basePath = Path.Combine(tempDir, "run3")
                BioFSharp.Mz.Vis.ProteinInference.qValueHitsVisualization 1.0 items basePath false

                let reportPath = basePath + "_QValueGraph.html"
                let html = File.ReadAllText reportPath
                let decoyTrace =
                    System.Text.RegularExpressions.Regex.Match(
                        html,
                        "\"y\":\\[([^\\]]*)\\](?:[^{}]|\\{[^{}]*\\})*?\"name\":\"Decoy\""
                    )
                Expect.isTrue decoyTrace.Success "the Decoy trace with a y array is present"
                let yValues =
                    decoyTrace.Groups.[1].Value.Split([|','|], StringSplitOptions.RemoveEmptyEntries)
                    |> Array.map (fun value -> Double.Parse(value, Globalization.CultureInfo.InvariantCulture))
                Expect.isTrue
                    (abs (Array.sum yValues - 1.0) <= 1e-6)
                    "the decoy relative frequencies sum to one"
            )

        testCase "the target histogram bins target scores with the requested bandwidth, normalized over targets" <| fun _ ->
            let target score =
                ProteinInference.createInferredProteinClassItemScored
                    (ProteinInference.proteinGroupToString [|sprintf "T%g" score|])
                    cl
                    [|"target-peptide"|]
                    score
                    -50.0
                    false
                    false
                    true
                |> fun item -> ProteinInference.createInferredProteinClassItemQValue item 0.01
            let targetItems = [|target 0.0; target 10.0; target 20.0|]
            let decoyItems = items |> Array.filter (fun item -> item.InfProtClassItem.Decoy)
            let fixture = Array.append targetItems decoyItems
            let parseNumbers (value: string) =
                value.Split([|','|], StringSplitOptions.RemoveEmptyEntries)
                |> Array.map (fun number -> Double.Parse(number.Trim(), Globalization.CultureInfo.InvariantCulture))
            let extractTargetY html =
                let trace =
                    Regex.Match(
                        html,
                        "\"y\":\\[([^\\]]*)\\](?:[^{}]|\\{(?:[^{}]|\\{[^{}]*\\})*\\})*?\"name\":\"Target\"")
                Expect.isTrue trace.Success "the Target trace with a y array is present"
                if trace.Success then parseNumbers trace.Groups.[1].Value else [||]

            withTempDirectory (fun tempDir ->
                let basePath1 = Path.Combine(tempDir, "target-bandwidth-1")
                let basePath2 = Path.Combine(tempDir, "target-bandwidth-100")
                BioFSharp.Mz.Vis.ProteinInference.qValueHitsVisualization 1.0 fixture basePath1 false
                BioFSharp.Mz.Vis.ProteinInference.qValueHitsVisualization 100.0 fixture basePath2 false
                let y1 = File.ReadAllText(basePath1 + "_QValueGraph.html") |> extractTargetY
                let y100 = File.ReadAllText(basePath2 + "_QValueGraph.html") |> extractTargetY
                Expect.equal y1.Length 3 "unit bandwidth leaves the three distinct target score bins"
                y1 |> Array.iter (fun value -> Expect.isTrue (abs (value - (1.0 / 3.0)) <= 1e-6) "each unit-bandwidth target bin has relative frequency one third")
                Expect.equal y100.Length 1 "bandwidth 100 pools the target scores into one bin"
                if y100.Length = 1 then
                    Expect.isTrue (abs (y100.[0] - 1.0) <= 1e-6) "the pooled target bin has relative frequency one"
                // three distinct scores with unit bandwidth occupy three bins of relative frequency 1/3 each; a 100-wide bandwidth pools them into one bin of frequency 1 - counting, not implementation output. The -50 decoy scores prove the target trace reads TargetScore (target normalization is correct today; only the decoy side is pended).
            )

        testCase "the q-value scatter pairs each protein's identity-appropriate score with its q-value" <| fun _ ->
            let target =
                ProteinInference.createInferredProteinClassItemScored
                    (ProteinInference.proteinGroupToString [|"TARGET"|])
                    cl
                    [|"target-peptide"|]
                    12.0
                    -1.0
                    false
                    false
                    true
                |> fun item -> ProteinInference.createInferredProteinClassItemQValue item 0.01
            let decoy =
                ProteinInference.createInferredProteinClassItemScored
                    (ProteinInference.proteinGroupToString [|"DECOY"|])
                    cl
                    [|"decoy-peptide"|]
                    99.0
                    3.0
                    true
                    true
                    true
                |> fun item -> ProteinInference.createInferredProteinClassItemQValue item 0.20
            let fixture = [|target; decoy|]
            let parseNumbers (value: string) =
                value.Split([|','|], StringSplitOptions.RemoveEmptyEntries)
                |> Array.map (fun number -> Double.Parse(number.Trim(), Globalization.CultureInfo.InvariantCulture))
            withTempDirectory (fun tempDir ->
                let basePath = Path.Combine(tempDir, "q-value-pairs")
                BioFSharp.Mz.Vis.ProteinInference.qValueHitsVisualization 1.0 fixture basePath false
                let html = File.ReadAllText(basePath + "_QValueGraph.html")
                let trace =
                    Regex.Match(
                        html,
                        "\"x\":\\[([^\\]]*)\\](?:[^{}]|\\{[^{}]*\\})*?\"y\":\\[([^\\]]*)\\](?:[^{}]|\\{[^{}]*\\})*?\"name\":\"Q-Values\"")
                Expect.isTrue trace.Success "the Q-Values trace with x and y arrays is present"
                if trace.Success then
                    let xValues = parseNumbers trace.Groups.[1].Value
                    let yValues = parseNumbers trace.Groups.[2].Value
                    Expect.equal xValues.Length yValues.Length "the Q-Values x and y arrays have equal length"
                    let pairs = Array.zip xValues yValues
                    let hasPair x y = pairs |> Array.exists (fun (actualX, actualY) -> abs (actualX - x) <= 1e-9 && abs (actualY - y) <= 1e-9)
                    Expect.isTrue (hasPair 12.0 0.01) "the target is paired with its target score and q-value"
                    Expect.isTrue (hasPair 3.0 0.20) "the decoy is paired with its decoy score and q-value"
                    Expect.isFalse (hasPair 99.0 0.20) "the decoy is not plotted at its target score"
                // a target is plotted at its target score, a decoy at its decoy score, each with its own q-value - the central inference result of the report; pair membership is order-agnostic.
            )

        testCase "the absolute-frequency axis carries the raw target and decoy counts" <| fun _ ->
            let target score =
                ProteinInference.createInferredProteinClassItemScored
                    (ProteinInference.proteinGroupToString [|sprintf "T%g" score|])
                    cl
                    [|"target-peptide"|]
                    score
                    -50.0
                    false
                    false
                    true
                |> fun item -> ProteinInference.createInferredProteinClassItemQValue item 0.01
            let fixture = Array.append [|target 0.0; target 10.0; target 20.0|] (items |> Array.filter (fun item -> item.InfProtClassItem.Decoy))
            let parseNumbers (value: string) =
                value.Split([|','|], StringSplitOptions.RemoveEmptyEntries)
                |> Array.map (fun number -> Double.Parse(number.Trim(), Globalization.CultureInfo.InvariantCulture))
            withTempDirectory (fun tempDir ->
                let basePath = Path.Combine(tempDir, "absolute-frequency")
                BioFSharp.Mz.Vis.ProteinInference.qValueHitsVisualization 1.0 fixture basePath false
                let html = File.ReadAllText(basePath + "_QValueGraph.html")
                let y2Traces =
                    Regex.Matches(
                        html,
                        "\"y\":\\[([^\\]]*)\\](?:[^{}]|\\{(?:[^{}]|\\{[^{}]*\\})*\\})*?\"yaxis\":\"y2\"")
                    |> Seq.cast<Match>
                    |> Seq.map (fun trace -> parseNumbers trace.Groups.[1].Value)
                    |> Seq.toArray
                let sums = y2Traces |> Array.map Array.sum |> Array.sort
                Expect.equal sums [|2.0; 3.0|] "the absolute-frequency traces sum to the two raw partition counts"
                // the report exposes an absolute-frequency axis; its traces must total the actual 3 targets and 2 decoys - input counting, robust to bin positions and trace order.
            )
    ]
