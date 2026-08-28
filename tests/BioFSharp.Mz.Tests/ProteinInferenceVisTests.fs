module ProteinInferenceVisTests

open System
open System.IO
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
    ]
