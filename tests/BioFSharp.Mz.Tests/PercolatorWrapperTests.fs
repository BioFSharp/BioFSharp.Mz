module PercolatorWrapperTests

open System.IO
open Expecto
open BioFSharp.Mz

[<Tests>]
let tests =
    testList "PercolatorWrapperTests" [
        testCase "fileInfoToWindowsPath returns the full native path" <| fun _ ->
            let fi = FileInfo(@"C:\Data\run.pin")
            Expect.equal
                (PercolatorWrapper.Parameters.fileInfoToWindowsPath fi)
                @"C:\Data\run.pin"
                "FileInfo.FullName is the native full path"
            // the input path is already canonical, so Windows conversion must be the identity on it - asserted against the literal, not against FileInfo.FullName (which would compare the implementation with itself). Windows-only fixture (drive paths resolve differently on Linux); tests do not run on the ubuntu docs workflow.

        testCase "fileInfoToLinuxPath maps a drive-letter path to its WSL mount" <| fun _ ->
            Expect.equal
                (PercolatorWrapper.Parameters.fileInfoToLinuxPath (FileInfo(@"C:\Data\run.pin")))
                "/mnt/c/Data/run.pin"
                "Windows drive-letter paths map to the documented WSL mount"
            // WSL mounts Windows drives at /mnt/<lowercase drive letter> with forward slashes - the documented WSL interop convention, independent of this codebase.

        testCase "stringOf serializes general options to their percolator CLI flags" <| fun _ ->
            let conv = fun (_: System.IO.FileInfo) -> "PATH"
            let options : PercolatorWrapper.Parameters.GeneralOptions list = [
                PercolatorWrapper.Parameters.Help
                PercolatorWrapper.Parameters.VerbosityOfOutput 3
            ]
            let s =
                PercolatorWrapper.Parameters.stringOf
                    conv
                    (PercolatorWrapper.Parameters.PercolatorParams.GeneralOptions options)
            Expect.stringContains s "--help" "the help flag is serialized"
            Expect.stringContains s "--verbose 3" "the verbosity flag and value are serialized"
            // --help and --verbose <n> are percolator's published command-line flags; the serializer's contract is emitting text the percolator binary accepts. Fragment containment (not exact whitespace) is asserted so incidental spacing is not pinned.

        testCase "stringOf routes file arguments through the supplied path converter" <| fun _ ->
            let conv = fun (_: System.IO.FileInfo) -> "CONVERTED"
            let rawPath = @"C:\x\y.tab"
            let output : PercolatorWrapper.Parameters.FileOutputOptions =
                PercolatorWrapper.Parameters.POUTXML (FileInfo rawPath)
            let s =
                PercolatorWrapper.Parameters.stringOf
                    conv
                    (PercolatorWrapper.Parameters.PercolatorParams.FileOutputOptions [output])
            Expect.stringContains s "CONVERTED" "the supplied path converter output is serialized"
            Expect.isFalse (s.Contains rawPath) "the raw file path is not serialized"
            // The wrapper must let the caller decide path dialect (Windows vs WSL); leaking the raw path would break cross-OS invocation. Observed via the stub converter - no filesystem or OS dependence.

        // PENDING: real proteomics paths commonly contain spaces ("C:\Proteomics Runs\..."), and the
        // serializer emits converted paths UNQUOTED, so standard Windows argument parsing splits them
        // into multiple operands - percolator receives a broken file reference. Argued-correct: the
        // serialized fragment must carry the full path as one operand (quoted).
        ptestCase "converted paths containing spaces remain single command-line operands" <| fun _ ->
            let s =
                PercolatorWrapper.Parameters.stringOf
                    PercolatorWrapper.Parameters.fileInfoToWindowsPath
                    (PercolatorWrapper.Parameters.PercolatorParams.FileOutputOptions [
                        PercolatorWrapper.Parameters.POUTTAB_PSMs (FileInfo @"C:\Proteomics Runs\run.psms")
                    ])
            Expect.stringContains s "\"C:\Proteomics Runs\run.psms\"" "a converted path containing spaces is quoted as one operand"

        testCase "each remaining option family serializes its published percolator flags" <| fun _ ->
            let conv = fun (_: System.IO.FileInfo) -> "CONVERTED"

            let fileInputValidation =
                PercolatorWrapper.Parameters.stringOf
                    conv
                    (PercolatorWrapper.Parameters.PercolatorParams.FileInputOptions [
                        PercolatorWrapper.Parameters.SkipSchemeValidation
                    ])
            Expect.stringContains fileInputValidation "--no-schema-validation" "schema validation can be disabled"

            let deprecatedXmlInput =
                PercolatorWrapper.Parameters.stringOf
                    conv
                    (PercolatorWrapper.Parameters.PercolatorParams.FileInputOptions [
                        PercolatorWrapper.Parameters.DeprecatedPINXML (FileInfo @"C:\x\in.xml")
                    ])
            Expect.stringContains deprecatedXmlInput "--xml-in" "deprecated XML input uses the XML input flag"
            Expect.stringContains deprecatedXmlInput "CONVERTED" "deprecated XML input routes its path through the converter"
            Expect.isFalse (deprecatedXmlInput.Contains @"C:\x\in.xml") "deprecated XML input does not serialize its raw path"

            let svmWeights =
                PercolatorWrapper.Parameters.stringOf
                    conv
                    (PercolatorWrapper.Parameters.PercolatorParams.SVMFeatureOptions [
                        PercolatorWrapper.Parameters.OUT_SVMWeights (FileInfo @"C:\x\w.txt")
                    ])
            Expect.stringContains svmWeights "--weights" "SVM weight output uses the weights flag"
            Expect.stringContains svmWeights "CONVERTED" "SVM weight output routes its path through the converter"

            let proteinInference =
                PercolatorWrapper.Parameters.stringOf
                    conv
                    (PercolatorWrapper.Parameters.PercolatorParams.ProteinInferenceOptions_Percolator [
                        PercolatorWrapper.Parameters.Fasta (FileInfo @"C:\x\db.fasta")
                    ])
            Expect.stringContains proteinInference "--picked-protein" "protein inference uses the picked-protein flag"
            Expect.stringContains proteinInference "CONVERTED" "protein inference routes its fasta path through the converter"

            let fido =
                PercolatorWrapper.Parameters.stringOf
                    conv
                    (PercolatorWrapper.Parameters.PercolatorParams.ProteinInferenceOptions_FIDO [
                        PercolatorWrapper.Parameters.Alpha 0.25
                    ])
            Expect.stringContains fido "--fido-alpha 0.25" "Fido alpha uses its flag and value"
            // one representative per option family: the serializer's contract is producing percolator's published long options with converter-routed paths; per-case coverage would pin fifty near-identical branches for no additional mechanism.

        testCase "an empty option list serializes to the empty string" <| fun _ ->
            let conv = fun (_: System.IO.FileInfo) -> "PATH"
            let actual =
                PercolatorWrapper.Parameters.stringOf
                    conv
                    (PercolatorWrapper.Parameters.PercolatorParams.GeneralOptions [])
            Expect.equal actual "" "no options produce no arguments"
            // No options, no arguments.

        testCase "target and decoy output flags are distinct for all three result types" <| fun _ ->
            let conv = fun (_: System.IO.FileInfo) -> "P"
            let serialize option =
                PercolatorWrapper.Parameters.stringOf
                    conv
                    (PercolatorWrapper.Parameters.PercolatorParams.FileOutputOptions [option])
            let assertPair name target targetFlag decoy decoyFlag =
                let targetOutput = serialize target
                let decoyOutput = serialize decoy
                Expect.stringContains targetOutput targetFlag (sprintf "%s target output uses its target flag" name)
                Expect.isFalse (targetOutput.Contains "decoy") (sprintf "%s target output does not contain a decoy flag" name)
                Expect.stringContains decoyOutput decoyFlag (sprintf "%s decoy output uses its decoy flag" name)
                Expect.isTrue (decoyOutput.Contains "decoy") (sprintf "%s decoy output is distinct from target output" name)
            // swapping a target/decoy output pair silently reverses target-decoy semantics downstream - the one confusion the per-family representative test cannot catch.
            assertPair
                "peptides"
                (PercolatorWrapper.Parameters.POUTTAB_Peptides (FileInfo @"C:\x\p.tsv"))
                "--results-peptides"
                (PercolatorWrapper.Parameters.POUTTAB_DecoyPeptides (FileInfo @"C:\x\p.decoy.tsv"))
                "--decoy-results-peptides"
            assertPair
                "PSMs"
                (PercolatorWrapper.Parameters.POUTTAB_PSMs (FileInfo @"C:\x\p.psms"))
                "--results-psms"
                (PercolatorWrapper.Parameters.POUTTAB_DecoyPSMs (FileInfo @"C:\x\p.decoy.psms"))
                "--decoy-results-psms"
            assertPair
                "proteins"
                (PercolatorWrapper.Parameters.POUTTAB_Proteins (FileInfo @"C:\x\p.proteins"))
                "--results-proteins"
                (PercolatorWrapper.Parameters.POUTTAB_DecoyProteins (FileInfo @"C:\x\p.decoy.proteins"))
                "--decoy-results-proteins"

        testCase "semantically confusable option siblings map to their own flags" <| fun _ ->
            let conv = fun (_: System.IO.FileInfo) -> "P"
            let serializeGeneral option =
                PercolatorWrapper.Parameters.stringOf
                    conv
                    (PercolatorWrapper.Parameters.PercolatorParams.GeneralOptions [option])
            let mixMax = serializeGeneral PercolatorWrapper.Parameters.PostProcessing_MIXMAX
            let tdc = serializeGeneral PercolatorWrapper.Parameters.PostProcessing_TargetDecoyCompetition
            Expect.stringContains mixMax "--post-processing-mix-max" "MIXMAX uses its own post-processing flag"
            Expect.isFalse (mixMax.Contains "--post-processing-tdc") "MIXMAX does not use the TDC flag"
            Expect.stringContains tdc "--post-processing-tdc" "TDC uses its own post-processing flag"
            Expect.isFalse (tdc.Contains "--post-processing-mix-max") "TDC does not use the MIXMAX flag"

            let crossValidation =
                PercolatorWrapper.Parameters.stringOf
                    conv
                    (PercolatorWrapper.Parameters.PercolatorParams.SVMTrainingOptions [
                        PercolatorWrapper.Parameters.FDR_CrossValidation 0.02
                    ])
            let positiveExamples =
                PercolatorWrapper.Parameters.stringOf
                    conv
                    (PercolatorWrapper.Parameters.PercolatorParams.SVMTrainingOptions [
                        PercolatorWrapper.Parameters.FDR_PositiveExamples 0.03
                    ])
            Expect.stringContains crossValidation "--testFDR 0.02" "cross-validation FDR uses testFDR"
            Expect.isFalse (crossValidation.Contains "--trainFDR") "cross-validation FDR does not use trainFDR"
            Expect.stringContains positiveExamples "--trainFDR 0.03" "positive-example FDR uses trainFDR"
            Expect.isFalse (positiveExamples.Contains "--testFDR") "positive-example FDR does not use testFDR"

            let fidoAlpha =
                PercolatorWrapper.Parameters.stringOf
                    conv
                    (PercolatorWrapper.Parameters.PercolatorParams.ProteinInferenceOptions_FIDO [
                        PercolatorWrapper.Parameters.Alpha 0.25
                    ])
            let fidoBeta =
                PercolatorWrapper.Parameters.stringOf
                    conv
                    (PercolatorWrapper.Parameters.PercolatorParams.ProteinInferenceOptions_FIDO [
                        PercolatorWrapper.Parameters.Beta 0.5
                    ])
            let fidoGamma =
                PercolatorWrapper.Parameters.stringOf
                    conv
                    (PercolatorWrapper.Parameters.PercolatorParams.ProteinInferenceOptions_FIDO [
                        PercolatorWrapper.Parameters.Gamma 0.75
                    ])
            Expect.stringContains fidoAlpha "--fido-alpha 0.25" "Fido alpha uses its own flag and value"
            Expect.isFalse (fidoAlpha.Contains "0.5") "Fido alpha does not use beta's value"
            Expect.isFalse (fidoAlpha.Contains "0.75") "Fido alpha does not use gamma's value"
            Expect.stringContains fidoBeta "--fido-beta 0.5" "Fido beta uses its own flag and value"
            Expect.isFalse (fidoBeta.Contains "0.25") "Fido beta does not use alpha's value"
            Expect.isFalse (fidoBeta.Contains "0.75") "Fido beta does not use gamma's value"
            Expect.stringContains fidoGamma "--fido-gamma 0.75" "Fido gamma uses its own flag and value"
            Expect.isFalse (fidoGamma.Contains "0.25") "Fido gamma does not use alpha's value"
            Expect.isFalse (fidoGamma.Contains "0.5") "Fido gamma does not use beta's value"
            // train/test FDR and the three Fido priors are semantically distinct knobs; distinct values make a swapped mapping observable.

        // PENDING: percolator's CLI requires whitespace (or =) between a flag and its value; the
        // implementation concatenates e.g. "--Cneg" directly with the number ("--Cneg0.75"), which the
        // percolator binary rejects as an unknown flag. Same pattern for --seed.
        ptestCase "numeric training options separate flag and value with whitespace" <| fun _ ->
            let conv = fun (_: System.IO.FileInfo) -> "PATH"
            let options : PercolatorWrapper.Parameters.SVMTrainingOptions list = [
                PercolatorWrapper.Parameters.Cpos 0.5
                PercolatorWrapper.Parameters.Cneg 0.75
                PercolatorWrapper.Parameters.SeedRndNumberGenerator 7.0
            ]
            let s =
                PercolatorWrapper.Parameters.stringOf
                    conv
                    (PercolatorWrapper.Parameters.PercolatorParams.SVMTrainingOptions options)
            Expect.stringContains s "--Cpos 0.5" "Cpos separates its flag and value with whitespace"
            Expect.stringContains s "--Cneg 0.75" "Cneg separates its flag and value with whitespace"
            Expect.stringContains s "--seed 7" "seed separates its flag and value with whitespace"
            // percolator also accepts '=' as separator; this test deliberately asserts the whitespace form - a correct '='-emitting fix should update the assertion, not be blocked by it.
]
