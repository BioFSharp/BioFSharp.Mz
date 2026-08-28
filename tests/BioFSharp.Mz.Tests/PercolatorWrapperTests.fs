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
