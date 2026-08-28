module ParametersTests

open System
open System.IO
open Expecto
open BioFSharp.Mz
open BioFSharp.WorkflowLanguage

let someWaveletParameters : SignalDetection.Wavelet.WaveletParameters = {
    NumberOfScales = 10
    YThreshold = 1.0
    MzTolerance = 0.05
    SNRS_Percentile = 95.0
    MinSNR = 1.0
    RefineMZ = true
    SumIntensities = false
}

let somePaddingParameters : SignalDetection.Padding.PaddingParameters = {
    PaddingYValue = 0.0
    MaximumPaddingPoints = Some 5
    MzTolerance = 0.05
    WindowSize = 3
    SpacingPerc = 50.0
}

let ppParams = MSProcessing.createPeakPickingParams true None None None

let ppParamsFull =
    MSProcessing.createPeakPickingParams false (Some someWaveletParameters) (Some somePaddingParameters) None

let op1 : Definition.Operation<MSProcessing.MSParameters> =
    {
        Id = Guid.Parse("00000000-0000-0000-0000-000000000001")
        Name = "peak-picking-minimal"
        Operator = "peak-picking-operator-1"
        Input = Definition.IOType.File (FileInfo("input-1.mzML"))
        Output = Definition.IOType.File (FileInfo("output-1.mzML"))
        Parameters = MSProcessing.MSParameters.PeakPicking ppParams
    }

let op2 : Definition.Operation<MSProcessing.MSParameters> =
    {
        Id = Guid.Parse("00000000-0000-0000-0000-000000000002")
        Name = "peak-picking-full"
        Operator = "peak-picking-operator-2"
        Input = Definition.IOType.File (FileInfo("input-2.mzML"))
        Output = Definition.IOType.File (FileInfo("output-2.mzML"))
        Parameters = MSProcessing.MSParameters.PeakPicking ppParamsFull
    }

[<Tests>]
let tests =
    testList "ParametersTests" [
        testCase "createPeakPickingParams preserves every field including option cases" <| fun _ ->
            Expect.equal ppParamsFull.CompressData false "CompressData preserves false"
            Expect.equal ppParamsFull.Ms1PeakPicking (Some someWaveletParameters) "Ms1PeakPicking preserves Some wavelet parameters"
            Expect.equal ppParamsFull.PaddingParameters (Some somePaddingParameters) "PaddingParameters preserves Some padding parameters"
            Expect.equal ppParamsFull.Ms2PeakPicking None "Ms2PeakPicking preserves None"
            Expect.equal ppParams.Ms1PeakPicking None "minimal Ms1PeakPicking preserves None"
            // a pure record builder; the asymmetric Some/None fixture exists because Ms1PeakPicking and Ms2PeakPicking share a type, so a positional transposition would compile silently and only this asymmetry catches it.

        testCase "the BioFSharp.Mz payload types round-trip through JSON on their own" <| fun _ ->
            let ppParamsRoundTripped =
                ppParamsFull
                |> Newtonsoft.Json.JsonConvert.SerializeObject
                |> Newtonsoft.Json.JsonConvert.DeserializeObject<MSProcessing.PeakPickingParams>

            let msParameters = MSProcessing.MSParameters.PeakPicking ppParamsFull
            let msParametersRoundTripped =
                msParameters
                |> Newtonsoft.Json.JsonConvert.SerializeObject
                |> Newtonsoft.Json.JsonConvert.DeserializeObject<MSProcessing.MSParameters>

            Expect.equal ppParamsRoundTripped ppParamsFull "PeakPickingParams JSON round trip preserves the payload"
            Expect.equal msParametersRoundTripped msParameters "MSParameters JSON round trip preserves the parameter DU"
            // isolation evidence: every BioFSharp.Mz-owned type in a workflow payload (records, options, the parameter DU) serializes and round-trips perfectly on this Newtonsoft version - the pended Operation round trip fails solely on the BioFSharp Operation type's embedded FileInfo.

        ptestCase "a single operation survives the JSON round trip" <| fun _ ->
            let roundTripped = MSProcessing.operationOfJson (MSProcessing.operationToJSon op1)
            let projectOperation (operation: Definition.Operation<MSProcessing.MSParameters>) =
                let projectIO = function
                    | Definition.IOType.File file -> Choice1Of2 file.FullName
                    | Definition.IOType.Files files -> Choice2Of2 (files |> List.map (fun file -> file.FullName))
                operation.Id, operation.Name, operation.Operator, operation.Parameters, projectIO operation.Input, projectIO operation.Output
            Expect.equal (projectOperation roundTripped) (projectOperation op1) "operation JSON round trip preserves the operation identity, parameters, and file paths"
            // The module exposes operationToJSon/operationOfJson as an explicit serialize/deserialize pair; a pair that does not invert is useless for its documented workflow-persistence purpose. No JSON shape is asserted - only the semantic round trip.
            // Observed behavior: Newtonsoft.Json.JsonSerializationException, "Unable to serialize instance of 'System.IO.FileInfo'."
            // SECOND obstacle (independent of serialization): FileInfo has only reference equality, so
            // whole-record equality on a round-tripped Operation can never hold - the projection compares
            // paths, which is what file identity means for a persisted workflow.

        ptestCase "a workflow of operations survives the JSON round trip preserving count and order" <| fun _ ->
            let rt = MSProcessing.workFlowToJson (MSProcessing.workFlowToJSon [op1; op2]) |> List.ofSeq
            let projectOperation (operation: Definition.Operation<MSProcessing.MSParameters>) =
                let projectIO = function
                    | Definition.IOType.File file -> Choice1Of2 file.FullName
                    | Definition.IOType.Files files -> Choice2Of2 (files |> List.map (fun file -> file.FullName))
                operation.Id, operation.Name, operation.Operator, operation.Parameters, projectIO operation.Input, projectIO operation.Output
            Expect.equal rt.Length 2 "workflow round trip preserves count"
            Expect.equal (projectOperation rt.[0]) (projectOperation op1) "workflow round trip preserves first operation"
            Expect.equal (projectOperation rt.[1]) (projectOperation op2) "workflow round trip preserves second operation"
            // The same pairing argument applies for sequences; order is part of a workflow's meaning.
            // Observed behavior: Newtonsoft.Json.JsonSerializationException, "Unable to serialize instance of 'System.IO.FileInfo'."
            // SECOND obstacle (independent of serialization): FileInfo has only reference equality, so
            // whole-record equality on a round-tripped Operation can never hold - the projection compares
            // paths, which is what file identity means for a persisted workflow.

        testCase "an empty workflow round-trips to an empty workflow" <| fun _ ->
            let rt =
                MSProcessing.workFlowToJson (
                    MSProcessing.workFlowToJSon ([]: BioFSharp.WorkflowLanguage.Definition.Operation<MSProcessing.MSParameters> list)
                )
                |> List.ofSeq
            Expect.equal rt [] "empty workflow round trip preserves emptiness"
    ]
