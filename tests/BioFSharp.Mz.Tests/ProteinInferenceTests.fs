module ProteinInferenceTests

open Expecto
open BioFSharp.Mz
open BioFSharp.PeptideClassification
open BioFSharp.FileFormats.GFF3

let private classC1a = BioFSharp.PeptideClassification.PeptideEvidenceClass.C1a

let private mkPsm pepSequenceID sequence score : ProteinInference.PSMInput =
    { PepSequenceID = pepSequenceID; Seq = sequence; Score = score }

let private mkScored decoyBetter =
    ProteinInference.createInferredProteinClassItemScored
        (ProteinInference.proteinGroupToString [|"P"|])
        classC1a
        [|"pep"|]
        10.0
        1.0
        false
        decoyBetter
        true

[<Tests>]
let tests =
    testList "ProteinInferenceTests" [
        testList "Helpers" [
            testCase "proteinGroupToString joins accessions with semicolons preserving order" <| fun _ ->
                Expect.equal
                    (ProteinInference.proteinGroupToString [|"P1";"P2";"P3"|])
                    "P1;P2;P3"
                    "protein group accessions are joined with semicolons in order"
                Expect.equal
                    (ProteinInference.proteinGroupToString [|"P9"|])
                    "P9"
                    "a single accession is preserved"
                // The semicolon-joined group string is the documented external representation of a protein group (order-preserving concatenation).

            testCase "removeModification strips lowercase modification codes and brackets down to the base sequence" <| fun _ ->
                Expect.equal
                    (ProteinInference.removeModification "[ox]MAGK")
                    "MAGK"
                    "bracketed lowercase modification codes are removed"
                Expect.equal
                    (ProteinInference.removeModification "AGSEK")
                    "AGSEK"
                    "an unmodified peptide is preserved"
                // Modified peptide strings carry bracketed lowercase codes; normalization must recover the plain uppercase residue sequence so modified and unmodified observations of the same peptide can be matched.

            testCase "createPeptideScoreMap keeps the best score per peptide sequence" <| fun _ ->
                let psms : ProteinInference.PSMInput list list =
                    [
                        [
                            { PepSequenceID = 1; Seq = "AA"; Score = 5.0 }
                            { PepSequenceID = 2; Seq = "AA"; Score = 9.0 }
                            { PepSequenceID = 3; Seq = "AA"; Score = 7.0 }
                        ]
                        [
                            { PepSequenceID = 4; Seq = "BB"; Score = 3.0 }
                        ]
                    ]
                let m = ProteinInference.createPeptideScoreMap psms
                Expect.equal m.["AA"] 9.0 "the best score is retained for AA"
                Expect.equal m.["BB"] 3.0 "the score for BB is retained"
                Expect.equal m.Count 2 "duplicate peptide sequences collapse to one map entry"
                // A peptide's evidence is its highest observed PSM score; the middle score 9.0 wins over the first 5.0, last 7.0, and sum 21.0 candidates.

            testCase "assignPeptideScores sums the member peptides' scores" <| fun _ ->
                Expect.floatClose
                    Accuracy.high
                    (ProteinInference.assignPeptideScores [|"AA";"BB"|] (Map.ofList ["AA",9.0;"BB",3.0]))
                    12.0
                    "member peptide scores are summed"
                Expect.floatClose
                    Accuracy.high
                    (ProteinInference.assignPeptideScores [||] (Map.ofList ["AA",9.0]))
                    0.0
                    "an empty peptide set has no evidence"
                // A protein's score aggregates its peptides' evidence additively (documented), and no peptides means no evidence.

            testCase "createReverseProteinScores totals decoy-protein evidence and drops proteins without any" <| fun _ ->
                let revProteins =
                    [| ("REV_P1", [|"AA";"BB"|]); ("REV_P2", [|"CC"|]) |]
                let scoreMap = Map.ofList ["AA",9.0;"BB",3.0]
                let m = ProteinInference.createReverseProteinScores revProteins scoreMap
                Expect.isTrue (m.ContainsKey "REV_P1") "REV_P1 has observed peptide evidence"
                let revP1Score, _ = m.["REV_P1"]
                Expect.floatClose Accuracy.high revP1Score 12.0 "REV_P1 receives the sum of its peptide scores"
                Expect.isFalse (m.ContainsKey "REV_P2") "REV_P2 has no observed peptide evidence"
                // A reversed (decoy) protein's score is the sum of its matched peptides' scores; an unobserved decoy protein is removed by the documented zero-total filter.

            testCase "assignDecoyScoreToTargetScore takes the best decoy score among the group's accessions" <| fun _ ->
                let decoyMap = Map.ofList ["P1", (5.0, [|"AA"|]); "P2", (8.0, [|"BB"|])]
                Expect.floatClose
                    Accuracy.high
                    (ProteinInference.assignDecoyScoreToTargetScore "P1;P2" decoyMap)
                    8.0
                    "the strongest decoy accession determines the group score"
                Expect.floatClose
                    Accuracy.high
                    (ProteinInference.assignDecoyScoreToTargetScore "P1;P9" decoyMap)
                    5.0
                    "a missing accession contributes no score"
                // NOTE: with any POSITIVE real score present the missing accession's 0.0 default cannot win; the negative-score hazard is covered by the pending test below.

            ptestCase "missing accessions do not outcompete negative decoy scores" <| fun _ ->
                // PENDING: Percolator-style scores can be negative. For group "P1;P9" where P1's decoy score
                // is -3.0 and P9 is absent, the best real decoy evidence is -3.0; the implementation's 0.0
                // missing-accession default competes in the max and wins, fabricating evidence.
                let decoyMap = Map.ofList ["P1", (-3.0, [|"AA"|])]
                Expect.floatClose
                    Accuracy.high
                    (ProteinInference.assignDecoyScoreToTargetScore "P1;P9" decoyMap)
                    -3.0
                    "a missing accession does not outcompete a negative real decoy score"
                // The decoy competitor of a target group is its strongest decoy counterpart.

            testCase "isGene and isRNA classify GFF lines by feature type" <| fun _ ->
                let geneEntry : GFFEntry =
                    {
                        Seqid = "chr1"
                        Source = "source"
                        Feature = "gene"
                        StartPos = 1
                        EndPos = 10
                        Score = 0.0
                        Strand = '+'
                        Phase = 0
                        Attributes = Map.empty
                        Supplement = [||]
                    }
                let rnaEntry : GFFEntry =
                    {
                        Seqid = "chr1"
                        Source = "source"
                        Feature = "mRNA"
                        StartPos = 1
                        EndPos = 10
                        Score = 0.0
                        Strand = '+'
                        Phase = 0
                        Attributes = Map.empty
                        Supplement = [||]
                    }
                let geneLine : GFFLine<seq<char>> = GFFEntryLine geneEntry
                let rnaLine : GFFLine<seq<char>> = GFFEntryLine rnaEntry
                Expect.isTrue (ProteinInference.isGene geneLine) "isGene recognizes the gene feature"
                Expect.isFalse (ProteinInference.isGene rnaLine) "isGene rejects the mRNA feature"
                Expect.isSome (ProteinInference.isRNA rnaLine) "isRNA recognizes the mRNA feature"
                Expect.isNone (ProteinInference.isRNA geneLine) "isRNA rejects the gene feature"
                // The documented feature-string predicates classify GFF lines by feature type.
        ]

        testList "Inference" [
            testCase "inferSequences respects the peptide-usage subset rules on nested protein groups" <| fun _ ->
                let items =
                    [
                        ProteinInference.createProteinClassItem [|"P1"|] classC1a "pa"
                        ProteinInference.createProteinClassItem [|"P1";"P2"|] classC1a "pb"
                        ProteinInference.createProteinClassItem [|"P1";"P2";"P3"|] classC1a "pc"
                    ]
                let rMax =
                    ProteinInference.inferSequences
                        ProteinInference.IntegrationStrictness.Maximal
                        ProteinInference.PeptideUsageForQuantification.Maximal
                        items
                    |> List.ofSeq
                let rInv =
                    ProteinInference.inferSequences
                        ProteinInference.IntegrationStrictness.Maximal
                        ProteinInference.PeptideUsageForQuantification.MaximalInverse
                        items
                    |> List.ofSeq

                let expectedPeptides = set ["pa"; "pb"; "pc"]
                let allPeptides (r: ProteinInference.InferredProteinClassItem<string> list) =
                    r
                    |> Seq.collect (fun item -> item.PeptideSequence)
                    |> Set.ofSeq
                Expect.equal
                    (allPeptides rMax)
                    expectedPeptides
                    "the Maximal inferred groups retain all peptide evidence"
                Expect.equal
                    (allPeptides rInv)
                    expectedPeptides
                    "the MaximalInverse inferred groups retain all peptide evidence"

                match rMax |> List.tryFind (fun item -> item.GroupOfProteinIDs = "P1") with
                | Some group ->
                    Expect.equal
                        (Set.ofArray group.PeptideSequence)
                        expectedPeptides
                        "Maximal usage assigns every peptide whose source protein set contains P1"
                | None ->
                    failtestf "Maximal inferred group structure did not contain P1; observed: %A" rMax

                match rInv |> List.tryFind (fun item -> item.GroupOfProteinIDs = "P1") with
                | Some group ->
                    Expect.equal
                        (Set.ofArray group.PeptideSequence)
                        (set ["pa"])
                        "MaximalInverse usage assigns only peptides whose source set is contained in P1"
                | None ->
                    failtestf "MaximalInverse inferred group structure did not contain P1; observed: %A" rInv
                // The two usage options are documented as opposite subset relations between a peptide's source proteins and the inferred group; the nested construction distinguishes them cleanly.

            testCase "inferSequences keeps unrelated evidence in separate groups" <| fun _ ->
                let items =
                    [
                        ProteinInference.createProteinClassItem [|"P1"|] classC1a "x"
                        ProteinInference.createProteinClassItem [|"P1";"P2"|] classC1a "y"
                        ProteinInference.createProteinClassItem [|"P3"|] classC1a "z"
                    ]
                let r =
                    ProteinInference.inferSequences
                        ProteinInference.IntegrationStrictness.Maximal
                        ProteinInference.PeptideUsageForQuantification.Maximal
                        items
                    |> List.ofSeq
                let groupProteins (item: ProteinInference.InferredProteinClassItem<string>) =
                    item.GroupOfProteinIDs.Split(';') |> Set.ofArray
                let p3Group = r |> List.tryFind (fun item -> item.GroupOfProteinIDs = "P3")
                Expect.isSome p3Group "P3 has its own inferred group"
                match p3Group with
                | Some group ->
                    Expect.equal
                        (Set.ofArray group.PeptideSequence)
                        (set ["z"])
                        "the unrelated P3 group contains only z"
                | None -> ()
                Expect.isTrue
                    (r |> List.forall (fun item ->
                        let proteins = groupProteins item
                        not (Set.contains "P3" proteins && (Set.contains "P1" proteins || Set.contains "P2" proteins))))
                    "unrelated P3 evidence is not merged with P1 or P2"
                let allPeptides = r |> Seq.collect (fun item -> item.PeptideSequence) |> Set.ofSeq
                Expect.equal allPeptides (set ["x"; "y"; "z"]) "all peptide evidence is retained"
                // Under IntegrationStrictness.Maximal the merge logic never runs (documented "groups stay intact"), so this test pins peptide assignment and evidence conservation, not merge behavior.

            testCase "minimal integration also keeps disjoint evidence separate" <| fun _ ->
                let items =
                    [
                        ProteinInference.createProteinClassItem [|"P1"|] classC1a "x"
                        ProteinInference.createProteinClassItem [|"P1";"P2"|] classC1a "y"
                        ProteinInference.createProteinClassItem [|"P3"|] classC1a "z"
                    ]
                let r =
                    ProteinInference.inferSequences
                        ProteinInference.IntegrationStrictness.Minimal
                        ProteinInference.PeptideUsageForQuantification.Maximal
                        items
                    |> List.ofSeq
                let groupProteins (item: ProteinInference.InferredProteinClassItem<string>) =
                    item.GroupOfProteinIDs.Split(';') |> Set.ofArray
                let p3Group = r |> List.tryFind (fun item -> item.GroupOfProteinIDs = "P3")
                Expect.isSome p3Group "P3 has its own inferred group"
                match p3Group with
                | Some group ->
                    Expect.equal
                        (Set.ofArray group.PeptideSequence)
                        (set ["z"])
                        "the disjoint P3 group contains only z"
                | None -> ()
                Expect.isTrue
                    (r |> List.forall (fun item ->
                        let proteins = groupProteins item
                        not (Set.contains "P3" proteins && (Set.contains "P1" proteins || Set.contains "P2" proteins))))
                    "disjoint P3 evidence is not merged with P1 or P2"
                let allPeptides = r |> Seq.collect (fun item -> item.PeptideSequence) |> Set.ofSeq
                Expect.equal allPeptides (set ["x"; "y"; "z"]) "all peptide evidence is retained"
                // Minimal strictness is the path where findAndIntegrate's subset/intersection merging actually executes; disjoint evidence must survive it untouched.

            ptestCase "evidence of every class survives inference" <| fun _ ->
                // PENDING: an item whose evidence class is Unknown silently vanishes from the output - the
                // class loop visits C1a..C3b and discards the remainder, losing peptide evidence without a
                // trace. Argued-correct: every input peptide appears in some output group.
                let items =
                    [
                        ProteinInference.createProteinClassItem [|"P1"|] classC1a "x"
                        ProteinInference.createProteinClassItem [|"P2"|] PeptideEvidenceClass.Unknown "u"
                    ]
                let r =
                    ProteinInference.inferSequences
                        ProteinInference.IntegrationStrictness.Maximal
                        ProteinInference.PeptideUsageForQuantification.Maximal
                        items
                    |> List.ofSeq
                let allPeptides = r |> Seq.collect (fun item -> item.PeptideSequence) |> Set.ofSeq
                Expect.isTrue
                    (Set.isSubset (set ["x"; "u"]) allPeptides)
                    "the union of inferred groups contains every input peptide"

            testCase "maximal integration preserves every nested group with superset peptide assignment" <| fun _ ->
                let items =
                    [
                        ProteinInference.createProteinClassItem [|"P1"|] classC1a "pa"
                        ProteinInference.createProteinClassItem [|"P1";"P2"|] classC1a "pb"
                        ProteinInference.createProteinClassItem [|"P1";"P2";"P3"|] classC1a "pc"
                    ]
                let result =
                    ProteinInference.inferSequences
                        ProteinInference.IntegrationStrictness.Maximal
                        ProteinInference.PeptideUsageForQuantification.Maximal
                        items
                    |> List.ofSeq
                let groupProteins (item: ProteinInference.InferredProteinClassItem<string>) =
                    item.GroupOfProteinIDs.Split(';') |> Set.ofArray
                let requireGroup expected =
                    match result |> List.tryFind (fun item -> groupProteins item = expected) with
                    | Some group -> group
                    | None -> failtestf "D18 observed group structure: %A" result
                let p1 = requireGroup (set ["P1"])
                let p1p2 = requireGroup (set ["P1";"P2"])
                let p1p2p3 = requireGroup (set ["P1";"P2";"P3"])
                Expect.equal (Set.ofArray p1.PeptideSequence) (set ["pa";"pb";"pc"]) "the singleton group receives every containing-source peptide"
                Expect.equal (Set.ofArray p1p2.PeptideSequence) (set ["pb";"pc"]) "the two-protein group receives every containing-source peptide"
                Expect.equal (Set.ofArray p1p2p3.PeptideSequence) (set ["pc"]) "the three-protein group receives its containing-source peptide"
                // documented: Maximal keeps all groups intact and usage-Maximal assigns each group every peptide whose source set CONTAINS it - the nested construction makes each expected set distinct.

            testCase "minimal integration collapses overlapping groups to their intersection" <| fun _ ->
                let items =
                    [
                        ProteinInference.createProteinClassItem [|"P1";"P2"|] classC1a "x"
                        ProteinInference.createProteinClassItem [|"P2";"P3"|] classC1a "y"
                    ]
                let run orderedItems =
                    let result =
                        ProteinInference.inferSequences
                            ProteinInference.IntegrationStrictness.Minimal
                            ProteinInference.PeptideUsageForQuantification.Maximal
                            orderedItems
                        |> List.ofSeq
                    let groupProteins (item: ProteinInference.InferredProteinClassItem<string>) =
                        item.GroupOfProteinIDs.Split(';') |> Set.ofArray
                    match result |> List.tryFind (fun item -> groupProteins item = set ["P2"]) with
                    | Some group -> result, group
                    | None -> failtestf "D19 observed group structure: %A" result
                let resultA, groupA = run items
                let resultB, groupB = run (List.rev items)
                Expect.equal (Set.ofArray groupA.PeptideSequence) (set ["x";"y"]) "the intersection group carries both overlapping peptides"
                Expect.equal (Set.ofArray groupB.PeptideSequence) (set ["x";"y"]) "reversing evidence order preserves the intersection group's peptides"
                Expect.equal (resultA |> Seq.collect (fun item -> item.PeptideSequence) |> Set.ofSeq) (set ["x";"y"]) "the first input order retains the union of peptides"
                Expect.equal (resultB |> Seq.collect (fun item -> item.PeptideSequence) |> Set.ofSeq) (set ["x";"y"]) "the reversed input order retains the union of peptides"
                // parsimony: P2 alone explains both peptides, so minimal integration must report the intersection - and inference must not depend on evidence order.

            testCase "usage Minimal assigns each group only its best-matching source evidence" <| fun _ ->
                let items =
                    [
                        ProteinInference.createProteinClassItem [|"P1"|] classC1a "pa"
                        ProteinInference.createProteinClassItem [|"P1";"P2"|] classC1a "pb"
                        ProteinInference.createProteinClassItem [|"P1";"P2";"P3"|] classC1a "pc"
                    ]
                let result =
                    ProteinInference.inferSequences
                        ProteinInference.IntegrationStrictness.Maximal
                        ProteinInference.PeptideUsageForQuantification.Minimal
                        items
                    |> List.ofSeq
                let groupProteins (item: ProteinInference.InferredProteinClassItem<string>) =
                    item.GroupOfProteinIDs.Split(';') |> Set.ofArray
                let requireGroup expected =
                    match result |> List.tryFind (fun item -> groupProteins item = expected) with
                    | Some group -> group
                    | None -> failtestf "D20 observed group structure: %A" result
                let p1 = requireGroup (set ["P1"])
                let p1p2 = requireGroup (set ["P1";"P2"])
                let p1p2p3 = requireGroup (set ["P1";"P2";"P3"])
                Expect.equal (Set.ofArray p1.PeptideSequence) (set ["pa"]) "the singleton group gets only its best-matching source evidence"
                Expect.equal (Set.ofArray p1p2.PeptideSequence) (set ["pb"]) "the two-protein group gets only its best-matching source evidence"
                Expect.equal (Set.ofArray p1p2p3.PeptideSequence) (set ["pc"]) "the three-protein group gets only its best-matching source evidence"
                // documented: usage-Minimal restricts each group's candidates (sources containing the group) to those with the SMALLEST source set - for nested sources that is exactly the group's own evidence.

            testCase "higher-confidence evidence classes are not displaced by overlapping lower classes" <| fun _ ->
                let items =
                    [
                        ProteinInference.createProteinClassItem [|"P1";"P2"|] classC1a "strong"
                        ProteinInference.createProteinClassItem [|"P2";"P3"|] PeptideEvidenceClass.C2a "weak"
                    ]
                let run orderedItems =
                    let result =
                        ProteinInference.inferSequences
                            ProteinInference.IntegrationStrictness.Minimal
                            ProteinInference.PeptideUsageForQuantification.Maximal
                            orderedItems
                        |> List.ofSeq
                    let candidate =
                        result
                        |> List.tryFind (fun item ->
                            let proteins = item.GroupOfProteinIDs.Split(';') |> Set.ofArray
                            Set.isSubset (set ["P1";"P2"]) proteins)
                    match candidate with
                    | Some group -> group
                    | None -> failtestf "D21 observed group structure: %A" result
                let groupA = run items
                let groupB = run (List.rev items)
                Expect.equal groupA.Class classC1a "the stronger C1a evidence class survives the first order"
                Expect.equal groupB.Class classC1a "the stronger C1a evidence class survives the reversed order"
                Expect.isTrue (Array.contains "strong" groupA.PeptideSequence) "the stronger peptide survives the first order"
                Expect.isTrue (Array.contains "strong" groupB.PeptideSequence) "the stronger peptide survives the reversed order"
                // evidence-class precedence: C1a outranks C2a, so overlapping weaker evidence must not overwrite or displace the stronger group's class or membership.
        ]

        testList "FDR" [
            testCase "calculateFDRwithDecoyTargetRatio is the decoy/target count ratio" <| fun _ ->
                let mkScored decoyBetter =
                    ProteinInference.createInferredProteinClassItemScored
                        (ProteinInference.proteinGroupToString [|"P"|])
                        classC1a
                        [|"pep"|]
                        10.0
                        1.0
                        false
                        decoyBetter
                        true
                Expect.floatClose
                    Accuracy.high
                    (ProteinInference.calculateFDRwithDecoyTargetRatio [| mkScored false; mkScored false; mkScored true |])
                    0.5
                    "one decoy-better item among two target-better items gives 1/2"
                Expect.floatClose
                    Accuracy.high
                    (ProteinInference.calculateFDRwithDecoyTargetRatio [| mkScored false; mkScored false |])
                    0.0
                    "no decoy wins gives zero estimated FDR"
                // The standard decoy/target FDR estimate is hand-counted here as 1/2 and 0.

            testCase "modified PSM sequences match unmodified reverse-digest peptides" <| fun _ ->
                let m =
                    ProteinInference.createReverseProteinScores
                        [| ("REV_P1", [|"MAGK"|]) |]
                        (Map.ofList ["M[ox]AGK", 9.0])
                Expect.isTrue (m.ContainsKey "REV_P1") "the modified PSM matches the reverse-digest peptide after normalization"
                Expect.equal (fst m.["REV_P1"]) 9.0 "the normalized reverse protein receives the modified PSM score"
                // the score map's keys are normalized with removeModification before matching, so a modified observation must still transfer its evidence to the base-sequence decoy peptide.

            testCase "the best PSM score survives when all scores are negative" <| fun _ ->
                let scoreMap =
                    ProteinInference.createPeptideScoreMap [
                        [ mkPsm 1 "AA" -5.0; mkPsm 2 "AA" -2.0; mkPsm 3 "AA" -7.0 ]
                    ]
                Expect.equal scoreMap.["AA"] -2.0 "the maximum negative PSM score is retained"
                // Percolator-style scores can be all-negative; a zero floor or first/last-wins would return the wrong evidence - the hand-evident maximum is -2.

            testCase "MAYU length bins stratify by sequence length" <| fun _ ->
                let inferred =
                    [|
                        for i in 1..4 do
                            yield
                                ProteinInference.createInferredProteinClassItemScored
                                    (ProteinInference.proteinGroupToString [|sprintf "P%d" i|])
                                    classC1a
                                    [|sprintf "pep%d" i|]
                                    10.0
                                    1.0
                                    false
                                    false
                                    true
                    |]
                let proteinsFromDB =
                    [|
                        ("P1", String.replicate 10 "A")
                        ("P2", String.replicate 20 "A")
                        ("P3", String.replicate 100 "A")
                        ("P4", String.replicate 200 "A")
                    |]
                let bins = ProteinInference.MAYU.binProteinsLength inferred proteinsFromDB 2.0
                let accessionSets =
                    bins
                    |> Array.map (fun bin ->
                        bin
                        |> Array.collect (fun item -> item.GroupOfProteinIDs.Split(';'))
                        |> Set.ofArray)
                Expect.equal accessionSets [|set ["P1";"P2"]; set ["P3";"P4"]|] "the two shortest proteins share MAYU's low length bin"
                // MAYU's model conditions FDR on protein length: the two shortest must share the low bin - directly from the documented length stratification, no incidental fields locked.

            // PENDING: with zero target wins, decoys/targets = 2/0 is IEEE Infinity; an FDR is a
            // proportion and its conservative bounded value is 1.0 (the module's MAYU path already applies
            // that bound).
            ptestCase "the decoy/target ratio is a bounded FDR" <| fun _ ->
                Expect.equal
                    (ProteinInference.calculateFDRwithDecoyTargetRatio [| mkScored true; mkScored true |])
                    1.0
                    "an all-decoy result is conservatively bounded at one"

            ptestCase "a decoy-free confident dataset has zero MAYU FDR" <| fun _ ->
                // PENDING: for a dataset of only confident targets (no decoy wins, all found in the DB) the
                // estimated number of false positives is 0, so the FDR must be 0. The implementation's
                // zero/non-finite fallback returns 1.0 instead - perfect data gets the worst possible FDR and
                // any downstream threshold rejects the entire dataset.
                let inferred =
                    [|
                        for i in 1..4 do
                            yield
                                ProteinInference.createInferredProteinClassItemScored
                                    (ProteinInference.proteinGroupToString [|sprintf "P%d" i|])
                                    classC1a
                                    [|sprintf "pep%d" i|]
                                    10.0
                                    0.0
                                    false
                                    false
                                    true
                    |]
                let proteinsFromDB =
                    [|
                        ("P1", "AAAA")
                        ("P2", "AAAAAA")
                        ("P3", "AAAAAAAA")
                        ("P4", "AAAAAAAAAA")
                    |]
                Expect.floatClose
                    Accuracy.high
                    (ProteinInference.calculateFDRwithMAYU inferred proteinsFromDB)
                    0.0
                    "a decoy-free confident dataset has zero MAYU FDR"

            testCase "MAYU.binProteinsLength conserves all inferred and database proteins across bins" <| fun _ ->
                let inferred =
                    [|
                        for i in 1..4 do
                            yield
                                ProteinInference.createInferredProteinClassItemScored
                                    (ProteinInference.proteinGroupToString [|sprintf "P%d" i|])
                                    classC1a
                                    [|"pep"|]
                                    1.0
                                    0.5
                                    false
                                    false
                                    true
                    |]
                let proteinsFromDB =
                    [|
                        for i in 1..6 do
                            yield (sprintf "P%d" i, String.replicate (10 * i) "A")
                    |]
                let bins = ProteinInference.MAYU.binProteinsLength inferred proteinsFromDB 2.0
                Expect.equal
                    (bins |> Array.sumBy Array.length)
                    6
                    "all inferred and unmatched database proteins occur exactly once"
                let allGroupOfProteinIDs =
                    bins
                    |> Array.collect (fun bin -> bin)
                    |> Array.map (fun item -> item.GroupOfProteinIDs)
                    |> Array.sort
                    |> Array.toList
                Expect.equal
                    allGroupOfProteinIDs
                    ["P1"; "P2"; "P3"; "P4"; "P5"; "P6"]
                    "all protein groups are conserved as an identity multiset"
                let syntheticEntries =
                    bins
                    |> Array.collect (fun bin -> bin)
                    |> Array.filter (fun item -> not item.FoundInDB)
                Expect.equal syntheticEntries.Length 2 "the two unmatched database proteins become synthetic entries"
                Expect.equal
                    (syntheticEntries |> Array.map (fun item -> item.GroupOfProteinIDs) |> Set.ofArray)
                    (set ["P5"; "P6"])
                    "entries with FoundInDB = false are exactly P5 and P6"
                // Conservation as identity, not count - a duplicate-P1-drop-P5 bug passes a count check but not this.

            testCase "assignQValueToIPCIS applies the q-value function to the appropriate score" <| fun _ ->
                let target =
                    ProteinInference.createInferredProteinClassItemScored
                        (ProteinInference.proteinGroupToString [|"P1"|])
                        classC1a
                        [|"pep"|]
                        5.0
                        1.0
                        false
                        false
                        true
                let decoy =
                    ProteinInference.createInferredProteinClassItemScored
                        (ProteinInference.proteinGroupToString [|"P2"|])
                        classC1a
                        [|"pep"|]
                        1.0
                        3.0
                        true
                        true
                        true
                let f = fun (s: float) -> s * 0.1
                Expect.floatClose
                    Accuracy.high
                    (ProteinInference.assignQValueToIPCIS f target).QValue
                    0.5
                    "a target item uses its target score"
                Expect.floatClose
                    Accuracy.high
                    (ProteinInference.assignQValueToIPCIS f decoy).QValue
                    0.3
                    "a decoy item uses its decoy score"
                // A target item is scored by its target score, a decoy item by its decoy score.
        ]
    ]
