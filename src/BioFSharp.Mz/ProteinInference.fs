namespace BioFSharp.Mz


open FSharpAux
open BioFSharp.PeptideClassification
open BioFSharp.IO.GFF3
open FSharpAux.IO.SchemaReader.Attribute
open System.Collections.Generic


module ProteinInference = 
    open BioFSharp.Elements

    type internal StructuralEqualityComparer<'a when 'a : equality >() = 
        interface IEqualityComparer<'a> with
            member x.Equals(a,b) = a = b
            member x.GetHashCode a = hash a 
    
    /// For a single peptide Sequence, contains information about all proteins it might originate from and its evidence class.
    type ProteinClassItem<'sequence> = 
        {
            GroupOfProteinIDs: string []
            PeptideSequence: 'sequence
            Class: PeptideEvidenceClass
        }

    /// For a group of proteins, contains information about all peptides that might be used for its quantification.
    type InferredProteinClassItem<'sequence> = 
        {
            GroupOfProteinIDs: string
            PeptideSequence: 'sequence []
            Class: PeptideEvidenceClass
        }    

    /// For a group of proteins, contains information about all peptides that might be used for its quantification and score calculated for it.
    type InferredProteinClassItemScored =
        {
            GroupOfProteinIDs: string
            PeptideSequence  : string []
            Class            : PeptideEvidenceClass
            TargetScore      : float
            DecoyScore       : float
            Decoy            : bool
            DecoyBigger      : bool
            FoundInDB        : bool
        }

    /// For a group of proteins, contains information about all peptides that might be used for its quantification and score / q-value calculated for it.
    type InferredProteinClassItemQValue =
        {
            GroupOfProteinIDs: string
            PeptideSequence  : string []
            Class            : PeptideEvidenceClass
            TargetScore      : float
            DecoyScore       : float
            QValue           : float
            Decoy            : bool
            DecoyBigger      : bool
            FoundInDB        : bool
        }

    /// For a group of proteins, contains information about all peptides that are put into the output file.
    type InferredProteinClassItemOut =
        {
            GroupOfProteinIDs: string
            PeptideSequence  : string
            Class            : PeptideEvidenceClass
            TargetScore      : float
            DecoyScore       : float
            QValue           : float
        }

    type PSMInput =
        {
            [<FieldAttribute("PepSequenceID")>]
            PepSequenceID   : int
            [<FieldAttribute("StringSequence")>]
            Seq             :string
            [<FieldAttribute("PercolatorScore")>]
            PercolatorScore : float
        }

    /// Input for QValue calulation
    type QValueInput =
        {
            Score    : float
            IsDecoy  : bool
        }

    type ScoreTargetDecoyCount =
        {
            Score      : float
            DecoyCount : float
            TargetCount: float
        }

    let private createInferredProteinClassItem proteinIDs evidenceClass peptideSequences = 
        {
            GroupOfProteinIDs = proteinIDs
            PeptideSequence = peptideSequences
            Class = evidenceClass
        }   

    let createProteinClassItem proteinIDs evidenceClass peptideSequence : ProteinClassItem<'sequence> = 
        {
            GroupOfProteinIDs = proteinIDs
            PeptideSequence = peptideSequence
            Class = evidenceClass
        }

    let createInferredProteinClassItemScored proteinIDs evidenceClass peptideSequences targetScore decoyScore isDecoy decoyBigger foundInDB =
        {
            GroupOfProteinIDs = proteinIDs
            PeptideSequence   = peptideSequences
            Class             = evidenceClass
            TargetScore       = targetScore
            DecoyScore        = decoyScore
            Decoy             = isDecoy
            DecoyBigger       = decoyBigger
            FoundInDB         = foundInDB
        }

    let createInferredProteinClassItemQValue proteinIDs evidenceClass peptideSequences targetScore decoyScore qValue isDecoy decoyBigger foundInDB =
        {
            GroupOfProteinIDs = proteinIDs
            PeptideSequence   = peptideSequences
            Class             = evidenceClass
            TargetScore       = targetScore
            DecoyScore        = decoyScore
            QValue            = qValue
            Decoy             = isDecoy
            DecoyBigger       = decoyBigger
            FoundInDB         = foundInDB
        }

    let createInferredProteinClassItemOut proteinIDs evidenceClass peptideSequences targetScore decoyScore qValue =
        {
            GroupOfProteinIDs = proteinIDs
            PeptideSequence   = peptideSequences
            Class             = evidenceClass
            TargetScore       = targetScore
            DecoyScore        = decoyScore
            QValue            = qValue
        }

    let createQValueInput score isDecoy =
        {
            Score     = score
            IsDecoy   = isDecoy
        }

    /// Gives the decoy and target count at a specific score
    let createScoreTargetDecoyCount score decoyCount targetCount =
        {
            Score       = score
            DecoyCount  = decoyCount
            TargetCount = targetCount
        }
 
    /// Used to decide wether overlapping groups of proteins should be kept or merged
    type IntegrationStrictness = 
        /// Results in the minimal set of proteins which still can explain all measured peptides
        | Minimal
        /// All protein groups stay intact
        | Maximal

    /// Used to decide which peptides should be used for quantification of protein groups
    type PeptideUsageForQuantification =
        /// Use only the best matching peptides
        | Minimal
        /// Use all available peptides which point to a protein group which includes the given protein group
        | Maximal
        /// Use all available peptides which point to a part of the protein group
        | MaximalInverse


    let private increaseClass c =
        match c with
        | PeptideEvidenceClass.C1a -> PeptideEvidenceClass.C1b
        | PeptideEvidenceClass.C1b -> PeptideEvidenceClass.C2a
        | PeptideEvidenceClass.C2a -> PeptideEvidenceClass.C2b
        | PeptideEvidenceClass.C2b -> PeptideEvidenceClass.C3a
        | PeptideEvidenceClass.C3a -> PeptideEvidenceClass.C3b

    let private mapClassName c =
        match c with
        | PeptideEvidenceClass.C1a -> "Class1A"
        | PeptideEvidenceClass.C1b -> "Class1B"
        | PeptideEvidenceClass.C2a -> "Class2A"
        | PeptideEvidenceClass.C2b -> "Class2B"
        | PeptideEvidenceClass.C3a -> "Class3A"
        | PeptideEvidenceClass.C3b -> "Class3B"

    /// Appends a collection of proteins to a single string
    let private proteinGroupToString (proteinGroup:string[]) =
        Array.reduce (fun x y ->  x + ";" + y) proteinGroup

    /// Contains all different classes with their associated proteins
    type private ClassMap =
        Map<PeptideEvidenceClass,BidirectionalDictionary<string,string>>

    module private ClassMap =
        ///Creates a new empty ClassMap
        let empty() : ClassMap =
            Map.empty
            |> Map.add PeptideEvidenceClass.C1a (BidirectionalDictionary<string,string>())
            |> Map.add PeptideEvidenceClass.C1b (BidirectionalDictionary<string,string>())
            |> Map.add PeptideEvidenceClass.C2a (BidirectionalDictionary<string,string>())
            |> Map.add PeptideEvidenceClass.C2b (BidirectionalDictionary<string,string>())
            |> Map.add PeptideEvidenceClass.C3a (BidirectionalDictionary<string,string>())
            |> Map.add PeptideEvidenceClass.C3b (BidirectionalDictionary<string,string>())

        ///Searches the Memory Map from highest ranking to lowest ranking class for any of the given proteinIds. If an overlap exists, returns its evidence class.
        let searchProts (m:ClassMap) (prots:string[]) =
            match Map.tryFindKey (fun c (d:BidirectionalDictionary<string,string>) ->
                    Array.exists (fun p -> d.ContainsKey p) prots
                ) m with
            | Some c -> c
            | None -> PeptideEvidenceClass.Unknown

        ///Adds the protein group to the given class in the memory map
        let addGroup (m:ClassMap) (c:PeptideEvidenceClass) prots =
            let name = proteinGroupToString prots
            Seq.iter (fun p -> m.[c].Add p name) prots
            m

    ///In a given class finds the overlap between the given and already added proteins and integrates them
    let private findAndIntegrate (m:ClassMap) (c:PeptideEvidenceClass) (proteins:string []) =
        let merge (dict:BidirectionalDictionary<string,string>) =
            let oldGroups =
                Array.choose (dict.TryGetByKey >> (Option.map Array.ofSeq)) proteins
                |> Array.concat
                |> Array.distinct
                |> Array.choose (fun pg -> Option.map (fun seqs -> seqs,pg) (dict.TryGetByValue pg))
            let newSet = Set.ofSeq proteins         
            oldGroups
            |> Array.iter (fun (x,y) -> printfn "%s" ((Seq.toArray >> proteinGroupToString) x) )
            oldGroups
            |> Array.iter (fun (seqs,group) ->
                let oldSet = Set.ofSeq seqs
                if Set.isProperSubset newSet oldSet then 
                    let newGroup = proteinGroupToString proteins
                    dict.RemoveValue group
                    Seq.iter (fun k -> dict.Add k newGroup) proteins
                elif Set.isProperSubset oldSet newSet then 
                    ()
                else
                    let protIntersect = Set.intersect newSet oldSet
                    let proteins = Set.toArray protIntersect
                    let newGroup = proteinGroupToString proteins
                    dict.RemoveValue group
                    Seq.iter (fun k -> dict.Add k newGroup) proteins
            ) 
            dict
        match ClassMap.searchProts m proteins with
        | PeptideEvidenceClass.Unknown ->
            printfn "Search did not return anything"
            ClassMap.addGroup m c proteins
        | foundC ->
            if foundC < c then
                m
            else
                Map.map (fun ca dict ->
                    if ca = c then merge dict
                    else dict) m
    
    let private createClassMap (integrationStrictness:IntegrationStrictness) (proteinClassItems: ProteinClassItem<'sequence> list ) =
        let rec loop (m:ClassMap) proteinClassItems remaining lookingForClass =
            match proteinClassItems with
            | (item:ProteinClassItem<'sequence>) :: l ->
                    if item.Class = lookingForClass then
                        let m' = 
                            if integrationStrictness = IntegrationStrictness.Maximal then 
                                ClassMap.addGroup m item.Class item.GroupOfProteinIDs
                            else  
                                (findAndIntegrate m item.Class item.GroupOfProteinIDs)
                        loop m' l remaining lookingForClass
                    else
                        loop m l (item::remaining) lookingForClass
            | [] ->
                if lookingForClass = PeptideEvidenceClass.C3b then
                    m
                else
                    loop m remaining [] (increaseClass lookingForClass)
        loop (ClassMap.empty()) proteinClassItems [] PeptideEvidenceClass.C1a

    /// Used to map the resulting protein groups to the sequences which are used for their quantification
    let createProteinToPepSequencesMap (proteinClassItems: ProteinClassItem<'sequence> list) =
        proteinClassItems
        |> List.collect (fun ci -> 
            let prots = ci.GroupOfProteinIDs
            List.init prots.Length (fun i -> prots.[i],ci)
            )
        |> List.groupBy fst
        |> List.map (fun (p,l) -> p, List.map snd l)
        |> Map.ofList
        //|> Map.ofList

    let inferSequences (integrationStrictness:IntegrationStrictness) (peptideUsageForQuantification:PeptideUsageForQuantification) (proteinClassItems: ProteinClassItem<'sequence> list) =
        ///Adds a value v to the set of values already associated with key k
        let classes : ClassMap = createClassMap integrationStrictness proteinClassItems
        /// Used to map the resulting protein groups to the sequences which are used for their quantification
        let seqDict = createProteinToPepSequencesMap proteinClassItems
        printfn "finish setup"
        classes
        |> Seq.collect (fun kv -> 
            let c = kv.Key
            let d = kv.Value
            d.GetArrayOfValues 
            |> Seq.map (fun group -> 
                let prots = Option.get (d.TryGetByValue group)
                let protSet = Set.ofSeq prots
                match peptideUsageForQuantification with
                | PeptideUsageForQuantification.Maximal ->                     
                    let potentialCandidates = 
                        prots
                        |> Seq.collect (fun prot -> 
                            seqDict.[prot]
                            |> Seq.choose (fun protClassItem -> 
                                let set = Set.ofArray protClassItem.GroupOfProteinIDs
                                if Set.isSubset protSet set then
                                    Some protClassItem
                                else 
                                    None
                                )
                            )                    
                    potentialCandidates
                    |> Seq.map (fun pci -> pci.PeptideSequence)
                    |> Seq.toArray
                    |> Array.distinct
                    |> createInferredProteinClassItem group c 
                | PeptideUsageForQuantification.Minimal -> 
                    let potentialCandidates = 
                        prots
                        |> Seq.collect (fun prot -> 
                            seqDict.[prot]
                            |> Seq.choose (fun protClassItem -> 
                                let set = Set.ofArray protClassItem.GroupOfProteinIDs
                                if Set.isSubset protSet set then
                                    Some protClassItem
                                else 
                                    None
                                )
                            )   
                    let withLength = 
                        potentialCandidates
                        |> Seq.map (fun pci -> pci.GroupOfProteinIDs.Length,pci)
                    let min = Seq.minBy fst withLength |> fst
                    Seq.filter (fst >> ((=) min)) withLength
                    |> Seq.map (snd >> (fun pci -> pci.PeptideSequence))
                    |> Seq.toArray
                    |> Array.distinct
                    |> createInferredProteinClassItem group c
                | PeptideUsageForQuantification.MaximalInverse ->    
                    let potentialCandidates = 
                        prots
                        |> Seq.collect (fun prot -> 
                            seqDict.[prot]
                            |> Seq.choose (fun protClassItem -> 
                                let set = Set.ofArray protClassItem.GroupOfProteinIDs
                                if Set.isSubset set protSet then
                                    Some protClassItem
                                else 
                                    None
                                )
                            )                    
                    potentialCandidates
                    |> Seq.map (fun pci -> pci.PeptideSequence)
                    |> Seq.toArray
                    |> Array.distinct
                    |> createInferredProteinClassItem group c                     
                )
            )

            