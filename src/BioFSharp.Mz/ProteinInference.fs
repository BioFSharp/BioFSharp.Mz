namespace BioFSharp.Mz


open FSharpAux
open BioFSharp.PeptideClassification
open BioFSharp.IO.GFF3
open FSharpAux.IO.SchemaReader.Attribute
open System.Collections.Generic
open FSharp.Stats
open FSharp.Stats.Fitting.NonLinearRegression
open FDRControl
open FDRControl.MAYU


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
            GroupOfProteinIDs  : string
            PeptideSequence    : string []
            Class              : PeptideEvidenceClass
            TargetScore        : float
            DecoyScore         : float
            Decoy              : bool
            DecoyHasBetterScore: bool
            FoundInDB          : bool
        }

    /// For a group of proteins, contains information about all peptides that might be used for its quantification and score / q-value calculated for it.
    type InferredProteinClassItemQValue =
        {
            InfProtClassItem : InferredProteinClassItemScored
            QValue           : float
        }

    /// For a group of proteins, contains information about all peptides that are put into the output file.
    type Result =
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
            Score : float
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

    let createInferredProteinClassItemScored proteinIDs evidenceClass peptideSequences targetScore decoyScore isDecoy decoyHasBetterScore foundInDB =
        {
            GroupOfProteinIDs   = proteinIDs
            PeptideSequence     = peptideSequences
            Class               = evidenceClass
            TargetScore         = targetScore
            DecoyScore          = decoyScore
            Decoy               = isDecoy
            DecoyHasBetterScore = decoyHasBetterScore
            FoundInDB           = foundInDB
        }

    let createInferredProteinClassItemQValue infProtClassItemScored qValue=
        {
            InfProtClassItem = infProtClassItemScored
            QValue           = qValue
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
    let proteinGroupToString (proteinGroup:string[]) =
        Array.reduce (fun x y ->  x + ";" + y) proteinGroup

    let removeModification pepSeq =
        String.filter (fun c -> System.Char.IsLower c |> not && c <> '[' && c <> ']') pepSeq

    /// Checks if GFF line describes gene
    let isGene (item: GFFLine<seq<char>>) =
        match item with
        | GFFEntryLine x -> x.Feature = "gene"
        | _ -> false

    /// Checks if GFF line describes rna
    let isRNA (item: GFFLine<seq<char>>) =
        match item with
        | GFFEntryLine x -> if x.Feature = "mRNA" then Some x else None
        | _ -> None

    /// Reads geographical information about protein from gff entry and builds the modelinfo of it
    /// This function takes an RNA gff3 entry and therefore will contain the splice variant id of the gene in its result.
    /// This splice variant id should be the same in the given FastA-file.
    let createProteinModelInfoFromEntry i locus (entry:GFFEntry) =
        let attributes = entry.Attributes
        /// Same as in FastA file
        let spliceVariantID =
            match Map.tryFind "Name" attributes with
            | Some res ->
                res.Head
            | None ->
                failwithf "could not find spliceVariantId for locus %s" locus
        let chromosomeID = entry.Seqid

        let direction =
            match entry.Strand with
            |'+' -> StrandDirection.Forward
            |'-' -> StrandDirection.Reverse

        createProteinModelInfo spliceVariantID chromosomeID direction locus i Seq.empty Seq.empty

    /// By reading GFF creates the protein models (relationships of proteins to each other) which basically means grouping the rnas over the gene loci
    /// TODO: Don't group over order but rather group over id
    let assignTranscriptsToGenes regexPattern (gffLines: seq<GFFLine<seq<char>>>)  =
        gffLines
        // transcripts are grouped by the gene they originate from
        |> Seq.groupWhen isGene
        |> Seq.map (fun group ->
            match Seq.head group with
            | GFFEntryLine x ->
                let locus = x.Attributes.["ID"].Head // =genename, this value is used to assign mRNAs of the same gene together
                group
                |> Seq.choose isRNA //the transcripts of the gene are chosen
                |> Seq.mapi (fun i element ->
                    // every transcript of gene gets its own number i and other info is collected from element and used for info of protein
                    let modelInfo = createProteinModelInfoFromEntry i locus element
                    let r = System.Text.RegularExpressions.Regex.Match(modelInfo.Id,regexPattern)
                    // the gff3 id has to be matched with the sequence in the fasta file. therefore the regexpattern is used
                    if r.Success then
                        r.Value,
                        modelInfo
                    else
                        failwithf "could not match gff3 entry id %s with regexpattern %s. Either gff3 file is corrupt or regexpattern is not correct" modelInfo.Id regexPattern
                )

            | _ -> Seq.empty
        )
        |> Seq.concat
        |> Map.ofSeq

    /// Creates a lookup data base to assign peptides to the proteins they are contained in
    let createPeptideProteinRelation (protModels:seq<ProteinModel<'id,'chromosomeId,'geneLocus,'sequence list> option>) =
        let ppRelation = BidirectionalDictionary<'sequence,ProteinModelInfo<'id,'chromosomeId,'geneLocus>>()
        protModels
        |> Seq.iter (fun prot ->
            // insert peptide-protein relationship
            // Todo: change type of proteinID in digest
            match prot with
            | Some proteinModel ->
                proteinModel.Sequence
                |> Seq.iter (fun pepSequence -> ppRelation.Add pepSequence proteinModel.ProteinModelInfo)
            | None -> ()
        )
        ppRelation

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

    // Creates a Map of peptides with their highest found score
    let createPeptideScoreMap (psmInputs: PSMInput list list) =
        psmInputs
        |> List.concat
        |> List.groupBy (fun psm -> psm.Seq)
        |> List.map (fun (sequence, psmList) ->
            // This sequence may contain modifications.
            // Depending on the type of lookup this map is used for, the modifications have to be removed.
            sequence,
            psmList
            |> List.maxBy (fun psm -> psm.Score)
            |> fun psm -> psm.Score
        )
        |> Map.ofList

    // Assigns a score to each protein with reverse digested peptides based on the peptides obtained in psm.
    let createReverseProteinScores (reverseProteins: (string*string[])[]) (peptideScoreMap: Map<string,float>) =
        // Remove modifications from map since protein inference was also done using unmodified peptides
        let scoreMapWithoutMods =
            peptideScoreMap
            |> Map.toArray
            |> Array.map (fun (seq, score) -> removeModification seq, score)
            |> Map.ofArray
        reverseProteins
        |> Array.map (fun (protein, peptides) ->
        protein,
            (
                peptides
                // looks wether the peptides resulting from the reverse digest appear in the peptides from the psm
                |> Array.map (fun pep ->
                    scoreMapWithoutMods.TryFind pep
                    |> (fun x ->
                        match x with
                        | Some score -> score
                        | None -> 0.
                    )
                )
                |> Array.sum,
                peptides
            )
        )
        // only hits are relevant
        |> Array.filter (fun (protein, (score, peptides)) -> score <> 0.)
        |> Map.ofArray

    /// Sums up score of all peptide sequences
    let assignPeptideScores (peptideSequences : string []) (peptideScoreMap : Map<string,float>) =
        peptideSequences
        |> Array.map (fun sequence -> peptideScoreMap.Item sequence)
        |> Array.sum

    /// Looks if the given protein accession is present in a map of identified decoy proteins and assigns its score when found.
    let assignDecoyScoreToTargetScore (proteins: string) (decoyScores: Map<string,(float * string[])>) =
        let prots = proteins |> String.split ';'
        prots
        |> Array.map (fun protein ->
            decoyScores.TryFind protein
            |> fun protOpt ->
                match protOpt with
                // peptides which pointed to the decoy version of this protein are discarded here, they can be included if needed
                | Some (score,peptides) -> score
                | None -> 0.
        )
        |> Array.max

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

    module MAYU =

        /// Returns proteins sorted into bins according to their size.
        /// proteins are the proteins which were found in either the reverse or forward proteininference, proteinsFromDB are the proteins with peptide sequence
        /// on which the inference was performed.
        let binProteinsLength (proteins: InferredProteinClassItemScored []) (proteinsFromDB: (string*string)[]) binCount =
            // need to bin all proteins in the db, not only those with hit?
            let proteinsNotMatchedDB =
                let proteinsMatched =
                    proteins
                    |> Array.collect (fun protein ->
                        protein.GroupOfProteinIDs
                        |> String.split ';'
                    )
                    |> Set.ofArray
                // Creates an entry for every protein that is present in the search db and wasn't inferred
                proteinsFromDB
                |> Array.choose (fun (proteinName, peptideSequence) ->
                    if Set.contains proteinName proteinsMatched then
                        None
                    else
                        Some ((createInferredProteinClassItemScored proteinName BioFSharp.PeptideClassification.PeptideEvidenceClass.Unknown [|peptideSequence|] (-1.) (-1.) false false false),
                             float peptideSequence.Length)
                )

            // Adds the length of the peptide sequence to every protein, since they should be binned according to it
            let addedSequenceLength =
                let proteinLengthMap =
                    proteinsFromDB
                    |> Array.map (fun (protein,sequence) -> protein, sequence.Length)
                    |> Map.ofArray
                proteins
                |> Array.map (fun ipcis ->
                    ipcis,
                    let groupOfProteins =
                        ipcis.GroupOfProteinIDs
                        |> String.split ';'
                    groupOfProteins
                    |> Array.map (fun protein ->
                        match proteinLengthMap.TryFind protein with
                        | None -> failwith "A protein where you are trying to get the length from isn't present in the database"
                        | Some length -> float length
                    )
                    // lengths are averaged for protein groups
                    |> Array.average
                )
            let combined = Array.append addedSequenceLength proteinsNotMatchedDB
            // sorting treats protein groups as one protein with average length. They are also treated as one protein for total and target counts.
            let sortLength =
                combined
                |> Array.sortBy snd
                |> Array.map fst
            let bins =
                let binSize =
                    ceil (float sortLength.Length / binCount)
                sortLength
                |> Array.chunkBySize (int binSize)
            bins

        // The original paper of Mayu describes a protein as:
        // FP = if only if all its peptides with q <= threshold are decoy
        // TP = at least one of its peptides with q <= threshold is target
        // However, the way mayu estimates it on the program is like this:
        // FP = any protein that contains a decoy psm with q <= threshold
        // TP = any protein that contains a target psm with q <= threshold
        // They do not consider protein containing both decoy and target psms.
        // ProteomIQon currently uses the picked target decoy approach with the following definitions:
        // FP = protein where the score of decoy hits is higher than score of target hits
        // TP = protein where the score of target hits is higher than score of decoy hits
        // Also, mayu estimates q as the empirical (target-decoy) q value.
        // Percolator estimates q as the empirical (target-decoy) q value and adjusted by pi0
        // Mayu extracts the list of TP and FP proteins from PSM level whereas percolator
        // extract the list of TP and FP proteins from peptide level, this avoids redundancy and
        // gives a better calibration since peptide level q values are re-adjusted in percolator.
        // ProteomIQon extracts TP and FP proteins from the result of the picked target decoy approach.
        // This creates sometimes a difference in the number of TP and FP proteins between percolator, Mayu, and ProteomIQon,
        // which causes a slight difference in the estimated protein FDR.

        /// Calculates the expected false positives for every protein bin and sums them up.
        let expectedFP (proteinBin: InferredProteinClassItemScored []) =
            let numberTarget =
                proteinBin
                |> Array.sumBy (fun protein ->
                    match not protein.DecoyHasBetterScore && protein.FoundInDB with
                    | true -> 1.
                    | false -> 0.
                )
            let numberDecoy =
                proteinBin
                 |> Array.sumBy (fun protein ->
                     match protein.DecoyHasBetterScore && protein.FoundInDB with
                     | true -> 1.
                     | false -> 0.
                 )
            let total =
                let notFound =
                    proteinBin
                    |> Array.sumBy (fun protein ->
                        match not protein.FoundInDB with
                        | true -> 1.
                        | false -> 0.
                    )
                notFound + numberTarget + numberDecoy
            // MAYU rounds the number of expected false positives for every bin to the first decimal
            let fpInBin = estimatePi0HG total numberTarget numberDecoy
            fpInBin

    
    /// Calculates the fdr of the data using the MAYU method. The proteinsFromDB is the DB that was used for the inference.
    let calculateFDRwithMAYU (data: InferredProteinClassItemScored []) (proteinsFromDB: (string*string)[]) =
        let proteinBins = MAYU.binProteinsLength data proteinsFromDB 10.
        let estimatedFP =
            proteinBins
            |> Array.fold (fun acc proteinBin -> acc + MAYU.expectedFP proteinBin) 0.
        let targetCount =
            data
            |> Array.sumBy (fun x ->
                match x.DecoyHasBetterScore with
                | true -> 0.
                | false -> 1.
            )
        let fdr =
            if (isNan estimatedFP) || (isInf estimatedFP) || estimatedFP = 0. then
                1.
            elif (estimatedFP / targetCount < 0.) || (estimatedFP / targetCount > 1.) then
                1.
            else
                estimatedFP / targetCount
        fdr

    /// Calculates Decoy/Target ratio
    let calculateFDRwithDecoyTargetRatio (data: InferredProteinClassItemScored []) =
        // Should decoy Hits be doubled?: Target-decoy search strategy for increasedconfidence in large-scale proteinidentifications by mass spectrometry
        let decoyCount  =
            data
            |> Array.sumBy (fun x ->
                match x.DecoyHasBetterScore with
                | true -> 1.
                | false -> 0.
            )
        let targetCount =
            data
            |> Array.sumBy (fun x ->
                match x.DecoyHasBetterScore with
                | true -> 0.
                | false -> 1.
            )
        decoyCount/targetCount

    /// Gives a function to calculate the q value for a score in a dataset using Lukas method and Levenberg Marguardt fitting
    let calculateQValueLogReg fdrEstimate bandwidth (data: 'a []) (isDecoy: 'a -> bool) (decoyScoreF: 'a -> float) (targetScoreF: 'a -> float) =
        // Input for q value calculation
        let createTargetDecoyInput =
            data
            |> Array.map (fun item ->
                if isDecoy item then
                    createQValueInput (decoyScoreF item) true
                else
                    createQValueInput (targetScoreF item) false
            )

        let scores,pep,qVal =
            binningFunction bandwidth fdrEstimate (fun (x: QValueInput) -> x.Score) (fun (x: QValueInput) -> x.IsDecoy) createTargetDecoyInput
            |> fun (scores,pep,qVal) -> scores.ToArray(), pep.ToArray(), qVal.ToArray()

        //Chart.Point (scores,qVal)
        //|> Chart.Show

        // gives a range of 1 to 30 for the steepness. This can be adjusted depending on the data, but normally it should lie in this range
        let initialGuess =
            Fitting.NonLinearRegression.LevenbergMarquardtConstrained.initialParamsOverRange scores qVal [|1. .. 30.|]
            |> Array.map (fun guess -> Table.lineSolverOptions guess)

        // performs Levenberg Marguardt Constrained algorithm on the data for every given initial estimate with different steepnesses and selects the one with the lowest RSS
        let estimate =
            initialGuess
            |> Array.map (fun initial ->
                if initial.InitialParamGuess.Length > 3 then failwith "Invalid initial param guess for Logistic Function"
                let lowerBound =
                    initial.InitialParamGuess
                    |> Array.map (fun param -> param - (abs param) * 0.1)
                    |> vector
                let upperBound =
                    initial.InitialParamGuess
                    |> Array.map (fun param -> param + (abs param) * 0.1)
                    |> vector
                LevenbergMarquardtConstrained.estimatedParamsWithRSS Table.LogisticFunctionDescending initial 0.001 10.0 lowerBound upperBound scores qVal
            )
            |> Array.filter (fun (param,rss) -> not (param |> Vector.exists System.Double.IsNaN))
            |> Array.minBy snd
            |> fst

        let logisticFunction = Table.LogisticFunctionDescending.GetFunctionValue estimate
        logisticFunction

    /// Gives a function to calculate the q value for a score in a dataset using Storeys method
    let calculateQValueStorey (data: 'a[]) (isDecoy: 'a -> bool) (decoyScoreF: 'a -> float) (targetScoreF: 'a -> float) =
        // Gives an array of scores with the frequency of decoy and target hits at that score
        let scoreFrequencies =
            data
            |> Array.map (fun x ->
                if isDecoy x then
                    decoyScoreF x, true
                else
                    targetScoreF x, false
            )
            // groups by score
            |> Array.groupBy fst
            // counts occurences of targets and decoys at that score
            |> Array.map (fun (score,scoreDecoyInfo) ->
                let decoyCount =
                    scoreDecoyInfo
                    |> Array.sumBy (fun (score, decoyInfo) ->
                        match decoyInfo with
                        | true -> 1.
                        | false -> 0.
                    )
                let targetCount =
                    scoreDecoyInfo
                    |> Array.sumBy (fun (score, decoyInfo) ->
                        match decoyInfo with
                        | true -> 0.
                        | false -> 1.
                    )
                createScoreTargetDecoyCount score decoyCount targetCount
            )
            |> Array.sortByDescending (fun x -> x.Score)

        // Goes through the list and assigns each protein a "q value" by dividing total decoy hits so far through total target hits so far
        let reverseQVal =
            scoreFrequencies
            |> Array.fold (fun (acc: (float*float*float*float) list) scoreCounts ->
                let _,_,decoyCount,targetCount = acc.Head
                // Should decoy hits be doubled?
                // accumulates decoy hits
                let newDecoyCount  = decoyCount + scoreCounts.DecoyCount(* * 2.*)
                // accumulates target hits
                let newTargetCount = targetCount + scoreCounts.TargetCount
                let newQVal =
                    let nominator =
                        if newTargetCount > 0. then newTargetCount
                        else 1.
                    newDecoyCount / nominator
                (scoreCounts.Score, newQVal, newDecoyCount, newTargetCount):: acc
            ) [0., 0., 0., 0.]
            // removes last part of the list which was the "empty" initial entry
            |> fun list -> list.[.. list.Length-2]
            |> List.map (fun (score, qVal, decoyC, targetC) -> score, qVal)

        //Assures monotonicity by going through the list from the bottom to top and assigning the previous q value if it is smaller than the current one
        let score, monotoneQVal =
            if reverseQVal.IsEmpty then
                failwith "Reverse qvalues in Storey calculation are empty"
            let head::tail = reverseQVal
            tail
            |> List.fold (fun (acc: (float*float) list) (score, newQValue) ->
                let _,qValue = acc.Head
                if newQValue > qValue then
                    (score, qValue)::acc
                else
                    (score, newQValue)::acc
            )[head]
            |> Array.ofList
            |> Array.sortBy fst
            |> Array.unzip
        // Linear Interpolation
        let linearSplineCoeff = Interpolation.LinearSpline.initInterpolateSorted score monotoneQVal
        // takes a score from the dataset and assigns it a q value
        let interpolation = Interpolation.LinearSpline.interpolate linearSplineCoeff
        interpolation

    // Assigns a q value to an InferredProteinClassItemScored
    let assignQValueToIPCIS (qValueF: float -> float) (item: InferredProteinClassItemScored) =
        if item.Decoy then
            createInferredProteinClassItemQValue item (qValueF item.DecoyScore)
        else
            createInferredProteinClassItemQValue item (qValueF item.TargetScore)
