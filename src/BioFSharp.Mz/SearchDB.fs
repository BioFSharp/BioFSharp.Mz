namespace BioFSharp.Mz


open BioFSharp
open BioFSharp.IO

open System
open System.Data
open FSharpAux
open AminoAcids 
open ModificationInfo
//open BioSequences

module SearchDB =
    
    open Digestion
    open System.Data.SQLite

    type SearchModType =
        | Minus
        | Plus 

    type SearchModSite =
        | Any      of ModificationInfo.ModLocation
        | Specific of AminoAcids.AminoAcid * ModificationInfo.ModLocation 
    
    type SearchModification = {
        Name        : string
        Accession   : string
        Description : string
        IsBiological: bool
        Composition : string
        Site        : SearchModSite list 
        MType       : SearchModType
        XModCode    : string
        }
    
    let createSearchModification name accession description isBiological composition site mType xModCode = {
        Name=name; Accession=accession; Description=description; IsBiological=isBiological ; Composition=composition; Site=site; MType=mType; XModCode=xModCode }
    
    type SearchInfoIsotopic = {
        Name        : string
        SourceEl    : Elements.Element
        TargetEl    : Elements.Element}

    let createSearchInfoIsotopic name sourceEl targetEl  = {
        Name=name; SourceEl=sourceEl; TargetEl=targetEl  }
                      
    let createIsotopicMod (infoIsotopic: SearchInfoIsotopic) = 
        createModification infoIsotopic.Name true ModLocation.Isotopic (fun f -> Formula.replaceElement f infoIsotopic.SourceEl infoIsotopic.TargetEl)
        
    let getModBy (smodi:SearchModification) = 
        let modLoc = 
              match smodi.Site.Head with
              | Any modLoc'           -> modLoc'
              | Specific (aa,modLoc') -> modLoc' 
        match smodi.MType with
        | SearchModType.Plus  -> createModificationWithAdd smodi.Name smodi.IsBiological modLoc smodi.Composition
        | SearchModType.Minus -> createModificationWithSubstract smodi.Name smodi.IsBiological modLoc smodi.Composition                             
    
    type MassMode = 
        | Average
        | Monoisotopic
        override this.ToString() = 
            match this with
            | Average -> "Average"
            | Monoisotopic -> "Monoisotopic"
    
    // fastaHeaderToName
    // Digestion filter params 
    type SearchDbParams = {
        // name of database i.e. Creinhardtii_236_protein_full_labeled
        Name            : string
        // path of db storage folder
        DbFolder            : string
        FastaPath           : string
        FastaHeaderToName   : string -> string
        Protease            : Protease
        MinMissedCleavages  : int
        MaxMissedCleavages  : int
        MaxMass             : float
        MinPepLength        : int
        MaxPepLength        : int
        // valid symbol name of isotopic label in label table i.e. #N15
        IsotopicMod         : SearchInfoIsotopic list 
        MassMode            : MassMode
        MassFunction        : IBioItem -> float  
        FixedMods           : SearchModification list            
        VariableMods        : SearchModification list
        VarModThreshold     : int
        }

    
    let createSearchDbParams name dbPath fastapath fastaHeaderToName protease minMissedCleavages maxMissedCleavages 
        maxmass minPepLength maxPepLength globalMod massMode massFunction fixedMods variableMods varModThreshold = {
         Name=name; 
         DbFolder=dbPath; 
         FastaPath=fastapath;
         FastaHeaderToName=fastaHeaderToName;  
         Protease=protease; 
         MinMissedCleavages=minMissedCleavages; 
         MaxMissedCleavages=maxMissedCleavages;
         MaxMass=maxmass;
         MinPepLength=minPepLength
         MaxPepLength=maxPepLength
         IsotopicMod=globalMod; 
         MassMode=massMode; 
         MassFunction=massFunction;
         FixedMods=List.sort fixedMods; 
         VariableMods=List.sort variableMods
         VarModThreshold=varModThreshold
         }

    let pepLengthLimitsBy maxExpCharge (mzAquisitionWindow: float*float) =
        let lowerborder = 4 
        let upperborder = ((float maxExpCharge / 2.)*(snd mzAquisitionWindow)) / 111. //TODO: (AminoAcids.monoisoMass AminoAcids.AminoAcid.Xaa)
        int lowerborder, int upperborder 

    let massFBy massMode = 
        match massMode with
            |Monoisotopic -> BioFSharp.BioItem.initMonoisoMassWithMemP
            |Average      -> BioFSharp.BioItem.initAverageMassWithMemP;
      
       

    ///needed as input if element of SearchModSite is of UnionCase | Any
    let listOfAA = [
        AminoAcid.Ala; 
        AminoAcid.Cys; 
        AminoAcid.Asp; 
        AminoAcid.Glu; 
        AminoAcid.Phe; 
        AminoAcid.Gly; 
        AminoAcid.His; 
        AminoAcid.Ile; 
        AminoAcid.Lys; 
        AminoAcid.Leu; 
        AminoAcid.Met; 
        AminoAcid.Asn; 
        AminoAcid.Pyl; 
        AminoAcid.Pro; 
        AminoAcid.Gln; 
        AminoAcid.Arg; 
        AminoAcid.Ser; 
        AminoAcid.Thr; 
        AminoAcid.Sel; 
        AminoAcid.Val; 
        AminoAcid.Trp; 
        AminoAcid.Tyr
        ]


    // Record for a peptide sequence and its precalculated mass calcualted mass
    type PeptideWithFeature<'a,'b> = {
        GlobalMod: int
        Sequence : 'a
        Feature  : 'b
    }

    /// Creates PeptideWithMass record
    let createPeptideWithFeature globalMod sequence feature = 
        {GlobalMod=globalMod; Sequence=sequence; Feature=feature}

    // Record for a peptide sequence and a container for all modified peptides
    type PeptideContainer = {
        PeptideId    : int    
        Sequence     : string
        MissCleavageStart : int
        MissCleavageEnd   : int
        MissCleavageCount : int      
        Container    : PeptideWithFeature<string,float> list 
    }

    let createPeptideContainer peptideId sequence  missCleavageStart missCleavageEnd missCleavageCount container =
        {PeptideId=peptideId; Sequence=sequence; MissCleavageStart=missCleavageStart; 
            MissCleavageEnd=missCleavageEnd; MissCleavageCount=missCleavageCount; Container=container;}


    type ProteinContainer = {
        ProteinId    : int
        DisplayID    : string
        Sequence     : string
        Container    : PeptideContainer list
    }

    let createProteinContainer proteinId displayID sequence container = {
        ProteinId=proteinId;DisplayID=displayID;Sequence=sequence;Container=container }

    type LookUpResult<'a when 'a :> IBioItem> = {
        ModSequenceID       : int
        PepSequenceID       : int        
        Mass                : float 
        RoundedMass         : int64
        StringSequence      : string
        BioSequence         : 'a list
        GlobalMod           : int                
    }

    let createLookUpResult modSequenceId pepSequenceId mass roundedMass stringSequence bioSequence globalMod =
        {ModSequenceID=modSequenceId ;PepSequenceID = pepSequenceId; Mass = mass;RoundedMass = roundedMass; StringSequence=stringSequence; BioSequence=bioSequence; GlobalMod=globalMod }
    

    module Db =

        type SqlErrorCode =
            | DbDataBaseNotFound
            | DbInternalLogicError
            | DbAccessDenied
            | DBGeneric of string * int
            | UnknownSqlException of SQLiteException
            | Unknown of Exception

        let sqlErrorCodeFromException (ex: Exception) =
            match ex with
            | :? SQLiteException  as ex ->
                    match ex.ErrorCode with
                    | 1 -> SqlErrorCode.DbDataBaseNotFound
                    | 2 -> SqlErrorCode.DbInternalLogicError
                    | 3 -> SqlErrorCode.DbAccessDenied
                    | _ -> SqlErrorCode.UnknownSqlException ex
            |  _ ->  SqlErrorCode.Unknown ex

        type SqlAction =
            | Create
            | Select
            | Insert
            | Delet
            | Update

        type PeptideLookUpError =
            | DbProtein of SqlAction * SqlErrorCode
            | DbProteinItemNotFound
            | DbCleavageIndex of SqlAction * SqlErrorCode 
            | DbCleavageIndexItemNotFound
            | DbPepSequence of SqlAction * SqlErrorCode
            | DbPepSequenceItemNotFound
            | DbModSequence of SqlAction * SqlErrorCode
            | DbModSequenceItemNotFound
            | DbSearchmodification of SqlAction * SqlErrorCode
            | DbSearchmodificationItemNotFound
            | DbSearchParams of SqlAction * SqlErrorCode
            | DbSearchParamsItemNotFound
            | DbInitialisation of SqlAction * SqlErrorCode
            | DbInitialisation_Database_Matching_The_Selected_Parameters_Already_Exists
            | DbInitialisation_Database_With_Identical_Name_But_Different_Parameters_Already_Exists


        ///  Prepared statements via Closure
        module  SQLiteQuery =
    
            open System.Data
            open Newtonsoft.Json
            //Create DB table statements
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /// Creates Table SearchDbParams
            let createTableSearchDbParams (cn:SQLiteConnection) = 
                let querystring = 
                    "CREATE TABLE SearchDbParams (ID INTEGER, 
                                                 Name TEXT NOT NULL, 
                                                 DbFolder TEXT NOT NULL, 
                                                 FastaPath  TEXT NOT NULL, 
                                                 Protease  TEXT NOT NULL, 
                                                 MinMissedCleavages INTEGER NOT NULL,
                                                 MaxMissedCleavages  INTEGER NOT NULL,
                                                 MaxMass  REAL NOT NULL, 
                                                 MinPepLength INTEGER NOT NULL,
                                                 MaxPepLength INTEGER NOT NULL,
                                                 IsotopicLabel  TEXT NOT NULL, 
                                                 MassMode  TEXT NOT NULL, 
                                                 FixedMods TEXT NOT NULL,   
                                                 VariableMods TEXT NOT NULL,   
                                                 VarModThreshold INTEGER NOT NULL, 
                                                 PRIMARY KEY (ID ASC)
                                                 )"
                let cmd  = new SQLiteCommand(querystring, cn)
                cmd.ExecuteNonQuery()

            /// Creates Table Protein
            let createTableProtein  (cn:SQLiteConnection) =
                let querystring = 
                    "CREATE TABLE Protein (ID  INTEGER,
                                           Accession  TEXT NOT NULL,
                                           Sequence  TEXT NOT NULL,
                                           PRIMARY KEY (ID ASC)
                                           )"
                let cmd = new SQLiteCommand(querystring, cn)
                cmd.ExecuteNonQuery()



            /// Creates Table CleavageIndex
            let createTableCleavageIndex  (cn:SQLiteConnection) = 
                let querystring = 
                    "CREATE TABLE CleavageIndex (ID INTEGER, 
                                                 ProteinID INTEGER NOT NULL, 
                                                 PepSequenceID INTEGER NOT NULL, 
                                                 CleavageStart  INTEGER NOT NULL, 
                                                 CleavageEnd  INTEGER NOT NULL, 
                                                 MissCleavages  INTEGER NOT NULL, 
                                                 PRIMARY KEY (ID ASC),
                                                 CONSTRAINT ProteinID FOREIGN KEY (ProteinID) REFERENCES Protein (ID),
                                                 CONSTRAINT PepSequenceID FOREIGN KEY (PepSequenceID) REFERENCES PepSequence (ID)
                                                 )"
                let cmd = new SQLiteCommand(querystring, cn)
                cmd.ExecuteNonQuery()



            /// Creates Table PepSequence
            let createTablePepSequence  (cn:SQLiteConnection) = 
                let querystring = 
                    "CREATE TABLE PepSequence (ID INTEGER,
                                               Sequence TEXT NOT NULL COLLATE NOCASE ,
                                               PRIMARY KEY (ID ASC),
                                               CONSTRAINT PepSequenceUnique UNIQUE (Sequence ASC) ON CONFLICT IGNORE
                                               )"
                let cmd = new SQLiteCommand(querystring, cn)
                cmd.ExecuteNonQuery()



            /// Creates Table ModSequence
            let createTableModSequence  (cn:SQLiteConnection) =
                let querystring = 
                    "CREATE TABLE ModSequence (ID	INTEGER,
	                                           PepSequenceID INTEGER NOT NULL,
	                                           RealMass REAL NOT NULL,
	                                           RoundedMass INTEGER NOT NULL,
	                                           Sequence TEXT NOT NULL,
	                                           GlobalMod INT NOT NULL,
	                                           PRIMARY KEY (ID),
	                                           FOREIGN KEY (PepSequenceID) REFERENCES PepSequence (ID) 
                                               )"
                let cmd = new SQLiteCommand(querystring, cn)
                cmd.ExecuteNonQuery()



            //Create Index Statements
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            let setMassIndexOnModSequence (cn:SQLiteConnection) = 
                let querystring = "CREATE INDEX RoundedMassIndex ON ModSequence (RoundedMass ASC) "
                let cmd = new SQLiteCommand(querystring, cn)    
                cmd.ExecuteNonQuery()


            let setSequenceIndexOnPepSequence (cn:SQLiteConnection) = 
                let querystring = "CREATE INDEX SequenceIndex ON PepSequence (Sequence ASC) "
                let cmd = new SQLiteCommand(querystring, cn)    
                cmd.ExecuteNonQuery()


            let setPepSequenceIDIndexOnCleavageIndex (cn:SQLiteConnection) = 
                let querystring = "CREATE INDEX PepSequenceIDIndex ON CleavageIndex (PepSequenceID) "
                let cmd = new SQLiteCommand(querystring, cn)    
                cmd.ExecuteNonQuery()

            ///
            let setIndexOnModSequenceAndGlobalMod (cn:SQLiteConnection) =
                let querystring = "CREATE INDEX IF NOT EXISTS SequenceAndGlobalModIndex ON ModSequence (Sequence,GlobalMod)"
                let cmd = new SQLiteCommand(querystring, cn)    
                cmd.ExecuteNonQuery()
            //let setSequenceIndexOnPepSequence (cn:SQLiteConnection) = 
            //    let querystring = "CREATE INDEX SequenceIndex ON PepSequence (Sequence ASC) "
            //    let cmd = new SQLiteCommand(querystring, cn)    
            //    try
            //        let exec = cmd.ExecuteNonQuery()
            //        if  exec < 1 then
            //            Either.succeed cn
            //        else 
            //            PeptideLookUpError.DbInitialisation (SqlAction.Create, DBGeneric ("INDEX SequenceIndex",exec)) 
            //            |> Either.fail
            //    with            
            //    | _ as ex -> 
            //            PeptideLookUpError.DbInitialisation (SqlAction.Create,sqlErrorCodeFromException ex) 
            //            |> Either.fail

            //Manipulate Pragma Statements
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            let pragmaSynchronousOFF (cn:SQLiteConnection) = 
                let querystring = "PRAGMA synchronous = 0 "
                let cmd = new SQLiteCommand(querystring, cn)
                // result equals number of affected rows
                cmd.ExecuteNonQuery()


         
            //Insert Statements
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /// Prepares statement to insert a Protein entry
            let prepareInsertProtein (cn:SQLiteConnection) (tr) =
                let querystring = "INSERT INTO Protein (ID, Accession, Sequence) VALUES (@id, @accession, @sequence)"
                let cmd = new SQLiteCommand(querystring, cn, tr)
                cmd.Parameters.Add("@id", Data.DbType.Int32) |> ignore
                cmd.Parameters.Add("@accession", Data.DbType.String) |> ignore
                cmd.Parameters.Add("@sequence", Data.DbType.String) |> ignore
    
                (fun (id:int32) (accession:string) (sequence:string)  ->  
                        cmd.Parameters.["@id"].Value        <- id
                        cmd.Parameters.["@accession"].Value <- accession
                        cmd.Parameters.["@sequence"].Value  <- sequence
                        // result equals number of affected rows
                        cmd.ExecuteNonQuery()
                        )

   
            /// Prepares statement to insert a CleavageIndex entry
            let prepareInsertCleavageIndex (cn:SQLiteConnection) (tr) =
                let querystring = "INSERT INTO CleavageIndex (ProteinID, 
                                                              PepSequenceID, 
                                                              CleavageStart, 
                                                              CleavageEnd, 
                                                              MissCleavages) 
                                                              VALUES (@proteinID, 
                                                                      @pepSequenceID, 
                                                                      @cleavageStart, 
                                                                      @cleavageEnd, 
                                                                      @missCleavages)"
                let cmd = new SQLiteCommand(querystring, cn, tr)
                //  cmd.Parameters.Add("@id", Data.DbType.Int64) |> ignore
                cmd.Parameters.Add("@proteinID", Data.DbType.Int32) |> ignore
                cmd.Parameters.Add("@pepSequenceID", Data.DbType.Int32) |> ignore
                cmd.Parameters.Add("@cleavageStart", Data.DbType.Int32) |> ignore
                cmd.Parameters.Add("@cleavageEnd", Data.DbType.Int32) |> ignore
                cmd.Parameters.Add("@missCleavages", Data.DbType.Int32) |> ignore
    
                (fun  (proteinID:int32) (pepSequenceID:int32) (cleavageStart:int32) (cleavageEnd:int32) (missCleavages:int32)  -> // (id:uint64)
                        // cmd.Parameters.["@id"].Value            <- id
                        cmd.Parameters.["@proteinID"].Value     <- proteinID
                        cmd.Parameters.["@pepSequenceID"].Value <- pepSequenceID
                        cmd.Parameters.["@cleavageStart"].Value <- cleavageStart
                        cmd.Parameters.["@cleavageEnd"].Value   <- cleavageEnd
                        cmd.Parameters.["@missCleavages"].Value <- missCleavages
                        // result equals number of affected rows
                        cmd.ExecuteNonQuery()
                        )


            /// Prepares statement to insert a PepSequence entry
            let prepareInsertPepSequence (cn:SQLiteConnection) (tr) =
                let querystring = "INSERT INTO PepSequence (ID, Sequence) VALUES (@id, @sequence)"
                let cmd = new SQLiteCommand(querystring, cn, tr)
                cmd.Parameters.Add("@id", Data.DbType.Int32) |> ignore
                cmd.Parameters.Add("@sequence", Data.DbType.String) |> ignore
    
                (fun (id:int32) (sequence:string)  ->
                        cmd.Parameters.["@id"].Value  <- id
                        cmd.Parameters.["@sequence"].Value  <- sequence
                        // result equals number of affected rows
                        cmd.ExecuteNonQuery()
                        )



            /// Prepares statement to insert a ModSequence entry
            let prepareInsertModSequence (cn:SQLiteConnection) (tr) =
                let querystring = "INSERT INTO ModSequence (PepSequenceID, 
                                                            RealMass, 
                                                            RoundedMass, 
                                                            Sequence, 
                                                            GlobalMod) 
                                                            VALUES (@pepSequenceID, 
                                                                    @realmass, 
                                                                    @roundedmass, 
                                                                    @sequence, 
                                                                    @globalMod)" //ID, @id
                let cmd = new SQLiteCommand(querystring, cn, tr)
                //cmd.Parameters.Add("@id", Data.DbType.Int64) |> ignore
                cmd.Parameters.Add("@pepSequenceID", Data.DbType.Int32) |> ignore
                cmd.Parameters.Add("@realmass", Data.DbType.Double) |> ignore
                cmd.Parameters.Add("@roundedmass", Data.DbType.Int64) |> ignore        
                cmd.Parameters.Add("@sequence", Data.DbType.String) |> ignore
                cmd.Parameters.Add("@globalMod", Data.DbType.Int32) |> ignore
    
                (fun  (pepSequenceID:int32) (realmass: float) (roundedmass: int64) (sequence: string) (globalMod:int)  -> //(id:uint64)
                        // cmd.Parameters.["@id"].Value            <- id            
                        cmd.Parameters.["@pepSequenceID"].Value     <- pepSequenceID
                        cmd.Parameters.["@realmass"].Value          <- realmass
                        cmd.Parameters.["@roundedmass"].Value       <- roundedmass
                        cmd.Parameters.["@sequence"].Value          <- sequence
                        cmd.Parameters.["@globalMod"].Value         <- globalMod
                        // result equals number of affected rows
                        cmd.ExecuteNonQuery()
                        )

            /// Prepares statement to insert a SearchDBParams entry
            let prepareInsertSearchDbParams (cn:SQLiteConnection) =
                let querystring = "INSERT INTO SearchDbParams (Name,
                                                               DbFolder, 
                                                               FastaPath,
                                                               Protease, 
                                                               MinMissedCleavages, 
                                                               MaxMissedCleavages,
                                                               MaxMass,
                                                               MinPepLength,
                                                               MaxPepLength, 
                                                               IsotopicLabel, 
                                                               MassMode, 
                                                               FixedMods, 
                                                               VariableMods,
                                                               VarModThreshold) 
                                                               VALUES (@name, 
                                                                       @dbFolder, 
                                                                       @fastaPath,
                                                                       @protease,
                                                                       @minMissedCleavages,
                                                                       @maxMissedCleavages, 
                                                                       @maxMass, 
                                                                       @minPepLength,
                                                                       @maxPepLength, 
                                                                       @isotopicLabel, 
                                                                       @massMode, 
                                                                       @fixedMods, 
                                                                       @variableMods,
                                                                       @varModThreshold)"
                let cmd = new SQLiteCommand(querystring, cn)
                cmd.Parameters.Add("@name", Data.DbType.String) |> ignore
                cmd.Parameters.Add("@dbFolder", Data.DbType.String) |> ignore
                cmd.Parameters.Add("@fastaPath", Data.DbType.String) |> ignore
                cmd.Parameters.Add("@protease", Data.DbType.String) |> ignore
                cmd.Parameters.Add("@minMissedCleavages", Data.DbType.Int32) |> ignore
                cmd.Parameters.Add("@maxMissedCleavages", Data.DbType.Int32) |> ignore
                cmd.Parameters.Add("@maxMass", Data.DbType.Double) |> ignore
                cmd.Parameters.Add("@minPepLength", Data.DbType.Int32) |> ignore
                cmd.Parameters.Add("@maxPepLength", Data.DbType.Int32) |> ignore
                cmd.Parameters.Add("@isotopicLabel", Data.DbType.String) |> ignore
                cmd.Parameters.Add("@massMode", Data.DbType.String) |> ignore
                cmd.Parameters.Add("@fixedMods", Data.DbType.String) |> ignore
                cmd.Parameters.Add("@variableMods", Data.DbType.String) |> ignore 
                cmd.Parameters.Add("@varModThreshold", Data.DbType.Int32) |> ignore 
                (fun (name:string) (dbFolder:string) (fastaPath:string) (protease:string) (minMissedCleavages:int32) (maxMissedCleavages:int32) (maxMass:double) (minPepLength:int32) (maxPepLength:int32) 
                    (isotopicLabel:string) (massMode:string) (fixedMods:string) (variableMods:string) (varModThreshold:int32)  ->  
                        cmd.Parameters.["@name"].Value                   <- name
                        cmd.Parameters.["@dbFolder"].Value               <- dbFolder
                        cmd.Parameters.["@fastaPath"].Value              <- fastaPath
                        cmd.Parameters.["@protease"].Value               <- protease
                        cmd.Parameters.["@minMissedCleavages"].Value     <- minMissedCleavages
                        cmd.Parameters.["@maxMissedCleavages"].Value     <- maxMissedCleavages
                        cmd.Parameters.["@maxMass"].Value                <- maxMass
                        cmd.Parameters.["@minPepLength"].Value           <- minPepLength
                        cmd.Parameters.["@maxPepLength"].Value           <- maxPepLength
                        cmd.Parameters.["@isotopicLabel"].Value          <- isotopicLabel
                        cmd.Parameters.["@massMode"].Value               <- massMode
                        cmd.Parameters.["@fixedMods"].Value              <- fixedMods
                        cmd.Parameters.["@variableMods"].Value           <- variableMods
                        cmd.Parameters.["@varModThreshold"].Value        <- varModThreshold
                        // result equals number of affected rows
                        cmd.ExecuteNonQuery()
                        )
                    
            //Select Statements
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /// Prepares statement to select all SearchDbParams 
            let selectSearchDbParams (cn:SQLiteConnection) =
                let querystring = "SELECT * FROM SearchDbParams"
                let cmd = new SQLiteCommand(querystring, cn)    
                use reader = cmd.ExecuteReader()            
                match reader.Read() with
                | true -> Some (reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetString(3), reader.GetString(4), 
                            reader.GetInt32(5), reader.GetInt32(6), reader.GetDouble(7),reader.GetInt32(8), reader.GetInt32(9), reader.GetString(10), reader.GetString(11), reader.GetString(12), reader.GetString(13),reader.GetInt32(14))
                | false -> None

            /// Prepares statement to select Accession and Sequence of Proteins from SearchDB
            let selectProteins (cn:SQLiteConnection) =
                let selectProteins =
                    let querystring = "SELECT Accession, Sequence FROM Protein"
                    let cmd = new SQLiteCommand(querystring, cn)
                    use reader = cmd.ExecuteReader()
                    (
                        [
                            while reader.Read() do
                                yield (reader.GetString(0), reader.GetString(1))
                        ]
                    )
                selectProteins
                
            /// Prepares statement to select all SearchDbParams entries by FastaPath, Protease, MinmissedCleavages, MaxmissedCleavages, MaxMass, MinPepLength, MaxPepLength, IsotopicLabel, MassMode, FixedMods, VariableMods, VarModsThreshold
            let prepareSelectSearchDbParamsbyParams (cn:SQLiteConnection) =
                let querystring = "SELECT * FROM SearchDbParams WHERE FastaPath=@fastaPath 
                                                                AND Protease=@protease 
                                                                AND MinMissedCleavages=@minMissedCleavages 
                                                                AND MaxMissedCleavages=@maxMissedCleavages
                                                                AND MaxMass=@maxMass
                                                                AND MinPepLength=@minPepLength
                                                                AND MaxPepLength=@maxPepLength
                                                                AND IsotopicLabel=@isotopicLabel 
                                                                AND MassMode=@massMode 
                                                                AND FixedMods=@fixedMods 
                                                                AND VariableMods=@variableMods
                                                                AND VarModThreshold=@varModThreshold"
                let cmd = new SQLiteCommand(querystring, cn) 
                cmd.Parameters.Add("@fastaPath", Data.DbType.String)         |> ignore
                cmd.Parameters.Add("@protease", Data.DbType.String)          |> ignore  
                cmd.Parameters.Add("@minMissedCleavages", Data.DbType.Int32) |> ignore  
                cmd.Parameters.Add("@maxMissedCleavages", Data.DbType.Int32) |> ignore  
                cmd.Parameters.Add("@maxMass", Data.DbType.Double)           |> ignore  
                cmd.Parameters.Add("@minPepLength", Data.DbType.Int32)       |> ignore  
                cmd.Parameters.Add("@maxPepLength", Data.DbType.Int32)       |> ignore  
                cmd.Parameters.Add("@isotopicLabel", Data.DbType.String)     |> ignore  
                cmd.Parameters.Add("@massMode", Data.DbType.String)          |> ignore  
                cmd.Parameters.Add("@fixedMods", Data.DbType.String)         |> ignore  
                cmd.Parameters.Add("@variableMods", Data.DbType.String)      |> ignore
                cmd.Parameters.Add("@varModThreshold", Data.DbType.Int32)    |> ignore       
                (fun (fastaPath:string) (protease:string) (minMissedCleavages:int32) (maxMissedCleavages:int32) (maxMass:float) 
                     (minPepLength:int32) (maxPepLength:int32) (isotopicLabel:string) (massMode:string) (fixedMods:string) (variableMods:string) (varModThreshold:int32) ->         
                    cmd.Parameters.["@fastaPath"].Value             <- fastaPath
                    cmd.Parameters.["@protease"].Value              <- protease
                    cmd.Parameters.["@minMissedCleavages"].Value    <- minMissedCleavages
                    cmd.Parameters.["@maxMissedCleavages"].Value    <- maxMissedCleavages
                    cmd.Parameters.["@maxMass"].Value               <- maxMass
                    cmd.Parameters.["@minPepLength"].Value          <- minPepLength
                    cmd.Parameters.["@maxPepLength"].Value          <- maxPepLength
                    cmd.Parameters.["@isotopicLabel"].Value         <- isotopicLabel
                    cmd.Parameters.["@massMode"].Value              <- massMode
                    cmd.Parameters.["@fixedMods"].Value             <- fixedMods
                    cmd.Parameters.["@variableMods"].Value          <- variableMods
                    cmd.Parameters.["@varModThreshold"].Value       <- varModThreshold
                    use reader = cmd.ExecuteReader()            
                    match reader.Read() with
                    | true -> Some (reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetString(3), reader.GetString(4), 
                                reader.GetInt32(5), reader.GetInt32(6), reader.GetDouble(7),reader.GetInt32(8), reader.GetInt32(9), reader.GetString(10), reader.GetString(11), reader.GetString(12), reader.GetString(13),reader.GetInt32(14))
                    | false -> None       
                )

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /// Prepares statement to select all SearchModifications by SearchDbParams
            let prepareSelectSearchModsbyParams (cn:SQLiteConnection) =
                let querystring = "SELECT * FROM SearchDbParams WHERE FastaPath=@fastaPath 
                                                                AND Protease=@protease 
                                                                AND MinMissedCleavages=@minMissedCleavages 
                                                                AND MaxMissedCleavages=@maxMissedCleavages
                                                                AND MaxMass=@maxMass
                                                                AND MinPepLength=@minPepLength
                                                                AND MaxPepLength=@maxPepLength
                                                                AND IsotopicLabel=@isotopicLabel 
                                                                AND MassMode=@massMode 
                                                                AND FixedMods=@fixedMods 
                                                                AND VariableMods=@variableMods
                                                                AND VarModThreshold=@varModThreshold"
                let cmd = new SQLiteCommand(querystring, cn) 
                cmd.Parameters.Add("@fastaPath", Data.DbType.String)         |> ignore
                cmd.Parameters.Add("@protease", Data.DbType.String)          |> ignore  
                cmd.Parameters.Add("@minMissedCleavages", Data.DbType.Int32) |> ignore  
                cmd.Parameters.Add("@maxMissedCleavages", Data.DbType.Int32) |> ignore  
                cmd.Parameters.Add("@maxMass", Data.DbType.Double)           |> ignore  
                cmd.Parameters.Add("@minPepLength", Data.DbType.Int32)       |> ignore  
                cmd.Parameters.Add("@maxPepLength", Data.DbType.Int32)       |> ignore  
                cmd.Parameters.Add("@isotopicLabel", Data.DbType.String)     |> ignore  
                cmd.Parameters.Add("@massMode", Data.DbType.String)          |> ignore  
                cmd.Parameters.Add("@fixedMods", Data.DbType.String)         |> ignore  
                cmd.Parameters.Add("@variableMods", Data.DbType.String)      |> ignore
                cmd.Parameters.Add("@varModThreshold", Data.DbType.Int32)    |> ignore          
                (fun (fastaPath:string) (protease:string) (minMissedCleavages:int32) (maxMissedCleavages:int32) (maxMass:float) 
                     (minPepLength:int32) (maxPepLength:int32) (isotopicLabel:string) (massMode:string) (fixedMods:string) (variableMods:string) (varModThreshold:int32) ->         
                    cmd.Parameters.["@fastaPath"].Value             <- fastaPath
                    cmd.Parameters.["@protease"].Value              <- protease
                    cmd.Parameters.["@minMissedCleavages"].Value    <- minMissedCleavages
                    cmd.Parameters.["@maxMissedCleavages"].Value    <- maxMissedCleavages
                    cmd.Parameters.["@maxMass"].Value               <- maxMass
                    cmd.Parameters.["@minPepLength"].Value          <- minPepLength
                    cmd.Parameters.["@maxPepLength"].Value          <- maxPepLength
                    cmd.Parameters.["@isotopicLabel"].Value         <- isotopicLabel
                    cmd.Parameters.["@massMode"].Value              <- massMode
                    cmd.Parameters.["@fixedMods"].Value             <- fixedMods
                    cmd.Parameters.["@variableMods"].Value          <- variableMods
                    cmd.Parameters.["@varModThreshold"].Value       <- varModThreshold   
                    use reader = cmd.ExecuteReader()            
                    match reader.Read() with
                    | true ->   (reader.GetInt32(0), reader.GetString(12), reader.GetString(13))
                                |> fun (id,fixMods,varMods) -> Newtonsoft.Json.JsonConvert.DeserializeObject<SearchModification list>(fixMods)@Newtonsoft.Json.JsonConvert.DeserializeObject<SearchModification list>(varMods) 
                    | false -> []
                                
                )
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            /// Prepares statement to select a Protein entry by ID        
            let prepareSelectProteinByID (cn:SQLiteConnection) (tr) =
                let querystring = "SELECT * FROM Protein WHERE ID=@id "
                let cmd = new SQLiteCommand(querystring, cn, tr) 
                cmd.Parameters.Add("@id", Data.DbType.Int32) |> ignore       
                (fun (id:int32)  ->         
                    cmd.Parameters.["@id"].Value <- id
                    use reader = cmd.ExecuteReader()            
                    match reader.Read() with
                    | true  -> (reader.GetInt32(0),reader.GetString(1), reader.GetString(2)) 
                    | false -> -1,"",""
                )   
                                    
            /// Prepares statement to select a Protein entry by Accession     
            let prepareSelectProteinByAccession (cn:SQLiteConnection) (tr) =
                let querystring = "SELECT * FROM Protein WHERE Accession=@accession"
                let cmd = new SQLiteCommand(querystring, cn, tr) 
                cmd.Parameters.Add("@accession", Data.DbType.String) |> ignore       
                (fun (accession:string)  ->         
                    cmd.Parameters.["@accession"].Value <- accession
                    use reader = cmd.ExecuteReader()            
                    match reader.Read() with
                    | true -> (reader.GetInt32(0),reader.GetString(1), reader.GetString(2))
                    | false -> -1,"",""
                )


            /// Prepares statement to select a Protein entry by Sequence     
            let prepareSelectProteinBySequence (cn:SQLiteConnection) (tr) =
                let querystring = "SELECT * FROM Protein WHERE Sequence=@sequence"
                let cmd = new SQLiteCommand(querystring, cn, tr) 
                cmd.Parameters.Add("@sequence", Data.DbType.String) |> ignore       
                (fun (sequence:string)  ->         
                    cmd.Parameters.["@sequence"].Value <- sequence

                    use reader = cmd.ExecuteReader()            
                    match reader.Read() with
                    | true -> (reader.GetInt32(0),reader.GetString(1), reader.GetString(2))
                    | false -> -1,"",""

                )

            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /// Prepares statement to select a CleavageIndex entry ProteinID 
            let prepareSelectCleavageIndexByProteinID (cn:SQLiteConnection) (tr) =
                let querystring = "SELECT * FROM CleavageIndex WHERE ProteinID=@proteinID"
                let cmd = new SQLiteCommand(querystring, cn, tr) 
                cmd.Parameters.Add("@proteinID", Data.DbType.Int64) |> ignore       
                (fun (proteinID:int32)  ->         
                    cmd.Parameters.["@proteinID"].Value <- proteinID     
                    use reader = cmd.ExecuteReader()            
                    match reader.Read() with
                    | true -> (reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetInt32(5)) 
                    | false -> -1,-1,-1,-1,-1,-1
                )

            /// Prepares statement to select a CleavageIndex entry PepSequenceID 
            let prepareSelectCleavageIndexByPepSequenceID  (cn:SQLiteConnection) (tr) =
                let querystring = "SELECT * FROM CleavageIndex WHERE PepSequenceID=@pepSequenceID"
                let cmd = new SQLiteCommand(querystring, cn, tr) 
                cmd.Parameters.Add("@pepSequenceID", Data.DbType.Int32) |> ignore
                let rec readerloop (reader:SQLiteDataReader) (acc:(int*int*int*int*int*int) list) =
                        match reader.Read() with 
                        | true  -> readerloop reader ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetInt32(5)) :: acc)
                        | false ->  acc 
                (fun (pepSequenceID:int32)  ->         
                    cmd.Parameters.["@pepSequenceID"].Value <- pepSequenceID
                    use reader = cmd.ExecuteReader()            
                    readerloop reader [] 
                )

            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /// Prepares statement to select a PepSequence entry by PepSequence 
            let prepareSelectPepSequenceBySequence' (cn:SQLiteConnection) (tr) =
                let querystring = "SELECT * FROM PepSequence WHERE Sequence=@sequence"
                let cmd = new SQLiteCommand(querystring, cn, tr) 
                cmd.Parameters.Add("@sequence", Data.DbType.String) |> ignore       
                (fun (sequence:string)  ->         
                    cmd.Parameters.["@sequence"].Value <- sequence
                    use reader = cmd.ExecuteReader()
                    match reader.Read() with 
                    | true ->  reader.GetInt32(0) 
                    | false -> -1
                )
            /// Prepares statement to select a PepSequence entry by PepSequence - Version without try.. with pattern to enhance the Select performance
            let prepareSelectPepSequenceBySequence (cn:SQLiteConnection) (tr) =
                let querystring = "SELECT * FROM PepSequence WHERE Sequence=@sequence"
                let cmd = new SQLiteCommand(querystring, cn, tr) 
                cmd.Parameters.Add("@sequence", Data.DbType.String) |> ignore       
                (fun (sequence:string)  ->         
                    cmd.Parameters.["@sequence"].Value <- sequence       
                    use reader = cmd.ExecuteReader()
                    reader.Read() |> ignore 
                    reader.GetInt32(0)           
                    )

            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /// Prepares statement to select a ModSequence entry by ModSequenceID
            let prepareSelectModsequenceByModSequenceID (cn:SQLiteConnection) tr =
                let querystring = "SELECT * FROM ModSequence WHERE ID=@id"
                let cmd = new SQLiteCommand(querystring, cn, tr) 
                cmd.Parameters.Add("@id", Data.DbType.Int32) |> ignore       
                (fun (id:int32)  ->        
                    cmd.Parameters.["@id"].Value <- id

                    use reader = cmd.ExecuteReader()            
                    match reader.Read() with
                    | true ->  (reader.GetInt32(0), reader.GetInt32(1),reader.GetDouble(2), reader.GetInt64(3), reader.GetString(4), reader.GetInt32(5))
                    | false -> -1,-1,nan,-1L,"",-1
                )

            /// Prepares statement to select a ModSequence entry by PepSequenceID
            let prepareSelectModsequenceByPepSequenceID (cn:SQLiteConnection) (tr) =
                let querystring = "SELECT * FROM ModSequence WHERE PepSequenceID=@pepSequenceID"
                let cmd = new SQLiteCommand(querystring, cn, tr) 
                cmd.Parameters.Add("@pepSequenceID", Data.DbType.Int32) |> ignore       
                (fun (pepSequenceID:int32)  ->        
                    cmd.Parameters.["@pepSequenceID"].Value <- pepSequenceID
                    use reader = cmd.ExecuteReader()            
                    match reader.Read() with
                    | true ->  (reader.GetInt32(0), reader.GetInt32(1),reader.GetDouble(2), reader.GetInt64(3), reader.GetString(4), reader.GetString(5)) 
                    | false -> -1,-1,nan,-1L,"",""

                )

            /// Prepares statement to select a ModSequence entry by Mass
            let prepareSelectModsequenceByMass (cn:SQLiteConnection) (tr) =
                let querystring = "SELECT * FROM ModSequence WHERE Mass=@mass"
                let cmd = new SQLiteCommand(querystring, cn, tr) 
                cmd.Parameters.Add("@mass", Data.DbType.Int32) |> ignore       
                (fun (mass: int)  ->        
                    cmd.Parameters.["@mass"].Value <- mass
                    use reader = cmd.ExecuteReader()            
                    match reader.Read() with
                    | true -> (reader.GetInt32(0), reader.GetInt32(1),reader.GetDouble(2), reader.GetInt64(3), reader.GetString(4), reader.GetString(5))  
                    | false -> -1,-1,nan,-1L,"",""
            )

            /// Prepares statement to select a ModSequence entry by Massrange (Between selected Mass -/+ the selected toleranceWidth)
            let prepareSelectModsequenceByMassRangeWith (cn:SQLiteConnection) tr =
                let querystring = "SELECT * FROM ModSequence WHERE RoundedMass BETWEEN @mass1 AND @mass2"
                let cmd = new SQLiteCommand(querystring, cn, tr) 
                cmd.Parameters.Add("@mass1", Data.DbType.Int64) |> ignore
                cmd.Parameters.Add("@mass2", Data.DbType.Int64) |> ignore
                let rec readerloop (reader:SQLiteDataReader) (acc:(int*int*float*int64*string*int) list) =
                        match reader.Read() with 
                        | true  -> readerloop reader (( reader.GetInt32(0), reader.GetInt32(1),reader.GetDouble(2), reader.GetInt64(3), reader.GetString(4), reader.GetInt32(5) ) :: acc)
                        | false ->  acc 
                fun (mass1:int64) (mass2:int64) ->
                cmd.Parameters.["@mass1"].Value <- mass1
                cmd.Parameters.["@mass2"].Value <- mass2
                use reader = cmd.ExecuteReader()            
                readerloop reader [] 

            /// Prepares statement to select a ModSequence entry by Massrange (Between selected Mass -/+ the selected toleranceWidth)
            let prepareSelectModsequenceByMassRange (cn:SQLiteConnection) =
                let querystring = "SELECT * FROM ModSequence WHERE RoundedMass BETWEEN @mass1 AND @mass2"
                let cmd = new SQLiteCommand(querystring, cn) 
                cmd.Parameters.Add("@mass1", Data.DbType.Int64) |> ignore
                cmd.Parameters.Add("@mass2", Data.DbType.Int64) |> ignore
                let rec readerloop (reader:SQLiteDataReader) (acc:(int*int*float*int64*string*int) list) =
                        match reader.Read() with 
                        | true  -> readerloop reader (( reader.GetInt32(0), reader.GetInt32(1),reader.GetDouble(2), reader.GetInt64(3), reader.GetString(4), reader.GetInt32(5) ) :: acc)
                        | false ->  acc 
                fun (mass1:int64) (mass2:int64) ->
                cmd.Parameters.["@mass1"].Value <- mass1
                cmd.Parameters.["@mass2"].Value <- mass2
                use reader = cmd.ExecuteReader()            
                readerloop reader [] 

            /// Prepares statement to select a ModSequence entry by Sequence
            let prepareSelectModsequenceBySequence (cn:SQLiteConnection) (tr) =
                let querystring = "SELECT * FROM ModSequence WHERE Sequence=@sequence"
                let cmd = new SQLiteCommand(querystring, cn, tr) 
                cmd.Parameters.Add("@sequence", Data.DbType.Double) |> ignore       
                (fun (sequence:string)  ->        
                    cmd.Parameters.["@sequence"].Value <- sequence         
                    use reader = cmd.ExecuteReader()            
                    match reader.Read() with
                    | true -> (reader.GetInt32(0), reader.GetInt32(1),reader.GetDouble(2), reader.GetInt64(3), reader.GetString(4), reader.GetString(5))  
                    | false -> -1,-1,nan,-1L,"",""
                )

            /// Prepares statement to select a ModSequence entry by Massrange (Between selected Mass -/+ the selected toleranceWidth)
            let prepareSelectMassByModSequenceAndGlobalMod (cn:SQLiteConnection) =
                let querystring = "SELECT RealMass FROM ModSequence WHERE Sequence=@sequence AND GlobalMod=@globalMod"
                let cmd = new SQLiteCommand(querystring, cn) 
                cmd.Parameters.Add("@sequence", Data.DbType.String) |> ignore
                cmd.Parameters.Add("@globalMod", Data.DbType.Int32) |> ignore
                (fun (sequence:string) (globalMod:int) ->
                cmd.Parameters.["@sequence"].Value  <- sequence
                cmd.Parameters.["@globalMod"].Value <- globalMod
                use reader = cmd.ExecuteReader()            
                match reader.Read() with 
                | true  -> Some (reader.GetDouble(0))
                | false -> Option.None)

            /// Prepares statement to select a Protein Accession entry by ID
            let prepareSelectProteinAccessionByID (cn:SQLiteConnection) (tr) =
                let querystring = "SELECT Accession FROM Protein WHERE ID=@id "
                let cmd = new SQLiteCommand(querystring, cn, tr)
                cmd.Parameters.Add("@id", DbType.Int32) |> ignore
                (fun (id:int32)  ->
                    cmd.Parameters.["@id"].Value <- id
                    use reader = cmd.ExecuteReader()
                    match reader.Read() with
                    | true  -> (reader.GetString(0))
                    | false -> ""
                )

            /// Prepares statement to select a Peptide Sequence entry by ID
            let prepareSelectPepSequenceByPepSequenceID (cn:SQLiteConnection) (tr) =
                let querystring = "SELECT Sequence FROM PepSequence WHERE ID=@pepSequenceID"
                let cmd = new SQLiteCommand(querystring, cn, tr)
                cmd.Parameters.Add("@pepSequenceID", DbType.Int32) |> ignore
                (fun (pepSequenceID:int)  ->
                    cmd.Parameters.["@pepSequenceID"].Value <- pepSequenceID
                    use reader = cmd.ExecuteReader()
                    reader.Read() |> ignore
                    reader.GetString(0)
                )

                //#define SQLITE_ERROR        1   /* SQL error or missing database */
                //#define SQLITE_INTERNAL     2   /* Internal logic error in SQLite */
                //#define SQLITE_PERM         3   /* Access permission denied */
                //#define SQLITE_ABORT        4   /* Callback routine requested an abort */
                //#define SQLITE_BUSY         5   /* The database file is locked */
                //#define SQLITE_LOCKED       6   /* A table in the database is locked */
                //#define SQLITE_NOMEM        7   /* A malloc() failed */
                //#define SQLITE_READONLY     8   /* Attempt to write a readonly database */
                //#define SQLITE_INTERRUPT    9   /* Operation terminated by sqlite3_interrupt()*/
                //#define SQLITE_IOERR       10   /* Some kind of disk I/O error occurred */
                //#define SQLITE_CORRUPT     11   /* The database disk image is malformed */
                //#define SQLITE_NOTFOUND    12   /* Unknown opcode in sqlite3_file_control() */
                //#define SQLITE_FULL        13   /* Insertion failed because database is full */
                //#define SQLITE_CANTOPEN    14   /* Unable to open the database file */
                //#define SQLITE_PROTOCOL    15   /* Database lock protocol error */
                //#define SQLITE_EMPTY       16   /* Database is empty */
                //#define SQLITE_SCHEMA      17   /* The database schema changed */
                //#define SQLITE_TOOBIG      18   /* String or BLOB exceeds size limit */
                //#define SQLITE_CONSTRAINT  19   /* Abort due to constraint violation */
                //#define SQLITE_MISMATCH    20   /* Data type mismatch */
                //#define SQLITE_MISUSE      21   /* Library used incorrectly */
                //#define SQLITE_NOLFS       22   /* Uses OS features not supported on host */
                //#define SQLITE_AUTH        23   /* Authorization denied */
                //#define SQLITE_FORMAT      24   /* Auxiliary database format error */
                //#define SQLITE_RANGE       25   /* 2nd parameter to sqlite3_bind out of range */
                //#define SQLITE_NOTADB      26   /* File opened that is not a database file */
                //#define SQLITE_NOTICE      27   /* Notifications from sqlite3_log() */
                //#define SQLITE_WARNING     28   /* Warnings from sqlite3_log() */
                //#define SQLITE_ROW         100  /* sqlite3_step() has another row ready */
                //#define SQLITE_DONE        101  /* sqlite3_step() has finished executing */





        /// Returns the database name given the SearchDbParams
        let getNameOf (sdbParams:SearchDbParams) =
            sprintf "%s\\%s.db" 
                (sdbParams.DbFolder |> FSharpAux.IO.PathFileName.normalizeFileName)
                (sdbParams.Name |> FSharpAux.IO.PathFileName.fileNameWithoutExtension )

        /// Returns a comma seperated string of given search modification list
        let getJsonStringOf item =
            Newtonsoft.Json.JsonConvert.SerializeObject(item)
                           
        /// Inserts SearchDbParams into DB
        let insertSdbParams cn (sdbParams:SearchDbParams) =         
            SQLiteQuery.prepareInsertSearchDbParams cn sdbParams.Name sdbParams.DbFolder sdbParams.FastaPath sdbParams.Protease.Name
                 sdbParams.MinMissedCleavages sdbParams.MaxMissedCleavages sdbParams.MaxMass sdbParams.MinPepLength sdbParams.MaxPepLength (getJsonStringOf sdbParams.IsotopicMod) (getJsonStringOf sdbParams.MassMode)
                (getJsonStringOf sdbParams.FixedMods) (getJsonStringOf sdbParams.VariableMods) sdbParams.VarModThreshold
        
        /// Select SearchDbParams entry from DB by given SearchDbParams
        let selectSdbParamsby cn (sdbParams:SearchDbParams) = 
            SQLiteQuery.prepareSelectSearchDbParamsbyParams cn sdbParams.FastaPath sdbParams.Protease.Name
             sdbParams.MinMissedCleavages sdbParams.MaxMissedCleavages sdbParams.MaxMass sdbParams.MinPepLength sdbParams.MaxPepLength (getJsonStringOf sdbParams.IsotopicMod) (getJsonStringOf sdbParams.MassMode)
                    (getJsonStringOf sdbParams.FixedMods) (getJsonStringOf sdbParams.VariableMods) sdbParams.VarModThreshold  
        

        /// 
        let xModToSearchModifications sMds = 
            match sMds with
            | [] -> 
                    sMds 
                    |> List.map (fun x -> x.XModCode , x)
                    |> Map.ofList 
            | _  -> Map.empty

        /// Builds xModToSearchMod Map from DB by given SearchDbParams
        let xModOfSearchDbParams cn (sdbParams:SearchDbParams) =       
            let allSMods =
                SQLiteQuery.prepareSelectSearchModsbyParams cn sdbParams.FastaPath sdbParams.Protease.Name sdbParams.MinMissedCleavages 
                    sdbParams.MaxMissedCleavages sdbParams.MaxMass sdbParams.MinPepLength sdbParams.MaxPepLength (getJsonStringOf sdbParams.IsotopicMod) (getJsonStringOf sdbParams.MassMode)
                    (getJsonStringOf sdbParams.FixedMods) (getJsonStringOf sdbParams.VariableMods) sdbParams.VarModThreshold 
            xModToSearchModifications allSMods

        /// Returns true if a db exists with the same parameter content
        let isExistsBy (sdbParams:SearchDbParams) =       
            let fileName = getNameOf sdbParams
            match FSharpAux.IO.FileIO.fileExists fileName with 
            | true  -> 
                let connectionString = sprintf "Data Source=%s;Version=3" fileName
                use cn = new SQLiteConnection(connectionString)
                cn.Open()
                match selectSdbParamsby cn sdbParams with
                | Some _   -> true
                | None   -> false                                                                           
            | false -> false

    
        /// Create a new file instance of the DB schema. Deletes already existing instance.
        let initDB fileName =
    
            let _ = FSharpAux.IO.FileIO.DeleteFile fileName 

            let connectionString = sprintf "Data Source=%s;Version=3" fileName
            use cn = new SQLiteConnection(connectionString)
  
            
            cn.Open()
            SQLiteQuery.createTableSearchDbParams       cn |> ignore
            SQLiteQuery.createTableProtein              cn |> ignore
            SQLiteQuery.createTableCleavageIndex        cn |> ignore
            SQLiteQuery.createTablePepSequence          cn |> ignore
            SQLiteQuery.createTableModSequence          cn |> ignore
            SQLiteQuery.setSequenceIndexOnPepSequence   cn |> ignore
            cn.Close()
        
            

        /// Bulk insert for a sequence of ProteinContainers
        let bulkInsert (cn:SQLiteConnection) (data:seq<ProteinContainer>) =
            cn.Open()
            let tr = cn.BeginTransaction()
            // Bind Insert and Select statements to connection and transaction / prepares statements
            let insertProtein = SQLiteQuery.prepareInsertProtein cn tr
            let insertCleavageIndex = SQLiteQuery.prepareInsertCleavageIndex cn tr
            let insertPepSequence = SQLiteQuery.prepareInsertPepSequence cn tr
            let selectPepSequenceBySequence = SQLiteQuery.prepareSelectPepSequenceBySequence cn tr
            let insertModSequence = SQLiteQuery.prepareInsertModSequence cn tr

            data
            |> Seq.iter 
                (fun (protContainer) ->
                    match insertProtein protContainer.ProteinId protContainer.DisplayID protContainer.Sequence with                                                                        
                    | 1 -> protContainer.Container
                           |> List.iter 
                                (fun pepContainer -> 
                                    match insertPepSequence pepContainer.PeptideId pepContainer.Sequence with                         
                                    | 1 ->  insertCleavageIndex protContainer.ProteinId  pepContainer.PeptideId pepContainer.MissCleavageStart
                                                pepContainer.MissCleavageEnd pepContainer.MissCleavageCount|>ignore
                                            pepContainer.Container
                                            |> List.iter (fun modPep -> 
                                                            insertModSequence (pepContainer.PeptideId) modPep.Feature ((Convert.ToInt64(modPep.Feature*1000000.)))
                                                                modPep.Sequence modPep.GlobalMod |> ignore
                                                            )  

                                    | _  -> insertCleavageIndex protContainer.ProteinId (selectPepSequenceBySequence pepContainer.Sequence) 
                                                pepContainer.MissCleavageStart pepContainer.MissCleavageEnd pepContainer.MissCleavageCount |>ignore                         
                                )   

                    | _ ->                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "Protein is already in the database" |>ignore
                )

            SQLiteQuery.setPepSequenceIDIndexOnCleavageIndex cn |> ignore 
            SQLiteQuery.setMassIndexOnModSequence cn |> ignore
            //Commit and dispose the transaction and close the SQLiteConnection
            tr.Commit()
            tr.Dispose()
            cn.Close()


    module ModCombinator =    
    
        open FSharpAux 
        open BioFSharp    
        open BioFSharp.AminoAcids
        open BioFSharp.ModificationInfo

        open Db

        /// Type abreviation
        type ModLookUpFunc = AminoAcid -> Modification list option

        type ModLookUp = {
            ResidualFixed                : ModLookUpFunc
            NTermAndResidualFixed        : ModLookUpFunc
            CTermAndResidualFixed        : ModLookUpFunc
            ProtNTermAndResidualFixed    : ModLookUpFunc
            ProtCTermAndResidualFixed    : ModLookUpFunc
            ResidualVariable             : ModLookUpFunc
            NTermAndResidualVariable     : ModLookUpFunc
            CTermAndResidualVariable     : ModLookUpFunc
            ProtNTermAndResidualVariable : ModLookUpFunc
            ProtCTermAndResidualVariable : ModLookUpFunc
            Total: string -> string option
            Isotopic : (GlobalModificationInfo.GlobalAAModificator*GlobalModificationInfo.GlobalModModificator) option
            }

        /// Flag indicates if potential modification is fixed
        [<Struct>]
        type AminoAcidWithFlag(flag:bool,amino:AminoAcid) =
            member this.IsFlaged  = flag
            member this.AminoAcid = amino
            new (amino) = AminoAcidWithFlag(false,amino)

        let isVarModified (a:AminoAcidWithFlag) =
            a.IsFlaged

        let setVarModifiedFlagOf (a:AminoAcid)  =
            AminoAcidWithFlag(true,a)

        let setFixedModifiedFlagOf (a:AminoAcid)  =
            AminoAcidWithFlag(false,a)
        
        /// Returns a list of all possible modified AminoAcids given the particular Searchmodification
        // param: aminoAcids is a list of all possible aminoacids
        let convertSearchModification (aminoAcids: AminoAcid list) (searchModification:SearchModification) =
            ///Creates AminoAcid Modification tuple. Concerns MType of the Searchmodification  
            let modificationOf (searchMod : SearchModification) modLocation =
                match searchMod.MType with
                | SearchModType.Plus  -> createModificationWithAdd searchMod.Name searchMod.IsBiological modLocation searchMod.Composition
                | SearchModType.Minus -> createModificationWithSubstract searchMod.Name searchMod.IsBiological modLocation searchMod.Composition
                       
            ///Matches elements of the SearchModification.Site with the SearchModSites Any or Specific; Returns tuple of AminoAcids and Modifications
            searchModification.Site 
            |> List.collect (fun site  ->                          
                                match site with
                                |Any(modLoc) -> 
                                    let tmpModification = modificationOf searchModification modLoc
                                    aminoAcids |> List.map (fun a -> (a,tmpModification))                                                                                                                                        
                                |Specific(aa, modLoc)->  
                                    let tmpModification = modificationOf searchModification modLoc
                                    [(aa,tmpModification)]
                            ) 
        

        /// Returns the ModLookup according to given SearchDbParams
        let modLookUpOf (dbParams:SearchDbParams) =         
        
            //Filters list of AminoAcid Modification tuples depending on the used logexp    
            let createAndFilterBy (logexp: AminoAcid*Modification -> bool) searchModifications = 
                searchModifications
                |> List.collect (convertSearchModification listOfAA)
                |> List.filter (logexp: _*Modification -> bool)
                |> Map.compose

            ///Logexp that returns true if ModLocation of modified AminoAcid equals Residual 
            let residual (_,modi) =  ModLocation.Residual.Equals(modi.Location)
            ///Logexp that returns true if ModLocation of modified AminoAcid equals Residual, Nterm or ProteinNterm 
            let nTermAndResidual (_,modi) = 
                ModLocation.Nterm.Equals(modi.Location) (*|| ModLocation.ProteinNterm.Equals(modi.Location)*) || ModLocation.Residual.Equals(modi.Location)
            ///Logexp that returns true if ModLocation of modified AminoAcid equals Residual, Cterm or ProteinCterm  
            let cTermAndResidual (_,modi) = 
                ModLocation.Cterm.Equals(modi.Location) (*|| ModLocation.ProteinCterm.Equals(modi.Location)*) || ModLocation.Residual.Equals(modi.Location)
            ///Logexp that returns true if ModLocation of modified AminoAcid equals Residual, Nterm or ProteinNterm 
            let protNTermAndResidual (_,modi) = 
                ModLocation.ProteinNterm.Equals(modi.Location) //|| ModLocation.Residual.Equals(modi.Location)
            ///Logexp that returns true if ModLocation of modified AminoAcid equals Residual, Cterm or ProteinCterm  
            let protCTermAndResidual (_,modi) = 
                ModLocation.ProteinCterm.Equals(modi.Location) //|| ModLocation.Residual.Equals(modi.Location)
                
            {                                                      
                ResidualFixed=
                    let lookUpR = createAndFilterBy residual dbParams.FixedMods
                    fun aa -> Map.tryFind aa lookUpR     
                NTermAndResidualFixed=
                    let lookUpNR = createAndFilterBy nTermAndResidual dbParams.FixedMods
                    fun aa -> Map.tryFind aa lookUpNR      
                CTermAndResidualFixed=
                    let lookUpCR = createAndFilterBy cTermAndResidual dbParams.FixedMods
                    fun aa -> Map.tryFind aa lookUpCR   
                ProtNTermAndResidualFixed=
                    let lookUpNR = createAndFilterBy protNTermAndResidual dbParams.FixedMods
                    fun aa -> Map.tryFind aa lookUpNR      
                ProtCTermAndResidualFixed=
                    let lookUpCR = createAndFilterBy protCTermAndResidual dbParams.FixedMods
                    fun aa -> Map.tryFind aa lookUpCR                         
                ResidualVariable =
                    let lookUpR = createAndFilterBy residual dbParams.VariableMods
                    fun aa -> Map.tryFind aa lookUpR     
                NTermAndResidualVariable =
                    let lookUpNR = createAndFilterBy nTermAndResidual dbParams.VariableMods
                    fun aa -> Map.tryFind aa lookUpNR      
                CTermAndResidualVariable =
                    let lookUpCR = createAndFilterBy cTermAndResidual dbParams.VariableMods
                    fun aa -> Map.tryFind aa lookUpCR   
                ProtNTermAndResidualVariable =
                    let lookUpNR = createAndFilterBy protNTermAndResidual dbParams.VariableMods
                    fun aa -> Map.tryFind aa lookUpNR      
                ProtCTermAndResidualVariable =
                    let lookUpCR = createAndFilterBy protCTermAndResidual dbParams.VariableMods
                    fun aa -> Map.tryFind aa lookUpCR   
                Total =
                    let lookupT = 
                        (dbParams.FixedMods@dbParams.VariableMods)
                        |> List.map (fun searchmod -> searchmod.Name, searchmod.XModCode)  
                        |> Map.ofList
                    fun aa -> Map.tryFind aa lookupT
                Isotopic = 
                    if dbParams.IsotopicMod <> [] then
                        let modiL = 
                            dbParams.IsotopicMod 
                            |> List.map createIsotopicMod 
                        Some (GlobalModificationInfo.initGlobalModificationDeltaOfAA dbParams.MassFunction modiL,GlobalModificationInfo.initGlobalModificationDeltaOfMod dbParams.MassFunction modiL )
                    else None
             }




        ///Returns modified or unmodified AminoAcid depending on the matching expression in a AminoAcidWithFlag struct
        ///The boolean value "false" is used to state that the Modification is fixed    
        let setFixModByLookUp (modLookUpFunc:ModLookUpFunc) (aa: AminoAcidWithFlag) =
            match modLookUpFunc aa.AminoAcid with
            | Some modiList -> setFixedModifiedFlagOf (AminoAcids.setModifications modiList aa.AminoAcid)
            | None -> aa     


        let mapModLocation loc = 
            match loc with 
            | ModLocation.Residual      -> ModLocation.Residual 
            | ModLocation.Nterm         -> ModLocation.Nterm
            | ModLocation.ProteinNterm  -> ModLocation.Nterm
            | ModLocation.Cterm         -> ModLocation.Cterm
            | ModLocation.ProteinCterm  -> ModLocation.Cterm
            | ModLocation.Isotopic      -> ModLocation.Isotopic
            

        ///Returns modified or unmodified AminoAcid depending on the matching expression in a AminoAcidWithFlag struct
        ///The boolean value "false" is used to state that the Modification is fixed    
        let setVarModByLookUp (modLookUpFunc:ModLookUpFunc) (aa: AminoAcidWithFlag) =
            match aa.AminoAcid with
            | Mod(aa',presentMods) -> 
                let occupiedModLocs =
                    presentMods
                    |> List.map (fun modi -> mapModLocation modi.Location )
                    |> Set.ofList
                match modLookUpFunc aa' with 
                | Some modiList -> 
                    let modsAtDifferentLocsCombined = 
                        modiList 
                        |> List.groupBy (fun x -> mapModLocation x.Location) 
                        |> List.map snd                        
                        |> FSharpAux.List.powerSetOf 
                        |> List.map FSharpAux.List.drawExaustively 
                        |> List.concat
                    let tmp = 
                        modsAtDifferentLocsCombined    
                        |> List.choose (fun mods -> 
                                            match mods |> List.exists (fun modi -> occupiedModLocs.Contains(mapModLocation modi.Location)) with 
                                            | true  -> None
                                            | false -> Some (setVarModifiedFlagOf (AminoAcids.setModifications mods aa.AminoAcid))
                                       )               
                    aa::tmp
                | None -> 
                    [aa]
            | _                  ->
                match modLookUpFunc aa.AminoAcid with 
                | Some modiList -> 
                    let modsAtDifferentLocsCombined = 
                        modiList 
                        |> List.groupBy (fun x -> mapModLocation x.Location) 
                        |> List.map snd                        
                        |> FSharpAux.List.powerSetOf 
                        |> List.map FSharpAux.List.drawExaustively  
                        |> List.concat
                    let tmp = 
                        modsAtDifferentLocsCombined    
                        |> List.map (fun mods -> setVarModifiedFlagOf (AminoAcids.setModifications mods aa.AminoAcid))
                    aa::tmp
                                                   
                | None -> 
                    [aa]
                                                              
        
        //let addProteinTerminalModifications minPos maxPos (modLookUp:ModLookUp) (digPep: DigestedPeptide<'a>) =    
        //    ///
        //    let modifyNTerminal (aal: AminoAcid list) = 
        //        match aal with 
        //        | []      -> []
        //        | h::tail -> 
        //            AminoAcidWithFlag h 
        //            |> setFixModByLookUp modLookUp.ProtNTermAndResidualFixed 
        //            |> setVarModByLookUp modLookUp.ProtNTermAndResidualVariable  
        //            |> List.map (fun item -> (item.AminoAcid)::tail)                              
        //    ///
        //    let modifyCTerminal (aal: AminoAcid list) = 
        //        let rec loop acc (aal: AminoAcid list) =
        //            match aal with 
        //            | h::[]   -> 
        //                AminoAcidWithFlag h 
        //                |> setFixModByLookUp modLookUp.ProtCTermAndResidualFixed 
        //                |> setVarModByLookUp modLookUp.ProtCTermAndResidualVariable
        //                |> List.map (fun item -> ((item.AminoAcid)::acc) |> List.rev )
        //            | h::tail -> loop (h::acc) tail

        //        loop [] aal             

        //    if digPep.CleavageStart = minPos && digPep.CleavageEnd = maxPos then
        //        let modifiedPeps = 
        //            digPep.PepSequence
        //            |> modifyNTerminal
        //            |> List.collect modifyCTerminal
        //            |> List.map (fun modS -> {digPep with PepSequence = modS})
        //        modifiedPeps
        //    elif digPep.CleavageStart = minPos then
        //        let modifiedPeps = 
        //            digPep.PepSequence
        //            |> modifyNTerminal
        //            |> List.map (fun modS -> {digPep with PepSequence = modS})
        //        modifiedPeps
        //    elif digPep.CleavageEnd = maxPos then
        //        let modifiedPeps = 
        //            digPep.PepSequence
        //            |> modifyCTerminal
        //            |> List.map (fun modS -> {digPep with PepSequence = modS})
        //        modifiedPeps 
        //    else
        //        [digPep]

        let addProteinTerminalModifications minPos maxPos (modLookUp:ModLookUp) (digPep: DigestedPeptide<'a>) =    
            ///
            let modifyNTerminal (aal: AminoAcid list) = 
                match aal with 
                | []      -> []
                | h::tail -> 
                    AminoAcidWithFlag h 
                    |> setFixModByLookUp modLookUp.ProtNTermAndResidualFixed 
                    |> setVarModByLookUp modLookUp.ProtNTermAndResidualVariable  
                    |> List.map (fun item -> (item.AminoAcid)::tail)                              
            ///
            let modifyCTerminal (aal: AminoAcid list) = 
                let rec loop acc (aal: AminoAcid list) =
                    match aal with 
                    | h::[]   -> 
                        AminoAcidWithFlag h 
                        |> setFixModByLookUp modLookUp.ProtCTermAndResidualFixed 
                        |> setVarModByLookUp modLookUp.ProtCTermAndResidualVariable
                        |> List.map (fun item -> ((item.AminoAcid)::acc) |> List.rev )
                    | h::tail -> loop (h::acc) tail

                loop [] aal             

            let modify peptide = 
                if digPep.CleavageStart = minPos && digPep.CleavageEnd = maxPos then
                    let modifiedPeps = 
                        peptide.PepSequence
                        |> modifyNTerminal
                        |> List.collect modifyCTerminal
                        |> List.map (fun modS -> {peptide with PepSequence = modS})
                    modifiedPeps
                elif digPep.CleavageStart = minPos then
                    let modifiedPeps = 
                        peptide.PepSequence
                        |> modifyNTerminal
                        |> List.map (fun modS -> {peptide with PepSequence = modS})
                    modifiedPeps
                elif digPep.CleavageEnd = maxPos then
                    let modifiedPeps = 
                        peptide.PepSequence
                        |> modifyCTerminal
                        |> List.map (fun modS -> {peptide with PepSequence = modS})
                    modifiedPeps 
                else
                    [peptide]

            match digPep.CleavageStart,digPep.PepSequence with 
            | 0, AminoAcid.Met::b::tail -> 
                (modify digPep)@(modify {digPep with CleavageStart=1; PepSequence = b::tail})     
            | _ -> 
                (modify digPep)

        /// Returns a list of all possible modified petide sequences and its masses according to the given modification-lookUp.
        /// It uses the given bioitem -> mass function and a function to aggregate the sequence.
        let combine (modLookUp:ModLookUp) threshold (massfunction:IBioItem -> float) seqfunction state (aal: AminoAcid list) =
            let rec loop c (modAcc:Modification list) seqAcc (aal: AminoAcid list) =
                match c with
                | c when c = threshold -> 
                    match aal with
                    | h::tail -> 
                        match tail with 
                        // Last amino acid => set NTermAndResidualFixed
                        | [] -> setFixModByLookUp modLookUp.CTermAndResidualFixed (AminoAcidWithFlag h) 
                        // all other amino acid => set ResidualFixed
                        | _  -> setFixModByLookUp modLookUp.ResidualFixed  (AminoAcidWithFlag h)
                        |> (fun item ->
                                match item.AminoAcid with
                                | Mod (a,m) -> loop c (m@modAcc) (seqfunction seqAcc item.AminoAcid) tail
                                | a         -> loop c modAcc  (seqfunction seqAcc item.AminoAcid) tail  
                            )
                
                    | [] -> [createPeptideWithFeature 0 seqAcc modAcc]
                | c -> 
                    match aal with
                    | h::tail -> 
                        match tail with
                        | [] -> 
                            AminoAcidWithFlag h 
                            |> setFixModByLookUp modLookUp.CTermAndResidualFixed 
                            |> setVarModByLookUp modLookUp.CTermAndResidualVariable 
                        | _  -> 
                            AminoAcidWithFlag h
                            |> setFixModByLookUp modLookUp.ResidualFixed  
                            |> setVarModByLookUp modLookUp.ResidualVariable
                        |> List.collect (fun item ->
                                                
                                                match (isVarModified item),item.AminoAcid with
                                                | true,Mod (a,m)   -> loop (c+1) (m@modAcc) (seqfunction seqAcc item.AminoAcid) tail                                       
                                                | false, Mod (a,m) -> loop c (m@modAcc) (seqfunction seqAcc item.AminoAcid) tail
                                                | false, a         -> loop c modAcc (seqfunction seqAcc item.AminoAcid) tail 
                                                | true,_ -> failwith "Matching case impossible: Check AminoAcidWithFlag"
                                            )                   
                                                   
                    | [] ->  [createPeptideWithFeature 0 seqAcc modAcc]   
            let massOfPeptide =  
                    aal
                    |> List.fold (fun acc x -> 
                                    match x with
                                    | Mod(a,_) ->
                                        acc + massfunction a
                                    | a -> 
                                        acc + massfunction a
                                 ) (massfunction ModificationInfo.Table.H2O) //add water
            let sequencesNoGlobalMod = 
                match aal with
                | [] -> [createPeptideWithFeature 0 "" []]   
                | h::tail -> 
                        AminoAcidWithFlag h 
                        |> setFixModByLookUp modLookUp.NTermAndResidualFixed 
                        |> setVarModByLookUp modLookUp.NTermAndResidualVariable 
                        |> List.collect (fun item ->
                                            match (isVarModified item),item.AminoAcid with
                                            | true,Mod (a,m)   -> 
                                                loop (1) (m) (seqfunction "" item.AminoAcid) tail                                       
                                            | false, Mod (a,m) -> 
                                                loop 0 (m) (seqfunction "" item.AminoAcid) tail
                                            | false, a         -> 
                                                loop 0 [] (seqfunction "" item.AminoAcid) tail 
                                            | true,_ -> failwith "Matching case impossible: Check AminoAcidWithFlag"
                                        ) 
            match modLookUp.Isotopic with
            | None   -> 
                sequencesNoGlobalMod
                |> List.map (fun x -> 
                                let mass = x.Feature |> List.fold (fun s x -> s + (massfunction x)) massOfPeptide
                                createPeptideWithFeature x.GlobalMod x.Sequence mass
                            )
            | Some (globalAAModifier,globalModModifier) ->  
                let massWithGlobalMod =     
                    aal
                    |> List.fold (fun acc x -> 
                                    match x with
                                    | Mod(a,_) ->
                                        acc + massfunction a + (globalAAModifier a)
                                    | a -> 
                                        acc + massfunction a + (globalAAModifier a)
                                 ) (massfunction ModificationInfo.Table.H2O) 
                let sequencesWithMass =
                    sequencesNoGlobalMod
                    |> List.map (fun x -> 
                                    let mass = x.Feature |> List.fold (fun s x -> s + (massfunction x)) massOfPeptide
                                    createPeptideWithFeature x.GlobalMod x.Sequence mass
                                )    
                let sequencesGlobalModWithMass =
                    sequencesNoGlobalMod
                    |> List.map (fun x -> 
                                    let mass = 
                                        x.Feature 
                                        |> List.fold (fun s x -> 
                                                        s + (massfunction x) + (globalModModifier x) 
                                                      ) massWithGlobalMod
                                    createPeptideWithFeature 1 x.Sequence mass
                                )
                List.append sequencesWithMass sequencesGlobalModWithMass

        /// Returns a ModString representation. 
        let ToModStringBy (xModLookUp:string->string option) (aa: AminoAcid) =
            match (aa:AminoAcid) with
            | AminoAcids.Mod (aa, mds) ->  
                                    mds
                                    |> List.fold (fun acc (md:Modification) -> match xModLookUp md.Name with
                                                                                | Some x -> "[" + x + "]" + acc 
                                                                                | None -> failwithf "Failed"                                                                                         
                                                                                ) "" 
                                    |> (fun x -> x + ((BioItem.symbol aa).ToString()))
                                                     
            | aa -> ((BioItem.symbol aa).ToString())  


        /// Returns a list of all possible modified petide sequences and its masses according to the given modification-lookUp.
        /// The peptide sequence representation is ModString.
        let combineToModString (modLookUp:ModLookUp) threshold (massfunction:IBioItem -> float) (aal: AminoAcid list) =
            let toString = ToModStringBy modLookUp.Total
            let seqfunction = (fun state amino -> state + (toString amino))
            combine modLookUp threshold massfunction seqfunction "" aal
    
        ///  
        let getModSequencesOf (peptideId:ref<int>) (modLookUp:ModLookUp) threshold (massfunction:IBioItem -> float) (digPeps:DigestedPeptide<'a> []) =
            if Array.isEmpty digPeps then [||]
            else
            let minPos = 0 
            let maxPos = (digPeps |> Array.maxBy (fun x -> x.CleavageEnd)).CleavageEnd
            digPeps
            |> Array.map (fun pep ->
                                let container = 
                                    pep
                                    |> addProteinTerminalModifications minPos maxPos modLookUp
                                    |> List.collect (fun pep -> combineToModString modLookUp threshold massfunction pep.PepSequence)
                                peptideId.Value <- peptideId.Value + 1 
                                createPeptideContainer (peptideId.Value) (BioList.toString pep.PepSequence) pep.CleavageStart pep.CleavageEnd pep.MissCleavages container
                            )
   
    // --------------------------------------------------------------------------------
    // PeptideLookUp continues
    
    /// Creates a mapping from the XModCode of a modification to the modification itself.
    let xModToMods (sMds:SearchModification list) = 
        List.map (fun md -> md.XModCode, getModBy md) sMds
        |> Map.ofList

    /// Generates amino acid sequence of one-letter-code string containing modified AminoAcids indicated by 2 lowercase digits per modification. 
    let initOfModAminoAcidString (isotopMods: SearchInfoIsotopic list) (fixedAndVarMods: SearchModification list) (globalMod:int) (aaStr: string) : BioFSharp.BioList.BioList<AminoAcid>  =
        let isotopMods = List.map createIsotopicMod isotopMods
        let xModToMod  = xModToMods fixedAndVarMods
        ///
        let setGlobalMods globalMod isotopMods aaAcc = 
            if globalMod = 0 then aaAcc
            else AminoAcids.setModifications isotopMods aaAcc
        ///
        let rec accumulateMods count stringAcc aaAcc (aaStr: string) = 
            if count < 0 then 
                aaAcc, count
            else
            let currentC = aaStr.[count]
            if Char.IsUpper currentC then 
                aaAcc, count
            elif currentC = ']' then
                accumulateMods (count-1) stringAcc aaAcc aaStr 
            elif currentC = '[' then
                match Map.tryFind stringAcc xModToMod with
                | Some md -> 
                    let tmp = AminoAcids.setModification md aaAcc 
                    accumulateMods (count-1) "" tmp aaStr
                | None ->
                    accumulateMods (count-1) "" aaAcc aaStr
            else
                accumulateMods (count-1) (currentC.ToString() + stringAcc) aaAcc aaStr
        ///
        let rec loop count (aaAcc:AminoAcid option) seqAcc (aaStr: string) = 
            match count with
            | c when c < 0 ->
                match aaAcc with
                | Some aaAcc -> 
                    let tmp = setGlobalMods globalMod isotopMods aaAcc
                    tmp::seqAcc
                | None -> 
                    seqAcc
            | _ -> 
                let currentC = aaStr.[count]
                match Char.IsUpper currentC with
                | true ->
                    let aa = BioItemsConverter.OptionConverter.charToOptionAminoAcid currentC 
                    match aaAcc with
                    | Some aaAcc -> 
                        let tmp = setGlobalMods globalMod isotopMods aaAcc 
                        loop (count-1) aa (tmp::seqAcc) aaStr            
                    | None ->
                        loop (count-1) aa (seqAcc) aaStr            
                | false -> 
                    match aaAcc with
                    | None    ->
                        loop (count-1) None seqAcc aaStr
                    | Some aa -> 
                        let mdAA, count = accumulateMods (count-1) "" aa aaStr
                        let tmp = setGlobalMods globalMod isotopMods mdAA
                        loop count None (tmp::seqAcc) aaStr
        loop (aaStr.Length-1) None [] aaStr    
        
    /// Creates a LookUpResult out of a entry in the ModSequence table
    let createLookUpResultBy parseAAString (modSID,pepID,realMass,roundMass,seqs,gMod)  =
            let bSeq = parseAAString gMod seqs
            createLookUpResult modSID pepID realMass roundMass seqs bSeq gMod   

    /// Opens a new connection to the database in given path
    let getDBConnection dbFilePath =
        let connectionString = sprintf "Data Source=%s;Version=3" dbFilePath
        let cn = new SQLiteConnection(connectionString)
        cn.Open()
        cn

    ///
    let getDBConnectionBy sdbParams =
        let dbFilePath = Db.getNameOf sdbParams
        let connectionString = sprintf "Data Source=%s;Version=3" dbFilePath
        let cn = new SQLiteConnection(connectionString)
        cn.Open()
        cn
        

    ///
    let connectOrCreateDB sdbParams =
        // Check existens by param
        if Db.isExistsBy sdbParams then            
            getDBConnectionBy sdbParams        
        else
            // Create db file
            let dbFileName = Db.getNameOf sdbParams
            Db.initDB dbFileName
            // prepares LookUpMaps of modLookUp based on the dbParams
            let modLookUp = ModCombinator.modLookUpOf sdbParams
            // Set name of global modification
            let connectionString = sprintf "Data Source=%s;Version=3" dbFileName
            let cn = new SQLiteConnection(connectionString)
            cn.Open()
            let _ = Db.insertSdbParams cn sdbParams
            cn.Close()
            // Read Fasta
            let fasta = 
                BioFSharp.IO.FastA.fromFile (BioArray.ofAminoAcidString) sdbParams.FastaPath                
            // Digest
            let peptideId = ref 1 
            fasta
            |> Seq.mapi 
                (fun i fastaItem ->
                    let proteinId = i // TODO                        
                    let peptideContainer =
                        Digestion.BioArray.digest sdbParams.Protease proteinId fastaItem.Sequence
                        |> Digestion.BioArray.concernMissCleavages sdbParams.MinMissedCleavages sdbParams.MaxMissedCleavages
                        |> Array.filter (fun x -> 
                                            let cleavageRange = x.CleavageEnd - x.CleavageStart
                                            cleavageRange > sdbParams.MinPepLength && cleavageRange < sdbParams.MaxPepLength 
                                        )
                        |> (ModCombinator.getModSequencesOf peptideId modLookUp sdbParams.VarModThreshold sdbParams.MassFunction)                     
                    createProteinContainer 
                        proteinId 
                            (sdbParams.FastaHeaderToName fastaItem.Header) 
                                (BioArray.toString fastaItem.Sequence) 
                                    (peptideContainer |> List.ofArray)
                ) 
            |> Db.bulkInsert cn
            |> ignore 
            cn.Open()
            cn

    /// Returns a list of proteins retrieved by PepsequenceID
    let getProteinLookUpFromFileBy sdbParams = 
        let dbFilePath = Db.getNameOf sdbParams
        let connectionString = sprintf "Data Source=%s;Version=3" dbFilePath
        let cn = getDBConnectionBy sdbParams
        cn.Open()
        let tr = cn.BeginTransaction()
        let selectCleavageIdxByPepSeqID   = Db.SQLiteQuery.prepareSelectCleavageIndexByPepSequenceID cn tr
        let selectProteinByProtID         = Db.SQLiteQuery.prepareSelectProteinByID cn tr        
        (fun pepSequenceID -> 
                selectCleavageIdxByPepSeqID pepSequenceID
                |> List.map (fun (_,protID,_,_,_,_) -> selectProteinByProtID protID )
        )

     /// Returns a LookUpResult list
    let getPeptideLookUpFromFileBy sdbParams = 
        let dbFilePath = Db.getNameOf sdbParams
        let connectionString = sprintf "Data Source=%s;Version=3" dbFilePath
        let cn = new SQLiteConnection(connectionString)
        cn.Open()
        let tr = cn.BeginTransaction()
        let parseAAString = initOfModAminoAcidString sdbParams.IsotopicMod (sdbParams.FixedMods@sdbParams.VariableMods)
        let selectModsequenceByMassRange = Db.SQLiteQuery.prepareSelectModsequenceByMassRangeWith cn tr
        (fun lowerMass upperMass  -> 
                let lowerMass' = Convert.ToInt64(lowerMass*1000000.)
                let upperMass' = Convert.ToInt64(upperMass*1000000.)
                selectModsequenceByMassRange lowerMass' upperMass'
                |> List.map (createLookUpResultBy parseAAString)
        )

    /// Prepares a function which returns a list of protein Accessions tupled with the peptide sequence whose ID they were retrieved by
    let getProteinPeptideLookUpFromFileBy (memoryDB: SQLiteConnection) =
        let tr = memoryDB.BeginTransaction()
        let selectCleavageIdxByPepSeqID   = Db.SQLiteQuery.prepareSelectCleavageIndexByPepSequenceID memoryDB tr
        let selectProteinByProtID         = Db.SQLiteQuery.prepareSelectProteinAccessionByID memoryDB tr
        let selectPeptideByPepSeqID       = Db.SQLiteQuery.prepareSelectPepSequenceByPepSequenceID memoryDB tr
        (
            fun pepSequenceID ->
                selectCleavageIdxByPepSeqID pepSequenceID
                |> List.map (fun (_,protID,pepID,_,_,_) -> selectProteinByProtID protID, selectPeptideByPepSeqID pepID )
        )

        /// Returns a LookUpResult list
    let getThreadSafePeptideLookUpFromFileBy (cn:SQLiteConnection) sdbParams = 
        let parseAAString = initOfModAminoAcidString sdbParams.IsotopicMod (sdbParams.FixedMods@sdbParams.VariableMods)
        let selectModsequenceByMassRange = Db.SQLiteQuery.prepareSelectModsequenceByMassRange cn 
        (fun lowerMass upperMass  -> 
                let lowerMass' = Convert.ToInt64(lowerMass*1000000.)
                let upperMass' = Convert.ToInt64(upperMass*1000000.)
                selectModsequenceByMassRange lowerMass' upperMass'
                |> List.map (createLookUpResultBy parseAAString)
        )   
   
    let copyDBIntoMemory (cn:SQLiteConnection) = 
        //cn.Open()
        let inMemoryDB = new SQLiteConnection("Data Source=:memory:;cache=shared;Version=3")
        inMemoryDB.Open()
        cn.BackupDatabase(inMemoryDB,"main","main",2000,null,10000)
        //cn.Close()
        inMemoryDB

    /// Returns a LookUpResult list 
    let getPeptideLookUpBy (sdbParams:SearchDbParams) =
        // Check existens by param
        if Db.isExistsBy sdbParams then            
            getPeptideLookUpFromFileBy sdbParams        
        else
            connectOrCreateDB sdbParams |> ignore
            getPeptideLookUpFromFileBy sdbParams
                
    ///    
    let getPeptideLookUpWithMemBy calcIonSeries massFunction (lookUpF: float -> float -> LookUpResult<AminoAcids.AminoAcid> list) (lookUpCache: Cache.Cache<int64,((LookUpResult<AminoAcids.AminoAcid>*Fragmentation.FragmentMasses) list)>) lowerMass upperMass = 
        let lowerMassAsInt = Convert.ToInt64(lowerMass * 1000000.)
        let upperMassAsInt = Convert.ToInt64(upperMass * 1000000.)
        match Cache.containsItemsBetween lookUpCache (lowerMassAsInt,upperMassAsInt) with
        // The cache contains items that lie within the search range
        | Some (lowerMassIdx,upperMassIdx)  -> 
            let lowermassAsInt' = lookUpCache.Keys.[upperMassIdx]
            if lowermassAsInt' < upperMassAsInt then
                let reducedLookUp = lookUpF (Convert.ToDouble(lowermassAsInt') / 1000000.) upperMass
                match reducedLookUp with
                | h::tail -> 
                    let updateCache =
                        reducedLookUp
                        |> List.map (fun lookUpResult -> 
                                        let ionSeries = calcIonSeries massFunction lookUpResult.BioSequence
                                        lookUpResult,ionSeries
                                    )                                       
                        |> List.groupBy (fun (lookUpResult,calcIonSeries) -> Convert.ToInt64(lookUpResult.Mass * 1000000.))
                        |> List.map (Cache.addItem lookUpCache)
                    let lookUpResults = 
                        match Cache.containsItemsBetween lookUpCache (lowerMassAsInt,upperMassAsInt) with 
                        | Some (lowerMassIdx',upperMassIdx') ->
                            Cache.getValuesByIdx lookUpCache (lowerMassIdx',upperMassIdx')
                            |> List.concat
                        | None -> 
                            []   
                    lookUpResults
                | []     ->
                    let lookUpResults = 
                        Cache.getValuesByIdx lookUpCache (lowerMassIdx, upperMassIdx)
                        |> List.concat
                    lookUpResults
            else
                let lookUpResults = 
                    Cache.getValuesByIdx lookUpCache (lowerMassIdx, upperMassIdx)
                    |> List.concat
                lookUpResults
        // The cache contains no items of this mass that lie within the search range.
        | None ->         
            let lookUpResults = 
                lookUpF lowerMass upperMass
                |> List.map (fun lookUpResult -> 
                                let ionSeries = calcIonSeries massFunction lookUpResult.BioSequence
                                lookUpResult,ionSeries
                            )   
            /// update the cache
            lookUpResults
            |> List.groupBy (fun (lookUpResult,calcIonSeries) -> Convert.ToInt64(lookUpResult.Mass * 1000000.))
            |> List.map (Cache.addItem lookUpCache)
            |> ignore 
            lookUpResults
                    
    /// Returns SearchDbParams of a existing database by filePath
    let getSDBParamsBy filePath = 
        let connectionString = sprintf "Data Source=%s;Version=3" filePath
        let cn = new SQLiteConnection(connectionString)
        cn.Open()
        match Db.SQLiteQuery.selectSearchDbParams cn with 
        | Some (iD,name,fo,fp,pr,minmscl,maxmscl,mass,minpL,maxpL,isoL,mMode,fMods,vMods,vThr) -> 
            createSearchDbParams 
                name fo fp id (Digestion.Table.getProteaseBy pr) minmscl maxmscl mass minpL maxpL 
                    (Newtonsoft.Json.JsonConvert.DeserializeObject<SearchInfoIsotopic list>(isoL)) (Newtonsoft.Json.JsonConvert.DeserializeObject<MassMode>(mMode)) (massFBy (Newtonsoft.Json.JsonConvert.DeserializeObject<MassMode>(mMode))) 
                        (Newtonsoft.Json.JsonConvert.DeserializeObject<SearchModification list>(fMods)) (Newtonsoft.Json.JsonConvert.DeserializeObject<SearchModification list>(vMods)) vThr
        | None ->
            failwith "This database does not contain any SearchParameters. It is not recommended to work with this file."

    /// Returns SearchDbParams of a existing database by SQLiteConnection
    let getSDBParamsByCn (cn :SQLiteConnection)=
        let cn =
            match cn.State with
            | ConnectionState.Open ->
                cn
            | ConnectionState.Closed ->
                cn.Open()
                cn
            | _ as x -> failwith "Data base is busy."
        match Db.SQLiteQuery.selectSearchDbParams cn with
        | Some (iD,name,fo,fp,pr,minmscl,maxmscl,mass,minpL,maxpL,isoL,mMode,fMods,vMods,vThr) ->
            createSearchDbParams
                name fo fp id (Digestion.Table.getProteaseBy pr) minmscl maxmscl mass minpL maxpL
                    (Newtonsoft.Json.JsonConvert.DeserializeObject<SearchInfoIsotopic list>(isoL)) (Newtonsoft.Json.JsonConvert.DeserializeObject<MassMode>(mMode)) (massFBy (Newtonsoft.Json.JsonConvert.DeserializeObject<MassMode>(mMode)))
                        (Newtonsoft.Json.JsonConvert.DeserializeObject<SearchModification list>(fMods)) (Newtonsoft.Json.JsonConvert.DeserializeObject<SearchModification list>(vMods)) vThr
        | None ->
            failwith "This database does not contain any SearchParameters. It is not recommended to work with this file."

    module Table =
        
        let phosphorylation'Ser'Thr'Tyr' =
            createSearchModification "Phosphorylation'Ser'Thr'Tyr'" "21" "Addition of a phosphate group to the AA residue" true "HO3P"
                [Specific(Ser,ModLocation.Residual); Specific(Thr,ModLocation.Residual); Specific(Tyr,ModLocation.Residual)] SearchModType.Plus "ph"
    
        let acetylation'ProtNTerm' =
            createSearchModification "Acetyl(Protein N-Term)" "1" "Acetylation of the protein N-terminus" true "C2H2O"
                [Any(ModLocation.ProteinNterm)] SearchModType.Plus "ac"

        let oxidation'Met' =
            createSearchModification "Oxidation'Met'" "35" "Oxidation" true "O"
                [Specific(Met,ModLocation.Residual);] SearchModType.Plus "ox"

        let carbamidomethyl'Cys' =
            createSearchModification "Carbamidomethyl'Cys'" "4" "Iodoacetamide derivative" false "C2H3NO"
                [Specific(Cys,ModLocation.Residual)] SearchModType.Plus "ca"
                                        