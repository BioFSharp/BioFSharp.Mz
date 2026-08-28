module SearchDBTests

open System
open System.IO
open System.Text.RegularExpressions
open System.Data.SQLite
open Expecto
open BioFSharp
open BioFSharp.Mz

let monoisotopicMass = SearchDB.massFBy SearchDB.MassMode.Monoisotopic

let expectWithin tolerance actual expected message =
    Expect.isTrue
        (abs (actual - expected) <= tolerance)
        (sprintf "%s; expected %g, got %g" message expected actual)

let withTemporaryDirectory action =
    let directory =
        Path.Combine(
            Path.GetTempPath(),
            "BioFSharpMzTests_" + Guid.NewGuid().ToString("N")
        )

    Directory.CreateDirectory(directory) |> ignore

    try
        action directory
    finally
        SQLiteConnection.ClearAllPools()
        GC.Collect()
        GC.WaitForPendingFinalizers()
        GC.Collect()

        try
            Directory.Delete(directory, true)
        with
        | _ -> ()

let stripBracketedModifications sequence =
    Regex.Replace(sequence, @"\[[^\]]*\]", "")

[<Tests>]
let tests =
    testList "SearchDBTests" [
        testList "PureHelpers" [
            testCase "pepLengthLimitsBy derives the length window from charge and upper acquisition bound" <| fun _ ->
                let actual = SearchDB.pepLengthLimitsBy 4 (300.0, 1665.0)
                Expect.equal actual (4, 30) "the peptide length window is derived from charge and the upper acquisition bound"
                // (maxExpCharge/2 * upperMz) / 111 = (2 * 1665)/111 = 30: the /2 is the IMPLEMENTATION'S
                // conservatism heuristic, not the physical maximum (a charge-4 ion at m/z 1665 could weigh
                // ~6656 Da, ~60 residues at the ~111 Da average residue mass). The lower bound is a hardcoded
                // practical minimum of 4. The lower acquisition-window bound is entirely unused.
                Expect.isTrue
                    (snd (SearchDB.pepLengthLimitsBy 6 (300.0, 1665.0)) > snd (SearchDB.pepLengthLimitsBy 4 (300.0, 1665.0)))
                    "higher allowed charge -> heavier observable peptides -> longer window (direction is domain truth even if the coefficient is a heuristic)"
                Expect.isTrue
                    (snd (SearchDB.pepLengthLimitsBy 4 (300.0, 2220.0)) > snd (SearchDB.pepLengthLimitsBy 4 (300.0, 1665.0)))
                    "wider acquisition -> longer window"

            testCase "createSearchDbParams stores modification lists sorted" <| fun _ ->
                let fx = [
                    SearchDB.Table.oxidation'Met'
                    SearchDB.Table.acetylation'ProtNTerm'
                ]

                let p =
                    SearchDB.createSearchDbParams
                        "n"
                        "f"
                        "fa"
                        id
                        (Digestion.Table.getProteaseBy "Trypsin")
                        0
                        1
                        3000.0
                        3
                        20
                        []
                        SearchDB.MassMode.Monoisotopic
                        (SearchDB.massFBy SearchDB.MassMode.Monoisotopic)
                        fx
                        []
                        3

                Expect.equal p.FixedMods (List.sort fx) "fixed modifications are stored in sorted order"
                // The parameter record is the database's identity for matching existing DBs; sorting makes that identity independent of caller-supplied list order (deterministic identity contract).

            testCase "getModBy realizes the predefined modifications with their published mass deltas" <| fun _ ->
                let modificationCases = [
                    SearchDB.Table.phosphorylation'Ser'Thr'Tyr', AminoAcids.Ser, 79.966331
                    SearchDB.Table.oxidation'Met', AminoAcids.Met, 15.994915
                    SearchDB.Table.carbamidomethyl'Cys', AminoAcids.Cys, 57.021464
                    SearchDB.Table.acetylation'ProtNTerm', AminoAcids.Ala, 42.010565
                ]

                modificationCases
                |> List.iter (fun (searchMod, residue, delta) ->
                    let md = SearchDB.getModBy searchMod
                    let modified = AminoAcids.setModification md residue
                    let actual =
                        monoisotopicMass (modified :> BioFSharp.IBioItem)
                        - monoisotopicMass (residue :> BioFSharp.IBioItem)

                    expectWithin 0.001 actual delta (sprintf "%s has its published monoisotopic delta" searchMod.Name))
                // These four deltas are published Unimod monoisotopic modification masses - fully independent of this codebase.

            testCase "createIsotopicMod applies a heavy-nitrogen label worth one N15-N14 difference per nitrogen" <| fun _ ->
                let isoMod =
                    SearchDB.createIsotopicMod
                        (SearchDB.createSearchInfoIsotopic "#N15" Elements.Table.N Elements.Table.Heavy.N15)

                let modifiedGlycine = AminoAcids.setModification isoMod AminoAcids.Gly
                let actual =
                    monoisotopicMass (modifiedGlycine :> BioFSharp.IBioItem)
                    - monoisotopicMass (AminoAcids.Gly :> BioFSharp.IBioItem)

                expectWithin 0.001 actual 0.997035 "one nitrogen labeled N15 has the published mass shift"
                // The published N15-N14 mass difference is 0.997035 Da; glycine contains exactly one nitrogen, so a full N15 label shifts its mass by exactly one such difference.

            testCase "initOfModAminoAcidString parses plain sequences and prefix-coded modifications" <| fun _ ->
                let plain = SearchDB.initOfModAminoAcidString [] [] 0 "AGS"
                Expect.equal plain [AminoAcids.Ala; AminoAcids.Gly; AminoAcids.Ser] "a plain sequence is parsed into its amino acids"

                let oxidation = SearchDB.Table.oxidation'Met'
                let code = oxidation.XModCode
                let bracketedCode =
                    if code.StartsWith("[") && code.EndsWith("]") then code
                    else "[" + code + "]"
                let s = bracketedCode + "M"
                let parsed = SearchDB.initOfModAminoAcidString [] [oxidation] 0 s

                Expect.equal (List.length parsed) 1 "the modified sequence contains exactly one residue"
                let modified = List.head parsed
                let actual =
                    monoisotopicMass (modified :> BioFSharp.IBioItem)
                    - monoisotopicMass (AminoAcids.Met :> BioFSharp.IBioItem)
                expectWithin 0.001 actual 15.994915 "the parsed residue carries the published oxidation delta"
                // The mod-string format places the bracketed modification code before its residue; the parsed residue must carry a modification worth the published oxidation delta.
        ]

        testList "SQLiteIntegration" [
            testCase "initDB creates the search database schema" <| fun _ ->
                withTemporaryDirectory (fun directory ->
                    let dbFile = Path.Combine(directory, "schema.db")
                    SearchDB.Db.initDB dbFile

                    use cn = new SQLiteConnection(sprintf "Data Source=%s;Version=3" dbFile)
                    cn.Open()
                    use command = new SQLiteCommand("SELECT name FROM sqlite_master WHERE type='table'", cn)
                    use reader = command.ExecuteReader()
                    let tableNames = [
                        while reader.Read() do
                            yield reader.GetString(0)
                    ]

                    Expect.isTrue (List.contains "SearchDbParams" tableNames) "SearchDbParams table exists"
                    Expect.isTrue (List.contains "Protein" tableNames) "Protein table exists"
                    Expect.isTrue (List.contains "CleavageIndex" tableNames) "CleavageIndex table exists"
                    Expect.isTrue (List.contains "PepSequence" tableNames) "PepSequence table exists"
                    Expect.isTrue (List.contains "ModSequence" tableNames) "ModSequence table exists"
                    // These five tables are the documented storage model of the peptide search database; without them no insert/lookup can work.
                )

            testCase "connectOrCreateDB digests a FASTA into a queryable peptide database" <| fun _ ->
                withTemporaryDirectory (fun directory ->
                    let fastaPath = Path.Combine(directory, "t7.fasta")
                    File.WriteAllText(fastaPath, ">T7\r\nMAGSTKLLVR\r\n")

                    let sdbParams =
                        SearchDB.createSearchDbParams
                            "t7"
                            directory
                            fastaPath
                            id
                            (Digestion.Table.getProteaseBy "Trypsin")
                            0
                            0
                            2000.0
                            2
                            20
                            []
                            SearchDB.MassMode.Monoisotopic
                            monoisotopicMass
                            []
                            [SearchDB.Table.oxidation'Met']
                            1

                    // MinPepLength 2 admits 4-mers: the digestion filter compares (CleavageEnd - CleavageStart),
                    // which is length - 1, strictly against the bound - so "MinPepLength n" actually means
                    // "length >= n + 2" (length-filter off-by-two).
                    use connection = SearchDB.connectOrCreateDB sdbParams
                    Expect.isTrue (SearchDB.Db.isExistsBy sdbParams) "the created database exists for its parameters"

                    let lookUp = SearchDB.getThreadSafePeptideLookUpFromFileBy connection sdbParams
                    let results = lookUp 0.0 2000.0
                    Expect.isTrue (not (List.isEmpty results)) "the peptide database returns query results"

                    let plainSequences =
                        results
                        |> List.map (fun result -> stripBracketedModifications result.StringSequence)

                    Expect.isTrue (List.contains "MAGSTK" plainSequences) "the tryptic MAGSTK peptide is present"
                    Expect.isTrue (List.contains "LLVR" plainSequences) "the tryptic LLVR peptide is present"
                    // initiator-Met excision: peptides starting at position 0 with Met also yield the Met-cleaved variant - documented protein N-terminal processing.
                    Expect.isTrue (List.contains "AGSTK" plainSequences) "the initiator-Met-cleaved tryptic peptide is present"

                    let llvrMass = 2.0 * 113.08406 + 99.06841 + 156.10111 + 18.010565  // Leu+Leu+Val+Arg residues + water, published monoisotopic values: 499.34821 Da
                    let narrowResults = lookUp (llvrMass - 0.005) (llvrMass + 0.005)
                    // the mass-range lookup must find exactly the peptide whose hand-computed neutral mass falls in the window - the one fully independent absolute-mass oracle of the pipeline.
                    Expect.equal narrowResults.Length 1 "the hand-computed LLVR mass window returns exactly one result"
                    Expect.equal
                        (narrowResults |> List.head |> fun result -> stripBracketedModifications result.StringSequence)
                        "LLVR"
                        "the hand-computed mass window returns LLVR"

                    plainSequences
                    |> List.iter (fun sequence ->
                        Expect.isTrue
                            ("MAGSTKLLVR".Contains(sequence))
                            (sprintf "the plain peptide sequence %s comes from the FASTA sequence" sequence))

                    results
                    |> List.iter (fun result ->
                        let residueMass =
                            result.BioSequence
                            |> List.sumBy (fun residue -> monoisotopicMass (residue :> BioFSharp.IBioItem))

                        // serialize->parse round trip: the stored ModString parses back to a residue list whose re-summed mass (same mass function) plus the published water mass equals the stored mass - catches format or mod-loss corruption, not absolute mass accuracy (that is the LLVR window probe's job).
                        expectWithin
                            0.001
                            result.Mass
                            (residueMass + 18.010565)
                            (sprintf "the mass of %s is the residue sum plus water" result.StringSequence))

                    let hasOxidationPair =
                        results
                        |> List.groupBy (fun result -> result.PepSequenceID)
                        |> List.exists (fun (_, peptideResults) ->
                            peptideResults
                            |> List.mapi (fun index result ->
                                peptideResults
                                |> List.skip (index + 1)
                                |> List.exists (fun other ->
                                    abs (abs (result.Mass - other.Mass) - 15.994915) <= 0.001))
                            |> List.exists id)

                    // 15.994915 is the independently known monoisotopic oxidation delta for methionine.
                    Expect.isTrue
                        hasOxidationPair
                        "two modified forms share a PepSequenceID and differ by one methionine oxidation"
                )

            // PENDING: MaxMass is recorded but never applied during candidate generation.
            // Observed sequences: ["WWWWK"; "LLVR"].
            ptestCase "the configured MaxMass excludes over-mass candidates" <| fun _ ->
                withTemporaryDirectory (fun directory ->
                    let fastaPath = Path.Combine(directory, "maxmass.fasta")
                    File.WriteAllText(fastaPath, ">sp|TESTPROT2|T2\r\nLLVRWWWWK\r\n")

                    let sdbParams =
                        SearchDB.createSearchDbParams
                            "testdb_maxmass"
                            directory
                            fastaPath
                            id
                            (Digestion.Table.getProteaseBy "Trypsin")
                            0
                            0
                            600.0
                            2
                            20
                            []
                            SearchDB.MassMode.Monoisotopic
                            monoisotopicMass
                            []
                            []
                            0

                    use connection = SearchDB.connectOrCreateDB sdbParams
                    let lookUp = SearchDB.getThreadSafePeptideLookUpFromFileBy connection sdbParams
                    let sequences =
                        lookUp 100.0 5000.0
                        |> List.map (fun result -> stripBracketedModifications result.StringSequence)

                    Expect.isTrue (List.contains "LLVR" sequences) "the under-MaxMass LLVR candidate is present"
                    Expect.isFalse (List.contains "WWWWK" sequences) "the over-MaxMass WWWWK candidate is excluded"
                    // MaxMass defines the search space; admitting over-mass candidates inflates the database and distorts multiple-testing statistics.
                )

            testCase "fixed modifications are mandatory and site-specific; variable modifications respect sites and the threshold" <| fun _ ->
                withTemporaryDirectory (fun directory ->
                    let fastaPath = Path.Combine(directory, "site-mods.fasta")
                    File.WriteAllText(fastaPath, ">sp|TESTPROT3|T3\r\nACDKMAMR\r\n")

                    let sdbParams =
                        SearchDB.createSearchDbParams
                            "testdb_site_mods"
                            directory
                            fastaPath
                            id
                            (Digestion.Table.getProteaseBy "Trypsin")
                            0
                            0
                            3000.0
                            2
                            20
                            []
                            SearchDB.MassMode.Monoisotopic
                            monoisotopicMass
                            [SearchDB.Table.carbamidomethyl'Cys'; SearchDB.Table.acetylation'ProtNTerm']
                            [SearchDB.Table.oxidation'Met']
                            1

                    let baseMass sequence =
                        SearchDB.initOfModAminoAcidString [] [] 0 sequence
                        |> List.sumBy (fun residue -> monoisotopicMass (residue :> BioFSharp.IBioItem))
                        |> (+) 18.010565

                    use connection = SearchDB.connectOrCreateDB sdbParams
                    let lookUp = SearchDB.getThreadSafePeptideLookUpFromFileBy connection sdbParams
                    let grouped =
                        lookUp 100.0 3000.0
                        |> List.groupBy (fun result -> stripBracketedModifications result.StringSequence)

                    let formsOf sequence =
                        grouped
                        |> List.tryFind (fun (plainSequence, _) -> plainSequence = sequence)
                        |> Option.map snd
                        |> Option.defaultValue []

                    let acdkForms = formsOf "ACDK"
                    let acdkBase = baseMass "ACDK"
                    let acdkExpectedMasses = [acdkBase + 57.021464; acdkBase + 57.021464 + 42.010565]
                    acdkForms
                    |> List.iter (fun result ->
                        Expect.isTrue
                            (acdkExpectedMasses |> List.exists (fun expected -> abs (result.Mass - expected) <= 0.001))
                            (sprintf "ACDK form %s carries mandatory Cys carbamidomethylation" result.StringSequence))
                    // fixed modifications are MANDATORY: with carbamidomethyl and protein-N-terminal
                    // acetylation both fixed and no variable mod applicable to ACDK, exactly one form
                    // exists - carrying both deltas (observed and now pinned).
                    Expect.equal acdkForms.Length 1 "fixed-only modification state yields exactly one ACDK form"
                    Expect.isTrue
                        (acdkForms |> List.exists (fun result -> abs (result.Mass - acdkExpectedMasses.[1]) <= 0.001))
                        "the protein-N-terminal ACDK form additionally carries acetylation"

                    let mamrForms = formsOf "MAMR"
                    let mamrBase = baseMass "MAMR"
                    Expect.equal mamrForms.Length 3 "the two Met sites produce base plus two singly oxidized MAMR forms"
                    if mamrForms.Length = 3 then
                        let actualMasses = mamrForms |> List.map (fun result -> result.Mass) |> List.sort
                        let expectedMasses = [mamrBase; mamrBase + 15.994915; mamrBase + 15.994915] |> List.sort
                        actualMasses
                        |> List.iter2 (fun actual expected -> expectWithin 0.001 actual expected "the MAMR form has the expected oxidation count") expectedMasses
                        Expect.isFalse
                            (mamrForms |> List.exists (fun result -> abs (result.Mass - (mamrBase + 2.0 * 15.994915)) <= 0.001))
                            "the variable-modification threshold excludes doubly oxidized MAMR"
                        Expect.isFalse
                            (mamrForms |> List.exists (fun result -> abs (result.Mass - (mamrBase + 42.010565)) <= 0.001))
                            "protein-terminal acetylation does not leak onto internal MAMR"
                    // fixed = mandatory on its site; variable = optional, one per threshold, Met-only; terminal mods bind to the protein terminus - each clause is a distinct search-space contract, all masses hand-composed from published deltas.
                )

            testCase "the missed-cleavage bounds select the stored candidate set" <| fun _ ->
                withTemporaryDirectory (fun directory ->
                    let fastaPath = Path.Combine(directory, "missed-cleavages.fasta")
                    File.WriteAllText(fastaPath, ">sp|TESTPROT4|T4\r\nAKLLVR\r\n")

                    let sdbParams =
                        SearchDB.createSearchDbParams
                            "testdb_missed_cleavages"
                            directory
                            fastaPath
                            id
                            (Digestion.Table.getProteaseBy "Trypsin")
                            1
                            1
                            3000.0
                            0
                            20
                            []
                            SearchDB.MassMode.Monoisotopic
                            monoisotopicMass
                            []
                            []
                            0

                    use connection = SearchDB.connectOrCreateDB sdbParams
                    let lookUp = SearchDB.getThreadSafePeptideLookUpFromFileBy connection sdbParams
                    let sequences =
                        lookUp 100.0 3000.0
                        |> List.map (fun result -> stripBracketedModifications result.StringSequence)

                    Expect.isTrue (List.contains "AKLLVR" sequences) "the one-missed-cleavage product is present"
                    Expect.isFalse (List.contains "AK" sequences) "the zero-missed-cleavage AK product is excluded"
                    Expect.isFalse (List.contains "LLVR" sequences) "the zero-missed-cleavage LLVR product is excluded"
                    // missed-cleavage bounds are a core digestion contract controlling database size and sensitivity; the products are hand-derived from trypsin specificity (cleave after K/R).
                )

            testCase "database identity is parameter-sensitive and parameters round-trip" <| fun _ ->
                withTemporaryDirectory (fun directory ->
                    let fastaPath = Path.Combine(directory, "identity.fasta")
                    File.WriteAllText(fastaPath, ">sp|TESTPROT5|T5\r\nAKLLVR\r\n")

                    let creatingParams =
                        SearchDB.createSearchDbParams
                            "testdb_identity"
                            directory
                            fastaPath
                            id
                            (Digestion.Table.getProteaseBy "Trypsin")
                            0
                            1
                            2750.0
                            2
                            20
                            []
                            SearchDB.MassMode.Monoisotopic
                            monoisotopicMass
                            [SearchDB.Table.carbamidomethyl'Cys']
                            [SearchDB.Table.oxidation'Met']
                            1

                    use connection = SearchDB.connectOrCreateDB creatingParams
                    let dbFile = SearchDB.Db.getNameOf creatingParams
                    let storedParams = SearchDB.getSDBParamsBy dbFile
                    let projection (p: SearchDB.SearchDbParams) =
                        p.Name,
                        p.Protease.Name,
                        p.MinMissedCleavages,
                        p.MaxMissedCleavages,
                        p.MaxMass,
                        p.MinPepLength,
                        p.MaxPepLength,
                        p.MassMode.ToString(),
                        (p.FixedMods |> List.map (fun modification -> modification.Name)),
                        (p.VariableMods |> List.map (fun modification -> modification.Name))

                    Expect.equal (projection storedParams) (projection creatingParams) "persisted search parameters round-trip through the database"
                    let changedParams = { creatingParams with MaxMissedCleavages = 2 }
                    Expect.isFalse
                        (SearchDB.Db.isExistsBy changedParams)
                        "changing a digestion parameter identifies a different database"
                    // a database is identified by its search parameters, not its file name: a changed digestion setting must not silently reuse an incompatible database, and stored parameters must survive persistence.
                )

            // PENDING: MinPepLength = 4 must admit the 4-residue tryptic peptide LLVR. The digestion
            // filter compares (CleavageEnd - CleavageStart) = length-1 strictly, so the effective minimum
            // is MinPepLength + 2 and LLVR is absent even at MinPepLength 3.
            ptestCase "the peptide length filter admits peptides of exactly the configured minimum length" <| fun _ ->
                withTemporaryDirectory (fun directory ->
                    let fastaPath = Path.Combine(directory, "t7.fasta")
                    File.WriteAllText(fastaPath, ">T7\r\nMAGSTKLLVR\r\n")

                    let sdbParams =
                        SearchDB.createSearchDbParams
                            "testdb2"
                            directory
                            fastaPath
                            id
                            (Digestion.Table.getProteaseBy "Trypsin")
                            0
                            0
                            2000.0
                            4
                            20
                            []
                            SearchDB.MassMode.Monoisotopic
                            monoisotopicMass
                            []
                            [SearchDB.Table.oxidation'Met']
                            1

                    use connection = SearchDB.connectOrCreateDB sdbParams
                    let lookUp = SearchDB.getThreadSafePeptideLookUpFromFileBy connection sdbParams
                    let plainSequences =
                        lookUp 0.0 2000.0
                        |> List.map (fun result -> stripBracketedModifications result.StringSequence)

                    Expect.isTrue (List.contains "LLVR" plainSequences) "the configured minimum-length LLVR peptide is present"
                )

            // PENDING: the function opens an already-open connection (getDBConnectionBy opens; the
            // function calls Open() again) and throws InvalidOperationException on every call.
            ptestCase "getProteinLookUpFromFileBy returns the protein for a peptide sequence ID" <| fun _ ->
                withTemporaryDirectory (fun directory ->
                    let fastaPath = Path.Combine(directory, "t7.fasta")
                    File.WriteAllText(fastaPath, ">sp|TESTPROT1|TEST\r\nMAGSTKLLVR\r\n")

                    let sdbParams =
                        SearchDB.createSearchDbParams
                            "t7"
                            directory
                            fastaPath
                            id
                            (Digestion.Table.getProteaseBy "Trypsin")
                            0
                            0
                            2000.0
                            2
                            20
                            []
                            SearchDB.MassMode.Monoisotopic
                            monoisotopicMass
                            []
                            [SearchDB.Table.oxidation'Met']
                            1

                    use connection = SearchDB.connectOrCreateDB sdbParams
                    let lookUp = SearchDB.getThreadSafePeptideLookUpFromFileBy connection sdbParams
                    let peptideResults = lookUp 0.0 2000.0
                    let pepSequenceID = peptideResults |> List.head |> fun result -> result.PepSequenceID
                    let protLookup = SearchDB.getProteinLookUpFromFileBy sdbParams
                    let proteins = protLookup pepSequenceID

                    Expect.isFalse (List.isEmpty proteins) "the peptide sequence ID maps to a protein"
                    let (_, accession, _) = List.head proteins
                    Expect.equal accession "sp|TESTPROT1|TEST" "the protein accession is returned"
                )

            // PENDING: the match branches are inverted - the empty-list case builds the map (of nothing)
            // and every non-empty list returns Map.empty, so the function ALWAYS returns an empty map.
            ptestCase "xModToSearchModifications maps modification names to their XMod codes" <| fun _ ->
                let m = SearchDB.Db.xModToSearchModifications [SearchDB.Table.oxidation'Met']
                Expect.equal m.Count 1 "the oxidation modification is mapped by its XMod code"

            // PENDING: the SQL references a nonexistent column `Mass` (the table has RealMass and RoundedMass)
            // and fails with "no such column" on every call.
            ptestCase "prepareSelectModsequenceByMass retrieves a stored mod sequence by its mass" <| fun _ ->
                withTemporaryDirectory (fun directory ->
                    let dbFile = Path.Combine(directory, "t6-mass.db")
                    SearchDB.Db.initDB dbFile

                    use connection = new SQLiteConnection(sprintf "Data Source=%s;Version=3" dbFile)
                    connection.Open()
                    use transaction = connection.BeginTransaction()
                    let insertModSequence = SearchDB.Db.SQLiteQuery.prepareInsertModSequence connection transaction
                    let insertedSequence = "LLVR"
                    let insertedMass = 499.348205
                    let insertedRoundedMass = 499348205L
                    insertModSequence 1 insertedMass insertedRoundedMass insertedSequence 0 |> ignore
                    transaction.Commit()

                    let selectModsequenceByMass = SearchDB.Db.SQLiteQuery.prepareSelectModsequenceByMass connection transaction
                    let selected = selectModsequenceByMass (int insertedRoundedMass)
                    let (_, _, _, _, sequence, _) = selected
                    Expect.equal sequence insertedSequence "the stored mod sequence is retrieved by its mass"
                )

            // PENDING: the sequence parameter is bound as DbType.Double, so any string argument raises FormatException before the query runs.
            ptestCase "prepareSelectModsequenceBySequence retrieves a stored mod sequence by its sequence" <| fun _ ->
                withTemporaryDirectory (fun directory ->
                    let dbFile = Path.Combine(directory, "t6-sequence.db")
                    SearchDB.Db.initDB dbFile

                    use connection = new SQLiteConnection(sprintf "Data Source=%s;Version=3" dbFile)
                    connection.Open()
                    use transaction = connection.BeginTransaction()
                    let insertModSequence = SearchDB.Db.SQLiteQuery.prepareInsertModSequence connection transaction
                    let insertedSequence = "LLVR"
                    let insertedMass = 499.348205
                    let insertedRoundedMass = 499348205L
                    insertModSequence 1 insertedMass insertedRoundedMass insertedSequence 0 |> ignore
                    transaction.Commit()

                    let selectModsequenceBySequence = SearchDB.Db.SQLiteQuery.prepareSelectModsequenceBySequence connection transaction
                    let selected = selectModsequenceBySequence insertedSequence
                    let (_, _, mass, _, _, _) = selected
                    Expect.equal mass insertedMass "the stored mod sequence mass is retrieved by its sequence"
                )

        ]
    ]
