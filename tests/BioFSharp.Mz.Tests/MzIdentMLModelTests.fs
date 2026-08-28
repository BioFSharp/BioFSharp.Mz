module MzIdentMLModelTests

open System
open System.Data
open System.Data.SQLite
open System.IO
open Expecto
open BioFSharp.Mz
open BioFSharp.Mz.MzIdentMLModel

let private parameterData () =
    let t1 = Term.create "MS:1" "MS" "scan time"
    let t2 = Term.create "MS:2" "MS" "charge"

    [
        CvParam.create (Guid.NewGuid()) t1 (box 12.5 :?> IConvertible)
        CvParam.create (Guid.NewGuid()) t2 (box 2 :?> IConvertible)
        CvParam.create (Guid.NewGuid()) t1 (box 99.0 :?> IConvertible)
    ]

let private withTemporaryDirectory action =
    let directory =
        Path.Combine(
            Path.GetTempPath(),
            "BioFSharpMzTests_" + Guid.NewGuid().ToString("N")
        )

    Directory.CreateDirectory(directory) |> ignore

    try
        action directory
    finally
        try SQLiteConnection.ClearAllPools() with _ -> ()
        try GC.Collect() with _ -> ()
        try GC.WaitForPendingFinalizers() with _ -> ()
        try GC.Collect() with _ -> ()
        try Directory.Delete(directory, true) with _ -> ()

[<Tests>]
let tests =
    testList "MzIdentMLModelTests" [
        testList "ParamContainer" [
            testCase "ofSeq stores parameters keyed by term id with last-one-wins on duplicates" <| fun _ ->
                let ps = parameterData ()
                let c = ParamContainer.ofSeq ps

                Expect.equal c.Count 2 "the container stores one value for each distinct term id"
                let actual =
                    ParamContainer.getCvParam "MS:1" c
                    |> fun param -> Convert.ToDouble(param.Value)
                Expect.floatClose Accuracy.high actual 99.0 "the later parameter replaces the earlier value"

            testCase "typed accessors convert stored values and return documented defaults when absent" <| fun _ ->
                let ps = parameterData ()
                let c = ParamContainer.ofSeq ps

                Expect.floatClose Accuracy.high (ParamContainer.getValueAsFloat "MS:1" c) 99.0 "the float accessor converts the stored value"
                Expect.equal (ParamContainer.getValueAsInt "MS:2" c) 2 "the int accessor converts the stored value"
                Expect.equal (ParamContainer.getValueAsString "MS:2" c) "2" "the string accessor converts the stored value"
                Expect.isTrue (Double.IsNaN (ParamContainer.getValueAsFloat "MS:404" c)) "the missing float accessor returns NaN"
                Expect.equal (ParamContainer.getValueAsInt "MS:404" c) -1 "the missing int accessor returns -1"
                Expect.equal (ParamContainer.getValueAsString "MS:404" c) "" "the missing string accessor returns an empty string"
                Expect.equal (ParamContainer.tryGetCvParam "MS:404" c) None "the missing parameter is absent"
                Expect.isSome (ParamContainer.tryGetCvParam "MS:1" c) "the stored parameter is present"

            testCase "addOrUpdateInPlace updates an existing term's value in the same container" <| fun _ ->
                let t1 = Term.create "MS:1" "MS" "scan time"
                let container =
                    ParamContainer.ofSeq [
                        CvParam.create (Guid.NewGuid()) t1 (box 12.5 :?> IConvertible)
                    ]

                ParamContainer.addOrUpdateInPlace
                    (CvParam.create (Guid.NewGuid()) t1 (box 55.0 :?> IConvertible))
                    container
                |> ignore

                Expect.floatClose Accuracy.high (ParamContainer.getValueAsFloat "MS:1" container) 55.0 "the existing term is updated in place"
                Expect.equal container.Count 1 "updating an existing term does not add a second entry"

            testCase "a replaced parameter keeps value and unit coupled" <| fun _ ->
                let fixedTs = DateTime(2026, 1, 15)
                let scanStartTime = Term.initOf "MS:1000016" "MS" "scan start time" fixedTs
                let secondTerm = Term.initOf "UO:0000010" "UO" "second" fixedTs
                let minuteTerm = Term.initOf "UO:0000031" "UO" "minute" fixedTs
                let container =
                    ParamContainer.ofSeq [
                        CvParam.createWithUnit
                            (Guid.NewGuid())
                            scanStartTime
                            (box 90.0 :?> IConvertible)
                            secondTerm
                    ]

                ParamContainer.addOrUpdateInPlace
                    (CvParam.createWithUnit
                        (Guid.NewGuid())
                        scanStartTime
                        (box 1.5 :?> IConvertible)
                        minuteTerm)
                    container
                |> ignore

                let stored = ParamContainer.getCvParam "MS:1000016" container
                Expect.equal container.Count 1 "replacing a parameter does not add a second term entry"
                Expect.floatClose Accuracy.high (Convert.ToDouble stored.Value) 1.5 "the replacement keeps its new value"
                Expect.equal (stored.Unit |> Option.map (fun unitTerm -> unitTerm.Id)) (Some minuteTerm.Id) "the replacement keeps its new unit"
                // 90 seconds = 1.5 minutes: a quantity is value AND unit together; replacement leaving a stale unit corrupts every consumer reading times or masses.
        ]

        testList "DataModel" [
            // PENDING: two Terms with equal Id but different RowVersion are unequal under Equals yet
            // compare as 0 under IComparable - so Set.ofList collapses them to one element while
            // List.distinct keeps two. The .NET framework contract requires Equals and CompareTo to agree
            // on identity; this assertion passes under ANY consistent fix (drop RowVersion from Equals, or
            // include it in CompareTo).
            ptestCase "equality and comparison agree on what makes two terms identical" <| fun _ ->
                let t1 = Term.initOf "MS:1" "MS" "name" (DateTime(2020, 1, 1))
                let t2 = Term.initOf "MS:1" "MS" "name" (DateTime(2021, 1, 1))

                Expect.equal
                    (Set.count (Set.ofList [t1; t2]))
                    (List.length (List.distinct [t1; t2]))
                    "Set and List.distinct agree on term identity"
        ]

        testList "SQLite" [
            testCase "entity table DDL creates queryable tables and a prepared insert round-trips a row" <| fun _ ->
                withTemporaryDirectory (fun directory ->
                    let file = Path.Combine(directory, "mzidentml.db")
                    let fixedTs = DateTime(2026, 1, 15)
                    use cn = new SQLiteConnection(sprintf "Data Source=%s;Version=3" file)
                    cn.Open()
                    use tr = cn.BeginTransaction()

                    MzIdentMLModel.Db.DBSequence.createDBSequenceTable cn |> ignore
                    MzIdentMLModel.Db.DBSequence.createDBSequenceParamTable cn |> ignore
                    let insert = MzIdentMLModel.Db.DBSequence.prepareInsertDBSequence cn tr
                    insert 17 "P01234" "ALBU_HUMAN" "Swiss-Prot-2026_01" fixedTs |> ignore
                    tr.Commit()

                    use command = new SQLiteCommand("SELECT ID, Accession, Name, SearchDBID, RowVersion FROM DBSequence WHERE ID = 17", cn)
                    use reader = command.ExecuteReader()
                    Expect.isTrue (reader.Read()) "the inserted row is returned by independent SQL"
                    if reader.Read() then
                        failtest "the query returned more than one row"
                    else
                        ()
                    // The first Read above established the row; re-querying the scalar fields keeps all five column assertions independent of the row-count check.
                    use rowCommand = new SQLiteCommand("SELECT ID, Accession, Name, SearchDBID, RowVersion FROM DBSequence WHERE ID = 17", cn)
                    use rowReader = rowCommand.ExecuteReader()
                    if rowReader.Read() then
                        Expect.equal (rowReader.GetInt32(0)) 17 "the ID round-trips"
                        Expect.equal (rowReader.GetString(1)) "P01234" "the accession round-trips"
                        Expect.equal (rowReader.GetString(2)) "ALBU_HUMAN" "the name round-trips"
                        Expect.equal (rowReader.GetString(3)) "Swiss-Prot-2026_01" "the SearchDBID round-trips"
                        let storedRowVersion = rowReader.GetDateTime(4)
                        if storedRowVersion <> fixedTs then
                            printfn "D7 RowVersion precision observed: expected %O, got %O; comparing to the second" fixedTs storedRowVersion
                        Expect.equal
                            (storedRowVersion.ToString("yyyy-MM-dd HH:mm:ss"))
                            (fixedTs.ToString("yyyy-MM-dd HH:mm:ss"))
                            "the RowVersion round-trips to the stored second"
                    else
                        failtest "the five-column row disappeared during the independent read"
                    Expect.isFalse (reader.Read()) "the query returns exactly one row"

                    use tablesCommand = new SQLiteCommand("SELECT name FROM sqlite_master WHERE type='table'", cn)
                    use tablesReader = tablesCommand.ExecuteReader()
                    let tableNames = [
                        while tablesReader.Read() do
                            yield tablesReader.GetString(0)
                    ]
                    Expect.isTrue (List.contains "DBSequenceParam" tableNames) "the DBSequenceParam table exists"
                )

            // PENDING: the SpectrumIdentificationItem prepared insert binds option values directly to
            // SQLite Int32 parameters and throws before the graph row can be persisted.
            ptestCase "a spectrum identification graph persists and joins back correctly" <| fun _ ->
                withTemporaryDirectory (fun directory ->
                    let file = Path.Combine(directory, "mzidentml-graph.db")
                    let fixedTs = DateTime(2026, 1, 15)
                    use cn = new SQLiteConnection(sprintf "Data Source=%s;Version=3" file)
                    cn.Open()

                    MzIdentMLModel.Db.DBSequence.createDBSequenceTable cn |> ignore
                    MzIdentMLModel.Db.Peptide.createPeptideTable cn |> ignore
                    MzIdentMLModel.Db.Modification.createModificationTable cn |> ignore
                    MzIdentMLModel.Db.ModLocation.createModLocationTable cn |> ignore
                    MzIdentMLModel.Db.PeptideEvidence.createPeptideEvidenceTable cn |> ignore
                    MzIdentMLModel.Db.SpectrumIdentificationResult.createSpectrumIdentificationResultTable cn |> ignore
                    MzIdentMLModel.Db.SpectrumIdentificationItem.createSpectrumIdentificationItemTable cn |> ignore

                    use tr = cn.BeginTransaction()
                    let insertDbSequence = MzIdentMLModel.Db.DBSequence.prepareInsertDBSequence cn tr
                    let insertPeptide = MzIdentMLModel.Db.Peptide.prepareInsertPeptide cn tr
                    let insertModification = MzIdentMLModel.Db.Modification.prepareInsertModification cn tr
                    let insertModLocation = MzIdentMLModel.Db.ModLocation.prepareInsertModLocation cn tr
                    let insertPeptideEvidence = MzIdentMLModel.Db.PeptideEvidence.prepareInsertPeptideEvidence cn tr
                    let insertSpectrumIdentificationResult = MzIdentMLModel.Db.SpectrumIdentificationResult.prepareInsertSpectrumIdentificationResult cn tr
                    let insertSpectrumIdentificationItem = MzIdentMLModel.Db.SpectrumIdentificationItem.prepareInsertSpectrumIdentificationItem cn tr

                    insertDbSequence 1 "P100" "protein" "sdb" fixedTs |> ignore
                    insertPeptide "pep1" "ACDMK" fixedTs |> ignore
                    insertModification 1 "Oxidation" "M" 15.994915 15.99 fixedTs |> ignore
                    insertModLocation 1 1 1 4 "M" fixedTs |> ignore
                    insertPeptideEvidence 1 1 1 None None None None None None None fixedTs |> ignore
                    insertSpectrumIdentificationResult 1 "scan=42" "sd1" None None fixedTs |> ignore
                    insertSpectrumIdentificationItem 1 (Some 1) None 1 None None "true" (Some 1) None 500.25 2 None None fixedTs |> ignore
                    tr.Commit()

                    use command = new SQLiteCommand(
                        "SELECT p.Sequence, m.MonoisotopicMassDelta, ml.Location, ml.Residue, dbs.Accession, sir.SpectrumID, sii.Rank, sii.ChargeState, sii.SampleID " +
                        "FROM SpectrumIdentificationItem AS sii " +
                        "JOIN SpectrumIdentificationResult AS sir ON sir.ID = sii.SpectrumIdentificationResultID " +
                        "JOIN Peptide AS p ON p.ID = printf('pep%d', sii.PeptideID) " +
                        "JOIN PeptideEvidence AS pe ON pe.PeptideID = sii.PeptideID " +
                        "JOIN DBSequence AS dbs ON dbs.ID = pe.DBSequenceID " +
                        "JOIN ModLocation AS ml ON ml.PeptideID = pe.PeptideID " +
                        "JOIN Modification AS m ON m.ID = ml.ModificationID " +
                        "WHERE sii.ID = 1",
                        cn)
                    use reader = command.ExecuteReader()
                    Expect.isTrue (reader.Read()) "the hand-written graph join returns the stored identification"
                    if reader.Read() then
                        failtest "the hand-written graph join returned more than one identification"
                    else
                        ()

                    use rowCommand = new SQLiteCommand(
                        "SELECT p.Sequence, m.MonoisotopicMassDelta, ml.Location, ml.Residue, dbs.Accession, sir.SpectrumID, sii.Rank, sii.ChargeState, sii.SampleID " +
                        "FROM SpectrumIdentificationItem AS sii " +
                        "JOIN SpectrumIdentificationResult AS sir ON sir.ID = sii.SpectrumIdentificationResultID " +
                        "JOIN Peptide AS p ON p.ID = printf('pep%d', sii.PeptideID) " +
                        "JOIN PeptideEvidence AS pe ON pe.PeptideID = sii.PeptideID " +
                        "JOIN DBSequence AS dbs ON dbs.ID = pe.DBSequenceID " +
                        "JOIN ModLocation AS ml ON ml.PeptideID = pe.PeptideID " +
                        "JOIN Modification AS m ON m.ID = ml.ModificationID " +
                        "WHERE sii.ID = 1",
                        cn)
                    use rowReader = rowCommand.ExecuteReader()
                    if rowReader.Read() then
                        Expect.equal (rowReader.GetString(0)) "ACDMK" "the peptide sequence joins back"
                        Expect.floatClose Accuracy.high (rowReader.GetDouble(1)) 15.994915 "the modification delta joins back"
                        Expect.equal (rowReader.GetInt32(2)) 4 "the modification location joins back"
                        Expect.equal (rowReader.GetString(3)) "M" "the modified residue joins back"
                        Expect.equal (rowReader.GetString(4)) "P100" "the DBSequence accession joins through peptide evidence"
                        Expect.equal (rowReader.GetString(5)) "scan=42" "the spectrum ID joins back"
                        Expect.equal (rowReader.GetInt32(6)) 1 "the rank joins back"
                        Expect.equal (rowReader.GetInt32(7)) 2 "the charge state joins back"
                        Expect.isTrue (rowReader.IsDBNull(8)) "the deliberately absent SampleID remains NULL"
                    else
                        failtest "the graph join had no row to inspect"
                    // the model's purpose is persisting identifications; a broken column mapping or optional binding silently corrupts every consumer - the join is hand-written so no broken select helper is exercised.
                )

            // PENDING: with PRAGMA foreign_keys = ON, inserting a ProteinDetectionHypothesis whose
            // DBSequenceID does not exist must fail - the table declares no FOREIGN KEY constraint, so the
            // orphan row inserts silently and downstream protein inference dereferences nothing.
            ptestCase "protein hypotheses cannot reference non-existent sequence evidence" <| fun _ ->
                withTemporaryDirectory (fun directory ->
                    let file = Path.Combine(directory, "mzidentml-orphan.db")
                    let fixedTs = DateTime(2026, 1, 15)
                    use cn = new SQLiteConnection(sprintf "Data Source=%s;Version=3" file)
                    cn.Open()
                    MzIdentMLModel.Db.DBSequence.createDBSequenceTable cn |> ignore
                    MzIdentMLModel.Db.ProteinDetectionList.createProteinDetectionListTable cn |> ignore
                    MzIdentMLModel.Db.ProteinAmbiguityGroup.createProteinAmbiguityGroupTable cn |> ignore
                    MzIdentMLModel.Db.ProteinDetectionHypothesis.createProteinDetectionHypothesisTable cn |> ignore
                    use pragma = new SQLiteCommand("PRAGMA foreign_keys=ON", cn)
                    pragma.ExecuteNonQuery() |> ignore
                    use tr = cn.BeginTransaction()
                    let insertList = MzIdentMLModel.Db.ProteinDetectionList.prepareInsertProteinDetectionList cn tr
                    let insertGroup = MzIdentMLModel.Db.ProteinAmbiguityGroup.prepareInsertProteinAmbiguityGroup cn tr
                    insertList 1 "pdl1" "protein detection list" "sdb" fixedTs |> ignore
                    insertGroup 1 1 None fixedTs |> ignore
                    let insertHypothesis = MzIdentMLModel.Db.ProteinDetectionHypothesis.prepareInsertProteinDetectionHypothesis cn tr
                    Expect.throws
                        (fun () -> insertHypothesis 1 999 1 None "true" fixedTs |> ignore)
                        "an orphan DBSequenceID is rejected when foreign keys are enabled"
                )

            // PENDING: the param-table insert SQL's column list is missing its closing ")" (and carries a
            // trailing comma before VALUES), so every call fails with "near VALUES: syntax error" - the
            // entire CvParam-persistence half of the model (insertDBSequenceToDb) can never complete.
            ptestCase "a prepared param insert persists a controlled-vocabulary parameter" <| fun _ ->
                withTemporaryDirectory (fun directory ->
                    let file = Path.Combine(directory, "mzidentml.db")
                    use cn = new SQLiteConnection(sprintf "Data Source=%s;Version=3" file)
                    cn.Open()
                    MzIdentMLModel.Db.DBSequence.createDBSequenceTable cn |> ignore
                    MzIdentMLModel.Db.DBSequence.createDBSequenceParamTable cn |> ignore
                    use tr = cn.BeginTransaction()
                    let insert = MzIdentMLModel.Db.DBSequence.prepareInsertDBSequenceParam cn tr
                    insert 1 7 "MS:1000001" "UO:0000000" "12.5" (DateTime(2020, 1, 1)) |> ignore
                    tr.Commit()

                    use command = new SQLiteCommand("SELECT Value FROM DBSequenceParam WHERE ID = 1 AND FKParamContainer = 7", cn)
                    use reader = command.ExecuteReader()
                    Expect.isTrue (reader.Read()) "the inserted parameter row is returned by independent SQL"
                    Expect.equal (reader.GetString(0)) "12.5" "the parameter value round-trips"
                )

            // PENDING: the select helpers' SQL is malformed - the FROM keyword is concatenated directly
            // to the table name ("FROMDBSequence") and parameter names are declared with stray spaces
            // ("@ id " vs "@id"), so selection fails at the SQLite layer although the documented purpose
            // is retrieving the stored entity.
            ptestCase "the prepared select helper retrieves an inserted row by ID" <| fun _ ->
                withTemporaryDirectory (fun directory ->
                    let file = Path.Combine(directory, "mzidentml.db")
                    use cn = new SQLiteConnection(sprintf "Data Source=%s;Version=3" file)
                    cn.Open()
                    use tr = cn.BeginTransaction()

                    MzIdentMLModel.Db.DBSequence.createDBSequenceTable cn |> ignore
                    let insert = MzIdentMLModel.Db.DBSequence.prepareInsertDBSequence cn tr
                    insert 1 "ACC1" "protein one" "sdb1" DateTime.Now |> ignore
                    tr.Commit()

                    use tr2 = cn.BeginTransaction()
                    let selected = MzIdentMLModel.Db.DBSequence.prepareSelectDBSequencebyID cn tr2 1
                    Expect.isFalse (List.isEmpty selected) "the prepared select returns the inserted row"
                    let _, accession, _, _, _ = List.head selected
                    Expect.equal accession "ACC1" "the selected row has the inserted accession"
                )

            // PENDING: initDB declares its connection with `use`, so the returned connection is already
            // disposed when the caller receives it - any command on it throws ObjectDisposedException.
            // A factory named initDB must hand back a usable connection.
            ptestCase "initDB returns a usable open connection" <| fun _ ->
                withTemporaryDirectory (fun directory ->
                    let tempFile = Path.Combine(directory, "mzidentml-init.db")
                    use cn = MzIdentMLModel.Db.initDB tempFile
                    Expect.equal cn.State ConnectionState.Open "initDB returns an open connection"
                    use command = new SQLiteCommand("SELECT COUNT(*) FROM DBSequence", cn)
                    try
                        command.ExecuteScalar() |> ignore
                    with ex ->
                        failtestf "an open initDB connection executes commands, but threw: %s" ex.Message
                )
        ]
    ]
