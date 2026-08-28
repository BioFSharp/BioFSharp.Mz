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
                    use cn = new SQLiteConnection(sprintf "Data Source=%s;Version=3" file)
                    cn.Open()
                    use tr = cn.BeginTransaction()

                    MzIdentMLModel.Db.DBSequence.createDBSequenceTable cn |> ignore
                    MzIdentMLModel.Db.DBSequence.createDBSequenceParamTable cn |> ignore
                    let insert = MzIdentMLModel.Db.DBSequence.prepareInsertDBSequence cn tr
                    insert 1 "ACC1" "protein one" "sdb1" DateTime.Now |> ignore
                    tr.Commit()

                    use command = new SQLiteCommand("SELECT Accession, Name FROM DBSequence WHERE ID = 1", cn)
                    use reader = command.ExecuteReader()
                    Expect.isTrue (reader.Read()) "the inserted row is returned by independent SQL"
                    Expect.equal (reader.GetString(0)) "ACC1" "the accession round-trips"
                    Expect.equal (reader.GetString(1)) "protein one" "the name round-trips"
                    Expect.isFalse (reader.Read()) "the query returns exactly one row"

                    use tablesCommand = new SQLiteCommand("SELECT name FROM sqlite_master WHERE type='table'", cn)
                    use tablesReader = tablesCommand.ExecuteReader()
                    let tableNames = [
                        while tablesReader.Read() do
                            yield tablesReader.GetString(0)
                    ]
                    Expect.isTrue (List.contains "DBSequenceParam" tableNames) "the DBSequenceParam table exists"
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
