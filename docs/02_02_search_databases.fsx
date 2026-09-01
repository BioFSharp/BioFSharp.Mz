(**
---
title: Peptide search databases
category: Peptide identification
categoryindex: 2
index: 2
---
*)

(*** hide ***)

(*** condition: prepare ***)
#r "nuget: FSharpAux, 2.1.0"
#r "nuget: FSharpAux.IO, 2.1.0"
#r "nuget: FSharp.Stats, 0.6.0"
#r "nuget: BioFSharp, 2.0.0"
#r "nuget: System.Data.SQLite.Core, 1.0.119"
#r "nuget: Newtonsoft.Json, 13.0.4"
#r "../src/BioFSharp.Mz/bin/Release/netstandard2.0/BioFSharp.Mz.dll"

(*** condition: ipynb ***)
#if IPYNB
#r "nuget: BioFSharp.Mz, {{fsdocs-package-version}}"
#endif // IPYNB

(*** hide ***)
let dbFolder = System.IO.Path.Combine(System.IO.Path.GetTempPath(), "BioFSharp.Mz_docs_searchdb")
System.IO.Directory.CreateDirectory dbFolder |> ignore

(**
[![Binder]({{root}}img/badge-binder.svg)](https://mybinder.org/v2/gh/BioFSharp/BioFSharp.Mz/gh-pages?filepath={{fsdocs-source-basename}}.ipynb)&emsp;
[![Script]({{root}}img/badge-script.svg)]({{root}}{{fsdocs-source-basename}}.fsx)&emsp;
[![Notebook]({{root}}img/badge-notebook.svg)]({{root}}{{fsdocs-source-basename}}.ipynb)

# Peptide search databases

Bottom-up proteomics identifies proteins through their peptides. The mass of an
intact protein is not specific enough to identify it, so the sample is digested
with a protease, almost always trypsin, and the mass spectrometer measures the
resulting peptides. Trypsin cuts C-terminal to lysine and arginine, except when
a proline follows, which yields peptides of around 14 residues on average that
ionize well. When the instrument selects a peptide ion for fragmentation, it
records the precursor m/z, and
[charge state determination]({{root}}01_03_charge_state_determination.html)
turns that into a neutral mass. Identification then becomes a lookup problem:
which peptides from the organism's proteome could have this mass?

A run holds tens of thousands of spectra, and each one needs its candidate
peptides with every modification variant and a precomputed mass. The `SearchDB`
module does this preparation once. It digests the FASTA proteome, generates the
modified variants of each peptide, precomputes their masses and stores
everything in a SQLite database file with an index on the mass column. A search
engine then asks for all candidates within a narrow mass window around a
measured precursor and gets them back from an indexed query.

This page builds such a database for the chloroplast proteome of
*Chlamydomonas reinhardtii* and queries it by mass. The same organism
accompanies the rest of these pages.

## Describing the search space

Everything the database will contain is decided up front by a `SearchDbParams`
record. The parameters are the identity of the database: two runs with the same
parameters mean the same database, and any change means a different one.

A fixed modification is applied to every occurrence of its target
residue, the standard example being carbamidomethylation of cysteine from the
sample preparation. A variable modification may or may not be present, so the
database stores each peptide with and without it. The most common variable
modification is oxidation of methionine, which adds one oxygen atom.
`createSearchModification` describes it: a name, a
[Unimod](https://www.unimod.org) accession, a description, a flag whether the
modification is biological in origin, the composition as an elemental formula,
the target sites, whether the composition is added or subtracted, and a short
code that will mark the modification inside stored sequence strings.
*)

open System.IO
open BioFSharp
open BioFSharp.Mz

let oxidationM =
    SearchDB.createSearchModification
        "Oxidation'Met'" "35" "Oxidation of methionine" true "O"
        [ SearchDB.Specific(AminoAcids.Met, ModificationInfo.ModLocation.Residual) ]
        SearchDB.SearchModType.Plus "ox"

let oxidationDelta =
    SearchDB.massFBy SearchDB.MassMode.Monoisotopic (SearchDB.getModBy oxidationM)

printfn "mass shift of %s: %.6f Da" oxidationM.Name oxidationDelta

(**
```text
mass shift of Oxidation'Met': 15.994915 Da
```

The composition "O" translates to the monoisotopic mass of one oxygen atom, the
expected +15.995 Da shift. With the modification in hand we can assemble the
full parameter record. `createSearchDbParams` takes every field as a positional
argument.
*)

let fastaPath = __SOURCE_DIRECTORY__ + "/data/Chlamy_Cp.fastA"

let fastaHeaderToName (header: string) = header.Split('|').[1].Trim()

let searchDbParams =
    SearchDB.createSearchDbParams
        "Chlamy_Cp_trypsin"                             // database name
        dbFolder                                        // folder for the db file
        fastaPath                                       // proteome FASTA
        fastaHeaderToName                               // header -> accession
        (Digestion.Table.getProteaseBy "Trypsin")       // protease
        0 2                                             // min/max missed cleavages
        15000.                                          // MaxMass
        4 40                                            // MinPepLength/MaxPepLength
        []                                              // isotopic mods
        SearchDB.MassMode.Monoisotopic                  // mass mode
        (SearchDB.massFBy SearchDB.MassMode.Monoisotopic) // matching mass function
        []                                              // fixed mods
        [ oxidationM ]                                  // variable mods
        2                                               // variable mod threshold

printfn "database file name: %s" (Path.GetFileName(SearchDB.Db.getNameOf searchDbParams))
printfn "protease: %s" searchDbParams.Protease.Name

(**
```text
database file name: Chlamy_Cp_trypsin.db
protease: Trypsin
```

The name and folder determine where the SQLite file lands, here a folder under the
system temp directory, and `SearchDB.Db.getNameOf` shows the resulting file name.
`fastaHeaderToName` extracts a protein accession from each FASTA header line: our
headers look like `sp|P19528| cytochrome b6/f complex subunit 4`, so splitting at `|`
and taking the second field yields the accession. The protease comes from BioFSharp's
`Digestion.Table`, either by name as shown or directly as `Digestion.Table.Trypsin`.

A missed cleavage is a lysine or arginine the protease failed to cut, and real
digests always contain such peptides, so search databases typically include peptides
with up to two or three missed cleavages, here up to two.

One implementation detail matters when choosing `MinPepLength` and `MaxPepLength`:
the digest filter compares
`CleavageEnd - CleavageStart`, which is the peptide length minus one, strictly
against both bounds. The shortest stored peptide therefore has
`MinPepLength + 2` residues and the longest has `MaxPepLength` residues. With
the values 4 and 40 above, the database contains peptides of 6 to 40 residues.
`MaxMass` is recorded with the database and participates in its identity.

Metabolic labeling such as full 15N would go into the isotopic modification list,
empty here, and the database would then store every peptide in a light and a heavy
form. `MassMode.Monoisotopic` together with `massFBy` selects the memoized
monoisotopic mass function used for all mass computations. Fixed modifications stay
empty. Oxidation of methionine goes in as a variable modification, and the threshold
of 2 caps how many variable modifications a single peptide may carry.

## Building and connecting to the database

`connectOrCreateDB` checks whether a database
file matching the parameters already exists. If it does, it simply reconnects.
If not, it reads the FASTA, digests every protein, generates the modified
variants, computes their masses and bulk inserts everything. The returned
value is an open `SQLiteConnection`. A few SQL counts show what the build
produced.
*)

open System.Data.SQLite

let cn = SearchDB.connectOrCreateDB searchDbParams

let countRows table =
    use cmd = new SQLiteCommand(sprintf "SELECT COUNT(*) FROM %s" table, cn)
    cmd.ExecuteScalar() :?> int64

printfn "database file exists: %b" (File.Exists(SearchDB.Db.getNameOf searchDbParams))
printfn "proteins:                  %i" (countRows "Protein")
printfn "distinct peptide sequences: %i" (countRows "PepSequence")
printfn "mass entries (ModSequence): %i" (countRows "ModSequence")

(**
```text
database file exists: true
proteins:                  74
distinct peptide sequences: 6415
mass entries (ModSequence): 9355
```

The 74 chloroplast proteins digest into several thousand distinct peptides,
and the variable methionine oxidation expands them into more rows in the
`ModSequence` table, which holds one row per modified variant with its
precomputed mass. Calling `connectOrCreateDB` a second time with the same
parameters finds the parameter record stored inside the file and reconnects
without rebuilding anything. Changing any parameter, even only `MaxMass`,
changes the identity and triggers a fresh build.

## Looking up candidates by precursor mass

`getThreadSafePeptideLookUpFromFileBy` prepares the mass window query against
the open connection. The result is a function taking a lower and an upper
neutral mass and returning every `ModSequence` row in between.

To query something realistic we first need a mass a spectrometer could have
measured. We digest the first protein of the FASTA in memory with the same
trypsin instance and pick the first tryptic peptide that fits the stored
length range and contains a methionine. Its neutral mass is the residue masses
summed up plus one water for the termini.
*)

open BioFSharp.IO

let firstProtein =
    Fasta.read BioArray.ofAminoAcidString fastaPath
    |> Seq.head

let targetPeptide =
    Digestion.BioArray.digest Digestion.Table.Trypsin 0 (firstProtein.Sequence |> Array.ofSeq)
    |> Array.find (fun p ->
        p.PepSequence.Length >= 6 && p.PepSequence.Length <= 40
        && List.contains AminoAcids.Met p.PepSequence)

let targetMass =
    targetPeptide.PepSequence
    |> List.sumBy BioItem.monoisoMass
    |> (+) (BioItem.monoisoMass ModificationInfo.Table.H2O)

printfn "protein:  %s" (fastaHeaderToName firstProtein.Header)
printfn "peptide:  %s" (BioList.toString targetPeptide.PepSequence)
printfn "neutral monoisotopic mass: %.5f Da" targetMass

(**
```text
protein:  P19528
peptide:  LLGVLLMAAVPAGLITVPFIESINK
neutral monoisotopic mass: 2578.51720 Da
```

A production pipeline queries the database with a 30 ppm window around the
measured precursor mass.
*)

let lookUpByMass = SearchDB.getThreadSafePeptideLookUpFromFileBy cn searchDbParams

let ppmToDalton ppm mass = mass * ppm / 1000000.

let showHits (hits: SearchDB.LookUpResult<AminoAcids.AminoAcid> list) =
    hits
    |> List.sortBy (fun h -> h.Mass)
    |> List.iter (fun h ->
        printfn "%-30s mass %11.5f  rounded %11i  pepSeqID %i  modSeqID %i  globalMod %i"
            h.StringSequence h.Mass h.RoundedMass h.PepSequenceID h.ModSequenceID h.GlobalMod)

let tolerance = ppmToDalton 30. targetMass

printfn "30 ppm at this mass: %.5f Da" tolerance

let hits = lookUpByMass (targetMass - tolerance) (targetMass + tolerance)

showHits hits

(**
```text
30 ppm at this mass: 0.07736 Da
LLGVLLMAAVPAGLITVPFIESINK      mass  2578.51720  rounded  2578517204  pepSeqID 4  modSeqID 3  globalMod 0
```

The window returns exactly the peptide we computed the mass for, as a
`LookUpResult`. `StringSequence` is the stored sequence, `Mass` the
precomputed neutral mass and `RoundedMass` that mass multiplied by one million
and stored as an integer, which is the indexed column the BETWEEN query runs
against. `PepSequenceID` identifies the plain peptide sequence and
`ModSequenceID` the specific modified variant. `GlobalMod` tells whether the
entry belongs to the isotopically labeled form of the database. We configured
no isotopic modification, so it is 0 for every entry. The `BioSequence` field,
not printed here, holds the sequence parsed back into a BioFSharp amino acid
list, ready for fragment prediction.

The oxidized form of the same peptide weighs one oxygen more and sits in its
own mass window. Querying 30 ppm around that mass returns the variant carrying
the modification.
*)

let oxidizedMass = targetMass + oxidationDelta

let oxHits = lookUpByMass (oxidizedMass - tolerance) (oxidizedMass + tolerance)

showHits oxHits

(**
```text
LLGVLL[ox]MAAVPAGLITVPFIESINK  mass  2594.51212  rounded  2594512119  pepSeqID 4  modSeqID 4  globalMod 0
```

The modification appears inside the sequence string as the `[ox]` code
directly in front of the modified residue, the same code we passed to
`createSearchModification`. Both variants share the `PepSequenceID` of the
plain sequence but have distinct `ModSequenceID`s and masses. A search engine
scoring this window would now treat the oxidized sequence as its own
candidate.

## Mapping peptides back to proteins

After spectra have been matched, the pipeline needs to know which proteins a
peptide belongs to, because the final goal is a protein list. The
`CleavageIndex` table stores this mapping, and
`getProteinPeptideLookUpFromFileBy` prepares a lookup from a `PepSequenceID`
to the accessions of all proteins containing that peptide, each accession
paired with the peptide sequence the ID stands for. It expects an
in-memory copy of the database, which `copyDBIntoMemory` produces from the
open file connection, so the many lookups of a full run avoid disk access.
*)

let memoryDB = SearchDB.copyDBIntoMemory cn

let proteinsOfPeptide = SearchDB.getProteinPeptideLookUpFromFileBy memoryDB

let firstHit = hits |> List.minBy (fun h -> h.Mass)

proteinsOfPeptide firstHit.PepSequenceID
|> List.iter (fun (accession, peptideSequence) ->
    printfn "peptide %s occurs in protein %s" peptideSequence accession)

(**
```text
peptide LLGVLLMAAVPAGLITVPFIESINK occurs in protein P19528
```

The peptide maps back to P19528, the cytochrome b6/f subunit the FASTA starts
with, and the returned accession is exactly what our `fastaHeaderToName`
function extracted during the build. This reverse mapping is the raw material
for [protein inference]({{root}}04_02_protein_inference.html), where shared
peptides make the mapping ambiguous and need to be resolved.

## Caching repeated lookups

Consecutive precursors in a run often have similar masses, so their candidate
windows overlap and the same peptides would be fetched and their fragments
predicted again. The `Cache` module addresses this with a thin wrapper around
`SortedList`, keyed by the same rounded integer masses the database uses.
`getPeptideLookUpWithMemBy` combines such a cache with the database lookup and
fragment prediction, so a peptide that already went through the pipeline is
served from memory. Creating a cache and storing a result under its rounded
mass looks like this.
*)

let lookUpCache = Cache.createCache<int64, SearchDB.LookUpResult<AminoAcids.AminoAcid> list>

Cache.addItem lookUpCache (firstHit.RoundedMass, hits)

printfn "cached entries: %i" lookUpCache.Count
printfn "contains key %i: %b" firstHit.RoundedMass (fst (Cache.getItemBy lookUpCache firstHit.RoundedMass))

(**
```text
cached entries: 1
contains key 2578517204: true
```

`addItem` inserts or replaces the value under a key and `getItemBy` retrieves it.
The memoized lookup itself is wired up by the search engine.

The next step is to predict each candidate's fragments and compare them to the
measured spectrum, the subject of
[SEQUEST-like scoring]({{root}}02_03_sequest_like_scoring.html).
*)
