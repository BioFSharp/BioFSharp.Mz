(**
---
title: Protein inference
category: Validation and inference
categoryindex: 4
index: 2
---
*)

(*** hide ***)

(*** condition: prepare ***)
#r "nuget: FSharpAux, 2.1.0"
#r "nuget: FSharpAux.IO, 2.1.0"
#r "nuget: FSharp.Stats, 0.6.0"
#r "nuget: BioFSharp, 2.0.0"
#r "nuget: Plotly.NET, 6.0.0-preview.2"
#r "../src/BioFSharp.Mz/bin/Release/netstandard2.0/BioFSharp.Mz.dll"
#r "../src/BioFSharp.Mz.Vis/bin/Release/netstandard2.0/BioFSharp.Mz.Vis.dll"

(*** condition: ipynb ***)
#if IPYNB
#r "nuget: BioFSharp.Mz, {{fsdocs-package-version}}"
#r "nuget: BioFSharp.Mz.Vis, {{fsdocs-package-version}}"
#endif // IPYNB

(**
[![Binder]({{root}}img/badge-binder.svg)](https://mybinder.org/v2/gh/BioFSharp/BioFSharp.Mz/gh-pages?filepath={{fsdocs-source-basename}}.ipynb)&emsp;
[![Script]({{root}}img/badge-script.svg)]({{root}}{{fsdocs-source-basename}}.fsx)&emsp;
[![Notebook]({{root}}img/badge-notebook.svg)]({{root}}{{fsdocs-source-basename}}.ipynb)

# Protein inference

The [FDR control]({{root}}04_01_fdr_control.html) page ended with a
statistically controlled list of accepted peptide identifications.
Biological questions are asked about proteins, and mapping the accepted
peptides back onto proteins runs into an ambiguity. A tryptic peptide often
occurs in more than one protein of the database. Splice variants of a gene
share most of their coding sequence, and gene families carry conserved
stretches that digest into identical peptides. When such a peptide is
observed, the evidence supports every protein that contains it, and nothing
in the measurement can tell those proteins apart. Protein inference resolves
this by reporting protein groups: sets of accessions that the observed
peptides cannot distinguish, carried through the rest of the pipeline as one
entry.

The `ProteinInference` module describes each observed peptide as a
`ProteinClassItem`, which couples the peptide sequence with the array of
protein accessions it maps to and a `PeptideEvidenceClass`. The class states
how specifically the peptide points at a gene model. `C1a` marks a peptide
that maps to exactly one splice variant of one gene. `C1b` marks one that
maps to several splice variants which translate to the same protein
sequence. `C2a` and `C2b` cover peptides shared between splice variants that
differ elsewhere, either a subset of them or all of them. `C3a` and `C3b`
cover peptides shared between different genes, where `C3a` is the case of
genes encoding an identical protein sequence and `C3b` the general case.
Inference processes the classes in that order, and in the output every group
carries the class of the evidence that created it.

In a full pipeline these classes come from genome annotation. The GFF3
helpers of the module read an annotation file: `assignTranscriptsToGenes`
collects the mRNA entries under their gene locus, with
`createProteinModelInfoFromEntry` extracting the splice variant identifier
of each transcript, and `createPeptideProteinRelation` maps the peptides of
an in silico digest onto the resulting gene models. From that relation the
evidence class of every peptide follows from how it distributes over loci
and splice variants. This page skips the annotation step and constructs the
class items directly, which keeps the example self-contained.

## Describing the peptide evidence

The example imitates a small slice of a Chlamydomonas reinhardtii
measurement, built so that every interesting mapping case occurs at least
once. Two chloroplast proteins have only unique peptides. A light-harvesting
gene has two annotated splice variants and every observed peptide sits in
their shared exons, so those peptides map to both variants. Two actin-like
paralogs each have a unique peptide and additionally share one conserved
peptide. One protein is a one-hit wonder with a single weakly scoring
peptide, which becomes important once decoy scores enter.
*)

open BioFSharp.Mz
open BioFSharp.PeptideClassification

let items =
    [
        // psbA (photosystem II D1) and rbcL (rubisco large subunit): two
        // unique peptides each, mapping to exactly one gene model (C1a).
        ProteinInference.createProteinClassItem [|"psbA"|] PeptideEvidenceClass.C1a "VLNTWADIINR"
        ProteinInference.createProteinClassItem [|"psbA"|] PeptideEvidenceClass.C1a "ETTENESANEGYR"
        ProteinInference.createProteinClassItem [|"rbcL"|] PeptideEvidenceClass.C1a "DTDILAAFR"
        ProteinInference.createProteinClassItem [|"rbcL"|] PeptideEvidenceClass.C1a "LTYYTPDYVVR"
        // A light-harvesting gene with two splice variants. Both observed
        // peptides lie in shared exons and map to all isoforms of the gene
        // (C2b), so the two variants are indistinguishable here.
        ProteinInference.createProteinClassItem
            [|"Cre01.g001000.t1"; "Cre01.g001000.t2"|] PeptideEvidenceClass.C2b "SSYLDGSAPGDFGFDPLGLGK"
        ProteinInference.createProteinClassItem
            [|"Cre01.g001000.t1"; "Cre01.g001000.t2"|] PeptideEvidenceClass.C2b "NLAGDVIGTR"
        // Two actin-like paralogs: each has one unique peptide (C1a) ...
        ProteinInference.createProteinClassItem [|"Cre02.g095150.t1"|] PeptideEvidenceClass.C1a "SYELPDGQVITIGNER"
        ProteinInference.createProteinClassItem [|"Cre03.g144807.t1"|] PeptideEvidenceClass.C1a "VAPEEHPVLLTEAPLNPK"
        // ... and they share one conserved peptide across the two gene loci
        // (C3b, different genes with different protein sequences).
        ProteinInference.createProteinClassItem
            [|"Cre02.g095150.t1"; "Cre03.g144807.t1"|] PeptideEvidenceClass.C3b "AGFAGDDAPR"
        // A one-hit wonder: a single peptide, and its PSM score will be low.
        ProteinInference.createProteinClassItem [|"Cre06.g250100.t1"|] PeptideEvidenceClass.C1a "LISWYDNEWGYSNR"
    ]

let distinctProteins =
    items
    |> List.collect (fun item -> List.ofArray item.GroupOfProteinIDs)
    |> List.distinct

printfn "%i peptides mapping onto %i proteins" items.Length distinctProteins.Length

(**
```text
10 peptides mapping onto 7 proteins
```

## Grouping the evidence into protein groups

`inferSequences` turns the class items into `InferredProteinClassItem` values, one
per protein group, each recording the group's evidence class and the peptides kept
for its quantification. The `IntegrationStrictness` decides what happens to
overlapping groups: `Minimal` merges them down to the smallest set of proteins that
still explains all observed peptides, while `Maximal` keeps every evidence group
exactly as the peptides stated it. The `PeptideUsageForQuantification` decides which
peptides a finished group keeps: under `Minimal` only those, among all peptides
consistent with it, whose own protein set is smallest and therefore points at the
group most specifically. Under `Maximal` it keeps every peptide whose protein set
contains the group, and under `MaximalInverse` every peptide whose protein set is
contained in it. The production settings of this module are `Maximal` strictness
with `Minimal` peptide usage, and the call writes a progress line to standard output
while it assembles its lookup structures.
*)

let inferred =
    ProteinInference.inferSequences
        ProteinInference.IntegrationStrictness.Maximal
        ProteinInference.PeptideUsageForQuantification.Minimal
        items
    |> List.ofSeq

(**
```text
finish setup
```

That progress marker is everything the call itself prints under `Maximal`
strictness. The result prints as one line per protein group.
`proteinGroupToString` joins the accessions of a group with semicolons, and
`inferSequences` already returns that joined representation.
*)

for g in inferred do
    printfn "%-33s %s  %s"
        g.GroupOfProteinIDs
        (sprintf "%A" g.Class)
        (String.concat ", " g.PeptideSequence)

(**
```text
psbA                              C1a  VLNTWADIINR, ETTENESANEGYR
rbcL                              C1a  DTDILAAFR, LTYYTPDYVVR
Cre02.g095150.t1                  C1a  SYELPDGQVITIGNER
Cre03.g144807.t1                  C1a  VAPEEHPVLLTEAPLNPK
Cre06.g250100.t1                  C1a  LISWYDNEWGYSNR
Cre01.g001000.t1;Cre01.g001000.t2 C2b  SSYLDGSAPGDFGFDPLGLGK, NLAGDVIGTR
Cre02.g095150.t1;Cre03.g144807.t1 C3b  AGFAGDDAPR
```

The two splice variants collapsed into a single group, because no observed
peptide separates them. The actin-like paralogs stayed apart, since each has
unique evidence, and their shared peptide formed a composite group of its
own. `Maximal` strictness meant that every evidence group was reported
unchanged, seven groups in total. The `Minimal` peptide usage shows in the
paralog groups: the shared peptide `AGFAGDDAPR` was consistent with the
`Cre02` group as well, but only the unique peptide points at that group
most specifically, so only the unique peptide is kept for quantifying it.

For comparison, the same evidence under `Minimal` strictness. On this path
the merge logic runs and traces its decisions to standard output, so the
call prints a diagnostic line for every piece of evidence that overlaps
nothing recorded so far and echoes the existing group when further evidence
arrives on the same proteins.
*)

let inferredMinimal =
    ProteinInference.inferSequences
        ProteinInference.IntegrationStrictness.Minimal
        ProteinInference.PeptideUsageForQuantification.Minimal
        items
    |> List.ofSeq

(**
```text
Search did not return anything
psbA
Search did not return anything
rbcL
Search did not return anything
Search did not return anything
Search did not return anything
Search did not return anything
Cre01.g001000.t1;Cre01.g001000.t2
finish setup
```

Every group is announced once when its first evidence arrives, and the
echoed group names accompany the integration of the second peptide on the
same proteins. The last line is the same progress marker as above, and the
comparison against the `Maximal` result comes from the returned values.
*)

let droppedGroups =
    inferred
    |> List.filter (fun g ->
        inferredMinimal
        |> List.forall (fun m -> m.GroupOfProteinIDs <> g.GroupOfProteinIDs))
    |> List.map (fun g -> g.GroupOfProteinIDs)

printfn "Maximal: %i groups, Minimal: %i groups" inferred.Length inferredMinimal.Length
printfn "dropped under Minimal: %s" (String.concat ", " droppedGroups)

(**
```text
Maximal: 7 groups, Minimal: 6 groups
dropped under Minimal: Cre02.g095150.t1;Cre03.g144807.t1
```

Under `Minimal` strictness the composite paralog group disappears, because
the two paralogs already explain the shared peptide on their own. The other
six groups come out identical under both settings.

## Scoring the groups against reversed proteins

A protein list needs the same error control as the peptide list before it,
again through a score compared against a decoy. The module
scores a protein group by summing the PSM scores of its peptides, and scores
its decoy counterpart from the peptides of the reversed protein sequences
that the search also matched. `PSMInput` is the row shape the module reads
from PSM result files: the peptide sequence and score of one match, plus
the database identifier of the peptide. Here the rows are built in memory,
one list per MS run, with deterministic scores on an Andromeda-like scale
where confident matches reach high double digits and random matches stay
low.
*)

let psm id sequence score : ProteinInference.PSMInput =
    { PepSequenceID = id; Seq = sequence; Score = score }

let run1 =
    [
        // Accepted target peptides. The one-hit wonder LISWYDNEWGYSNR
        // scores weakly by design.
        "VLNTWADIINR",           88.7
        "ETTENESANEGYR",         74.5
        "DTDILAAFR",            102.8
        "LTYYTPDYVVR",           68.9
        "SSYLDGSAPGDFGFDPLGLGK", 83.6
        "NLAGDVIGTR",            57.2
        "SYELPDGQVITIGNER",      79.4
        "VAPEEHPVLLTEAPLNPK",    88.1
        "AGFAGDDAPR",            61.7
        "LISWYDNEWGYSNR",        45.3
        // Decoy matches: some spectra matched peptides of the reversed
        // protein sequences. Their scores sit in the random-match range.
        "NIIDAWTNLVR",           18.3
        "VVYDPTYYTLR",           24.6
        "TGIVDGALNR",            21.4
        "ENGITIVQGDPLEYSR",      15.8
        "PNLPAETLLVPHEEPAVK",    27.5
        "NSYGWENDYWSILR",        26.7
        "GLTSAFVEK",             21.5
    ]
    |> List.mapi (fun i (sequence, score) -> psm i sequence score)

let run2 =
    [
        // A second injection rescored two peptides.
        "VLNTWADIINR", 96.2
        "DTDILAAFR",   95.1
    ]
    |> List.mapi (fun i (sequence, score) -> psm (100 + i) sequence score)

let peptideScores = ProteinInference.createPeptideScoreMap [run1; run2]

printfn "best score for VLNTWADIINR: %.1f" peptideScores.["VLNTWADIINR"]

(**
```text
best score for VLNTWADIINR: 96.2
```

`createPeptideScoreMap` flattens the runs and keeps the best score seen for
each peptide sequence, so the rescored `VLNTWADIINR` carries 96.2 from the
second run while `DTDILAAFR` keeps its 102.8 from the first.

The decoy side starts from the reverse digest: for every accession, the
peptides that its reversed sequence produces. The peptides follow the
reverse-and-keep-cleavage-site convention of common search engines, and the
score map above already contains the low scores they earned.
`createReverseProteinScores` sums, per reversed protein, the scores of those
of its peptides that were actually observed. The reversed one-hit wonder
picked up two such matches, and their sum is about to matter.
*)

let reversedProteins =
    [|
        "psbA",             [|"NIIDAWTNLVR"|]
        "rbcL",             [|"VVYDPTYYTLR"|]
        // Both splice variants produce the same reversed peptide.
        "Cre01.g001000.t1", [|"TGIVDGALNR"|]
        "Cre01.g001000.t2", [|"TGIVDGALNR"|]
        "Cre02.g095150.t1", [|"ENGITIVQGDPLEYSR"|]
        "Cre03.g144807.t1", [|"PNLPAETLLVPHEEPAVK"|]
        "Cre06.g250100.t1", [|"NSYGWENDYWSILR"; "GLTSAFVEK"|]
    |]

let decoyScores = ProteinInference.createReverseProteinScores reversedProteins peptideScores

(**
With both maps in place, each inferred group receives its two scores.
`assignPeptideScores` sums the scores of the peptides the group kept for
quantification, and `assignDecoyScoreToTargetScore` takes the best decoy
score among the group's accessions. `createInferredProteinClassItemScored`
assembles the record the downstream estimators work with: these entries
describe target proteins found in the search database, and the
`DecoyHasBetterScore` flag records the picked target-decoy comparison of
the two scores.
*)

let scored =
    inferred
    |> List.map (fun g ->
        let targetScore = ProteinInference.assignPeptideScores g.PeptideSequence peptideScores
        let decoyScore  = ProteinInference.assignDecoyScoreToTargetScore g.GroupOfProteinIDs decoyScores
        ProteinInference.createInferredProteinClassItemScored
            g.GroupOfProteinIDs
            g.Class
            g.PeptideSequence
            targetScore
            decoyScore
            false                       // these are target-side entries
            (decoyScore > targetScore)  // the picked target-decoy comparison
            true)                       // every accession is in the search DB
    |> Array.ofList

for s in scored do
    printfn "%-33s target %5.1f  decoy %4.1f  decoy wins %b"
        s.GroupOfProteinIDs s.TargetScore s.DecoyScore s.DecoyHasBetterScore

(**
```text
psbA                              target 170.7  decoy 18.3  decoy wins false
rbcL                              target 171.7  decoy 24.6  decoy wins false
Cre02.g095150.t1                  target  79.4  decoy 15.8  decoy wins false
Cre03.g144807.t1                  target  88.1  decoy 27.5  decoy wins false
Cre06.g250100.t1                  target  45.3  decoy 48.2  decoy wins true
Cre01.g001000.t1;Cre01.g001000.t2 target 140.8  decoy 21.4  decoy wins false
Cre02.g095150.t1;Cre03.g144807.t1 target  61.7  decoy 27.5  decoy wins false
```

Six groups beat their decoy counterpart clearly. The one-hit wonder loses:
its single peptide scored 45.3 while the two random matches on its reversed
sequence sum to 48.2, so the group is flagged as a decoy win. This is the
protein-level analogue of a wrong PSM, and having such a case in the set is
what gives the error estimators of the next section something to measure.

## Estimating the protein-level FDR

`calculateFDRwithDecoyTargetRatio` is the direct estimate: the number of
decoy wins divided by the number of target wins, the protein-level version
of the decoy counting from the [FDR control]({{root}}04_01_fdr_control.html)
page.
*)

let ratioFDR = ProteinInference.calculateFDRwithDecoyTargetRatio scored

printfn "decoy/target ratio FDR: %.4f" ratioFDR

(**
```text
decoy/target ratio FDR: 0.1667
```

One decoy win against six target wins gives an estimated FDR of about 17
percent. The number is high because the set is tiny, and it says that
roughly one of the six reported groups could be wrong at this acceptance.

`calculateFDRwithMAYU` refines the estimate with the MAYU model introduced
in the MAYU section of the [FDR control]({{root}}04_01_fdr_control.html)
page. It needs the search database as accession and sequence pairs, because
the model stratifies proteins by sequence length before applying the
hypergeometric expectation per length bin, and unidentified database
proteins enter those bins too. Only the length of each sequence is used, so
placeholder residue strings of realistic lengths stand in for the real
sequences here.
*)

let proteinsFromDB =
    [|
        // The seven accessions of the example, with realistic lengths.
        "psbA",             String.replicate 353 "A"
        "rbcL",             String.replicate 475 "A"
        "Cre01.g001000.t1", String.replicate 254 "A"
        "Cre01.g001000.t2", String.replicate 231 "A"
        "Cre02.g095150.t1", String.replicate 380 "A"
        "Cre03.g144807.t1", String.replicate 377 "A"
        "Cre06.g250100.t1", String.replicate 335 "A"
        // Database proteins that were never identified in this run.
        "Cre09.g396920.t1", String.replicate 497 "A"
        "Cre12.g554400.t1", String.replicate 155 "A"
        "Cre16.g678000.t1", String.replicate 612 "A"
    |]

let mayuFDR = ProteinInference.calculateFDRwithMAYU scored proteinsFromDB

printfn "MAYU FDR: %.4f" mayuFDR

(**
```text
MAYU FDR: 0.1111
```

MAYU lands somewhat below the raw ratio. It converts the decoy observation
in its length bin into an expected number of false positives among the
targets, and that expectation stays below one here, so the decoy win enters
the estimate with less than the full weight the ratio gave it.

## Attaching protein q-values

For the final list every group receives a q-value, computed from the same
score axis the estimators used. `FDRControl.calculateQValueStorey` from the
[FDR control]({{root}}04_01_fdr_control.html) page works on any record type
once it is told how to read the data, so it applies to the scored groups
directly. The decoy flag is `DecoyHasBetterScore`, and with the two score
selectors a decoy win enters the estimate at its decoy score while a
target win enters at its target score. The returned function
maps a score to its q-value, and `assignQValueToIPCIS` applies it to each
record, wrapping the result in an `InferredProteinClassItemQValue`.
*)

let proteinQValue =
    FDRControl.calculateQValueStorey
        scored
        (fun s -> s.DecoyHasBetterScore)
        (fun s -> s.DecoyScore)
        (fun s -> s.TargetScore)

let qValued =
    scored
    |> Array.map (ProteinInference.assignQValueToIPCIS proteinQValue)

for q in qValued do
    printfn "%-33s score %5.1f  q %6.4f"
        q.InfProtClassItem.GroupOfProteinIDs
        q.InfProtClassItem.TargetScore
        q.QValue

(**
```text
psbA                              score 170.7  q 0.0000
rbcL                              score 171.7  q 0.0000
Cre02.g095150.t1                  score  79.4  q 0.0000
Cre03.g144807.t1                  score  88.1  q 0.0000
Cre06.g250100.t1                  score  45.3  q 0.2025
Cre01.g001000.t1;Cre01.g001000.t2 score 140.8  q 0.0000
Cre02.g095150.t1;Cre03.g144807.t1 score  61.7  q 0.0000
```

Every group that beat its decoy counterpart sits above the score of the
single decoy win, so accepting any of them accumulates no decoy evidence
and their q-values are zero. The one-hit wonder sits at the bottom of the
score axis and collects the worst q-value of the run. Filtering this list
at the usual one percent threshold keeps the six confident groups and
removes the one-hit wonder.

## Rendering the inference report with BioFSharp.Mz.Vis

`BioFSharp.Mz.Vis` is the second library of this repository. It holds the
visual reporting that pipeline runs produce alongside their result tables,
built on Plotly.NET, and it ships as its own assembly so that the core
library carries no charting dependency.

Its `ProteinInference.qValueHitsVisualization` renders the q-value results
of this page into an HTML report. The arguments are the histogram bandwidth
on the score axis, the array of q-value records, the output path stem and a
flag for grouped output. With `groupFiles = false` the function appends
`_QValueGraph` to the path stem and Plotly adds the `.html` extension. A
bandwidth of 10 suits the score axis of this example, which spans from the
forties up to about 170.
*)

let reportStem =
    System.IO.Path.Combine(System.IO.Path.GetTempPath(), "chlamy_protein_inference")

BioFSharp.Mz.Vis.ProteinInference.qValueHitsVisualization 10.0 qValued reportStem false

let reportPath = reportStem + "_QValueGraph.html"

printfn "report exists: %b" (System.IO.File.Exists reportPath)
printfn "report size: %i KB" (System.IO.FileInfo(reportPath).Length / 1024L)

(**
```text
report exists: true
report size: 11 KB
```

The report is a single Plotly figure. Column traces show the score
histograms of the target groups and of the decoy wins, drawn against a
relative frequency axis on the left and accompanied by an absolute count
axis on the right. A scatter trace on the same score axis draws each
group's q-value over its score, which on real data traces the rising error
curve toward low scores. In this small run the separation is easy to read:
the decoy-win bar sits at the low end of the score axis below every
confident target bar, and the q-value points are flat at zero everywhere
above it. The HTML loads plotly.js from a CDN, so viewing the report needs
network access.

## The end of the identification pipeline

The q-valued protein list is the final product of the identification side of the
pipeline. What remains is connecting the identities to amounts: the ion areas
determined by
[quantification]({{root}}03_01_quantification.html) are aggregated over exactly
these protein groups, using the peptides each group kept for quantification.
*)
