(**
---
title: False discovery rate control
category: Validation and inference
categoryindex: 4
index: 1
---
*)

(*** hide ***)

(*** condition: prepare ***)
#r "nuget: FSharpAux, 2.1.0"
#r "nuget: FSharpAux.IO, 2.1.0"
#r "nuget: FSharp.Stats, 0.6.0"
#r "nuget: BioFSharp, 2.0.0"
#r "../src/BioFSharp.Mz/bin/Release/netstandard2.0/BioFSharp.Mz.dll"

(*** condition: ipynb ***)
#if IPYNB
#r "nuget: BioFSharp.Mz, {{fsdocs-package-version}}"
#endif // IPYNB

(**
[![Binder]({{root}}img/badge-binder.svg)](https://mybinder.org/v2/gh/BioFSharp/BioFSharp.Mz/gh-pages?filepath={{fsdocs-source-basename}}.ipynb)&emsp;
[![Script]({{root}}img/badge-script.svg)]({{root}}{{fsdocs-source-basename}}.fsx)&emsp;
[![Notebook]({{root}}img/badge-notebook.svg)]({{root}}{{fsdocs-source-basename}}.ipynb)

# False discovery rate control

Over a whole run, the per-spectrum rankings of the scoring pages
([SEQUEST-like]({{root}}02_03_sequest_like_scoring.html) and
[Andromeda-like]({{root}}02_04_andromeda_like_scoring.html)) accumulate into a long
list of peptide spectrum matches (PSMs), each carrying a score.
The score states how well a spectrum and a candidate agree. It says nothing
about how often agreement this good arises by chance, so a score alone
cannot separate true identifications from false ones. The search database
never covers every peptide actually present in the sample, so for some
spectra the true answer was not among the candidates and the best available
match is wrong. The instrument also sometimes selects precursor ions that
were never peptides. And in a large search space, some wrong candidate will
score well on some spectrum simply because so many candidates were tried.

The standard way to put numbers on this problem is the decoy database.
Reversed or shuffled sequences are searched alongside the real ones, and any
PSM that lands on a decoy is wrong by construction. Counting decoy hits
above a score threshold therefore estimates how many target hits above the
same threshold are wrong, because a wrong target assignment behaves like a
decoy assignment.

Three quantities are built on that count. The global false discovery rate
(FDR) of an accepted PSM list is the fraction of the list that is wrong. The
q-value attached to a score is the FDR of the list you get by accepting
everything at that score or better, made monotone so that relaxing the
threshold never claims a lower error. It is the number used for cutoffs, and
"filter at 1 percent FDR" means keeping PSMs with q at or below 0.01. The
posterior error probability (PEP) asks the question for a single PSM: how
likely is this one match wrong. Tools in this space include Qvality, which
estimates PEPs and q-values by regression on decoy counts, and Percolator,
which rescores PSMs with an iteratively trained SVM and appears again at the
end of this page. The `FDRControl` module implements the counting and
fitting machinery for q-values.

## Simulating a run's worth of PSM scores

The scoring pages worked on one spectrum, so for this page we simulate the
score list of a whole run. A seeded generator and a Box-Muller helper draw
normally distributed scores on an Andromeda-like scale from 0 to about 120.
*)

open BioFSharp.Mz

let rng = System.Random 42

// Box-Muller transform: two uniform draws become one normal draw.
let gaussian mean sd =
    let u1 = 1.0 - rng.NextDouble()
    let u2 = rng.NextDouble()
    mean + sd * (sqrt (-2.0 * log u1) * cos (2.0 * System.Math.PI * u2))

// Decoy PSMs are wrong by construction and score low.
let decoyScores = Array.init 5000 (fun _ -> max 0.0 (gaussian 25.0 8.0))

// Target PSMs are a mixture. 60 percent are correct identifications and
// score high, 40 percent are wrong assignments. The wrong targets draw from
// the same low distribution as the decoys, because that equivalence is the
// premise of the method: the decoy score distribution estimates how wrong
// target assignments score.
let correctScores = Array.init 3000 (fun _ -> max 0.0 (gaussian 70.0 15.0))
let wrongScores   = Array.init 2000 (fun _ -> max 0.0 (gaussian 25.0 8.0))
let targetScores  = Array.append correctScores wrongScores

let psms =
    Array.append
        (targetScores |> Array.map (fun s -> FDRControl.createQValueInput s false))
        (decoyScores  |> Array.map (fun s -> FDRControl.createQValueInput s true))

let median xs =
    let sorted = Array.sort xs
    sorted.[Array.length sorted / 2]

printfn "PSMs: %i targets, %i decoys" targetScores.Length decoyScores.Length
printfn "median target score: %.1f  median decoy score: %.1f"
    (median targetScores) (median decoyScores)

(**
```text
PSMs: 5000 targets, 5000 decoys
median target score: 56.4  median decoy score: 25.1
```

`createQValueInput` builds the `QValueInput` record the module works with, a
score plus an `IsDecoy` flag. The functions below take the data together
with three accessors, one for the decoy flag and one score projection each
for decoys and targets. Here both projections read the same field.

## Computing q-values by counting decoys

`calculateQValueStorey` implements the direct counting estimate and returns
a function from score to q-value. Internally it walks the scores from best
to worst and accumulates decoy and target counts. At each score the raw
estimate is accumulated decoys divided by accumulated targets. A second pass
makes the sequence monotone, so q never decreases as the score drops, and a
linear interpolation over the resulting support turns the table into a
function usable at any score.
*)

let storeyQ =
    FDRControl.calculateQValueStorey
        psms
        (fun x -> x.IsDecoy)
        (fun x -> x.Score)
        (fun x -> x.Score)

let probes = [90.; 70.; 55.; 40.]

for s in probes do
    printfn "score %5.0f  q = %.5f" s (storeyQ s)

(**
```text
score    90  q = 0.00000
score    70  q = 0.00000
score    55  q = 0.00000
score    40  q = 0.05299
```

No decoy in this simulated run scores above 55, so down to there the
accumulated decoy count is zero and so is q. By 40 the upper tail of the
decoy distribution is contributing steadily and q has passed 5 percent.
Accepting at a q threshold is a score cutoff in disguise, so counting the
targets that pass is one filter away.
*)

let acceptedTargets threshold =
    targetScores |> Array.filter (fun s -> storeyQ s <= threshold) |> Array.length

printfn "targets accepted at q <= 0.01: %i" (acceptedTargets 0.01)
printfn "targets accepted at q <= 0.05: %i" (acceptedTargets 0.05)
printfn "wrong targets among the q <= 0.01 list: %i"
    (wrongScores |> Array.filter (fun s -> storeyQ s <= 0.01) |> Array.length)

(**
```text
targets accepted at q <= 0.01: 2897
targets accepted at q <= 0.05: 3003
wrong targets among the q <= 0.01 list: 6
```

Against ground truth, the 1 percent list holds 2897 PSMs of which 6 are actually
wrong, a true error rate of 0.2 percent. The estimate is conservative here because all 5000
decoys stand in for the 2000 wrong targets, so the decoy count overstates
the false hits. The assumed fraction of wrong targets is exactly what the
pi0 parameter of the next method makes explicit.

## Smoothing the estimate with a logistic fit

The counting estimate is a step function shaped by every random fluctuation
in the decoy tail. `calculateQValueLogReg` trades the steps for a smooth
curve. First `binningFunction` slices the score axis into bandwidth-wide
bins and computes raw estimates per bin. The decoy count is normalized by
the overall decoy share, then doubled and scaled by pi0. Dividing that by the
total count in the bin gives a PEP-like value, and dividing it by the total
count at this bin and above gives a q-like value. These raw bin values serve
purely as input points for the fit. A descending logistic curve is then fitted through the
(score, q) points with constrained Levenberg-Marquardt over a range of
steepness guesses, and the fitted curve is the deliverable, a smooth
function from score to q.

The first argument is pi0, the assumed fraction of wrong assignments among
the targets. We pass 0.4 because the generated mixture is 40 percent wrong,
which is what a real analysis would estimate from the score overlap. The
bandwidth of 3 gives bins a few score units wide, enough entries per bin on
a 0 to 120 scale to keep the raw points stable.
*)

let logRegQ =
    FDRControl.calculateQValueLogReg
        0.4
        3.0
        psms
        (fun x -> x.IsDecoy)
        (fun x -> x.Score)
        (fun x -> x.Score)

printfn "%5s %10s %10s" "score" "Storey q" "logReg q"
for s in probes do
    printfn "%5.0f %10.5f %10.5f" s (storeyQ s) (logRegQ s)

(**
```text
score   Storey q   logReg q
   90    0.00000    0.00000
   70    0.00000    0.00000
   55    0.00000    0.00000
   40    0.05299    0.00065
```

The methods agree over most of the range. Everywhere above 55 both report
zero, and both place the onset of error below 45. They part ways in the
transition zone: at a score of 40 the counting estimate already sees the
159 decoys above that score, while the fitted sigmoid centers its drop on
the bulk of the decoy distribution near 30 and has fallen close to zero by
40. The pi0 scaling also lowers the fitted curve as a whole, since the fit
treats only 40 percent of the targets as potentially wrong where Storey's
count treats every decoy as evidence against a target. Which output feeds
the pipeline is a choice between the assumption-free step function and the
smooth curve that needs a pi0.

## From list-level q to single-PSM PEP

The q-value describes a whole accepted list. The PEP asks how likely one particular
PSM is wrong. A PSM
sitting exactly at a 1 percent q cutoff typically carries a PEP well above 1
percent, because the list average is pulled down by all the confident hits
above it. The library derives PEP-like values inside the binning and fitting
machinery you just saw, where each bin's decoy share is a local error
estimate, and it also exposes dedicated PEP functions: `calculatePEPValues`
computes the per-bin decoy shares directly and `initCalculateLin` returns a
fitted score-to-PEP mapping. The pipeline-facing outputs on this page are the
two q-value functions.

## Estimating protein-level false positives with MAYU

Once PSMs are aggregated to proteins, the error question returns one level
up: how many of the target proteins with at least one accepted hit are
false? `FDRControl.MAYU.estimatePi0HG` answers with the MAYU model,
translated from Percolator's ProteinFDREstimator. It treats the decoy
protein hits as draws from a hypergeometric distribution over the database
and returns the expected number of false positives among the target hits.
The arguments are the total number of candidate proteins in the database,
the number of target protein hits and the number of decoy protein hits.
*)

let expectedFP = FDRControl.MAYU.estimatePi0HG 2000. 800. 30.

printfn "expected false positive protein hits: %.2f" expectedFP
printfn "implied protein-level FDR: %.4f" (expectedFP / 800.)

(**
```text
expected false positive protein hits: 18.27
implied protein-level FDR: 0.0228
```

With 800 of 2000 database proteins hit and 30 decoy hits observed, the model
expects about 18 of the 800 target hits to be false, a protein-level FDR
just over 2 percent. The expectation is a weighted average over all feasible
false positive counts from zero up to the decoy count, so it can land below
the raw decoy count, as it does here. The [protein inference]({{root}}04_02_protein_inference.html)
page applies this estimator when assembling peptide evidence into protein
lists.

## Rescoring with Percolator

The methods above interpret one score column. Percolator improves the
separation before the interpretation: it trains a support vector machine on
the PSMs themselves, using the decoy PSMs as the negative training set and
high-confidence target PSMs as positives, on features such as the raw
score, the delta to the next-best candidate, the precursor mass error and
the precursor charge. The learned score replaces the raw one and the
positive training set is re-derived from it, and this classification round
repeats until the score stabilizes. The output is again q-values and PEPs,
computed on a score axis where targets and decoys separate better.

Percolator ships as an external program, and the `PercolatorWrapper` module
wraps its command line. The option groups of the CLI are mirrored by
discriminated unions under `PercolatorWrapper.Parameters`, so a parameter
set is ordinary typed F# data, and the `PercolatorParams` cases collect the
groups. Building a set looks like this.
*)

open BioFSharp.Mz.PercolatorWrapper

let percolatorParams : Parameters.PercolatorParams list =
    [ Parameters.GeneralOptions
        [ Parameters.VerbosityOfOutput 2
          Parameters.PostProcessing_TargetDecoyCompetition ]
      Parameters.FileInputOptions
        [ Parameters.PINTAB (System.IO.FileInfo "run1.pin") ]
      Parameters.FileOutputOptions
        [ Parameters.POUTTAB_PSMs (System.IO.FileInfo "run1_psms.tsv")
          Parameters.POUTTAB_DecoyPSMs (System.IO.FileInfo "run1_decoy_psms.tsv") ]
      Parameters.SVMTrainingOptions
        [ Parameters.FDR_CrossValidation 0.01
          Parameters.MaxIterations 10 ] ]

let groupSummary (p: Parameters.PercolatorParams) =
    match p with
    | Parameters.GeneralOptions s                     -> "general options", Seq.length s
    | Parameters.FileInputOptions s                   -> "file input", Seq.length s
    | Parameters.FileOutputOptions s                  -> "file output", Seq.length s
    | Parameters.SVMFeatureOptions s                  -> "SVM features", Seq.length s
    | Parameters.SVMTrainingOptions s                 -> "SVM training", Seq.length s
    | Parameters.ProteinInferenceOptions_FIDO s       -> "Fido inference", Seq.length s
    | Parameters.ProteinInferenceOptions_Percolator s -> "protein inference", Seq.length s

for p in percolatorParams do
    let name, count = groupSummary p
    printfn "%-16s %i option(s)" name count

(**
```text
general options  2 option(s)
file input       1 option(s)
file output      2 option(s)
SVM training     2 option(s)
```

The set reads a pin-tab input file, applies target-decoy competition,
evaluates cross validation at 1 percent FDR with at most ten SVM iterations,
and writes separate tab-delimited result files for target and decoy PSMs.
Running the set goes through `PercolatorWrapper(os).Percolate`, which
renders the options to a command line for the chosen operating system and
launches the percolator process. Executing it requires an external
percolator installation, so this page stops at the constructed parameters.

## Where the identifications go next

Whichever q-value function the pipeline uses, the PSMs that pass the chosen
threshold form the statistically controlled identification list. Those
identifications carry their peptide evidence into
[protein inference]({{root}}04_02_protein_inference.html), and the ion areas
determined by [quantification]({{root}}03_01_quantification.html) are
aggregated over the same accepted set.
*)
