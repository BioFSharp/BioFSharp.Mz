module XScoringTests

open System
open Expecto
open BioFSharp
open BioFSharp.Mz

let mf = SearchDB.massFBy SearchDB.MassMode.Monoisotopic

let expectWithin tolerance actual expected message =
    Expect.isTrue
        (abs (actual - expected) <= tolerance)
        (sprintf "%s; expected %g, got %g" message expected actual)

let peptide =
    [AminoAcids.Ala; AminoAcids.Gly; AminoAcids.Ser; AminoAcids.Glu; AminoAcids.Lys]

let pepMass =
    (peptide |> List.sumBy (fun aa -> mf (aa :> BioFSharp.IBioItem))) + 18.010565 // published monoisotopic water mass

let lookup =
    SearchDB.createLookUpResult 1 1 pepMass (int64 (pepMass * 1e6)) "AGSEK" peptide 0

let fragMasses =
    Fragmentation.Series.fragmentMasses
        Fragmentation.Series.bOfBioList
        Fragmentation.Series.yOfBioList
        mf
        peptide

let scanlimits = (100.0, 2000.0)
let chargeState = 2
let precursorMz = Mass.toMZ pepMass 2.0

let theoSpecsX = XScoring.getTheoSpecs scanlimits chargeState [(lookup, fragMasses)]

let spectrum =
    fragMasses.TargetMasses
    |> List.map (fun family ->
        let m = family.MainPeak.Mass
        Mass.toMZ m 1.0, 100.0)
    |> List.sortBy fst
    |> List.toArray
    |> PeakArray.zipMzInt

let qMinAndMax = (1, 10)
let matchingTolPPM = 20.0

let targetScore (results: SearchEngineResult.SearchEngineResult<float> list) =
    results |> List.find (fun result -> result.IsTarget) |> fun result -> result.Score

let taggedFamilyAtMz mz =
    Peaks.createPeakFamily
        (TaggedPeak.TaggedPeak(Ions.IonTypeFlag.B, mz, nan))
        []

let lookupB =
    SearchDB.createLookUpResult 2 2 (pepMass + 1.0) (int64 ((pepMass + 1.0) * 1e6)) "B" peptide 0

[<Tests>]
let tests =
    testList "XScoringTests" [
        testList "Counting" [
            testCase "findMaxIterator returns the first q index at which a peak's local rank qualifies" <| fun _ ->
                Expect.equal
                    (XScoring.findMaxIterator 1 10 (XScoring.createRatedPeak (500.0, 100.0) 0))
                    0
                    "a peak with no more-intense neighbors qualifies at every q >= 1"
                Expect.equal
                    (XScoring.findMaxIterator 2 10 (XScoring.createRatedPeak (500.0, 100.0) 3))
                    2
                    "a peak with rank 3 first qualifies at q = 4 when qMin = 2"
                // the countArrayL parameter is dead in the implementation (no clamping); out-of-range results are harmless because the crediting loop never reaches them.
                // q is the "top q most abundant peaks per window" threshold; a peak ranked r qualifies exactly when q > r - the Andromeda peak-depth definition, hand-computed.

            testCase "a matched B peak raises N everywhere but K and matched intensity only at qualifying depths" <| fun _ ->
                let countArray =
                    Array.init 3 (fun i -> XScoring.createCountedMatches (i + 1) 0 0 0 0 0 0 0.0)
                XScoring.raiseCountsAndHits
                    1
                    countArray
                    (TaggedPeak.TaggedPeak(Ions.IonTypeFlag.B, 500.0, nan))
                    (XScoring.createRatedPeak (500.0, 80.0) 1)
                countArray
                |> Array.iteri (fun i entry ->
                    Expect.equal entry.N_BIons 1 "the theoretical B peak is offered at every depth"
                    if i >= 1 then
                        Expect.equal entry.K_BIons 1 "the matched B peak survives at qualifying depths"
                        Expect.equal entry.MatchedSum 80.0 "only surviving matches contribute intensity"
                    else
                        Expect.equal entry.K_BIons 0 "the rank-1 B peak does not qualify at q = 1"
                        Expect.equal entry.MatchedSum 0.0 "the filtered match contributes no intensity"
                    Expect.equal entry.N_YIons 0 "a B peak does not raise Y offers"
                    Expect.equal entry.K_YIons 0 "a B peak does not raise Y hits"
                    Expect.equal entry.N_Neutral 0 "a B peak does not raise neutral offers"
                    Expect.equal entry.K_Neutral 0 "a B peak does not raise neutral hits")
                // N counts offered theoretical peaks, K counts matches surviving the local-abundance threshold, and only surviving matches contribute intensity - the Andromeda/X!Tandem counting semantics, hand-computed for rank 1.

            testCase "composite loss tags are classified as their backbone ion, leaving the Neutral machinery untouched" <| fun _ ->
                let countArray =
                    Array.init 3 (fun i -> XScoring.createCountedMatches (i + 1) 0 0 0 0 0 0 0.0)
                XScoring.raiseCountsAndHits
                    1
                    countArray
                    (TaggedPeak.TaggedPeak(Ions.IonTypeFlag.B ||| Ions.IonTypeFlag.lossH2O, 500.0, nan))
                    (XScoring.createRatedPeak (500.0, 60.0) 0)
                countArray
                |> Array.iter (fun entry ->
                    Expect.equal entry.N_BIons 1 "a composite B/loss peak raises B offers"
                    Expect.equal entry.K_BIons 1 "a composite B/loss peak raises B hits"
                    Expect.equal entry.MatchedSum 60.0 "a composite B/loss peak contributes matched intensity"
                    Expect.equal entry.N_Neutral 0 "a composite B/loss peak does not raise neutral offers"
                    Expect.equal entry.K_Neutral 0 "a composite B/loss peak does not raise neutral hits"
                    Expect.equal entry.N_YIons 0 "a composite B/loss peak does not raise Y offers")
                // every loss dependent the fragment model produces carries a composite backbone|loss flag (never the separate Neutral flag), and classification is first-match-wins - so loss peaks are counted as their backbone ion and the Neutral counters (and the hyperscore's neutral factorial, and Andromeda's with/without-neutral split) are unreachable from production input. Characterization of live behavior.

            testCase "non-backbone ion tags contribute nothing to the counts" <| fun _ ->
                let countArray =
                    Array.init 3 (fun i -> XScoring.createCountedMatches (i + 1) 0 0 0 0 0 0 0.0)
                XScoring.raiseCountsAndHits
                    1
                    countArray
                    (TaggedPeak.TaggedPeak(Ions.IonTypeFlag.X, 500.0, nan))
                    (XScoring.createRatedPeak (500.0, 60.0) 0)
                XScoring.raiseOnlyCounts countArray (TaggedPeak.TaggedPeak(Ions.IonTypeFlag.C, 600.0, nan))
                countArray
                |> Array.iter (fun entry ->
                    Expect.equal entry.N_BIons 0 "an X/C peak does not raise B offers"
                    Expect.equal entry.N_YIons 0 "an X/C peak does not raise Y offers"
                    Expect.equal entry.N_Neutral 0 "an X/C peak does not raise neutral offers"
                    Expect.equal entry.K_BIons 0 "an X/C peak does not raise B hits"
                    Expect.equal entry.K_YIons 0 "an X/C peak does not raise Y hits"
                    Expect.equal entry.K_Neutral 0 "an X/C peak does not raise neutral hits"
                    Expect.equal entry.MatchedSum 0.0 "an X/C peak does not contribute matched intensity")
                // the scoring model counts only b and y ions (plus the dead Neutral category); matched or offered a/c/x/z peaks fall through the default branch untouched - a B/Y-only model, defensible for CID.

            testCase "an unmatched theoretical Y peak raises only the Y offer count" <| fun _ ->
                let countArray =
                    Array.init 3 (fun i -> XScoring.createCountedMatches (i + 1) 0 0 0 0 0 0 0.0)
                XScoring.raiseOnlyCounts countArray (TaggedPeak.TaggedPeak(Ions.IonTypeFlag.Y, 600.0, nan))
                countArray
                |> Array.iter (fun entry ->
                    Expect.equal entry.N_YIons 1 "the unmatched theoretical Y peak raises every Y offer count"
                    Expect.equal entry.K_YIons 0 "an unmatched theoretical Y peak raises no Y hits"
                    Expect.equal entry.MatchedSum 0.0 "an unmatched theoretical peak contributes no intensity"
                    Expect.equal entry.N_BIons 0 "an unmatched Y peak does not raise B offers"
                    Expect.equal entry.N_Neutral 0 "an unmatched Y peak does not raise neutral offers"
                    Expect.equal entry.K_Neutral 0 "an unmatched Y peak does not raise neutral hits")
                // an unmatched theoretical peak enlarges the trial count but can neither be a hit nor contribute intensity.
        ]

        testList "HyperScore" [
            testCase "calcHyperScore is the log X!Tandem hyperscore" <| fun _ ->
                let s = XScoring.calcHyperScore (XScoring.createCountedMatches 4 5 5 0 2 3 0 1000.0)
                expectWithin
                    1e-9
                    s
                    (log 1000.0 + log 2.0 + log 6.0)
                    "ln(matchedSum) + ln(Kb!) + ln(Ky!) + ln(Kneutral!)"
                Expect.equal
                    (XScoring.calcHyperScore (XScoring.createCountedMatches 4 0 0 0 0 0 0 0.0))
                    0.0
                    "a no-match hyperscore is clamped to zero"
                // no matches means no explained intensity: ln(0) = -infinity is clamped to the floor 0 - the production no-match score.
                // the published X!Tandem hyperscore is (matched intensity sum) * b! * y!; in log space that is ln(sum) + ln(Kb!) + ln(Ky!) - hand-computed with 2! = 2, 3! = 6, 0! = 1.

            testCase "the hyperscore grows with matched intensity and with matched ion counts" <| fun _ ->
                let higherIntensity =
                    XScoring.calcHyperScore (XScoring.createCountedMatches 4 5 5 0 2 3 0 2000.0)
                let lowerIntensity =
                    XScoring.calcHyperScore (XScoring.createCountedMatches 4 5 5 0 2 3 0 1000.0)
                Expect.isTrue
                    (higherIntensity > lowerIntensity)
                    "more explained intensity produces a better score"
                let moreBMatches =
                    XScoring.calcHyperScore (XScoring.createCountedMatches 4 5 5 0 3 3 0 1000.0)
                let fewerBMatches =
                    XScoring.calcHyperScore (XScoring.createCountedMatches 4 5 5 0 2 3 0 1000.0)
                Expect.isTrue
                    (moreBMatches > fewerBMatches)
                    "an additional matched B ion increases b! and the score"
                Expect.floatClose
                    Accuracy.high
                    (XScoring.calcHyperScore (XScoring.createCountedMatches 4 5 5 0 1 0 0 1000.0))
                    (XScoring.calcHyperScore (XScoring.createCountedMatches 4 5 5 0 0 0 0 1000.0))
                    "0! = 1! so counts 0 and 1 contribute equally"
                // monotonicity follows from log/factorial monotonicity - mathematical properties, no constants pinned.
        ]

        testList "EndToEnd" [
            testCase "combined scoring ranks the true peptide above its decoy in both engines" <| fun _ ->
                let andro, xtandem =
                    XScoring.calcAndromedaAndXTandemScore
                        qMinAndMax
                        scanlimits
                        matchingTolPPM
                        spectrum
                        30.0
                        chargeState
                        precursorMz
                        theoSpecsX
                        "spec1"
                Expect.equal andro.Length 2 "target and decoy are returned for AndromedaLike"
                Expect.equal xtandem.Length 2 "target and decoy are returned for XTandemLike"
                Expect.isTrue
                    (andro |> List.forall (fun result -> result.SearchEngine = SearchEngineResult.SearchEngine.AndromedaLike))
                    "every Andromeda result is tagged AndromedaLike"
                let androTarget = andro |> List.find (fun result -> result.IsTarget)
                let androDecoy = andro |> List.find (fun result -> not result.IsTarget)
                Expect.isTrue (androTarget.Score > androDecoy.Score) "the true peptide scores above the decoy in AndromedaLike"
                Expect.isTrue
                    (xtandem |> List.forall (fun result -> result.SearchEngine = SearchEngineResult.SearchEngine.XTandemLike))
                    "every XTandem result is tagged XTandemLike"
                let xtandemTarget = xtandem |> List.find (fun result -> result.IsTarget)
                let xtandemDecoy = xtandem |> List.find (fun result -> not result.IsTarget)
                Expect.isTrue (xtandemTarget.Score > xtandemDecoy.Score) "the true peptide scores above the decoy in XTandemLike"
                // the decoy shares only the sequence-invariant full-length b_n/y_n families with the spectrum, so its K counts and matched intensity are strictly smaller - strict inequality is guaranteed by construction.
                Expect.isTrue (xtandemTarget.Score > 0.0) "the true peptide has a positive XTandemLike score"
                Expect.floatClose Accuracy.high (List.head andro).NormDeltaBestToRest 0.0 "the top Andromeda result has zero normalized delta"
                Expect.floatClose Accuracy.high (List.head xtandem).NormDeltaBestToRest 0.0 "the top XTandem result has zero normalized delta"
                // a spectrum made of the true peptide's fragments must favor the true peptide under both scoring models; engine tags keep the two result streams distinguishable for downstream FDR.

            testCase "XScoring's Andromeda results agree with the AndromedaLike module" <| fun _ ->
                let theoSpecsA = AndromedaLike.getTheoSpecs scanlimits chargeState [(lookup, fragMasses)]
                let reference =
                    AndromedaLike.calcAndromedaScore
                        qMinAndMax
                        scanlimits
                        matchingTolPPM
                        spectrum
                        30.0
                        chargeState
                        precursorMz
                        theoSpecsA
                        "spec1"
                let andro, _ =
                    XScoring.calcAndromedaAndXTandemScore
                        qMinAndMax
                        scanlimits
                        matchingTolPPM
                        spectrum
                        30.0
                        chargeState
                        precursorMz
                        theoSpecsX
                        "spec1"
                let actual = andro |> List.map (fun x -> x.IsTarget, x.Score)
                let expected = reference |> List.map (fun x -> x.IsTarget, x.Score)
                Expect.equal actual.Length expected.Length "the two Andromeda implementations return the same number of results"
                List.iter2
                    (fun (actualTarget, actualScore) (expectedTarget, expectedScore) ->
                        Expect.equal actualTarget expectedTarget "the target flags agree"
                        expectWithin 1e-9 actualScore expectedScore "the Andromeda scores agree")
                        actual
                    expected
                let unequalSpectrum =
                    fragMasses.TargetMasses
                    |> List.map (fun family -> Mass.toMZ family.MainPeak.Mass 1.0)
                    |> List.sort
                    |> List.toArray
                    |> Array.mapi (fun i mz -> mz, if i % 2 = 0 then 100.0 else 40.0)
                    |> PeakArray.zipMzInt
                let referenceUnequal =
                    AndromedaLike.calcAndromedaScore
                        qMinAndMax
                        scanlimits
                        matchingTolPPM
                        unequalSpectrum
                        30.0
                        chargeState
                        precursorMz
                        theoSpecsA
                        "spec1"
                let androUnequal, _ =
                    XScoring.calcAndromedaAndXTandemScore
                        qMinAndMax
                        scanlimits
                        matchingTolPPM
                        unequalSpectrum
                        30.0
                        chargeState
                        precursorMz
                        theoSpecsX
                        "spec1"
                let actualUnequal = androUnequal |> List.map (fun x -> x.IsTarget, x.Score)
                let expectedUnequal = referenceUnequal |> List.map (fun x -> x.IsTarget, x.Score)
                Expect.equal actualUnequal.Length expectedUnequal.Length "the two Andromeda implementations return the same number of results for unequal intensities"
                List.iter2
                    (fun (actualTarget, actualScore) (expectedTarget, expectedScore) ->
                        Expect.equal actualTarget expectedTarget "the target flags agree for unequal intensities"
                        expectWithin 1e-9 actualScore expectedScore "the Andromeda scores agree for unequal intensities")
                    actualUnequal
                    expectedUnequal
                // unequal intensities produce nonzero local ranks, exercising the duplicated peak-rating and depth machinery through the agreement oracle.
                // XScoring duplicates the Andromeda-like model for B/Y fragment spectra; the two public implementations must agree on identical inputs - cross-module consistency, the preferred oracle class. If they diverge, one of them is wrong.
                // agreement guards against divergence of the duplicated code paths, not against their shared defects (hardcoded corrections, dead scanLimits, charge cap); this oracle cannot certify correctness, only consistency.

            testCase "the ppm tolerance bounds XTandem matching" <| fun _ ->
                let theoSpec =
                    TheoreticalSpectra.createTheoreticalSpectrum
                        lookup
                        [|taggedFamilyAtMz 500.0|]
                        [|taggedFamilyAtMz 1500.0|]
                let scoreFor tolerance =
                    let _, xtandem =
                        XScoring.calcAndromedaAndXTandemScore
                            (1, 1)
                            scanlimits
                            tolerance
                            (PeakArray.zipMzInt [|(500.004, 100.0)|])
                            30.0
                            2
                            500.0
                            [theoSpec]
                            "s"
                    xtandem |> targetScore
                let inside = scoreFor 10.0
                let outside = scoreFor 5.0
                // hyperscore for one matched B ion: ln(matchedSum) + ln(1!) = ln 100 - hand-derived; the 0.004 Da offset sits between the two allowances.
                expectWithin 1e-6 inside (log 100.0) "a 0.004 Da offset is inside the 10 ppm tolerance"
                Expect.equal outside 0.0 "a 0.004 Da offset is outside the 5 ppm tolerance"

            testCase "XTandem matching honors the local rating window and peak depth" <| fun _ ->
                let theoSpec =
                    TheoreticalSpectra.createTheoreticalSpectrum
                        lookup
                        [|taggedFamilyAtMz 400.0|]
                        [|taggedFamilyAtMz 1500.0|]
                let scoreFor spectrum =
                    let _, xtandem =
                        XScoring.calcAndromedaAndXTandemScore
                            (1, 1)
                            scanlimits
                            20.0
                            (PeakArray.zipMzInt spectrum)
                            30.0
                            2
                            500.0
                            [theoSpec]
                            "s"
                    xtandem |> targetScore
                let insideWindowCompetitor = scoreFor [|(400.0, 50.0); (440.0, 100.0)|]
                let outsideWindowCompetitor = scoreFor [|(400.0, 50.0); (460.0, 100.0)|]
                // same 100-Da window and depth semantics as the Andromeda side, observed through the hyperscore: ln(50) for the sole matched ion.
                Expect.equal insideWindowCompetitor 0.0 "a more intense peak 40 Da away excludes the lower-ranked match"
                expectWithin 1e-6 outsideWindowCompetitor (log 50.0) "a more intense peak 60 Da away leaves the match eligible"

            testCase "every candidate is scored under its own identity" <| fun _ ->
                let candidateA =
                    TheoreticalSpectra.createTheoreticalSpectrum
                        lookup
                        [|taggedFamilyAtMz 400.0|]
                        [|taggedFamilyAtMz 1500.0|]
                let candidateB =
                    TheoreticalSpectra.createTheoreticalSpectrum
                        lookupB
                        [|taggedFamilyAtMz 800.0|]
                        [|taggedFamilyAtMz 1500.0|]
                let andro, xtandem =
                    XScoring.calcAndromedaAndXTandemScore
                        (1, 10)
                        scanlimits
                        20.0
                        (PeakArray.zipMzInt [|(400.0, 100.0)|])
                        30.0
                        2
                        500.0
                        [candidateA; candidateB]
                        "s"
                let assertStream name (stream: SearchEngineResult.SearchEngineResult<float> list) =
                    Expect.equal stream.Length 4 (sprintf "%s returns target and decoy for both candidates" name)
                    Expect.equal
                        (stream |> List.countBy (fun result -> result.PepSequenceID) |> List.sort)
                        [(1, 2); (2, 2)]
                        (sprintf "%s contains exactly two records per candidate ID" name)
                    Expect.equal
                        (stream |> List.map (fun result -> result.Score))
                        (stream |> List.map (fun result -> result.Score) |> List.sortDescending)
                        (sprintf "%s is sorted by descending score" name)
                    let head = List.head stream
                    Expect.equal head.PepSequenceID 1 (sprintf "%s ranks candidate A first" name)
                    Expect.isTrue head.IsTarget (sprintf "%s ranks candidate A's target first" name)
                // a scorer that processes only the first candidate, or reuses one candidate's spectra, fails the per-ID accounting; A matches and B does not, so A-target must lead both streams.
                assertStream "Andromeda" andro
                assertStream "XTandem" xtandem
        ]
    ]
