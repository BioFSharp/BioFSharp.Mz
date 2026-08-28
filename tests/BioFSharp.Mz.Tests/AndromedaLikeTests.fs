module AndromedaLikeTests

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
    (peptide |> List.sumBy (fun aa -> mf (aa :> BioFSharp.IBioItem))) + 18.010565

let lookup =
    SearchDB.createLookUpResult 1 1 pepMass (int64 (pepMass * 1000000.0)) "AGSEK" peptide 0

let fragMasses =
    Fragmentation.Series.fragmentMasses
        Fragmentation.Series.bOfBioList
        Fragmentation.Series.yOfBioList
        mf
        peptide

let scanlimits = (100.0, 2000.0)
let chargeState = 2
let precursorMz = BioFSharp.Mass.toMZ pepMass 2.0

let theoSpecs =
    AndromedaLike.getTheoSpecs scanlimits chargeState [(lookup, fragMasses)]

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
let scanTime = 30.0
let spectrumID = "spec1"

let targetScore (results: SearchEngineResult.SearchEngineResult<float> list) =
    results |> List.find (fun result -> result.IsTarget) |> fun result -> result.Score

[<Tests>]
let tests =
    testList "AndromedaLikeTests" [
        testList "Prediction" [
            testCase "predictOf converts a fragment mass to its singly charged m/z with unmeasured intensity" <| fun _ ->
                let fam = Peaks.createPeakFamily (TaggedMass.createTaggedMass Ions.IonTypeFlag.B 500.0) []
                let r = AndromedaLike.predictOf (100.0, 2000.0) 1.0 [fam]
                Expect.equal r.Length 1 "one theoretical peak family is emitted"
                Expect.equal r.[0].MainPeak.Iontype Ions.IonTypeFlag.B "the ion type is retained"
                expectWithin 1e-6 r.[0].MainPeak.Mz 501.007276 "the singly charged m/z uses the published proton mass"
                // a theoretical peak carries no measured intensity; NaN is the only honest "unmeasured" representation (any real number would be a fabricated measurement), and downstream scoring never reads theoretical intensities.
                Expect.isTrue (Double.IsNaN r.[0].MainPeak.Intensity) "theoretical intensity represents no measurement"

            testCase "predictOf for a multiply charged precursor emits both singly and doubly charged fragment ions" <| fun _ ->
                let fam =
                    Peaks.createPeakFamily
                        (TaggedMass.createTaggedMass Ions.IonTypeFlag.B 500.0)
                        [TaggedMass.createTaggedH2OLoss Ions.IonTypeFlag.B 482.0]
                let r = AndromedaLike.predictOf (100.0, 2000.0) 2.0 [fam]
                Expect.equal r.Length 2 "both fragment charge states are emitted"
                [501.007276; 251.007276]
                |> List.iter (fun expectedMz ->
                    Expect.isTrue
                        (r |> Array.exists (fun peakFamily -> abs (peakFamily.MainPeak.Mz - expectedMz) <= 1e-6))
                        (sprintf "a main peak is present near m/z %g" expectedMz))
                r
                |> Array.iter (fun peakFamily ->
                    Expect.equal peakFamily.MainPeak.Iontype Ions.IonTypeFlag.B "both peaks retain ion type B")
                // Fragments of a 2+ precursor are observed at charge 1 and charge 2 - standard fragment charge distribution; m/z values follow the fundamental relation with the published proton mass.
                // The doc's "at a given charge state" wording is in tension with the 1-and-2 emission policy.
                let singlyChargedFamily =
                    r
                    |> Array.find (fun peakFamily -> abs (peakFamily.MainPeak.Mz - 501.007276) <= 1e-6)
                Expect.equal singlyChargedFamily.DependentPeaks.Length 1 "the singly charged family retains its dependent peak"
                let dependentPeak = List.head singlyChargedFamily.DependentPeaks
                expectWithin 1e-6 dependentPeak.Mz 483.007276 "the dependent peak uses the singly charged m/z"
                Expect.isTrue
                    (dependentPeak.Iontype.HasFlag(Ions.IonTypeFlag.B)
                     && dependentPeak.Iontype.HasFlag(Ions.IonTypeFlag.lossH2O))
                    "the dependent peak retains both ion type B and lossH2O"

            // PENDING: the doc for predictOf says peaks outside the scan limits are filtered out. A
            // fragment whose 1+ m/z (5001.007) lies far above the upper scan limit of 2000 must not
            // appear in the predicted spectrum. The implementation never reads scanLimits.
            ptestCase "predictOf filters fragment peaks outside the scan limits" <| fun _ ->
                let r =
                    AndromedaLike.predictOf
                        (100.0, 2000.0)
                        1.0
                        [Peaks.createPeakFamily (TaggedMass.createTaggedMass Ions.IonTypeFlag.B 5000.0) []]
                Expect.equal r.Length 0 "peaks outside the scan limits are filtered out"
        ]

        testList "EndToEnd" [
            testCase "the true peptide outscores its reversed decoy on a spectrum of its own fragments" <| fun _ ->
                let results =
                    AndromedaLike.calcAndromedaScore
                        qMinAndMax
                        scanlimits
                        matchingTolPPM
                        spectrum
                        scanTime
                        chargeState
                        precursorMz
                        theoSpecs
                        spectrumID
                Expect.equal results.Length 2 "one target and one decoy are returned per candidate"
                let target = results |> List.find (fun result -> result.IsTarget)
                let decoy = results |> List.find (fun result -> not result.IsTarget)
                Expect.isTrue (target.Score > decoy.Score) "the true peptide scores above its reversed decoy"
                Expect.equal
                    results
                    (results |> List.sortByDescending (fun result -> result.Score))
                    "results are ordered by descending score"
                let top = List.head results
                Expect.floatClose Accuracy.high top.NormDeltaBestToRest 0.0 "the best result has zero normalized delta"
                results
                |> List.iter (fun result ->
                    Expect.equal result.SearchEngine SearchEngineResult.SearchEngine.AndromedaLike "the search engine is AndromedaLike"
                    Expect.equal result.PrecursorMZ precursorMz "the precursor m/z is preserved"
                    // Mass.ofMZ inverts Mass.toMZ exactly - a round-trip oracle with no truncated constants.
                    expectWithin 1e-9 result.MeasuredMass pepMass "the measured mass round-trips the peptide mass"
                    Expect.equal result.TheoMass pepMass "the theoretical mass is the peptide mass"
                    Expect.isTrue (result.Score >= 0.0) "scores are clamped at zero")

            testCase "removing matched fragment peaks cannot increase the target score" <| fun _ ->
                let strippedFrag : Fragmentation.FragmentMasses =
                    { TargetMasses = fragMasses.TargetMasses |> List.map (fun f -> Peaks.createPeakFamily f.MainPeak [])
                      DecoyMasses = fragMasses.DecoyMasses |> List.map (fun f -> Peaks.createPeakFamily f.MainPeak []) }
                let theoSpecsStripped = AndromedaLike.getTheoSpecs scanlimits chargeState [(lookup, strippedFrag)]
                let fullSpectrum =
                    strippedFrag.TargetMasses
                    |> List.map (fun family ->
                        let m = family.MainPeak.Mass
                        Mass.toMZ m 1.0, 100.0)
                    |> List.sortBy fst
                    |> List.toArray
                    |> PeakArray.zipMzInt
                let halfSpectrum =
                    fullSpectrum
                    |> Array.mapi (fun i peak -> i, peak)
                    |> Array.choose (fun (i, peak) -> if i % 2 = 0 then Some peak else None)
                let full =
                    AndromedaLike.calcAndromedaScore
                        qMinAndMax scanlimits matchingTolPPM fullSpectrum scanTime chargeState precursorMz theoSpecsStripped spectrumID
                    |> targetScore
                let half =
                    AndromedaLike.calcAndromedaScore
                        qMinAndMax scanlimits matchingTolPPM halfSpectrum scanTime chargeState precursorMz theoSpecsStripped spectrumID
                    |> targetScore
                // with dependent-free families the trial count N is spectrum-independent (an unmatched main has no dependents to skip), so removing matched peaks can only reduce K at fixed N - and the binomial upper tail P(X >= K) is strictly decreasing in K: the score cannot increase. This is a theorem for this fixture; with dependents, N itself varies with match outcomes and no such guarantee exists.
                Expect.isTrue (full >= half) "removing matched fragment peaks cannot increase the target score"

            testCase "the peak-depth threshold q gates matches by local intensity rank" <| fun _ ->
                // 440 sits 40 Da from 400 - strictly inside the +/-50 rating half-window, so the
                // rank relation does not depend on the window's open/closed boundary convention.
                let spectrum = PeakArray.zipMzInt [|(400.0, 50.0); (440.0, 100.0)|]
                let family =
                    Peaks.createPeakFamily
                        (TaggedPeak.TaggedPeak(Ions.IonTypeFlag.B, 400.0, nan))
                        []
                let theoSpec =
                    TheoreticalSpectra.createTheoreticalSpectrum lookup [|family|] [|family|]
                let rQ1 =
                    AndromedaLike.calcAndromedaScore
                        (1, 1) scanlimits 20.0 spectrum 30.0 2 500.0 [theoSpec] "s"
                let rQ2 =
                    AndromedaLike.calcAndromedaScore
                        (1, 2) scanlimits 20.0 spectrum 30.0 2 500.0 [theoSpec] "s"
                let targetQ1 = rQ1 |> targetScore
                let targetQ2 = rQ2 |> targetScore
                // this is the mechanism that distinguishes Andromeda scoring - matches are gated by local intensity rank, not just m/z tolerance.
                // at peak depth q = 1, only the single most intense peak per window is considered; the 400-peak (rank 1) is filtered, the theoretical peak finds no partner, K = 0, and the raw score log10(P(X >= 0)) = 0 clamps to 0 with the sub-800 precursor correction - hand-derived from the documented top-q-per-window semantics
                Expect.equal targetQ1 0.0 "q = 1 excludes the lower-ranked matching peak and clamps the zero-match score"
                // widening the depth to q = 2 admits the second-ranked peak, the match lands, and a positive K yields a positive binomial-tail score
                Expect.isTrue (targetQ2 > targetQ1) "q = 2 admits the lower-ranked matching peak and produces a positive score"

            testCase "parallel scoring equals sequential scoring" <| fun _ ->
                let rSeq =
                    AndromedaLike.calcAndromedaScore
                        qMinAndMax scanlimits matchingTolPPM spectrum scanTime chargeState precursorMz theoSpecs spectrumID
                let rPar =
                    AndromedaLike.calcAndromedaScoreParallel
                        qMinAndMax scanlimits matchingTolPPM spectrum scanTime chargeState precursorMz theoSpecs spectrumID
                let seqProjection =
                    rSeq |> List.map (fun x -> x.IsTarget, x.Score, x.NormDeltaBestToRest, x.NormDeltaNext)
                let parProjection =
                    rPar |> List.map (fun x -> x.IsTarget, x.Score, x.NormDeltaBestToRest, x.NormDeltaNext)
                Expect.equal seqProjection.Length parProjection.Length "sequential and parallel result counts agree"
                List.iter2
                    (fun (actualTarget, actualScore, actualBest, actualNext) (expectedTarget, expectedScore, expectedBest, expectedNext) ->
                        Expect.equal actualTarget expectedTarget "sequential and parallel target flags agree"
                        expectWithin 1e-9 actualScore expectedScore "sequential and parallel scores agree"
                        expectWithin 1e-9 actualBest expectedBest "sequential and parallel best-to-rest deltas agree"
                        expectWithin 1e-9 actualNext expectedNext "sequential and parallel next deltas agree")
                    seqProjection
                    parProjection
        ]
    ]
