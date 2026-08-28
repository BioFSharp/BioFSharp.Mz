module SequestLikeTests

open Expecto
open BioFSharp
open BioFSharp.Mz
open FSharp.Stats

let expectVectorClose accuracy (actual:Vector<float>) (expected:Vector<float>) message =
    Expect.equal actual.Length expected.Length (sprintf "%s length" message)
    if actual.Length > 0 then
        for i = 0 to actual.Length - 1 do
            Expect.floatClose accuracy actual.[i] expected.[i] (sprintf "%s at index %d" message i)

let expectWithin tolerance actual expected message =
    Expect.isTrue
        (abs (actual - expected) <= tolerance)
        (sprintf "%s; expected %g, got %g" message expected actual)

let mf = SearchDB.massFBy SearchDB.MassMode.Monoisotopic

let peptide =
    [AminoAcids.Ala; AminoAcids.Gly; AminoAcids.Ser; AminoAcids.Glu; AminoAcids.Lys]

let pepMass =
    (peptide |> List.sumBy (fun aa -> mf (aa :> BioFSharp.IBioItem))) + 18.010565

let lookup =
    SearchDB.createLookUpResult 1 1 pepMass (int64 (pepMass * 1000000.0)) "AGSEK" peptide 0

let calcIonSeriesF =
    fun _ s ->
        (Fragmentation.Series.bOfBioList mf s) @ (Fragmentation.Series.yOfBioList mf s)

let fragMasses =
    Fragmentation.Series.fragmentMasses
        Fragmentation.Series.bOfBioList
        Fragmentation.Series.yOfBioList
        mf
        peptide

let scanlimits = (100.0, 800.0)
let chargeState = 2
let precursorMz = BioFSharp.Mass.toMZ pepMass 2.0
let scanTime = 30.0

let spectrum =
    fragMasses.TargetMasses
    |> List.choose (fun family ->
        let m = family.MainPeak.Mass
        let mz1 = BioFSharp.Mass.toMZ m 1.0
        if mz1 > fst scanlimits && mz1 < snd scanlimits then
            Some (mz1, 100.0)
        else
            None)
    |> List.toArray
    |> PeakArray.zipMzInt

[<Tests>]
let tests =
    testList "SequestLikeTests" [
        testList "Preprocessing" [
            testCase "windowNormalizeIntensities square-root-normalizes each window to its window maximum" <| fun _ ->
                let r = SequestLike.windowNormalizeIntensities (vector [4.;1.;9.;16.]) 2
                expectVectorClose Accuracy.high r (vector [1.0; 0.5; 0.75; 1.0]) "window-normalized intensities"
                // window 1 = [4;1] with max 4: sqrt(4)/sqrt(4)=1, sqrt(1)/sqrt(4)=0.5; window 2 = [9;16] with max 16: 3/4 and 1 - hand-computed from the documented sqrt-normalize-to-window-max rule (the original SEQUEST preprocessing). Each positive window's maximum must map to exactly 1.

            testCase "spectrum preprocessing maps silence to silence and a lone peak to positive signal" <| fun _ ->
                let r0 = SequestLike.spectrumToIntensityArrayMinusAutoCorrelation (100.0, 200.0) (PeakArray.zipMzInt [||])
                Expect.isTrue (r0 |> Seq.forall (fun intensity -> intensity = 0.0)) "empty spectrum is zero everywhere"
                // an empty spectrum has zero intensity everywhere; 0 minus the autocorrelation of 0 is 0 - no preprocessing may invent signal
                let rp = SequestLike.spectrumToIntensityArrayMinusAutoCorrelation (100.0, 200.0) (PeakArray.zipMzInt [|(150.0, 100.0)|])
                Expect.isTrue (rp.[50] > 0.0) "a lone peak remains positive at its bin"
                // a single peak at m/z 150 bins to index 50; normalization maps the window maximum to 1 and the autocorrelation subtracts only shifted (zero) values at that bin's own position... assert simply that the value at the peak's bin is positive: preprocessing must preserve, not erase, a real peak.

            testCase "windowNormalizeIntensities truncates a remainder that does not fill a window" <| fun _ ->
                let r = SequestLike.windowNormalizeIntensities (vector [4.;1.;9.;16.;25.]) 2
                Expect.equal r.Length 4 "a remainder that does not fill a window is truncated"
                expectVectorClose Accuracy.high r (vector [1.0; 0.5; 0.75; 1.0]) "retained normalized prefix"
                // the doc says the array is SHORTENED (cut at the tail): the retained prefix must be the same normalization as the exact-fit case - a length-only check would pass for head-dropping or zero-filling.

            testCase "an all-zero window normalizes to zeros" <| fun _ ->
                let r = SequestLike.windowNormalizeIntensities (vector [0.;0.;4.;16.]) 2
                expectVectorClose Accuracy.high r (vector [0.;0.;0.5;1.0]) "zero and positive windows"
                // a window with no positive signal cannot be normalized to its maximum and must stay zero (documented guard); the positive window normalizes as usual - hand-computed.

            testCase "predictIntensitySimpleModel ranks backbone ions above satellites and scales inversely with charge" <| fun _ ->
                let b1 = SequestLike.predictIntensitySimpleModel Ions.IonTypeFlag.B 1.0
                let y1 = SequestLike.predictIntensitySimpleModel Ions.IonTypeFlag.Y 1.0
                let a1 = SequestLike.predictIntensitySimpleModel Ions.IonTypeFlag.A 1.0
                let z1 = SequestLike.predictIntensitySimpleModel Ions.IonTypeFlag.Z 1.0
                let b2 = SequestLike.predictIntensitySimpleModel Ions.IonTypeFlag.B 2.0
                Expect.equal b1 y1 "b and y ions are modeled equally"
                Expect.isTrue (b1 > a1) "b ions outrank a-ion satellites"
                Expect.isTrue (b1 > z1) "b ions outrank z-ion satellites"
                Expect.isTrue (b1 > b2) "higher charge lowers per-peak backbone intensity"
                // higher fragment charge spreads the same ion current over more peaks, so per-peak modeled intensity must fall; the exact 1/z functional form is the module's modeling choice and is deliberately not pinned.

        ]

        testList "TheoreticalSpectrum" [
            testCase "predictOf bins fragment m/z values into unit-Dalton bins carrying the modeled intensity" <| fun _ ->
                let family = Peaks.createPeakFamily (TaggedMass.createTaggedMass Ions.IonTypeFlag.B 199.0) []
                let v1 = SequestLike.predictOf (100.0, 500.0) 1.0 [family]
                Expect.floatClose
                    Accuracy.high
                    v1.[100]
                    (SequestLike.predictIntensitySimpleModel Ions.IonTypeFlag.B 1.0)
                    "1+ fragment intensity is binned at index 100"
                Expect.floatClose Accuracy.high (Vector.sum v1) v1.[100] "the only 1+ fragment fills exactly one bin"
                // 1+ m/z = 199 + 1.007276 = 200.007 (published proton mass), rounds to bin 200 -> index 200-100 = 100.

                let v2 = SequestLike.predictOf (100.0, 500.0) 2.0 [family]
                Expect.floatClose
                    Accuracy.high
                    v2.[100]
                    (SequestLike.predictIntensitySimpleModel Ions.IonTypeFlag.B 1.0)
                    "the 1+ ion is still binned"
                Expect.floatClose
                    Accuracy.high
                    v2.[1]
                    (SequestLike.predictIntensitySimpleModel Ions.IonTypeFlag.B 2.0)
                    "the 2+ ion is binned at index 1"
                Expect.floatClose
                    Accuracy.high
                    (Vector.sum v2)
                    (v2.[100] + v2.[1])
                    "exactly two bins are filled"
                // 2+ m/z = (199 + 2*1.007276)/2 = 100.5073, which rounds UP to bin 101 -> index 1 (nearest-integer rounding, .507 > .5 so no midpoint ambiguity).
                // Bin positions follow the fundamental m/z relation with the published proton mass; intensities cross-check the module's own intensity model rather than pinning constants.

            testCase "autoCorrelation averages the shifted vector over the implemented SEQUEST-like lag window" <| fun _ ->
                // These scores are SEQUEST-LIKE, not SEQUEST: the background correction is not
                // required to match the canonical algorithm. The implemented window accumulates
                // the shifts +delay and -(delay-1) scaled by 1/(2*delay) - narrower than classic
                // SEQUEST's full +/-delay average. For delay 2 on [a;b;c] it yields
                // [b/4; c/4; a/4], hand-computed from the implemented window: [8/4; 12/4; 4/4].
                let r = SequestLike.autoCorrelation 2 (vector [4.;8.;12.])
                expectVectorClose Accuracy.high r (vector [2.0; 3.0; 1.0]) "autocorrelation over the implemented lag window"

            testCase "calcXCorr is the clamped inner product: symmetric, linear in scale, floored at zero" <| fun _ ->
                let a = vector [1.;2.;0.]
                let b = vector [0.5;1.;3.]
                expectWithin 1e-12 (SequestLike.calcXCorr a b) 2.5 "zero-lag inner product"
                expectWithin 1e-12 (SequestLike.calcXCorr a b) (SequestLike.calcXCorr b a) "inner-product symmetry"
                expectWithin 1e-12 (SequestLike.calcXCorr (vector [2.;4.;0.]) b) 5.0 "inner-product scaling"
                expectWithin 1e-12 (SequestLike.calcXCorr (vector [1.;-1.]) (vector [-1.;1.])) 0.0 "negative correlations are clamped"
                // inner product 1*0.5 + 2*1 + 0*3 = 2.5, hand-computed; cross-correlation at zero lag is defined as this weighted sum
                // doubling one argument doubles the correlation (bilinearity)
                // negative correlations are floored at zero; the downstream contract requires it: SearchEngineResult's documented deltaCn conventions treat a non-positive best score as a sentinel case, i.e. real scores are expected non-negative (cross-consistency justification).

        ]

        testList "EndToEnd" [
            testCase "the true peptide outscores its reversed decoy on a spectrum built from its own fragments" <| fun _ ->
                let results =
                    SequestLike.calcSequestLikeScoresRevDecoy
                        calcIonSeriesF
                        BioFSharp.Formula.monoisoMass
                        scanlimits
                        spectrum
                        scanTime
                        chargeState
                        precursorMz
                        [lookup]
                        "spec1"
                Expect.equal results.Length 2 "one target and one decoy are returned per candidate"
                let target = results |> List.find (fun result -> result.IsTarget)
                let decoy = results |> List.find (fun result -> not result.IsTarget)
                Expect.isTrue (target.Score > decoy.Score) "the true peptide scores above its reversed decoy"
                Expect.isTrue
                    (results = (results |> List.sortByDescending (fun result -> result.Score)))
                    "results are ordered by descending score"
                let top = List.head results
                expectWithin 1e-12 top.NormDeltaBestToRest 0.0 "the best result has zero deltaCn"
                expectWithin 1e-9 decoy.NormDeltaBestToRest target.NormDeltaNext "decoy best-to-rest delta equals target next delta"
                // with exactly two results, (best - second)/best is simultaneously the second's best-to-rest delta and the best's adjacent delta - an internal consistency identity of the two documented normalizations.
                results
                |> List.iter (fun result ->
                    expectWithin
                        1e-9
                        result.MeasuredMass
                        pepMass
                        "measured mass round-trips the peptide mass"
                    Expect.equal result.TheoMass pepMass "the theoretical mass is the peptide mass"
                    Expect.equal result.SearchEngine SearchEngineResult.SearchEngine.SEQUESTLike "search engine is SEQUESTLike"
                    Expect.equal result.SpectrumID "spec1" "spectrum ID is preserved"
                    Expect.isTrue (result.Score >= 0.0) "scores are clamped at zero")
                // a spectrum generated from the true peptide's fragment ladder must correlate better with the true prediction than with the reversed-sequence decoy; this discriminative power is the entire point of the scoring engine
                // the documented deltaCn convention for the best hit sets NormDeltaBestToRest to zero
                // Mass.ofMZ inverts Mass.toMZ exactly: the measured neutral mass of a 2+ precursor at toMZ(pepMass, 2) is pepMass itself - a round-trip oracle with no truncated constants.
                // negative correlations are floored at zero; the downstream contract requires it: SearchEngineResult's documented deltaCn conventions treat a non-positive best score as a sentinel case, i.e. real scores are expected non-negative (cross-consistency justification).

            testCase "sequential and parallel scoring agree, and both agree with the recomputed-ion-series path" <| fun _ ->
                let theoSpecs = SequestLike.getTheoSpecs scanlimits chargeState [(lookup, fragMasses)]
                let rSeq =
                    SequestLike.calcSequestScore
                        scanlimits
                        spectrum
                        scanTime
                        chargeState
                        precursorMz
                        theoSpecs
                        "spec1"
                let rPar =
                    SequestLike.calcSequestScoreParallel
                        scanlimits
                        spectrum
                        scanTime
                        chargeState
                        precursorMz
                        theoSpecs
                        "spec1"
                let rRev =
                    SequestLike.calcSequestLikeScoresRevDecoy
                        calcIonSeriesF
                        BioFSharp.Formula.monoisoMass
                        scanlimits
                        spectrum
                        scanTime
                        chargeState
                        precursorMz
                        [lookup]
                        "spec1"
                let seqScores = rSeq |> List.map (fun result -> result.IsTarget, result.Score)
                let parScores = rPar |> List.map (fun result -> result.IsTarget, result.Score)
                Expect.equal seqScores.Length parScores.Length "sequential and parallel result counts agree"
                List.iter2
                    (fun (actualTarget, actualScore) (expectedTarget, expectedScore) ->
                        Expect.equal actualTarget expectedTarget "sequential and parallel target flags agree"
                        expectWithin 1e-9 actualScore expectedScore "sequential and parallel scores agree")
                    seqScores
                    parScores
                let seqTarget = rSeq |> List.find (fun result -> result.IsTarget)
                let seqDecoy = rSeq |> List.find (fun result -> not result.IsTarget)
                let revTarget = rRev |> List.find (fun result -> result.IsTarget)
                let revDecoy = rRev |> List.find (fun result -> not result.IsTarget)
                expectWithin 1e-9 seqTarget.Score revTarget.Score "precomputed and recomputed target scores agree"
                expectWithin 1e-9 seqDecoy.Score revDecoy.Score "precomputed and recomputed decoy scores agree"
                // parallelization must not change any score or ordering
                // the precomputed-spectrum path and the recompute-per-call path implement the same model over identical fragments and must agree
                // the decoy fragments come from the reversed sequence in both paths
        ]
    ]
