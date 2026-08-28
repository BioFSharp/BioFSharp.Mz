module FragmentationTests

open Expecto
open BioFSharp
open BioFSharp.Mz

let mf = BioFSharp.BioItem.monoisoMass
let ags = [AminoAcids.Ala; AminoAcids.Gly; AminoAcids.Ser]

let mainMasses (fams: PeakFamily<TaggedMass.TaggedMass> list) = fams |> List.map (fun f -> f.MainPeak.Mass)

let containsMass (tol: float) (m: float) (masses: float list) =
    masses |> List.exists (fun x -> abs (x - m) <= tol)

let expectWithin (tol: float) (actual: float) (expected: float) message =
    Expect.isTrue
        (abs (actual - expected) <= tol)
        (sprintf "%s; expected %g, got %g" message expected actual)

let familyWithMainMass (tol: float) (mass: float) (fams: PeakFamily<TaggedMass.TaggedMass> list) =
    fams
    |> List.tryFind (fun fam -> abs (fam.MainPeak.Mass - mass) <= tol)
    |> Option.defaultWith (fun () -> failtestf "No family found near mass %f" mass)

[<Tests>]
let tests =
    testList "FragmentationTests" [
        testList "IonSeriesMasses" [
            testCase "b-series masses are cumulative N-terminal residue sums" <| fun _ ->
                let masses = mainMasses (Fragmentation.Series.bOfBioList mf ags)
                Expect.isTrue
                    (containsMass 0.001 71.03711 masses)
                    (sprintf "b1 is the Ala residue mass; expected a mass within %g of %g; got %A" 0.001 71.03711 masses)
                Expect.isTrue
                    (containsMass 0.001 128.05857 masses)
                    (sprintf "b2 is the Ala + Gly residue mass; expected a mass within %g of %g; got %A" 0.001 128.05857 masses)
                // The b fragment ion of a peptide is the neutral cumulative sum of its N-terminal residue masses - textbook fragment chemistry, values from published residue mass tables.

            testCase "series functions tag every main peak with their own ion type, and combined series are the union of the singles" <| fun _ ->
                let seriesCases =
                    [ (Fragmentation.Series.aOfBioList, Ions.IonTypeFlag.A)
                      (Fragmentation.Series.bOfBioList, Ions.IonTypeFlag.B)
                      (Fragmentation.Series.cOfBioList, Ions.IonTypeFlag.C)
                      (Fragmentation.Series.xOfBioList, Ions.IonTypeFlag.X)
                      (Fragmentation.Series.yOfBioList, Ions.IonTypeFlag.Y)
                      (Fragmentation.Series.zOfBioList, Ions.IonTypeFlag.Z) ]

                seriesCases
                |> List.iter (fun (seriesF, expectedFlag) ->
                    seriesF mf ags
                    |> List.iter (fun family ->
                        Expect.equal
                            family.MainPeak.Iontype
                            expectedFlag
                            (sprintf "every main peak from %A carries its own ion flag" expectedFlag)))

                let mainProjection (families: PeakFamily<TaggedMass.TaggedMass> list) =
                    families
                    |> List.map (fun family ->
                        family.MainPeak.Iontype, System.Math.Round(family.MainPeak.Mass, 6))
                    |> List.sort

                let abcSingles =
                    [ Fragmentation.Series.aOfBioList
                      Fragmentation.Series.bOfBioList
                      Fragmentation.Series.cOfBioList ]
                    |> List.collect (fun seriesF -> mainProjection (seriesF mf ags))
                    |> List.sort
                let xyzSingles =
                    [ Fragmentation.Series.xOfBioList
                      Fragmentation.Series.yOfBioList
                      Fragmentation.Series.zOfBioList ]
                    |> List.collect (fun seriesF -> mainProjection (seriesF mf ags))
                    |> List.sort

                Expect.equal
                    (mainProjection (Fragmentation.Series.abcOfBioList mf ags))
                    abcSingles
                    "the abc wrapper is the union of a, b, and c main peaks"
                Expect.equal
                    (mainProjection (Fragmentation.Series.xyzOfBioList mf ags))
                    xyzSingles
                    "the xyz wrapper is the union of x, y, and z main peaks"
                // downstream scoring dispatches on ion flags (intensity models, B/Y counting) - chemically correct masses with wrong tags corrupt every score; the union property pins that combined wrappers select exactly the requested series.

            testCase "fragment masses are computed with the SUPPLIED mass function" <| fun _ ->
                let mfAvg = SearchDB.massFBy SearchDB.MassMode.Average
                let bAvg = mainMasses (Fragmentation.Series.bOfBioList mfAvg ags)
                Expect.isTrue
                    (containsMass 0.001 (mfAvg (AminoAcids.Ala :> BioFSharp.IBioItem)) bAvg)
                    "b1 equals the supplied function's Ala mass"
                Expect.isFalse
                    (containsMass 0.0005 71.03711 bAvg)
                    "the monoisotopic b1 is absent from the average-mass series"
                // parametricity: the ladder must consume the caller's mass model; a hard-coded monoisotopic table would pass every existing test and fail this one. Average - monoisotopic Ala differ by ~0.041 Da, far beyond tolerance.

            // PENDING: backbone cleavage of a 3-residue peptide produces exactly 2 fragment positions per
            // series; the full-length b3/y3 "fragment" is the intact peptide (y_n equals the precursor
            // neutral mass) and is not a cleavage product - its presence inflates every theoretical
            // spectrum. The implementation emits n families per series.
            ptestCase "a fragment ladder of an n-residue peptide has n-1 cleavage fragments" <| fun _ ->
                Expect.equal (List.length (Fragmentation.Series.bOfBioList mf ags)) 2 "a 3-mer b series has two cleavage fragments"

            testCase "y-series masses are C-terminal residue sums plus water" <| fun _ ->
                let masses = mainMasses (Fragmentation.Series.yOfBioList mf ags)
                Expect.isTrue
                    (containsMass 0.001 105.04260 masses)
                    (sprintf "y1 is Ser + H2O; expected a mass within %g of %g; got %A" 0.001 105.04260 masses)
                Expect.isTrue
                    (containsMass 0.001 162.06406 masses)
                    (sprintf "y2 is Ser + Gly + H2O; expected a mass within %g of %g; got %A" 0.001 162.06406 masses)
                // y ions retain the C-terminus and the full water of the intact peptide - textbook.

            testCase "a and c ions differ from b by CO and NH3" <| fun _ ->
                let aMasses = mainMasses (Fragmentation.Series.aOfBioList mf ags)
                let cMasses = mainMasses (Fragmentation.Series.cOfBioList mf ags)
                Expect.isTrue
                    (containsMass 0.001 (71.03711 - 27.994915) aMasses)
                    (sprintf "a1 is b1 - CO; expected a mass within %g of %g; got %A" 0.001 (71.03711 - 27.994915) aMasses)
                Expect.isTrue
                    (containsMass 0.001 (71.03711 + 17.026549) cMasses)
                    (sprintf "c1 is b1 + NH3; expected a mass within %g of %g; got %A" 0.001 (71.03711 + 17.026549) cMasses)
                Expect.isTrue
                    (containsMass 0.001 100.063655 aMasses)
                    (sprintf "a2 is b2 - CO; expected a mass within %g of %g; got %A" 0.001 100.063655 aMasses)
                Expect.isTrue
                    (containsMass 0.001 145.085119 cMasses)
                    (sprintf "c2 is b2 + NH3; expected a mass within %g of %g; got %A" 0.001 145.085119 cMasses)
                // a and c ion relations use the published carbon monoxide and ammonia masses.

            testCase "an oxidized residue shifts exactly the fragments containing it" <| fun _ ->
                let oxMet = AminoAcids.setModification (SearchDB.getModBy SearchDB.Table.oxidation'Met') AminoAcids.Met
                let plain = [AminoAcids.Ala; AminoAcids.Met; AminoAcids.Gly]
                let modified = [AminoAcids.Ala; oxMet; AminoAcids.Gly]
                let bPlain = mainMasses (Fragmentation.Series.bOfBioList mf plain) |> List.sort
                let bModified = mainMasses (Fragmentation.Series.bOfBioList mf modified) |> List.sort
                let yPlain = mainMasses (Fragmentation.Series.yOfBioList mf plain) |> List.sort
                let yModified = mainMasses (Fragmentation.Series.yOfBioList mf modified) |> List.sort
                let b1Plain, b2Plain = bPlain.[0], bPlain.[1]
                let b1Modified, b2Modified = bModified.[0], bModified.[1]
                let y1Plain, y2Plain = yPlain.[0], yPlain.[1]
                let y1Modified, y2Modified = yModified.[0], yModified.[1]

                expectWithin 1e-9 b1Modified b1Plain "b1 does not contain the modified Met"
                expectWithin 1e-9 y1Modified y1Plain "y1 does not contain the modified Met"
                expectWithin 0.001 (b2Modified - b2Plain) 15.994915 "b2 carries the oxidation delta"
                expectWithin 0.001 (y2Modified - y2Plain) 15.994915 "y2 carries the oxidation delta"
                // published Unimod oxidation delta; a PTM must shift precisely the fragments that contain the modified residue - the positional signature every search engine relies on.

            testCase "z ions differ from y by ammonia" <| fun _ ->
                let zMasses = mainMasses (Fragmentation.Series.zOfBioList mf ags)
                Expect.isTrue
                    (containsMass 0.001 (105.04260 - 17.026549) zMasses)
                    (sprintf "z1 is y1 - NH3; expected a mass within %g of %g; got %A" 0.001 (105.04260 - 17.026549) zMasses)
                Expect.isTrue
                    (containsMass 0.001 (162.06406 - 17.026549) zMasses)
                    (sprintf "z2 is y2 - NH3; expected a mass within %g of %g; got %A" 0.001 (162.06406 - 17.026549) zMasses)
                // z ions differ from the corresponding y ion by the published ammonia mass.
                // this is the classic Biemann z ion (y - NH3), deliberately NOT the ETD z-dot (z + 1.00783); the API claims no ETD semantics, so the Biemann convention is the correct passing expectation.

            // PENDING: textbook x-ion mass is y + CO - H2 (= y + 25.979265); the implementation applies a
            // plain +CO modification (y + 27.994915), 2.01565 Da (one H2) too heavy.
            ptestCase "x ions differ from y by CO minus H2" <| fun _ ->
                let xMasses = mainMasses (Fragmentation.Series.xOfBioList mf ags)
                Expect.isTrue
                    (containsMass 0.001 (105.04260 + 25.979265) xMasses)
                    (sprintf "x1 is y1 + CO - H2; expected a mass within %g of %g; got %A" 0.001 (105.04260 + 25.979265) xMasses)
                // x = y + CO - H2, using published CO and H2 masses: 27.994915 - 2.015650 = 25.979265.
        ]

        testList "NeutralLosses" [
            testCase "a water-loss residue triggers an H2O-loss satellite from its first fragment onwards" <| fun _ ->
                let famsSA = Fragmentation.Series.bOfBioList mf [AminoAcids.Ser; AminoAcids.Ala]
                let famsAS = Fragmentation.Series.bOfBioList mf [AminoAcids.Ala; AminoAcids.Ser]
                let serFirst = familyWithMainMass 0.001 87.03203 famsSA
                let serAla = familyWithMainMass 0.001 158.06914 famsSA
                let alaFirst = familyWithMainMass 0.001 71.03711 famsAS
                let alaSer = familyWithMainMass 0.001 158.06914 famsAS
                Expect.isTrue
                    (serAla.DependentPeaks
                     |> List.exists (fun peak ->
                         Ions.hasFlag peak.Iontype Ions.IonTypeFlag.lossH2O
                         && abs (peak.Mass - 140.058575) <= 0.001))
                    (sprintf "the Ser + Ala fragment has a water-loss satellite; expected a mass within %g of %g; got %A" 0.001 140.058575 (serAla.DependentPeaks |> List.map (fun peak -> peak.Mass)))
                // the loss channel PERSISTS: Ser is no longer the current residue at b2, but once a loss-prone residue is inside the fragment every longer fragment keeps the channel.
                Expect.isTrue
                    (alaSer.DependentPeaks
                     |> List.exists (fun peak ->
                         Ions.hasFlag peak.Iontype Ions.IonTypeFlag.lossH2O
                         && abs (peak.Mass - (158.06914 - 18.010565)) <= 0.001))
                    "the Ala + Ser fragment has a water-loss satellite"
                Expect.isTrue
                    (serFirst.DependentPeaks
                     |> List.exists (fun peak ->
                         Ions.hasFlag peak.Iontype Ions.IonTypeFlag.lossH2O
                         && abs (peak.Mass - (87.03203 - 18.010565)) <= 0.001))
                    "the Ser-first fragment has a water-loss satellite"
                Expect.isEmpty alaFirst.DependentPeaks "the Ala-first fragment has no loss-prone residue yet"
                // Ser side chains are prone to water loss; once a loss-prone residue is part of the fragment, every longer fragment retains the loss channel - documented residue-loss convention plus published masses.

            testCase "an ammonia-loss residue triggers an NH3-loss satellite" <| fun _ ->
                let fams = Fragmentation.Series.bOfBioList mf [AminoAcids.Lys; AminoAcids.Ala]
                let lysFirst = familyWithMainMass 0.001 128.09496 fams
                Expect.isTrue
                    (lysFirst.DependentPeaks
                     |> List.exists (fun peak ->
                         Ions.hasFlag peak.Iontype Ions.IonTypeFlag.lossNH3
                         && abs (peak.Mass - (128.09496 - 17.026549)) <= 0.001))
                    "the Lys fragment has an ammonia-loss satellite"
                // Lys is in the documented ammonia-loss set; published NH3 mass.

            testCase "the residue loss table is complete and exclusive" <| fun _ ->
                let assertOneLoss residue lossFlag lossMass =
                    let family =
                        Fragmentation.Series.bOfBioList mf [residue; AminoAcids.Ala]
                        |> List.minBy (fun family -> family.MainPeak.Mass)
                    Expect.equal family.DependentPeaks.Length 1 "the b1 family has exactly one dependent"
                    let dependent = List.head family.DependentPeaks
                    Expect.isTrue
                        (Ions.hasFlag dependent.Iontype lossFlag)
                        (sprintf "the dependent carries the expected loss flag %A" lossFlag)
                    expectWithin
                        0.001
                        dependent.Mass
                        (family.MainPeak.Mass - lossMass)
                        "the dependent has the published neutral-loss mass"

                [ (AminoAcids.Ser, Ions.IonTypeFlag.lossH2O, 18.010565)
                  (AminoAcids.Thr, Ions.IonTypeFlag.lossH2O, 18.010565)
                  (AminoAcids.Glu, Ions.IonTypeFlag.lossH2O, 18.010565)
                  (AminoAcids.Asp, Ions.IonTypeFlag.lossH2O, 18.010565)
                  (AminoAcids.Arg, Ions.IonTypeFlag.lossNH3, 17.026549)
                  (AminoAcids.Lys, Ions.IonTypeFlag.lossNH3, 17.026549)
                  (AminoAcids.Gln, Ions.IonTypeFlag.lossNH3, 17.026549)
                  (AminoAcids.Asn, Ions.IonTypeFlag.lossNH3, 17.026549) ]
                |> List.iter (fun (residue, lossFlag, lossMass) -> assertOneLoss residue lossFlag lossMass)

                [AminoAcids.Ala; AminoAcids.Gly]
                |> List.iter (fun residue ->
                    let family =
                        Fragmentation.Series.bOfBioList mf [residue; AminoAcids.Ala]
                        |> List.minBy (fun family -> family.MainPeak.Mass)
                    Expect.isEmpty family.DependentPeaks "a residue outside both loss tables has no dependents")
                // the documented residue-loss convention, table-driven: missing residues, spurious channels, and misassigned channels all fail; published loss masses.

            testCase "water and ammonia loss channels coexist on one fragment" <| fun _ ->
                let fams = Fragmentation.Series.bOfBioList mf [AminoAcids.Ser; AminoAcids.Lys; AminoAcids.Ala]
                let b2 = fams |> List.sortBy (fun family -> family.MainPeak.Mass) |> List.item 1
                expectWithin 0.001 b2.MainPeak.Mass (87.03203 + 128.09496) "the b2 main mass is Ser plus Lys"
                Expect.equal b2.DependentPeaks.Length 2 "the b2 family has both loss channels"
                let ammoniaLoss =
                    b2.DependentPeaks
                    |> List.find (fun peak -> Ions.hasFlag peak.Iontype Ions.IonTypeFlag.lossNH3)
                let waterLoss =
                    b2.DependentPeaks
                    |> List.find (fun peak -> Ions.hasFlag peak.Iontype Ions.IonTypeFlag.lossH2O)
                expectWithin 0.001 ammoniaLoss.Mass (b2.MainPeak.Mass - 17.026549) "the ammonia-loss satellite is present"
                expectWithin 0.001 waterLoss.Mass (b2.MainPeak.Mass - 18.010565) "the water-loss satellite is present"
                // once both a water-prone and an ammonia-prone residue are inside the fragment, both channels must be offered - published loss masses, hand-summed main.

            testCase "y-series loss channels follow C-terminal suffix membership" <| fun _ ->
                let yASG = Fragmentation.Series.yOfBioList mf [AminoAcids.Ala; AminoAcids.Ser; AminoAcids.Gly]
                let ySAG = Fragmentation.Series.yOfBioList mf [AminoAcids.Ser; AminoAcids.Ala; AminoAcids.Gly]
                let y1ASG = yASG |> List.sortBy (fun family -> family.MainPeak.Mass) |> List.item 0
                let y2ASG = yASG |> List.sortBy (fun family -> family.MainPeak.Mass) |> List.item 1
                let y1SAG = ySAG |> List.sortBy (fun family -> family.MainPeak.Mass) |> List.item 0
                let y2SAG = ySAG |> List.sortBy (fun family -> family.MainPeak.Mass) |> List.item 1
                expectWithin 0.001 y1ASG.MainPeak.Mass 75.03203 "y1 is Gly plus water"
                expectWithin 0.001 y1SAG.MainPeak.Mass 75.03203 "y1 is Gly plus water"
                Expect.isEmpty y1ASG.DependentPeaks "y1 has no loss-prone residue in its suffix"
                Expect.isEmpty y1SAG.DependentPeaks "y1 has no loss-prone residue in its suffix"
                expectWithin 0.001 y2ASG.MainPeak.Mass 162.06406 "y2 is Ser plus Gly plus water"
                Expect.isTrue
                    (y2ASG.DependentPeaks |> List.exists (fun peak -> Ions.hasFlag peak.Iontype Ions.IonTypeFlag.lossH2O))
                    "y2 has a water-loss channel when Ser enters the suffix"
                expectWithin 0.001 y2SAG.MainPeak.Mass 146.06914 "y2 is Ala plus Gly plus water"
                Expect.isEmpty y2SAG.DependentPeaks "y2 has no loss channel when Ser is outside the suffix"
                // y ions grow from the C-terminus: the loss channel activates exactly when the loss-prone residue enters the suffix - hand-placed by construction.
        ]

        testList "TargetDecoy" [
            testCase "decoy fragment masses equal target fragment masses of the reversed sequence" <| fun _ ->
                let fm = Fragmentation.Series.fragmentMasses Fragmentation.Series.bOfBioList Fragmentation.Series.yOfBioList mf ags
                let fmRev = Fragmentation.Series.fragmentMasses Fragmentation.Series.bOfBioList Fragmentation.Series.yOfBioList mf (List.rev ags)
                Expect.equal fm.DecoyMasses.Length fmRev.TargetMasses.Length "target and decoy fragment counts match"
                List.iter2
                    (fun (actual: PeakFamily<TaggedMass.TaggedMass>) (expected: PeakFamily<TaggedMass.TaggedMass>) ->
                        expectWithin 1e-9 actual.MainPeak.Mass expected.MainPeak.Mass "decoy and reversed-target masses match")
                    fm.DecoyMasses
                    fmRev.TargetMasses
                List.iter2
                    (fun (actual: PeakFamily<TaggedMass.TaggedMass>) (expected: PeakFamily<TaggedMass.TaggedMass>) ->
                        Expect.equal actual.DependentPeaks.Length expected.DependentPeaks.Length "decoy and reversed-target dependent peak counts match"
                        List.iter2
                            (fun (actualPeak: TaggedMass.TaggedMass) (expectedPeak: TaggedMass.TaggedMass) ->
                                expectWithin 1e-9 actualPeak.Mass expectedPeak.Mass "decoy and reversed-target dependent masses match")
                            actual.DependentPeaks
                            expected.DependentPeaks)
                    fm.DecoyMasses
                    fmRev.TargetMasses
                // The decoy is the standard reversed-decoy convention; its fragment set must be exactly the target fragment set of the reversed input - cross-consistency, no mass constants involved.
        ]

        testList "Laddering" [
            testCase "ladderAndChargeElement converts neutral masses to m/z per charge" <| fun _ ->
                let fam = Peaks.createPeakFamily (TaggedMass.createTaggedMass Ions.IonTypeFlag.B 200.0) [TaggedMass.createTaggedH2OLoss Ions.IonTypeFlag.B 182.0]
                let r = Fragmentation.ladderAndChargeElement [1.0; 2.0] [fam]
                Expect.equal r.Length 2 "one family is returned for each charge"
                let chargeOne = r |> List.find (fun family -> family.MainPeak.Charge = 1.0)
                let chargeTwo = r |> List.find (fun family -> family.MainPeak.Charge = 2.0)
                let chargeOneDependent = chargeOne.DependentPeaks |> List.head
                let chargeTwoDependent = chargeTwo.DependentPeaks |> List.head
                expectWithin 0.001 chargeOne.MainPeak.MassOverCharge 201.00728 "charge-1 m/z is 200 + proton"
                Expect.equal chargeOne.MainPeak.Charge 1.0 "charge-1 charge is preserved"
                Expect.equal chargeOne.MainPeak.Number 1 "charge-1 number is one"
                Expect.equal chargeOne.MainPeak.Iontype Ions.IonTypeFlag.B "charge-1 ion type is B"
                expectWithin 0.001 chargeTwo.MainPeak.MassOverCharge 101.00728 "charge-2 m/z is 200 / 2 + proton"
                Expect.equal chargeTwo.MainPeak.Charge 2.0 "charge-2 charge is preserved"
                Expect.equal chargeTwo.MainPeak.Number 1 "charge-2 number is one"
                Expect.equal chargeOne.DependentPeaks.Length 1 "charge-1 has one dependent peak"
                expectWithin 0.001 chargeOneDependent.MassOverCharge 183.007276 "charge-1 dependent m/z is 182 + proton"
                Expect.equal chargeOneDependent.Number chargeOne.MainPeak.Number "charge-1 dependent number matches the main peak"
                Expect.isTrue
                    (Ions.hasFlag chargeOneDependent.Iontype Ions.IonTypeFlag.B
                     && Ions.hasFlag chargeOneDependent.Iontype Ions.IonTypeFlag.lossH2O)
                    (sprintf "charge-1 dependent ion type carries B and lossH2O; got %A" chargeOneDependent.Iontype)
                Expect.equal chargeTwo.DependentPeaks.Length 1 "charge-2 has one dependent peak"
                expectWithin 0.001 chargeTwoDependent.MassOverCharge 92.007276 "charge-2 dependent m/z is 182 / 2 + proton"
                Expect.equal chargeTwoDependent.Number chargeTwo.MainPeak.Number "charge-2 dependent number matches the main peak"
                Expect.isTrue
                    (Ions.hasFlag chargeTwoDependent.Iontype Ions.IonTypeFlag.B
                     && Ions.hasFlag chargeTwoDependent.Iontype Ions.IonTypeFlag.lossH2O)
                    (sprintf "charge-2 dependent ion type carries B and lossH2O; got %A" chargeTwoDependent.Iontype)
                // The m/z of an ion is (neutral mass + z protons)/z - fundamental MS relation with the published proton mass.

            testCase "ladderElement numbers each ion series independently starting at one" <| fun _ ->
                let bFam1 = Peaks.createPeakFamily (TaggedMass.createTaggedMass Ions.IonTypeFlag.B 100.0) []
                let bFam2 = Peaks.createPeakFamily (TaggedMass.createTaggedMass Ions.IonTypeFlag.B 250.0) []
                let yFam = Peaks.createPeakFamily (TaggedMass.createTaggedMass Ions.IonTypeFlag.Y 150.0) []
                let r = Fragmentation.ladderElement [bFam2; yFam; bFam1] [1.0]
                let bFamilies = r |> List.filter (fun family -> family.MainPeak.Iontype = Ions.IonTypeFlag.B)
                let yFamily = r |> List.find (fun family -> family.MainPeak.Iontype = Ions.IonTypeFlag.Y)
                let bLighter = bFamilies |> List.minBy (fun family -> family.MainPeak.MassOverCharge)
                let bHeavier = bFamilies |> List.maxBy (fun family -> family.MainPeak.MassOverCharge)
                Expect.equal bLighter.MainPeak.Number 1 "the lighter B ion is b1"
                Expect.equal bHeavier.MainPeak.Number 2 "the heavier B ion is b2"
                Expect.equal yFamily.MainPeak.Number 1 "the Y ion series starts at y1"
                // Ladder numbering (b1, b2, ..., y1, ...) is per ion series in ascending mass order, starting at 1 - the standard fragment-ladder naming convention.
        ]

        testList "TheoreticalSpectra" [
            testCase "getTheoSpec feeds target and decoy fragment masses through the predictor and preserves the lookup" <| fun _ ->
                let lookup = SearchDB.createLookUpResult 7 42 328.15 328150000L "AGS" ags 0
                let fm = Fragmentation.Series.fragmentMasses Fragmentation.Series.bOfBioList Fragmentation.Series.yOfBioList mf ags
                let predictor = fun scanlimits charge masses -> (scanlimits, charge, masses |> List.length, masses)
                let spec = TheoreticalSpectra.getTheoSpec predictor (100.0, 1000.0) 2 (lookup, fm)
                Expect.equal spec.LookUpResult lookup "the lookup result is preserved"
                let sl, ch, n, ms = spec.TheoSpec
                Expect.equal sl (100.0, 1000.0) "target scan limits are forwarded"
                Expect.equal ch 2.0 "target charge is forwarded as a float"
                Expect.equal n fm.TargetMasses.Length "target mass count is forwarded"
                Expect.equal ms fm.TargetMasses "target masses are forwarded"
                let decoyScanLimits, decoyCharge, decoyN, decoyMs = spec.DecoyTheoSpec
                Expect.equal decoyScanLimits (100.0, 1000.0) "decoy scan limits are forwarded"
                Expect.equal decoyCharge 2.0 "decoy charge is forwarded as a float"
                Expect.equal decoyN fm.DecoyMasses.Length "decoy mass count is forwarded"
                Expect.equal decoyMs fm.DecoyMasses "decoy masses are forwarded"
                // the decoy prediction must receive the SAME acquisition context as the target - a decoy scored under different limits or charge silently biases FDR.

            testCase "getTheoSpecs yields one spectrum per candidate" <| fun _ ->
                let lookup = SearchDB.createLookUpResult 7 42 328.15 328150000L "AGS" ags 0
                let fm = Fragmentation.Series.fragmentMasses Fragmentation.Series.bOfBioList Fragmentation.Series.yOfBioList mf ags
                let lookup2 = SearchDB.createLookUpResult 8 43 203.12 203120000L "AG" [AminoAcids.Ala; AminoAcids.Gly] 0
                let fm2 = Fragmentation.Series.fragmentMasses Fragmentation.Series.bOfBioList Fragmentation.Series.yOfBioList mf [AminoAcids.Ala; AminoAcids.Gly]
                let specs =
                    TheoreticalSpectra.getTheoSpecs
                        (fun _ _ ms -> ms)
                        (100.0, 1000.0)
                        2
                        [(lookup, fm); (lookup2, fm2)]
                Expect.equal specs.Length 2 "one spectrum is returned per candidate"
                let ids = specs |> List.map (fun spec -> spec.LookUpResult.PepSequenceID) |> Set.ofList
                Expect.equal ids (set [42; 43]) "candidate membership is preserved"
                specs
                |> List.iter (fun spec ->
                    let candidate =
                        [(lookup, fm); (lookup2, fm2)]
                        |> List.find (fun (candidateLookup, _) -> candidateLookup.PepSequenceID = spec.LookUpResult.PepSequenceID)
                    let _, candidateMasses = candidate
                    let actualTargetMasses = mainMasses spec.TheoSpec
                    let expectedTargetMasses = mainMasses candidateMasses.TargetMasses
                    let actualDecoyMasses = mainMasses spec.DecoyTheoSpec
                    let expectedDecoyMasses = mainMasses candidateMasses.DecoyMasses
                    Expect.equal actualTargetMasses.Length expectedTargetMasses.Length "target payload lengths match the candidate"
                    List.iter2
                        (fun actual expected -> expectWithin 1e-9 actual expected "target payload belongs to the candidate")
                        actualTargetMasses
                        expectedTargetMasses
                    Expect.equal actualDecoyMasses.Length expectedDecoyMasses.Length "decoy payload lengths match the candidate"
                    List.iter2
                        (fun actual expected -> expectWithin 1e-9 actual expected "decoy payload belongs to the candidate")
                        actualDecoyMasses
                        expectedDecoyMasses)
                Expect.equal
                    (TheoreticalSpectra.getTheoSpecs (fun _ _ ms -> 0) (100.0, 1000.0) 2 [])
                    []
                    "an empty candidate list returns no spectra"
                // count and ID membership alone would pass if payloads were swapped between candidates or one candidate's spectra were reused; per-ID payload identity is the actual association contract.
        ]
    ]
