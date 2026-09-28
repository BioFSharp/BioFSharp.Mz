namespace BioFSharp.Mz

open BioFSharp
open SearchDB
open Fragmentation
open TheoreticalSpectra
open FSharpAux
open SearchEngineResult
open Ions
open FSharp.Stats
//open MathNet.Numerics
//open MathNet.Numerics.LinearAlgebra.Double

module SequestLike =

    /// normalize the intensities within a window to maximum of the window
    /// Attention shortens the array  (cuts)
    let windowNormalizeIntensities (intensities:Vector<float>) (numberOfWindows:int) =
        // finds max within range
        let rec findMax (array:Vector<float>) (cMax:float) (lowerLimit:int) (counter:int)  =
            if counter < lowerLimit then
                cMax
            else
                let nMax = max array.[counter] cMax
                findMax array nMax lowerLimit (counter - 1)
        // sqrt intensities and normalize to maximum within lower-upper limit
        let normToMaxSqrtInPlace (array:Vector<float>) (lowerLimit:int) (upperLimit:int)  =
            let nMax = sqrt (findMax array 0.0 lowerLimit upperLimit)
            for i = lowerLimit to upperLimit do
                if nMax > 0. then                
                    array.[i] <- (sqrt array.[i]) / nMax 
                else
                    array.[i] <- 0.
        let windowSize =  (intensities.Length / numberOfWindows)    
        let tmpIntensities =            
            Vector.init (windowSize*numberOfWindows) (fun i -> intensities.[i])            
        for i = 1 to numberOfWindows do
            //printfn "window: %i lower: %i counter: %i " i (windowSize * (i - 1)) (windowSize * i - 1)
            normToMaxSqrtInPlace tmpIntensities (windowSize * (i - 1)) (windowSize * i - 1) 
        tmpIntensities
    
    /// Predicts the intensity of a theoretical peak based on the given charge and iontype.
    let predictIntensitySimpleModel (ionType:Ions.IonTypeFlag) (charge:float) =
        match ionType with 
        | ionT when hasFlag Ions.IonTypeFlag.Unknown    ionT -> 1.  / charge
        | ionT when hasFlag Ions.IonTypeFlag.Diagnostic ionT -> 1.  / charge
        | ionT when hasFlag Ions.IonTypeFlag.Neutral    ionT -> 0.6 / charge
        | ionT when hasFlag Ions.IonTypeFlag.Immonium   ionT -> 0.6 / charge
        | ionT when hasFlag Ions.IonTypeFlag.lossNH3    ionT -> 0.2 / charge
        | ionT when hasFlag Ions.IonTypeFlag.lossH2O    ionT -> 0.2 / charge
        | ionT when hasFlag Ions.IonTypeFlag.A          ionT -> 0.2 / charge   
        | ionT when hasFlag Ions.IonTypeFlag.B          ionT -> 1.  / charge
        | ionT when hasFlag Ions.IonTypeFlag.C          ionT -> 0.2 / charge
        | ionT when hasFlag Ions.IonTypeFlag.X          ionT -> 0.2 / charge
        | ionT when hasFlag Ions.IonTypeFlag.Y          ionT -> 1.  / charge
        | ionT when hasFlag Ions.IonTypeFlag.Z          ionT -> 0.2 / charge
        | _                                                  -> 0.2 / charge

    let setVectorInplace (vector:Vector<float>) charge maxIndex lowerScanLimit (taggedMass:TaggedMass.TaggedMass) = 
        let index = int(System.Math.Round (Mass.toMZ taggedMass.Mass charge) ) - lowerScanLimit 
        if index < maxIndex-1 && index > -1 then
            vector.[index] <- max vector.[index] (predictIntensitySimpleModel (taggedMass.Iontype) charge) 

    /// Converts the fragment ion ladder to a theoretical Sequestlike spectrum at a given charge state. 
    let theoSpecOf (lowerScanLimit,upperScanLimit) (maxcharge:float) (fragments:PeakFamily<TaggedMass.TaggedMass> list) =
        let lowerScanLimit = int lowerScanLimit
        let upperScanLimit = int upperScanLimit
        let maxIndex = upperScanLimit - lowerScanLimit + 1        
        let vector = Vector.create (maxIndex-1) 0.
        fragments 
        |> List.iter (fun p ->  
            let peaks = p.MainPeak::p.DependentPeaks
            [1. .. maxcharge]  
            |> List.iter (fun ch -> peaks |> List.iter (setVectorInplace vector ch maxIndex lowerScanLimit) ) 
            )
        vector

    /// Computes the autocorrelation of the vector +/- plusMinusMaxDelay.
    let autoCorrelation (plusMinusMaxDelay:int) (vector:Vector<float>) =
        // Element by element the same sums as the former recursion, which stopped one shift short
        // on both sides. The positive shifts run from +delay down to +2. The negative shifts run
        // from -1 down to -(delay-1). The two sums are added and divided by twice the delay. A
        // shift outside the vector added a zero, and adding zero changes nothing, so those shifts
        // are skipped. A delay below 1 made the former recursion run forever, so it is rejected.
        if plusMinusMaxDelay < 1 then invalidArg "plusMinusMaxDelay" "Must be at least 1."
        let n = vector.Length
        let result = Vector.zero n
        let scale = float plusMinusMaxDelay * 2.
        for i = 0 to n - 1 do
            let mutable plus = 0.
            let mutable tau = min plusMinusMaxDelay i
            while tau >= 2 do
                plus <- plus + vector.[i - tau]
                tau <- tau - 1
            let mutable minus = 0.
            let mutable k = 1
            let kMax = min (plusMinusMaxDelay - 1) (n - 1 - i)
            while k <= kMax do
                minus <- minus + vector.[i + k]
                k <- k + 1
            result.[i] <- (plus + minus) / scale
        result

    /// Converts the fragment ion ladder to a theoretical Sequestlike spectrum at a given charge state. 
    /// Subsequently, the spectrum is binned to the nearest mz bin (binwidth = 1 Da). Filters out peaks 
    /// that are not within the scanLimits.
    let predictOf (lowerScanLimit,upperScanLimit) charge ionSeries =         
        theoSpecOf (lowerScanLimit,upperScanLimit) charge ionSeries      
        

//    /// Measured spectrum to sequest-like normalized intensity array
//    /// ! Uses 10 as number of windows for window normalization (like in original sequest algorithm)    
//    let spectrumToNormalizedIntensityArray (scanlimits:int*int) (spectrum:PeakArray<_>) =
//        let lowerScanLimit,upperScanLimit = scanlimits    
//        let si  = PeakArray.peaksToNearestUnitDaltonBin spectrum lowerScanLimit upperScanLimit
//        let nsi = windowNormalizeIntensities si 10 // |> Seq.toArray
//        nsi

    /// Measured spectrum to sequest-like normalized intensity array
    /// minus auto-correlation (delay 75 -> like in original sequest algorithm)
    /// ! Uses 10 as number of windows for window normalization (like in original sequest algorithm)    
    let spectrumToIntensityArrayMinusAutoCorrelation (lowerScanLimit,upperScanLimit) (spectrum:PeakArray<_>) =
        let si  = PeakArray.peaksToNearestUnitDaltonBinVector spectrum lowerScanLimit upperScanLimit
        let nsi = windowNormalizeIntensities si 10    
        let nsi' = autoCorrelation 75 nsi
        (nsi - nsi')  


    /// Calculates the Cross-Correlation of p_nis and ms_nis
    let calcXCorr (p_nis:Vector<float>) (ms_nis:Vector<float>) =
        let tmp = Vector.dot p_nis ms_nis 
        if tmp < 0. then 0. else tmp 

    /// Substracts Auto Correlation of measured normalized intensity spectrum (ms_nis). Subsequently, 
    /// the Cross-Correlation of the predicted normalized intensity spectrum (p_nis) and ms_nis is calculated. 
    let scoreSingle (p_nis:Vector<float>) (ms_nis:Vector<float>)  =
        let ms_nis' = autoCorrelation 75 ms_nis
        (ms_nis - ms_nis') 
        |> calcXCorr p_nis 
        
    /// Converts the fragment ion ladders to a theoretical Sequestlike spectrum at a given charge state. 
    /// Subsequently, the spectrum is binned to the nearest mz bin (binwidth = 1 Da). Filters out peaks 
    /// that are not within the scanLimits.
    let getTheoSpecs scanlimits chargeState (possiblePeptideInfos:list<LookUpResult<AminoAcids.AminoAcid>*FragmentMasses>) =
        TheoreticalSpectra.getTheoSpecs predictOf scanlimits chargeState possiblePeptideInfos

    /// Calculates the SequestLike Scores for all peptides (possiblePeptideInfos).
    let calcSequestLikeScoresRevDecoy calcIonSeries (massfunction:Formula.Formula -> float) (scanlimits) (spectrum:PeakArray<_>) scanTime chargeState isolationWindowTargetMz (possiblePeptideInfos:list<LookUpResult<AminoAcids.AminoAcid>>) spectrumID =                             
        // measured normailzed intensity array (spectrum) minus auto-correlation
        let ms_nis =  spectrumToIntensityArrayMinusAutoCorrelation scanlimits spectrum
        // float charge
        let fCharge = float chargeState
        // measured mass
        let ms_mass = Mass.ofMZ isolationWindowTargetMz fCharge
        //
        let ides = 
            possiblePeptideInfos 
            |> List.fold (fun acc lookUpResult -> 
                              let sequence = lookUpResult.BioSequence        
                              let ionSeries = calcIonSeries massfunction sequence    
                              //predicted  normalized intensity array (spectrum) 
                              let p_nis = predictOf scanlimits fCharge ionSeries 
                              let xcorr = calcXCorr p_nis ms_nis
                              let targetScore = createSearchEngineResult SearchEngineResult.SearchEngine.SEQUESTLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod true scanTime lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass  (List.length sequence) xcorr nan nan
                              
                              let revPeptide_decoy = sequence |> List.rev
                              let ionSeries_decoy  = calcIonSeries massfunction revPeptide_decoy
                              let p_nis_decoy      = predictOf scanlimits fCharge ionSeries_decoy 
                              let xcorr_decoy      = calcXCorr p_nis_decoy ms_nis
                              let decoyScore = createSearchEngineResult SearchEngineResult.SearchEngine.SEQUESTLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod false scanTime lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass revPeptide_decoy.Length xcorr_decoy nan nan
                              targetScore::decoyScore::acc  
                     ) []
        calcNormDeltaBestToRest (ides |> List.sortBy (fun sls -> - sls.Score))        
        |> calcNormDeltaNext

    /// Calculates the SequestLike Scores for all theoretical spectra.    
    let calcSequestScore scanlimits (spectrum:PeakArray<_>) scanTime chargeState isolationWindowTargetMz (theoreticalSpectra:list<TheoreticalSpectrum<Vector<float>>>) spectrumID = 
        let fCharge = float chargeState
        // measured normailzed intensity array (spectrum) minus auto-correlation
        let ms_nis =  spectrumToIntensityArrayMinusAutoCorrelation scanlimits spectrum
        // measured mass
        let ms_mass = Mass.ofMZ isolationWindowTargetMz fCharge        

        let ides = 
            theoreticalSpectra 
            |> List.fold (fun acc theoreticalSpectrum -> 
                              let lookUpResult = theoreticalSpectrum.LookUpResult
                              let sequence = lookUpResult.BioSequence
                              let p_nis = theoreticalSpectrum.TheoSpec 
                              let xcorr = calcXCorr p_nis ms_nis
                              let targetScore = createSearchEngineResult SearchEngineResult.SearchEngine.SEQUESTLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod true scanTime lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass sequence.Length xcorr nan nan
                              
                              let p_nis_decoy = theoreticalSpectrum.DecoyTheoSpec
                              let xcorr_decoy = calcXCorr p_nis_decoy  ms_nis
                              let decoyScore = createSearchEngineResult SearchEngineResult.SearchEngine.SEQUESTLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod false scanTime lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass sequence.Length xcorr_decoy nan nan
                              targetScore::decoyScore::acc  
                     ) []

        calcNormDeltaBestToRest (ides  |> List.sortBy (fun sls -> - sls.Score))        
        |> calcNormDeltaNext
        
    /// Calculates the sequestLike Scores for all theoretical spectra. Implemented using Async parallel. 
    let calcSequestScoreParallel scanlimits (spectrum:PeakArray<_>) scanTime chargeState isolationWindowTargetMz (theoreticalSpectra:list<TheoreticalSpectrum<Vector<float>>>) spectrumID = 
        let fCharge = float chargeState
        // measured normailzed intensity array (spectrum) minus auto-correlation
        let ms_nis =  spectrumToIntensityArrayMinusAutoCorrelation scanlimits spectrum
        // measured mass
        let ms_mass = Mass.ofMZ isolationWindowTargetMz fCharge     
        let computeScores theoreticalSpectrum =
            let lookUpResult = theoreticalSpectrum.LookUpResult
            let sequence = lookUpResult.BioSequence
            let p_nis = theoreticalSpectrum.TheoSpec 
            let xcorr = calcXCorr p_nis ms_nis
            let targetScore = createSearchEngineResult SearchEngineResult.SearchEngine.SEQUESTLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod true scanTime lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass sequence.Length xcorr nan nan
                              
            let p_nis_decoy      = theoreticalSpectrum.DecoyTheoSpec
            let xcorr_decoy      = calcXCorr p_nis_decoy ms_nis
            let decoyScore = createSearchEngineResult SearchEngineResult.SearchEngine.SEQUESTLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod false scanTime lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass sequence.Length xcorr_decoy nan nan
            [targetScore;decoyScore]
        
        let ides = 
            [for i in theoreticalSpectra ->  async {return computeScores i}]
            |> Async.Parallel 
            |> Async.RunSynchronously
            |> List.concat    
        calcNormDeltaBestToRest (ides  |> List.sortBy (fun sls -> - sls.Score))        
        |> calcNormDeltaNext