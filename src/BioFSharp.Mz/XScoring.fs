namespace BioFSharp.Mz

open System
open BioFSharp
open BioFSharp.Mz.Peaks
open SearchDB
open Fragmentation
open TheoreticalSpectra
open SearchEngineResult
open FSharpAux.Array
open FSharp.Stats.SpecialFunctions.Factorial

module XScoring =
    
    type MatchingScore = {
        Score   : float 
        N       : int 
        K       : int
        Q       : int
        //Matches : PeakAnnotation list 
        }
    
    let createMatchingScore score n k q = {
     Score = score; N = n ; K = k; Q = q ;//Matches = matches
    }

    [<Struct>]
    type RatedPeak = {
        Peak: float*float
        MoreIntensePeaksWithinWindow: int
        }

    let createRatedPeak peak moreIntensePeaks = {
        Peak = peak; MoreIntensePeaksWithinWindow= moreIntensePeaks }

    let getIntensity rPeak =
        snd rPeak.Peak 
    /// Counts how many peaks present in spectrum data within a a given mzWindow (mzwindow/2 +/- ) have a higher 
    /// intensity then the peak at position spectrumData.[peakIdx]. 
    let private moreIntensePeaksWithin mzWindow (spectrumData:(float*float) []) peakIdx =
        let halfWindowSize = mzWindow / 2. 
        let currentPeakMz,currentPeakIntensity = spectrumData.[peakIdx] 
        let upperBorder = currentPeakMz + halfWindowSize
        let lowerBorder = currentPeakMz - halfWindowSize
        let rightSidePeaks =
            let rec loop currentIdx moreIntensePeaks = 
                if currentIdx = spectrumData.Length then moreIntensePeaks
                else 
                    let tempPeakMz,tempPeakInt = spectrumData.[currentIdx] 
                    if tempPeakMz > upperBorder then moreIntensePeaks
                    elif tempPeakInt > currentPeakIntensity then 
                         loop (currentIdx+1) (moreIntensePeaks+1) 
                    else loop (currentIdx+1) moreIntensePeaks 
            loop peakIdx 0 
        let leftSidePeaks = 
            let rec loop currentIdx moreIntensePeaks = 
                if currentIdx < 0 then moreIntensePeaks
                else 
                    let tempPeakMz,tempPeakInt = spectrumData.[currentIdx] 
                    if tempPeakMz < lowerBorder then moreIntensePeaks
                    elif tempPeakInt > currentPeakIntensity then 
                         loop (currentIdx-1) (moreIntensePeaks+1) 
                    else loop (currentIdx-1) moreIntensePeaks 
            loop peakIdx 0         
        createRatedPeak (currentPeakMz,currentPeakIntensity) (leftSidePeaks+rightSidePeaks) 
    
    /// Applies the function "moreIntensePeaksWithin" to every peak in spectrumData, thus rating every peak by counting how many more abundant peaks
    /// are present within the given mzWindow (mzwindow/2 +/- ). Filters out peaks that are not within the scanLimits and show to have more abundand peaks than 
    /// defined by qMostAbundandpepsMax in their neighbourhood.
    let private ratedSpectrum mzWindow (qMostAbundandPepsMin, qMostAbundandPepsMax) (lowerScanLimit,upperScanLimit) (spectrumData:(float*float) []) = 
        if spectrumData.Length = 0 then [||]
        else
        let filteredAndRanked = 
            spectrumData
            |> Array.mapi   (fun i (mzData,intensityData) -> 
                                moreIntensePeaksWithin mzWindow spectrumData i      
                            )   
            |> Array.filter (fun rPeak -> 
                                let mzData = fst rPeak.Peak  
                                // scanLimits vllt raus
                                mzData >= lowerScanLimit && mzData <= upperScanLimit && rPeak.MoreIntensePeaksWithinWindow < qMostAbundandPepsMax
                            )                        
        filteredAndRanked
          
    /// Converts the fragmentIon ladder to a theoretical spectrum at a given charge state. Filters out all theoretical peaks, that lie
    /// outside the given lower and upperScanLimits.
    let predictOf (lowerScanLimit,upperScanLimit) chargeState (fragments:PeakFamily<TaggedMass.TaggedMass> list)  =
        let predictPeak charge (taggedMass: TaggedMass.TaggedMass) =
            TaggedPeak.TaggedPeak(taggedMass.Iontype, (Mass.toMZ taggedMass.Mass charge), nan)
        let computePeakFamily charge fragments =
            let mainPeak = predictPeak charge fragments.MainPeak
            let dependentPeaks =
                fragments.DependentPeaks
                |> List.fold (fun acc dependent ->
                                 if charge <= 1. then
                                     predictPeak charge dependent
                                     :: acc
                                 else
                                     acc
                              ) []
            createPeakFamily mainPeak dependentPeaks
        // The families come out in reversed fragment order, and with two charge states the
        // singly charged family of a fragment precedes its doubly charged one, which is the order
        // the former list building produced. Filling the array from its end gives it directly.
        let fragments = List.toArray fragments
        let perFragment = if chargeState <= 1. then 1 else 2
        let result = Array.zeroCreate (fragments.Length * perFragment)
        for f = 0 to fragments.Length - 1 do
            let slot = (fragments.Length - 1 - f) * perFragment
            if perFragment = 1 then
                result.[slot] <- computePeakFamily 1. fragments.[f]
            else
                result.[slot] <- computePeakFamily 1. fragments.[f]
                result.[slot + 1] <- computePeakFamily 2. fragments.[f]
        result

    /// Searches for a Peak in ratedSpectrum that matches the target mz (tarmz) within a certain tolerance defined by mzmatchingTolerance.
    let private hasMatchingPeakMZTol (mzMatchingTolerance: float) (tarmz: float) (ratedSpectrum: RatedPeak []) =        
        match binarySearchIndexBy  (fun (rPeak: RatedPeak) ->
                                let mz = fst rPeak.Peak  
                                if ((mz-tarmz) |> abs) <= mzMatchingTolerance then 0 
                                elif tarmz < mz then -1
                                else 1
                            ) ratedSpectrum with
        | x when x >= 0 -> Some ratedSpectrum.[x] 
        | x             ->  None 

    /// Searches for a Peak in ratedSpectrum that matches the target mz (tarmz) within a certain tolerance expressed in parts per million of the tarmz, defined by mzMatchingTolerancePPM.
    let private hasMatchingPeakPPM (mzMatchingTolerancePPM: float) (tarmz: float) (ratedSpectrum: RatedPeak []) = 
        match binarySearchIndexBy (fun (rPeak: RatedPeak) ->
                                let mz = fst rPeak.Peak  
                                if ((mz-tarmz) |> abs) <= Mass.deltaMassByPpm mzMatchingTolerancePPM tarmz then 0 
                                elif tarmz < mz then -1
                                else 1
                            ) ratedSpectrum with
        | x when x >= 0 -> Some ratedSpectrum.[x] 
        | _             -> None 


    [<Struct>]
    type CountedMatches = {
        Q                       : int
        N_BIons                 : int
        N_YIons                 : int
        N_Neutral               : int
        K_BIons                 : int
        K_YIons                 : int
        K_Neutral               : int
        MatchedSum              : float
        }

    let createCountedMatches q n_BIons n_YIons n_Neutral k_BIons k_YIons k_Neutral matchedSum = {
        Q           = q                       
        N_BIons     = n_BIons                 
        N_YIons     = n_YIons                 
        N_Neutral   = n_Neutral               
        K_BIons     = k_BIons                 
        K_YIons     = k_YIons                 
        K_Neutral   = k_Neutral               
        MatchedSum  = matchedSum  
        }

    let findMaxIterator qMostAbundandPepsMin countArrayL (rPeak:RatedPeak) = 
        match (1+rPeak.MoreIntensePeaksWithinWindow)-qMostAbundandPepsMin with 
        | x when x <= 0 -> 0
        | x when x >= countArrayL-1 -> x
        | x -> x      

    let raiseCountsAndHits qMostAbundandPepsMin (countArray: CountedMatches []) (theoPeak:TaggedPeak.TaggedPeak) (rPeak:RatedPeak) =
        let maxIterator = findMaxIterator qMostAbundandPepsMin countArray.Length rPeak
        match theoPeak.Iontype with
        | iontype when iontype.HasFlag(Ions.IonTypeFlag.B) ->
            for i = 0 to countArray.Length-1 do
                let currentMatchingCount = countArray.[i] 
                if i >= maxIterator then
                    countArray.[i] <- createCountedMatches 
                                        currentMatchingCount.Q                       
                                        (currentMatchingCount.N_BIons+1)                 
                                        (currentMatchingCount.N_YIons)                 
                                        (currentMatchingCount.N_Neutral)               
                                        (currentMatchingCount.K_BIons+1)                 
                                        (currentMatchingCount.K_YIons)                 
                                        (currentMatchingCount.K_Neutral)               
                                        (currentMatchingCount.MatchedSum+getIntensity rPeak)
                else 
                    countArray.[i] <- createCountedMatches 
                                        currentMatchingCount.Q                       
                                       (currentMatchingCount.N_BIons+1)                 
                                       (currentMatchingCount.N_YIons)                 
                                       (currentMatchingCount.N_Neutral)               
                                       (currentMatchingCount.K_BIons)                 
                                       (currentMatchingCount.K_YIons)                 
                                       (currentMatchingCount.K_Neutral)               
                                       (currentMatchingCount.MatchedSum) 
        | iontype when iontype.HasFlag(Ions.IonTypeFlag.Y) ->
            for i = 0 to countArray.Length-1 do
                let currentMatchingCount = countArray.[i] 
                if i >= maxIterator then
                    countArray.[i] <- createCountedMatches 
                                        currentMatchingCount.Q                       
                                        (currentMatchingCount.N_BIons)                 
                                        (currentMatchingCount.N_YIons+1)                 
                                        (currentMatchingCount.N_Neutral)               
                                        (currentMatchingCount.K_BIons)                 
                                        (currentMatchingCount.K_YIons+1)                 
                                        (currentMatchingCount.K_Neutral)               
                                        (currentMatchingCount.MatchedSum+getIntensity rPeak)
                else 
                    countArray.[i] <- createCountedMatches 
                                        currentMatchingCount.Q                       
                                       (currentMatchingCount.N_BIons)                 
                                       (currentMatchingCount.N_YIons+1)                 
                                       (currentMatchingCount.N_Neutral)               
                                       (currentMatchingCount.K_BIons)                 
                                       (currentMatchingCount.K_YIons)                 
                                       (currentMatchingCount.K_Neutral)               
                                       (currentMatchingCount.MatchedSum) 
        | iontype when iontype.HasFlag(Ions.IonTypeFlag.Neutral) ->
            for i = 0 to countArray.Length-1 do
                let currentMatchingCount = countArray.[i] 
                if i >= maxIterator then
                    countArray.[i] <- createCountedMatches 
                                        currentMatchingCount.Q                       
                                        (currentMatchingCount.N_BIons)                 
                                        (currentMatchingCount.N_YIons)                 
                                        (currentMatchingCount.N_Neutral+1)               
                                        (currentMatchingCount.K_BIons)                 
                                        (currentMatchingCount.K_YIons)                 
                                        (currentMatchingCount.K_Neutral+1)               
                                        (currentMatchingCount.MatchedSum+getIntensity rPeak)   
                else 
                    countArray.[i] <- createCountedMatches 
                                        currentMatchingCount.Q                       
                                       (currentMatchingCount.N_BIons)                 
                                       (currentMatchingCount.N_YIons)                 
                                       (currentMatchingCount.N_Neutral+1)               
                                       (currentMatchingCount.K_BIons)                 
                                       (currentMatchingCount.K_YIons)                 
                                       (currentMatchingCount.K_Neutral)               
                                       (currentMatchingCount.MatchedSum)    
        | _ -> 
            for i = 0 to countArray.Length-1 do
                let currentMatchingCount = countArray.[i] 
                if i >= maxIterator then
                    countArray.[i] <- createCountedMatches 
                                        currentMatchingCount.Q                       
                                        (currentMatchingCount.N_BIons)                 
                                        (currentMatchingCount.N_YIons)                 
                                        (currentMatchingCount.N_Neutral)               
                                        (currentMatchingCount.K_BIons)                 
                                        (currentMatchingCount.K_YIons)                 
                                        (currentMatchingCount.K_Neutral)               
                                        (currentMatchingCount.MatchedSum)                
                else 
                    countArray.[i] <- createCountedMatches 
                                        currentMatchingCount.Q                       
                                       (currentMatchingCount.N_BIons)                 
                                       (currentMatchingCount.N_YIons)                 
                                       (currentMatchingCount.N_Neutral)               
                                       (currentMatchingCount.K_BIons)                 
                                       (currentMatchingCount.K_YIons)                 
                                       (currentMatchingCount.K_Neutral)               
                                       (currentMatchingCount.MatchedSum) 

    let raiseOnlyCounts (countArray: CountedMatches []) (theoPeak:TaggedPeak.TaggedPeak) =
        match theoPeak.Iontype with
        | iontype when iontype.HasFlag(Ions.IonTypeFlag.B) ->
            for i = 0 to countArray.Length-1 do
                let currentMatchingCount = countArray.[i] 
                countArray.[i] <- createCountedMatches 
                                    currentMatchingCount.Q                       
                                    (currentMatchingCount.N_BIons+1)                 
                                    (currentMatchingCount.N_YIons)                 
                                    (currentMatchingCount.N_Neutral)               
                                    (currentMatchingCount.K_BIons)                 
                                    (currentMatchingCount.K_YIons)                 
                                    (currentMatchingCount.K_Neutral)               
                                    (currentMatchingCount.MatchedSum) 
        | iontype when iontype.HasFlag(Ions.IonTypeFlag.Y) ->
            for i = 0 to countArray.Length-1 do
                let currentMatchingCount = countArray.[i] 
                countArray.[i] <- createCountedMatches 
                                    currentMatchingCount.Q                       
                                    (currentMatchingCount.N_BIons)                 
                                    (currentMatchingCount.N_YIons+1)                 
                                    (currentMatchingCount.N_Neutral)               
                                    (currentMatchingCount.K_BIons)                 
                                    (currentMatchingCount.K_YIons)                 
                                    (currentMatchingCount.K_Neutral)               
                                    (currentMatchingCount.MatchedSum) 
        | iontype when iontype.HasFlag(Ions.IonTypeFlag.Neutral) ->
            for i = 0 to countArray.Length-1 do
                let currentMatchingCount = countArray.[i] 
                countArray.[i] <- createCountedMatches 
                                    currentMatchingCount.Q                       
                                    (currentMatchingCount.N_BIons)                 
                                    (currentMatchingCount.N_YIons)                 
                                    (currentMatchingCount.N_Neutral+1)               
                                    (currentMatchingCount.K_BIons)                 
                                    (currentMatchingCount.K_YIons)                 
                                    (currentMatchingCount.K_Neutral)               
                                    (currentMatchingCount.MatchedSum)   
        | _ -> 
            for i = 0 to countArray.Length-1 do
                let currentMatchingCount = countArray.[i] 
                countArray.[i] <- createCountedMatches 
                                    currentMatchingCount.Q                       
                                    (currentMatchingCount.N_BIons)                 
                                    (currentMatchingCount.N_YIons)                 
                                    (currentMatchingCount.N_Neutral)               
                                    (currentMatchingCount.K_BIons)                 
                                    (currentMatchingCount.K_YIons)                 
                                    (currentMatchingCount.K_Neutral)               
                                    (currentMatchingCount.MatchedSum)  
        
    /// Iterates theoSpect and counts the number of common peaks in ratedspectrum. A peak is considered as common if the mz value lies within a user defined matching tolerance. This procedure
    /// is carried out (qMostAbundandPepsMax-qMostAbundandPepsMin+1) times every time with a different qThreshold to account for peak intensity. A peak is counted as matching if it has
    /// less more abundant peaks in the neighbourhood than the current qThreshold.
    let private countMatches (lowerScanLimit,upperScanLimit) (qMostAbundandPepsMin, qMostAbundandPepsMax) (mzMatchingTolerance: float) (theoSpect: PeakFamily<TaggedPeak.TaggedPeak> []) (ratedSpectrum: RatedPeak []) =    
        let countArray = Array.init (qMostAbundandPepsMax-qMostAbundandPepsMin+1) (fun i -> createCountedMatches (qMostAbundandPepsMin+i) 0 0 0 0 0 0 0.)   
        let rec findMatches (mzMatchingTolerance: float) (theoSpect: PeakFamily<TaggedPeak.TaggedPeak> []) (ratedSpectrum: RatedPeak []) idx =
            if idx = theoSpect.Length then
                countArray
            else
                let currentPeakFam = theoSpect.[idx]
                if currentPeakFam.MainPeak.Mz >= lowerScanLimit && currentPeakFam.MainPeak.Mz <= upperScanLimit then                
                    match hasMatchingPeakPPM mzMatchingTolerance currentPeakFam.MainPeak.Mz ratedSpectrum with 
                    | Some rPeak ->
                        raiseCountsAndHits qMostAbundandPepsMin countArray currentPeakFam.MainPeak rPeak
                        currentPeakFam.DependentPeaks 
                        |> List.iter ( fun dependentAnnotPeak ->
                                        // Peak not within bounds loop further. 
                                        if dependentAnnotPeak.Mz >= lowerScanLimit && dependentAnnotPeak.Mz <= upperScanLimit then 
                                             match hasMatchingPeakPPM mzMatchingTolerance dependentAnnotPeak.Mz ratedSpectrum with 
                                             //Successful match: raise n and k
                                             | Some rPeak -> raiseCountsAndHits qMostAbundandPepsMin countArray dependentAnnotPeak rPeak
                                             | None       -> raiseOnlyCounts countArray dependentAnnotPeak 
                                        )  
                        findMatches mzMatchingTolerance theoSpect ratedSpectrum (idx+1) 
                    | None ->
                        raiseOnlyCounts countArray currentPeakFam.MainPeak
                        findMatches mzMatchingTolerance theoSpect ratedSpectrum (idx+1)
                else 
                    // Peak not within bounds 
                    findMatches mzMatchingTolerance theoSpect ratedSpectrum (idx+1)
        findMatches mzMatchingTolerance theoSpect ratedSpectrum 0                          

    /// gammaLn of the integers 0 .. n, filled with FSharp.Stats gammaLn on demand.
    type private GammaLnTable() =
        let gate = obj()
        let mutable values : float[] = Array.init 4096 (fun i -> FSharp.Stats.SpecialFunctions.Gamma.gammaLn (float i))
        /// Reads through a local reference, so a table that grows on another thread at the same
        /// time cannot shrink under this reader. Growth happens under the lock and only ever
        /// replaces the array with a longer one.
        member _.At (n: int) =
            let current = values
            if n < current.Length then current.[n]
            else
                lock gate (fun () ->
                    if n >= values.Length then
                        let old = values
                        values <- Array.init (max (n + 1) (old.Length * 2)) (fun i -> if i < old.Length then old.[i] else FSharp.Stats.SpecialFunctions.Gamma.gammaLn (float i))
                    values.[n])

    /// gammaLn at the integers 0 .. n, computed once with the same function the score used
    /// directly, so the table holds the same values. Grows when a larger argument shows up.
    let private lnGammaTable : GammaLnTable = GammaLnTable()

    /// Helper function to calculate the AndromedaScore
    let private lnProb (n: int) (k: int) lnp lnq =
        let fn = float n
        let fk = float k
        let s1 = -fk * lnp - (fn - fk) * lnq - lnGammaTable.At(n + 1)
                    + lnGammaTable.At(k + 1) + lnGammaTable.At(n - k + 1)
        s1

    /// Calculates the andromedaLike Score based on offered peaks n and matched peaks k, as well as the q value threshold topx
    let private scoreFuncCompute n k topx =
        let fTopx = float topx
        let p1 = Math.Min(fTopx/100.0,0.5)
        let lnp = Math.Log(p1)
        let lnq = Math.Log(1.0-p1)
        let mutable acc = 0.
        for i = k to n do
            acc <- acc + Math.Exp(-lnProb n i lnp lnq)
        acc <- -Math.Log(acc)
        10. * acc / Math.Log(10.)

    /// The score depends on nothing but the three integers, so every computed value is kept.
    let private scoreFuncMemo = System.Collections.Concurrent.ConcurrentDictionary<struct (int * int * int), float>()

    let private scoreFuncImpl n k topx =
        let key = struct (n, k, topx)
        match scoreFuncMemo.TryGetValue key with
        | true, v -> v
        | _ ->
            let v = scoreFuncCompute n k topx
            scoreFuncMemo.[key] <- v
            v

    /// Correction factor applided to the andromeda score based on the observed precursor Mz.
    /// Values are taken from the original andromeda release.
    let private mzCorrection mass = 
        0.024*(mass - 600.)

    /// Correction factor applided to the andromeda score based on the amount modifications.
    /// Values are taken from the original andromeda release.
    let private modCorrection nMod = 
        match nMod with 
        | 0  -> 42.
        | 1  -> 28.
        | 2  -> 22.
        | 3  -> 16.
        | 4  -> 9.
        | 5  -> 5.
        | 6  -> 2.
        | _  -> 0.
    
    /// Correction factor applided to the andromeda score based on the amount of missCleavages.
    /// Values are taken from the original andromeda release.
    let private cleavageCorrection nCleave consecutiveCleavage =
        match nCleave with 
        | 0 -> 53.2
        | 1 -> 
            if consecutiveCleavage then 
                42.1
            else 
                31.1
        | 2 -> 17. 
        |_  -> 0.

    let optimizeAndromedaScore precursorMZ countedMatches = 
        let massCorr = mzCorrection precursorMZ //(Mass.toMZ lookUpResult.Mass charge)
        let modCorr  = modCorrection 0
        let cleavageCorr = cleavageCorrection 0 false
        let nWithoutNeutral = countedMatches.N_BIons+countedMatches.N_YIons
        let kWithoutNeutral = countedMatches.K_BIons+countedMatches.K_YIons
        let nWithNeutral = countedMatches.N_BIons+countedMatches.N_YIons+countedMatches.N_Neutral
        let kWithNeutral = countedMatches.K_BIons+countedMatches.K_YIons+countedMatches.K_Neutral
        let rawScore1 = scoreFuncImpl nWithoutNeutral kWithoutNeutral countedMatches.Q
        let rawScore2 = scoreFuncImpl nWithNeutral kWithNeutral countedMatches.Q
        if rawScore1 > rawScore2 then 
            let finalScore = 
                let tmp = (rawScore1 + massCorr + modCorr + cleavageCorr - 100.)
                if tmp < 0. then 0. else tmp 
            createMatchingScore finalScore nWithoutNeutral kWithoutNeutral countedMatches.Q
        else 
            let finalScore = 
                let tmp = (rawScore2 + massCorr + modCorr + cleavageCorr - 100.)
                if tmp < 0. then 0. else tmp 
            createMatchingScore finalScore nWithNeutral kWithNeutral countedMatches.Q 
    
    /// Calculates the X!TandemHyperScore
    let calcHyperScore (countedMatches:CountedMatches) =
        let res = log countedMatches.MatchedSum + factorialLn countedMatches.K_BIons + factorialLn countedMatches.K_YIons + factorialLn countedMatches.K_Neutral
        Math.Max(0.,res)

    /// Computes the andromedalike score of a theoretical spectrum vs. the ratedSpectrum. Subsequently applies mzCorrection - modCorrection and cleavageCorrection to the score. 
    let private scoreTheoVsRecordedSpec (lowerScanLimit,upperScanLimit) (qMostAbundandPepsMin, qMostAbundandPepsMax)  matchingTolPPM  precursorMZ (*charge*) theoSpec (ratedSpectrum: RatedPeak []) =     
        let matches = countMatches (lowerScanLimit,upperScanLimit) (qMostAbundandPepsMin, qMostAbundandPepsMax) matchingTolPPM theoSpec ratedSpectrum
        let androScore = 
            matches
            |> Array.map (optimizeAndromedaScore precursorMZ)
            |> Array.maxBy (fun x -> x.Score) 
        let bestXTandem = 
            matches
            |> Array.maxBy (fun x -> x.Q) 
        let xTandemScore = 
            bestXTandem
            |> calcHyperScore
        let matchSum = bestXTandem.MatchedSum
        androScore,xTandemScore,matchSum

    /// Converts the fragment ion ladders to a theoretical Sequestlike spectrum at a given charge state. 
    /// Subsequently, the spectrum is binned to the nearest mz bin (binwidth = 1 Da). Filters out peaks 
    /// that are not within the scanLimits.
    let getTheoSpecs scanlimits chargeState (possiblePeptideInfos:list<LookUpResult<AminoAcids.AminoAcid>*FragmentMasses>) =
        TheoreticalSpectra.getTheoSpecs predictOf scanlimits chargeState possiblePeptideInfos
 
    /// Calculates the AndromedaLike scores for all theoretical spectra.
    let calcAndromedaAndXTandemScore qMinAndMax scanlimits matchingTolPPM (spectrum:PeakArray<_>) scanTime chargeState isolationWindowTargetMz (theoreticalSpectra:list<TheoreticalSpectrum<PeakFamily<TaggedPeak.TaggedPeak>[]>> ) spectrumID =
        
        let fCharge = float chargeState
        //
        let binnedRecSpec = 
            spectrum 
            |> Array.map (fun peak -> peak.Mz,peak.Intensity)
            |> ratedSpectrum 100. qMinAndMax scanlimits
        //
        let ms_mass = Mass.ofMZ isolationWindowTargetMz fCharge
        let androScores,xTandemScores = 
            theoreticalSpectra 
            |> List.fold (fun (andro,xTandem) theoreticalSpectrum -> 
                                let lookUpResult = theoreticalSpectrum.LookUpResult
                                let seqL = lookUpResult.BioSequence.Length
 
                                // Calculates score for target sequence
                                let theoSpec = theoreticalSpectrum.TheoSpec
                                let targetScoreAndro,targetScoreXTandem,matchSum = scoreTheoVsRecordedSpec scanlimits qMinAndMax matchingTolPPM isolationWindowTargetMz theoSpec binnedRecSpec
                                let targetResultAndro = createSearchEngineResult SearchEngineResult.SearchEngine.AndromedaLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod true scanTime lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass seqL targetScoreAndro.Score nan nan
                                let targetResultXTandem = createSearchEngineResult SearchEngineResult.SearchEngine.XTandemLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod true scanTime lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass seqL targetScoreXTandem nan nan
                                
                                let theoSpec_Decoy = theoreticalSpectrum.DecoyTheoSpec
                                let decoyScoreAndro,decoyScoreXTandem,matchSum = scoreTheoVsRecordedSpec scanlimits qMinAndMax matchingTolPPM isolationWindowTargetMz theoSpec_Decoy binnedRecSpec
                                let decoyResultAndro = createSearchEngineResult SearchEngineResult.SearchEngine.AndromedaLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod false scanTime lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass seqL decoyScoreAndro.Score nan nan
                                let decoyResultXTandem = createSearchEngineResult SearchEngineResult.SearchEngine.XTandemLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod false scanTime lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass seqL decoyScoreXTandem nan nan
                                
                                targetResultAndro :: decoyResultAndro :: andro,targetResultXTandem :: decoyResultXTandem :: xTandem      
                     ) ([],[])
        let androScoresFinal =
            calcNormDeltaBestToRest (androScores |> List.sortBy (fun sls -> - sls.Score))       
            |> calcNormDeltaNext
        let xTandemScoresFinal =
            calcNormDeltaBestToRest (xTandemScores |> List.sortBy (fun sls -> - sls.Score))       
            |> calcNormDeltaNext
        androScoresFinal,xTandemScoresFinal 
   
   