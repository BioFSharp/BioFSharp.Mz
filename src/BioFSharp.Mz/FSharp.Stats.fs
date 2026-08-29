namespace BioFSharp.Mz

open FSharp.Stats
open FSharp.Stats.Signal
open System


module NonLinearRegression' = 
    ///
    let standardErrorOfPrediction dOF (predicted:float []) (actual:float [])  =
        let n = actual.Length-1 |> float 
        match n with
        | x when x > dOF -> 
            let sumOfResSq = Array.fold2 (fun acc yReal yPred  -> acc + ((yPred-yReal)**2.) ) 0.0  actual predicted
            sqrt( (sumOfResSq / (n - dOF)))
        | _             -> nan

module PeakDetection = 

    module Wavelet = 

        type WaveletPeak = {
            XAxisIndex  : int 
            ScaleIndex  : int
            Scale       : float
            Correlation : float
            }

        type Gaussian = {
            Amplitude   : float
            XLoc        : float
            Stdev       : float
            Trace       : (float*float) []
            Function    : float -> float
            }

        type PeakGroup = {   
            Start       : float
            End         : float
            Data        : (float*float)[]
            Fits        : Gaussian list
            SumTrace    : float []
            Convolved   : bool
            }

        type Parameters = {
            Borderpadding           : int 
            BorderPadMethod         : Padding.BorderPaddingMethod 
            InternalPaddingMethod   : Padding.InternalPaddingMethod 
            HugeGapPaddingMethod    : Padding.HugeGapPaddingMethod 
            HugeGapPaddingDistance  : float
            MaxPeakLength           : float 
            NoiseQuantile           : float 
            MinSNR                  : float
            }

        ///
        let private getColLowBound  (xData: float[]) centerIdx mzTol = 
            let rec loop  (xData: float []) center i mzTol =
                if i <= 0 then 0
                else                                              
                    match (xData.[center] - xData.[i]) with 
                    | x when x >= mzTol -> i+1   
                    | x when x <= mzTol -> loop xData center (i-1) mzTol
            loop xData centerIdx (centerIdx-1) mzTol

        /// 
        let private getColHighBound (xData: float []) centerIdx mzTol = 
            let rec loop (xData: float []) center i mzTol = 
                if i >= xData.Length-1 then xData.Length-1
                else                                             
                    match xData.[i] - xData.[center] with 
                    | x when x >= mzTol -> i-1
                    | x when x <= mzTol -> loop xData center (i+1) mzTol
            loop xData centerIdx (centerIdx+1) mzTol

        ///
        let private getScales maxScale averageDistance =
            let tmp =
                [|averageDistance*2. .. Ops.roundTo 2 (averageDistance/2.) .. maxScale|]
            tmp

        ///
        let private getMinPeakDistance averageDistance = averageDistance*3.

        ///
        let private aic n nParams likelihood = 
            (float n * Math.Log(likelihood / float n)) + 2. * (float nParams) 

        ///
        let private optimizeWaveletFits (origData:(float*float) []) (groupedPeaks:(WaveletPeak*Gaussian)list) = 
            match groupedPeaks with 
            | (wp,g)::[] ->
                let min = g.XLoc - 2. * g.Stdev
                let max = g.XLoc + 2. * g.Stdev
                let data = origData |> Array.filter (fun (x,y) -> x >= min  && x <= max )
                let yHat = origData |> Array.map (fst >> g.Function)
                {Start=min;End=max;Data=data;Fits=[g];SumTrace=yHat;Convolved=false}
            | _ -> 
                let min  = groupedPeaks |> List.map (fun (wp,g) -> g.XLoc - 2. * g.Stdev) |> List.min
                let max  = groupedPeaks |> List.map (fun (wp,g) -> g.XLoc + 2. * g.Stdev) |> List.max
                let xy = origData |> Array.filter (fun (x,y) -> x >= min && x <= max) 
                let x,y = xy |> Array.unzip        
                let bestGroup = 
                    groupedPeaks
                    |> List.map (fun (wp,g) -> 
                        let yHat = x |> Array.map g.Function
                        g,yHat
                    )
                    |> FSharpAux.List.powerSetOf
                    |> List.minBy (fun groupedPeaks ->    
                        let m = Matrix.ofJaggedArray (groupedPeaks |> List.map snd |> Array.ofList)
                        let yHat'= m |> Matrix.mapiCols (fun i v -> v |> Vector.sum)
                        let rss = Seq.map2 (fun x y -> (x-y)**2. ) y yHat' |> Seq.sum
                        aic m.NumCols m.NumRows rss
                    )
                    |> List.map fst
                let sumTrace = 
                    let m = Matrix.ofJaggedArray (bestGroup |> List.map (fun g ->  origData |> Array.map (fst >> g.Function)) |> Array.ofList )
                    let yHat' = m |> Matrix.mapiCols (fun i v -> v |> Vector.sum)
                    yHat'
                    |> Array.ofSeq            
                {Start=min;End=max;Data=xy;Fits=bestGroup;SumTrace=sumTrace;Convolved=true}
             
        ///        
        let identifyPeaksBy (borderpadding:int) (borderPadMethod:Padding.BorderPaddingMethod) (internalPaddingMethod:Padding.InternalPaddingMethod) (hugeGapPaddingMethod:Padding.HugeGapPaddingMethod) 
            (maxDistance:float) maxPeakLength noiseQuantile minSNR (trace:(float*float)[]) =
            let maxScale = maxPeakLength / 6.
            let averageDistance = Padding.HelperFunctions.getMedianSpacing trace (-)
            let minDistancePadding = averageDistance / 2. 
            let paddedData = Padding.pad trace minDistancePadding maxDistance (-) (+) borderpadding borderPadMethod internalPaddingMethod hugeGapPaddingMethod
            let minPeakDistance = getMinPeakDistance minDistancePadding
            let scales = getScales maxScale minDistancePadding
            let transformations = 
                scales
                |> Array.map Wavelet.createRicker
                |> Array.map (fun x -> ContinuousWavelet.transform paddedData (-) borderpadding x |> Array.unzip)
            let xVals = fst transformations.[0]
            let corrMatrix = 
                transformations
                |> Array.map snd
                |> Matrix.ofJaggedArray
            let noiseLevel = (corrMatrix.Row 0).ToArray() |> FSharp.Stats.Quantile.compute noiseQuantile
            let peakMatrix = 
                corrMatrix
                |> Matrix.mapi (fun m n x -> 
                    if n < 2 || n >= corrMatrix.NumCols-3 then 0. 
                    elif x > corrMatrix.[m,n-1] && x > corrMatrix.[m,n-2] && x > corrMatrix.[m,n+1] && x > corrMatrix.[m,n+2] && x >= (minSNR*noiseLevel) then 
                        x
                    else 0.
                )
            let maxCorrScale = 
                peakMatrix 
                |> Matrix.mapiCols (fun i x ->
                    let mutable maxScaleIdx,maxScale,maxCorr = 0,0.,0.
                    for scale = 1 to x.Length-1 do 
                        if x.[scale] > maxCorr then 
                            maxScaleIdx <- scale 
                            maxScale    <- scales.[scale]
                            maxCorr     <- x.[scale]
                    {XAxisIndex=i ;ScaleIndex=maxScaleIdx; Scale=maxScale; Correlation=maxCorr}   
                )
                |> Array.ofSeq
            let mergedPeaks = 
                let finalPeaks = ResizeArray<WaveletPeak>(100)
                for i = 2 to maxCorrScale.Length-3 do    
                    let currentPeak = maxCorrScale.[i]
                    let currentCorr = currentPeak.Correlation
                    if currentCorr < maxCorrScale.[i-1].Correlation ||
                        currentCorr < maxCorrScale.[i-2].Correlation ||
                        currentCorr < maxCorrScale.[i+1].Correlation ||
                        currentCorr < maxCorrScale.[i+2].Correlation ||
                        currentCorr = 0. then ()
                    else 
                        if currentPeak.ScaleIndex > 0 (*nouse*) then 
                            let lowBound  = getColLowBound xVals i (minPeakDistance/2.)
                            let highBound = getColHighBound xVals i (minPeakDistance/2.)
                            let mutable maxCorr = 0.0
                            let mutable maxCol = 0
                            for j = lowBound to highBound do 
                                let scale = maxCorrScale.[j].ScaleIndex
                                if corrMatrix.[scale,j] > maxCorr && scale > 0  then
                                    maxCorr <- corrMatrix.[scale,j]
                                    maxCol <- j 
                            let refinedPeak = maxCorrScale.[maxCol] 
                            if finalPeaks.Count = 0 then 
                                finalPeaks.Add refinedPeak
                            else 
                                let prevPeak = finalPeaks.[finalPeaks.Count-1]
                                let mzDiff = xVals.[refinedPeak.XAxisIndex] - xVals.[prevPeak.XAxisIndex]
                                if mzDiff > minPeakDistance then
                                    finalPeaks.Add refinedPeak
                                elif refinedPeak.Correlation > prevPeak.Correlation then
                                    finalPeaks.[finalPeaks.Count-1] <- refinedPeak
                finalPeaks.ToArray()
            let refinedPeaks = 
                mergedPeaks 
                |> Array.map (fun wp -> 
                    let loc,apex = paddedData |> Array.minBy (fun x -> abs(fst x - (xVals.[wp.XAxisIndex])))
                    let f = FSharp.Stats.Fitting.NonLinearRegression.Table.gaussModel.GetFunctionValue ([|apex;loc;wp.Scale|] |> vector)
                    let peakSim = 
                        trace 
                        |> Array.map fst 
                        |> Array.map (fun x -> x,f x)
                    let gaussian = 
                        {
                            Amplitude   = apex
                            XLoc        = loc
                            Stdev       = wp.Scale
                            Trace       = peakSim
                            Function = f
                        }
                    wp,gaussian 
                )
                |> Array.sortBy (fun (wp,g) -> wp.XAxisIndex)
                |> List.ofArray
                |> List.filter (fun (wp,g) -> g.Stdev > 0. && g.Amplitude > 0.)
            let groupedPeaks = 
                let isOverlapping (currentPeak:(WaveletPeak*Gaussian)) (candidate:(WaveletPeak*Gaussian)) = 
                    (((snd candidate).XLoc < ((snd currentPeak).XLoc + ((snd currentPeak).Stdev * 2.))) && ((snd candidate).XLoc > ((snd currentPeak).XLoc - ((snd currentPeak).Stdev * 2.))) ) ||
                    ( ((snd currentPeak).XLoc < ((snd candidate).XLoc + ((snd candidate).Stdev * 2.))) && ((snd currentPeak).XLoc > ((snd candidate).XLoc - ((snd candidate).Stdev * 2.)))     ) 
                let overlappingPeaks currentFamily peaksLeft = 
                    let rec findOverlappingPeaks  (overlappingPeaks:(WaveletPeak*Gaussian) list) (nonOverlappingPeaks:(WaveletPeak*Gaussian) list) (currentPeak:(WaveletPeak*Gaussian)) (peaksLeft:(WaveletPeak*Gaussian) list) =
                        match peaksLeft with 
                        | [] -> overlappingPeaks,nonOverlappingPeaks
                        | h::t ->
                            if isOverlapping currentPeak h then 
                                findOverlappingPeaks (h::overlappingPeaks) nonOverlappingPeaks currentPeak t 
                            else 
                                findOverlappingPeaks overlappingPeaks (h::nonOverlappingPeaks) currentPeak t 
                              
                    let rec loop newPeakFamily currentFamily peaksLeft =
                        match currentFamily with 
                        | [] -> newPeakFamily,peaksLeft
                        | cp::t -> 
                            let overlappingPeaks,nonOverlappingPeaks = findOverlappingPeaks [] [] cp peaksLeft
                            loop (overlappingPeaks@newPeakFamily) t nonOverlappingPeaks
                    loop [] currentFamily peaksLeft             
                let acc = ResizeArray<(WaveletPeak*Gaussian)list>()
                let rec loop (acc:ResizeArray<(WaveletPeak*Gaussian)list>) currentFam peaksLeft =
                    match peaksLeft with 
                    | [] -> 
                        currentFam 
                        |> acc.Add 
                        acc |> Seq.toList
                    | cP::t ->
                        match currentFam with 
                        | [] -> 
                            loop acc (cP::currentFam) (t)                     
                        | cf  ->
                            let overlapping = overlappingPeaks cf (cP::t) 
                            match overlapping with 
                            | [],rest   -> 
                                cf 
                                |> acc.Add 
                                loop acc [] rest
                            | newPeaks,rest -> 
                                let cf' = newPeaks@cf
                                loop acc cf' rest 
                loop acc [] refinedPeaks                     
            let origData = 
                xVals 
                |> Array.map (fun x -> 
                    paddedData |> Array.minBy (fun (xx,yy) -> abs(xx-x) )
                )

            let fits = 
                groupedPeaks
                |> List.filter  (List.isEmpty >> not)
                |> List.map (optimizeWaveletFits origData)
            let final =
                let m =
                    fits
                    |> List.map (fun x -> x.SumTrace)
                    |> Array.ofList
                // JaggedArray.transpose requires a non-empty outer array; no fits means no peaks
                if Array.isEmpty m then [] else
                m
                |> JaggedArray.transpose
                |> Array.map (fun v -> 
                    let max = v |> Seq.max
                    v 
                    |> Array.map (fun x -> if x < max then 0. else x)
                )
                |> JaggedArray.transpose
                |> Array.mapi (fun i trace -> 
                    let g = fits.[i]
                    let y' = trace
                    let max = y' |> Array.max
                    let sumTrace' = 
                        Array.map2 (fun (x,y) y' -> x,y,y') origData y'
                        |> Array.choose (fun (x,y,y') -> if y' > 0. && y' > max * 0.01 then Some (x,y) else None)           
                    if Array.isEmpty sumTrace' then 
                        None
                    else
                        let start' =  sumTrace' |> Array.minBy fst |> fst
                        let end'   =  sumTrace' |> Array.maxBy fst |> fst
                        Some {g with Start=start';End=end';Data=sumTrace' }
                )
                |> Array.choose id
                |> Array.map (fun pk -> 
                                let x,y = pk.Data |> Array.unzip
                                let traces = 
                                    Array.map2 (fun x y -> 
                                        let yHats = 
                                            pk.Fits 
                                            |> List.mapi (fun i g -> i,g.Function x)
                                        let maxIndex,value = List.maxBy snd yHats
                                        yHats 
                                        |> List.map (fun (i,value) -> if i = maxIndex then x,y else x,0.)
                                        |> Array.ofList
                                    ) x y
                                    |> JaggedArray.transpose
                            
                                let refined = 
                                    pk.Fits
                                    |> List.mapi (fun i g -> 
                                        let trace = 
                                            let tmp = 
                                                traces.[i] |> Array.filter (fun (x,y) -> y > 0.)
                                            if Array.isEmpty tmp then 
                                                None
                                            else    
                                                let maxX,maxY = tmp |> Array.maxBy snd 
                                                tmp 
                                                |> Array.filter (fun x -> snd x > maxY * 0.01)
                                                |> Some 
                                        if trace.IsNone then None 
                                        else
                                        Some {g with Trace = trace.Value}
                                    )
                                    |> List.choose id
                                {pk with Fits = refined}
                            )
                |> List.ofSeq
            final
           
        ///        
        let identifyPeaks (parameters:Parameters) (trace:(float*float)[]) =
            identifyPeaksBy parameters.Borderpadding parameters.BorderPadMethod parameters.InternalPaddingMethod parameters.HugeGapPaddingMethod
                parameters.HugeGapPaddingDistance parameters.MaxPeakLength parameters.NoiseQuantile parameters.MinSNR trace

        //
        let substractBaseLine (yData:float []) =
            let baseLine = FSharp.Stats.Signal.Baseline.baselineAls 10 6 0.05 yData |> Array.ofSeq
            Array.map2 (fun y b ->
                           let c = y - b
                           if c < 0. then 0. else c
                       ) yData baseLine
