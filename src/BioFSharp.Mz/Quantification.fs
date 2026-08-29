namespace BioFSharp.Mz

module Quantification =

    open System
    open SignalDetection
    open FSharp.Stats
    open Fitting.NonLinearRegression
    open Signal.PeakDetection


    module Integration = 

        /// Returns the estimated area beneath the data using the trapezoidal rule. Requires uniform spaced datapoints. 
        let trapezEstAreaOfUniform (xData: float[]) (yData: float[]) =
            if xData.Length <> yData.Length then
                failwith "Both input collections must be of the same size"
            // lower border of area estimation
            let a = xData.[0]
            // upper border of area estimation
            let b = xData.[xData.Length-1]
            let mutable trapezSum = 0.0
            for i = 1 to xData.Length-2 do
                trapezSum <- trapezSum + 2. * yData.[i]
            trapezSum <- trapezSum + yData.[0]
            trapezSum <- trapezSum + yData.[xData.Length-1]
            let area = ( (b-a) / (2. * float xData.Length) ) * trapezSum
            area    
         
        /// Returns the estimated area beneath the data using the trapezoidal rule.
        let trapezEstAreaOf (xData: float[]) (yData: float[]) =
            if xData.Length <> yData.Length then
                failwith "Both input collections must be of the same size"
            let mutable trapezSum = 0.0
            for i = 0 to xData.Length-2 do
                trapezSum <- trapezSum + (xData.[i+1] - xData.[i]) * (yData.[i+1] + yData.[i])
            let area = 0.5 * trapezSum
            area    

    module ParameterEstimation =

        ///
        type CaruanaEstimates = {
            Amplitude: float
            MeanX    : float
            STD      : float
            FWHM     : float
            }

        /// 
        let createCaruanaEstimates amplitude meanX std fwhm = {
            Amplitude=amplitude; MeanX=meanX; STD=std; FWHM=fwhm }

        /// Returns the value of the standard deviation
        let toSTD fwhm = fwhm / ( 2. * sqrt(2. * log 2.) )

        /// Returns the full width at half maximum  
        let toFWHM std = 2. * sqrt(2. * log 2.) * std

        /// Returns the yValue of a gauss function at a given position x.
        let gaussFunc amplitude meanX std x = 
            amplitude * exp((-1.)*((((x-meanX)**2.)/(2.*std**2.))))
  
        /// Estimates the Parameters of a Gaussian function
        /// Warning: This method is sensitive to noisy data. If the noise level of the input parameters is high, smoothing of 
        /// the data is strongly recommended. 
        let caruanaAlgorithm (mzData:float []) (intensityData:float []) =
            let mzData,intensityData = 
                Array.zip mzData intensityData
                |> Array.filter (fun (_,intensity) -> intensity <> 0.)
                |> Array.unzip
                |> fun (x,y) -> Vector.ofArray x, y
            if mzData.Length < 3 || intensityData.Length < 3 then Option.None 
            else 
            let logTransIntensityData = 
                intensityData
                |> Array.map (fun x -> 
                                if x <= 0. then log 1. 
                                else log x
                             )
                |> Vector.ofArray                          
            let polCoeff = FSharp.Stats.Fitting.LinearRegression.OLS.Polynomial.fit 2 mzData logTransIntensityData 
            // f(x) = a1 + a2 * x + a3 * x**2
            let a = polCoeff.[0]
            let b = polCoeff.[1]
            let c = polCoeff.[2]
            let amplitude = exp(a-((b**2.)/(4.*c)))
            let meanX = -b/(2.*c) 
            let fwhm  = sqrt(-1./(2.*c)) * sqrt(2.*log(2.))*2.
            let std   = toSTD fwhm
            Some (createCaruanaEstimates amplitude meanX std fwhm)


        ///
        type EstimatedMoments = {
            ModeY : float 
            MeanX : float 
            Var   : float 
            Std   : float 
            Skew  : float 
            Tau   : float
            }

        ///
        let createEstimatedMoments modeY meanX var std skew tau = {
            ModeY = modeY 
            MeanX = meanX 
            Var   = var   
            Std   = std   
            Skew  = skew  
            Tau   = tau      
            }

        ///
        let meanOfGaussian xData yData = 
            Array.fold2 (fun (sum,acc) x y -> sum+y, acc + (y*x) ) (0.,0.)  xData yData 
            |> fun (sum,acc) -> acc / sum

        ///
        let varianceBy mean xData yData =     
            Array.fold2 (fun (sum,acc) x y -> sum+y, acc + (y * (x - mean)**2.)  ) (0.,0.)  xData yData 
            |> fun (sum,acc) -> acc/sum

        ///
        let varianceOf xData yData =     
            let mean = meanOfGaussian xData yData
            varianceBy mean xData yData 

        ///
        let skewBy mean var xData yData = 
            let m3   =  
                Array.fold2 (fun (sum,acc) x y -> sum+y, acc + (y * (x - mean)**3.) ) (0.,0.)  xData yData 
                |> fun (sum,acc) -> acc/sum
            m3 / (var**1.5)

        ///
        let skewOf xData yData = 
            let mean = meanOfGaussian xData yData
            let var  = varianceBy mean xData yData
            skewBy mean var xData yData 
    
        ///    
        let estTau stdev skew  =
            stdev * (skew/2.) ** 0.333

        ///
        let methodOfMoments (p: IdentifiedPeak) =
            try
            let modeY   = p.Apex.YVal
            let meanX   = meanOfGaussian p.XData p.YData
            let var     = varianceOf p.XData p.YData
            let std     = sqrt(var)
            let skew    = skewBy meanX var p.XData p.YData
            let tau     = estTau std skew
            Some (createEstimatedMoments modeY meanX var std skew tau)
            with 
            | _ -> Option.None
        ///
        let estimatePeakIntegrity (p: IdentifiedPeak) = 
            let zipped = Array.zip p.XData p.YData 
            let left  = zipped |> Array.exists (fun (x,y) -> x < p.Apex.XVal && y > 0.2 * p.Apex.YVal) 
            let right = zipped |> Array.exists (fun (x,y) -> x > p.Apex.XVal && y > 0.2 * p.Apex.YVal)
            left && right
            
        ///
        let estimateMoments (p: IdentifiedPeak ) =
            if estimatePeakIntegrity p then 
                methodOfMoments p 
            else 
                let xData,yData = 
                    Array.zip p.XData p.YData 
                    |> Array.filter (fun (x,y) -> y > 0.3 * p.Apex.YVal)
                    |> Array.unzip
                match caruanaAlgorithm xData yData with 
                | Some caruanaEst -> 
                    let skew = skewBy caruanaEst.MeanX  (caruanaEst.STD**2.) p.XData p.YData
                    let tau  = estTau caruanaEst.STD skew
                    Some (createEstimatedMoments p.Apex.YVal caruanaEst.MeanX (caruanaEst.STD**2.) caruanaEst.STD skew tau)
                | Option.None            -> 
                    methodOfMoments p 
      
 
    module HULQ = 

        open Integration
        open ParameterEstimation
        open LevenbergMarquardtConstrained

        type PeakModel = 
            | Gaussian of Model
            | EMG of Model 

        type FittedPeak = {
            Model                       : PeakModel 
            EstimatedParams             : float []
            StandardErrorOfPrediction   : float
            YPredicted                  : float []
            Area                        : float
            }

        ///
        let createFittedPeak model estimatedParams standardErrorOfPrediction yPredicted area = {
            Model                       = model                       
            EstimatedParams             = estimatedParams            
            StandardErrorOfPrediction   = standardErrorOfPrediction  
            YPredicted                  = yPredicted                 
            Area                        = area
            }

        type QuantifiedPeak = {
            Model                       : PeakModel option 
            YPredicted                  : float []
            EstimatedParams             : float []
            StandardErrorOfPrediction   : float
            Area                        : float
            MeasuredApexIntensity       : float
            }

        ///
        let createQuantifiedPeak model yPredicted estimatedParams standardErrorOfPrediction area measuredApexIntensity = {
            Model                       = model                       
            YPredicted                  = yPredicted                  
            EstimatedParams             = estimatedParams             
            StandardErrorOfPrediction   = standardErrorOfPrediction   
            Area                        = area      
            MeasuredApexIntensity       = measuredApexIntensity
            }

        /// Return Option
        let getPeakBy (peaks:IdentifiedPeak []) x =
            let isPartOfPeak x (p:IdentifiedPeak) = p.LeftEnd.XVal < x && p.RightEnd.XVal > x 
            match Array.tryFind (isPartOfPeak x) peaks with
            | Some p        ->  p
            | Option.None   -> Array.minBy (fun p -> abs(p.Apex.XVal - x)) peaks 

        ///
        let integralOfGaussian (coeffs:vector ) =
            sqrt(2.) * coeffs.[0] * coeffs.[2] * sqrt(Math.PI)

        ///
        let integralOfEMGBy meanX sigma tau intensityAtx xPos =    
            let fst = intensityAtx * (sqrt(Math.PI)*tau)
            let intg =         
                let firstT = sqrt(Math.PI) / 2.
                let errorF = 
                    ((xPos-meanX)/(sqrt(2.)*sigma)) - (sigma/(sqrt(2.)*tau))
                    |> FSharp.Stats.SpecialFunctions.Errorfunction.Erf
                firstT * (1.+errorF)
            let exp = exp ( 0.5*((sigma/tau)**2.) - ((xPos-meanX)/tau) )
            fst / (intg * exp)

        ///
        let calcArea model (estParams:vector) =
            match model with 
            | Gaussian f -> 
                integralOfGaussian estParams
            | EMG      f -> 
                let intensityAtx = f.GetFunctionValue estParams estParams.[1]
                integralOfEMGBy estParams.[1] estParams.[2] estParams.[3] intensityAtx estParams.[1]

        ///
        let tryFit model solverOptions lowerBound upperBound xData yData =
            let modelF = 
                match model with 
                | Gaussian f -> f
                | EMG f      -> f       
            try
                let estParams = estimatedParamsVerbose modelF solverOptions 0.001 10.0 lowerBound upperBound xData yData |> Seq.last                                  
                let y' = Array.map (fun xValue -> modelF.GetFunctionValue estParams xValue) xData
                let sEoE_Gauss = NonLinearRegression'.standardErrorOfPrediction (float solverOptions.InitialParamGuess.Length) y' yData
                let area = calcArea model estParams
                if nan.Equals area || infinity.Equals area || infinity.Equals (area * -1.) || nan.Equals sEoE_Gauss then 
                    Option.None 
                else
                    Option.Some (createFittedPeak model (estParams.ToArray()) sEoE_Gauss y' area)
            with 
            | ex -> Option.None      

        ///
        let tryFitGaussian initAmp initMeanX initStdev xData yData =
            let deltaAmp = initAmp * 0.3
            let deltaX = Math.Min((Array.max xData - Array.min xData) / 2.,10.) 
            let deltaStabw = initStdev * 0.5
            let lowerBound = [|Math.Max(initAmp-deltaAmp,0.); Math.Max(initMeanX-deltaX,0.); Math.Max(initStdev-deltaStabw,0.)|] |> vector
            let upperBound = [|initAmp+deltaAmp; initMeanX+deltaX; initStdev+deltaStabw|] |> vector
            let solverOptions = createSolverOption 0.00001 0.001 1000 [|initAmp;initMeanX;initStdev|]
            tryFit (Gaussian Table.gaussModel) solverOptions lowerBound upperBound xData yData

        ///
        let tryFitEMG initAmp initMeanX initStdev initTau xData yData =
            let deltaAmp = initAmp * 0.3
            let deltaX = Math.Min((Array.max xData - Array.min xData) / 2.,10.) 
            let deltaStabw = initStdev * 0.5
            let deltaTau = initTau
            let lowerBound = [|Math.Max(initAmp-deltaAmp,0.); Math.Max(initMeanX-deltaX,0.); Math.Max(initStdev-deltaStabw,0.);Math.Max(initTau-deltaTau,0.00001)|] |> vector
            let upperBound = [|initAmp+deltaAmp; initMeanX+deltaX; initStdev+deltaStabw; initTau+deltaTau|] |> vector
            let solverOptions = createSolverOption 0.00001 0.001 1000 [|initAmp;initMeanX;initStdev;initTau|]
            tryFit (EMG Table.emgModel) solverOptions lowerBound upperBound xData yData

        ///
        let selectModel (models: FittedPeak []) = 
            Array.minBy (fun (x:FittedPeak) -> x.StandardErrorOfPrediction) models


        ///
        let quantifyPeak (p: IdentifiedPeak ) =
            match estimateMoments p with 
            | Some moments -> 
                let gaussian = tryFitGaussian moments.ModeY moments.MeanX moments.Std p.XData p.YData
                let emg      = tryFitEMG moments.ModeY moments.MeanX moments.Std moments.Tau p.XData p.YData
                match gaussian, emg with 
                | Some g, Some emg -> 
                    let finalFit   = selectModel [|g;emg|]
                    createQuantifiedPeak  (Some finalFit.Model) finalFit.YPredicted finalFit.EstimatedParams finalFit.StandardErrorOfPrediction finalFit.Area p.Apex.YVal
                | Some finalFit, Option.None  | Option.None, Some finalFit -> 
                    createQuantifiedPeak  (Some finalFit.Model) finalFit.YPredicted finalFit.EstimatedParams finalFit.StandardErrorOfPrediction finalFit.Area  p.Apex.YVal
                | _                       ->
                    let area = trapezEstAreaOf p.XData p.YData
                    createQuantifiedPeak Option.None [||] [||] nan area p.Apex.YVal
            | Option.None   -> 
                let area = trapezEstAreaOf p.XData p.YData
                createQuantifiedPeak Option.None [||] [||] nan area p.Apex.YVal              

          
