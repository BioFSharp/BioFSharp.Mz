namespace BioFSharp.Mz

open FSharp.Stats
open FSharpAux

module FDRControl = 

    module MAYU =

        // FDR estimation using MAYU
        // Code form 'stirlingLogFactorial' to 'estimatePi0HG' translated from percolator 'ProteinFDREstimator.cpp'

        let private stirlingLogFacorial (n: float) =
            log(sqrt(2. * pi  *n)) + n * log(n) - n

        let private exactLogFactorial (n: float) =
            let rec loop i log_fact =
                if i > n then
                    log_fact
                else
                    let new_log_fact = log_fact + (log i)
                    loop (i + 1.) new_log_fact
            loop 2. 0.

        let private logFactorial (n: float) =
            if n < 1000. then
                exactLogFactorial n
            else
                stirlingLogFacorial n

        let private logBinomial (n: float) (k: float) =
            (logFactorial n) - (logFactorial k) - (logFactorial (n - k))

        let private hypergeometric (x: float) (n: float) (w: float) (d: float) =
            //natural logarithm of the probability
            if (d > 0.) then
                exp((logBinomial w x) + (logBinomial (n - w) (d - x)) - (logBinomial n d))
            else 0.

        /// Estimates the false positives given the total number of entries, the number of target hits and the number of decoy hits
        let estimatePi0HG (n: float) (targets: float) (cf: float) =
            let rec loop (fp: float) (logprob: float list) =
                if fp > cf then
                    logprob |> List.rev
                else
                    let tp = targets - fp
                    let w = n - tp
                    let prob = hypergeometric fp n w cf
                    loop (fp + 1.) (prob::logprob)
            let logprob = loop 0. []
            let sum = logprob |> List.sum
            let logprob_Norm =
                logprob
                |> List.map (fun x ->
                    x / sum
                )
            // MAYU rounds here to first decimal
            let expectation_value_FP_PID =
                logprob_Norm
                |> List.foldi (fun i acc x ->
                    acc + x * (float i)
                ) 0.
            if (isNan expectation_value_FP_PID) || (isInf expectation_value_FP_PID) then
                0.
            else
                expectation_value_FP_PID

    /// for given data, creates a logistic regression model and returns a mapping function for this model
    let getLogisticRegressionFunction (x:vector) (y:vector) epsilon = 
        let alpha = 
            match FSharp.Stats.Fitting.LogisticRegression.Univariable.estimateAlpha epsilon x y with 
            | Some a -> a
            | None -> failwith "Could not find an alpha for logistic regression of fdr data"
        let weight = FSharp.Stats.Fitting.LogisticRegression.Univariable.coefficient epsilon alpha x y
        FSharp.Stats.Fitting.LogisticRegression.Univariable.fit weight

    /// returns scores, pep, q
    let binningFunction bandwidth pi0 (scoreF: 'A -> float) (isDecoyF: 'A -> bool) (data:'A[])  = 
        let totalDecoyProportion = 
            let decoyCount = Array.filter isDecoyF data |> Array.length |> float
            let totalCount = data |> Array.length  |> float
            1. / (2. * decoyCount / totalCount)
        data
        |> Array.groupBy (fun s -> floor (scoreF s / bandwidth))
        |> Array.sortBy fst
        |> Array.map (fun (k,values)->
            let median     = values |> Array.map scoreF |> Array.average
            let totalCount = values |> Array.length |> float
            let decoyCount = values |> Array.filter isDecoyF |> Array.length |> float |> (*) totalDecoyProportion
            //(median |> float,(decoyCount * pi0  / totalCount))
            median,totalCount,decoyCount
                //(median, totalCount )
        )
        |> fun a ->
            a
            |> Array.mapi (fun i (median,totalCountBin,decoyCountBin) ->
                            /// TODO: Accumulate totalCount + totalDecoyCount beforeHand and skip the time intensive mapping accross the array in each iteration.
                            let _,totalCountRight,decoyCountRight = a.[i..a.Length-1] |> Array.reduce (fun (x,y,z) (x',y',z') -> x+x',y+y',z+z')
                            (median,(pi0 * 2. * decoyCountBin / totalCountBin),(pi0 * 2. * decoyCountRight / totalCountRight))
                          )
        |> Array.sortBy (fun (score,pep,q) -> score) 
        |> Array.unzip3 
        |> fun (score,pep,q) -> vector score, vector pep, vector q
    
    /// Calculates q value mapping funtion for target/decoy dataset
    let getQValueFunc pi0 (scoreF: 'A -> float) (isDecoyF: 'A -> bool) (data:'A[]) = 
        let bw = 
            data 
            |> Array.map scoreF
            |> FSharp.Stats.Distributions.Bandwidth.nrd0
        let (scores,_,q) = binningFunction bw pi0 scoreF isDecoyF data
        getLogisticRegressionFunction scores q 0.0000001

    /// Calculates q values for target/decoy dataset
    let getQValues pi0 (scoreF: 'A -> float) (isDecoyF: 'A -> bool) (data:'A[]) = 
        let f = getQValueFunc pi0 scoreF isDecoyF data
        Array.map (scoreF >> f) data

    /// Calculates pep value mapping funtion for target/decoy dataset
    let getPEPValueFunc pi0 (scoreF: 'A -> float) (isDecoyF: 'A -> bool) (data:'A[]) = 
        let bw = 
            data 
            |> Array.map scoreF
            |> FSharp.Stats.Distributions.Bandwidth.nrd0
        let (scores,pep,_) = binningFunction bw pi0 scoreF isDecoyF data
        getLogisticRegressionFunction scores pep 0.0000001

    /// Calculates pep values for target/decoy dataset
    let getPEPValues pi0 (scoreF: 'A -> float) (isDecoyF: 'A -> bool) (data:'A[]) = 
        let f = getPEPValueFunc pi0 scoreF isDecoyF data 
        Array.map (scoreF >> f) data 
    
    /// Calculates Q-Values from pep-values
    let getQValuesFromPEPValues (pepValues : float []) = 
        let q : float [] = Array.zeroCreate pepValues.Length
        let pep = Array.sort pepValues
        let rec loop sum i =
            if i = pepValues.Length then 
                q
            else 
                let newSum  = pep.[i] + sum
                q.[i] <- newSum / (float (i+1))
                loop newSum (i+1)
        loop 0. 0