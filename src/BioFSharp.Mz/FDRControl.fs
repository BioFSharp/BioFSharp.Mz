namespace BioFSharp.Mz

open FSharp.Stats
open FSharpAux
open FSharp.Stats.Fitting
open FSharp.Stats.Fitting.NonLinearRegression

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
            let logprob = 
                [0. .. cf]
                |> List.fold (fun acc fp ->
                    let tp = targets - fp
                    let w = n - tp
                    let prob = hypergeometric fp n w cf
                    prob::acc
                ) []
                |> List.rev
            let sum = logprob |> List.sum
            let logprob_Norm =
                logprob
                |> List.map (fun x -> x / sum)
            // MAYU rounds here to first decimal
            let expectation_value_FP_PID =
                logprob_Norm
                |> List.foldi (fun i acc x -> acc + x * (float i)) 0.
            if (isNan expectation_value_FP_PID) || (isInf expectation_value_FP_PID) then
                0.
            else
                expectation_value_FP_PID

    type ScoreTargetDecoyCount =
        {
            Score      : float
            DecoyCount : float
            TargetCount: float
        }
     
    /// Gives the decoy and target count at a specific score
    let createScoreTargetDecoyCount score decoyCount targetCount =
        {
            Score       = score
            DecoyCount  = decoyCount
            TargetCount = targetCount
        }
    
    /// returns scores, pep, q
    let binningFunction bandwidth pi0 (scoreF: 'A -> float) (isDecoyF: 'A -> bool) (data:'A[])  = 
        let totalDecoyProportion = 
            let decoyCount = 
                Array.filter isDecoyF data
                |> Array.length
                |> float
            let totalCount =
                data
                |> Array.length
                |> float
            1. / (2. * decoyCount / totalCount)
        data
        |> Array.groupBy (fun s -> floor (scoreF s / bandwidth))
        |> Array.sortBy fst
        |> Array.map (fun (k,values)->
            let median     = values |> Array.map scoreF |> Array.average
            let totalCount = values |> Array.length |> float
            let decoyCount = values |> Array.filter isDecoyF |> Array.length |> float |> (*) totalDecoyProportion
            median,totalCount,decoyCount
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
    


    /// Input for QValue calulation
    type QValueInput =
        {
            Score    : float
            IsDecoy  : bool
        }

    let createQValueInput score isDecoy =
        {
            Score     = score
            IsDecoy   = isDecoy
        }

    /// Gives a function to calculate the q value for a score in a dataset using Lukas method and Levenberg Marguardt fitting
    let calculateQValueLogReg fdrEstimate bandwidth (data: 'a []) (isDecoy: 'a -> bool) (decoyScoreF: 'a -> float) (targetScoreF: 'a -> float) =
        // Input for q value calculation
        let createTargetDecoyInput =
            data
            |> Array.map (fun item ->
                if isDecoy item then
                    createQValueInput (decoyScoreF item) true
                else
                    createQValueInput (targetScoreF item) false
            )

        let scores,pep,qVal =
            binningFunction bandwidth fdrEstimate (fun (x: QValueInput) -> x.Score) (fun (x: QValueInput) -> x.IsDecoy) createTargetDecoyInput
            |> fun (scores,pep,qVal) -> scores.ToArray(), pep.ToArray(), qVal.ToArray()

        //Chart.Point (scores,qVal)
        //|> Chart.Show

        // gives a range of 1 to 30 for the steepness. This can be adjusted depending on the data, but normally it should lie in this range
        let initialGuess =
            Fitting.NonLinearRegression.LevenbergMarquardtConstrained.initialParamsOverRange scores qVal [|1. .. 30.|]
            |> Array.map (fun guess -> Table.lineSolverOptions guess)

        // performs Levenberg Marguardt Constrained algorithm on the data for every given initial estimate with different steepnesses and selects the one with the lowest RSS
        let estimate =
            initialGuess
            |> Array.map (fun initial ->
                if initial.InitialParamGuess.Length > 3 then failwith "Invalid initial param guess for Logistic Function"
                let lowerBound =
                    initial.InitialParamGuess
                    |> Array.map (fun param -> param - (abs param) * 0.1)
                    |> vector
                let upperBound =
                    initial.InitialParamGuess
                    |> Array.map (fun param -> param + (abs param) * 0.1)
                    |> vector
                LevenbergMarquardtConstrained.estimatedParamsWithRSS Table.LogisticFunctionDescending initial 0.001 10.0 lowerBound upperBound scores qVal
            )
            |> Array.filter (fun (param,rss) -> not (param |> Vector.exists System.Double.IsNaN))
            |> Array.minBy snd
            |> fst

        let logisticFunction = Table.LogisticFunctionDescending.GetFunctionValue estimate
        logisticFunction



    /// Gives a function to calculate the q value for a score in a dataset using Storeys method
    let calculateQValueStorey (data: 'a[]) (isDecoy: 'a -> bool) (decoyScoreF: 'a -> float) (targetScoreF: 'a -> float) =
        // Gives an array of scores with the frequency of decoy and target hits at that score
        let scoreFrequencies =
            data
            |> Array.map (fun x ->
                if isDecoy x then
                    decoyScoreF x, true
                else
                    targetScoreF x, false
            )
            // groups by score
            |> Array.groupBy fst
            // counts occurences of targets and decoys at that score
            |> Array.map (fun (score,scoreDecoyInfo) ->
                let decoyCount =
                    scoreDecoyInfo
                    |> Array.sumBy (fun (score, decoyInfo) ->
                        match decoyInfo with
                        | true -> 1.
                        | false -> 0.
                    )
                let targetCount =
                    scoreDecoyInfo
                    |> Array.sumBy (fun (score, decoyInfo) ->
                        match decoyInfo with
                        | true -> 0.
                        | false -> 1.
                    )
                createScoreTargetDecoyCount score decoyCount targetCount
            )
            |> Array.sortByDescending (fun x -> x.Score)

        // Goes through the list and assigns each protein a "q value" by dividing total decoy hits so far through total target hits so far
        let reverseQVal =
            scoreFrequencies
            |> Array.fold (fun (acc: (float*float*float*float) list) scoreCounts ->
                let _,_,decoyCount,targetCount = acc.Head
                // Should decoy hits be doubled?
                // accumulates decoy hits
                let newDecoyCount  = decoyCount + scoreCounts.DecoyCount(* * 2.*)
                // accumulates target hits
                let newTargetCount = targetCount + scoreCounts.TargetCount
                let newQVal =
                    let nominator =
                        if newTargetCount > 0. then newTargetCount
                        else 1.
                    newDecoyCount / nominator
                (scoreCounts.Score, newQVal, newDecoyCount, newTargetCount):: acc
            ) [0., 0., 0., 0.]
            // removes last part of the list which was the "empty" initial entry
            |> fun list -> list.[.. list.Length-2]
            |> List.map (fun (score, qVal, decoyC, targetC) -> score, qVal)

        //Assures monotonicity by going through the list from the bottom to top and assigning the previous q value if it is smaller than the current one
        let score, monotoneQVal =
            if reverseQVal.IsEmpty then
                failwith "Reverse qvalues in Storey calculation are empty"
            let head::tail = reverseQVal
            tail
            |> List.fold (fun (acc: (float*float) list) (score, newQValue) ->
                let _,qValue = acc.Head
                if newQValue > qValue then
                    (score, qValue)::acc
                else
                    (score, newQValue)::acc
            )[head]
            |> Array.ofList
            |> Array.sortBy fst
            |> Array.unzip
        // Linear Interpolation
        let linearSplineCoeff = Interpolation.LinearSpline.initInterpolateSorted score monotoneQVal
        // takes a score from the dataset and assigns it a q value
        let interpolation = Interpolation.LinearSpline.interpolate linearSplineCoeff
        interpolation