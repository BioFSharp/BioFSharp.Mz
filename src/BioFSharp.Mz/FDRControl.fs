namespace BioFSharp.Mz

open System
open FSharp.Stats
open FSharpAux
open FSharp.Stats.Fitting
open FSharp.Stats.Fitting.NonLinearRegression

module FDRControl = 

    module MAYU =

        // FDR estimation using MAYU
        // Code form 'stirlingLogFactorial' to 'estimatePi0HG' translated from percolator 'ProteinFDREstimator.cpp'
        let private stirlingLogFacorial (n: float) =
            log(sqrt(2. * Ops.pi  *n)) + n * log(n) - n

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
            if (Ops.isNan expectation_value_FP_PID) || (Ops.isInf expectation_value_FP_PID) then
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
        let linearSplineCoeff = Interpolation.LinearSpline.interpolateSorted score monotoneQVal
        // takes a score from the dataset and assigns it a q value
        let interpolation = Interpolation.LinearSpline.predict linearSplineCoeff
        interpolation

    /// for given data, creates a logistic regression model and returns a mapping function for this model
    let getLogisticRegressionFunction (x:vector) (y:vector) epsilon =
        let alpha =
            match FSharp.Stats.Fitting.LogisticRegression.Univariable.estimateAlpha epsilon x y with
            | Some a -> a
            | None -> failwith "Could not find an alpha for logistic regression of fdr data"
        let weight = FSharp.Stats.Fitting.LogisticRegression.Univariable.coefficient epsilon alpha x y
        FSharp.Stats.Fitting.LogisticRegression.Univariable.predict weight

    /// Creates a Histogram based on a given score of a target/decoy dataset. Each bin contains the information of the total count, the decoy count and the median score.
    /// (Bin, Count, DecoyCount, Median Score)
    let createTargetDecoyHis bandwidth (isDecoy: 'a -> bool) (decoyScoreF: 'a -> float) (targetScoreF: 'a -> float) (data: 'a[]) =
        let halfBw = bandwidth / 2.0
        let scoreDecoyInfo =
            data
            |> Array.map (fun x ->
                if isDecoy x then
                    {|Score = decoyScoreF x; Decoy = true|}
                else
                    {|Score = targetScoreF x; Decoy = false|}
            )
        scoreDecoyInfo
        |> Array.groupBy (fun x ->
            floor (x.Score / bandwidth))
        |> Array.map (fun (k,values) ->
            let count = (Array.length(values))
            let decoyCount = (values |> Array.filter (fun x -> x.Decoy = true) |> Array.length)
            let medianScore = values |> Array.map (fun x -> x.Score) |> Array.median
            // first part of the tuple only needed for debugging
            if k < 0. then
                ((k  * bandwidth) + halfBw, count, decoyCount, medianScore)
            else
                ((k + 1.) * bandwidth - halfBw, count, decoyCount, medianScore)
        )

    /// Calculates the PEP value based on the ratio of Decoys to targets at a given score
    let calculatePEPValues (totalCountF: 'a -> float) (decoyCountF: 'a -> float) (scoreF: 'a -> float) (dataFreq: 'a[]) =
        dataFreq
        |> Array.map (fun x ->
            scoreF x,(decoyCountF x)/(totalCountF x)
        )
        |> Array.sortBy fst
        |> Array.toList

    /// Logit transforms pep values (log10)
    let logitTransformPepValues score pepVal  =
        Array.zip score pepVal
        // 0 and 1 are + and - infinity
        |> Array.filter (fun (y,x) -> x <> 0. && x <> 1.)
        |> Array.map (fun (score,pep) ->
            score,
            log10 (pep/(1.-pep))
        )
        |> Array.unzip

    /// Calculates monotonized PEP values for a target/decoy dataset based on the decoy/target ratio. Entries are binned with a given bandwidth as intital estiamtor based on the scores.
    /// Returns a function which maps from score to PEP value based on a fit of a linear function using linear regression. The linear regression is performed on the logit transformed
    /// pep values. The fit focuses on the pep values centered aound the middle of the score distribution
    let initCalculateLin (trace: string -> unit) bandwidth (isDecoy: 'a -> bool) (decoyScoreF: 'a -> float) (targetScoreF: 'a -> float) (data: 'a[]) =
        let lowerScore, upperScore =
            let decoy =
                data
                |> Array.filter isDecoy
                |> Array.map decoyScoreF
                |> Array.filter (fun x -> x < 0.)
                |> Array.median
            let target =
                data
                |> Array.filter (isDecoy >> not)
                |> Array.map targetScoreF
                |> Array.filter (fun x -> x > 0.)
                |> Array.median
            decoy, target
        trace (sprintf "Lower Score: %f; Upper Score: %f" lowerScore upperScore)
        let filteredData =
            data
            |> Array.filter (fun entry ->
                if isDecoy entry then
                    let score = decoyScoreF entry
                    score >= lowerScore && score <= upperScore
                else
                    let score = targetScoreF entry
                    score >= lowerScore && score <= upperScore
            )
        trace (sprintf "Initial Bandwidth: %f" bandwidth)
        let fittingFunction, score, pep =
            let xPointRange =
                let min = Math.Min((Array.minBy targetScoreF filteredData) |> targetScoreF, (Array.minBy decoyScoreF filteredData) |> decoyScoreF)
                let max = Math.Max((Array.maxBy targetScoreF filteredData) |> targetScoreF, (Array.maxBy decoyScoreF filteredData) |> decoyScoreF)
                max-min
            let upperBW = Math.Min(10., xPointRange/10.)
            [|bandwidth .. 0.1 .. upperBW|]
            |> Array.choose (fun bw ->
                let targetDecoyHis = createTargetDecoyHis bw (isDecoy: 'a -> bool) (decoyScoreF: 'a -> float) (targetScoreF: 'a -> float) (filteredData: 'a[])
                let score',pep' =
                    calculatePEPValues (fun (_,count,_,_) -> float count) (fun (_,_,decoyCount,_) -> float decoyCount) (fun (_,_,_,medianScore) -> medianScore) targetDecoyHis
                    |> Array.ofList
                    |> Array.unzip
                let logitScore, logitPEPVal = logitTransformPepValues score' pep'
                let coeff = Fitting.LinearRegression.OLS.Linear.Univariable.fit (vector logitScore) (vector logitPEPVal)
                let fittingFunction' = (Fitting.LinearRegression.OLS.Linear.Univariable.predict coeff) >> (fun x -> 10.**(x)/(1.+10.**(x)))
                let sos = FSharp.Stats.Fitting.GoodnessOfFit.calculateSumOfSquares fittingFunction' score' pep'
                if coeff.[1] < 0. then
                    Some (sos.Error/sos.Count, fittingFunction', score', pep', bw)
                else
                    None
            )
            |> Array.minBy (fun (error,_,_,_,_) -> error)
            |> fun (error, fit,s,p,bw) -> trace (sprintf "Chosen Bandwidth: %f" bw); fit,s,p
        fittingFunction
