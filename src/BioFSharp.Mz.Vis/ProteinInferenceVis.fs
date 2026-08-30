namespace BioFSharp.Mz.Vis

open Plotly.NET
open BioFSharp.Mz.ProteinInference

module ProteinInference =

    let qValueHitsVisualization bandwidth (inferredProteinClassItemQValue: InferredProteinClassItemQValue[]) path (groupFiles: bool) =
        let decoy, target =
            inferredProteinClassItemQValue
            |> Array.partition (fun x -> x.InfProtClassItem.DecoyHasBetterScore)
        // Histogram with relative abundance
        let relFreqTarget = 
            FSharp.Stats.Distributions.Frequency.create bandwidth (target |> Array.map (fun x -> x.InfProtClassItem.TargetScore))
            |> Map.toArray
            |> Array.map (fun x -> fst x, (float (snd x)) / (float target.Length))
        let relFreqDecoy  =
            FSharp.Stats.Distributions.Frequency.create bandwidth (decoy |> Array.map (fun x -> x.InfProtClassItem.DecoyScore))
            |> Map.toArray
            |> Array.map (fun x -> fst x, (float (snd x)) / (float target.Length))
        // Histogram with absolute values
        let absFreqTarget =
            FSharp.Stats.Distributions.Frequency.create bandwidth (target |> Array.map (fun x -> x.InfProtClassItem.TargetScore))
            |> Map.toArray
        let absFreqDecoy  =
            FSharp.Stats.Distributions.Frequency.create bandwidth (decoy |> Array.map (fun x -> x.InfProtClassItem.DecoyScore))
            |> Map.toArray
        let histogram =
            [
                Chart.Column relFreqTarget
                |> Chart.withTraceInfo (Name = "Target")
                |> Chart.withAxisAnchor(Y=1);
                Chart.Column relFreqDecoy
                |> Chart.withTraceInfo (Name = "Decoy")
                |> Chart.withAxisAnchor(Y=1);
                Chart.Column absFreqTarget
                |> Chart.withAxisAnchor(Y=2)
                |> Chart.withMarkerStyle (Opacity = 0.)
                |> Chart.withTraceInfo (ShowLegend = false);
                Chart.Column absFreqDecoy
                |> Chart.withAxisAnchor(Y=2)
                |> Chart.withMarkerStyle (Opacity = 0.)
                |> Chart.withTraceInfo (ShowLegend = false)
            ]
            |> Chart.combine

        let sortedQValues =
            inferredProteinClassItemQValue
            |> Array.map (fun x ->
                if x.InfProtClassItem.Decoy then
                    x.InfProtClassItem.DecoyScore, x.QValue
                else
                    x.InfProtClassItem.TargetScore, x.QValue
            )
            |> Array.sortBy fst

        [
            Chart.Point sortedQValues
            |> Chart.withTraceInfo (Name = "Q-Values");
            histogram
        ]
        |> Chart.combine
        |> Chart.withYAxisStyle("Relative Frequency / Q-Value",Side=StyleParam.Side.Left,Id=StyleParam.SubPlotId.YAxis 1, MinMax = (0., 1.))
        |> Chart.withYAxisStyle("Absolute Frequency",Side=StyleParam.Side.Right,Id=StyleParam.SubPlotId.YAxis 2,Overlaying=StyleParam.LinearAxisId.Y 1, MinMax = (0., float target.Length))
        |> Chart.withXAxisStyle "Score"
        |> Chart.withSize (900, 900)
        |> if groupFiles then
               Chart.saveHtml (path + @"\QValueGraph")
           else
               Chart.saveHtml (path + @"_QValueGraph")
