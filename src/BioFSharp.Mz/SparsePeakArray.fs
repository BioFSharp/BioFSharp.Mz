namespace BioFSharp.Mz

module SparsePeakArray =

    ///
    type SparsePeakArray = {
        Data            : System.Collections.Generic.IDictionary<int,float>
        MzToBinIdx      : float -> int
        BinIdxToMz      : int -> float
        }

    ///
    let dot (x:SparsePeakArray) (y:SparsePeakArray) =
        x.Data
        |> Seq.fold (fun (acc:float) xi ->
            let present,yi = y.Data.TryGetValue xi.Key
            if present then acc + (yi * xi.Value) else acc
            ) 0.

    ///
    let initMzToBinIdx width offset x = int ((x / width) + offset)

    let initBinIdxToMz width offset x =
        ((float x) - offset) * width


    ///
    let peaksToNearestBinVector binWidth offset (minMassBoarder:float) (maxMassBoarder:float) (pkarr:PeakArray<_>) =
        let mzToBinIdx = initMzToBinIdx binWidth offset
        let binIdxToMz = initBinIdxToMz binWidth offset
        let keyValues =
            pkarr
            |> Array.choose (fun p ->
                if p.Mz < maxMassBoarder && p.Mz > minMassBoarder then
                    let index = (mzToBinIdx p.Mz)
                    Some (index, p.Intensity)
                else
                    None
                )
            |> Array.groupBy fst
            |> Array.map (fun (idx,data) -> idx, data |> Array.sumBy snd)
        {
            Data        = keyValues |> dict
            MzToBinIdx  = mzToBinIdx
            BinIdxToMz  = binIdxToMz
        }
