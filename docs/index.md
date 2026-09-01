# BioFSharp.Mz

BioFSharp.Mz is an F# library for computational proteomics. It covers the steps between a
raw mass spectrum and a quantified, statistically controlled protein list: peak
processing, signal detection and centroiding, charge state determination, in silico
peptide fragmentation, peptide search databases, spectrum scoring, quantification, false
discovery rate control, and protein inference. The companion package BioFSharp.Mz.Vis
adds report visualizations for protein inference results.

The library builds on [BioFSharp](https://biofsharp.com/BioFSharp/) for biological data
structures and on [FSharp.Stats](https://fslab.org/FSharp.Stats/) for the numerical
routines.

## Installation

Get the latest package from NuGet:

```shell
dotnet add package BioFSharp.Mz
```

Or reference it in an F# script:

```fsharp
#r "nuget: BioFSharp.Mz"
```

## How this documentation is organized

The pages in the sidebar walk through the library one processing stage at a time, in the
order the stages occur in a proteomics workflow. Each page is a literate F# script that
you can download and run as it stands, and the outputs shown come from runs of the
examples against the library before publication.

Signatures and doc comments for every public function are in the
[API Reference](reference/index.html).

## Contributing and copyright

The project is hosted on [GitHub](https://github.com/BioFSharp/BioFSharp.Mz) where you
can report issues and submit pull requests. The library is available
under the MIT license. For details see the
[License file](https://github.com/BioFSharp/BioFSharp.Mz/blob/developer/LICENSE) in the
GitHub repository.
