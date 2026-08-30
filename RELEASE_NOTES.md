#### 0.2.0 (Released 2026-08-29)

* Upgrade to BioFSharp 2.0.0 (bundles the former BioFSharp.IO package) and migrate the codebase over the breaking API changes
* Align the dependency stack: FSharpAux / FSharpAux.IO 2.1.0, FSharp.Stats 0.6.0, Newtonsoft.Json 13.0.4, Plotly.NET 6.0.0-preview.2
* Replace the legacy FAKE 5 build script with the FAKE 6 build project ported from BioFSharp (build.cmd / build.sh wrappers, GitHub Actions build-test and gh-pages docs deployment, fsdocs 20 documentation pipeline)

#### 0.1.5 - Friday, February 19, 2021
* SignalDetection: add option for peak summation and intensity weighted based mz refinement

#### 0.1.4 - Friday, Oktober 9, 2020
* Minor changes 
* SearchDB: added default inclusion of peptides resulting of protein n terminal methionin cleavages.
* FDRControl: added more stable qValue calculation, added storeys qvalue method
* ProteinInference: added FDR calculations

#### 0.1.3 - Monday, March 16, 2020
* Minor changes 
* Quantification: improved handling of unsuccessful fits

#### 0.1.2 - Thursday, December 12, 2019
* Minor changes 

#### 0.1.1 - Thursday, December 12, 2019
* Minor changes 

#### 0.1.0 - Thursday, December 12, 2019
*
#### 0.0.117 - Friday, September 13, 2019
* include faster hasFlag function
* speed improvements in SequestLike: faster computation of theoretical spectra. 

#### 0.0.116 - Saturday, September 7, 2019
* refactor according to changes latest BioFSharp Release

#### 0.0.115 - Friday, September 6, 2019
* add XScoring, a combination of XTandemLike and AndromedaLike Scoring

#### 0.0.114 - Wednesday, September 4, 2019
* add ProteinInference

#### 0.0.113 - Thursday, August 1, 2019
* refactor Fragmentation module according to changes in BioFSharp.Mz release 0.1.1

#### 0.0.112 - Thursday, April 25, 2019
* paket template fix

#### 0.0.111 - Thursday, April 25, 2019
* Initial release

#### 0.0.1 - Thursday, April 25, 2019
* Initial release
