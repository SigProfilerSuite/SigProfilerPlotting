# SigProfilerPlotting #
----------

SigProfilerPlotting provides a standard tool for displaying mutational
signatures and mutational patterns in cancer genomes. It integrates with the
other SigProfiler tools and accepts their tabular matrix outputs.

## Current capabilities

- Static SBS, DBS, ID, CNV, and SV plots, including SBS4608 and the ID89 and
  ID476 indel schemes.
- A dedicated SBS4608 interactive HTML report and portable multi-page PDF.
- A local interactive dashboard with search, sample/signature selection,
  percentage/proportion/counts controls, SVG/PDF export, a choice of four
  themes, and optional COSMIC comparison using a user-supplied authorized
  catalogue.
- One-click dashboard launchers for macOS, Linux, and Windows.
- Combined sample portraits assembled from SBS, DBS, and ID matrices.

Start with the [Python installation guide](Installation-Python.md), then see
[Plotting Substitutions](Plotting-Substitutions.md),
[Plotting ID89/ID476](Plotting-ID89-ID476.md), or the
[Interactive Dashboard](Plotting-Dashboard.md).

The SigProfilerPlotting library can be found on GitHub [here][2].
The [README.md][3] file provides a great overview if you already have python and pip installed on your computer. For users that prefer working in an R environment, we provide an R wrapper (SigProfilerPlottingR) that can be found on github [here][12].

## Citation ##
Bergstrom EN, Huang MN, Mahto U, Barnes M, Stratton MR, Rozen SG, and Alexandrov LB (2019) SigProfilerMatrixGenerator: a tool for visualizing and exploring patterns of small mutational events. **BMC Genomics** 20, Article number: 685.
https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-6041-2
## License ## 

Copyright (c) 2019, Erik Bergstrom [Alexandrov Lab]
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## Contact ##

Please address any queries or bug reports to S M Ashiqul Islam at
sislam6@albany.edu.

  [1]: https://files.osf.io/v1/resources/2aj6t/providers/osfstorage/5c73d34f8d5d98001a393840?mode=render
  [2]: https://github.com/SigProfilerSuite/SigProfilerPlotting
  [3]: https://github.com/SigProfilerSuite/SigProfilerPlotting/blob/master/README.md
  [4]: https://realpython.com/installing-python/
  [5]: https://docs.brew.sh/Homebrew-and-Python
  [6]: https://files.osf.io/v1/resources/2aj6t/providers/osfstorage/5c6b0c28737a510019aa3e15?mode=render
  [7]: https://github.com/SigProfilerSuite/SigProfilerPlotting
  [8]: https://files.osf.io/v1/resources/2aj6t/providers/osfstorage/5c5de3d676653c0018270138?mode=render
  [9]: https://files.osf.io/v1/resources/2aj6t/providers/osfstorage/5c73d2c482a3950017d09b92?mode=render
  [10]: https://files.osf.io/v1/resources/2aj6t/providers/osfstorage/5c73c72862c82a0019dae02f?mode=render
  [11]: https://files.osf.io/v1/resources/2aj6t/providers/osfstorage/5c73c7eb82a3950018cbfddc?mode=render
  [12]: https://github.com/SigProfilerSuite/SigProfilerPlottingR
