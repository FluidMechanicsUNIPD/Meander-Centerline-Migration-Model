This repository contains the source codes of the numerical model developed by [Bogoni et al. WRR 2017](http://onlinelibrary.wiley.com/doi/10.1002/2017WR020726/abstract) and featured as AGU Research Spotlight by [EOS](https://doi.org/10.1029/2017EO078667).

[![DOI](https://zenodo.org/badge/95902287.svg)](https://zenodo.org/badge/latestdoi/95902287)

CODE PROPERTIES
* The source code is written in Fortran 77/90 language, thus a Fortran compiler is required (e.g. gfortran).
* The code can be compiled through the script _compile_unix_ (for UNIX/macOS) or _compile_windows.cmd_ (for Windows).
* Information about the model and the sctructure of input and output files are reported in the file _README.pdf_
* The folder _input_ contains a couple of templates.
* Input files must be into the _input_ folder, located in the same directory of the executable file.
* Folders _temp_ and _output_ are required in the same directory of the executable file.

See:

[MCMM website](https://fluidmechanicsunipd.github.io/Meander-Centerline-Migration-Model)

[MCMM GitHub](https://github.com/FluidMechanicsUNIPD/Meander-Centerline-Migration-Model)

[CSDMS website](http://csdms.colorado.edu/wiki/Model:Meander_Centerline_Migration_Model)

[CSDMS GitHub](https://github.com/csdms-contrib/Meander-Centerline-Migration-Model)

Please cite:
Bogoni, M., Putti M., and Lanzoni S. (2017), Modeling meandermorphodynamics over self-formed heterogeneous floodplains, Water Resour. Res. 53, 5137–5157, doi:10.1002/2017WR020726
