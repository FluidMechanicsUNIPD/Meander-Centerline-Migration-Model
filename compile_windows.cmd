::---------------------------------------------------------------------
:: COMMAND FILE TO COMPILE THE SOURCE CODE OF MEANDER MIGRATION MODEL
:: last update: 2017/01/02 by Manuel Bogoni
::---------------------------------------------------------------------

:: set the directory with the source code
set "dir=fortran"

::  set the name of the executable
set "exe=mcmm.exe"

::  set the Fortran compiler
set "fc90=gfortran"

::  clean the dashboard
cls

::  enter the directory
cd %dir%

::  compile the modules
%fc90% -c module_global.f90 module_geometry.f90 module_zs.f90

::  compile the code
%fc90% -w -static -fcheck=all  ^
main.f90 ^
angle.f  ^
banks.f ^
coefZS.f90 ^
curv.f ^
dashline.f90 ^
derivative.f ^
dUIPS.f ^
dUZS.f ^
findcutoff.f ^
initialize.f90 ^
integrals_cavsimp.f ^
integrals_semiana_functions.f ^
k0134.f ^
load_floodplain.f90 ^
load_river.f90 ^
load_sim.f90 ^
module_global.f90 ^
module_geometry.f90 ^
module_zs.f90 ^
neckcutoff.f ^
print_floodplain.f90 ^
printcutoff.f90 ^
print_lambda.f90 ^
printscreen.f90 ^
printsim.f90 ^
printtime.f90 ^
randomconfig.f90 ^
randomnormal.f90 ^
resampling.f ^
resistance.f ^
rungeg.f90 ^
savgolfilter_new.f ^
searchpoint_new.f90 ^
smoothing.f ^
storeoxbow.f ^
storepointbar.f ^
time2sec.f90 ^
timestep.f ^
windingnumber.f90 ^
zroots.f  ^
-o %exe%

::  copy executable code
copy %exe% ..

::  exit directory
cd ..
pause
