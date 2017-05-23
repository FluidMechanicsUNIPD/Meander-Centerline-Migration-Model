!-----------------------------------------------------------------------
!  Module containing all the global variables of the code
!-----------------------------------------------------------------------
      module module_global
      implicit none
      
! file names
      character*20 dirIN      ! input directory
      character*20 dirOUT     ! output directory
      character*20 dirTMP     ! temporary directory
      character*20 simname    ! name of the simulation
      character*20 filesim    ! name of the file sim (into dirIN)
      character*20 filexy     ! name of the file with the initial coordinates
      character*100 filename  ! auxiliary file name
      
! flow and sediment
      integer flagbed         ! bed setting (1:from Rp, 2:from indbed)
      integer typebed         ! bed configuration (1:flat bed, 2:dune-covered)
      real*8 beta, beta0, betaold     ! aspect ratio B0/D0 (current, initial, previous)
      real*8 theta, theta0, thetaold  ! Shield parameter (current, initial, previous)
      real*8 ds, ds0, dsold           ! dimensionless grain size ds/D0 (current, initial, previous)
      real*8 Rp                       ! Reynolds particle number Rp=(Delta g ds^3)^0.5/ni
      real*8 rpic, rpic0              ! transverse transport parameter (Talmon)
      integer jmodel          ! flag for flow field model (1:ZS, 2:IPS)
      real*8 Cf, Cf0, Cfold   ! friction coefficient (current, initial, previous)
      real*8 CD, CT           ! derivatives of friction coefficient
      real*8 phi, phiD, phiT  ! sediment transport intensity and derivatives
      real*8 F0               ! Froude number
      
! floodplain structure
      real*8 Ef               ! erodibility coefficient of the floodplain
      real*8 Eb               ! erodibility coefficient of the point bars
      real*8 Eo               ! erodibility coefficient of the oxbow lakes
      
! river configuration
      integer flagxy0         ! initial path configuration (1:random, 2:given)
      integer flag_ox         ! flag for existing floodplain structure (0:no, 1:yes)) 
      integer N, Nold, N0     ! number of points (current, previous, initial)
      integer Nrand           ! number of points interested by a perturbation if flagxy0 = 1
      integer jre             ! removed points before and after a cutoff
      integer jnco            ! minimum threshold points for neck cutoff searching (>2)
      integer ksavgol         ! Savitzky–Golay flag: if <0 number of smoothings per iter, if >0 number of iters between smoothing, if 0 no smoothing
      real*8 deltas, deltas0  ! grid size (current, initial)
      real*8 dsmin, dsmax     ! extreme threshold for grid size (times deltas0)
      real*8 tollc            ! minimum threshold before neck-cutoff
      real*8 stdv             ! standard deviation of initial perturbation flagxy0 = 1
     
! simulation parameters
      integer flag_time       ! final time assignment (1:time, 2:coefficient, 3:iteration number, 4:cutoff number, 5:printed conf number)
      integer flag_dt         ! flag for time marching (1:fixed and equal to dt0, 2:dynamic <dt0)
      integer ivideo          ! number of iterations between two video prints
      integer ifile           ! number of iterations between two files prints
      integer nend	          ! item number (iterations, cutoffs, printed confs)
      real*8 TT, TTs          ! simulation time (dimensionless years) (initial)
      real*8 kTTfco	          ! coefficient for time of first cutoff occurrence
      real*8 tt0              ! starting time of simulation (dimensionless years)
      real*8 dt0              ! fixed time step (dimensionless years)
      real*8 cstab            ! stability coefficient for time marching

!-----------------------------------------------------------------------
      end module module_global
!-----------------------------------------------------------------------
