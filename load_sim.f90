!-----------------------------------------------------------------------
! LOAD FILE SIM: river and simulation parameters
!----------------------------------------------------------------------- 
      subroutine load_sim

      use module_global
      use module_zs
      implicit none

! open file
      filesim = ADJUSTL(filesim)
      filename = TRIM(dirIN)//TRIM(filesim)
      open(unit=1, file=filename, status='old')
      
! intro
      read(1,*)
      read(1,*) simname    ! name of the simulation
!      read(1,*) dirOUT    ! output directory

! morphodynamics
      read(1,*)
      read(1,*) beta0      ! aspect ratio
      read(1,*) theta0     ! Shields number
      read(1,*) ds0        ! grain roughness
      read(1,*) flagbed    ! bedset    (1:from Rp, 2:from bed input)              ex IndBed
      read(1,*) Rp         ! particle Reynolds Rp=(Delta g ds^3)^0.5/ni
      read(1,*) typebed    ! bed configuration (1:flat bed, 2:dune covered)       ex IndC
      read(1,*) rpic0      ! transverse transport parameter (Talmon)
      read(1,*) jmodel     ! flag for flow field model (1:ZS, 2:IPS)
      read(1,*) Nz         ! number of points for the vertical integration of the flow field
      read(1,*) Mdat       ! order of Fourier expansion

! geomorphology
      read(1,*)
      read(1,*) Ef         ! erodibility coefficient of the floodplain
      read(1,*) Eb         ! erodibility coefficient of the point bars
      read(1,*) Eo         ! erodibility coefficient of the oxbow lakes
      read(1,*) flag_ox    ! flag for existing floodplain structure (0:no, 1:yes)) 

! river configuration
      read(1,*)
      read(1,*) N0	   ! initial number of points
      read(1,*) flagxy0    ! initial path configuration (1:random, 2:given)
      read(1,*) filexy     ! name of geometry file
      read(1,*) deltas0    ! distance between axis points
      read(1,*) dsmin      ! min value of grid size (times deltas0)
      read(1,*) dsmax      ! max value of grid size (times deltas0)
      read(1,*) Nrand      ! number of point interested by a perturbation if flagxy0 = 1
      read(1,*) stdv       ! standard deviation of initial perturbation if flagxy0 = 1
	read(1,*) tollc      ! minimum threshold before neck-cutoff
      read(1,*) jre        ! removed points before and after a cutoff
      read(1,*) jnco       ! minimum threshold points for neck cutoff searching (>2)   
      read(1,*) ksavgol    ! Savitzky–Golay flag: if <0 number of smoothings per iter, if >0 number of iters between smoothing, if 0 no smoothing
      
! time step and printing
      read(1,*)
      read(1,*) flag_time  ! final time assignment (1:time, 2:coefficient, 3:iteration number, 4:cutoff number, 5:printed conf number)
      read(1,*) TTs        ! simulation time if flag_time = 1 (dimensionless years)
      read(1,*) kTTfco	   ! coefficient for first cutoff time if flag_time = 2
      read(1,*) nend	   ! item number (iterations, cutoffs, printed confs)
      read(1,*) tt0        ! starting time of simulation (dimensionless years)
      read(1,*) flag_dt    ! flag for time marching (1:fixed and equal to dt0, 2:dynamic <dt0)
      read(1,*) dt0        ! fixed time step (dimensionless years)
      read(1,*) cstab      ! coefficient for time marching
      read(1,*) ivideo     ! number of iterations between two video prints
      read(1,*) ifile      ! number of iterations between two files prints
      
! valley boundaries
      read(1,*)
      read(1,*) jbound           ! transition shape (0: no boundary, 1:sharp,2:spline,3:exp-decaying,4:exp-reversed)
      read(1,*) Ebound           ! erodibility coefficient of the valley boundaries
      read(1,*) Lhalfvalley      ! transverse half-width of the usual floodplain
      read(1,*) Ltransvalley     ! thickness of the transition layer
       
! close file sim
      close(1)

! print summary on screen
      call dashline(6)
      write(6,*) 'INTRO'
      write(6,*) 'simulation name : ', simname
      select case(flagxy0)
      case(1)
        write(6,*) 'straight river with slight random perturbation'
      case(2)
        write(6,*) 'coordinate file : ', filexy
      case default
        stop 'ERROR! Wrong flag for initial river nfiguration'
      end select
      
! end of the subroutine
      return
      end
