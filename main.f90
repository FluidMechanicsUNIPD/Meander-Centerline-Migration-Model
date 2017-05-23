!-----------------------------------------------------------------------
! MEANDER CENTERLINE MIGRATION MODEL - MCMM
! Simulation of the planform evolution of freely-evolving meanders with
! three environments (floodplain, scroll bars, oxbow lakes).
! Last update: 2017/01/02 by Manuel Bogoni
!-----------------------------------------------------------------------
! FLOW FIELD:
! ZS approach (Zolezzi and Seminara, JFM 2001)
! IPS approach (Ikeda, Parker and Sawai, JFM 1981)
! FLOODPLAIN:
! Camporeale et al. JGR 2005, M.Bogoni (PhD Thesis 2016)
!-----------------------------------------------------------------------
      program mcmm

      use module_global
      use module_zs
      use module_geometry
      implicit none
      integer DateTimeI(8), DateTimeE(8)
      integer flagco, ierr, trasp, fondo, RUNSIM, j, jt, jprint
      real*8 tfc, time, year2sec, Tend, Tstart, dt, timestep, dCf, dLs
      parameter ( year2sec = 86400.d0 * 365.d0 )
      real*8, allocatable :: dUx(:), dUy(:), dUxold(:), dUyold(:)

      write(6,*) '#####################################################'
      write(6,*) '##       Meander Centerline Migration Model        ##'
      write(6,*) '#####################################################'

!-----------------------------------------------------------------------
! INPUT DATA
!-----------------------------------------------------------------------

! set directories
      dirIN = ADJUSTL('input/')
      dirOUT = ADJUSTL('output/')
      dirTMP = ADJUSTL('temp/')

! set the name of the file sim
      write(6,*)
      write(6,*) 'Insert file sim:'
      read(5,*) filesim
!      filesim = 'sim.dat'
!      write(6,*)

! load file sim with the characteristic parameters
      call load_sim
      simname = ADJUSTL(simname)
      
! load existing structure of the floodplain
      call load_floodplain

! load or build the initial configuration of the river 
      call load_river
      
! print summary on screen
      call dashline(6)

      write(6,*) 'SIMULATION TIMES'
      write(6,101) 'dt0 =', dt0
      write(6,101) 'TTi =', tt0 
      
      select case (flag_time)
        case (1)
          write(6,101) 'TTe =', TTs
        case (2)
          write(6,101) 'TTe =', kTTfco, '* tfc'
        case (3)
          write(6,100) 'TTe -->', nend, 'iters'
        case (4)
          write(6,100) 'TTe -->', nend, 'cutoffs'
        case (5)
          write(6,100) 'TTe -->', nend, 'prints'
        case default
          stop 'ERROR! Flag for simulation end'
      end select  
      
      call dashline(6)
      write(6,*) 'VALLEY STRUCTURE'
      write(6,102) 'Ef =', Ef
      write(6,102) 'Eb =', Eb
      write(6,102) 'Eo =', Eo
      call dashline(6)
      write(6,*) 'LEADING PARAMETERS'
      write(6,103) 'beta = ', beta0, 'theta =', theta0
      write(6,103) '  ds = ', ds0,   '   Rp =', Rp
      call dashline(6)

      select case (flagbed)
      case(1)
      case(2)
      case default
        stop 'ERROR! Flag for bed set'
      end select

      select case (typebed)
      case(1)
      case(2)
      case default
        stop 'ERROR! Flag for bed configuration'
      end select

100   format(a7, 1x, i6, 1x, a7)
101   format(a5, 1x, f12.3, 1x, a5)
102   format(a4, 1x, e12.6)
103   format(a6, 1x, f12.6, 3x, a7, 1x, f12.6)

!-----------------------------------------------------------------------
! open output files
!-----------------------------------------------------------------------

      filename = TRIM(dirOUT)//TRIM(simname)//'_timesim.dat'
      open(unit=1, file=filename, status='replace')

      filename = TRIM(dirOUT)//TRIM(simname)//'_simulation.dat'
      open(unit=2, file=filename, status='replace')
      
      filename = TRIM(dirOUT)//TRIM(simname)//'_parameters.dat'
      open(unit=3, file=filename, status='replace')

      filename = TRIM(dirOUT)//TRIM(simname)//'_cutoffs.dat'
      open(unit=4, file=filename, status='replace')

      filename = TRIM(dirOUT)//TRIM(simname)//'_lambda.dat'
      open(unit=7, file=filename, status='replace')

!-----------------------------------------------------------------------
! INITIALIZE SIMULATION
!-----------------------------------------------------------------------

! initialize time
      time = tt0 
      
! initialize time of first cutoff occurrence and flag
      tfc = 0.d0
      flagco = 0

! initialize iteration counter
      jt = 0

! initialize printing counter
      jprint = 0

! initialize the system
      call initialize(jt, jprint, time, Tstart, DateTimeI)

! initialize simulation switch
      if ((flag_time.eq.1).and.(time.ge.TTs)) then
        RUNSIM = 0
      else
        RUNSIM = 1
      end if

! initialize arrays
      allocate ( dUx(N), dUy(N), dUxold(N), dUyold(N), stat=ierr )
      if (ierr.ne.0) go to 998

!#######################################################################
!#######################################################################
      do while (RUNSIM.gt.0)
!#######################################################################
!#######################################################################

!-----------------------------------------------------------------------
! TIME MARCHING
!-----------------------------------------------------------------------

! update timestep
        dt = timestep(N, Ef, deltas, dU, cstab, tollc, dt0, flag_dt)
        if (dt.lt.TINY(dt)) stop 'ERROR! Time step tends to zero'

! select case for continuing/stopping the iterations
        select case (flag_time)

! case 1: given final time
        case(1)

! check whether the updated time will be over total period of simulation
          if ((time+dt).ge.TTs) then
            dt = TTs - time
            RUNSIM = 0
          end if

! case 2: final time proportional to the time of first cutoff occurrence
        case(2)

! if first cutoff occurs, then fix the relative occurrence time
          if (flagco.eq.0) then
            if (Nnco.ne.Nncoold) then
              tfc = kTTfco * time
              flagco = 1
              write(6,'(a16,1x,f12.5 )') '## Final time = ', tfc
            end if

! check if current time is larger than the proportional final time
          else if (flagco.eq.1) then
            if ((time+dt).ge.tfc) then
              dt = tfc - time
              RUNSIM = 0
            end if
          end if

! case 3: given number of iterations
        case (3)

          if (jt.ge.nend) then
            RUNSIM = 0
          end if

! case 4: given number of cutoffs
        case (4)

          if (Nnco.ge.nend) then
            RUNSIM = 0
          end if

! case 5: given number of printed configurations
        case (5)

          if (jprint.ge.nend) then
            RUNSIM = 0
          end if

! case default: input flag error
        case default
          stop 'ERROR! Flag for simulation period'

! end select for final time of simulation
        end select

! update simulation time
        time = time + dt

! update iteration counter
        jt = jt + 1

!-----------------------------------------------------------------------
! lateral migration of the river
!-----------------------------------------------------------------------

! if number of node changed, then restore arrays
        if (N.ne.Nold) then
          deallocate (dUx,dUy,dUxold,dUyold, E, whereisnode, stat=ierr)
          if (ierr.ne.0) go to 998
          allocate ( dUx(N), dUy(N), dUxold(N), dUyold(N), stat=ierr)
          if (ierr.ne.0) go to 998
          allocate ( E(N), whereisnode(N), stat=ierr)
          if (ierr.ne.0) go to 998
          whereisnode(:) = 0

! update previous projected velocity
        else
          dUxold(:) = dUx(:)
          dUyold(:) = dUy(:)
        end if

! initialize node counters
        Npf = 0
        Npo = 0
        Npb = 0
        
! loop over points
        do j = 1, N

! project curvature-driven excess velocity onto the reference system
          dUx(j) = dU(j) * DSIN(th(j))
          dUy(j) = dU(j) * DCOS(th(j))

! #### ERODIBILITY BASED ON OXBOW / POINT BAR FORMATION ####
! assign erodibility to the current point: if there is no oxbows, then
! assign the floodplain erodibility
          if ((Nnco.eq.0).or.       &
     &      ( (ABS(Ef-Eo).lt.1.D-12).and.(ABS(Ef-Eb).lt.1.D-12)) ) then
            E(j) = Ef
            Npf = Npf + 1

! assign erodibility to the current point: if there are oxbows, then
! check the position of the point
          else
            call searchpointnew(j)

! node in floodplain
            if (whereisnode(j).eq.0) then
              E(j) = Ef
              Npf = Npf + 1

! node in point bar
            else if (whereisnode(j).gt.0) then
              E(j) = Eb
              Npb = Npb + 1

! node in oxbow lake
            else if (whereisnode(j).lt.0) then
              E(j) = Eo
              Npo = Npo + 1
            end if

! end if for erodibility assigning
          end if

! lateral migration of the channel axis
          if ((N.eq.Nold).and.(jt.gt.1)) then
            x(j) = x(j) + dt*E(j) * 0.5d0 * (dUx(j)+dUxold(j)) * year2sec
            y(j) = y(j) + dt*E(j) * 0.5d0 * (dUy(j)+dUyold(j)) * year2sec
          else
            x(j) = x(j) + dt * E(j) * dUx(j) * year2sec
            y(j) = y(j) + dt * E(j) * dUy(j) * year2sec
          end if

! end loop over points
        end do

! check for points
        if ((Npf+Npb+Npo).ne.N) then
          stop 'ERROR! Point missing in floodplain structure seeking'
        end if

!-----------------------------------------------------------------------
! cutoff seeking
!-----------------------------------------------------------------------

! update point number
        Nold = N

! update cutoff counter
        Nncoold = Nnco
              
! look for neck cutoffs
        call neckcutoff(time,4,jt)

! print oxbow and point bar arrays
        if (Nnco.ne.Nncoold) then
          call print_floodplain
        end if

!-----------------------------------------------------------------------
! river geometry
!-----------------------------------------------------------------------

! point resampling
        call resampling
        Nsampold = Nsamp

! if number of node changed, then restore arrays
        if (N.ne.Nold) then
          deallocate (th, c, dU, stat=ierr)
          if (ierr.ne.0) go to 998
          allocate (th(N), c(N), dU(N), stat=ierr)
          if (ierr.ne.0) go to 998
        end if
        
! angle with longitudinal direction
        call angle(N, x, y, th)

! axis curvature (multiplied by v0)
        call curv(N, s, th, c, ksavgol, jt)

! update intrinsit length of the river
        Lsold = Ls
        Ls = s(N)
        
! update Cartesian longitudinal length of the river
        Lxold = Lx
        Lx = x(N)-x(1)
!        Lx = DSQRT((x(N0)-x(1))**2.d0 + (y(N0)-y(1))**2.d0)
        
!-----------------------------------------------------------------------
! flow resistance and sediment transport
!-----------------------------------------------------------------------
        Cfold = Cf
        call resistance(1, 2, 0, trasp, fondo)
        
!-----------------------------------------------------------------------                      
! update leading parameters
!-----------------------------------------------------------------------
        dCf = Cfold / Cf
        dLs = (Lx/Ls) / (Lxold/Lsold)
        beta  = beta  * dLs**(1.d0/3.d0) * dCf**( 1.d0/3.d0)
        theta = theta * dLs**(2.d0/3.d0) * dCf**(-1.d0/3.d0)
        ds    = ds    * dLs**(1.d0/3.d0) * dCf**( 1.d0/3.d0)

!-----------------------------------------------------------------------
! curvature-driven flow field
!-----------------------------------------------------------------------

! initialize array of excess near-bank longitudinal velocity
        dU(:) = 0.d0

! set flow field model
        select case(jmodel)
! ZS model      
        case(1)
          call coefZS(jt)
          call dUZS
! IPS model
        case(2)
          call dUIPS
        end select
      		  	                            
!-----------------------------------------------------------------------
! print results
!-----------------------------------------------------------------------

! print on file
        if ((MOD(jt,ifile).eq.0).or.(RUNSIM.eq.0)) then
          jprint = jprint + 1
          call printsim(jt, jprint, time, 2, 3, trasp, fondo)
          call printlambda(jt, jprint, time, 7)
        end if

! print on screen
        if ((MOD(jt,ivideo).eq.0).or.(RUNSIM.eq.0)) then
          call printscreen(jt, jprint, time, Nnco)
        end if

! print backup
        filename = TRIM(dirTMP)//'xy.dat'
        filename = ADJUSTL(filename)
        open(unit=501, file=filename, status='replace')
        
        write(501,88) jt, 'jt'
        write(501,89) time, 'time'
        write(501,89) beta, 'beta'
        write(501,89) theta, 'theta'
        write(501,89) ds, 'ds'
        write(501,89) Cf, 'Cf'
        write(501,88) jprint, 'jp'
        write(501,88) N, 'N'
        call dashline(501)
        do j = 1, N
          write(501,*) x(j), y(j)
        end do
        close(501)
                 
!#######################################################################
!#######################################################################
! end loop over time
      end do
!#######################################################################
!#######################################################################

! close files
      close(2) ! simulation
      close(3) ! parameters
      close(4) ! cutoffs
      close(7) ! lambda

! reset arrays
      deallocate (x0, y0, x, y, s, th, c, dU, dUx, dUy, E, whereisnode)
      if (ALLOCATED(dUxold)) deallocate (dUxold, dUyold)
      if (ALLOCATED(Am)) deallocate (Am, lamb1, lamb2, lamb3, lamb4)
      if (ALLOCATED(Am)) deallocate (g10,g20,g30,g40,g11,g21,g31,g41)

!-----------------------------------------------------------------------
! compute CPU time
!-----------------------------------------------------------------------
      call time2sec(DateTimeE, Tend)
      call printtime(DateTimeE, 1, 1)
      time = (Tend - Tstart) / 1000.d0

! print on file
      write(1,85)'Computation time    : ', time,         ' seconds'
      write(1,86)'Computation time    : ', time/3600.d0, '   hours'
      close(1)

! print on screen
      call dashline(6)
      write(6,82)'Computation time : ', time,         ' seconds'
      write(6,83)'Computation time : ', time/3600.d0, '   hours'   
      call dashline(6)           
      write(6,*) '...THE END!'

!-----------------------------------------------------------------------
! printing formats
!-----------------------------------------------------------------------
82    format(a19,1x,e12.5,2x,a8)
83    format(a19,1x,f12.5,2x,a8)
85    format(a22,1x,e12.5,2x,a8)
86    format(a22,1x,f12.5,2x,a8)
88    format(5x,i10,5x,a2)
89    format(f15.7,2x,a5)

!-----------------------------------------------------------------------
! end of code
!-----------------------------------------------------------------------
      stop
998   write(6,*) 'ERROR! Allocation in main'
      stop
      end
