!------------------------------------------------------------------------
! SYSTEM INITIALIZATION
!------------------------------------------------------------------------
      subroutine initialize(jt, jprint, time, Tstart, DateTimeI)
      use module_global
      use module_zs
      use module_geometry
      implicit none
      integer j, jt, jprint, ierr, trasp, fondo
      integer DateTimeI(8)
      real*8 time, Tstart

!-----------------------------------------------------------------------
! river geometry
!-----------------------------------------------------------------------

! number of axis points
      N = N0
      Nold = N0

! initialize arrays for the geometry of the river
      allocate ( x(N), y(N), s(N), stat=ierr )
      if (ierr.ne.0) go to 997
      x(:) = x0(:)
      y(:) = y0(:)

! point resampling
      Nsamp = 0
      call resampling
      Nsampold = Nsamp

! initialize arrays for the geometry of the river
      allocate ( th(N), c(N), stat=ierr )
      if (ierr.ne.0) go to 997

! angle with longitudinal direction
      call angle(N, x, y, th)

! axis curvature (multiplied by v0)
      call curv(N, s, th, c, ksavgol, jt)

! update intrinsit length of the river
      Ls0 = s(N)
      Lsold = s(N)
      Ls = s(N)

! update Cartesian longitudinal length of the river
      Lx0 = x(N)-x(1)
!      Lx0 = DSQRT((x(N0)-x(1))**2.d0 + (y(N0)-y(1))**2.d0)
      Lxold = x(N)-x(1)
      Lx = x(N)-x(1)
!      Lx = DSQRT((x(N0)-x(1))**2.d0 + (y(N0)-y(1))**2.d0)

! number of cutoffs
      Nnco = Nnco0

! initialize arrays for the floodplain structure
      allocate ( E(N), whereisnode(N), stat=ierr )
      if (ierr.ne.0) go to 997
      E(:) = 0.d0
      whereisnode(:) = 0

!-----------------------------------------------------------------------
! flow resistance and sediment transport
!-----------------------------------------------------------------------

! initialize leading parameters
      beta = beta0
      theta = theta0
      ds = ds0

      call resistance(1, 2, 1, trasp, fondo)
      Cf0 = Cf
      Cfold = Cf

!-----------------------------------------------------------------------
! curvature-driven flow field
!-----------------------------------------------------------------------

! initialize arrays for the flow field
      allocate ( dU(N), stat=ierr )
      if (ierr.ne.0) go to 997
      dU(:) = 0.d0

! initialize ZS flow field
      if (jmodel.eq.1) then
        allocate (Am(Mdat))
        allocate (lamb1(Mdat), lamb2(Mdat), lamb3(Mdat), lamb4(Mdat))
        allocate (g10(Mdat), g20(Mdat), g30(Mdat), g40(Mdat))
        allocate (g11(Mdat), g21(Mdat), g31(Mdat), g41(Mdat))
      end if

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
      call findcutoff(time,0,0,N,4,jt)
      call printsim(jt, jprint, time, 2, 3, trasp, fondo)
      call printlambda(jt, jprint, time, 7)

! flag for continuing the simulation
      write(6,*) 'System initialized. Press any key to continue.'
      read(5,*)

! start and print CPU time
      call time2sec(DateTimeI, Tstart)
      call printtime(DateTimeI, 1, 0)

! print on screen
      call printscreen(jt, jprint, time, Nnco)

! print backup
      filename='temp/xy.dat'
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

!-----------------------------------------------------------------------
! end of the subroutine
!-----------------------------------------------------------------------

! printing formats
88    format(5x,i10,5x,a2)
89    format(f15.7,2x,a5)

      return
997   write(6,*) 'ERROR! Allocation in initial configuration'
      stop
      end
