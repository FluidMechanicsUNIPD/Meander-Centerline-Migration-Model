!-----------------------------------------------------------------------
! LOCAL ERODIBILITY VALUE OF A FLOODPLAIN WITH GIVEN SYMMETRIC BUNDARIES
!-----------------------------------------------------------------------
! jshape : switch for the shape of the floodplain
!        (0:no boundaries, 1: straight boundaries, 2:bend floodplain, 3:specified shape)
! jbound : switch for the kind of boundary transition
!        (1:sharp, 2:spline, 3:exp-decaying, 4:exp-reversed)
! i    : index of the current point
! N    : number of points of axis river
! x, y : coordinates of the axis river
!-----------------------------------------------------------------------
       real*8 function boundaries (jshape, i, N, x, y, E)
       implicit none
       integer jshape, jbound, i, N
       real*8 a, b, c, d, E, Ef, Eb, Lf, Lb, phi
       real*8 x(N), y(N)
       parameter ( phi = 3.1415927d0 )

!-----------------------------------------------------------------------
! input data
!-----------------------------------------------------------------------

! open file and read boundary limit and erodibilities
       open(unit=72, file='input/boundary_transition.dat', status='old')

! read transition shape (1:sharp,2:spline,3:exp-decaying,4:exp-reversed)
       read(72,*) jbound

! read boundary limits
       read(72,*) Lf     ! transverse half-width of the usual floodplain
       read(72,*) Lb     ! thickness of the transtion layer

! read erodibilities
       read(72,*) Ef     ! erodibility of the usual floodplain
       read(72,*) Eb     ! erodibility of the lateral boundary

! close file
       close(72)

!-----------------------------------------------------------------------
! erodibility for the current point
!-----------------------------------------------------------------------

! set boundary transition shape
       select case (jbound)

!-----------------------------------------------------------------------
! no boundaries
         case (0)
!-----------------------------------------------------------------------
           boundaries = Ef

!-----------------------------------------------------------------------
! sharp transition
         case (1)
!-----------------------------------------------------------------------

           if ( ABS(y(i)).le.(Lf+Lb*0.5d0) ) then
             boundaries = Ef
           else
             boundaries = Eb
           end if

!-----------------------------------------------------------------------
! cubic spline
         case (2)
!-----------------------------------------------------------------------

           if ( ABS(y(i)).le.Lf ) then
             boundaries = Ef
           else if ( ABS(y(i)).gt.(Lf+Lb) ) then
             boundaries = Eb
           else
             call cubecoeff(Lf, Ef, 0.d0, Lf+Lb, Eb, 0.d0, a, b, c, d)
             boundaries = a * (ABS(y(i)))**3.d0 +  &
     &                     b * (ABS(y(i)))**2.d0 +  &
     &                     c * ABS(y(i)) + d
           end if

!-----------------------------------------------------------------------
! exponential decaying
         case (3)
!-----------------------------------------------------------------------

           if ( ABS(y(i)).le.Lf ) then
             boundaries = Ef
           else
             boundaries = Eb + (Ef-Eb) * EXP(-phi/Lb * (ABS(y(i))-Lf))
           end if

!-----------------------------------------------------------------------
! reversed exponential
         case (4)
!-----------------------------------------------------------------------

           if ( ABS(y(i)).le.Lf ) then
             boundaries = Ef
           else if ( ABS(y(i)).gt.(Lf+Lb) ) then
             boundaries = Eb
           else
             boundaries = MAX(Eb, 1.d0+Ef-EXP(Ef/Lb*(ABS(y(i))-Lf)) )
           end if

!-----------------------------------------------------------------------
! error
         case default
!-----------------------------------------------------------------------
           stop 'ERROR! Wrong flag for boundary transition'

!-----------------------------------------------------------------------
! end transition setting
       end select

! end of function
       return
       end