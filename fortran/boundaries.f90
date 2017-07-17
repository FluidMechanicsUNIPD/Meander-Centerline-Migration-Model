!-----------------------------------------------------------------------
! LOCAL ERODIBILITY VALUE OF A FLOODPLAIN WITH GIVEN SYMMETRIC BUNDARIES
!-----------------------------------------------------------------------
! jshape : switch for the shape of the floodplain
!        (0:no boundaries, 1: straight boundaries, 2:bend floodplain, 3:specified shape)
! jbound : switch for the kind of boundary transition
!        (0:no boundary, 1:sharp, 2:spline, 3:exp-decaying, 4:exp-reversed)
! i    : index of the current point
! N    : number of points of axis river
! x, y : coordinates of the axis river
!-----------------------------------------------------------------------
       real*8 function boundaries(i)
       use module_global
       use module_geometry
       implicit none
       integer i
       real*8 c1, c2, c3, c4, Eaux, pi
       parameter ( pi = 3.1415927d0 )

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

           if ( ABS(y(i)).le.Lhalfvalley ) then
             boundaries = Ef
           else
             boundaries = Ebound
           end if

!-----------------------------------------------------------------------
! cubic spline
         case (2)
!-----------------------------------------------------------------------

           if ( ABS(y(i)).le.Lhalfvalley ) then
             boundaries = Ef
           else if ( ABS(y(i)).gt.(Lhalfvalley+Ltransvalley) ) then
             boundaries = Ebound
           else
             call cubecoeff(Lhalfvalley, Ef, 0.d0, &
     &              Lhalfvalley+Ltransvalley, Ebound, 0.d0, c1,c2,c3,c4)
             boundaries = c1 * (ABS(y(i)))**3.d0 +  &
     &                    c2 * (ABS(y(i)))**2.d0 +  &
     &                    c3 * ABS(y(i)) + &
     &                    c4
           end if

!-----------------------------------------------------------------------
! exponential decaying
         case (3)
!-----------------------------------------------------------------------

           if ( ABS(y(i)).le.Lhalfvalley ) then
             boundaries = Ef
           else
             boundaries = Ebound + (Ef-Ebound) * &
     &              EXP(-pi/Ltransvalley * (ABS(y(i))-Lhalfvalley))
           end if

!-----------------------------------------------------------------------
! reversed exponential
         case (4)
!-----------------------------------------------------------------------

           if ( ABS(y(i)).le.Lhalfvalley ) then
             boundaries = Ef
           else if ( ABS(y(i)).gt.(Lhalfvalley+Ltransvalley) ) then
             boundaries = Ebound
           else
             boundaries = MAX(Ebound, &
     &          1.d0+Ef-EXP(Ef/Ltransvalley*(ABS(y(i))-Lhalfvalley)) )
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
