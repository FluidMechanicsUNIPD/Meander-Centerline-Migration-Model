!-----------------------------------------------------------------------
! POINT INSIDE POLYGON WITH WINDING NUMBER PROCEDURE
! [Hormann & Agathos 2001]
!-----------------------------------------------------------------------
! Given a loop polyline C(x,y) made by N points and the point R(xR,yR), 
! the question is whether the point lies inside or outside the polygon. 
! The winding number wn is the number of revolution made by C around P:
! if P is inside C then W is nonzero.
!-----------------------------------------------------------------------
      integer function windingnumber(xR, yR, N, x, y)
      implicit none
      integer i, N, wn, dq, q(N)
      real*8 det, xp, xp1, yp, yp1
      real*8 xR, yR, x(N), y(N)
      
! loop over the polygon nodes
      do i = 1, N-1   

        if ((x(i).gt.xR).and.(y(i).ge.yR)) then 
          q(i) = 0
        else if ((x(i).le.xR).and.(y(i).gt.yR)) then
          q(i) = 1
        else if ((x(i).lt.xR).and.(y(i).le.yR)) then
          q(i) = 2
        else if ((x(i).ge.xR).and.(y(i).lt.yR)) then
          q(i) = 3
        end if

      end do
      q(N) = q(1)

! initialize wn
      wn = 0  

! loop for the polygon nodes
      do i = 1, N-1  
        dq = q(i+1) - q(i)

        if ((dq.eq.(1)).or.(dq.eq.(-3))) then  
          wn = wn + 1
        else if ((dq.eq.(-1)).or.(dq.eq.(3))) then
          wn = wn - 1
        else if ((dq.eq.(-2)).or.(dq.eq.(2))) then

! evaluate orientation of triangle (R, P, P1)
          xp = x(i)
          xp1 = x(i+1)
          yp = y(i)
          yp1 = y(i+1) 
          det = (xp-xR)*(yp1-yR) - (xp1-xR)*(yp-yR)

          if (det.gt.0.d0) then
            wn = wn + 2
          else if (det.lt.0.d0) then
            wn = wn - 2
          end if
        end if

      end do

      windingnumber = wn

!-----------------------------------------------------------------------
! end of subroutine
!-----------------------------------------------------------------------
      return
      end
