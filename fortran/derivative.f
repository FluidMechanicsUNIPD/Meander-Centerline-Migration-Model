c-----------------------------------------------------------------------
c DISCRETE DERIVATIVE OF A FUNCTION (x,y) OF N POINTS
c-----------------------------------------------------------------------
c Flag is a switch for the numerical differentiation method:
c 1 = derivative in a point is computed by averaging the slopes of the 
c     two segments sharing that points
c 2 = derivative is computed with second-order central finite difference
c 3 = derivative is computed with quasi fourth-order central finite difference
c 4 = derivative is computed with fourth-order central finite difference
c     (last solutions requires N >=6 )
c NB requires the function slope (reported below), to check eventual 
c denominator approaching zero
c-----------------------------------------------------------------------
      subroutine derivative(N, x, y, dydx, flag)
      implicit none
      integer i, N, flag
      real*8 dy, dx, aux1, aux2, slope
      real*8 x(N), y(N), dydx(N)

c-----------------------------------------------------------------------
c BACK-AND-FORTH SLOPE WEIGHTING
      if (flag.eq.1) then
c-----------------------------------------------------------------------

c forward finite difference for first point
        i = 1
        dy = y(i+1) - y(i)
        dx = x(i+1) - x(i)
        aux1 = slope(dx,dy)
        dydx(i) = aux1

c back and forth slope averaging for central points
        do i = 2, N-1
          dy = y(i+1) - y(i)
          dx = x(i+1) - x(i)
          aux2 = slope(dx,dy)

          dydx(i) = 0.5d0 * (aux1 + aux2)
          aux1 = aux2
        end do

c backward finite difference for first point
        i = N
        dy = y(i) - y(i-1)
        dx = x(i) - x(i-1)
        dydx(i) = slope(dx,dy)

c-----------------------------------------------------------------------
c SECOND-ORDER CENTRAL FINITE DIFFERENCE
      else if (flag.eq.2) then
c-----------------------------------------------------------------------

c forward finite difference for first point
        i = 1
        dy = y(i+1) - y(i)
        dx = x(i+1) - x(i)
        dydx(i) = slope(dx,dy)

c central finite difference for central points
        do i = 2, N-1
          dy = y(i+1)-y(i-1)
          dx = x(i+1)-x(i-1)
          dydx(i) = slope(dx,dy)
        end do

c backward finite difference for last point
        i = N
        dy = y(i) - y(i-1)
        dx = x(i) - x(i-1)
        dydx(i) = slope(dx,dy)

c-----------------------------------------------------------------------
c QUASI FOURTH-ORDER FINITE DIFFERENCE
      else if (flag.eq.3) then
c-----------------------------------------------------------------------

c first point
        i = 1
        dy = y(i+1) - y(i)
        dx = x(i+1) - x(i)
        dydx(i) = slope(dx,dy)

c second point
        i = 2
        dy = y(i+1) - y(i-1)
        dx = x(i+1) - x(i-1)
        dydx(i) = slope(dx,dy)
      
c central finite difference for central points
        do i = 3, N-2
          dy = -1.d0/12.d0*y(i+2) + 2.d0/3.d0*y(i+1) 
     2         -2.d0/3.d0*y(i-1) + 1.d0/12.d0*y(i-2)
          dx = (x(i+1)-x(i-1))/2.d0
          dydx(i) = slope(dx,dy)
        end do

c second to last point
        i = N-1
        dy = y(i+1) - y(i-1)
        dx = x(i+1) - x(i-1)
        dydx(i) = slope(dx,dy)

c last point
        i = N
        dy = y(i) - y(i-1)
        dx = x(i) - x(i-1)
        dydx(i) = slope(dx,dy)

c-----------------------------------------------------------------------
c FOURTH-ORDER FINITE DIFFERENCE
      else if (flag.eq.4) then
c-----------------------------------------------------------------------

c forward finite difference for first and second points
        do i = 1, 2
          dy = - 3.d0/2.d0*y(i+2) + 2.d0*y(i+1) - 1.d0/2.d0*y(i)
          dx = (x(i+1)-x(i))
          dydx(i) = slope(dx,dy)
        end do

c central finite difference for central points
        do i = 3, N-2
          dy = -1.d0/12.d0*y(i+2) + 2.d0/3.d0*y(i+1)
     2         -2.d0/3.d0*y(i-1) + 1.d0/12.d0*y(i-2)
          dx = (x(i+1)-x(i-1))/2.d0
          dydx(i) = slope(dx,dy)
        end do

c backward finite difference for second-last and last points
        do i = N-1,N
          dy = 1.d0/2.d0*y(i-2) - 2.d0*y(i-1) + 3.d0/2.d0*y(i)
          dx = (x(i)-x(i-1))/2.d0
          dydx(i) = slope(dx,dy)
        end do

      else
        stop 'ERROR! Wrong flag for derivative subroutine'
      end if

c-----------------------------------------------------------------------
c end fo subroutine
c-----------------------------------------------------------------------
      return
      end

c-----------------------------------------------------------------------
c SLOPE CHECKING FOR EVENTUAL ZERO DENOMINATOR
c-----------------------------------------------------------------------
      real*8 function slope(dx,dy)
      implicit none
      real*8 dx, dy, toll
      parameter (toll=1.d-06)

c if the denominator approaches zero, then set the ratio equal to the
c largest available finite number
      if (DABS(dx).lt.toll) then
        slope = HUGE(dy) * SIGN(1.d0,dy) * SIGN(1.d0,dx)
      else
        slope = dy/dx
      end if

      return
      end
