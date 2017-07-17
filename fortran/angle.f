c-----------------------------------------------------------------------
c ANGLE THETA OF THE RIVER AXIS WITH RESPECT TO THE CARTESIAN REFERENCE
c-----------------------------------------------------------------------
c Given a plane curve specified by its cartesian coordinate (x,y), 
c calculates the angle theta between the local tangent and the x-axis
c Setting th(j)=-datan2(dy,dx) implies that :
c     c(js) = ( th(js+1)-th(js-1) )/2*ds
c     xnew(js) = x(js) + dt * deltaU(js) * dsin(th(js))
c     ynew(js) = y(js) + dt * deltaU(js) * dcos(th(js))
c-----------------------------------------------------------------------
      subroutine angle(N, x, y, th)
      implicit none
      integer j, N
      real*8 dx, dy, a1, a2, pi
      real*8 aux0, aux1
      real*8 th(N), x(N), y(N)
      parameter (pi = DATAN(1.d0)*4.d0)	! pi = 3.1415
     
c first point
      j = 1
      dy = y(j+1)-y(j)
      dx = x(j+1)-x(j)
      th(j) = -DATAN2(dy,dx)
      aux0 = 0.5d0 * th(j)
      
c central points      
      do j = 2, N-1
        dy = y(j+1) - y(j)
        dx = x(j+1) - x(j)
        aux1 = -0.5D0 * DATAN2(dy,dx)
        th(j) = aux1 + aux0
        aux0 = aux1
      end do

c last point      
      j = N
      dy = y(j)-y(j-1)
      dx = x(j)-x(j-1)
      th(j) = -DATAN2(dy,dx)

c set critical angles
      a1 = pi/2.0d0
      a2 = pi*3.d0/2.d0

c check if function DATAN2 provided sharp gaps crossing pi angles
      j = 1
      if ((dabs(th(j+1)-th(j)).gt.a1).and.
     &                (dabs(th(j+1)-th(j)).lt.a2)) then    
        th(j) = th(j) - pi * dsign(1.d0,(th(j+1)-th(j)))
       else if (dabs(th(j+1)-th(j)).gt.pi) then
        th(j) = th(j) - 2.d0 * pi * dsign(1.d0,(th(j+1)-th(j)))
      end if

      do j = 2, N
        if ((dabs(th(j)-th(j-1)).gt.a1).and.
     &                (dabs(th(j)-th(j-1)).lt.a2)) then    
          th(j) = th(j) - pi * dsign(1.d0,(th(j)-th(j-1)))
         else if (dabs(th(j)-th(j-1)).gt.pi) then
          th(j) = th(j) - 2.d0 * pi * dsign(1.d0,(th(j)-th(j-1)))
        end if
      end do

c-----------------------------------------------------------------------
c end of the subroutine
      return
      end
