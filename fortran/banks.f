c-----------------------------------------------------------------------
c RIVER BANKS FROM AXIS AND LOCAL WIDTH
c-----------------------------------------------------------------------
c Input:
c N   : numer of axis points
c x,y : axis coordinates
c b   : local width of the cross section
c Outputs
c xbr, ybr : coordinates of the first bank
c xbl, ybl : coordinates of the second bank
c-----------------------------------------------------------------------

      subroutine banks (N, x, y, b, xbr, ybr, xbl, ybl)
      implicit none
      integer i, j, k1, k2, N, nt
      parameter (nt=2)
      real*8 toll, dx, dy
      parameter ( toll = 1.D-09 )
      real*8 x(N), y(N), b(N), halfb(N), xbr(N), ybr(N), xbl(N), ybl(N)
      real*8 m(N), vnorm(nt), xplot(N, nt), yplot(N, nt)

c-----------------------------------------------------------------------
c loop over axis points
      do i = 1, N
c-----------------------------------------------------------------------

c set half-width
        halfb(i) = b(i) * 0.5d0

c-----------------------------------------------------------------------
c inverse slope of the river axis (i.e. slope of the transverse
c direction with respect to the longitudinal river axis)
c-----------------------------------------------------------------------

c first point
        if (i.eq.1) then
          k2 = i+1
          k1 = i
c last point
        else if (i.eq.N) then
          k2 = i
          k1 = i-1
c central points
        else
          k2 = i+1
          k1 = i-1
        end if

c coordinate differences
        dx = x(k2) - x(k1)
        dy = y(k2) - y(k1)

c local slope
        if ( ABS(dy).gt.toll ) then
          m(i) = - dx/dy
        else
          m(i) = - HUGE(1.d0) * SIGN(1.d0,dx) * SIGN(1.d0,dy)
        end if

c-----------------------------------------------------------------------
c loop over transverse direction
        do j = 1, nt
c-----------------------------------------------------------------------

c transverse coordinate
          vnorm(j) = -1.d0 + 2.d0 * DBLE(j-1)/DBLE(nt-1)

c-----------------------------------------------------------------------
c left and right banks
c-----------------------------------------------------------------------

c first upstream point
          if (i.eq.1) then
            k2 = i+1
            k1 = i

c last downstream point
          else if (i.eq.N) then
            k2 = i
            k1 = i-1

c points from second to second-last
          else if (i.gt.1) then
            k2 = i+1
            k1 = i-1
          end if

c horizontal segment
          if ( ABS(y(k2)-y(k1)) .lt. toll ) then
            xplot(i,j) = x(i)

c rightward segment
            if ( x(k2) .gt. x(k1) ) then
              yplot(i,j) = y(i) + vnorm(j)*halfb(i)

c leftward segment
            else
              yplot(i,j) = y(i) - vnorm(j)*halfb(i)

            end if

c upper-quadrant segment
          else if ( y(k2) .gt. y(k1) ) then

c vertical upward segment
            if ( ABS(x(k2)-x(k1)) .lt. toll ) then
              xplot(i,j) = x(i) - vnorm(j) * halfb(i)
              yplot(i,j) = y(i)

c rightward segment
            elseif ( x(k2) .gt. x(k1) ) then
              xplot(i,j) = x(i) -
     1          vnorm(j) * halfb(i) / (1+m(i)**2.d0)**0.5d0
              yplot(i,j) = y(i) +
     1          vnorm(j) * abs(m(i)) * halfb(i) / (1+m(i)**2.d0)**0.5d0

c leftward segment
            else if ( x(k2) .lt. x(k1) ) then
              xplot(i,j) = x(i) -
     2          vnorm(j) * halfb(i) / (1+m(i)**2.d0)**0.5d0
              yplot(i,j) = y(i) -
     2          vnorm(j) * abs(m(i)) * halfb(i) / (1+m(i)**2.d0)**0.5d0

            end if

c lower-quadrant segment
          else if ( y(k2) .lt. y(k1) ) then

c vertical downward segment
            if ( ABS(x(k2)-x(k1)) .lt. toll ) then
              xplot(i,j) = x(i) + vnorm(j) * halfb(i)
              yplot(i,j) = y(i)

c leftward segment
            else if ( x(k2) .lt. x(k1) ) then
              xplot(i,j) = x(i) +
     3          vnorm(j) * halfb(i) / (1+m(i)**2.d0)**0.5d0
              yplot(i,j) = y(i) -
     3          vnorm(j) * abs(m(i)) * halfb(i) / (1+m(i)**2.d0)**0.5d0

c rightward segment
            else if ( x(k2) .gt. x(k1) ) then
              xplot(i,j) = x(i) +
     4          vnorm(j) * halfb(i) / (1+m(i)**2.d0)**0.5d0
              yplot(i,j) = y(i) +
     4          vnorm(j) * abs(m(i)) * halfb(i) / (1+m(i)**2.d0)**0.5d0
            end if

c end if for segment orientation
          end if

c define right bank
          if (j.eq.1) then
            xbr(i) = xplot(i,j)
            ybr(i) = yplot(i,j)

c define left bank
          else if (j.eq.nt) then
            xbl(i) = xplot(i,j)
            ybl(i) = yplot(i,j)
          end if

c-----------------------------------------------------------------------
c end loop over transverse direction
        end do
c end loop over axis points
      end do
c-----------------------------------------------------------------------

c end of subroutine
      return
      end
