c-----------------------------------------------------------------------
c DIMENSIONLESS CURVATURE OF THE RIVER AXIS
c-----------------------------------------------------------------------
c The dimensionless curvature c of river axis corresponds to the 
c negative derivative of the axis angle theta with respect to the
c curvilinear coordinate, i.e. c = - nu0 * d(th)/ds
c NB angle is already computed as -th
c-----------------------------------------------------------------------
      subroutine curv(N, s, th, c, ksavgol, jt)
      implicit none
      integer ksavgol, N, j, jt, np
      parameter (np = 11)
      real*8 s(N), th(N), c(N), caux(N)

c curvature
      call derivative(N, s, th, c, 3)
      
c numerical smoothing 
      call smoothing(N, c)

c repeated Savitzsky-Golay filtering at each iterations
      if (ksavgol.lt.0) then
        do j = 1, ABS(ksavgol)
          call savgolfilter_new (np, N, c, caux)
          c(:) = caux(:)
        end do

c periodic Savitzsky-Golay filtering
      else if (ksavgol.gt.0) then
        if (MOD(jt,ksavgol).eq.0) then
          call savgolfilter_new (np, N, c, caux)
          c(:) = caux(:)
        end if
      end if
     
c-----------------------------------------------------------------------
c end of the subroutine
      return
      end
