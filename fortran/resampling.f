c-----------------------------------------------------------------------
c RESAMPLE THE RIVER AXIS THROUGH A CUBIC SPLINE
c-----------------------------------------------------------------------
      subroutine resampling
      use module_global
      use module_geometry
      implicit none

      integer j, flag, Nnew, ierr
      real*8 dist
      real*8, allocatable :: saux(:), xaux(:), yaux(:)

c initialize flag
      flag = 0
      
c initialize longitudinal coordinate of axis
      s(1) = 0.d0

c loop over points, from second to last
      do j = 2, N
      
c evaluate distance between two consecutive points
        dist = DSQRT((x(j)-x(j-1))**2.d0 + (y(j)-y(j-1))**2.d0)

c check whether distance is into the allowed range
        if ((dist.gt.(deltas0*dsmax)).or.(dist.lt.(deltas0*dsmin))) then
          flag = flag + 1
        end if
        
c update coordinate
        s(j) = s(j-1) + dist

c end loop over poitns
      end do

c mean distance
      deltas = s(N) / DBLE(N-1)

c if flag is not zero, then resample points
      if ((flag.gt.0).or.(Nnco.ne.Nncoold)) then

c update resampling counter
        Nsamp = Nsamp + 1
        
c update point number
        Nnew = INT(s(N)/deltas0) + 1
        deltas = s(N) / DBLE(Nnew-1)
   
c initialize auxiliary arrays
        allocate ( saux(Nnew), xaux(Nnew), yaux(Nnew), stat = ierr )
        if(ierr.ne.0) go to 999

c new longitudinal coordinate 
        saux(1) = 0.d0
        do j = 2, Nnew
          saux(j) = saux(j-1) + deltas
        end do

c resampling x- and y-coordinates
        call Cspline(s,x,N,saux,xaux,Nnew)
        call Cspline(s,y,N,saux,yaux,Nnew)

c restore pristine arrays
        deallocate (s, x, y)
        allocate ( s(Nnew), x(Nnew), y(Nnew), stat = ierr )
        if(ierr.gt.0) go to 999

        s(:) = saux(:)
        x(:) = xaux(:)
        y(:) = yaux(:)

        deallocate (saux, xaux, yaux)

c update point number
        N = Nnew
        
c end if for point resampling
      end if

c-----------------------------------------------------------------------
c end of subroutine
      return
999   write(6,*) 'ERROR! Allocation in point resampling'
      stop
      end

c-----------------------------------------------------------------------
c Given the arrays xa(1:N) and ya(1:N), containing a tabulated function 
c (i.e., yi=f(xi)), with x1<x2<...<xn, given the array y2a(1:N), which 
c is the output of spline, this routine returns the cubic spline 
c interpolated values y(1:NN) for x(1:NN)
c-----------------------------------------------------------------------
      subroutine Cspline(xa,ya,N,xx,yy,NN)
      implicit none
      integer N, NN, j
      real*8 xa(N), ya(N), y2a(N), xx(NN), yy(NN)
      real*8 x, y
c-----------------------------------------------------------------------
      call spline(xa,ya,N,y2a)

      yy(1) = ya(1)
      do j = 2, NN-1
        x = xx(j)
        call splint(xa,ya,y2a,N,x,y)
        yy(j) = y
      end do
      yy(NN) = ya(N)
      return
      end

c-----------------------------------------------------------------------
c Given the arrays xa(1:N) and ya(1:N), containing a tabulated function 
c (i.e., yi=f(xi)), with x1<x2<...<xn, given the array y2a(1:N), which 
c is the output of spline, and given a value of x, this routine returns 
c the cubic spline interpolated value y
c-----------------------------------------------------------------------
      subroutine splint(xa,ya,y2a,N,x,y)
      implicit none
      integer N, k, khi, klo
      real*8 x, y, a, b, h
      real*8 xa(N), y2a(N), ya(N)
c-----------------------------------------------------------------------
      klo = 1
      khi = N
      do while ((khi-klo).gt.1)
        k = INT( DBLE((khi+klo)/2.d0) )
        if(xa(k).gt.x) then
          khi = k
        else
          klo = k
        end if
      end do
        
      h = xa(khi) - xa(klo)
      if (h.eq.0.d0) then
        write(6,*) 'khi= ',khi, '  klo = ', klo, '  N = ', N
        stop 'ERROR! Bad xa in splint'
      end if
      
      a = (xa(khi)-x) / h
      b = (x-xa(klo)) / h
      y = a * ya(klo) + b * ya(khi) +
     5  ( (a**3.d0-a)*y2a(klo) + (b**3.d0-b)*y2a(khi))*(h**2.d0)/6.d0
      return
      end

c-----------------------------------------------------------------------
c Given arrays x(1:N) and y(1:N) containing a tabulated function 
c (i.e., yi=f(xi)), with x1<x2<...<xn, returns an array y2(1:N) of 
c length N which contains the second derivatives of the interpolating 
c function at the tabulated points. The boundary conditions for a 
c natural spline are adopted (i.e., a zero second derivative is assumed 
c at the boundaries) Nmax: largest anticipated value of N
c-----------------------------------------------------------------------
        subroutine spline(x,y,N,y2)
        implicit none
        integer i, k, N, Nmax
        real*8 p, qn, sig, un
        real*8 x(N), y(N), y2(N), u(N)
c-----------------------------------------------------------------------
c lower boundary condition is set to be natural
        y2(1) = 0.d0
        u(1) = 0.d0

c decomposition loop of the tridiagonal algorithm
c y2 and u are used for temporary storage of decomposed factors
        do i = 2, N-1
          sig = (x(i)-x(i-1)) / (x(i+1)-x(i-1))
          p = sig*y2(i-1) + 2.d0
          y2(i) = (sig - 1.d0) / p
          u(i) = (6.d0*((y(i+1)-y(i))/(x(i+1)-x(i)) -(y(i)-y(i-1))
     &          /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
        end do

c upper boundary condition is set to be natural
        qn = 0.d0
        un = 0.d0
        y2(N) = (un-qn*u(N-1))/(qn*y2(N-1) + 1.d0)

c backsostitution loop of the tridiagonal algorithm
        do k = N-1, 1, -1
          y2(k) = y2(k) * y2(k+1) + u(k)
        end do

        return
        end
