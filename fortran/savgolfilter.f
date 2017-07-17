      subroutine savgolfilter(n, c, cfilt)
      implicit none
      integer i, j, jn, n, nr, nl, np ,m, ld
      parameter (nr=5, nl=5, np=nl+nr+1, m=2, ld=0)
      integer indexSAVGOL(np), indexPOS(N)
      real*8 c(n), aux(n), cfilt(n), coeff(np)

      call savgol(coeff,np,nl,nr,ld,m)
      aux(:) = c(:)
      
! initial points
      do j = 1, 1+nl
      
! weighting the current point
        aux(j) = coeff(1) * c(j)

! weighting through the next points (right points)
        do jn = 1, nr
          aux(j) = aux(j) + coeff(np-jn+1) * c(j+jn)
        end do
        
      end do
      
! central points     
      do j = 1+nl , N-nr
      
! weighting the current point
        aux(j) = coeff(1) * c(j)

! weighting through the previous points (left points)
        do jn = 1, nl
          aux(j) = aux(j) + coeff(jn+1) * c(j-jn)
        end do

! weighting through the next points (right points)
        do jn = 1, nr
          aux(j) = aux(j) + coeff(np-jn+1) * c(j+jn)
        end do
        
      end do
      
! end points
      do j = N-nr, N
      
! weighting the current point
        aux(j) = coeff(1) * c(j)

! weighting through the previous points (left points)
        do jn = 1, nl
          aux(j) = aux(j) + coeff(jn+1) * c(j-jn)
        end do
        
      end do
      
      cfilt(:) = aux(:)
      
c-----------------------------------------------------------------------
      return
999     format(i3,1x,*(f10.6,1x))
      end
      
!-------------------------------------------------------------------------------------------- 
!USES lubksb,ludcmp given below. 
!Returns in c(1:np), in wrap-around order (see reference) consistent with the argument respns 
!in routine convlv, a set of Savitzky-Golay filter coefficients. nl is the number of leftward 
!(past) data points used, while nr is the number of rightward (future) data points, making 
!the total number of data points used nl +nr+1. ld is the order of the derivative desired 
!(e.g., ld = 0 for smoothed function). m is the order of the smoothing polynomial, also 
!equal to the highest conserved moment; usual values are m = 2 or m = 4. 
!--------------------------------------------------------------------------------------------
      SUBROUTINE savgol(c,np,nl,nr,ld,m)
      INTEGER ld,m,nl,np,nr,MMAX
      REAL*8 c(np)
      PARAMETER (MMAX=6)
CU    USES lubksb,ludcmp
      INTEGER imj,ipj,j,k,kk,mm,indx(MMAX+1)
      REAL*8 d,fac,sum,a(MMAX+1,MMAX+1),b(MMAX+1)
      if(np.lt.nl+nr+
     *1.or.nl.lt.0.or.nr.lt.0.or.ld.gt.m.or.m.gt.MMAX.or.nl+nr.lt.m)then
	  write(6,*) 'bad args in savgol'
	  stop
      end if  
      do 14 ipj=0,2*m
        sum=0.D0
        if(ipj.eq.0)sum=1.D0
        do 11 k=1,nr
          sum=sum+dble(k)**ipj
11      continue
        do 12 k=1,nl
          sum=sum+dble(-k)**ipj

12      continue
        mm=min(ipj,2*m-ipj)
        do 13 imj=-mm,mm,2
          a(1+(ipj+imj)/2,1+(ipj-imj)/2)=sum
13      continue
14    continue
      call ludcmp(a,m+1,MMAX+1,indx,d)
      do 15 j=1,m+1
        b(j)=0.D0
15    continue
      b(ld+1)=1.D0
      call lubksb(a,m+1,MMAX+1,indx,b)
      do 16 kk=1,np
        c(kk)=0.D0
16    continue
      do 18 k=-nl,nr
        sum=b(1)
        fac=1.D0
        do 17 mm=1,m
          fac=fac*k
          sum=sum+b(mm+1)*fac
17      continue
        kk=mod(np-k,np)+1
        c(kk)=sum
18    continue
      return
      end

!***************************************************************
!* Given an N x N matrix A, this routine replaces it by the LU *
!* decomposition of a rowwise permutation of itself. A and N   *
!* are input. INDX is an output vector which records the row   *
!* permutation effected by the partial pivoting; D is output   *
!* as -1 or 1, depending on whether the number of row inter-   *
!* changes was even or odd, respectively. This routine is used *
!* in combination with LUBKSB to solve linear equations or to  *
!* invert a matrix. Return code is 1, if matrix is singular.   *
!***************************************************************
      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      REAL*8 d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.e-20)
      INTEGER i,imax,j,k
      REAL*8 aamax,dum,sum,vv(NMAX)
      d=1.D0
      do 12 i=1,n
        aamax=0.D0
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.D0) pause 'singular matrix in ludcmp'
        vv(i)=1.D0/aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.D0
        do 16 i=j,n

          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.D0)a(j,j)=TINY
        if(j.ne.n)then
          dum=1.D0/a(j,j)

          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      end

!******************************************************************
!* Solves the set of N linear equations A . X = B.  Here A is     *
!* input, not as the matrix A but rather as its LU decomposition, *
!* determined by the routine LUDCMP. INDX is input as the permuta-*
!* tion vector returned by LUDCMP. B is input as the right-hand   *
!* side vector B, and returns with the solution vector X. A, N and*
!* INDX are not modified by this routine and can be used for suc- *
!* cessive calls with different right-hand sides. This routine is *
!* also efficient for plain matrix inversion.                     *
!******************************************************************
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL*8 sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.D0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      end
