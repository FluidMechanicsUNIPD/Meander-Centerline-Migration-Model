c-----------------------------------------------------------------------
c Calcolo radici polinomio di grado m a(i)x**(i-1) [i=1:m+1]
c-----------------------------------------------------------------------
c Numerical Recipes in Fortran 77 Ch. 9.5, p. 366-367
c Given the degree m and the coefficient a(1:m+1) of the polinomial 
c a(i)x**(i-1)  this routine calls laguer and finds all m complex roots
c the logical variable polish should be input as .true. if polishing is 
c desired, .false. if the roots will be not polished
c-----------------------------------------------------------------------
      subroutine zroots(a,roots)
      implicit none
      logical polish
      integer m,maxm,i,j,jj,its
      real*8 eps
      parameter (m=4, polish=.true.)
      parameter (eps=dble(1.e-8), maxm=101) !Desired accuracy and max anticipated value of m+1
      complex*16 a(m+1), roots(m)
      complex*16 ad(maxm),x,b,c

!copy of coefficient for successive deflation
      do j = 1, m+1                      
        ad(j)=a(j)
      end do
        
!loop over each root to be found
      do j = m, 1, -1                     
        x = dcmplx(1.D0,0.D0)

!find the root
        call laguer(ad,j,x,its)
        if( abs(dimag(x)) .le. (2.D0*eps**2*abs(dreal(x))) )
     &            x = dcmplx(dreal(x),0.D0)
        roots(j) = x

!forward deflation
        b = ad(j+1)

        do jj = j, 1, -1
          c = ad(jj)
          ad(jj) = b
          b = x*b + c
        end do
      end do

!polish the roots using the undeflated coefficients
      if (polish) then
        do j = 1, m
          call laguer(a,m,roots(j),its)
          if(dabs(dimag(roots(j))).lt.eps) then
            roots(j) = dcmplx(dreal(roots(j)),0.D0)
          end if
        end do
      end if

!      call ordina(roots, m)
      do i = 1, m-1
        do j = i+1, m
          if (DBLE(roots(j)).gt.DBLE(roots(i))) then
            x = roots(i)
            roots(i) = roots(j)
            roots(j) = x
          end if
        end do
      end do

      return
      end

c-----------------------------------------------------------------------
c Sort the four complex roots of the 4th order polynomial
c-----------------------------------------------------------------------
      subroutine ordina(roots, m)
      implicit none
      integer m, i, j
      real*8 toll
      parameter (toll=dble(1.e-6))
      complex*16 roots(m), x

      do i = 1,m
        if ( dreal(roots(i)).gt.dreal(roots(1)) ) then
          if ( ABS(imag(roots(i))).lt.toll ) then
            x = roots(1)
            roots(1) = roots(i)
            roots(i) = x
          end if
        end if
      end do
        
      do i = 2,m
        if ( dreal(roots(i)).lt.dreal(roots(m)) ) then
          if ( ABS(imag(roots(i))).lt.toll ) then
            x = roots(m)
            roots(m) = roots(i)
            roots(i) = x       
          end if
        end if
      end do         

      return
      end             

c-----------------------------------------------------------------------
        subroutine laguer(a,m,x,its)
c-----------------------------------------------------------------------
      implicit none
        integer m, its, maxit, mr, mt
        real*8 eps,epss
        complex*16 a(m+1), x
        parameter(epss=dble(2.e-7), mr=8, mt=10, maxit=mt*mr)
        integer iter, j
        real*8 abx, abp, abm, err, frac(mr)
        complex*16 dx, x1, b, d, f, g, h, sq, gp, gm, g2
        save frac
        data frac/.5d0,.25d0,.75d0,.13d0,.38d0,.62d0,.88d0,1.d0/
        do iter=1,maxit
                 its=iter
                 b=a(m+1)
                 err=cdabs(b)
                 d=dcmplx(0.d0,0.d0)
                 f=dcmplx(0.d0,0.d0)
                 abx=cdabs(x)
                 do j=m,1,-1
                          f=f*x+d
                          d=x*d+b
                          b=x*b+a(j)
                          err=cdabs(b)+abx*err
                 end do
                 err=epss*err
                 if(cdabs(b).le.err) then
                          return
                 else
                          g=d/b
                          g2=g*g
                          h=g2-2.d0*f/b
                          sq=sqrt((m-1)*(m*h-g2))
                          gp=g+sq
                          gm=g-sq
                          abp=cdabs(gp)
                          abm=cdabs(gm)
                          if(abp.lt.abm) gp=gm
                          if(max(abp,abm).gt.0.D0) then
                                        dx=m/gp
                          else
                             dx=zexp(dcmplx(dlog(1.d0+abx),dble(iter)))
                          end if
                 end if
                 x1=x-dx
                 if(x.eq.x1) return
                 if(mod(iter,mt).ne.0) then
                          x=x1
                 else
                          x=x-dx*frac(iter/mt)
                 end if
        end do
        stop 'too many interactions in laguer'
        return
        end
