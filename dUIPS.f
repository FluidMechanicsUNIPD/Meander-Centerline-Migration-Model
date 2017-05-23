c-----------------------------------------------------------------------
c EXCESS NEAR-BANK LONGITUDINAL VELOCITY
c Ikeda et al. (1981), with free boundary conditions
c-----------------------------------------------------------------------
      subroutine dUIPS
      use module_global
      use module_geometry
      implicit none

      integer js, jc
      real*8 A, alfa, conv, lambda
      
c characteristic exponent
      lambda = -2.d0*beta*Cf

c scour factor
      alfa = 3.d0
      A = 1.d0 + alfa

c-----------------------------------------------------------------------
c loop over points
      do js = 1, N
c-----------------------------------------------------------------------
      
c convolution term
        conv = 0.d0
        do jc = 1, js
          conv = conv + c(jc) * DEXP(lambda * ( s(js)-s(jc) )) * deltas
        end do

c excess near-bank velocity
        dU(js) = 2.d0 * (- c(js) + beta * Cf * (A+F0**2.d0) * conv )

c-----------------------------------------------------------------------        
c end loop over points
      end do
c-----------------------------------------------------------------------

c end of subroutine
      return
      end
