c-----------------------------------------------------------------------
c SEMI-ANALYTIC COMPUTATION OF CONVOLUTION INTEGRALS
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c complex exponent > 0
c-----------------------------------------------------------------------
      complex*16 function SEMIANA1(js,N,c,jsend,deltas,lm)
      implicit none
      integer j, js, jsend, N
      real*8 deltas, c(N)
      complex*16 lm, aux1, lmds, lm2ds

      SEMIANA1 = dcmplx(0.D0)

      lmds = lm * deltas 
      lm2ds= lm * lmds
      aux1 = -2.D0 + zexp(lmds) + zexp(-lmds)

      j = js
      SEMIANA1 = c(j) / lm2ds * (zexp(-lmds)+lmds-1)    

      do j = js+1, jsend-1
        SEMIANA1 = SEMIANA1 + (c(j)*aux1/lm2ds) * zexp(lmds*dble(js-j))
      enddo
      
      j = jsend
      SEMIANA1 = SEMIANA1 - 
     6       c(j)/lm2ds*(1.D0+lmds-zexp(lmds))*zexp(lmds*dble(js-j))

      return
      end function SEMIANA1

c-----------------------------------------------------------------------
c complex exponent < 0
c-----------------------------------------------------------------------
      complex*16 function SEMIANA2(js,N,c,jsend,deltas,lm)
      implicit none
      integer j, js, jsend, N
      real*8 deltas, c(N)
      complex*16 lm, aux1, lmds, lm2ds

      SEMIANA2=dcmplx(0.D0)

      lmds = lm * deltas 
      lm2ds= lm * lmds
      aux1 = -2.D0 + zexp(lmds) + zexp(-lmds)

      j = jsend
      SEMIANA2 = c(j) / 
     6       lm2ds*(-1.D0+lmds+zexp(-lmds))*zexp(lmds*dble(js-j))   

      do j = jsend+1, js-1
        SEMIANA2 = SEMIANA2 + (c(j)*aux1/lm2ds) * zexp(lmds*dble(js-j))
      enddo
      
      j = js
      SEMIANA2 = SEMIANA2 - c(j) / lm2ds * (1.D0+lmds-zexp(lmds))

      return
      end function SEMIANA2
