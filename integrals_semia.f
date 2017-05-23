c-----------------------------------------------------------------------
c SUBROUTINES PER IL CALCOLO SEMIANALITICO DI INTEGRALI
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c Integrazione lm>0; esponenziale reale
c-----------------------------------------------------------------------
      subroutine SEMIANA1(js,N,c,jsend,deltas,lm,x)
      implicit none
      integer j, js, jsend, N
      real*8 lm, lmds, lm2ds, aux1, deltas, c(N)
      complex*16 x

      x=dcmplx(0.D0)

      lmds = lm * deltas 
      lm2ds= lm * lmds
      aux1 = -2.D0 + dexp(lmds) + dexp(-lmds)

      j = js
      x = c(j) / lm2ds * (dexp(-lmds)+lmds-1)    

      do j = js+1, jsend-1
        x = x + ( c(j)*aux1/lm2ds) * dexp(lmds*dble(js-j) )
      enddo
      
      j = jsend
      x = x-c(j)/lm2ds*(1.D0+lmds-dexp(lmds))*dexp(lmds*dble(js-j))

      return
      end

c-----------------------------------------------------------------------
c Integrazione lm>0; esponenziale complesso
c-----------------------------------------------------------------------
      subroutine SEMIANA2(js,N,c,jsend,deltas,lm,x)
      implicit none
      integer j, js, jsend, N
      real*8 deltas, c(N)
      complex*16 x, lm, aux1, lmds, lm2ds

      x = dcmplx(0.D0)

      lmds = lm * deltas 
      lm2ds= lm * lmds
      aux1 = -2.D0 + zexp(lmds) + zexp(-lmds)

      j = js
      x = c(j) / lm2ds * (zexp(-lmds)+lmds-1)    

      do j = js+1, jsend-1
        x = x + (c(j)*aux1/lm2ds) * zexp(lmds*dble(js-j))
      enddo
      
      j = jsend
      x = x-c(j)/lm2ds*(1.D0+lmds-zexp(lmds))*zexp(lmds*dble(js-j))

      return
      end

c-----------------------------------------------------------------------
c Integrazione lm<0; esponenziale reale
c-----------------------------------------------------------------------
      subroutine SEMIANA3(js,N,c,jsend,deltas,lm,x)
      implicit none
      integer j, js, jsend, N
      real*8 lm, lmds, lm2ds, aux1, deltas, c(N)
      complex*16 x

      x = dcmplx(0.D0)

      lmds = lm * deltas 
      lm2ds= lm * lmds
      aux1 = -2.D0 + dexp(lmds) + dexp(-lmds)

      j = jsend
      x = c(j)/lm2ds*(-1.D0+lmds+dexp(-lmds))*dexp(lmds*dble(js-j))   

      do j = jsend+1, js-1
        x = x + (c(j)*aux1/lm2ds) * dexp(lmds*dble(js-j))
      enddo
      
      j = js
      x = x - c(j) / lm2ds * (1.D0+lmds-dexp(lmds))

      return
      end

c-----------------------------------------------------------------------
c Integrazione lm<0; esponenziale complesso
c-----------------------------------------------------------------------
      subroutine SEMIANA4(js,N,c,jsend,deltas,lm,x)
      implicit none
      integer j, js, jsend, N
      real*8 deltas, c(N)
      complex*16 x, lm, aux1, lmds, lm2ds

      x=dcmplx(0.D0)

      lmds = lm * deltas 
      lm2ds= lm * lmds
      aux1 = -2.D0 + zexp(lmds) + zexp(-lmds)

      j = jsend
      x = c(j) / lm2ds*(-1.D0+lmds+zexp(-lmds))*zexp(lmds*dble(js-j))   

      do j = jsend+1, js-1
        x = x + (c(j)*aux1/lm2ds) * zexp(lmds*dble(js-j))
      enddo
      
      j = js
      x = x - c(j) / lm2ds * (1.D0+lmds-zexp(lmds))

      return
      end
