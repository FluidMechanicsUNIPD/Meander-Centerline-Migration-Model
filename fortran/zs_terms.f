c#######################################################################
c TERMS FOR ZS MODEL (Zolezzi and Seminara 2001)
c#######################################################################
      subroutine zs_terms(tau,phi,Cf,dCD,dCT,dphD,dphT,CD,CT,phiD,phiT)
      implicit none
      real*8 tau, phi, Cf, dCD, dCT, dphD, dphT
      real*8 CD, CT, phiD, phiT

c terms related to the friction coefficient
      CD = dCD / Cf
      CT = dCT * tau / Cf
      
c terms related to the sediment intensity
      if (phi.gt.0.d0) then
        phiD = dphD / phi
        phiT = dphT * tau / phi
      else
        phiD = 0.d0
        phiT = 0.d0
      end if
      
c end of subroutine
      return
      end
