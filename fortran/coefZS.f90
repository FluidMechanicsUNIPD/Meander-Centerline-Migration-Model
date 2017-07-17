!-----------------------------------------------------------------------
! COEFFICIENTI PER IL CAMPO DI MOTO A CURVATURA VARIABILE PER CONDIZIONI
! AL CONTORNO PERIODICHE (Zolezzi e Seminara 2001)
!-----------------------------------------------------------------------
      subroutine coefZS(jt)
      
      use module_global
      use module_zs
      implicit none

      integer jm, jt

      real*8 k0, k1, k3, k4
      real*8 a1, a2, a3, a4, a5, a6
      real*8 b1, b2, b3, b4, b5, b6
      real*8 et3, et4, et5, et6
      real*8 T0, T1, T2, T3, Tc1, Tc2, Tc3, Tc4, Tc5
      real*8 s1, s2, f1, f2, M
      real*8 h1bar, h2bar, h3bar
      real*8 d1bar, d2bar, d3bar
      real*8 alf0, alf1, alf2
      real*8 bet2, bet3, bet4
      real*8 csi1, csi2, csi3, csi4
      real*8 del0, del1, del2
      real*8 Delta, Delta0, Delta1, Delta2
      real*8 eps3, eps4
      real*8 gam2, gam3, gam4, gam5, gam6
      real*8 mu1, mu2, mu3, mu4, mu5, mu6
      real*8 rho0, rho1, rho2, rho3, rho4, rho5, rho6
      real*8 sigma0, sigma1, sigma2, sigma3
      real*8 pi
      parameter (pi = DATAN(1.d0)*4.d0)	! pi = 3.1415

      complex*16 lam1, lam2, lam3, lam4
      complex*16 waux1, waux2, waux3
      complex*16 wdet, W1det, W2det, W3det, W4det
      complex*16 coef(5), roots(4), wj(4)

!-----------------------------------------------------------------------
! coefficients k0, k1, k2, k3, k4
!-----------------------------------------------------------------------
	  call k0134(Cf,Nz,k0,k1,k3,k4)

!-----------------------------------------------------------------------
! ZS coefficients
!-----------------------------------------------------------------------

      s1 = 2.d0 / (1.d0-CT)
      s2 = CD / (1.d0-CT)
      f1 = 2.d0 * phiT / (1.d0-CT)
      f2 = phiD + CD * phiT / (1.d0-CT)

! coefficients a1 - a6
      a1 = beta * Cf * s1
      a2 = beta * Cf * (s2 - 1.d0)
      a3 = beta * Cf
      a4 = f1
      a5 = f2
      a6 = rpic / (beta*theta0**0.5d0)

! coefficients b1 - b6
      b1 = - beta*Cf
      b2 = 1.d0 - dsqrt(Cf) * k3
      b3 = -k0 / (beta*Cf**0.5d0) - k4/beta
      b4 = k3 * theta0**0.5d0 / (rpic*Cf**0.5d0)
      b5 = - k1 / (Cf*beta**2.d0)
      b6 = k4 * theta0**0.5d0 / (beta*Cf*rpic)

! coefficients hbar and dbar
      h1bar = b2
      h2bar = b3
      h3bar = b5
      d1bar = (F0**2.d0) * h1bar - b4
      d2bar = (F0**2.d0) * h2bar - b6
      d3bar = (F0**2.d0) * h3bar

      alf0 = a2
      alf1 = 1.d0 / (F0**2.d0)

      bet2 = a1
      bet3 = 1.d0

      gam2 = b1 - a2 * d1bar
      gam3 = -h1bar - a2 * d2bar

      del1 = a5 - 1.d0 - (F0**2.d0) * a3 * a6;
      del2 = - (F0**2.d0) * a6;

      eps3 = a4 - 1.d0 - (F0**2.d0) * a3 * a6
      eps4 = del2

      et3 = - del1 * d1bar
      et4 = - del1 * d2bar + (F0**2.d0) * a6 * d1bar
      et5 = - del1 * d3bar + (F0**2.d0) * a6 * d2bar
      et6 = (F0**2.d0) * a6 * d3bar

!-----------------------------------------------------------------------
      do jm = 1, Mdat
        M = ( 2.d0*DBLE(jm-1) + 1.d0 ) * pi / 2.d0
        Am(jm) = ((-1.d0)**(jm-1)) * 2.d0/(M**2.d0)
   
        alf2 = (1.d0-a5) / ( (M**2.d0)*(F0**2.d0)*a6 )
 
        bet4  =(1.d0-a4) / ( (M**2.d0)*(F0**2.d0)*a6 )
   
        gam4 = - alf2 * d1bar - h2bar - a2 * d3bar
        gam5 = - alf2 * d2bar - h3bar
        gam6 = - alf2 * d3bar
   
        del0 = - (M**2.d0) * a6;
   
        Delta = del2 * alf1 - del1 * alf2
        Delta0 = (del2 * alf0 - del0 * alf2) / Delta
        Delta1 = del2 * Delta0 - del1
        Delta2 = Delta1 * Delta0 + del0
 
        T0 = - Delta0
        T1 = - del2 * bet2 / Delta
        T2 = - (del2 * bet3 - alf2 * eps3) / Delta
        T3 = - (del2 * bet4 - alf2 * eps4) / Delta

        Tc1 = del2 * gam2 / Delta
        Tc2 = (del2 * gam3 - alf2 * et3) / Delta
        Tc3 = (del2 * gam4 - alf2 * et4) / Delta
        Tc4 = (del2 * gam5 - alf2 * et5) / Delta
!        Tc5 = (del2 * gam6 - alf2 * et6) / Delta         !Tc5=0
         
        csi1 = -Delta * Delta1 * T1
        csi2 = Delta * (-Delta1 * T2 + del2 * T1 + eps3)
        csi3 = Delta * (-Delta1 * T3 + del2 * T2 + eps4)
        csi4 = Delta * del2 * T3

        mu1 = Delta * Delta1 * Tc1
        mu2 = Delta * (Delta1 * Tc2 - del2 * Tc1 + et3)
        mu3 = Delta * (Delta1 * Tc3 - del2 * Tc2 + et4)
        mu4 = Delta * (Delta1 * Tc4 - del2 * Tc3 + et5) 
!        mu5 = Delta * (Delta1 * Tc5 - del2 * Tc4 + et6)  !mu5=0   
!        mu6 = -Delta * del2 * Tc5                        !mu6=0

        sigma0 = (Delta0 * csi1 + Delta * Delta2 * T1) / csi4
        sigma1 = (csi1 + Delta0 * csi2 + Delta * Delta2 * T2) / csi4
        sigma2 = (csi2 + Delta0 * csi3 + Delta * Delta2 * T3) / csi4
        sigma3 = (csi3 + Delta0 * csi4) / csi4

        rho0 = (Delta0 * mu1 - Delta * Delta2 * Tc1) / csi4
        rho1 = (mu1 + Delta0 * mu2 - Delta * Delta2 * Tc2) / csi4
        rho2 = (mu2 + Delta0 * mu3 - Delta * Delta2 * Tc3) / csi4
        rho3 = (mu3 + Delta0 * mu4 - Delta * Delta2 * Tc4) / csi4
        rho4 =  mu4 /csi4
!        rho4 = (mu4 + Delta0 * mu5 - Delta * Delta2 * Tc5) / csi4
!        rho5 = (mu5 + Delta0 * mu6 ) / csi4              !rho5=0
!        rho6 =  mu6 / csi4                               !rho6=0

!-----------------------------------------------------------------------
! quartic roots
!-----------------------------------------------------------------------
        coef(1) = dcmplx(sigma0,0.d0)
        coef(2) = dcmplx(sigma1,0.d0)
        coef(3) = dcmplx(sigma2,0.d0)
        coef(4) = dcmplx(sigma3,0.d0)
        coef(5) = (1.d0,0.d0)
               
        call zroots(coef,roots)

        lamb1(jm) = roots(1)
        lamb2(jm) = roots(2)
        lamb3(jm) = roots(3)
        lamb4(jm) = roots(4)

!-----------------------------------------------------------------------
! Wronskians
!-----------------------------------------------------------------------
        lam1 = roots(1)
        lam2 = roots(2)
        lam3 = roots(3)
        lam4 = roots(4)

        waux1 = lam3 * lam4 * (lam4-lam3)
        waux2 =-lam2 * lam4 * (lam4-lam2)
        waux3 = lam2 * lam3 * (lam3-lam2)
        W1det =-(waux1 + waux2 + waux3)

        wdet = -lam2 * lam3 * lam4 * W1det

        waux1 = lam3 * lam4 * (lam4-lam3)
        waux2 =-lam1 * lam4 * (lam4-lam1)
        waux3 = lam1 * lam3 * (lam3-lam1)
        W2det = (waux1 + waux2 + waux3)
                
        wdet = wdet - lam1 * lam3 * lam4 * W2det

        waux1 = lam2 * lam4 * (lam4-lam2)
        waux2 =-lam1 * lam4 * (lam4-lam1)
        waux3 = lam1 * lam2 * (lam2-lam1)
        W3det =-(waux1 + waux2 + waux3)

        wdet = wdet - lam1 * lam2 * lam4 * W3det

        waux1 = lam2 * lam3 * (lam3-lam2)
        waux2 =-lam1 * lam3 * (lam3-lam1)
        waux3 = lam1 * lam2 * (lam2-lam1)
        W4det = (waux1 + waux2 + waux3)

        wdet = wdet - lam1 * lam2 * lam3 * W4det

        wj(1) = W1det / wdet
        wj(2) = W2det / wdet
        wj(3) = W3det / wdet
        wj(4) = W4det / wdet

!-----------------------------------------------------------------------
! gj0
!-----------------------------------------------------------------------
        g10(jm) = wj(1) * ( rho0 + rho1*lam1 &
     &          + rho2*lam1**2.d0 + rho3*lam1**3.d0 + rho4*lam1**4.d0 )

        g20(jm) = wj(2) * ( rho0 + rho1*lam2 &
     &          + rho2*lam2**2.d0 + rho3*lam2**3.d0 + rho4*lam2**4.d0 )

        g30(jm) = wj(3) * ( rho0 + rho1*lam3 &
     &          + rho2*lam3**2.d0 + rho3*lam3**3.d0 + rho4*lam3**4.d0 )

        g40(jm) = wj(4) * ( rho0 + rho1*lam4 &
     &          + rho2*lam4**2.d0 + rho3*lam4**3.d0 + rho4*lam4**4.d0 )

!-----------------------------------------------------------------------
! gj1
!-----------------------------------------------------------------------
        g11(jm) = wj(1) * ( rho1 + rho2*lam1 &
     &          + rho3*lam1**2.d0 + rho4*lam1**3.d0 )

        g21(jm) = wj(2) * ( rho1 + rho2*lam2 &
     &          + rho3*lam2**2.d0 + rho4*lam2**3.d0 )

        g31(jm) = wj(3) * ( rho1 + rho2*lam3 &
     &          + rho3*lam3**2.d0 + rho4*lam3**3.d0 )

        g41(jm) = wj(4) * ( rho1 + rho2*lam4 &
     &          + rho3*lam4**2.d0 + rho4*lam4**3.d0 )

      end do
      
!-----------------------------------------------------------------------
! end of subroutine
!-----------------------------------------------------------------------
      return
      end
