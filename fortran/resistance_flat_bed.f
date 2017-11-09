c#######################################################################
c FLOW RESISTANCE AND SEDIMENT TRANSPORT INTENSITY - FLAT BED
c#######################################################################
      subroutine resistance_flat_bed(j_print)
      use module_global
      implicit none
      integer j_bedload, j_print
      real*8 dCD, dCT, dphiD, dphiT

c bed load formula (1=MPM, 2=Parker)
      j_bedload = 2

c print on screen
      if (j_print.eq.1) then
        write(6,*) 'RIVER MORPHOLOGY FROM INPUT'
        write(6,*) 'Flat bed'
      end if

c flow resistance and derivatives
      call resistance_keulegan(ds, Cf, dCD, dCT, rpic0, rpic)

c sediment transport intensity and derivatives

c Meyer-Peter & Muller 1948
      if (j_bedload.eq.1) then
        call seditrans_mpm(theta, phi, dphiD, dphiT)
c Parker (1982, 1990)
      else if (j_bedload.eq.2) then
        call seditrans_parker(theta, phi, dphiD, dphiT)
c flag error
      else
        stop 'ERROR! Wrong flag for bedload transport'
      end if

c parameters for ZS model
      if (jmodel.eq.1) then
        call zs_terms(theta,phi,Cf,dCD,dCT,dphiD,dphiT,CD,CT,phiD,phiT)
      end if
      
c Froude number
      F0 = DSQRT(1.65d0 * ds * theta / Cf)

c print on screen
      if (j_print.eq.1) then
        write(6,'(a3,3(1x,f12.5))') 'Cf', Cf, CD, CT
        write(6,'(a3,3(1x,f12.5))') 'phi', phi, phiD, phiT
        write(6,'(a3,3(1x,f12.5))') 'Fr0', F0
        call dashline(6)
      end if
      
c end of subroutine
      return
      stop
      end
      
c#######################################################################
c RESISTANCE FOR PLANE BED (Meyer-Peter & Muller, Keulegan 1938)
c#######################################################################
      subroutine resistance_keulegan(ds, Cf, dCD, dCT, rpic0, rpic)
      implicit none
      real*8 ds, Cf, dCD, dCT, rpic0, rpic

c friction coeffient
      Cf = ( 6.d0 + 2.5d0 * log(1.d0/(2.5d0*ds)))**(-2.d0)

c derivatives of friction coefficients
      dCD = - 5.d0 * ( 6.d0 + 2.5d0 * log(1.d0/(2.5d0*ds)))**(-3.d0)
      dCT = 0.d0

      rpic = rpic0

c end of subroutine
      return
      end

c#######################################################################
c SEDIMENT TRANSPORT INTENSITY (Meyer-Peter & Muller 1948)
c#######################################################################
      subroutine seditrans_mpm(theta, phi, dphiD, dphiT)
      implicit none
      real*8 theta, thetacr, phi, dphiD, dphiT

c original formula

c critical threshold
c      thetacr = 0.047d0

c sediment transport intensity
c      phi = MAX(0.d0, 8.d0*(theta-thetacr)**(1.5d0))

c derivatives of sediment transport intensity
c      dphiD = 0.d0
c      dphiT = MAX(0.d0, 8.d0 * 1.5d0 *(theta-thetacr)**(1.5d0-1.d0))

c modified formula

c critical threshold
      thetacr = 0.0495d0

c sediment transport intensity
      phi = MAX(0.d0, 3.97d0*(theta-thetacr)**(1.5d0))

c derivatives of sediment transport intensity
      dphiD = 0.d0
      dphiT = MAX(0.d0, 3.97d0 * 1.5d0 *(theta-thetacr)**(1.5d0-1.d0))

c end of subroutine
      return
      end

c#######################################################################
c SEDIMENT TRANSPORT INTENSITY (Parker 1982, 1990)
c#######################################################################
      subroutine seditrans_parker(theta, phi, dphiD, dphiT)
      implicit none
      real*8 theta, thetar, phi, dphiD, dphiT
      real*8 csi, A, B, D, F, G0, dG

      thetar = 0.0386d0
      csi = theta / thetar
      A = 0.00218d0

      if (csi.lt.1.d0) then
        B = 14.2d0
        G0 = A * csi**B
        dG = B * G0 / theta
      else if ((csi.ge.1.d0).and.(csi.le.1.65d0)) then
        B = 14.2d0
        F = 9.28d0
        G0 = A * DEXP(B * (csi-1.d0) - F*(csi-1.d0)**2.d0)
        dG = G0/thetar * (B - 2.d0 * F *(csi-1.d0) )
      elseif (csi.gt.1.65d0) then
        B = 5474.d0
        F = 0.853d0
        D = 4.5d0
        G0 = A * B * (1.d0 - F/csi)**D
        dG = A*B*D*F*thetar/theta**2.d0 * (1.d0-F/csi)**(D-1.d0)
      end if

c sediment transport intensity
      phi = G0 * theta**1.5d0

c derivatives of sediment transport intensity
      dphiD = 0.d0
      dphiT = theta**1.5d0 * dG + 1.50 * phi/theta

c end of the subroutine
      return
      end
