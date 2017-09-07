c#######################################################################
c FLOW RESISTANCE AND SEDIMENT TRANSPORT INTENSITY
c#######################################################################
      subroutine resistance(j_thetac, j_bedload, j_print, trasp, fondo)
      use module_global
      implicit none
      integer fondo, j_thetac, j_bedload, j_print, trasp
      real*8 thetac, dCD, dCT, dphiD, dphiT

c j_thetac  : thetac               (1=Brownlie, 2=Van Rijn)
c j_bedload : bed loaf formula     (1=MPM, 2=Parker)
c j_print   : print on screen      (0=no, 1=yes)

c-----------------------------------------------------------------------
c BED TYPE CHARACTERIZATION
c-----------------------------------------------------------------------
      call bed_characterization(Rp, j_thetac, theta, thetac, fondo)

c print on screen
      if (j_print.eq.1) then
      
        write(6,*) 'RIVER MORPHOLOGY FROM INPUT'
        if (flagbed.eq.0) then
          write(6,*) 'given by Rp (see below)'
        else if (flagbed.eq.1) then
          write(6,*) 'Flat bed'
        else if (flagbed.eq.2) then
          write(6,*) 'Dune-covered bed'
        end if
      
        write(6,*) 'RIVER MORPHOLOGY FROM DATA'
        if (fondo.eq.1) then
          write(6,*) 'Flat bed'
        else if (fondo.eq.2) then
          write(6,*) 'Dune-covered bed'
        end if

      end if

c-----------------------------------------------------------------------
c BED AND TRANSPORT TYPES FROM Rp + VAN RIJN METHOD(1984)
      if (flagbed.eq.0) then
c-----------------------------------------------------------------------

c FRICTION COEFFICIENT

c plane bed
        if (fondo.eq.1) then
          call resistance_planebed(ds, Cf, dCD, dCT, rpic0, rpic)
c dune-covered bed
        else if (fondo.eq.2) then
          call resistance_dunebed(ds, theta, Cf, dCD, dCT, rpic0, rpic)
        end if

c SEDIMENT TRANSPORT CHARACTERIZATION

c no sediment transport transport
        if (theta.lt.thetac) then
          trasp = 0
          write(6,*) 'ERROR! No sediment transport, theta = ', theta, 
     6       ' < thetac = ', thetac
          stop 
c sediment transport
        else
          call seditrans_characterization(theta, ds, Rp, Cf, F0, trasp)
        end if

c SEDIMENT TRANSPORT INTENSITY
        
c bedload
        if (trasp.eq.1) then

c stampa a video
          if (j_print.eq.1) then
            write(6,*) 'Bedload transport'
          end if

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

c total load = bedload transport + suspended transport
        else if (trasp.eq.2) then

c stampa a video
          if (j_print.eq.1) then
            write(6,*) 'Suspended transport'
          end if

c Engelund and Hansen (1967)
          call seditrans_EH(theta, Cf, dCD, dCT, phi, dphiD, dphiT)

c end if for sediment transport intensity
        end if

c-----------------------------------------------------------------------
c BED AND TRANSPORT TYPES FROM FROM INPUT SETTING - FLAT BED
      else if (flagbed.eq.1) then
c-----------------------------------------------------------------------

c flow resistance and derivatives
        call resistance_planebed(ds, Cf, dCD, dCT, rpic0, rpic)

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
        
c bed characterization
        call seditrans_characterization(theta, ds, Rp, Cf, F0, trasp)
        
c-----------------------------------------------------------------------
c BED AND TRANSPORT TYPES FROM FROM INPUT SETTING - DUNE COVERED BED
      else if (flagbed.eq.2) then
c-----------------------------------------------------------------------

c flow resistance and derivatives
        call resistance_dunebed(ds, theta, Cf, dCD, dCT, rpic0, rpic)

c sediment transport intesity and derivatives
c Engelund and Hansen (1967)
        call seditrans_EH(theta, Cf, dCD, dCT, phi, dphiD, dphiT)
     
c bed characterization
        call seditrans_characterization(theta, ds, Rp, Cf, F0, trasp)

c-----------------------------------------------------------------------
c end if for resistance and intensity
      else
        stop 'ERROR! Wrong flag for bed configuration'
      end if
c-----------------------------------------------------------------------

c parameters for ZS model
      if (jmodel.eq.1) then
        call zs_terms(theta,phi,Cf,dCD,dCT,dphiD,dphiT,CD,CT,phiD,phiT)
      end if

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
c BED TYPE CHARACTERIZATION (Van Rijn 1984)
c#######################################################################
      subroutine bed_characterization(Rp,j_thetac,theta,thetac,fondo)
      implicit none
      integer j_thetac, fondo
      real*8 Rp, dsdim, dsmm, dstar, thetac, theta, Tvr

c sediment size (meters)
      dsdim = ( Rp*Rp / (10.d0**12.d0*1.65d0*9.81d0) )**(1.d0/3.d0)

c sediment size (millimeters)
      dsmm = dsdim*1000.d0

c dimensionless diameter
      dstar = Rp**(2.d0/3.d0)

c critical Shield number (Brownlie 1981)
      if (j_thetac.eq.1) then
        thetac = 0.22d0*Rp**(-0.6d0)+0.06d0*DEXP(-17.77d0*Rp**(-0.6d0))

c critical Shield number (Van Rijn 1989)
      else if (j_thetac.eq.2) then

        if (dstar.le.4.d0) then
          thetac = 0.24d0 / dstar
        else if ((dstar.gt.4.d0).and.(dstar.le.10.d0)) then
          thetac = 0.14d0 * dstar**(-0.64d0)
        else if ((dstar.gt.10.d0).and.(dstar.le.20.d0)) then
          thetac = 0.04d0 * dstar**(-0.10d0)
        else if ((dstar.gt.20.d0).and.(dstar.le.150.d0)) then
          thetac = 0.013d0 * dstar**0.29d0
        else if (dstar.gt.150.d0) then
          thetac = 0.055d0
        end if

c flag error
      else
        stop 'ERROR! Wrong flag for critical Shield stress'
      end if

c transport parameter
      Tvr = (theta-thetac) / thetac

c type classification (Van Rijn 1984)
      if (dstar.le.10.d0) then

        if (Tvr.le.3.d0) then
          fondo = 1     ! plane bed
        else if ((Tvr.gt.3.d0).and.(Tvr.le.15.d0)) then
          fondo = 2     ! dune-covered bed
        else if (Tvr.gt.15.d0) then
          fondo = 1     ! plane bed
        end if

      else

        if (Tvr.gt.15.d0) then
          fondo = 1     ! plane bed
        else
          fondo = 2     ! dune-covered bed
        end if

      end if

c end of subroutine
      return
      end

c#######################################################################
c SEDIMENT TRANSPORT TYPE CHARACTERIZATION (Van Rijn 1984)
c#######################################################################
      subroutine seditrans_characterization(theta,ds,Rp,Cf,F0,trasp)
      implicit none
      integer trasp
      real*8 theta, ds, Cf, Rp, F0, F02
      real*8 U0, ustar, dsdim, CC, DD, ws0, wsdim, Frstar

c sediment size (meters)
      dsdim = ( Rp*Rp / (10.d0**12.d0*1.65d0*9.81d0) )**(1.d0/3.d0)

c Froude number
      F0 = DSQRT(1.65d0 * ds * theta / Cf)
      F02 = F0*F0

c uniform flow velocity
      U0 = DSQRT(F02 * 9.81d0 * dsdim / ds)

c shear velocity
      ustar = U0 * DSQRT(Cf)

c dimensionless settling velocity (Parker 1978, Van Oyen et al. 2011)
      CC = log10(Rp)
      DD = -1.181d0 + 0.966d0*CC - 0.1804d0*CC**2.d0 +
     8      0.003746d0*CC**3.d0 + 0.0008782d0*CC**4.d0
      ws0 = 10.d0**DD

c dimensional settling velocity
      wsdim = ws0 * DSQRT(1.65d0 * 9.81d0 * dsdim)

c threshold for suspended load
      Frstar = ustar / wsdim

c sediment flux characterization (Van Rjin 1984)
      if (Rp.gt.31.62d0) then

c bedload transport
        if (Frstar.lt.0.4d0) then
          trasp = 1
c bed load + suspended load
        else
          trasp = 2
        end if

      else

c bed load
        if (Frstar.lt.(4.d0*Rp**(-2.d0/3.d0))) then
          trasp = 1
c bed load + suspended load
        else
          trasp = 2
        end if

      end if

c end of subroutine
      return
      end

c#######################################################################
c RESISTANCE FOR PLANE BED (Meyer-Peter & Muller, Keulegan 1938)
c#######################################################################
      subroutine resistance_planebed(ds, Cf, dCD, dCT, rpic0, rpic)
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
c RESISTANCE FOR DUNE-COVERED BED (Engelund-Hansen 1967)
c#######################################################################
      subroutine resistance_dunebed(ds, theta, Cf, dCD, dCT, rpic0,rpic)
      implicit none
      real*8 ds, theta, Cf, dCD, dCT, rpic0, rpic
      real*8 theta1, dtheta1, A, B, D, F

c experimental function
      if (theta.lt.0.06d0) then
        theta1 = theta
        dtheta1 = 1.d0
      else if ((theta.ge.0.06d0).and.(theta.le.0.55d0)) then
        theta1 = 0.06d0 + 0.4d0 * theta**2.0d0
        dtheta1 = 0.8d0 * theta
c        theta1 = 0.06d0 + 0.3d0 * theta**1.5
c        dtheta1 = 0.45d0 * theta**0.5
      else if ((theta.gt.0.55d0).and.(theta.le.0.80d0)) then
        theta1 = 0.06d0 + 0.4d0 * theta**2.0d0
        dtheta1 = 0.8d0 * theta
c        theta1 = 0.06d0 + 0.3d0 * theta**1.5
c        dtheta1 = 0.45d0 * theta**0.5
      else if ((theta.gt.0.80d0).and.(theta.lt.1.1068d0)) then
        A = 0.439218723d0
        B = 1.072070082d0 - A
        D = 0.85d0
        F = 1.1068d0 - D
        theta1 = A + B/F * (theta-D)
        dtheta1 = B/F
      else if (theta.ge.1.1068d0) then
        theta1 = (0.3d0 + 0.7d0*theta**(-1.8d0))**(-0.56d0)
        dtheta1 = 0.56d0 * (0.3d0 + 0.7d0*theta**(-1.8d0))**(-1.56d0)
     8          * 0.7d0*1.8d0*theta**(-2.8d0)
      end if

c friction resistance
      Cf = ( 6.d0 + 2.5d0 * log((theta1/theta)/(2.5d0*ds)))**(-2.d0) * 
     8              (theta1/theta)**(-1.d0)

c derivatives of friction coefficient
      dCD = - 5.d0 * (theta1/theta)**(-1.d0) * 
     7      ( 6.d0 + 2.5d0 * log((theta1/theta)/(2.5d0*ds)))**(-3.d0)
          
      dCT = - (theta1/theta)**(-2.d0) * 
     1      (1.d0/theta*dtheta1 - theta1/theta**2.d0) * 
     2      ( 6.d0 + 2.5d0 * log((theta1/theta)/(2.5d0*ds)))**(-2.d0) *
     3      (1.d0 + 5.d0 * 
     4      ( 6.d0 + 2.5d0 * log((theta1/theta)/(2.5d0*ds)))**(-1.d0))
     
      rpic = rpic0/DSQRT(theta1/theta)
c      rpic = rpic0

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

c#######################################################################
c SEDIMENT TRANSPORT INTENSITY (Engelund and Hansen 1967)
c#######################################################################
      subroutine seditrans_EH(theta, Cf, dCD, dCT, phi, dphiD, dphiT)
      implicit none
      real*8 theta, Cf, dCD, dCT, phi, dphiD, dphiT, A, B, F

      A = 0.05d0
      B = 2.50d0

c sediment transport intensity
      phi  = A/Cf * theta**B

c derivatives of sediment transport intensity
      dphiD = -phi/Cf * dCD
      dphiT = -phi/Cf * dCT + B*phi/theta

c end of subroutine
      return
      end
      
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
