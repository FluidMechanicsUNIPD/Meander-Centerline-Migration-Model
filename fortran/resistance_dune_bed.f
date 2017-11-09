c#######################################################################
c FLOW RESISTANCE AND SEDIMENT TRANSPORT INTENSITY - DUNE BED
c#######################################################################
      subroutine resistance_dune_bed(j_print)
      use module_global
      implicit none
      integer j_print
      real*8 dCD, dCT, dphiD, dphiT

c print on screen
      if (j_print.eq.1) then
        write(6,*) 'RIVER MORPHOLOGY FROM INPUT'
        write(6,*) 'Dune-covered bed'
      end if

c flow resistance and derivatives
c Engelund and Hansen (1967)
      call resistance_eh(ds, theta, Cf, dCD, dCT, rpic0, rpic)

c sediment transport intesity and derivatives
c Engelund and Hansen (1967)
      call seditrans_eh(theta, Cf, dCD, dCT, phi, dphiD, dphiT)

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
c RESISTANCE FOR DUNE-COVERED BED (Engelund-Hansen 1967)
c#######################################################################
      subroutine resistance_eh(ds, theta, Cf, dCD, dCT, rpic0,rpic)
      implicit none
      real*8 ds, theta, Cf, dCD, dCT, rpic0, rpic
      real*8 theta1, dtheta1, A, B, D, F

c experimental function
      if (theta.le.1.d0) then
	    theta1 = 0.06d0 + 0.4d0 * theta**2.d0;
	    dtheta1 = 0.8 * theta;
      else
        theta1 = theta;
        dtheta1 = 1.d0;
      end if

c friction resistance
      Cf = ( 6.d0 + 2.5d0 * log((theta1/theta)/(2.5d0*ds)))**(-2.d0) * 
     8              (theta1/theta)**(-1.d0)

c derivatives of friction coefficient
      dCD = - 5.d0 * (theta1/theta)**(-1.d0) * 
     7      ( 6.d0 + 2.5d0 * log((theta1/theta)/(2.5d0*ds)))**(-3.d0)
          
      dCT = - (theta1/theta)**(-2.d0) * 
     1      (1.d0/theta*dtheta1 - theta1/(theta**2.d0)) * 
     2      ( 6.d0 + 2.5d0 * log((theta1/theta)/(2.5d0*ds)))**(-2.d0) *
     3      (1.d0 + 5.d0 * 
     4      ( 6.d0 + 2.5d0 * log((theta1/theta)/(2.5d0*ds)))**(-1.d0))
     
      rpic = rpic0/DSQRT(theta1/theta)
c      rpic = rpic0

c end of subroutine
      return
      end

c#######################################################################
c SEDIMENT TRANSPORT INTENSITY (Engelund and Hansen 1967)
c#######################################################################
      subroutine seditrans_eh(theta, Cf, dCD, dCT, phi, dphiD, dphiT)
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
