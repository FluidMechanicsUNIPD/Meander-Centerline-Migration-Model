c-----------------------------------------------------------------------
c Coefficients k0, k1, k3, k4
c-----------------------------------------------------------------------
      subroutine k0134(Cf,Nz,k0,k1,k3,k4)

      implicit none
      integer Nz,ig, j

      real*8 k0, k1, k3, k4, kvk, A, B, Cf
      real*8 I10, I20, I30, z0, deltaz, const, xiini, xiend, deltaxi, xi
      real*8 zz, yj1, yj2, yjj1, yjj2, q1, q2, ha0, q0
      real*8 aux1, aux2, aux3, Sumg0, Sumg1, Sumg2, Du0
      real*8 g0(Nz), dg0(Nz), g1(Nz), dg1(Nz), g2(Nz), dg2(Nz)
      real*8 zg(Nz), u0(Nz), ha1
      real*8 GG0(Nz), dGG0(Nz), GG1(Nz), dGG1(Nz), GG2(Nz), dGG2(Nz)
      real*8 GG0aux, k0z(Nz), k1z(Nz)

c-----------------------------------------------------------------------
c dati in ingresso:  Nz, A, B, kvk, Cf 
      kvk = 0.41d0 !von karman constant
      A = 1.84d0   !coeff. distribuzione niT
      B = -1.56d0  !coeff. distribuzione niT

      z0 = dexp(-kvk/sqrt(Cf) - 0.777d0)
      deltaz = (1.d0-z0)/dble(Nz-1)

c-----------------------------------------------------------------------
c turbulent viscosity nuT0(z) and velocity u0(z)
      I10 = dlog(z0)-dlog(1.d0-z0)
      I20 = -(z0+dlog(1.d0-z0))
      I30 = -(0.5d0*z0**2.d0 + z0 + dlog(1.d0-z0))
      const = -(I10 + 2.d0*A*I20 + 3.d0*B*I30)

      xiini = 0.d0
      xiend = -dlog(z0)
      deltaxi = (xiend - xiini)/dble(Nz-1)

c-----------------------------------------------------------------------
c Calcolo G0
c-----------------------------------------------------------------------
      g0(1)  = 0.d0
      dg0(1) = 1.d0
      g1(1)  = 0.d0
      dg1(1) = 1.d0
      g2(1)  = 0.d0
      dg2(1) = 1.d0

      do j = 1, Nz-1
        xi = xiini + dble(j-1)*deltaxi
        zg(j) = z0*dexp(xi)
        deltaz = zg(j)*(dexp(deltaxi)-1.d0)

c calcolo di u0(z)
        aux1 = dlog(zg(j)/z0)
        aux2 = A*(zg(j)**2.d0-z0**2.d0)
        aux3 = B*(zg(j)**3.d0 -z0**3.d0)
        u0(j)= (aux1 + aux2 + aux3)*sqrt(Cf)/kvk   
   
c calcolo g0
        ig=0
        zz = zg(j)
        yj1 = g0(j)
        yj2 = dg0(j)
        call Rungeg(A,B,Cf,z0,kvk,GG0aux,zz,deltaz,yj1,yj2,yjj1,yjj2,ig)
        g0(j+1) = yjj1  
        dg0(j+1) = yjj2

c calcolo g1
        ig = 1
        zz = zg(j)
        yj1 = g1(j)
        yj2 = dg1(j)
        call Rungeg(A,B,Cf,z0,kvk,GG0aux,zz,deltaz,yj1,yj2,yjj1,yjj2,ig)
        g1(j+1) = yjj1  
        dg1(j+1) = yjj2

c calcolo g2
        ig = 2
        zz = zg(j)
        yj1 = g2(j)
        yj2 = dg2(j)
        call Rungeg(A,B,Cf,z0,kvk,GG0aux,zz,deltaz,yj1,yj2,yjj1,yjj2,ig)
        g2(j+1) = yjj1 
        dg2(j+1) = yjj2
      end do

      xi = xiini+(Nz-1)*deltaxi
      zg(Nz) = z0*dexp(xi)
      q1 = -dg1(Nz)/dg0(Nz)
      q2 = -dg2(Nz)/dg0(Nz)

c calcolo di u0(N)
      aux1 = dlog(zg(Nz)/z0)
      aux2 = A*(zg(Nz)**2.d0 -z0**2.d0)
      aux3 = B*(zg(Nz)**3.d0 -z0**3.d0)
      u0(Nz) = (aux1 + aux2 + aux3)*sqrt(Cf)/kvk
        
c calcolo integrali funzioni g1,g2,g3
      call CavSimp3(zg,g0,g1,g2,Sumg0,Sumg1,Sumg2,Nz)

      ha0 = -(Sumg2 + q2*Sumg0)/(Sumg1 + q1*Sumg0)
      q0 = q1*ha0+q2

      do j = 1, Nz
        GG0(j) = q0*g0(j) + ha0*g1(j) + g2(j)
        dGG0(j) = q0*dg0(j) + ha0*dg1(j) + dg2(j)
      end do

c-----------------------------------------------------------------------
c Calcolo G1
c-----------------------------------------------------------------------
      g2(1) = 0.D0
      dg2(1) = 1.D0
      do j = 1, Nz-1
        deltaz = zg(j)*(dexp(deltaxi)-1.d0)

c calcolo g21
        ig = 21
        zz = zg(j)
        yj1 = g2(j)
        yj2 = dg2(j)
        GG0aux = (GG0(j)+GG0(j+1))/2.d0
        call Rungeg(A,B,Cf,z0,kvk,GG0aux,
     5                   zz,deltaz,yj1,yj2,yjj1,yjj2,ig)
        g2(j+1) = yjj1  
        dg2(j+1) = yjj2
      end do

      q2 = -dg2(Nz)/dg0(Nz)

c calcolo integrali funzioni g1,g2,g3
      call CavSimp3(zg,g0,g1,g2,Sumg0,Sumg1,Sumg2,Nz)

      ha1 = -(Sumg2 + q2*Sumg0)/(Sumg1 + q1*Sumg0)
      q0 = q1*ha1+q2

      do j = 1, Nz
        GG1(j) = q0*g0(j) + ha1*g1(j) + g2(j)
        dGG1(j) = q0*dg0(j) + ha1*dg1(j) + dg2(j)
      end do

c-----------------------------------------------------------------------
c calcolo k3,k4,k5
c-----------------------------------------------------------------------
      Du0 = sqrt(Cf)*(1.d0/z0 + 2.d0*A*z0 + 3.d0*B*z0**2.d0)/kvk
      k3 = dGG0(1)/Du0
      k4 = dGG1(1)/Du0

c-----------------------------------------------------------------------
c calcolo k0,k1,k2
c-----------------------------------------------------------------------
      do j = 1, Nz
        k0z(j) = GG0(j)*u0(j)
        k1z(j) = GG1(j)*u0(j)
      end do

      call CavSimp2(zg,k0z,k1z,k0,k1,Nz)

c-----------------------------------------------------------------------
c end of subroutines
c-----------------------------------------------------------------------
      return
      end
