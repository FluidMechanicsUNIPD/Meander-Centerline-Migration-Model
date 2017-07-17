!---------------------------------------------------------------------
! Runge-Kutta for second order ordinary differential equation
! input: zz,yj1,yj2       output yjj1,yjj2
!---------------------------------------------------------------------
      subroutine Rungeg(A,B,Cf,z0,kvk,GG0aux, &
      &            zz,deltaz,yj1,yj2,yjj1,yjj2,ig)
      implicit none
      integer ig
      real*8 A,B,Cf,kvk,GG0aux
      real*8 zz,deltaz,yj1,yj2,yjj1,yjj2
      real*8 zrk,yrk1,yrk2,fff,z0
      real*8 k11,k12,k21,k22,k31,k32,k41,k42

!---------------------------------------------------------------------
      zrk=zz
      yrk1=yj1
      yrk2=yj2
      call fgj(A,B,Cf,z0,kvk,GG0aux,zrk,deltaz,yrk1,yrk2,fff,ig)
      k11=deltaz*yrk2
      k12=deltaz*fff
  
      zrk=zz+deltaz/2.D0
      yrk1=yj1+k11/2.D0
      yrk2=yj2+k12/2.D0
      call fgj(A,B,Cf,z0,kvk,GG0aux,zrk,deltaz,yrk1,yrk2,fff,ig)
      k21=deltaz*yrk2
      k22=deltaz*fff
   
      zrk=zz+deltaz/2.D0
      yrk1=yj1+k21/2.D0
      yrk2=yj2+k22/2.D0
      call fgj(A,B,Cf,z0,kvk,GG0aux,zrk,deltaz,yrk1,yrk2,fff,ig)
      k31=deltaz*yrk2
      k32=deltaz*fff
   
      zrk=zz+deltaz
      yrk1=yj1+k31
      yrk2=yj2+k32
      call fgj(A,B,Cf,z0,kvk,GG0aux,zrk,deltaz,yrk1,yrk2,fff,ig)
      k41=deltaz*yrk2
      k42=deltaz*fff

      yjj1=yj1+(k11 + 2.D0*(k21+k31)+k41)/6.D0
      yjj2=yj2+(k12 + 2.D0*(k22+k32)+k42)/6.D0

!---------------------------------------------------------------------
      return
      end

!--------------------------------------------------------
! Subroutine for calculating the function fff
!--------------------------------------------------------
      subroutine fgj(A,B,Cf,z0,kvk,GG0aux,zrk,deltaz,yrk1,yrk2,fff,ig)
      implicit none
      integer ig
      real*8 zrk,deltaz,yrk1,yrk2,fff,pf,p0,p1
      real*8 A,B,Cf,kvk,z0,GG0aux,dGG0aux
      real*8 aux1,aux2,aux3,uzero,auxnu,nuTzero,Dauxnu,DnuTzero

!--------------------------------------------------------       
      aux1 = dlog(zrk/z0)
      aux2 = A * (zrk**2.D0 - z0**2.D0)
      aux3 = B * (zrk**3.D0 - z0**3.D0)
      uzero = (aux1 + aux2 +aux3) * dsqrt(Cf) / kvk

      if (zrk.lt.(1.D0 - deltaz/4.D0)) then
        auxnu = 1.D0 + 2.D0*A*(zrk**2.D0) + 3.D0*B*(zrk**3D0)
        nuTzero = kvk*zrk*(1.D0-zrk)/auxnu
        Dauxnu = 4.D0*A*zrk+9.D0*B*zrk**2.D0
        DnuTzero = kvk*((1.D0-2.D0*zrk)*auxnu- &
     &            zrk*(1.D0-zrk)*Dauxnu)/(auxnu**2.D0)
       else
        nuTzero =   0.0614D0
        DnuTzero = -0.0338D0
      end if

!--------------------------------------------------------       
      select case(ig)

        case(0)
        pf = 0.d0
        p0 = 0.d0
        p1 = DnuTzero/nuTzero

        case(1)
        pf = 1.d0/nuTzero
        p0 = 0.d0
        p1 = DnuTzero/nuTzero

        case(2)
        pf = -(uzero**2.D0)/nuTzero
        p0 = 0.d0
        p1 = DnuTzero/nuTzero

        case(12)
        pf = -1.d0/nuTzero
        p0 = 0.d0
        p1 = DnuTzero/nuTzero

        case(21)
        pf = (uzero*GG0aux)/nuTzero
        p0 = 0.d0
        p1 = DnuTzero/nuTzero

!        case(22)
!        pf = (1.d0-zrk)*uzero*dGG0aux/nuTzero
!        p0 = 0.d0
!        p1 = DnuTzero/nuTzero

      end select
      fff = -p1*yrk2 - p0*yrk1 + pf

!--------------------------------------------------------       
      return
      end
