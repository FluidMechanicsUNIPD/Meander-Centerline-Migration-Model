c-----------------------------------------------------------------------
c EXCESS NEAR-BANK LONGITUDINAL VELOCITY
c Zolezzi & Seminara model (2001), with free boundary conditions,
c i.e., cm1=cm2=cm3=cm4=0.
c-----------------------------------------------------------------------
      subroutine dUZS
      use module_global
      use module_zs
      use module_geometry
      implicit none
      
      integer j, jm, jev, jsend, jtoll
      real*8 TOLL, c_bc
      parameter(TOLL=1.D-04)
      real*8, allocatable :: um(:)
      complex*16 lambda, gj0, CONV
      complex*16 SEMIANA1, SEMIANA2
      complex*16, allocatable :: UPSTR(:), DWSTR(:), LOCAL(:)
      complex*16, allocatable :: DWSTR_BC(:), UPSTR_BC(:)
c UPSTR : upstream propagating influence (i.e. negative-exp convolution)
c DWSTR : downstream propagating influence (i.e. positive-exp convolution)
c DWSTR_BC : upstream-propagating donwstream boundary conditions (positive exp) 
c UPSTR_BC : downstream-propagating upstream boundary conditions (negative exp)
c LOCAL : local effect of the curvature     
       
c initialize the excess near-bank velocity array
      dU(:) = 0.d0
      
c initialize the five components
      allocate ( UPSTR(N),DWSTR(N),UPSTR_BC(N),DWSTR_BC(N),LOCAL(N) )
      allocate ( um(N) )

c-----------------------------------------------------------------------
c loop over Fourier components
      do jm = 1, Mdat
c-----------------------------------------------------------------------

c reset the five components
        UPSTR(:) = 0.d0
        DWSTR(:) = 0.d0
        UPSTR_BC(:) = 0.d0
        DWSTR_BC(:) = 0.d0
        LOCAL(:) = 0.d0

c-----------------------------------------------------------------------
c loop over the eigenvalues
        do jev = 1, 4
c-----------------------------------------------------------------------
        
c set current eigenvalues
          select case (jev)
            case(1)
              lambda = lamb1(jm)
              gj0 = g10(jm)
              c_bc = 0.d0
            case(2)
              lambda = lamb2(jm)
              gj0 = g20(jm)
              c_bc = 0.d0
            case(3)
              lambda = lamb3(jm)
              gj0 = g30(jm)
              c_bc = 0.d0
            case(4)
              lambda = lamb4(jm)
              gj0 = g40(jm)
              c_bc = 0.d0
          end select
        
c initialize pointer
          jtoll = 1
          
c UPSTREAM PROPAGATING INFLUENCE (+ DOWNSTREAM BOUNDARY CONDITIONS)
          if ( DBLE(lambda).gt.0.d0) then
        
            do j = 1, N
              if ( DEXP(-DBLE(lambda)*(s(j)-s(1))).gt.TOLL ) then
                jtoll = jtoll + 1
              end if
            end do              

            do j = 1, N, +1
              if (j.le.(N-1)) then
                jsend = MIN(j+jtoll-1,N)
                CONV = SEMIANA1(j,N,c,jsend,deltas,lambda)
                UPSTR(j) = UPSTR(j) - Am(jm) * gj0 * CONV
              end if
              DWSTR_BC(j) = DWSTR_BC(j) + c_bc * 
     6                      DEXP(-DBLE(lambda)*(s(N)-s(j)))  
            end do
            
c DOWNSTREAM PROPAGATING INFLUENCE (+ UPSTREAM BOUNDARY CONDITIONS)
          else if ( DBLE(lambda).lt.0.d0) then
            
            do j = 1, N
              if ( DEXP(DBLE(lambda)*(s(j)-s(1))).gt.TOLL ) then
                jtoll = jtoll + 1
              end if
            end do 
            
            do j = N, 1, -1
              if (j.ge.2) then
                jsend = MAX(j-jtoll+1,1)
                CONV = SEMIANA2(j,N,c,jsend,deltas,lambda)
                DWSTR(j) = DWSTR(j) + Am(jm) * gj0 * CONV
              end if
              UPSTR_BC(j) = UPSTR_BC(j) + c_bc * 
     6                      DEXP(DBLE(lambda)*(s(j)-s(1)))         
            end do
            
          end if
          
c-----------------------------------------------------------------------
c end loop over the eigenvalues
        end do
c-----------------------------------------------------------------------
            
c local effect of the curvature, longitudinal velocity and excess near-
c bank velocity
        do j = 1, N
          LOCAL(j) = Am(jm) * c(j) * (g11(jm)+g21(jm)+g31(jm)+g41(jm))
          um(j) = DBLE( UPSTR(j) + UPSTR_BC(j) +
     6                  DWSTR(j) + DWSTR_BC(j) + LOCAL(j) )
          dU(j) = dU(j) + 2.d0 * um(j) * (-1.d0)**(jm-1)
        end do

c-----------------------------------------------------------------------
c end loop over the Fourier components
      end do
c-----------------------------------------------------------------------

c end of the subroutine
      return
      end
