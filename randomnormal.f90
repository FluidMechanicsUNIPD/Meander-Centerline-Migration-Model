!-----------------------------------------------------------------------
! PSEUDO-RANDOM NUMBER GENERATOR FROM NORMAL DISTRIBUTION N(A,B)
!-----------------------------------------------------------------------
! A : mean of the probability density function
! B : standard deviation of the probability density function
!-----------------------------------------------------------------------
! If the generator is used into a loop, before the loop you have to 
! call randinitialize() and then call *loop subroutine into the loop.
! If the generator is used for a one-time call, then you have to call
! the subroutine without *loop, which already contain randinitialize()
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! SINGLE NUMBER FOR ONE-TIME CALL
! If the random generator is one-time called, call randnormnumber(x)
! which contains randinitialize() 
!-----------------------------------------------------------------------
      subroutine randnormnumber(A, B, x)
      implicit none
      real*8 pi, y1, y2, z, A, B, x

! calculate value of pi = 3.141592653589793
      pi = 4.d0 * DATAN(1.d0)

! get two independent random numbers distributed in the interval (0, 1]
      call randinitialize()
      call RANDOM_NUMBER(y1)
      call RANDOM_NUMBER(y2)

! find a random value with a standard normal pdf N(0,1)
! [Box–Muller transform, 1958]
      z = SQRT(-2.0D0 * LOG(y1)) * COS(2.0D0*pi*y2)

! find a random value with a scaled normal pdf N(A,B)
      x = A + B * z
    
      return
      end subroutine
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
! SINGLE NUMBER INTO A LOOP
! If the random generator is into a loop, before it call randinitialize() 
! and then into the loop call randnormnumberloop(x)
!-----------------------------------------------------------------------
      subroutine randnormnumberloop(A, B, x)
      implicit none
      real*8 pi, y1, y2, z, A, B, x

! calculate value of pi = 3.141592653589793
      pi = 4.d0 * DATAN(1.d0)

! get two independent random numbers distributed in the interval (0, 1]
      call RANDOM_NUMBER(y1)
      call RANDOM_NUMBER(y2)

! find a random value with a standard normal pdf N(0,1)
! [Box–Muller transform, 1958]
      z = SQRT(-2.0D0 * LOG(y1)) * COS(2.0D0*pi*y2)

! find a random value with a scaled normal pdf N(A,B)
      x = A + B * z
    
      return
      end subroutine
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
! INITIALIZATION OF PSEUDO-RANDOM NUMBER GENERATOR
!-----------------------------------------------------------------------
      subroutine randinitialize()
      implicit none
      integer j, time(8), timems, k
      integer, allocatable :: seed(:)
! VALUE(1):	year 
! VALUE(2):	month 
! VALUE(3):	day of the month 
! VALUE(4):	difference with UTC in minutes 
! VALUE(5):	hour of the day 
! VALUE(6):	minutes of the hour 
! VALUE(7):	seconds of the minute 
! VALUE(8):	milliseconds of the second 

      call RANDOM_SEED(SIZE=k)
      allocate (seed(k))

! get the integer vector of system time
    1 call DATE_AND_TIME(VALUES=time)
      time(1)=time(1)*365*30*24*3600*1000
      time(2)=time(2)*30*24*3600*1000
      time(3)=time(3)*24*3600*1000
!      time(4)=time(4)*60*1000
      time(5)=time(5)*3600*1000
      time(6)=time(6)*60*1000
      time(7)=time(7)*1000
      
! find the value of system time in milliseconds
      timems = 0
      do j = 1,8
        timems = timems + time(j)
      end do

! assign the values for seed
      seed(:) = timems

! check the vaue of seed
      if(seed(1).eq.0) then
        write(6,*) 'Warning: pseudo-random number generator!'
        go to 1
      end if

! restart the pseudorandom number generator
      call RANDOM_SEED(PUT=seed)
      return
      end subroutine
!-----------------------------------------------------------------------
