!----------------------------------------------------------------------
! SAVITZKY-GOLAY FILTERING
! A Savitzky-Golay filter smooth the input signal x of size N into the
! output signal y, through the following relation:
! y(i) = 1/h * sum _{j = -(np-1)/2}^{(np-1)/2} aj * x_{i+j}
! which involves np terms (np must be odd).
!----------------------------------------------------------------------
      subroutine savgolfilter_new (np, n, x, y)
      implicit none
      integer i, j, n, np, nh
      real*8 x(n), y(n), h
      real*8 a(-INT(REAL(np-1)/2.0):INT(REAL(np-1)/2.0))

! assign the half-number of points
      nh = INT(REAL(np-1)/2.0)
      
! set the Savitzky-Golay coefficients and the normalization term h
      select case (np)
      case(5)
        h = 35.d0
        a(0:nh) = (/ 17.d0, 12.d0, -3.d0 /)
      case(7)
        h = 21.d0
        a(0:nh) = (/ 7.d0, 6.d0, 3.d0, -2.d0 /)
      case(9)
        h = 231.d0
        a(0:nh) = (/ 59.d0, 54.d0, 39.d0, 14.d0, -21.d0 /)
      case(11)
        h = 429.d0
        a(0:nh) = (/ 89.d0, 84.d0, 69.d0, 44.d0, 9.d0, -36.d0 /)
      case(13)
        h = 143.d0
        a(0:nh) = (/ 25.d0, 24.d0, 21.d0, 16.d0, 9.d0, 0.d0, -11.d0 /)
      case(15)
        h = 1105.d0
        a(0:nh) = (/ 167.d0, 162.d0, 147.d0, 122.d0, 87.d0,
     1                42.d0, -13.d0, -78.d0 /)
      case(17)
        h = 323.d0
        a(0:nh) = (/ 43.d0, 42.d0, 39.d0, 34.d0, 27.d0,
     2               18.d0, 7.d0, -6.d0, -21.d0 /)
      case(19)
        h = 2261.d0
        a(0:nh) = (/ 269.d0, 264.d0, 249.d0, 224.d0,
     3               189.d0, 144.d0, 89.d0, 24.d0, -21.d0, -136.d0 /)
      case(21)
        h = 3059.d0
        a(0:nh) = (/ 329.d0, 324.d0, 309.d0, 284.d0, 249.d0, 204.d0,
     4               149.d0, 84.d0, 9.d0, -76.d0, -171.d0 /)
      case(23)
        h = 805.d0
        a(0:nh) = (/ 79.d0, 78.d0, 75.d0, 70.d0, 63.d0, 54.d0,
     5               43.d0, 30.d0, 15.d0, -2.d0, -21.d0, -42.d0 /)
      case(25)
        h = 5175.d0
        a(0:nh) = (/ 467.d0, 462.d0, 447.d0, 422.d0, 387.d0, 343.d0,
     6        287.d0, 222.d0, 147.d0, 62.d0, -33.d0, -138.d0, -253.d0 /)
      case default
        stop 'ERROR! Wrong value for Savitzky-Golay filtering'
      end select
      
! set the mirror part of the coefficient array
      do j = 1, nh
        a(-j) = a(j)
      end do

! set the initial part of the signal equal to the input
      do i = 1, nh
        y(i) = x(i)
      end do

! filter the central part of the input signal
      do i = nh+1, n-nh
        y(i) = 0.d0
        do j = -nh, +nh
          y(i) = y(i) + a(j) * x(i+j) / h
        end do
      end do
      
! set the final part of the signal equal to the input
      do i = n-nh+1, n
        y(i) = x(i)
      end do
      
! end of the subroutine 
      return
      end
