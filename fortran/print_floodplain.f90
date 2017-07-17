!-----------------------------------------------------------------------
! PRINT FILES OF OXBOW LAKES AND POINT BARS
!-----------------------------------------------------------------------
      subroutine print_floodplain
      use module_global
      use module_geometry
      implicit none
      integer i

!-----------------------------------------------------------------------
! print oxbow arrays
!-----------------------------------------------------------------------

      filename = TRIM(dirTMP)//TRIM(simname)//'_xo.dat'
      open(unit=51, file=filename, status='replace')
      filename = TRIM(dirTMP)//TRIM(simname)//'_yo.dat'
      open(unit=52, file=filename, status='replace')
      do i = 1, oxrow(Nnco+1)
        write(51,*) xo(i)
        write(52,*) yo(i)
      end do
      close(51)
      close(52)

      filename = TRIM(dirTMP)//TRIM(simname)//'_ox_sum.dat'
      open(unit=53, file=filename, status='replace')
      filename = TRIM(dirTMP)//TRIM(simname)//'_ox_row.dat'
      open(unit=54, file=filename, status='replace')
      do i = 1, Nnco+1
        write(53,*) oxsum(i)
        write(54,*) oxrow(i)
      end do
      close(53)
      close(54)

!-----------------------------------------------------------------------
! print point bar arrays
!-----------------------------------------------------------------------

      filename = TRIM(dirTMP)//TRIM(simname)//'_xb.dat'
      open(unit=51, file=filename, status='replace')
      filename = TRIM(dirTMP)//TRIM(simname)//'_yb.dat'
      open(unit=52, file=filename, status='replace')
      do i = 1, pbrow(Nnco+1)
        write(51,*) xb(i)
        write(52,*) yb(i)
      end do
      close(51)
      close(52)

      filename = TRIM(dirTMP)//TRIM(simname)//'_pb_sum.dat'
      open(unit=53, file=filename, status='replace')
      filename = TRIM(dirTMP)//TRIM(simname)//'_pb_row.dat'
      open(unit=54, file=filename, status='replace')
      do i = 1, Nnco+1
        write(53,*) pbsum(i)
        write(54,*) pbrow(i)
      end do
      close(53)
      close(54)

!-----------------------------------------------------------------------
! print cutoff number
!-----------------------------------------------------------------------

      filename = TRIM(dirTMP)//TRIM(simname)//'_nnco.dat'
      open(unit=55, file=filename, status='replace')
      write(55,*) Nnco
      close(55)

!-----------------------------------------------------------------------
! end of subroutine
!-----------------------------------------------------------------------
      return
      end
