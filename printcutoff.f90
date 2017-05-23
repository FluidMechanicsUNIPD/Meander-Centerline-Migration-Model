!-----------------------------------------------------------------------
! PRINT OUTPUT FILES
!-----------------------------------------------------------------------
      subroutine printcutoff(Ncb, xcb, ycb, Nco, xco, yco)
      use module_global
      use module_geometry
      implicit none
      character*20 p1
      integer j, Ncb, Nco
      real*8 xco(Nco), yco(Nco)     ! oxbow coordinates
      real*8 xcb(Ncb), ycb(Ncb)     ! point bar coordinates

!-----------------------------------------------------------------------
! print current oxbow lake
!-----------------------------------------------------------------------

! file name
      write( p1, '(i10)' )  Nnco
      p1 = adjustl(p1)
      filename = TRIM(dirOUT)//'oxbow_'//TRIM(p1)//'.out'
      
! open file
      open(unit=22, file=filename, status='replace')
      
! print current oxbow
      do j = 1, Nco
        write(22,*) xco(j), yco(j)
      end do
      
! close file
      close(22)

!-----------------------------------------------------------------------
! print current point bar
!-----------------------------------------------------------------------

! file name
      write( p1, '(i10)' )  Nnco
      p1 = adjustl(p1)
      filename = TRIM(dirOUT)//'point_bar_'//TRIM(p1)//'.out'

! open file
      open(unit=22, file=filename, status='replace')

! print current oxbow
      do j = 1, Ncb
        write(22,*) xcb(j), ycb(j)
      end do

! close file
      close(22)

!-----------------------------------------------------------------------
! end of subroutine
!-----------------------------------------------------------------------
      return
      end