!-----------------------------------------------------------------------
! LOAD THE INITIAL CONFIGURATION OF THE RIVER
!----------------------------------------------------------------------- 
      subroutine load_river
      use module_global
      use module_geometry
      implicit none
      integer j, ierr

! allocate coordinate arrays
      allocate ( x0(N0), y0(N0), stat=ierr )
      if(ierr.gt.0) go to 999
      
      select case (flagxy0)
      
! generate a straight path slightly perturbed in the transverse direction
        case(1)
          call randomconfig2 (N, Nrand, deltas0, stdv, x0, y0)

! import path from file
        case(2)
        
! open file
          filexy = ADJUSTL(filexy)
          filename = TRIM(dirIN)//TRIM(filexy)
          open(unit=2, file=filename, status='old')

! read coordinates
          do j = 1, N0
            read(2,*) x0(j), y0(j)
          end do
          
! close file
          close(2)

! warning
        case default
          stop 'ERROR! Flag for geometry file'
      end select

! end of the subroutine
      return
999   write(6,*) 'ERROR! Allocation in loading river data'
      stop
      end
