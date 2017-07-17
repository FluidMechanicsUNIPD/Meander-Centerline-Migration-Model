!-----------------------------------------------------------------------
! LOAD FLOODPLAIN STRUCTURES: oxbow lakes and scroll bars
!----------------------------------------------------------------------- 
      subroutine load_floodplain
      use module_global
      use module_geometry
      implicit none
      integer i, ierr

      select case(flag_ox)
      
! no existing structure of the floodplain
      case(0)
        write(6,*) 'existing structure of the floodplain : no'
        Nnco0 = 0
        
! existing structure of the floodplain
      case(1)
        write(6,*) 'existing structure of the floodplain : yes'

! cutoff counter
        filename = TRIM(dirIN)//TRIM(simname)//'_nnco.dat'
        open(unit=1, file=filename, status='old')
        read(1,*) Nnco0
        close(1)

!-----------------------------------------------------------------------
! existing oxbow lakes
!-----------------------------------------------------------------------
        
! oxbow pointers
        allocate ( oxrow(Nnco0+1), oxsum(Nnco0+1), stat=ierr )
        if(ierr.gt.0) go to 999    

        filename = TRIM(dirIN)//TRIM(simname)//'_ox_row.dat'
        open(unit=12, file=filename, status='old')
        filename = TRIM(dirIN)//TRIM(simname)//'_ox_sum.dat'
        open(unit=13, file=filename, status='old')

        do i = 1, Nnco0+1
          read(12,*) oxrow(i)
          read(13,*) oxsum(i)
        end do

        close(12)
        close(13) 

! oxbow coordinates
        allocate ( xo(oxrow(Nnco0+1)), yo(oxrow(Nnco0+1)), stat=ierr )
        if(ierr.gt.0) go to 999     
 
        filename = TRIM(dirIN)//TRIM(simname)//'_xo.dat'
        open(unit=14, file=filename, status='old')
        filename = TRIM(dirIN)//TRIM(simname)//'_yo.dat'
        open(unit=15, file=filename, status='old')

        do i = 1, oxrow(Nnco0+1)
          read(14,*) xo(i)
          read(15,*) yo(i)
        end do

        close(14)
        close(15)
        
!-----------------------------------------------------------------------
! existing scroll bars
!-----------------------------------------------------------------------

! scroll bar pointers
        allocate ( pbrow(Nnco0+1), pbsum(Nnco0+1), stat=ierr )
        if(ierr.gt.0) go to 999    

        filename = TRIM(dirIN)//TRIM(simname)//'_pb_row.dat'
        open(unit=12, file=filename, status='old')
        filename = TRIM(dirIN)//TRIM(simname)//'_pb_sum.dat'
        open(unit=13, file=filename, status='old')

        do i = 1, Nnco0+1
          read(12,*) pbrow(i)
          read(13,*) pbsum(i)
        end do

        close(12)
        close(13) 

! scroll bar coordinates
        allocate ( xb(pbrow(Nnco0+1)), yb(pbrow(Nnco0+1)), stat=ierr )
        if(ierr.gt.0) go to 999     
 
        filename = TRIM(dirIN)//TRIM(simname)//'_xb.dat'
        open(unit=14, file=filename, status='old')
        filename = TRIM(dirIN)//TRIM(simname)//'_yb.dat'
        open(unit=15, file=filename, status='old')

        do i = 1, pbrow(Nnco0+1)
          read(14,*) xb(i)
          read(15,*) yb(i)
        end do

        close(14)
        close(15)

        write(6,*) 'number of scroll bars/oxbow lakes = ', Nnco0
        
!-----------------------------------------------------------------------
! flag error
      case default
        stop 'ERROR! Flag for existing floodplain structure'
      end select
      
! end of the subroutine
      return
999   write(6,*) 'ERROR! Allocation in loading floodplain structure'
      stop
      end
