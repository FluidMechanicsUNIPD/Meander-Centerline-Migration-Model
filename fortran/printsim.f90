!-----------------------------------------------------------------------
! PRINT OUTPUT FILES
!-----------------------------------------------------------------------
      subroutine printsim(jt, jp, time, jfile1, jfile2, trasp, fondo)
      use module_global
      use module_geometry
      implicit none
      character*20 print1
      integer j, jt, jp, jfile1, jfile2, trasp, fondo
      real*8 time
      
!-----------------------------------------------------------------------
! print current configuration
!-----------------------------------------------------------------------

! file name
      write( print1, '(i10)' )  jp
      print1 = ADJUSTL(print1)
      filename = TRIM(dirOUT)//'configuration_'//TRIM(print1)//'.out'
      
! open file
      open(unit=21, file=filename, status='replace')
      
! print current configuration
      do j = 1, N
        write(21,101) x(j), y(j), s(j), th(j), c(j), dU(j)
      end do
      
! close file
      close(21) 

!-----------------------------------------------------------------------
! summary file
!-----------------------------------------------------------------------  

! print head (only first iteration)
      if (jt.eq.0) then
        write(jfile1,102) 'jp', 'jt', 'time', 'Np', 'Npb', 'Npo', 'Nnco'
        call dashline(jfile1)
      end if
      
! print quantities of current iteration
      write(jfile1,103) jp, jt, time, N, Npb, Npo, Nnco

! print head (only first iteration)
      if (jt.eq.0) then
        write(jfile2,104) 'jp','jt','time','L','beta','theta','ds', &
     &      'Cf','bed','load'
        call dashline(jfile2)
      end if
      
! print quantities of current iteration
      write(jfile2,105) jp,jt,time,Ls,beta,theta,ds,Cf,fondo,trasp
     
!-----------------------------------------------------------------------
! end of subroutine
!-----------------------------------------------------------------------
101   format(*(f16.8,1x))
102   format(4x,a2,7x,a2,11x,a4,5x,a2,2(4x,a3),3x,a4)
103   format(i6,1x,i8,1x,f14.4,4(1x,i6))
104   format(4x,a2,7x,a2,11x,a4,10x,a1,7x,a4,6x,a5,9x,a2,9x,a2,2x, &
     &                      a3,1x,a4)
105   format(i6,1x,i8,1x,f14.4,1x,f10.3,4(1x,f10.6),2(4x,i1))

      return
      end
