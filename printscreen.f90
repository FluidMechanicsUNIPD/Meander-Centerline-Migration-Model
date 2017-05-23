!-----------------------------------------------------------------------
! PRINT OUTPUT FILES
!-----------------------------------------------------------------------
      subroutine printscreen(jt, jp, time, Nnco)
      use module_global
      implicit none
      integer jt, jp, Nnco
      real*8 time
      
!----------------------------------------------------------------------- 
! print on screen
!-----------------------------------------------------------------------  

! print head (only first iteration)
      if (jt.eq.0) then
        call dashline(6)
        write(6,102) 'jt', 'time', 'N','Nnco','Cf', 'jp'
        call dashline(6)
      end if
 
! print quantities of current iteration
      write(6,103) jt, time, N, Nnco, Cf, jp!, beta
     
!-----------------------------------------------------------------------
! end of subroutine
!-----------------------------------------------------------------------
102   format(7x,a2,7x,a4,4x,a2,2x,a4,8x,a2,5x,a2)
103   format(1x,i8,1x,f10.2,2(1x,i5),1x,f9.5,1x,i6,f12.5)

      return
      end
