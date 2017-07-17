      subroutine dashline(i)
      implicit none
      integer i
      write(i,*) '-----------------------------------------------------'
      return
      end

!      subroutine dashline_title(i, title)
!      implicit none
!      integer i
!      character*15 title

!      title = ADJUSTL(title)
!      write(6,*) '- '
!      write(6,'(a10)',ADVANCE = "NO") title
!      write(6,*) ' -'
!      return
!      end
