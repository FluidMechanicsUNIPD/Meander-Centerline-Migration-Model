!-----------------------------------------------------------------------
! PRINT CHARACTERISTIC ENVALUES
!-----------------------------------------------------------------------
      subroutine printlambda(jt, jp, time, jfile)
      use module_global
      use module_zs
      implicit none
      integer jt, jp, jm, jfile
      real*8 time

! print head (only first iteration)
      if (jt.eq.0) then
        if (jmodel.eq.1) then
          write(jfile,101) 'jp', 'jt', 'jm', 'time', &
&                  'lambda1', 'lambda2', 'lambda3', 'lambda4'
        else if (jmodel.eq.2) then
          write(jfile,102) 'jp', 'jt', 'time', 'lambda'
        end if
        call dashline(jfile)
      end if

! print quantities of current iteration
      if (jmodel.eq.1) then
        do jm = 1, Mdat
          write(jfile,103) jp, jt, jm, time, &
&            DBLE(lamb1(jm)), IMAG(lamb1(jm)), &
&            DBLE(lamb2(jm)), IMAG(lamb2(jm)), &
&            DBLE(lamb3(jm)), IMAG(lamb3(jm)), &
&            DBLE(lamb4(jm)), IMAG(lamb4(jm))
          end do
      else if (jmodel.eq.2) then
        write(jfile,104) jp, jt, time, -2.d0*beta*Cf
      end if

!-----------------------------------------------------------------------
! end of subroutine
!-----------------------------------------------------------------------
101   format(4x,a2,7x,a2,1x,a2,11x,a4,4x,*(a7,15x))
102   format(4x,a2,7x,a2,11x,a4,5x,a6)
103   format(i6,1x,i8,1x,i2,1x,f14.4,*(1x,f10.6))
104   format(i6,1x,i8,1x,f14.4,*(1x,f10.6))

      return
      end
