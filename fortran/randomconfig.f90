!-----------------------------------------------------------------------
! SUBROUTINES PER LA GENERAZIONE DELLA CONFIGURAZIONE INIZIALE
!----------------------------------------------------------------------- 

!----------------------------------------------------------------------- 
! Manuel
!----------------------------------------------------------------------- 
      subroutine randomconfig2 (N, Nrand, deltas, stdev, x, y)
      
      implicit none
      integer j, N, Nrand, Nini, Nend
      real*8 deltas, stdev, yy
      real*8 x(N), y(N)

! inizializzo i vettori delle coordinate
      x(:) = 0.d0
      y(:) = 0.d0

! costruisco il vettore delle coordinate longitudinali
      do j = 1, N
        x(j) = deltas * dble(j-1)
      end do

! calcolo l'intervallo dei punti perturbati
      Nini = INT(REAL(N)/2.d0 - REAL(Nrand)/2.d0 + 1.d0)
      Nend = INT(REAL(N)/2.d0 + REAL(Nrand)/2.d0)

! genero il vettore contenente i numeri estratti da una distribuzione
! normale con media 0 e deviazione standard dev N(0,dev)
      call randinitialize()
      do j = Nini, Nend
        call randnormnumberloop(0.d0, stdev, yy)
        y(j) = yy
      end do

      return
      end
