c-----------------------------------------------------------------------

	subroutine kappa_alpha(Cf,k0,k1,k3,k4)

      implicit none
	character*50 nome
      integer j, NCf, jdown, jup, ierr
      real*8 Cf, k0, k1, k3, k4, maxCf, minCf, deltaCf, resto
      real*8 Cfupdown, wk0updown, wk1updown, wk3updown, wk4updown
      real*8 wmk0, wmk1, wmk3, wmk4
      real*8, allocatable :: a1(:), a2(:), a3(:), a4(:), a5(:)
      real*8, allocatable :: a6(:), a7(:), a8(:), a9(:), Cfv(:)
      real*8, allocatable :: k0v(:),k1v(:),k2v(:),k3v(:),k4v(:),k5v(:)

c-----------------------------------------------------------------------
c apertura file di libreria con i coefficienti
c-----------------------------------------------------------------------
      nome="input/kappa_alpha.dat"
      open(unit=8, file=nome, status='old')
      read(8,*) NCf

      allocate ( a1(Ncf), a2(Ncf), a3(Ncf), a4(Ncf), a5(Ncf), 
     &           a6(Ncf), a7(Ncf), a8(Ncf), a9(Ncf), stat=ierr )
      if(ierr.gt.0) stop 'allocation error: ZS coefficients'

      do j = 1, NCf
	  read(8,*) a1(j),a2(j),a3(j),a4(j),a5(j),a6(j),a7(j),a8(j),a9(j)
      end do
      close(8)

c-----------------------------------------------------------------------
      allocate (Cfv(NCf),k0v(Ncf),k1v(Ncf),k3v(Ncf),k4v(NCf),stat=ierr)
      if(ierr.gt.0) stop 'allocation error: ZS coefficients'

	Cfv = a1
	k0v = a2
	k1v = a3
	k3v = a6
	k4v = a7

	maxCf = MAXVAL(Cfv)
	minCf = MINVAL(Cfv)
	deltaCf = (maxCf-minCf) / (NCf-1.d0)

	resto = DMOD((Cf-minCf),deltaCf)

c-----------------------------------------------------------------------
c calcolo i coefficienti k0, k1, k2, k3, k4
c-----------------------------------------------------------------------
	if (resto.eq.0.d0) then
        j = (Cf - minCf) / deltaCf + 1
	  k0 = k0v(j)
	  k1 = k1v(j)
	  k3 = k3v(j)
	  k4 = k4v(j)
	else 
	  jdown = FLOOR((Cf-minCf)/deltaCf) + 1
c      write(6,'(4(e12.6,1x),i2)')Cf,minCf,deltaCf,(Cf-minCf)/deltaCf

	  jup = jdown + 1
	  wk0updown = k0v(jup) - k0v(jdown)
	  wk1updown = k1v(jup) - k1v(jdown)
     	  wk3updown = k3v(jup) - k3v(jdown)
	  wk4updown = k4v(jup) - k4v(jdown)  
      
	  Cfupdown=Cfv(jup)-Cfv(jdown)
	  wmk0 = wk0updown/Cfupdown
	  wmk1 = wk1updown/Cfupdown
	  wmk3 = wk3updown/Cfupdown
	  wmk4 = wk4updown/Cfupdown
    
	  k0 = k0v(jdown) + wmk0 * (Cf-Cfv(jdown))
        k1 = k1v(jdown) + wmk1 * (Cf-Cfv(jdown))
	  k3 = k3v(jdown) + wmk3 * (Cf-Cfv(jdown))
	  k4 = k4v(jdown) + wmk4 * (Cf-Cfv(jdown))
	end if
c-----------------------------------------------------------------------
      deallocate ( a1, a2, a3, a4, a5, a6, a7, a8, a9, stat=ierr )
      if(ierr.gt.0) stop 'allocation error: ZS coefficients'
      deallocate ( Cfv, k0v, k1v, k3v, k4v, stat=ierr)
      if(ierr.gt.0) stop 'allocation error: ZS coefficients'

c-----------------------------------------------------------------------
	return
	end
