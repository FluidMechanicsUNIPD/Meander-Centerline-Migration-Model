c-----------------------------------------------------------------------
c STORE THE NEW OXBOW INTO GLOBAL STORAGE ARRAYS
c-----------------------------------------------------------------------
      subroutine storeoxbow(Nco, xco, yco)
      use module_geometry
      implicit none
      integer i, j, k, l, Nco, ierr
      real*8 xco(Nco), yco(Nco)
      integer, allocatable :: aux3(:), aux4(:)
      real*8, allocatable :: aux1(:), aux2(:)

c-----------------------------------------------------------------------
c pointer arrays
c-----------------------------------------------------------------------

c if pointer vectors oxsum e oxrow do not exist, then they are created
      if (Nnco.eq.1) then

        allocate (oxsum(Nnco+1), oxrow(Nnco+1) , stat=ierr )
        if(ierr.gt.0) go to 999

        oxsum(1) = 0
        oxrow(1) = 0
		  
c if there are already oxbow, then update vector size
      else if (Nnco.gt.1) then

c transef existing vectors on auxiliary arrays, reset and resize them,
c and finally restore them
        allocate (aux3(Nnco), aux4(Nnco))

        aux3(:) = oxsum(:)
        aux4(:) = oxrow(:)

        deallocate (oxsum, oxrow)
        allocate (oxsum(Nnco+1), oxrow(Nnco+1) , stat=ierr )
        if(ierr.gt.0) go to 999

        do j = 1, Nnco
          oxsum(j) = aux3(j)
          oxrow(j) = aux4(j)
        end do

        deallocate (aux3, aux4) 

      end if
  
c store the point number of current oxbow
      oxsum(Nnco+1) = Nco
      oxrow(Nnco+1) = oxrow(Nnco) + Nco

c-----------------------------------------------------------------------
c update storage vectors for coordinates
c-----------------------------------------------------------------------

c if coordinate vectors xo e yo do not exist, then they are created
      if (Nnco.eq.1) then
        k = 0

        allocate(xo(Nco), yo(Nco), stat=ierr)
        if(ierr.gt.0) go to 999 

c if there are already oxbows, then update vector size
      else if (Nnco.gt.1) then

c oxrow(Nnco) = number of the row where last cutoff ends
        k = oxrow(Nnco)

c transfer existing vectors on auxiliary arrays, reset and resize them,
c and finally restore them
        allocate (aux1(k), aux2(k))

        aux1(:) = xo(:)
        aux2(:) = yo(:)

        deallocate (xo, yo)
        allocate (xo(k+Nco),yo(k+Nco), stat=ierr)
        if(ierr.gt.0) go to 999 

        do l = 1, k
          xo(l) = aux1(l)
          yo(l) = aux2(l)
        end do

        deallocate (aux1, aux2)

      end if 

c-----------------------------------------------------------------------
c storage of oxbow coordinates
c-----------------------------------------------------------------------

c loop over new oxbow points
      do j = 1, Nco

c store coordinates
        xo(k+j) = xco(j)
        yo(k+j) = yco(j)     
	
c end loop over new oxbow points
      end do

c-----------------------------------------------------------------------
c end of subroutine
c-----------------------------------------------------------------------
      return
 999  write(6,*) 'ERROR! Allocation in storing oxbow lakes'
      stop
      end