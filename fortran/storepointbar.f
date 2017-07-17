c-----------------------------------------------------------------------
c STORE THE NEW SCROLL BAR INTO GLOBAL STORAGE ARRAYS
c-----------------------------------------------------------------------
      subroutine storepointbar(Ncb, xcb, ycb)
      use module_geometry
      implicit none
      integer i, j, k, l, Ncb, ierr
      real*8 xcb(Ncb), ycb(Ncb)
      integer, allocatable :: aux3(:), aux4(:)
      real*8, allocatable :: aux1(:), aux2(:)

c-----------------------------------------------------------------------
c pointer arrays
c-----------------------------------------------------------------------

c if pointer vectors pbsum e pbrow do not exist, then they are created
      if (Nnco.eq.1) then

        allocate (pbsum(Nnco+1), pbrow(Nnco+1) , stat=ierr )
        if(ierr.gt.0) go to 999

        pbsum(1) = 0
        pbrow(1) = 0

c if there are already point bars, then update vector size
      else if (Nnco.gt.1) then

c transef existing vectors on auxiliary arrays, reset and resize them,
c and finally restore them
        allocate (aux3(Nnco), aux4(Nnco))

        aux3(:) = pbsum(:)
        aux4(:) = pbrow(:)

        deallocate (pbsum, pbrow)
        allocate (pbsum(Nnco+1), pbrow(Nnco+1) , stat=ierr )
        if(ierr.gt.0) go to 999

        do j = 1, Nnco
          pbsum(j) = aux3(j)
          pbrow(j) = aux4(j)
        end do

        deallocate (aux3, aux4)

      end if

c store the point number of current point bar
      pbsum(Nnco+1) = Ncb
      pbrow(Nnco+1) = pbrow(Nnco) + Ncb

c-----------------------------------------------------------------------
c update storage vectors for coordinates
c-----------------------------------------------------------------------

c if coordinate vectors xb e yb do not exist, then they are created
      if (Nnco.eq.1) then
        k = 0

        allocate(xb(Ncb), yb(Ncb), stat=ierr)
        if(ierr.gt.0) go to 999

c if there are already point bars, then update vector size
      else if (Nnco.gt.1) then

c pbrow(Nnco) = number of the row where last cutoff ends
        k = pbrow(Nnco)

c transfer existing vectors on auxiliary arrays, reset and resize them,
c and finally restore them
        allocate (aux1(k), aux2(k))

        aux1(:) = xb(:)
        aux2(:) = yb(:)

        deallocate (xb, yb)
        allocate (xb(k+Ncb),yb(k+Ncb), stat=ierr)
        if(ierr.gt.0) go to 999

        do l = 1, k
          xb(l) = aux1(l)
          yb(l) = aux2(l)
        end do

        deallocate (aux1, aux2)

      end if

c-----------------------------------------------------------------------
c storage of point bar coordinates
c-----------------------------------------------------------------------

c loop over new point bar points
      do j = 1, Ncb

c store coordinates
        xb(k+j) = xcb(j)
        yb(k+j) = ycb(j)

c end loop over new point bar points
      end do

c-----------------------------------------------------------------------
c end of subroutine
c-----------------------------------------------------------------------
      return
 999  write(6,*) 'ERROR! Allocation in storing point bars'
      stop
      end
