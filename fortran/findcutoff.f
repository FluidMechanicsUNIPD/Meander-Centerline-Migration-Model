c-----------------------------------------------------------------------
c REMOVE ABANDONED REACH AND STORE ITS COORDINATES
c-----------------------------------------------------------------------
      subroutine findcutoff(time,jmin,jmax,Nnew,jfile,jt)
      use module_global
      use module_geometry
      implicit none
      integer jmin, jmax, jco, k, j, jpup, jpdo, jt
      integer Nnew, Nrem, Ncb, Nco, jfile, ierr
      real*8 time, s1, s2
      real*8, allocatable :: xrem(:), yrem(:), xaux(:), yaux(:), baux(:)
      real*8, allocatable :: xco(:), yco(:), xcb(:), ycb(:) ! oxbow, point bar
      real*8, allocatable :: x1(:), y1(:), x2(:), y2(:) ! oxbow boundaries

!-----------------------------------------------------------------------
! summary file
!-----------------------------------------------------------------------

! print head (only first iteration)
      if (jt.eq.0) then
        write(jfile,102) 'j','Np','N','jur','juo','jdo','jdr',
     5                          'Nnew','time','jt'
        call dashline(jfile)
        write(jfile,103) 0,0,N,0,jmin,jmax,0,N,time,jt
        return
      end if

102   format(7x,a1,6x,a2,7x,a1,4(5x,a3),4x,a4,11x,a4,7x,a2)
103   format(8(1x,i7),1x,f14.5,1x,i8)

c-----------------------------------------------------------------------
c find abandoned reach
c-----------------------------------------------------------------------

c update neck cutoff counter
      Nnco = Nnco + 1

c number of points forming the abandoner reach
      Nrem = jmax - jmin + 1
                               
c coordinated of the abandoned reach
      allocate ( xrem(Nrem), yrem(Nrem) )
      do j = 1, Nrem
        xrem(j) = x(jmin+j-1)
        yrem(j) = y(jmin+j-1)
      end do

c-----------------------------------------------------------------------
c current oxbow lake
c-----------------------------------------------------------------------

c find oxbow left and right boundaries (i.e. banks of abandoned reach)
      allocate ( baux(Nrem), x1(Nrem), y1(Nrem), x2(Nrem), y2(Nrem) )

c set total width
      baux(:) = 2.d0

      call banks (Nrem, xrem, yrem, baux, x1, y1, x2, y2)
      deallocate (baux)

c number of points in the current owbow lake
      Nco = Nrem + Nrem + 1
      allocate (xco(Nco), yco(Nco))

c assemble oxbow coordinates and compute bank lenghts
      s1 = 0.d0
      s2 = 0.d0

      do j = 1, Nrem
        xco(j) = x1(j)
        yco(j) = y1(j)
        xco(Nrem + j) = x2(Nrem - j + 1)
        yco(Nrem + j) = y2(Nrem - j + 1)

        if (j.ge.2) then
          s1 = s1 + DSQRT((x1(j)-x1(j-1))**2.d0 + (y1(j)-y1(j-1))**2.d0)
          s2 = s2 + DSQRT((x2(j)-x2(j-1))**2.d0 + (y2(j)-y2(j-1))**2.d0)
        end if

      end do

c repeat first point
      xco(Nco) = x1(1)
      yco(Nco) = y1(1)

c-----------------------------------------------------------------------
c current point bar
c-----------------------------------------------------------------------

c number of points in the current point bar
      Ncb = Nrem + 1

c consider the shorter oxbow bank as point bar boundary
      allocate ( xcb(Ncb), ycb(Ncb) )
      do j = 1, Nrem
        if (s1.lt.s2) then
          xcb(j) = x1(j)
          ycb(j) = y1(j)
        else
          xcb(j) = x2(j)
          ycb(j) = y2(j)
        end if
c        xcb(j) = xrem(j)
c        ycb(j) = yrem(j)
      end do

c repeat first point
c      xcb(Ncb) = xrem(1)
c      ycb(Ncb) = yrem(1)
      if (s1.lt.s2) then
        xcb(Ncb) = x1(1)
        ycb(Ncb) = y1(1)
      else
        xcb(Ncb) = x2(1)
        ycb(Ncb) = y2(1)
      end if

c reset arrays
      deallocate (x1, y1, x2, y2, xrem, yrem)

c-----------------------------------------------------------------------
c print and store current cutoff (oxbow + point bar)
c-----------------------------------------------------------------------

c print on file the current cutoff (oxbow + point bar)
      if (jt.gt.0) then
        call printcutoff (Ncb, xcb, ycb, Nco, xco, yco, jt)
      end if

c store oxbow
      call storeoxbow(Nco, xco, yco)

c store point bar
      call storepointbar(Ncb, xcb, ycb)

c reset auxiliary arrays
      deallocate (xco, yco, xcb, ycb)

c-----------------------------------------------------------------------
c remove points from the active channel
c-----------------------------------------------------------------------

c last point upstream of the cutoff
      jpup = MAX(1,jmin - jre)

c first point downstream of the cutoff
      jpdo = MIN(N, jmax + jre)

c total number of removed points
      k = jpdo - jpup - 1

c new number of points forming the active river
      Nnew = N-k

c initialize auxiliary arrays
      allocate(xaux(Nnew), yaux(Nnew), stat=ierr )

c cutoff upstream points do not move
      do j = 1, jpup
        xaux(j) = x(j)
        yaux(j) = y(j)
      end do

c cutoff-downstream points are shifted backward
      do j = jpup+1, Nnew
        xaux(j) = x(j+k)
        yaux(j) = y(j+k)
      end do

c initialize river coordinates
      deallocate(x,y)
      allocate(x(Nnew),y(Nnew))

c restore pristine river coordinates
      x(:) = xaux(:)
      y(:) = yaux(:)

c reset auxiliary arrays
      deallocate(xaux,yaux)

! print quantities of current iteration
      write(jfile,103) Nnco,Nrem,N,jpup,jmin,jmax,jpdo,Nnew,time,jt

c-----------------------------------------------------------------------
c end of subroutine
c-----------------------------------------------------------------------
      return
      end
