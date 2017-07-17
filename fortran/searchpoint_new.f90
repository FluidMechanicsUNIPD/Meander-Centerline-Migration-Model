!-----------------------------------------------------------------------
! POSITION OF THE AXIS POINTS WITH RESPECT TO THE GEOMORPHIC UNITS
!-----------------------------------------------------------------------
! whereisnode indicates the position of the current point (xp,yp), and is
! both an input and and output. First it is used as input to determine 
! whether the point lies on the same geomorphic unit (floodplain, point 
! bar, oxbow lake) as previously. If not, then it is updated.
! whereisnode=0 means the point lies into the floodplain.
! whereisnode>0 means the point lies in the point bar with index whereisnode
! whereisnode<0 means the point lies in the oxbow lake with index -whereisnode
!-----------------------------------------------------------------------

      subroutine searchpointnew (js)
      use module_global
      use module_geometry
      implicit none
      integer i, j, js, jup, jdw, naux, wn, windingnumber
      real*8 xp, yp
      real*8, allocatable :: xaux(:), yaux(:)

! set current point
      xp = x(js)
      yp = y(js)

!-----------------------------------------------------------------------
! check if the point lies in the same point bar as previously
      if (whereisnode(js).gt.0) then
!-----------------------------------------------------------------------

! set initial and last point of the point bar
        jup = pbrow(whereisnode(js))
        jdw = pbrow(whereisnode(js)+1)
        naux = jdw - jup

! allocate auxiliary arrays for point bar coordinates
        allocate ( xaux(naux), yaux(naux) )

! set current point bar
        do i = 1, naux
          xaux(i) = xb(jup+i)
          yaux(i) = yb(jup+i)
        end do

! check for point position, if wn=/=0 the point is inside the point bar
        if ( (xp.ge.MINVAL(xaux)).and.(xp.le.MAXVAL(xaux)).and. &
     &       (yp.ge.MINVAL(yaux)).and.(yp.le.MAXVAL(yaux)) ) then
          wn = windingnumber(xp, yp, naux, xaux, yaux)
        end if

! reset auxiliary coordinates arrays
        deallocate (xaux, yaux)

! if the winding number is not zero, then the point remains into the
! point bar which was hosting it
        if (wn.ne.0) then
          return
        end if

!-----------------------------------------------------------------------
! check if the point lies in the same oxbow lake as previously
      else if (whereisnode(js).lt.0) then
!-----------------------------------------------------------------------

! set initial and last point of the oxbow lake
        jup = oxrow(ABS(whereisnode(js)))
        jdw = oxrow(ABS(whereisnode(js))+1)
        naux = jdw - jup

! allocate auxiliary arrays for point bar coordinates
        allocate ( xaux(naux), yaux(naux) )

! set current oxbow lake
        do i = 1, naux
          xaux(i) = xo(jup+i)
          yaux(i) = yo(jup+i)
        end do

! check for point position, if wn=/=0 the point is inside the oxbow lake
        if ( (xp.ge.MINVAL(xaux)).and.(xp.le.MAXVAL(xaux)).and. &
     &       (yp.ge.MINVAL(yaux)).and.(yp.le.MAXVAL(yaux)) ) then
          wn = windingnumber(xp, yp, naux, xaux, yaux)
        end if

! reset auxiliary coordinates arrays
        deallocate (xaux, yaux)

! if the winding number is not zero, then the point remains into the
! oxbow lake which was hosting it
        if (wn.ne.0) then
          return
        end if

! end if for existing position of the point
      end if

!-----------------------------------------------------------------------
! search the eventual oxbow lake or point bar where the point lies
!-----------------------------------------------------------------------

! initialize position pointer
      whereisnode(js) = 0

! loop for all oxbow lakes and point bars alternatively, starting from 
! last ones, i.e. more recent
      j = Nnco

      do while ((j.ge.1).and.(whereisnode(js).eq.0))

! set initial and last point of the current oxbow lake
        jup = oxrow(j)
        jdw = oxrow(j+1)
        naux = jdw - jup

! allocate arrays for oxbow coordinates
        allocate ( xaux(naux), yaux(naux) )

! set current oxbow lake
        do i = 1, naux
          xaux(i) = xo(jup+i)
          yaux(i) = yo(jup+i)
        end do

! check for point position, if wn=/=0 the point is inside the oxbow lake
        if ( (xp.ge.MINVAL(xaux)).and.(xp.le.MAXVAL(xaux)).and. &
     &       (yp.ge.MINVAL(yaux)).and.(yp.le.MAXVAL(yaux)) ) then
          wn = windingnumber(xp, yp, naux, xaux, yaux)
          if (wn.ne.0) then
            whereisnode(js) = - j  ! NB whereisnode<0 for oxbow lake lying
          end if
        end if

! reset auxiliary arrays
        deallocate (xaux, yaux)

        if (whereisnode(js).eq.0) then

! set initial and last point of the current point bar
          jup = pbrow(j)
          jdw = pbrow(j+1)
          naux = jdw - jup

! allocate arrays for point bar coordinates
          allocate ( xaux(naux), yaux(naux) )

! set current point bar
          do i = 1, naux
            xaux(i) = xb(jup+i)
            yaux(i) = yb(jup+i)
          end do

! check for point position, if wn=/=0 the point is inside the point bar
          if ( (xp.ge.MINVAL(xaux)).and.(xp.le.MAXVAL(xaux)).and. &
     &         (yp.ge.MINVAL(yaux)).and.(yp.le.MAXVAL(yaux)) ) then
            wn = windingnumber(xp, yp, naux, xaux, yaux)
            if (wn.ne.0) then
              whereisnode(js) = + j    ! NB whereisnode>0 for point bar lying
            end if
          end if

! reset auxiliary arrays
          deallocate (xaux, yaux)

        end if

! update oxbow & point bar pointer
        j = j - 1

! end loop over oxbow lakes and point bars
      end do

!-----------------------------------------------------------------------
! end of subroutine
!-----------------------------------------------------------------------
      return
      end
