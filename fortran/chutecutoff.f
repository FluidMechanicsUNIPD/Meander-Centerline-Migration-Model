!-----------------------------------------------------------------------
! FIND CHUTE CUTOFFS
!-----------------------------------------------------------------------
      subroutine chutecutoffs
      use variables
      implicit none

      integer i, j, k, jup, jdw, Nnew, counter, ierr
      integer flag(N)
      real*8 dist
      real*8, allocatable :: vec1(:), vec2(:)

!-----------------------------------------------------------------------
! initialize pointer
      flag(:) = 1
      counter = 0
      
! loop from the first point to the (n-jcco)th one
      do while (i.lt.(N - jcco))
        i = i + 1
        jdw = i

! loop for downstream closed points
        do j = i+jcco, MIN(i+jnco-1,N)

! evaluate distance between points
          dist = ((x(i)-x(j))**2.d0 + (y(i)-y(j))**2.d0)**0.5d0

          if (dist.lt.tollc) then
            jup = jdw+1
            jdw = j-1
            i = j

! delete all points upstream
            do k = jup, jdw
              if (flag(k).ne.0) then
                flag(k) = 0
                counter = counter + 1
              end if
            end do

          end if

! end loop over points
        end do
      end do

!-----------------------------------------------------------------------
! delete chute cutoffs
      if (counter.gt.0) then

! update point number
        Nnew = N - counter

! set auxiliary arrays
        allocate( vec1(N), vec2(N), stat=ierr )
        if(ierr.gt.0) stop 'ERROR! Allocation in chute cutoffs'
        do i = 1, N
          vec1(i) = x(i)
          vec2(i) = y(i)
        end do

        deallocate(x,y)
        allocate(x(Nnew),y(Nnew))   

! transfer coordinates in pristine arrays
        j = 0
        do i = 1, N
          if(flag(i).ne.0) then
            j = j+1
            x(j) = vec1(i)
            y(j) = vec2(i)
            
! update chute cutoff counter
            if ((i.gt.1).and.(flag(i).eq.0).and.(flag(i-1).ne.1)) then
              Ncco = Ncco + 1
            end if
            
          end if
        end do

! reset auxiliary arrays
        deallocate(vec1,vec2)

! restore new point number
        N = Nnew
      end if
         
c-----------------------------------------------------------------------
c end of subroutine
c-----------------------------------------------------------------------
      return
      end          
