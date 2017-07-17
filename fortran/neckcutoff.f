c-----------------------------------------------------------------------
c NECK CUTOFF SEARCHING
c-----------------------------------------------------------------------
      subroutine neckcutoff(time,jfile,jt)

      use module_global
      use module_geometry
      implicit none

      integer i, j, k, jt, jb, jk, js, jmin, jmax, jpp, ko, jfile, ierr
      integer nxmin, nxmax, nymin, nymax, Nr, Nc, NrB, nn, mm, ncell
      integer Ngmax, ncellini, ncellfin, resto, Nnew
      real*8 time, xmin, xmax, ymin, ymax
      real*8 xgridmin, xgridmax, ygridmin, ygridmax, dxy
      integer, allocatable :: GRID(:,:), indpos(:)
      
c-----------------------------------------------------------------------
c structure of floodplain grid
c-----------------------------------------------------------------------

c extreme coordinates
60    xmin = MINVAL(x)
      xmax = MAXVAL(x)
      ymin = MINVAL(y)
      ymax = MAXVAL(y)

c set grid cells
      nxmin = FLOOR(xmin/tollc) - 1
      nxmax = CEILING(xmax/tollc) + 1
      nymin = FLOOR(ymin/tollc) - 1
      nymax = CEILING(ymax/tollc) + 1

c row and column numbers
      Nr = ABS(nymax-nymin)
      Nc = ABS(nxmax-nxmin)

c coordinates of grid edges
      xgridmin = DBLE(nxmin) * tollc
      xgridmax = DBLE(nxmax) * tollc
      ygridmin = DBLE(nymin) * tollc
      ygridmax = DBLE(nymax) * tollc

c total number of grid cells
      NrB = Nr * Nc

c maximum number of point in a grid cell
      Ngmax = INT( SQRT(2.D0) * tollc/deltas) + 1

c build grid matrix
      if (ALLOCATED(GRID)) deallocate(GRID, indpos)
      allocate ( GRID(NrB,Ngmax+6) , indpos(NrB), stat=ierr )
      if (ierr.gt.0) stop 'ERROR! Allocation in neck cutoff'
      GRID(:,:) = 0
      indpos(:) = 0

c-----------------------------------------------------------------------
c store river point in floodplain grid
c-----------------------------------------------------------------------

c loop over points
      do i = 1, N
        nn = FLOOR( (ABS(x(i)-xgridmin)) / tollc) + 1 ! column index of the ith point from left
        mm = FLOOR( (ABS(ygridmax-y(i))) / tollc) + 1 ! row index of the ith point from above
        ncell = Nc * (mm-1) + nn                     ! row of GRID cointaining this point  

c the kth row of GRID contains the indices of the points within the 
c kth cell
        indpos(ncell) = indpos(ncell) + 1
        GRID(ncell,indpos(ncell)) = i

c end loop over ponts
      end do

c-----------------------------------------------------------------------
c look for neck cutoffs
c-----------------------------------------------------------------------

c the first row (1<ncell<Nc) as well as the last row (NrB-Nc+1<ncell<NrB)
c of the grid are not considered; the first (Nc+1) and the last (NrB-Nc)
c columns of the second and the second-last rows are neglected; summarizing, 
c the first element is the second of the second row while the last is the
c second-last of the second-last row.
      ncellini = Nc + 2
      ncellfin = NrB-Nc-1 

      do ncell = ncellini, ncellfin

c resto = 1 : first column of GRID
c resto > 1 : from second column of GRID
        resto = mod(ncell,Nc)

c first and last column of GRID are not considered
        if (resto.gt.1) then 
        
c check if the current cell cointains nodes
          if (GRID(ncell,1).ne.0) then
          
c-----------------------------------------------------------------------
c 1) looking for cutoffs into the current cell
c-----------------------------------------------------------------------
            jb = 1

            do while ((GRID(ncell,jb+1).ne.0).and.(jb.le.(Ngmax+4)))

c indices of the two nodes
              jk = GRID(ncell,jb+1)
              js = GRID(ncell,jb)

              if (abs(js-jk).ge.jnco) then

c distance between the two nodes
                dxy = DSQRT((x(jk)-x(js))**2.d0 + (y(jk)-y(js))**2.d0)
                         
c if the distance if lower than the threshold, then a cutoff occurs    
                if (dxy.lt.tollc) then
                  jmin = min(js,jk)
                  jmax = max(js,jk)
                  call findcutoff(time,jmin,jmax,Nnew,jfile,jt)

c update the point number                  
                  N = Nnew 
                  ncellini = ncell+1 
                  goto 60
                end if

              end if

c go to the next point within the considered cell
              jb = jb + 1
            end do

c-----------------------------------------------------------------------
c 2) looking for cutoffs into the next cell ncell+1
c-----------------------------------------------------------------------
            jpp = 1

            if (GRID(ncell+1,jpp).ne.0) then
              jb = 1  
              do while ((GRID(ncell,jb).ne.0).and.(jb.le.(Ngmax+4)))

c indices of the two nodes
                js = GRID(ncell+1,jpp)
                jk = GRID(ncell,jb)               

                if (abs(js-jk).ge.jnco) then

c distance between the two nodes
                  dxy = DSQRT((x(jk)-x(js))**2.d0 + (y(jk)-y(js))**2.d0)
              
c if the distance if lower than the threshold, then a cutoff occurs   
                  if (dxy.lt.tollc) then
                    jmin = min(js,jk)
                    jmax = max(js,jk)
                    call findcutoff(time,jmin,jmax,Nnew,jfile,jt)            

c update the point number                    
                    N = Nnew 
                    ncellini = ncell+1 
                    goto 60
                  end if
                end if

c go to the next point within the considered cell
                jb = jb + 1
              end do
            end if 

c-----------------------------------------------------------------------
c 3) looking for cutoffs into the previous cell ncell-1
c-----------------------------------------------------------------------
            jpp = 1

            if (GRID(ncell-1,jpp).ne.0) then
              jb = 1                            
              do while ((GRID(ncell,jb).ne.0).and.(jb.le.(Ngmax+4)))

c indices of the two nodes
                js = GRID(ncell-1,jpp)
                jk = GRID(ncell,jb)    

                if (abs(js-jk).ge.jnco) then

c distance between the two nodes
                  dxy = DSQRT((x(jk)-x(js))**2.d0 + (y(jk)-y(js))**2.d0)

c if the distance if lower than the threshold, then a cutoff occurs                              
                  if (dxy.lt.tollc) then
                    jmin = min(js,jk)
                    jmax = max(js,jk)
                    call findcutoff(time,jmin,jmax,Nnew,jfile,jt)           

c update the point number                    
                    N = Nnew 
                    ncellini = ncell+1 
                    goto 60
                  end if
                end if

c go to the next point within the considered cell
                jb = jb + 1
              end do
            end if
 
c-----------------------------------------------------------------------
c 4) looking for cutoffs into the three cells above the current one
c-----------------------------------------------------------------------                     
            do ko = -1, 1                     
              jpp = 1

              if (GRID(ncell-Nc+ko,jpp).ne.0) then
                jb = 1  
                do while ((GRID(ncell,jb).ne.0).and.(jb.le.(Ngmax+4)))

c indices of the two nodes
                  js = GRID(ncell-Nc+ko,jpp)
                  jk = GRID(ncell,jb)

                  if(abs(js-jk).ge.jnco) then

c distance between the two nodes
                    dxy = DSQRT((x(jk)-x(js))**2.d0+(y(jk)-y(js))**2.d0)

c if the distane if lower than the threshold, then a cutoff occurs                             
                    if (dxy.lt.tollc) then
                      jmin = min(js,jk)
                      jmax = max(js,jk)
                      call findcutoff(time,jmin,jmax,Nnew,jfile,jt)            

c update the point number                      
                      N = Nnew 
                      ncellini = ncell+1 
                      goto 60
                    end if
                  end if

c go to the next point within the considered cell
                  jb = jb + 1
                end do
              end if
            end do

c-----------------------------------------------------------------------
c 5) looking for cutoffs into the three cells below the current one
c-----------------------------------------------------------------------
            do ko = -1, 1
              jpp = 1
                   
              if (GRID(ncell+Nc+ko,jpp).ne.0) then
                jb = 1  
                do while ((GRID(ncell,jb).ne.0).and.(jb.le.(Ngmax+4))) !jb+1

c indices of the two nodes
                  js = GRID(ncell+Nc+ko,jpp)
                  jk = GRID(ncell,jb)

                  if(abs(js-jk).ge.jnco) then

c distance between the two nodes
                    dxy = DSQRT((x(jk)-x(js))**2.d0+(y(jk)-y(js))**2.d0) 

c if the distance if lower than the threshold, then a cutoff occurs                            
                    if (dxy.lt.tollc) then
                      jmin = min(js,jk)
                      jmax = max(js,jk)
                      call findcutoff(time,jmin,jmax,Nnew,jfile,jt)

c update the point number                      
                      N = Nnew 
                      ncellini = ncell+1 
                    goto 60
                    end if
                  end if

c go to the next point within the considered cell
                  jb = jb + 1
                end do
              end if
            end do

c-----------------------------------------------------------------------
          end if ! end of checking points into a cellfine controlli presenza elementi in una cella
        end if   ! end of checking inside GRID
      end do     ! end loop over cells
c-----------------------------------------------------------------------

c reset grid
      deallocate (GRID, indpos)
      
c end of subroutine
      return
28    format(1x,i9,1x,e12.5,1x,f12.5,4(4x,i9),2(1x,f12.5))
      end
