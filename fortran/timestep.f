!-----------------------------------------------------------------------
! TIME STEP
!-----------------------------------------------------------------------
      real*8 function timestep(N,Ef,deltas,dU,cstab,tollc,dt0,flag_dt)
      implicit none
      integer N, flag_dt
      real*8 maxCSI, Ef, deltas, dt0, tollc, cstab, dtmax, year2sec
      real*8 dU(N)
      parameter ( year2sec = 86400.d0 * 365.d0 )

! fixed time step
      if (flag_dt.eq.1) then
        timestep = dt0
        
! dynamic time step
      else if (flag_dt.eq.2) then

! find maximim value of migration rate
        maxCSI = Ef * MAXVAL(dU)
      
! find time step threshold (dimensionless years)
        dtmax = cstab * MIN(deltas,tollc)/(maxCSI) / year2sec
      
! set current time step
        timestep = MIN(dt0, dtmax)
        
! wrong flag
      else
        stop 'ERROR! Wrong flag for time step'
      end if

! end of function
      return
      end
