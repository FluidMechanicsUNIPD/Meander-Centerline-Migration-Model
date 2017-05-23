!-----------------------------------------------------------------------
!  Module containing the variables related to geometry
!-----------------------------------------------------------------------
      module module_geometry
      implicit none
      
      integer Nnco, Nnco0, Nncoold     ! number of neck cutoff (current, initial, previous)
      integer Npf, Npb, Npo   ! number of points in floodplain, point bars, and oxbow lakes
      integer Nsamp, Nsampold ! number of resampling processes (current, previous)
      real*8 Ls, Ls0, Lsold   ! length of the river (current, initial, previous)   
      real*8 Lx, Lx0, Lxold   ! projected length of the river (current, initial, previous)

! oxbow lakes and scroll bars
      integer, allocatable :: oxrow(:), oxsum(:)      ! oxbow pointers
      integer, allocatable :: pbrow(:), pbsum(:)      ! point bar pointers
      integer, allocatable :: whereisnode(:)          ! position of node with respect to oxbows or point bars
      
! river
      real*8, allocatable :: x(:), y(:), x0(:), y0(:) ! coordinates of axis grid (initial)
      real*8, allocatable :: xo(:), yo(:)             ! coordinates of oxbow lakes
      real*8, allocatable :: xb(:), yb(:)             ! coordinates of point bars
      real*8, allocatable :: s(:)                     ! longitudinal axis coordinate
      real*8, allocatable :: th(:)                    ! axis angle with respect to the longitudinal direction
      real*8, allocatable :: c(:)                     ! axis curvature (multiplied by v0)
      real*8, allocatable :: dU(:)                    ! excess near-bank longitudinal velocity      ex deltaU
      real*8, allocatable :: E(:)                     ! erodibility arrays
          
!-----------------------------------------------------------------------
      end module module_geometry
!-----------------------------------------------------------------------
