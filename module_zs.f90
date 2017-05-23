!-----------------------------------------------------------------------
!  Module containing the variables related to the ZS model
!-----------------------------------------------------------------------
      module module_zs
      implicit none
      
      integer Nz     ! number of vertical points for the integration of the flow field
      integer Mdat   ! order of Fourier expansion
      complex*16, allocatable :: Am(:)
      complex*16, allocatable :: lamb1(:), lamb2(:), lamb3(:), lamb4(:)
      complex*16, allocatable :: g10(:), g20(:), g30(:), g40(:)
      complex*16, allocatable :: g11(:), g21(:), g31(:), g41(:)
      
!-----------------------------------------------------------------------
      end module module_zs
!-----------------------------------------------------------------------