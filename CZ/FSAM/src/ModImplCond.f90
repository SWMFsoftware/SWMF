module ModImplCond
  use ModPar,  ONLY: in, jn, kn
  implicit none
  
  public :: init_impl_cond
  
  real, allocatable :: MAT(:,:,:,:), res(:,:,:), rhs(:), dw(:), HeatCond(:,:,:,:)
  
contains
  
  subroutine init_impl_cond
    
    if(.not.allocated(MAT)) allocate(MAT(in,jn,kn,7))
    if(.not.allocated(res)) allocate(res(in,jn,kn))
    if(.not.allocated(rhs)) allocate(rhs((in-5)*(jn-5)*(kn-5)))
    if(.not.allocated(dw))  allocate(dw((in-5)*(jn-5)*(kn-5)))
    if(.not.allocated(HeatCond)) allocate(HeatCond(in,jn,kn,3))
    
  end subroutine init_impl_cond
  
  subroutine clean_impl_cond
    
    if(allocated(MAT)) deallocate(MAT)
    if(allocated(res)) deallocate(res)
    if(allocated(rhs)) deallocate(rhs)
    if(allocated(dw))  deallocate(dw)
    if(allocated(HeatCond)) deallocate(HeatCond)
    
  end subroutine clean_impl_cond
  
end module ModImplCond
