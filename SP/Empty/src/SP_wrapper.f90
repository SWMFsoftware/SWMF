subroutine SP_run
  implicit none
end subroutine SP_run
subroutine write_ihdata
  implicit none
end subroutine write_ihdata
subroutine sp_smooth_ihdata
  implicit none
end subroutine sp_smooth_ihdata
subroutine sp_set_ihdata(nPointIn,XyzIn_DI)
  implicit none
  integer,intent(in)::nPointIn
  real,intent(in)::XyzIn_DI
end subroutine sp_set_ihdata
subroutine SP_put_from_ih(nPartial,iPutStart,Put,W,DoAdd,Buff_I,nVar)
  use CON_coupler, ONLY: IndexPtrType, WeightPtrType
  implicit none
  integer,intent(in)::nPartial,iPutStart,nVar
  type(IndexPtrType),intent(in)::Put
  type(WeightPtrType),intent(in)::W
  logical,intent(in)::DoAdd
  real,dimension(nVar),intent(in)::Buff_I
end subroutine SP_put_from_ih
