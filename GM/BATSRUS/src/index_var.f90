!^CFG COPYRIGHT UM
!If the Name coincides with one of the elements of NameVar_V array,
!then index_var is the number of this element of the array, 
!otherwise NameVar_V is 0 
integer function index_var(Name)
  use ModVarIndexes,ONLY:nVar,NameVar_V
  implicit none
  character(LEN=*),intent(in)::Name
  character(LEN=len_trim(Name))::NameTrim
  index_var=nVar;NameTrim=trim(Name)
  do while(index_var>0)
     if(trim(NameVar_V(index_var))==NameTrim)return
     index_var=index_var-1
  end do
end function index_var
