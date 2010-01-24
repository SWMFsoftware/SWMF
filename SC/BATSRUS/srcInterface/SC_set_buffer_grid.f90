!^CFG COPYRIGHT UM
subroutine SC_set_buffer_grid(DD)
  use SC_ModBuffer,ONLY:&
       set_spher_buffer_grid,set_buffer_name,&
       DomainDecompositionType,&
       LocalBufferDD
  use CON_coupler,ONLY:SC_,LC_,is_proc
  implicit none
  type(DomainDecompositionType),&
       intent(out)::DD

  call set_spher_buffer_grid(&
       DD,SC_,IsLocal=.false.)
  if(.not.is_proc(SC_))return

  call set_spher_buffer_grid(&
       LocalBufferDD,SC_,IsLocal=.true.)
  call set_buffer_name('SC_from_lc',LC_)

end subroutine SC_set_buffer_grid
