!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module CON_buffer_grid
  use CON_global_message_pass
  implicit none
  private!Except
  public::init_buffer_grid_couple,couple_buffer_grid
  real,dimension(:,:),pointer::State_VI
contains
  subroutine init_buffer_grid_couple(&
       SourceGD,&
       TargetGD,     &
       RouterToBuffer,&
       nVar,NameBuffer)
    type(GridDescriptorType),intent(in)::SourceGD,TargetGD
    type(RouterType),intent(out)::RouterToBuffer
    character(LEN=*),intent(in)::NameBuffer
    integer,intent(in)::nVar
    call init_router(&
         SourceGD,&
         TargetGD,&
         RouterToBuffer,&
         nIndexesTarget=1)
    if(is_proc(compid_grid(TargetGD%DD%Ptr)))&
         call allocate_vector(NameBuffer,nVar,TargetGD)
    if(is_proc0(compid_grid(TargetGD%DD%Ptr)))&
         write(*,*)'Allocated '//NameBuffer//' bounds=',&
         ubound_vector(NameBuffer)
  end subroutine init_buffer_grid_couple
  !=============================================================!
  subroutine couple_buffer_grid(&
!Router: should be preset. See CON_router for the instructions. !
      Router,&
      nVar,  &
!You should use one external subroutines. It                    !
!is important to mention that in the procedure itself there is  !
!absolutely no names of the data sets to be sent-received. You  !
!can send whatever you want but you should form the procedure   !
!which GET the data you want from the points defined by router  !
!(the cell and block indexes for these points are submitted by  !
!global_message_pass itself, and you should place               !
!point-by-point, whatever you  want to a short (of the nVar     !
!length) piece of the buffer which then will be placed to a     !
!particular place of the global buffer and proberly             !
!sent-received. A name of your  procedure which gets your data  !
!and place them to a small buffer should be used here. Use      !
!fill_buffer=your_procedure_name best of all...
!
      fill_buffer,&
      NameBuffer,TargetID_)
!INPUT ARGUMENTS:         
      type(RouterType),intent(in)::Router
    interface
!---------------------------------------------------------------!
!Should get nVar values from some arrays using indexes          !
!                                                               !
!Get%Index_II(:,iGetStart:iGetStart+nPartial-1), take the sum(s)!             
!with weight coefficients                                       !
!                                                               !
!Weight_I(iGetStart:iGetStart+nPartial-1), and put the nVar     !
!results to the buffer array VSend_I, starting from iBuffStart  !  
!value of indexes.                                              !
!                                                               !
!Pseudocode:
!\end{verbatim}  
!$$                                                  
!VSend_I(iBuffStart:iBuffStart+nVar-1)=
!$$
!$$                                 
!\sum_{iGet=iGetStart}^              
!{iGetStart+nPartial-1}                      
!V(1:nVar)(Get\  Index(:,iGet)),     
!$$            
!\begin{verbatim}                      
!where                                                          !
!V(1:nVar)(IndexVector)      is a vector of nVar real variables,!
!which are indexed by vector index Get%Index_II(:,iGet). As an  !
!example, the vector index can be                               !
!                                                               !
!   iCell,jCell,kCell,iBlock,                                   !
!                                                               !
!and the vector V can be the set of the conserved variables in  !
!the thus numbered control volumes                              !
       subroutine fill_buffer(nPartial,&
                              iGetStart,&
                              Get,&
                              Weight,&
                              Buff_I,nVar)
         use CON_router
         implicit none
         integer,intent(in)::nPartial,iGetStart,nVar
         type(IndexPtrType),intent(in)::Get
         type(WeightPtrType),intent(in)::Weight
         real,dimension(nVar),intent(out)::Buff_I
       end subroutine fill_buffer
    end interface
    character(LEN=*),intent(in)::NameBuffer
    integer,intent(in)::nVar,TargetID_
    if(is_proc(TargetID_))&
         call associate_with_global_vector(&
         State_VI,NameBuffer)
    call global_message_pass(Router,&
         nVar,&
         fill_buffer,&
         put_to_global_vector)
    if(is_proc(TargetID_))&
         call bcast_global_vector(NameBuffer,0,i_comm(TargetID_))
  end subroutine couple_buffer_grid
  !=============================================================!
  subroutine put_to_global_vector(nPartial,&
       iPutStart,&
       Put,&
       Weight,&
       DoAdd,&
       State_V,nVar)
    integer,intent(in)::nPartial,iPutStart,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::Weight
    logical,intent(in)::DoAdd
    real,dimension(nVar),intent(in)::State_V
    integer::iPoint
    iPoint=Put%iCB_II(1,iPutStart)
    if(DoAdd)then
       State_VI(:,iPoint)=State_VI(:,iPoint)+State_V
    else
       State_VI(:,iPoint)=State_V
    end if
  end subroutine put_to_global_vector
end module CON_buffer_grid
