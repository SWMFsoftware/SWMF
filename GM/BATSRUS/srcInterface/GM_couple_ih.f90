!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE IH

module GM_couple_ih

  ! Couple with IH component

  implicit none
  save

  private ! except

  public:: GM_put_from_ih         ! coupling toolkit based coupler
  public:: GM_put_from_ih_buffer  ! buffer grid based coupler

  character(len=3),  public:: NameCoord
  integer,           public:: nY, nZ
  real,              public:: yMin, yMax, zMin, zMax
  real, allocatable, public:: State_VII(:,:,:)

contains
  !============================================================================
  !BOP
  !ROUTINE: GM_put_from_ih - transform and put the data got from IH_
  !INTERFACE:
  subroutine GM_put_from_ih(nPartial,&
       iPutStart,&
       Put,& 
       Weight,&
       DoAdd,&
       StateSI_V,&
       nVar)
    !USES:
    use CON_router, ONLY: IndexPtrType, WeightPtrType
    use ModAdvance, ONLY: State_VGB,rho_,rhoUx_,rhoUz_,Bx_,Bz_,P_
    use ModB0,      ONLY: B0_DGB
    use ModPhysics, ONLY: Si2No_V, UnitRho_, UnitRhoU_, UnitP_, UnitB_

    !INPUT ARGUMENTS:
    integer,intent(in)::nPartial,iPutStart,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::Weight
    logical,intent(in)::DoAdd
    real,dimension(nVar),intent(in)::StateSI_V

    !REVISION HISTORY:
    !18JUL03     I.Sokolov <igorsok@umich.edu> - intial prototype/code
    !23AUG03                                     prolog
    !03SEP03     G.Toth    <gtoth@umich.edu>   - simplified
    !19JUL04     I.Sokolov <igorsok@umich.edu> - sophisticated back 
    !                  (this is what we refer to as a development)
    !EOP

    character (len=*), parameter :: NameSub='GM_put_from_ih.f90'

    real,dimension(nVar)::State_V
    integer:: i, j, k, iBlock

    !The meaning of state intdex in buffer and in model can be 
    !different. Below are the conventions for buffer:
    integer,parameter::&
         BuffRho_  =1,&
         BuffRhoUx_=2,&
         BuffRhoUz_=4,&
         BuffBx_   =5,&
         BuffBy_   =6,&
         BuffBz_   =7,&
         BuffP_    =8

    !-----------------------------------------------------------------------

    State_V(BuffRho_)              = StateSI_V(BuffRho_) *Si2No_V(UnitRho_)
    State_V(BuffRhoUx_:BuffRhoUz_) = StateSI_V(BuffRhoUx_:BuffRhoUz_) &
         *Si2No_V(UnitRhoU_)
    State_V(BuffBx_:BuffBz_)       = StateSI_V(BuffBx_:BuffBz_)* Si2No_V(UnitB_)
    State_V(BuffP_)                = StateSI_V(BuffP_)         * Si2No_V(UnitP_)

    i      = Put%iCB_II(1,iPutStart)
    j      = Put%iCB_II(2,iPutStart)
    k      = Put%iCB_II(3,iPutStart)
    iBlock = Put%iCB_II(4,iPutStart)

    if(DoAdd)then
       State_VGB(rho_,i,j,k,iBlock) = State_VGB(rho_,i,j,k,iBlock) + &
            State_V(BuffRho_)
       State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock) = &
            State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock) + &
            State_V(BuffRhoUx_:BuffRhoUz_)
       State_VGB(Bx_:Bz_,i,j,k,iBlock) = &
            State_VGB(Bx_:Bz_,i,j,k,iBlock) + &
            State_V(BuffBx_:BuffBz_)
       State_VGB(P_,i,j,k,iBlock) = State_VGB(P_,i,j,k,iBlock) + &
            State_V(BuffP_)

    else
       State_VGB(rho_,i,j,k,iBlock)= State_V(BuffRho_)
       State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock) =  State_V(BuffRhoUx_:BuffRhoUz_)
       State_VGB(Bx_:Bz_,i,j,k,iBlock) = State_V(BuffBx_:BuffBz_) - &
            B0_DGB(:,i,j,k,iBlock)
       State_VGB(P_,i,j,k,iBlock)  = State_V(BuffP_)
    end if
  end subroutine GM_put_from_ih

  !============================================================================

  subroutine GM_put_from_ih_buffer( &
       NameCoordIn, nYIn, nZIn, yMinIn, yMaxIn, zMinIn, zMaxIn, Buffer_VII)

    use ModVarIndexes
    use ModPhysics, ONLY: Si2No_V, UnitX_,UnitRho_,UnitU_,UnitB_,UnitP_
    use ModMain, ONLY: TypeCellBc_I

    character(len=*), intent(in) :: NameCoordIn
    integer,          intent(in) :: nYIn, nZIn
    real,             intent(in) :: yMinIn, yMaxIn, zMinIn, zMaxIn
    real,             intent(in) :: Buffer_VII(nVar, nYIn, nZIn)

    integer                      :: j, k
    character(len=*), parameter  :: NameSub = 'GM_put_from_ih_buffer.f90'
    !--------------------------------------------------------------------------
    if(.not.allocated(State_VII)) then
       ! Check coordinate system. Only GSM and GSE make sense.
       if(NameCoordIn /= 'GSM' .and. NameCoord /= 'GSE') &
            call CON_stop(NameSub//': cannot handle coord system=' &
            //NameCoordIn)

       ! Store grid information
       NameCoord = NameCoordIn
       yMin = yMinIn * Si2No_V(UnitX_)
       yMax = yMaxIn * Si2No_V(UnitX_)
       zMin = zMinIn * Si2No_V(UnitX_)
       zMax = zMaxIn * Si2No_V(UnitX_)
       nY   = nYIn
       nZ   = nZIn
       ! Allocate space
       allocate(State_VII(nVar,nY,nZ))

       ! Make sure that GM uses the IH buffer
       TypeCellBc_I(2) = 'ihbuffer'

       ! Debugging
       !write(*,*)'!!! NameCoord, nY, nZ=',NameCoord,nY,nZ
       !write(*,*)'!!! yMin, yMax, zMin, zMax=',yMin, yMax, zMin, zMax
    end if

    ! Store input data
    State_VII = Buffer_VII

    ! Convert units and velocity to momentum
    do k=1,nZ; do j=1,nY
       State_VII(Rho_,j,k)          = State_VII(Rho_,j,k)    *Si2No_V(UnitRho_)
       State_VII(RhoUx_:RhoUz_,j,k) = &
            State_VII(Rho_,j,k)*State_VII(Rhoux_:RhoUz_,j,k) *Si2No_V(UnitU_)
       State_VII(Bx_:Bz_,j,k)       = State_VII(Bx_:Bz_,j,k) *Si2No_V(UnitB_)
       State_VII(P_,j,k)            = State_VII(P_,j,k)      *Si2No_V(UnitP_)
    end do; end do

  end subroutine GM_put_from_ih_buffer

end module GM_couple_ih
