!^CFG COPYRIGHT UM
!========================================================================
module ModUser

  use ModSize
  use ModUserEmpty,                        &
       IMPLEMENTED1 => user_set_ics,       &
       IMPLEMENTED2 => user_read_inputs,   &
       IMPLEMENTED3 => user_calc_sources,  &
       IMPLEMENTED4 => user_update_states
       

  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 0.9
  character (len=*), parameter :: NameUserModule = &
       'Rubin single species Europa MHD module, Sep. 2009'

  real, public, dimension(1:nI, 1:nJ, 1:nK, nBLK) :: Neutral_BLK
  real :: n0, H, v, alpha, mi_mean, kin
  real :: vNorm, alphaNorm, kinNorm, nNorm

contains

  !========================================================================
  !  SUBROUTINE user_read_inputs
  !========================================================================
  subroutine user_read_inputs

    use ModMain
    use ModProcMH,    ONLY: iProc
    use ModReadParam
    use ModIO,        ONLY: write_prefix, write_myname, iUnitOut

    character (len=100) :: NameCommand

   !-------------------------------------------------------------------------

    if(iProc==0.and.lVerbose > 0)then
       call write_prefix; write(iUnitOut,*)'User read_input Europa starts'
    endif
    
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
          
         case("#EUROPA")
            call read_var('n0' , n0)           !! Neutral surface density [1/cm^3]
            call read_var('H' , H)             !! Neutral scale height [km]
            call read_var('v' , v)             !! Ionization rate [1/s]
            call read_var('alpha' , alpha)     !! Recombination rate [cm^3/s]
            call read_var('kin' , kin)         !! ion neutral friction [cm^3/s]
            call read_var('mi_mean' , mi_mean) !! mean ion mass [amu]
            H=H*1E3                            !! conversion to SI  
            n0=n0*1E6                          !! conversion to SI
         case('#USERINPUTEND')
            if(iProc==0.and.lVerbose > 0)then
               call write_prefix;
             write(iUnitOut,*)'User read_input EUROPA ends'
          endif
          EXIT
      case default
         if(iProc==0) then
            call write_myname; write(*,*) &
                  'ERROR: Invalid user defined #COMMAND in user_read_inputs. '
            write(*,*) '--Check user_read_inputs for errors'
            write(*,*) '--Check to make sure a #USERINPUTEND command was used'
            write(*,*) '  *Unrecognized command was: '//NameCommand
             call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
          end if
       end select
    end do

  end subroutine user_read_inputs
 

  !========================================================================
  !  SUBROUTINE user_set_ICs
  !========================================================================
  subroutine user_set_ICs

    use ModMain, ONLY: globalBLK
    use ModPhysics
    use ModNumConst
    use ModGeometry,ONLY: R_BLK
    integer :: i,j,k

    ! neutral density in SI units
    do k=1,nK; do j=1,nJ; do i=1,nI
       Neutral_BLK(i,j,k,globalBLK) = n0*exp(-(R_BLK(i,j,k,globalBLK) - Rbody)&
            /(H/rPlanetSi))

    end do;  end do;  end do

    vNorm=1/Si2No_V(UnitT_)                           !! conversion to unitless
    alphaNorm=1E-6/Si2No_V(UnitN_)/Si2No_V(UnitT_)    !! conversion to SI to unitless
    kinNorm=1E-6/Si2No_V(UnitN_)/Si2No_V(UnitT_)      !! conversion to SI to unitless
    nNorm=Si2No_V(UnitN_)                             !! conversion to unitless

  end subroutine user_set_ICs


  !========================================================================
  !  SUBROUTINE user_calc_sources
  !========================================================================
  subroutine user_calc_sources

    use ModMain,    ONLY: GlobalBlk, nI, nJ, nK
    use ModAdvance, ONLY: State_VGB, Source_VC, &
         Rho_, RhoUx_, RhoUy_, RhoUz_, Bx_,By_,Bz_, p_, Energy_
    use ModGeometry,ONLY: x_BLK,y_BLK,z_BLK,R_BLK
    use ModPhysics
    use ModProcMH

    real, dimension(1:nI,1:nJ,1:nK) :: ux, uy, uz, uxyz, ne !!, Te
    real, dimension(1:nI,1:nJ,1:nK) :: Srho,SrhoUx,SrhoUy,SrhoUz,SBx,SBy,SBz,Sp,SE 

    integer :: i,j,k

    ux=State_VGB(rhoUx_,1:nI,1:nJ,1:nK,globalBLK) / &
          State_VGB(rho_,1:nI,1:nJ,1:nK,globalBLK)
    uy=State_VGB(rhoUy_,1:nI,1:nJ,1:nK,globalBLK) / &
          State_VGB(rho_,1:nI,1:nJ,1:nK,globalBLK)
    uz=State_VGB(rhoUz_,1:nI,1:nJ,1:nK,globalBLK) / &
          State_VGB(rho_,1:nI,1:nJ,1:nK,globalBLK)

    uxyz = ux*ux+uy*uy+uz*uz

    ! ne is the electron density in SI units
    ne = State_VGB(rho_,1:nI,1:nJ,1:nK,globalBLK)*No2SI_V(UnitN_)/mi_mean

    ! Electron temperature calculated from pressure assuming Te=Ti to calculate a more 
    ! appropriate ion-electron recombination rate. p=nkT with n=ne+ni and ne=ni (quasi-neutrality)
    !Te=State_VGB(p_,1:nI,1:nJ,1:nK,globalBLK) * NO2SI_V(UnitP_) * mi_mean * cProtonMass / &
    !     ( 2.0 * NO2SI_V(UnitRho_) * cBoltzmann * State_VGB(rho_,1:nI,1:nJ,1:nK,globalBLK) )

    ! Set the source arrays for this block to zero
    Srho   = 0.0
    SrhoUx = 0.0
    SrhoUy = 0.0
    SrhoUz = 0.0
    SBx    = 0.0
    SBy    = 0.0
    SBz    = 0.0
    SP     = 0.0
    SE     = 0.0

    do k=1,nK; do j=1,nJ; do i=1,nI
       Srho(i,j,k) = Neutral_BLK(i,j,k,globalBLK)*nNorm*mi_mean*v*vNorm &   !! newly ionized neutrals
            - alpha*alphaNorm*State_VGB(rho_,i,j,k,globalBLK)*ne(i,j,k)*Si2No_V(UnitN_) !! loss due to recombination
       
       SrhoUx(i,j,k) = - State_VGB(rho_,i,j,k,globalBLK)*( &
            Neutral_BLK(i,j,k,globalBLK)*nNorm*kin*kinNorm  &               !! loss due to charge exchange
            + alpha*alphaNorm*ne(i,j,k)*Si2No_V(UnitN_))*ux(i,j,k)          !! loss due to recombination

       SrhoUy(i,j,k) = - State_VGB(rho_,i,j,k,globalBLK)*( &
            Neutral_BLK(i,j,k,globalBLK)*nNorm*kin*kinNorm  &               !! loss due to charge exchange
            + alpha*alphaNorm*ne(i,j,k)*Si2No_V(UnitN_))*uy(i,j,k)          !! loss due to recombination

       SrhoUz(i,j,k) = - State_VGB(rho_,i,j,k,globalBLK)*( &
            Neutral_BLK(i,j,k,globalBLK)*nNorm*kin*kinNorm  &               !! loss due to charge exchange
            + alpha*alphaNorm*ne(i,j,k)*Si2No_V(UnitN_))*uz(i,j,k)          !! loss due to recombination

       SP(i,j,k) = 1/3*(v*vNorm*mi_mean + kin*kinNorm*State_VGB(rho_,i,j,k,globalBLK))* &
            Neutral_BLK(i,j,k,globalBLK)*nNorm*uxyz(i,j,k)  &               !! newly generated ions
            - State_VGB(p_,i,j,k,globalBLK)*kin *kinNorm* &
            Neutral_BLK(i,j,k,globalBLK)*nNorm &                            !! loss due to charge exchange
            - State_VGB(p_,i,j,k,globalBLK)*alpha*alphaNorm*ne(i,j,k)*Si2No_V(UnitN_) !! loss due to recombination

       SE(i,j,k) = - 0.5*State_VGB(rho_,i,j,k,globalBLK)*( &
            kin *kinNorm*Neutral_BLK(i,j,k,globalBLK)*nNorm &               !! loss due to charge exchange
            + alpha *alphaNorm*ne(i,j,k)*Si2No_V(UnitN_))*uxyz(i,j,k) &     !! loss due to recombination
            - inv_gm1*(kin*kinNorm - alpha *alphaNorm)*State_VGB(p_,i,j,k,globalBLK)
    end do;  end do;  end do

    Source_VC(rho_   ,:,:,:) = Srho   + Source_VC(rho_   ,:,:,:)
    Source_VC(rhoUx_ ,:,:,:) = SrhoUx + Source_VC(rhoUx_ ,:,:,:)
    Source_VC(rhoUy_ ,:,:,:) = SrhoUy + Source_VC(rhoUy_ ,:,:,:)
    Source_VC(rhoUz_ ,:,:,:) = SrhoUz + Source_VC(rhoUz_ ,:,:,:)
    Source_VC(Bx_    ,:,:,:) = SBx    + Source_VC(Bx_    ,:,:,:)
    Source_VC(By_    ,:,:,:) = SBy    + Source_VC(By_    ,:,:,:)
    Source_VC(Bz_    ,:,:,:) = SBz    + Source_VC(Bz_    ,:,:,:)
    Source_VC(P_     ,:,:,:) = SP     + Source_VC(P_     ,:,:,:)
    Source_VC(Energy_,:,:,:) = SE     + Source_VC(Energy_,:,:,:)

  end subroutine user_calc_sources


  !========================================================================
  !  SUBROUTINE user_update_states(iStage,iBlock)
  !========================================================================
  subroutine user_update_states(iStage,iBlock)
    use ModVarIndexes
    use ModSize
    use ModAdvance, ONLY: State_VGB
    use ModPhysics
    use ModEnergy
    integer,intent(in):: iStage,iBlock
    integer:: i,j,k

    real :: Tmin = 50.0 !! Minimum ion temperature (Europa's nightside surface temperature)

    call update_states_MHD(iStage,iBlock)

    ! Force minimum temperature:
    ! If the temperature is less than the prescribed minimum 'Tmin',
    ! set it to Tmin.

    where( State_VGB(p_,1:nI,1:nJ,1:nK,iBlock)*NO2SI_V(UnitP_) < &
           (State_VGB(rho_,1:nI,1:nJ,1:nK,iBlock)*NO2SI_V(UnitN_)/mi_mean)*cBoltzmann*Tmin )
       State_VGB(p_,1:nI,1:nJ,1:nK,iBlock) = &
           (State_VGB(rho_,1:nI,1:nJ,1:nK,iBlock)*NO2SI_V(UnitN_)/mi_mean)*cBoltzmann*Tmin/NO2SI_V(UnitP_)
    end where

    call calc_energy_cell(iBlock)

  end subroutine user_update_states

end module ModUser
