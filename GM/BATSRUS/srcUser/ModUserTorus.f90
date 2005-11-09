!^CFG COPYRIGHT UM
!========================================================================
module ModUser

  use ModUserEmpty, ONLY:               &
       user_read_inputs,                &
       ! user_set_ics,                  &
       ! user_init_session,             &
       user_initial_perturbation,       &
       user_set_boundary_cells,         &
       user_face_bcs,                   &
       ! user_set_outerbcs,             &
       user_specify_initial_refinement, &
       user_amr_criteria,               &
       user_write_progress,             &
       user_get_log_var,                &
       user_calc_sources,               &
       user_heat_source,                &
       ! user_get_b0,                   &
       user_update_states

  use ModCovariant

  include 'user_module.h' !list of public methods

  !\
  ! Here you must define a user routine Version number and a 
  ! descriptive string.
  !/
  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = &
       'Magnetic fusion device with a planar magnetic axis, I.V. Sokolov'
  real::GradShafranovAZ,GradShafranovAR
  real::PoloidalField=0.5,ToroidalField=4.0 !T
  real::WallDensity=cHundred*cE12           !cm^{-3}
  real::WallTemperature=cOne                !eV
  real::FullPoloidalFlux

contains
  !===========================================================================
  subroutine user_init_session
    use ModProcMH
    use ModMain
    use ModPhysics
    use ModVarIndexes
    implicit none

    real :: Qqpmag, Oopmag, Gsun
    real :: Mbody_dim
    real :: MBody2Dim                 !^CFG IF SECONDBODY
    real :: cosTheta, sinTheta, cosPhi, sinPhi, & 
         xx, yy, zz

    integer :: i, iBoundary

    logical :: oktest, oktest_me
    real, external :: energy_in
    !-------------------------------------------------------------------------
    call set_oktest('user_init_session',oktest, oktest_me)

    unitSI_x       = 1.0           ! Radius - NOTE - MUST BE IN meters
    Mbody_dim      = 0.00          ! Mass of body in kg
    rot_period_dim = 0.00          ! rotation period in hours
    

    ! Second body mass is set to zero by default   !^CFG IF SECONDBODY
    MBody2Dim = 0.0                                !^CFG IF SECONDBODY

    !\
    ! Call set_dimensions, which set the quantities for converting from
    ! normalized to  dimensional quantities and vice versa.  Also
    ! sets input and output strings for the conversion quantities
    !/
    call user_set_dimensions

    if(oktest .and. iProc==0) then
       write(*,'(E15.6,11X,E15.6)') unitUSER_x, unitSI_x
       write(*,'(E15.6,11X,E15.6)') unitUSER_t, unitSI_t
       write(*,'(E15.6,11X,E15.6)') unitUSER_angle, unitSI_angle       
       write(*,'(E15.6,11X,E15.6)') unitUSER_rho,  unitSI_rho
       write(*,'(E15.6,11X,E15.6)') unitUSER_n, unitSI_n 
       write(*,'(E15.6,11X,E15.6)') unitUSER_U, unitSI_U          
       write(*,'(E15.6,11X,E15.6)') unitUSER_p, unitSI_p 
       write(*,'(E15.6,11X,E15.6)') unitUSER_B, unitSI_B                        
       write(*,'(E15.6,11X,E15.6)') unitUSER_rhoU, unitSI_rhoU
       write(*,'(E15.6,11X,E15.6)') unitUSER_energydens, unitSI_energydens
       write(*,'(E15.6,11X,E15.6)') unitUSER_J, unitSI_J
       write(*,'(E15.6,11X,E15.6)') unitUSER_electric, unitSI_electric
       write(*,'(E15.6,11X,E15.6)') unitUSER_DivB, unitSI_DivB
       write(*,'(E15.6,11X,E15.6)') unitUSER_temperature, unitSI_temperature
    end if
    !Dimensionless poloidal and toroidal field
    !Both are calculated at r=rTorusLarge,z=rTorusSmall
    PoloidalField=PoloidalField/UnitUser_B
    ToroidalField=ToroidalField/UnitUser_B
    
    !Poloidal magnetic field is B_r=-a_z*z/(cPi*r), hence:
    GradShafranovAZ=PoloidalField*cPi*rTorusLarge/rTorusSmall

    !To ensure purely round magnetic field line near the magnetic axis

    GradShafranovAR=GradShafranovAZ/(cFour*rTorusLarge**2)

    !Poloidal flux through the last magnetic surface

    FullPoloidalFlux=GradShafranovAZ*rTorusSmall**2

    call set_toroidal_surface

    !Here the arrays of the FACE VALUE are formed
    !Initialization
    do iBoundary=body2_,Top_
       FaceState_VI(:,iBoundary)=DefaultState_V
    end do

    !Cell State is used for filling the ghostcells
    CellState_VI=FaceState_VI

    do iBoundary=body2_,Top_  
       CellState_VI(rhoUx_:rhoUz_,iBoundary) = &
            FaceState_VI(Ux_:Uz_,iBoundary)*FaceState_VI(rho_,iBoundary)
    end do
  contains

    !==========================================================================

    subroutine user_set_dimensions
      use ModProcMH, ONLY:iProc
      use ModMain
      use ModPhysics
      use ModVarIndexes
      implicit none

      logical :: oktest, oktest_me
      !-----------------------------------------------------------------------
      call set_oktest('user_set_dimensions',oktest, oktest_me)

      !\
      ! Load variables used for converting from dimensional to non-dimensional 
      ! units and back.  Also load the name variables for each of the units for
      ! use in writing output.
      !/

      unitSI_angle = 180/cPi   
      unitUSER_angle = unitSI_angle
      unitstr_TEC_angle = '[degree]'            
      unitstr_IDL_angle = '--'            

      !\
      ! set independent normalizing SI variables
      !/
      unitSI_x    = cOne                                         ! m
      unitSI_rho  = cProtonMass*(cMillion)*WallDensity           ! kg/m^3
      unitSI_U    = sqrt(energy_in('eV')*WallTemperature/&
           cProtonMass)                                          ! m/s

      !\
      ! set other normalizing SI variables from the independent ones
      !/
      unitSI_t           = unitSI_x/unitSI_U                     ! s
      unitSI_n           = unitSI_rho/cProtonMass                ! (#/m^3)
      unitSI_p           = unitSI_rho*unitSI_U**2                ! Pa
      unitSI_B           = unitSI_U*sqrt(cMu*unitSI_rho)         ! T
      unitSI_rhoU        = unitSI_rho*unitSI_U                   ! kg/m^2/s
      unitSI_energydens  = unitSI_p                              ! J/m^3
      unitSI_Poynting    = unitSI_energydens*unitSI_U            ! J/m^2/s
      unitSI_J           = unitSI_B/(unitSI_x*cMu)               ! A/m^2
      unitSI_electric    = unitSI_U*unitSI_B                     ! V/m
      unitSI_DivB        = unitSI_B/unitSI_x                     ! T/m
      ! set temperature - note that the below is only strictly true for a
      ! pure proton plasma.  If the are heavy ions of mass #*mp then you could
      ! be off in temperature by as much as a factor of #.  
      ! There is no way around this in MHD.
      unitSI_temperature = (unitSI_p/unitSI_rho)*(cProtonMass/cBoltzmann) !K 

      !\
      ! set USER variables used for normalization, input, and output
      ! They all have user defined coordinates and are not consistent 
      ! in a unit sense.
      ! Every variable is case dependent.
      ! 
      ! Also load the string variables associated with the USER variables - 
      ! note that they are loaded differently for IDL and TEC output
      !/
      unitUSER_x           = unitSI_x                         ! m
      unitUSER_rho         = 1.0E-6*unitSI_n                  ! (amu/cm^3)
      unitUSER_U           = unitSI_U                         ! km/s
      unitUSER_t           = unitSI_t                         ! s
      unitUSER_n           = 1.0E-6*unitSI_n                  ! (#/cm^3)
      unitUSER_p           = unitSI_p                         ! Pa
      unitUSER_B           = unitSI_B                         ! T
      unitUSER_rhoU        = unitSI_rhoU                      ! kg/m^2/s
      unitUSER_energydens  = unitSI_energydens                ! J/m^3
      unitUSER_Poynting    = unitSI_Poynting                  ! J/m^2/s
      unitUSER_J           = unitSI_J                         ! A/m^2
      unitUSER_electric    = unitSI_electric                  ! V/m
      unitUSER_DivB        = unitSI_DivB                      ! T/m
      ! set temperature - note that the below is  only strictly true for a
      ! pure proton plasma.  If the are heavy ions of mass #*mp then you could
      ! be off in temperature by as much as a factor of #.  
      ! There is no way around this in MHD.
      unitUSER_temperature = UnitSI_temperature*cBoltzmann/&
           energy_in('eV')                                    ! eV 

      !\
      ! set string variables used for writing output - TECPLOT
      !/
      unitstr_TEC_x           = '[R]'            
      unitstr_TEC_rho         = '[amu/cm^3]'       
      unitstr_TEC_U           = '[m/s]'          
      unitstr_TEC_t           = '[s]'             
      unitstr_TEC_n           = '[amu/cm^3]'        
      unitstr_TEC_p           = '[Pa]'           
      unitstr_TEC_B           = '[T]'            
      unitstr_TEC_rhoU        = '[kg m^-^2 s^-^2]'
      unitstr_TEC_energydens  = '[J/m^3]'             
      unitstr_TEC_Poynting    = '[J m^-^2 s^-^1]'
      unitstr_TEC_J           = '[`A/m^2]'       
      unitstr_TEC_electric    = '[V/m]'          
      unitstr_TEC_DivB        = '[nT/R]'           
      unitstr_TEC_temperature = '[eV]'             
      !\
      ! set string variables used for writing output - IDL
      !/
      unitstr_IDL_x           = 'R'            
      unitstr_IDL_rho         = 'amu/cm3'       
      unitstr_IDL_U           = 'm/s'          
      unitstr_IDL_t           = 's'             
      unitstr_IDL_n           = 'amu/cm3'        
      unitstr_IDL_p           = 'Pa'           
      unitstr_IDL_B           = 'T'            
      unitstr_IDL_rhoU        = 'kg/m2s2'
      unitstr_IDL_energydens  = 'J/m3'           
      unitstr_IDL_Poynting    = 'J/m^2s'
      unitstr_IDL_J           = 'A/m2'       
      unitstr_IDL_electric    = 'V/m'          
      unitstr_IDL_DivB        = 'T/R'           
      unitstr_IDL_temperature = 'eV'             

      unitUSERVars_V(rho_)     = unitUSER_rho
      unitUSERVars_V(rhoUx_)   = unitUSER_rhoU
      unitUSERVars_V(rhoUy_)   = unitUSER_rhoU
      unitUSERVars_V(rhoUz_)   = unitUSER_rhoU
      unitUSERVars_V(Bx_)      = unitUSER_B
      unitUSERVars_V(By_)      = unitUSER_B
      unitUSERVars_V(Bz_)      = unitUSER_B
      unitUSERVars_V(P_)       = unitUSER_p

      TypeUnitVarsTec_V(rho_)    = unitstr_TEC_rho
      TypeUnitVarsTec_V(rhoUx_)  = unitstr_TEC_rhoU
      TypeUnitVarsTec_V(rhoUy_)  = unitstr_TEC_rhoU
      TypeUnitVarsTec_V(rhoUz_)  = unitstr_TEC_rhoU
      TypeUnitVarsTec_V(Bx_)     = unitstr_TEC_B
      TypeUnitVarsTec_V(By_)     = unitstr_TEC_B
      TypeUnitVarsTec_V(Bz_)     = unitstr_TEC_B
      TypeUnitVarsTec_V(P_)  = unitstr_TEC_p

      TypeUnitVarsIdl_V(rho_)    = unitstr_IDL_rho
      TypeUnitVarsIdl_V(rhoUx_)  = unitstr_IDL_rhoU
      TypeUnitVarsIdl_V(rhoUy_)  = unitstr_IDL_rhoU
      TypeUnitVarsIdl_V(rhoUz_)  = unitstr_IDL_rhoU
      TypeUnitVarsIdl_V(Bx_)     = unitstr_IDL_B
      TypeUnitVarsIdl_V(By_)     = unitstr_IDL_B
      TypeUnitVarsIdl_V(Bz_)     = unitstr_IDL_B
      TypeUnitVarsIdl_V(P_)  = unitstr_IDL_p
    end subroutine user_set_dimensions

    !==========================================================================

    subroutine set_toroidal_surface
      use ModIO,ONLY:NamePlotDir
      use ModIoUnit,ONLY:io_unit_new
      use ModProcMH
      real::PoloidalAngle,rOld,rNew,fNew,CosAngle,SinAngle
      real::R1,R2,F1,F2
      integer::iPoint,iFile
      real,parameter::dAngle=cTwoPi/nToroidalBoundaryPoints
      !-----------------------------------------------------------------------
      IsInitializedTorusGeometry=.true.
      do iPoint=0,nToroidalBoundaryPoints-1
         PoloidalAngle=iPoint*dAngle
         CosAngle=cos(PoloidalAngle)
         SinAngle=sin(PoloidalAngle)
         R1=rTorusSmall*cHalf
         F1=poloidal_flux(R1*CosAngle+rTorusLarge,R1*SinAngle)-FullPoloidalFlux
         R2=rTorusSmall*cTwo
         F2=poloidal_flux(R2*CosAngle+rTorusLarge,R2*SinAngle)-FullPoloidalFlux
         rOld=R1
         TorusSurface_I(iPoint)=R2
         if(F1*F2>cZero.or.F1==F2)call stop_mpi('Chord method failed')
         CHORD:do while (abs(TorusSurface_I(iPoint)-rOld)>cTolerance)
            rOld=TorusSurface_I(iPoint)
            TorusSurface_I(iPoint)=(F1*R2-F2*R1)/(F1-F2)
            fNew=poloidal_flux(TorusSurface_I(iPoint)*CosAngle+rTorusLarge,&
                 TorusSurface_I(iPoint)*SinAngle)-FullPoloidalFlux
            if(fNew==cZero)EXIT CHORD
            if(fNew*F1>cZero)then
               F1=fNew
               R1=TorusSurface_I(iPoint)
            else
               F2=fNew
               R2=TorusSurface_I(iPoint)
            end if
         end do CHORD
      end do
      TorusSurface_I(nToroidalBoundaryPoints)= TorusSurface_I(0)
      if(iProc/=0)return
      iFile=io_unit_new()
      open(iFile,file=trim(NamePlotDir)//'torus.dat',status='replace')
      write(iFile,*)nToroidalBoundaryPoints,rTorusSmall,rTorusLarge
      do iPoint=0,nToroidalBoundaryPoints
         write(iFile,*)iPoint,TorusSurface_I(iPoint)
      end do
      close(iFile)
    end subroutine set_toroidal_surface
  end subroutine user_init_session

  !==========================================================================

  real function poloidal_flux(R,Z)

    real,intent(in)::R,Z
    !-------------------------------------------------------------------------
    poloidal_flux=GradShafranovAZ*Z**2+&
         GradShafranovAR*(R**2-rTorusLarge**2)**2
  end function poloidal_flux

  !==========================================================================

  subroutine get_grad_shafranov(X_D,B_D,P)

    use ModMain,ONLY:x_,y_,z_
    real,dimension(nDim),intent(in)::X_D
    real,dimension(nDim),intent(out)::B_D
    real,intent(out)::P
    real::R,Z,SinPhi,CosPhi,PoloidalFlux
    !-------------------------------------------------------------------------
    R=sqrt(X_D(x_)**2+X_D(y_)**2)
    SinPhi=X_D(y_)/R;CosPhi=X_D(x_)/R
    Z=X_D(z_)
    PoloidalFlux=poloidal_flux(R,Z)
    !Poloidal field:
    B_D(z_)=cTwo*GradShafranovAR*(R**2-rTorusLarge**2)/cPi
    B_D(x_:y_)=-(/CosPhi,SinPhi/)*GradShafranovAZ*Z/(cPi*R)+&
                (/-SinPhi,CosPhi/)*&                 !ToroidalField
                (sign(cOne,ToroidalField)*sqrt(&
                (ToroidalField*rTorusLarge/R)**2+&
                (FullPoloidalFlux-PoloidalFlux)*GradShafranovAZ/&
                (cPi*R)**2)-ToroidalField*rTorusLarge/R)
    P=cOne+cTwo*GradShafranovAR*(FullPoloidalFlux-PoloidalFlux)/cPi**2
  end subroutine get_grad_shafranov

  !==========================================================================

  subroutine user_set_ICs
    use ModMain,ONLY:x_,y_,z_
    use ModMain,ONLY:globalBLK,nI,nJ,nK,gcn
    use ModIO,ONLY:restart
    use ModAdvance,ONLY:State_VGB
    use ModVarIndexes
    use ModGeometry,ONLY:x_BLK,y_BLK,z_BLK
    integer::i,j,k
    !-------------------------------------------------------------------------
    call set_b0(globalBLK)
    if(restart)return
    do k=1-gcn,nK+gcn;do j=1-gcn,nJ+gcn;do i=1-gcn,nI+gcn
       State_VGB(rho_,i,j,k,globalBLK)=cOne
       State_VGB(rhoUx_:rhoUz_,i,j,k,globalBLK)=cZero
       call get_grad_shafranov(&
            (/x_BLK(i,j,k,globalBLK),&
              y_BLK(i,j,k,globalBLK),&
              z_BLK(i,j,k,globalBLK)/),&
              State_VGB(Bx_:Bz_,i,j,k,globalBLK),&
              State_VGB(P_,i,j,k,globalBLK))
      State_VGB(P_,i,j,k,globalBLK)=max(&
            State_VGB(P_,i,j,k,globalBLK),cOne)
   end do;end do;end do
       
  end subroutine user_set_ICs

  !==========================================================================

  subroutine user_set_outerbcs(iBLK,iSide,TypeBc,IsFound)
    use ModAdvance,ONLY:State_VGB
    use ModVarIndexes
    use ModGeometry,ONLY:x_BLK,y_BLK,z_BLK
    use ModSize
    integer,intent(in)::iBLK,iSide
    character(LEN=*),intent(in)::TypeBc
    logical,intent(out)::IsFound
    integer::iStart,iFinal,jStart,jFinal,kStart,kFinal,i,j,k
    !-------------------------------------------------------------------------
    IsFound=.true.
    iStart=1-gcn;iFinal=nI+gcn
    jStart=1-gcn;jFinal=nJ+gcn
    kStart=1-gcn;kFinal=nK+gcn
    select case(iSide)
    case(East_)
       iFinal=0
    case(Bot_)
       kFinal=0
    case(West_)
       iStart=nI+1
    case(Top_)
       kStart=nK+1
    case default
       call stop_mpi('Wrong iSide in user_set_outerBCs')
    end select
    do k=kStart,kFinal;do j=jStart,jFinal;do i=iStart,iFinal
       State_VGB(rho_,i,j,k,iBLK)=cOne
       State_VGB(rhoUx_:rhoUz_,i,j,k,iBLK)=cZero
       call get_grad_shafranov(&
            (/x_BLK(i,j,k,iBLK),&
              y_BLK(i,j,k,iBLK),&
              z_BLK(i,j,k,iBLK)/),&
              State_VGB(Bx_:Bz_,i,j,k,iBLK),&
              State_VGB(P_,i,j,k,iBLK))
      State_VGB(P_,i,j,k,iBLK)=max(&
            State_VGB(P_,i,j,k,iBLK),cOne)
   end do;end do;end do

  end subroutine user_set_outerbcs

  !==========================================================================

  subroutine user_get_b0(X0,Y0,Z0,B0_D)
    use ModMain,ONLY:x_,y_,z_
    real,intent(in)::X0,Y0,Z0
    real,intent(out),dimension(nDim)::B0_D
    B0_D(z_)=cZero
    B0_D(x_:y_)=(/-Y0,X0/)/(X0**2+Y0**2)*&
         ToroidalField*rTorusLarge
  end subroutine user_get_b0

  !==========================================================================

end module ModUser
