!^CFG COPYRIGHT UM
!========================================================================
module ModUser

  use ModUserEmpty,               &
       IMPLEMENTED1 => user_set_ics,                  &
       IMPLEMENTED2 => user_init_session,             &
       IMPLEMENTED3 => user_set_outerbcs,             &
       IMPLEMENTED4 => user_get_b0,                   &
       IMPLEMENTED5 => user_normalization,            &
       IMPLEMENTED6 => user_io_units

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
  real::WallDensity=100.0*1.0E12           !cm^{-3}
  real::WallTemperature=1.                !eV
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
    !-------------------------------------------------------------------------
    call set_oktest('user_init_session',oktest, oktest_me)

    !Dimensionless poloidal and toroidal field
    !Both are calculated at r=rTorusLarge,z=rTorusSmall
    PoloidalField=PoloidalField*Io2No_V(UnitB_)
    ToroidalField=ToroidalField*Io2No_V(UnitB_)

    !Poloidal magnetic field is B_r=-a_z*z/(cPi*r), hence:
    GradShafranovAZ=PoloidalField*cPi*rTorusLarge/rTorusSmall

    !To ensure purely round magnetic field line near the magnetic axis

    GradShafranovAR=GradShafranovAZ/(4.0*rTorusLarge**2)

    !Poloidal flux through the last magnetic surface

    FullPoloidalFlux=GradShafranovAZ*rTorusSmall**2

    call set_toroidal_surface

    !Here the arrays of the FACE VALUE are formed
    !Initialization
    do iBoundary=body2_,Top_
       FaceState_VI(:,iBoundary)=DefaultState_V(1:nVar)
    end do

    !Cell State is used for filling the ghostcells
    CellState_VI=FaceState_VI

    do iBoundary=body2_,Top_  
       CellState_VI(rhoUx_:rhoUz_,iBoundary) = &
            FaceState_VI(Ux_:Uz_,iBoundary)*FaceState_VI(rho_,iBoundary)
    end do
  contains
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
         R1=rTorusSmall*0.5
         F1=poloidal_flux(R1*CosAngle+rTorusLarge,R1*SinAngle)-FullPoloidalFlux
         R2=rTorusSmall*2.0
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

  subroutine user_normalization

    use ModPhysics, ONLY: No2Si_V, UnitX_, UnitRho_, UnitU_
    use ModConst, ONLY: cProtonMass
    real, external :: energy_in

    No2Si_V(UnitX_)    = cOne                            
    No2Si_V(UnitRho_)  = 1000000*cProtonMass*WallDensity 
    No2Si_V(UnitU_)    = &
         sqrt(energy_in('eV')*WallTemperature/cProtonMass)

  end subroutine user_normalization

  !==========================================================================

  subroutine user_io_units

    use ModPhysics
    real, external :: energy_in
    !-----------------------------------------------------------------------
    ! Set I/O units in terms of SI units if they differ

    Io2Si_V(UnitRho_)         = cProtonMass/1.0E-6             ! (amu/cm^3)
    Io2Si_V(UnitU_)           = 1000.0                         ! km/s
    Io2Si_V(UnitN_)           = 1./1.0E-6                      ! (#/cm^3)
    Io2Si_V(UnitTemperature_) = cBoltzmann/energy_in('eV')     ! eV 

    !\
    ! set strings for writing Tecplot output
    !/
    NameTecUnit_V(UnitRho_)         = '[amu/cm^3]'
    NameTecUnit_V(UnitU_)           = '[km/s]'          
    NameTecUnit_V(UnitN_)           = '[cm^-^3]'        
    NameTecUnit_V(UnitTemperature_) = '[eV]'
    !\
    ! set strings for writing IDL output
    !/
    NameTecUnit_V(UnitRho_)         = 'mp/cc'
    NameTecUnit_V(UnitU_)           = 'km/s'          
    NameTecUnit_V(UnitN_)           = '/cc'        
    NameTecUnit_V(UnitTemperature_) = 'eV'

  end subroutine user_io_units

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
    B_D(z_)=2.0*GradShafranovAR*(R**2-rTorusLarge**2)/cPi
    B_D(x_:y_)=-(/CosPhi,SinPhi/)*GradShafranovAZ*Z/(cPi*R)+&
                (/-SinPhi,CosPhi/)*&                 !ToroidalField
                (sign(cOne,ToroidalField)*sqrt(&
                (ToroidalField*rTorusLarge/R)**2+&
                (FullPoloidalFlux-PoloidalFlux)*GradShafranovAZ/&
                (cPi*R)**2)-ToroidalField*rTorusLarge/R)
    P=cOne+2.0*GradShafranovAR*(FullPoloidalFlux-PoloidalFlux)/cPi**2
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
       State_VGB(rhoUx_:rhoUz_,i,j,k,iBLK)=0.0
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
