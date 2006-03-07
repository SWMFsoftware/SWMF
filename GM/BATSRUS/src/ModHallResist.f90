!^CFG FILE HALLRESFLUX
!^CFG COPYRIGHT UM
module ModHallResist

  use ModPhysics,  ONLY: Eta0Resist 
  implicit none

  real, allocatable :: HallJ_CD(:,:,:,:)
  real, allocatable :: BxPerRho_G(:,:,:),ByPerRho_G(:,:,:),BzPerRho_G(:,:,:)
  real              :: HallFactor
  real              :: Eta0Resist_ND,  IonMassPerCharge 

  ! Logical for adding hall resistivity
  logical:: UseHallResistFlux=.false.

  ! Logical for taking whistler wave speed into account
  logical:: UseHallCmax=.true.

contains
  !============================================================================
  subroutine init_hall_resist
    use ModSize,    ONLY: nI, nJ, nK, nDim
    use ModConst,   ONLY: cProtonMass, cElectronCharge
    use ModPhysics, ONLY: UnitSI_B, UnitSI_T, UnitSI_X

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='init_hall_resist'
    !--------------------------------------------------------------------------

    call set_oktest(NameSub, DoTest, DoTestMe)

    allocate(&
         HallJ_CD(nI,nJ,nK,nDim), &
         BxPerRho_G(0:nI+1,0:nJ+1,0:nK+1),&
         ByPerRho_G(0:nI+1,0:nJ+1,0:nK+1),&
         BzPerRho_G(0:nI+1,0:nJ+1,0:nK+1) )
    IonMassPerCharge =HallFactor*(cProtonMass/cElectronCharge) &
         * (UnitSI_B*UnitSI_T/UnitSI_X**2)
    Eta0Resist_ND =Eta0Resist*(UnitSI_T/UnitSI_X**2) 
    if (DoTestMe) then
       write(*,*) ''
       write(*,*) '>>>>>>>>>>>>>>>>> HALL Resistivity Parameters <<<<<<<<<<'
       write(*,*)
       write(*,*) 'Eta0             = ',Eta0Resist
       write(*,*) 'Eta0_ND          = ',Eta0Resist_ND
       write(*,*) 'Hallfactor       = ',Hallfactor
       write(*,*) 'IonMassPerCharge = ',IonMassPerCharge
       write(*,*) 'Omega_Bi=B0/IonMassPerCharge'
       write(*,*)
       write(*,*) '>>>>>>>>>>>>>>>>>                       <<<<<<<<<<<<<<<<<'
       write(*,*) ''
    end if

  end subroutine init_hall_resist

  !============================================================================

  subroutine add_hall_resist_flux(DoResChangeOnly)
    use ModProcMH,   ONLY:iProc
    use ModSize,     ONLY:nI,nJ,nK,gcn,nBLK
    use ModMain,     ONLY:nIFace,nJFace,nKFace,&
         iMinFaceY,iMaxFaceY,iMinFaceZ,iMaxFaceZ, &
         jMinFaceX,jMaxFaceX,jMinFaceZ,jMaxFaceZ, &
         kMinFaceX,kMaxFaceX,kMinFaceY,kMaxFaceY, &
         globalBLK, &
         x_, y_, z_, &
         iTest,jTest,kTest,VarTest,BlkTest,ProcTest
    use ModVarIndexes,ONLY:Bx_,By_,Bz_,&
         rho_,Energy_
    use ModGeometry, ONLY:x_BLK,y_BLK,z_BLK,dx_BLK, &
         dy_BLK,dz_BLK,fAx_BLK,fAy_BLK,fAz_BLK
    use ModParallel, ONLY:neiLeast,neiLwest,neiLsouth, &
         neiLnorth,neiLtop,neiLbot
    use ModAdvance,  ONLY:State_VGB, &
         B0xCell_BLK,B0yCell_BLK,B0zCell_BLK, & 
         VdtFace_x,VdtFace_y,VdtFace_z, &
         Flux_VX,Flux_VY,Flux_VZ
    use ModNumConst, ONLY:cOne,cTwo,cFour,cHalf, &
         cZero,cTiny,cHundred,cHundredth,cPi
    use ModPhysics,  ONLY: gm1
    use ModMpi

    logical, intent(in):: DoResChangeOnly
    integer:: i,j,k
    real:: AreaHalf
    real:: Jx, Jy, Jz, HallCoeff, EtaJx, EtaJy, EtaJz, DiffVDt
    real,dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn):: &
         Bx_G,By_G,Bz_G,EtaJx_G,EtaJy_G,EtaJz_G

    character(len=*), parameter :: NameSub = 'add_hall_resist_flux'
    logical :: DoTest, DoTestMe
    !---------------------------------------------------------------------------
    if(globalBLK == BlkTest .and. iProc == ProcTest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest = .false.; DoTestMe = .false.
    end if

    if(.not.allocated(HallJ_CD)) call init_hall_resist

    if(DoTestMe)then
       write(*,*)NameSub,' starting'
       write(*,*)NameSub,':Left  Flux_VX=',Flux_VX(VarTest,iTest,jTest,kTest)
       write(*,*)NameSub,':Right Flux_VX=',Flux_VX(VarTest,iTest+1,jTest,kTest)
    end if

    Bx_G = State_VGB(Bx_,:,:,:,globalBLK)+B0xCell_BLK(:,:,:,globalBLK)
    By_G = State_VGB(By_,:,:,:,globalBLK)+B0yCell_BLK(:,:,:,globalBLK)
    Bz_G = State_VGB(Bz_,:,:,:,globalBLK)+B0zCell_BLK(:,:,:,globalBLK)

    !\
    ! Compute and add the x_resistive_flux to the x-face fluxes 
    !/
    AreaHalf = cHalf*fAx_BLK(globalBLK)
    DiffVDt  = ctwo*Eta0Resist_ND*fAx_BLK(globalBLK)/dx_BLK(globalBLK) 
    if (.not.DoResChangeOnly) then
       do k=kMinFaceX,kMaxFaceX
          do j=jMinFaceX,jMaxFaceX
             do i=1,nIFace
                call add_hallresistive_flux_x
             end do
          end do
       end do
    else if (neiLeast(globalBLK)==+1) then
       i=1
       do k=1,nK
          do j=1,nJ
             call add_hallresistive_flux_x
          end do
       end do
    else if ( neiLwest(globalBLK)==+1) then
       i=nIFace
       do k=1,nK
          do j=1,nJ
             call add_hallresistive_flux_x
          end do
       end do
    end if
    !\  
    ! Compute and add the y_resistive_flux to the y-face fluxes 
    !/
    AreaHalf = cHalf*fAy_BLK(globalBLK)
    DiffVDt  = ctwo*Eta0Resist_ND*fAy_BLK(globalBLK)/dy_BLK(globalBLK) 
    if (.not.DoResChangeOnly) then
       do k=kMinFaceY,kMaxFaceY
          do j=1,nJFace
             do i=iMinFaceY,iMaxFaceY
                call add_hallresistive_flux_y
             end do
          end do
       end do
    else if(neiLsouth(globalBLK)==+1)then
       j=1
       do k=1,nK
          do i=1,nI
             call add_hallresistive_flux_y
          end do
       end do
    else if (neiLnorth(globalBLK)==+1) then
       j=nJFace 
       do k=1,nK
          do i=1,nI
             call add_hallresistive_flux_y
          end do
       end do
    end if
    !\
    ! Compute and add the z_resistive_flux to the z-face fluxes  
    !/
    AreaHalf = cHalf*fAz_BLK(globalBLK)
    DiffVDt  = ctwo*Eta0Resist_ND*fAz_BLK(globalBLK)/dz_BLK(globalBLK) 
    if (.not.DoResChangeOnly) then
       do k=1,nKFace
          do j=jMinFaceZ,jMaxFaceZ
             do i=iMinFaceZ,iMaxFaceZ
                call add_hallresistive_flux_z
             end do
          end do
       end do
    else if (neiLbot(globalBLK)==+1) then
       k=1
       do j=1,nJ
          do i=1,nI
             call add_hallresistive_flux_z
          end do
       end do
    else if ( neiLtop(globalBLK)==+1) then
       k=nKFace            
       do j=1,nJ
          do i=1,nI
             call add_hallresistive_flux_z
          end do
       end do
    end if

    if(DoTestMe)then
       write(*,*)NameSub,' finished'
       write(*,*)NameSub,': Left  Flux_VX=',Flux_VX(VarTest,iTest,jTest,kTest)
       write(*,*)NameSub,': Right Flux_VX=',Flux_VX(VarTest,iTest+1,jTest,kTest)
    end if

  contains
    !==========================================================================
    subroutine add_hallresistive_flux_x

      real :: EtaJXBx

      ! curl E = curl (Eta.J) = curl(Eta0*J + HallFactor*J x B/(n*e))
      HallCoeff = IonMassPerCharge/&
           (State_VGB(Rho_,i,j,k,globalBLK)+State_VGB(Rho_,i-1,j,k,globalBLK))

      call get_face_current(x_,i,j,k,globalBLK, Jx, Jy, Jz)

      EtaJy = fAx_BLK(globalBLK)*( Eta0Resist_ND*Jy + HallCoeff* &
           ((Bx_G(i,j,k)+Bx_G(i-1,j,k))*Jz  &
           -(Bz_G(i,j,k)+Bz_G(i-1,j,k))*Jx))

      EtaJz = fAx_BLK(globalBLK)*( Eta0Resist_ND*Jz + HallCoeff* &
           ((By_G(i,j,k)+By_G(i-1,j,k))*Jx  &
           -(Bx_G(i,j,k)+Bx_G(i-1,j,k))*Jy ))

      ! Poynting flux: (Eta.J x B)_x
      EtaJXBx = cHalf * &
           (EtaJy * (Bz_G(i,j,k) + Bz_G(i-1,j,k))   &  
           -EtaJz * (By_G(i,j,k) + By_G(i-1,j,k)))
      !\
      ! Update the `resistive flux' in the induction equation
      !/
      ! Flux_x(By) = -(eta.J)_z
      Flux_VX(By_,i,j,k) = Flux_VX(By_,i,j,k) - EtaJz

      ! Flux_x(Bz) =  (eta.J)_y
      Flux_VX(Bz_,i,j,k) = Flux_VX(Bz_,i,j,k) + EtaJy

      ! Update the `resistive flux' in the energy equation
      Flux_VX(Energy_,i,j,k) = Flux_VX(Energy_,i,j,k) + EtaJXBx

      ! Update the CFL condition
      VdtFace_x(i,j,k) = VdtFace_x(i,j,k) + DiffVDt

    end subroutine add_hallresistive_flux_x

    !==========================================================================

    subroutine add_hallresistive_flux_y
      implicit none
      real :: EtaJXBy

      logical :: DoTestCell

      DoTestCell = DoTestMe .and. i==iTest.and.j==jTest .and. k==kTest

      ! curl E = curl (Eta.J) = curl(Eta0*J + HallFactor*J x B/(n*e))
      HallCoeff = IonMassPerCharge/&
           (State_VGB(Rho_,i,j,k,globalBLK)+State_VGB(Rho_,i,j-1,k,globalBLK))

      call get_face_current(y_,i,j,k,globalBLK, Jx, Jy, Jz)

      EtaJx = fAy_BLK(globalBLK)*( Eta0Resist_ND*Jx + HallCoeff* &
           ((Bz_G(i,j,k)+Bz_G(i,j-1,k))*Jy  &
           -(By_G(i,j,k)+By_G(i,j-1,k))*Jz))

      EtaJz = fAy_BLK(globalBLK)*( Eta0Resist_ND*Jz + HallCoeff* &
           ((By_G(i,j,k)+By_G(i,j-1,k))*Jx  &
           -(Bx_G(i,j,k)+Bx_G(i,j-1,k))*Jy ))

      ! Poynting flux: (Eta.J x B)_y
      EtaJXBy = cHalf * &
           (EtaJz * (Bx_G(i,j,k) + Bx_G(i,j-1,k))   &  
           -EtaJx * (Bz_G(i,j,k) + Bz_G(i,j-1,k)))
      !\
      ! Update the `resistive flux' in the induction equation
      !/
      ! Flux_x(Bx) =  (eta.J)_z
      Flux_VY(Bx_,i,j,k) = Flux_VY(Bx_,i,j,k) + EtaJz

      ! Flux_x(Bz) = -(eta.J)_x
      Flux_VY(Bz_,i,j,k) = Flux_VY(Bz_,i,j,k) - EtaJx

      ! Update the `resistive flux' in the energy equation
      Flux_VY(Energy_,i,j,k) = Flux_VY(Energy_,i,j,k) + EtaJXBy

      ! Update the CFL condition
      VdtFace_y(i,j,k) = VdtFace_y(i,j,k) + DiffVDt

    end subroutine add_hallresistive_flux_y

    !==========================================================================

    subroutine add_hallresistive_flux_z
      implicit none
      real :: EtaJXBz

      ! curl E = curl (Eta.J) = curl(Eta0*J + HallFactor*J x B/(n*e))
      HallCoeff = IonMassPerCharge/&
           (State_VGB(Rho_,i,j,k,globalBLK)+State_VGB(Rho_,i,j,k-1,globalBLK))

      call get_face_current(z_,i,j,k,globalBLK, Jx, Jy, Jz)

      EtaJx = fAz_BLK(globalBLK)*( Eta0Resist_ND*Jx + HallCoeff* &
           ((Bz_G(i,j,k)+Bz_G(i,j,k-1))*Jy  &
           -(By_G(i,j,k)+By_G(i,j,k-1))*Jz))

      EtaJy = fAz_BLK(globalBLK)*( Eta0Resist_ND*Jy + HallCoeff* &
           ((Bx_G(i,j,k)+Bx_G(i,j,k-1))*Jz  &
           -(Bz_G(i,j,k)+Bz_G(i,j,k-1))*Jx ))

      ! Poynting flux: (Eta.J x B)_z
      EtaJXBz = cHalf * &
           (EtaJx * (By_G(i,j,k) + By_G(i,j,k-1))   &  
           -EtaJy * (Bx_G(i,j,k) + Bx_G(i,j,k-1)))
      !\
      ! Update the `resistive flux' in the induction equation
      !/
      ! Flux_z(Bx) =  -(eta.J)_y
      Flux_VZ(Bx_,i,j,k) = Flux_VZ(Bx_,i,j,k) - EtaJy

      ! Flux_x(By) =  (eta.J)_x
      Flux_VZ(By_,i,j,k) = Flux_VZ(By_,i,j,k) + EtaJx

      ! Update the `resistive flux' in the energy equation
      Flux_VZ(Energy_,i,j,k) = Flux_VZ(Energy_,i,j,k) + EtaJXBz

      ! Update the CFL condition
      VdtFace_z(i,j,k) = VdtFace_z(i,j,k) + DiffVDt

    end subroutine add_hallresistive_flux_z

  end subroutine add_hall_resist_flux

  !============================================================================

  subroutine get_face_current(iDir, i, j, k, iBlock, Jx, Jy, Jz)

    use ModAdvance, ONLY: State_VGB, Bx_, By_, Bz_
    use ModMain,    ONLY: x_, y_, z_
    use ModGeometry,ONLY: Dx_BLK, Dy_BLK, Dz_BLK

    implicit none
    !  logical, intent(in):: DoTest
    integer, intent(in):: iDir, i, j, k, iBlock
    real, intent(out)  :: Jx, Jy, Jz

    real :: InvDx, InvDy, InvDz, InvDx4, InvDy4, InvDz4


    InvDx = 1.0/dx_Blk(iBlock)
    InvDy = 1.0/dy_Blk(iBlock)
    InvDz = 1.0/dz_Blk(iBlock)
    InvDx4 = 0.25*InvDx
    InvDy4 = 0.25*InvDy
    InvDz4 = 0.25*InvDz

    select case(iDir)
    case(x_)
       Jx = +InvDy4*(State_VGB(Bz_,i,j+1,k  ,iBlock)+State_VGB(Bz_,i-1,j+1,k  ,iBlock)  &
            -        State_VGB(Bz_,i,j-1,k  ,iBlock)-State_VGB(Bz_,i-1,j-1,k  ,iBlock)) &
            -InvDz4*(State_VGB(By_,i,j  ,k+1,iBlock)+State_VGB(By_,i-1,j  ,k+1,iBlock)  &
            -        State_VGB(By_,i,j  ,k-1,iBlock)-State_VGB(By_,i-1,j  ,k-1,iBlock))

       Jy = +InvDz4*(State_VGB(Bx_,i,j,k+1,iBlock)+State_VGB(Bx_,i-1,j,k+1,iBlock)  &
            -        State_VGB(Bx_,i,j,k-1,iBlock)-State_VGB(Bx_,i-1,j,k-1,iBlock)) &
            -InvDx* (State_VGB(Bz_,i,j,k  ,iBlock)-State_VGB(Bz_,i-1,j,k  ,iBlock))

       Jz = +InvDx* (State_VGB(By_,i,j  ,k,iBlock)-State_VGB(By_,i-1,j  ,k,iBlock)) &
            -InvDy4*(State_VGB(Bx_,i,j+1,k,iBlock)+State_VGB(Bx_,i-1,j+1,k,iBlock)  &
            -        State_VGB(Bx_,i,j-1,k,iBlock)-State_VGB(Bx_,i-1,j-1,k,iBlock))

    case(y_)
       Jx = +InvDy* (State_VGB(Bz_,i,j,k  ,iBlock)-State_VGB(Bz_,i,j-1,k  ,iBlock)) &
            -InvDz4*(State_VGB(By_,i,j,k+1,iBlock)+State_VGB(By_,i,j-1,k+1,iBlock)  &
            -        State_VGB(By_,i,j,k-1,iBlock)-State_VGB(By_,i,j-1,k-1,iBlock))

       Jy = +InvDz4*(State_VGB(Bx_,i  ,j,k+1,iBlock)+State_VGB(Bx_,i  ,j-1,k+1,iBlock)  &
            -        State_VGB(Bx_,i  ,j,k-1,iBlock)-State_VGB(Bx_,i  ,j-1,k-1,iBlock)) &
            -InvDx4*(State_VGB(Bz_,i+1,j,k  ,iBlock)+State_VGB(Bz_,i+1,j-1,k  ,iBlock)  &
            -        State_VGB(Bz_,i-1,j,k  ,iBlock)-State_VGB(Bz_,i-1,j-1,k  ,iBlock))

       Jz = +InvDx4*(State_VGB(By_,i+1,j,k,iBlock)+State_VGB(By_,i+1,j-1,k,iBlock)  &
            -        State_VGB(By_,i-1,j,k,iBlock)-State_VGB(By_,i-1,j-1,k,iBlock)) &
            -InvDy* (State_VGB(Bx_,i  ,j,k,iBlock)-State_VGB(Bx_,i  ,j-1,k,iBlock)) 

    case(z_)
       Jx = +InvDy4*(State_VGB(Bz_,i,j+1,k,iBlock)+State_VGB(Bz_,i,j+1,k-1,iBlock)  &
            -        State_VGB(Bz_,i,j-1,k,iBlock)-State_VGB(Bz_,i,j-1,k-1,iBlock)) &
            -InvDz* (State_VGB(By_,i,j  ,k,iBlock)-State_VGB(By_,i,j  ,k-1,iBlock))  

       Jy = +InvDz* (State_VGB(Bx_,i  ,j,k,iBlock)-State_VGB(Bx_,i  ,j,k-1,iBlock)) &
            -InvDx4*(State_VGB(Bz_,i+1,j,k,iBlock)+State_VGB(Bz_,i+1,j,k-1,iBlock)  &
            -        State_VGB(Bz_,i-1,j,k,iBlock)-State_VGB(Bz_,i-1,j,k-1,iBlock))

       Jz = +InvDx4*(State_VGB(By_,i+1,j  ,k,iBlock)+State_VGB(By_,i+1,j  ,k-1,iBlock)  &
            -        State_VGB(By_,i-1,j  ,k,iBlock)-State_VGB(By_,i-1,j  ,k-1,iBlock)) &
            -InvDy4*(State_VGB(Bx_,i  ,j+1,k,iBlock)+State_VGB(Bx_,i  ,j+1,k-1,iBlock)  &
            -        State_VGB(Bx_,i  ,j-1,k,iBlock)-State_VGB(Bx_,i  ,j-1,k-1,iBlock))

       !if(DoTest)then
       !   write(*,*)'i,j,k,iBlock=',i,j,k,iBlock
       !   write(*,*)'InvDx4,InvDy4,InvDz4=',InvDx4,InvDy4,InvDz4
       !   write(*,*)'State_VGB(Bz)=',State_VGB(Bz_,i,j+1,k,iBlock),
       !end if
    case default
       write(*,*)'Error in get_face_current: iDir=',iDir
       call stop_mpi('DEBUG')
    end select

  end subroutine get_face_current

end module ModHallResist
