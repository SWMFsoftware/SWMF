!^CFG COPYRIGHT UM

module ModUser

  use ModSize,     ONLY: nI,nJ,nK,gcn,nBLK
  use ModPhysics,  ONLY: ElectronCharge

  use ModUserEmpty, &
       IMPLEMENTED1 => user_calc_sources, &
       IMPLEMENTED2 => user_read_inputs

  include 'user_module.h' !list of public methods

  real, parameter :: VersionUserModule = 1.3
  character (len=*), parameter :: NameUserModule = &
       'Non-Gyrotropic Reconnection Model, Kuznetsova, Toth'

  real, parameter :: UnsetX_ = -100.0

  ! Reconnection region where source terms are applied
  real :: xRecoMin  = -99.0
  real :: xRecoMax  =   0.0

  ! Search region for reconnection line
  real :: xLineMin  = -87.0
  real :: xLineMax  =  -3.0 ! for first line
  real :: xLineMax1 =  -5.0 ! for second line

  real :: yLineMax = 23.38  ! maximum distance in Y direction
  real :: DyLine   = 0.125  ! Y resolution of reconnection line search
  integer :: nJJ   = 187    ! = floor(yLineMax/DyLine)

  integer, parameter :: MaxJJ=1000 ! should exceed nJJ

  ! By default assume hydrogene for the effective ion mass
  real :: IonMassReconnect = 1.0, IonMassPerCharge

  real, dimension(MaxJJ) :: YYR, &
       XX0, DWZ_0, DWX_0, EF0_0, BZPR_0, VXPR_0, Cos_Thet_0, Sin_Thet_0, &
       XX1, DWZ_1, DWX_1, EF0_1, BZPR_1, VXPR_1, Cos_Thet_1, Sin_Thet_1

  integer, parameter :: nRecParam=7, dBzDx_=1, dVxDx_=2, WidthX_=3, &
       Efield_=4, WidthZ_=5, CosTheta_=6, SinTheta_=7

  real :: RecParam_I(nRecParam)

contains

  subroutine user_read_inputs
    use ModProcMH,    ONLY: iProc
    use ModReadParam
    implicit none

    character (len=100) :: NameCommand
    !-------------------------------------------------------------------------
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case('#LIMITX')
          call read_var('xRecoMin',xRecoMin)
          call read_var('xRecoMax',xRecoMax)
          call read_var('xLineMin',xLineMin)
          call read_var('xLineMax',xLineMax)
          call read_var('xLineMax1',xLineMax1)

       case('#LIMITY')
          call read_var('yLineMax',yLineMax)
          call read_var('DyLine',DyLine)
          nJJ = nint(yLineMax / DyLine)

       case('#IONMASS')
          call read_var('IonMass', IonMassReconnect)

       case('#USERINPUTEND')
          if(iProc==0) write(*,*)'USERINPUTEND'
          EXIT

       case default
          if(iProc==0) call stop_mpi( &
               'user_read_inputs: unrecognized command: '//NameCommand)
       end select
    end do
  end subroutine user_read_inputs

  !============================================================================

  subroutine user_calc_sources
    use ModMain
    use ModVarIndexes
    use ModAdvance
    use ModGeometry
    use ModConst
    use ModNumConst
    use ModPhysics
    use ModProcMH
    use ModPhysics
    use ModParallel

    ! User declared local variables go here
    integer :: i,j,k, jj
    real :: xt, yt, zt
    real :: xxx, zzz, yyy
    real :: InvWidthX, InvWidthZ, Ey, dBxDt, dBzDt

    integer :: nStepOld = -1

    character(len=*), parameter :: NameSub='user_calc_sources'
    logical :: DoTest, DoTestMe
    !-------------------------------------------------------------------------

    IonMassPerCharge = IonMassReconnect / ElectronCharge

    ! find reconnection lines at the beginning of the time step
    if(n_Step /= nStepOld)then
       nStepOld = n_Step
       call set_rec_line
       call set_rec_line1

       call set_oktest(NameSub, DoTest, DoTestMe)
       if (DoTestMe) then
          do jj = 1, nJJ/2+1, nJJ/4
             write(*,*)'jj=',jj,' Y=',YYR(jj)
             if(XX0(jj) <= UnsetX_) CYCLE
             write(*,*)'XX0=',XX0(jj), &
                  'BZPR=',BZPR_0(jj), 'VXPR=',VXPR_0(jj), &
                  'DWX=', DWX_0(jj), 'DWZ=',DWZ_0(jj), &
                  'EF0=', EF0_0(jj)
             !    'CosThet=',Cos_Thet_0(1), 'SinThet=', Sin_Thet_0(jj)

             if(XX1(jj) <= UnsetX_) CYCLE
             write(*,*)'XX1=',XX1(jj), &
                  'BZPR=',BZPR_1(jj), 'VXPR=',VXPR_1(jj), &
                  'DWX=', DWX_1(jj), 'DWZ=',DWZ_1(jj), &
                  'EF0=', EF0_1(jj)
             !    'CosThet=',Cos_Thet_1(1), 'SinThet=', Sin_Thet_1(jj)
          end do
       endif
    end if

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       xt = x_BLK(i,j,k,globalBLK)
       yt = abs(y_BLK(i,j,k,globalBLK))
       zt = z_BLK(i,j,k,globalBLK)

       if( yt > yLineMax) CYCLE

       jj = 1 + (yt - DyLine/2)/DyLine

       if(jj < 1 .or. jj > nJJ)then
          write(*,*)'ERROR: yt, jj, nJJ=',yt,jj,nJJ
          call stop_mpi('jj index out of 1..nJJ range')
       end if

       if(xt > xRecoMin .and. xt < xRecoMax .and. XX0(jj) > UnsetX_)then

          InvWidthX = 2/DWX_0(jj)
          InvWidthZ = 2/DWZ_0(jj)
          ! InvWidthZ = 2/(IonMassPerCharge*2.*sqrt(2.)) !!! Masha's experiment

          xxx = InvWidthX*(xt - XX0(jj))
          yyy = yt/yLineMax
          zzz = InvWidthZ*zt
          Ey  = EF0_0(jj)/(cosh(zzz)*cosh(xxx)*cosh(yyy))

          ! B' = B - curl Ey: Bx' = Bx + dEy/dZ  and Bz' = Bz - dEy/dX
          ! d(1/cosh(xxx))/dx = -sinh(xxx)/cosh(xxx)^2*InvWidthX

          dBxDt =  - Ey*tanh(zzz)*InvWidthZ
          dBzDt =  + Ey*tanh(xxx)*InvWidthX
          
          Source_VC(Bx_,i,j,k)     = Source_VC(Bx_,i,j,k) + dBxDt
          Source_VC(Bz_,i,j,k)     = Source_VC(Bz_,i,j,k) + dBzDt
          Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) &
               + State_VGB(Bx_,i,j,k,globalBlk)*dBxDt &
               + State_VGB(Bz_,i,j,k,globalBlk)*dBzDt

          !Source_VC(Bx_,i,j,k) = Source_VC(Bx_,i,j,k) &
          !     - Ey*tanh(zzz)*InvWidthZ*Cos_Thet_0(jj)
          !Source_VC(Bz_,i,j,k) = Source_VC(Bz_,i,j,k) &
          !     + Ey*tanh(xxx)*InvWidthX*Sin_Thet_0(jj)
       endif

       ! Add term from second reconnection site
       ! Make sure that the second line is separated from the first one
       if(xt > xRecoMin .and. xt < XX0(jj) .and. XX1(jj) > UnsetX_ &
            .and. XX1(jj) < XX0(jj)-DWX_0(jj)*0.5) then

          InvWidthX = 2/DWX_1(jj)
          InvWidthZ = 2/DWZ_1(jj)

          xxx = InvWidthX*(xt - XX1(jj))
          yyy = yt/yLineMax
          zzz = InvWidthZ*zt
          Ey  = EF0_1(jj)/(cosh(zzz)*cosh(xxx)*cosh(yyy))

          ! B' = B - curl Ey: Bx' = Bx + dEy/dZ  and Bz' = Bz - dEy/dX
          ! d(1/cosh(xxx))/dx = -sinh(xxx)/cosh(xxx)^2*InvWidthX

          dBxDt =  - Ey*tanh(zzz)*InvWidthZ
          dBzDt =  + Ey*tanh(xxx)*InvWidthX
          
          Source_VC(Bx_,i,j,k)     = Source_VC(Bx_,i,j,k) + dBxDt
          Source_VC(Bz_,i,j,k)     = Source_VC(Bz_,i,j,k) + dBzDt
          Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) &
               + State_VGB(Bx_,i,j,k,globalBlk)*dBxDt &
               + State_VGB(Bz_,i,j,k,globalBlk)*dBzDt

          !Source_VC(Bx_,i,j,k) = Source_VC(Bx_,i,j,k) &
          !     - Ey*tanh(zzz)*InvWidthZ*Cos_Thet_1(jj)
          !Source_VC(Bz_,i,j,k) = Source_VC(Bz_,i,j,k) &
          !     + Ey*tanh(xxx)*InvWidthX*Sin_Thet_1(jj)
       endif

    end do; end do; end do     ! end the k,j,i loops

  end subroutine user_calc_sources

  !========================================================================
  subroutine set_rec_line
    use ModMain
    use ModVarIndexes
    use ModSize
    use ModIO
    use ModPhysics
    use ModGeometry
    use ModAdvance
    use ModParallel
    use ModProcMH
    use ModConst
    use ModMpi

    implicit none

    integer :: iPE, iPEmax, iBLK, iBLKmax

    integer :: i,j,k
    integer :: ii,jj
    integer :: j1, iX, jY
    integer :: iError
    real :: XX0max, XX0max_all, yj1
    !--------------------------------------------------------------------------
    do jj=1, njj
       !   YYR(jj)=0.03125+0.0625*(jj -1)
       YYR(jj) = DyLine/2 + (jj-1)*DyLine

       XX0(jj) = UnsetX_
       XX0max  = UnsetX_

       BZPR_0(jj) = -1.
       VXPR_0(jj) = 0.
       DWX_0(jj) = -1.
       DWZ_0(jj) = -1.
       EF0_0(jj) = 0.
       Cos_Thet_0(jj) = 1.
       Sin_Thet_0(jj) = 0.
       RecParam_I(1) = -1.
       RecParam_I(2) = 0.
       RecParam_I(3) = -1.
       RecParam_I(4) = 0.
       RecParam_I(5) = -1.
       RecParam_I(6) = 1.
       RecParam_I(7) = 0.
       !
       !
       !
       iPE = -1
       iPEmax = -1
       iBLKmax = -1

       BLOCKS: do iBLK = 1,nBlockMax

          if (unusedBLK(iBLK)) CYCLE BLOCKS

          ! Check if block is within the search region
          if ( x_BLK(nI,1,1,iBLK) > xLineMin .and. &
               x_BLK(nI,1,1,iBLK) < xLineMax .and. &
               YYR(jj) >= y_BLK(1,1,1,iBLK)  - 0.5*dy_BLK(iBLK) .and. &
               YYR(jj) <  y_BLK(1,nJ,1,iBLK) + 0.5*dy_BLK(iBLK) .and. &
               0.25*dz_BLK(iBLK) >=z_BLK(1,1,1,iBLK)  - 0.5*dz_BLK(iBLK) .and.&
               0.25*dz_BLK(iBLK) < z_BLK(1,1,nK,iBLK) + 0.5*dz_BLK(iBLK) ) then

             ! Find the j index corresponding to YYR
             yj1 = (YYR(jj)-y_BLK(1,1,1,iBLK))/dy_BLK(iBLK)  + 1
             j1  = max(1,floor(yj1))

             do i=nI,1,-1
                ! Find reversal in Bz
                if ( B0_DGB(z_,i+1,j1,1,iBLK) + &
                     State_VGB(Bz_,i+1,j1,1,iBLK) > 0. .and. &
                     B0_DGB(z_,i-1,j1,1,iBLK) + &
                     State_VGB(Bz_,i-1,j1,1,iBLK) <=0.) then
                   
                   XX0max  = x_BLK(i,j1,1,iBLK)
                   iX      = i
                   iBLKmax = iBLK
                   jY      = j1

                   EXIT BLOCKS
                endif
             enddo

          end if
       enddo BLOCKS

       ! Get the X coordinate of the first reconnection line
       call MPI_ALLREDUCE(XX0max, XX0max_all, 1, MPI_REAL,MPI_MAX,iComm,iError)
       XX0(jj) = XX0max_all

       ! Calculate the reconnection box size and electric field
       if (XX0max == XX0max_all .and. XX0max_all > xRecoMin) then
          iPE = iProc
          call set_rec_parameters(iX, jY, 1, iBlkMax)
       else
          iPE = -1
       endif

       ! Tell everyone which processor contains iBLKmax (if any)
       call MPI_ALLREDUCE(iPE, iPEmax, 1, MPI_INTEGER, MPI_MAX, iComm, iError)

       ! If there is a reconnection line, broadcast its parameters
       if (iPEmax >= 0) then
          call MPI_Bcast(RecParam_I,nRecParam,MPI_REAL,iPEmax,iComm,iError)
          BZPR_0(jj)     = RecParam_I(dBzDx_)
          VXPR_0(jj)     = RecParam_I(dVxDx_)
          DWX_0(jj)      = RecParam_I(WidthX_)
          EF0_0(jj)      = RecParam_I(Efield_) ! multiplied by 2 for run hr4
          DWZ_0(jj)      = RecParam_I(WidthZ_)
          Cos_Thet_0(jj) = RecParam_I(CosTheta_)
          Sin_Thet_0(jj) = RecParam_I(SinTheta_)
       endif

    enddo

  end subroutine set_rec_line

  !======================================================================

  subroutine set_rec_line1
    use ModMain
    use ModVarIndexes
    use ModSize
    use ModIO
    use ModPhysics
    use ModGeometry
    use ModAdvance
    use ModParallel
    use ModProcMH
    use ModConst
    use ModMpi

    implicit none

    integer :: iPE, iPEmax, iBLK, iBLKtemp, iBLKmax
    integer :: iBLKtemp1, iBLKtemp2

    integer :: i,j,k
    integer :: ii,jj
    integer :: j1, iX, jY
    integer :: iError
    real    :: XX0max, XX0max_all, yj1
    !--------------------------------------------------------------------------
    do jj=1, njj
       YYR(jj) = DyLine/2 + (jj-1)*DyLine
       XX1(jj) = UnsetX_
       XX0max = UnsetX_

       BZPR_1(jj) = -1.
       VXPR_1(jj) = 0.
       DWX_1(jj) = -1.
       DWZ_1(jj) = -1.
       EF0_1(jj) = 0.
       Cos_Thet_1(jj) = 1.
       Sin_Thet_1(jj) = 0.
       RecParam_I(1) = -1.
       RecParam_I(2) = 0.
       RecParam_I(3) = -1.
       RecParam_I(4) = 0.
       RecParam_I(5) = -1.
       RecParam_I(6) = 1.
       RecParam_I(7) = 0.
       !
       iPE = -1
       iPEmax = -1
       iBLKmax = -1

       BLOCKS: do iBLK = 1,nBlock

          if (unusedBLK(iBLK)) CYCLE

          if ( x_BLK(nI,1,1,iBLK) > xLineMin .and. &
               x_BLK(nI,1,1,iBLK) < xLineMax1 .and. &
               x_BLK(nI,1,1,iBLK) < XX0(jj) - 2. .and. &
               YYR(jj) >= y_BLK(1,1,1,iBLK) - 0.5*dy_BLK(iBLK) .and. &
               YYR(jj) < y_BLK(1,nJ,1,iBLK) + 0.5*dy_BLK(iBLK) .and. &
               0.25*dz_BLK(iBLK) >=z_BLK(1,1,1,iBLK) - 0.5*dz_BLK(iBLK) .and. &
               0.25*dz_BLK(iBLK) < z_BLK(1,1,nK,iBLK) + 0.5*dz_BLK(iBLK) ) then

             yj1=(YYR(jj)-y_BLK(1,1,1,iBLK))/dy_BLK(iBLK)+1.
             j1=max(1,floor(yj1))

             do i=nI,1,-1
                !
                if ( B0_DGB(z_,i+1,j1,1,iBLK) + &
                     State_VGB(Bz_,i+1,j1,1,iBLK) > 0 .and. &
                     B0_DGB(z_,i-1,j1,1,iBLK) + &
                     State_VGB(Bz_,i-1,j1,1,iBLK) <=0) then

                   XX0max  = x_BLK(i,j1,1,iBLK)
                   iX      = i
                   iBLKmax = iBLK
                   jY      = j1
                   EXIT BLOCKS
                endif
             enddo

          end if
       enddo BLOCKS

       ! Get the X coordinate of the second reconnection line
       call MPI_ALLREDUCE(XX0max, XX0max_all, 1, MPI_REAL,MPI_MAX,iComm,iError)
       XX1(jj) = XX0max_all

       ! Calculate the reconnection box size and electric field
       if (XX0max == XX0max_all .and. XX0max_all > xRecoMin) then
          iPE = iproc
          call set_rec_parameters(iX, jY, 1, iBlkMax)
       else
          iPE = -1
       endif
       ! Tell everyone which processor contains iBLKmax (if any)
       call MPI_ALLREDUCE(iPE, iPEmax, 1, MPI_INTEGER, MPI_MAX, iComm,iError)

       ! If there is a reconnection line, broadcast its parameters
       if (iPEmax >= 0) then
          call MPI_Bcast(RecParam_I,nRecParam,MPI_REAL,iPEmax,iComm,iError)
          BZPR_1(jj)     = RecParam_I(dBzDx_)
          VXPR_1(jj)     = RecParam_I(dVxDx_)
          DWX_1(jj)      = RecParam_I(WidthX_)
          EF0_1(jj)      = RecParam_I(Efield_) ! multiplied by 2 for run hr4
          DWZ_1(jj)      = RecParam_I(WidthZ_)
          Cos_Thet_1(jj) = RecParam_I(CosTheta_)
          Sin_Thet_1(jj) = RecParam_I(SinTheta_)
       endif

    enddo

  end subroutine set_rec_line1

  !========================================================================
  subroutine set_rec_parameters(i, j, k, iBlock)

    use ModAdvance, ONLY: State_VGB, Bx_, Bz_, Rho_, RhoUx_, RhoUz_, p_, &
         B0_DGB
    use ModMain,ONLY:x_,y_,z_
    use ModGeometry, ONLY: Dx_Blk, Dy_Blk, Dz_Blk

    integer, intent(in) :: i, j, k, iBlock

    real:: dBzDx, dBxDz, dVxDx, WidthX, WidthZ, Efield, IonSpeed
    !-----------------------------------------------------------------------
    ! dBz/dx
    dBzDx = ((State_VGB(Bz_,i+1,j,k,iBlock) &
         +      B0_DGB(z_,i+1,j,k,iBlock))  &
         -   (State_VGB(Bz_,i-1,j,k,iBlock) &
         +      B0_DGB(z_,i-1,j,k,iBlock))) / (2*Dx_BLK(iBlock))

    ! dBx/dz
    dBxDz = ((State_VGB(Bx_,i,j,k+1,iBlock) &
         +      B0_DGB(x_,i,j,k+1,iBlock))  &
         -   (State_VGB(Bx_,i,j,k-1,iBlock) &
         +      B0_DGB(x_,i,j,k-1,iBlock))) / (2*Dz_BLK(iBlock))

    ! d(Vx)/dx
    dVxDx = (State_VGB(RhoUx_,i+1,j,1,iBlock) &
         /   State_VGB(Rho_  ,i+1,j,1,iBlock) &
         -   State_VGB(RhoUx_,i-1,j,1,iBlock) &
         /   State_VGB(Rho_  ,i-1,j,1,iBlock)) / (2*Dx_BLK(iBlock))

    ! IonSpeed = sqrt(2*P/Rho)
    IonSpeed = sqrt(2*State_VGB(p_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock))

    Efield = IonMassPerCharge*IonSpeed*dVxDx

    ! WidthX = sqrt(IonMassPerCharge*sqrt(2*p/rho) / (dBz/dx))
    WidthX = 4*dx_BLK(iBlock)
    if (dBzDx > 0) WidthX = max(WidthX, sqrt(IonMassPerCharge*IonSpeed/dBzDx))

    ! WidthZ = max(4*Dz, sqrt(IonMassPerCharge/(dBx/dz) * sqrt(2*p/rho)))
    WidthZ = 4*dz_BLK(iBlock)
    if (dBxDz > 0) WidthZ = max(WidthZ, sqrt(IonMassPerCharge*IonSpeed/dBxDz))

    RecParam_I(dBzDx_)    = dBzDx
    RecParam_I(dVxDx_)    = dVxDx
    RecParam_I(Efield_)   = Efield
    RecParam_I(WidthX_)   = WidthX
    RecParam_I(WidthZ_)   = WidthZ
    RecParam_I(CosTheta_) = 1.0
    RecParam_I(SinTheta_) = 0.0

    ! The code below takes the tilt in the X-Y plane into account

    !BX = State_VGB(Bx_,:,:,:,iBlock)+B0xCell_BLK(:,:,:,iBlock)
    !BY = State_VGB(By_,:,:,:,iBlock)+B0yCell_BLK(:,:,:,iBlock)
    !BZ = State_VGB(Bz_,:,:,:,iBlock)+B0zCell_BLK(:,:,:,iBlock)
    !        B2(:,:,:) = BX**2+BY**2+BZ**2
    !Jyy=0.5*((BX(i,j,2)-BX(i,j,0))/dz_BLK(iBlock) - &
    !     (BZ(i+1,j,1)-BZ(i-1,j,1))/dx_BLK(iBlock))
    !Jxx=0.5*((BZ(i,j+1,1)-BZ(i,j-1,1))/dy_BLK(iBlock)     - &
    !     (BY(i,j,2)-BY(i,j,0))/dz_BLK(iBlock))
    !Jxy=sqrt(Jyy**2+Jxx**2)
    !CosThet=Jyy/Jxy
    !SinThet=Jxx/Jxy
    !BzPr=(BZ(i+1,j,1)-BZ(i-1,j,1))/(2.*dx_BLK(iBlock))*CosThet-&
    !        (BZ(i,j+1,1)-BZ(i,j-1,1))/(2.*dy_BLK(iBlock))*SinThet
    !
    !Rho3 = State_VGB(rho_,i,j-1,1,iBlock)
    !Rho4 = State_VGB(rho_,i,j+1,1,iBlock)
    !VYPRX=(State_VGB(rhoUy_,i+1,j,1,iBlock)/Rho2 - &
    !     State_VGB(rhoUy_,i-1,j,1,iBlock)/Rho1)/(2.*dx_BLK(iBlock)) 
    !VXPRY=(State_VGB(rhoUx_,i,j+1,1,iBlock)/Rho4 - &
    !     State_VGB(rhoUx_,i,j-1,1,iBlock)/Rho3)/(2.*dy_BLK(iBlock))
    !VYPRY=(State_VGB(rhoUy_,i,j+1,1,iBlock)/Rho4 - &
    !     State_VGB(rhoUy_,i,j-1,1,iBlock)/Rho3)/(2.*dy_BLK(iBlock))
    !RecParam_I(2)=VXPRX*CosThet**2-CosThet*SinThet*(VYPRX+VXPRY) &
    !       +VYPRY*SinThet**2
    !write(*,*)'VXPRX=',VXPRX,'VYPRX=',VYPRY,'VXPRY=',VXPRY, &
    !       'VYPRY=',VYPRY

  end subroutine set_rec_parameters

end module ModUser

