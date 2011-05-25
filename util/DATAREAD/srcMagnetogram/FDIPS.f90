module ModPotentialField
  
  use ModMpi

  implicit none

  ! input parameter
  logical:: DoReadMagnetogram = .true.

  ! grid and domain parameters
  integer:: nR = 150, nThetaAll = 180, nPhiAll = 360
  real   :: rMin = 1.0, rMax = 2.5

  ! domain decomposition
  integer:: nProcTheta = 2, nProcPhi = 2

  ! solver parameters
  logical:: UsePreconditioner = .true.
  real   :: Tolerance         = 1e-10

  ! magnetogram parameters
  character(len=100):: NameFileIn = 'fitsfile.dat'  ! filename
  logical           :: UseCosTheta = .true. 
  real              :: BrMax = 3500.0               ! Saturation level of MDI

  ! output paramters
  logical           :: DoSaveField   = .true.
  character(len=100):: NameFileField = 'potentialfield'
  character(len=5)  :: TypeFileField = 'real8'

  logical           :: DoSavePotential   = .true.
  character(len=100):: NameFilePotential = 'potentialtest'
  character(len=5)  :: TypeFilePotential = 'real8'
  
  logical           :: DoSaveTecplot   = .false.
  character(len=100):: NameFileTecplot = 'potentialfield.dat'

  ! testing parameters
  logical :: UseTiming = .true.
  real    :: TimeStart, TimeEnd

  logical :: DoTestMe  = .true.
  integer :: iRTest = 1, iPhiTest = 1, iThetaTest = 2

  ! local variables --------------------
  character(len=100):: NameFile

  logical :: UseBr = .true.

  real, dimension(:), allocatable :: &
       Radius_I, Theta_I, Phi_I, SinTheta_I, &
       dRadius_I, dPhi_I, dCosTheta_I, &
       RadiusNode_I, ThetaNode_I, PhiNode_I, SinThetaNode_I, &
       dRadiusNode_I, dTheta_I, dThetaNode_I, dPhiNode_I, dCosThetaNode_I

  real, allocatable:: Br_II(:,:), Potential_C(:,:,:), Rhs_C(:,:,:), &
       B0_DF(:,:,:,:), DivB_C(:,:,:), PlotVar_VC(:,:,:,:), BrLocal_II(:,:)

  ! Variables for hepta preconditioner
  real, parameter:: PrecondParam = 1.0 ! see ModLinearSolver

  ! Seven diagonals for the preconditioner
  real, dimension(:), allocatable :: &
       d_I, e_I, e1_I, e2_I, f_I, f1_I, f2_I

  integer, parameter :: iComm = MPI_COMM_WORLD
  integer :: iProc, nProc, iProcTheta, iProcPhi
  integer :: iTheta0, iPhi0
  integer :: nTheta, nPhi
  real,  allocatable :: tmpXPhi0_II(:,:),tmpXPhipi_II(:,:)
  integer :: nThetaLgr,nThetaSml,nPhiLgr,nPhiSml
  integer :: nProcThetaLgr,nProcThetaSml,nProcPhiLgr,nProcPhiSml

  real, allocatable :: &
       sendBC010_II(:,:), sendBC180_II(:,:), sendBC12_II(:,:), &
       sendBC21_II(:,:),  sendBC34_II(:,:), sendBC43_II(:,:), &
       sendBC020_II(:,:), sendBC360_II(:,:), &
       recvBCLgr010_II(:,:), recvBCSml010_II(:,:), &
       recvBCLgr180_II(:,:), recvBCSml180_II(:,:), &
       recvBC12_II(:,:), recvBC21_II(:,:), &
       recvBC34_II(:,:), recvBC43_II(:,:), &
       recvBC020_II(:,:), recvBC360_II(:,:)

contains

  !===========================================================================
  subroutine read_fdips_param

    use ModReadParam

    character(len=lStringLine) :: NameCommand
    character(len=10):: TypeOutput
    integer:: i

    character(len=*), parameter:: NameSub = 'read_fdips_param'
    !-----------------------------------------------------------------------
    call read_file('FDIPS.in')
    call read_init

    ! Default decomposition
    nProcTheta = nProc
    nProcPhi   = 1

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case("#DOMAIN")
          call read_var('rMin', rMin)
          call read_var('rMax', rMax)
       case("#GRID")
          call read_var('nR    ', nR)
          call read_var('nThetaAll', nThetaAll)
          call read_var('nPhiAll  ', nPhiAll)
       case("#PARALLEL")
          call read_var('nProcTheta', nProcTheta)
          call read_var('nProcPhi'  , nProcPhi)
       case("#MAGNETOGRAM")
          call read_var('NameFileIn' , NameFileIn)
          call read_var('UseCosTheta', UseCosTheta)
          call read_var('BrMax'      , BrMax)
       case("#TIMING")
          call read_var('UseTiming', UseTiming)
       case("#TEST")
          call read_var('DoTestMe', DoTestMe)
       case("#TESTIJK")
          call read_var('iRTest'    , iRTest)
          call read_var('iPhiTest'  , iPhiTest)
          call read_var('iThetaTest', iThetaTest)
       case("#SOLVER")
          call read_var('UsePreconditioner', UsePreconditioner)
          call read_var('Tolerance',         Tolerance)
       case("#OUTPUT")
          call read_var('TypeOutput', TypeOutput, IsLowerCase=.true.)
          select case(TypeOutput)
          case('field')
             DoSaveField = .true.
             call read_var('NameFileField', NameFileField)
             call read_var('TypeFileField', TypeFileField)
             ! remove .out extension if present
             i = index(NameFileField,'.out')
             if(i>0) NameFileField = NameFileField(1:i-1)
          case('potential')
             DoSavePotential = .true.
             call read_var('NameFilePotential', NameFilePotential)
             call read_var('TypeFilePotential', TypeFilePotential)
             ! remove .out extension if present
             i = index(NameFilePotential,'.out')
             if(i>0) NameFilePotential = NameFilePotential(1:i-1)
          case('tecplot')
             if(nProc > 1)call CON_stop(NameSub// &
                  ': TypeOutput=tecplot works for serial runs only')
             DoSaveTecplot = .true.
             call read_var('NameFileTecplot', NameFileTecplot)
          case default
             call CON_stop(NameSub//': unknown TypeOutput='//trim(TypeOutput))
          end select
       case default
          call CON_stop(NameSub//': unknown command='//trim(NameCommand))
       end select
    end do

    if ( nProcTheta*nProcPhi /= nProc .and. iProc==0) then
       write(*,*)NameSub,': nProcTheta, nProcPhi, nProc=', &
            nProcTheta, nProcPhi, nProc
       call CON_stop(NameSub//': nProc should be nProcTheta*nProcPhi')
    end if

    ! Do timing on proc 0 only, if at all
    if(iProc > 0) UseTiming = .false.
    
    ! Test on proc 0 only, if at all
    if(iProc > 0) DoTestMe = .false.

  end subroutine read_fdips_param
  !===========================================================================
  subroutine read_magnetogram

    use ModIoUnit, ONLY: UnitTmp_

    ! Read the raw magnetogram file into a 2d array

    integer:: iError
    integer:: nCarringtonRotation
    integer:: nTheta0, nPhi0, nThetaRatio, nPhiRatio
    integer:: iTheta, iPhi, iTheta0, iTheta1, jPhi0, jPhi1, jPhi, kPhi
    real :: BrAverage, Weight
    character (len=100) :: String

    real, allocatable:: Br0_II(:,:)

    character(len=*), parameter:: NameSub = 'read_magnetogram'
    !------------------------------------------------------------------------
    open(UnitTmp_, file=NameFileIn, status='old', iostat=iError)
    if(iError /= 0) call CON_stop(NameSub// &
         ': could not open input file'//trim(NameFileIn))
    do 
       read(UnitTmp_,'(a)', iostat = iError ) String
       if(index(String,'#CR')>0)then
          read(UnitTmp_,*) nCarringtonRotation
       endif
       if(index(String,'#ARRAYSIZE')>0)then
          read(UnitTmp_,*) nPhi0
          read(UnitTmp_,*) nTheta0
       endif
       if(index(String,'#START')>0) EXIT
    end do
    
    write(*,*)'nCarringtonRotation, nTheta0, nPhi0: ',&
         nCarringtonRotation, nTheta0, nPhi0

    allocate(Br0_II(nTheta0,nPhi0))

    ! input file is in longitude, latitude
    do iTheta = nTheta0, 1, -1
       do iPhi = 1, nPhi0
          read(UnitTmp_,*) Br0_II(iTheta,iPhi)
          if (abs(Br0_II(iTheta,iPhi)) > BrMax) &
               Br0_II(iTheta,iPhi) = sign(BrMax, Br0_II(iTheta,iPhi))
       end do
    end do

    if(nTheta0 > nThetaAll)then
       ! Set integer coarsening ratio
       nThetaRatio = nTheta0 / nThetaAll
       nThetaAll      = nTheta0 / nThetaRatio
    else
       nThetaRatio = 1
       nThetaAll      = nTheta0
    end if

    if(nPhi0 > nPhiAll)then
       nPhiRatio = nPhi0 / nPhiAll
       nPhiAll      = nPhi0 / nPhiRatio
    else
       nPhiRatio = 1
       nPhiAll      = nPhi0
    end if

    allocate(Br_II(nThetaAll,nPhiAll))
    Br_II = 0.0

    do iPhi = 1, nPhiAll
       jPhi0 = nPhiRatio*(iPhi-1) - nPhiRatio/2 + 1
       jPhi1 = nPhiRatio*(iPhi-1) + nPhiRatio/2 + 1

       do jPhi = jPhi0, jPhi1

          if( modulo(nPhiRatio,2) == 0 .and. &
               (jPhi == jPhi0 .or. jPhi == jPhi1) )then
             ! For even coarsening ratio use 0.5 weight at the two ends
             Weight = 0.5
          else
             Weight = 1.0
          end if

          ! Apply periodicity
          kPhi = modulo(jPhi-1,nPhi0) + 1

          do iTheta = 1, nThetaAll
             iTheta0 = nThetaRatio*(iTheta-1) + 1
             iTheta1 = iTheta0 + nThetaRatio - 1
          
             Br_II(iTheta,iPhi) = Br_II(iTheta,iPhi) &
                  + Weight * sum( Br0_II(iTheta0:iTheta1, kPhi))
          end do
       end do
    end do

    Br_II = Br_II / (nThetaRatio*nPhiRatio)

    ! remove monopole
    BrAverage = sum(Br_II)/(nThetaAll*nPhiAll)
    Br_II = Br_II - BrAverage

    deallocate(Br0_II)

    close(UnitTmp_)

  end subroutine read_magnetogram

  !============================================================================

  subroutine init_potential_field

    use ModConst, ONLY: cPi, cTwoPi

    integer :: iR, iTheta, iPhi
    real:: dR, dTheta, dPhi, dZ, z
    !--------------------------------------------------------------------------

    ! The processor coordinate
    iProcTheta = floor(real(iProc)/nProcPhi)
    iProcPhi   = iProc - iProcTheta*nProcPhi

    ! Calculate the nTheta, nPhi. To distribute as even as possible,
    ! there will be two different nTheta and nPhi. 
    nThetaLgr  = ceiling(real(nThetaAll)/nProcTheta)
    nThetaSml  = floor(  real(nThetaAll)/nProcTheta)
    nPhiLgr    = ceiling(real(nPhiAll)/nProcPhi)
    nPhiSml    = floor(  real(nPhiAll)/nProcPhi)

    ! Calculate the number of processors which has large/small 
    ! local number of Theta/Phi.
    nProcThetaLgr = mod(nThetaAll,nProcTheta)
    nProcThetaSml = nProcTheta - nProcThetaLgr
    nProcPhiLgr = mod(nPhiAll,nProcPhi)
    nProcPhiSml = nProcPhi - nProcPhiLgr

    ! Test if the partitioning works
    if (iProc == 0) then
       write(*,*) 'nThetaLgr = ',nThetaLgr, 'nThetaSml = ', nThetaSml
       write(*,*) 'nPhiLgr   = ', nPhiLgr,  'nPhiSml   = ', nPhiSml
       write(*,*) 'Partitioning in nThetaAll gives: ', nThetaLgr*nProcThetaLgr + &
                                                   nThetaSml*nProcThetaSml, &
                  'Actual nThetaAll is: ', nThetaAll
       write(*,*) 'Partitioning in nPhiAll gives:   ', nPhiLgr*nProcPhiLgr +  &
                                                 nPhiSml*nProcPhiSml, &
                  'Actual nPhiAll is:   ', nPhiAll
    end if

    !Both iProcTheta and iProcPhi in the large region
    if (iProcTheta < nProcThetaLgr .and. iProcPhi < nProcPhiLgr) then
       nTheta = nThetaLgr
       nPhi   = nPhiLgr
       iTheta0 = iProcTheta* nThetaLgr
       iPhi0   = iProcPhi  * nPhiLgr
    end if

    !Only iProcTheta in the large region
    if (iProcTheta < nProcThetaLgr .and. iProcPhi >= nProcPhiLgr) then
       nTheta = nThetaLgr
       nPhi   = nPhiSml
       iTheta0 = iProcTheta  * nThetaLgr
       iPhi0   = nProcPhiLgr * nPhiLgr + (iProcPhi - nProcPhiLgr)*nPhiSml
    end if

    !Only iProcPhi in the large region
    if (iProcTheta >= nProcThetaLgr .and. iProcPhi < nProcPhiLgr) then
       nTheta = nThetaSml
       nPhi   = nPhiLgr
       iTheta0 = nProcThetaLgr * nThetaLgr + (iProcTheta - nProcThetaLgr)*nThetaSml
       iPhi0   = iProcPhi      * nPhiLgr
    end if

    !Both iProcTheta and iProcPhi in the small region
    if (iProcTheta >= nProcThetaLgr .and. iProcPhi >= nProcPhiLgr) then
       nTheta = nThetaSml
       nPhi   = nPhiSml
       iTheta0 = nProcThetaLgr*nThetaLgr &
            + (iProcTheta - nProcThetaLgr)*nThetaSml
       iPhi0   = nProcPhiLgr  *nPhiLgr   &
            + (iProcPhi   - nProcPhiLgr)  *nPhiSml
    end if

     allocate( BrLocal_II(nTheta,nPhi), &
          Radius_I(0:nR+1), Theta_I(0:nTheta+1), Phi_I(0:nPhi+1), &
          dRadius_I(nR), dPhi_I(nPhi), &
          SinTheta_I(0:nTheta+1), dTheta_I(nTheta), dCosTheta_I(nTheta), &
          SinThetaNode_I(nTheta+1), dCosThetaNode_I(nTheta+1), &
          RadiusNode_I(nR+1), ThetaNode_I(nTheta+1), PhiNode_I(nPhi+1), &
          dRadiusNode_I(nR+1), dThetaNode_I(nTheta+1), dPhiNode_I(nPhi+1) , &
          Potential_C(nR,nTheta,nPhi), &
          Rhs_C(nR,nTheta,nPhi), &
          B0_DF(3,nR+1,nTheta+1,nPhi+1), &
          DivB_C(nR,nTheta,nPhi))

     ! Set BrLocal_II, this is used in set_boundary when UseBr is true
     BrLocal_II(:,:) = Br_II(iTheta0 + 1: iTheta0 + nTheta, &
                             iPhi0   + 1: iPhi0   + nPhi)
     
    ! nR is the number of mesh cells in radial direction
    ! cell centered radial coordinate
    dR = (rMax - rMin)/nR
    do iR = 0, nR+1
       Radius_I(iR) = rMin + (iR - 0.5)*dR
    end do
    ! node based radial coordinate
    do iR = 1, nR+1
       RadiusNode_I(iR) = rMin + (iR - 1)*dR
    end do
    dRadius_I = RadiusNode_I(2:nR+1) - RadiusNode_I(1:nR)
    dRadiusNode_I = Radius_I(1:nR+1) - Radius_I(0:nR)

    if(UseCosTheta)then
       dZ = 2.0/nThetaAll

       !Set Theta_I
       do iTheta = 0, nTheta+1
          z = max(-1.0, min(1.0, 1 - (iTheta + iTheta0 - 0.5)*dZ))
          Theta_I(iTheta) = acos(z)
       end do
        
       !Set the boundary condition of Theta_I
       if (iProcTheta == 0) &
            Theta_I(0)             = -Theta_I(1)
       if (iProcTheta == nProcTheta-1) &
            Theta_I(nTheta+1) = cTwoPi - Theta_I(nTheta)

       !Set ThetaNode_I
       do iTheta = 1, nTheta + 1
          z = max(-1.0, min(1.0, 1 - (iTheta + iTheta0 -1)*dZ))
          ThetaNode_I(iTheta) = acos(z)
       end do
    else
       
       dTheta = cPi/nThetaAll
       
       !Set Theta_I
       do iTheta = 0, nTheta+1
          Theta_I(iTheta) = (iTheta  + iTheta0 - 0.5)*dTheta
       end do
 
       !Set ThetaNode_I
       do iTheta = 1, nTheta+1
          ThetaNode_I(iTheta) = (iTheta + iTheta0 - 1)*dTheta
       end do
    end if

    dTheta_I = ThetaNode_I(2:nTheta+1) - ThetaNode_I(1:nTheta)
    SinTheta_I = sin(Theta_I)
    SinThetaNode_I = sin(ThetaNode_I)

    if(UseCosTheta)then
       dCosTheta_I = dZ
       dCosThetaNode_I = dZ

       ! The definitions below work better for l=m=1 harmonics test
       !dCosTheta_I(1:nThetaAll) = SinTheta_I(1:nThetaAll)*dTheta_I
       !dCosThetaNode_I(2:nThetaAll) = SinThetaNode_I(2:nThetaAll)* &
       !     (Theta_I(2:nThetaAll)-Theta_I(1:nThetaAll-1))

    else
       dCosTheta_I(1:nTheta) = SinTheta_I(1:nTheta)*dTheta
       dThetaNode_I = dTheta
    end if

    dPhi = cTwoPi/nPhiAll
    !Set Phi_I
    do iPhi = 0, nPhi+1
       Phi_I(iPhi) = (iPhi + iPhi0 - 1)*dPhi
    end do
    
    PhiNode_I = Phi_I(1:nPhi+1) - 0.5*dPhi
    dPhi_I = PhiNode_I(2:nPhi+1) - PhiNode_I(1:nPhi)
    dPhiNode_I = Phi_I(1:nPhi+1) - Phi_I(0:nPhi)

    Potential_C       =   0.0
    Rhs_C             =   0.0

  end subroutine init_potential_field

  !============================================================================

  subroutine save_potential_field

    use ModIoUnit,      ONLY: UnitTmp_
    use ModNumConst,    ONLY: cHalfPi
    use ModPlotFile,    ONLY: save_plot_file

    integer :: iR, jR, iTheta, iPhi
    real    :: r, CosTheta, SinTheta, CosPhi, SinPhi
    real    :: Br, Btheta, Bphi
    real    :: rI, rJ, rInv
    real, allocatable :: B_DX(:,:,:,:), B_DII(:,:,:)
    integer :: iError
    integer :: iStatus_I(mpi_status_size)
    integer:: nPhiOut
    !-------------------------------------------------------------------------

    ! Only the last processors in the phi direction write out the ghost cell
    if(iProcPhi == nProcPhi - 1)then
       nPhiOut = nPhi + 1
    else
       nPhiOut = nPhi
    end if

    allocate(B_DX(3,nR+1,nPhiOut,nTheta), B_DII(3,nR+1,nTheta))

    ! Average the magnetic field to the R face centers

    ! For the radial component only the theta index changes
    do iPhi = 1, nPhi; do iTheta = 1, nTheta; do iR = 1, nR+1
       B_DX(1,iR,iPhi,nTheta+1-iTheta) = B0_DF(1,iR,iTheta,iPhi)
    end do; end do; end do

    ! Use radius as weights to average Bphi and Btheta 
    ! Also swap phi and theta components (2,3) -> (3,2)
    do iPhi = 1, nPhi; do iTheta = 1, nTheta; do iR = 1, nR

       ! Use first order approximation at lower boundary. Reduces noise.
       jR = max(1,iR-1)

       rI   = Radius_I(iR)
       rJ   = Radius_I(jR)
       rInv = 0.25/RadiusNode_I(iR)

       B_DX(2,iR,iPhi,nTheta+1-iTheta) = rInv* &
            ( rI*(B0_DF(3,iR,iTheta,iPhi) + B0_DF(3,iR,iTheta,iPhi+1)) &
            + rJ*(B0_DF(3,jR,iTheta,iPhi) + B0_DF(3,jR,iTheta,iPhi+1)) )

       B_DX(3,iR,iPhi,nTheta+1-iTheta) = rInv* &
            ( rI*(B0_DF(2,iR,iTheta,iPhi) + B0_DF(2,iR,iTheta+1,iPhi)) &
            + rJ*(B0_DF(2,jR,iTheta,iPhi) + B0_DF(2,jR,iTheta+1,iPhi)) )

    end do; end do; end do

    ! set tangential components to zero at the top
    B_DX(2:3,nR+1,:,:) = 0.0
    
    ! Apply periodicity in Phi to fill the nPhi+1 ghost cell
    if (nProcPhi > 1) then
       if (iProcPhi ==0 ) then 
          b_DII = B_DX(:,:,1,:)
          call mpi_send(b_DII, 3*(nR+1)*nTheta, MPI_REAL, &
               iProcTheta*nProcPhi + nProcPhi-1, 21, iComm,  iError)
       end if
       if (iProcPhi == nProcPhi -1) then
          call mpi_recv(b_DII, 3*(nR+1)*nTheta, MPI_REAL, &
               iProcTheta*nProcPhi , 21, iComm, iStatus_I, iError)
          B_DX(:,:,nPhiOut,:) = b_DII
       end if
    else
       B_DX(:,:,nPhiOut,:) = B_DX(:,:,1,:)
    end if

    if(DoSaveField)then
       ! Note the fake processor index to be used by redistribute.pl
       write(NameFile,'(a,2i2.2,a,i3.3,a)') &
            trim(NameFileField)//'_np01', nProcPhi, nProcTheta, '_', &
            iProcPhi + (nProcTheta - 1 - iProcTheta)*nProcPhi, '.out'

       call save_plot_file(NameFile, TypeFileIn=TypeFileField, &
            StringHeaderIn = &
            'Radius [Rs] Longitude [Rad] Latitude [Rad] B [G]', &
            nameVarIn = 'Radius Longitude Latitude Br Bphi Btheta' &
            //' Ro_PFSSM Rs_PFSSM PhiShift nRExt', &
            ParamIn_I = (/ rMin, rMax, 0.0, 0.0 /), &
            nDimIn=3, VarIn_VIII=B_DX, &
            Coord1In_I=RadiusNode_I, &
            Coord2In_I=Phi_I(1:nPhiOut), &
            Coord3In_I=cHalfPi-Theta_I(nTheta:1:-1))
    end if

    if(DoSaveTecplot)then
       open(unit = UnitTmp_, file=NameFileTecplot, status='replace')

       write (UnitTmp_, '(a)') 'Title = "'     // 'PFSSM' // '"'
       write (UnitTmp_, '(a)') &
         'Variables = ' // trim ('"X [Rs]", "Y [Rs]", "Z [Rs]","Bx [G]",'// &
	' "By [G]", "Bz [G]"')
       write(UnitTmp_, '(a)') 'ZONE T="Rectangular zone"'
       write(UnitTmp_, '(a,i6,a,i6,a,i6,a)') &
            ' I = ', nR+1, ', J=', nThetaAll, ', K=', nPhiAll+1, ', ZONETYPE=Ordered'
       write(UnitTmp_, '(a)')' DATAPACKING=POINT'
       write(UnitTmp_, '(a)')' DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE )'

       do iPhi = 1, nPhiAll+1; do iTheta = 1, nThetaAll; do iR = 1, nR+1
          Br     = B_DX(1,iR,iPhi,nThetaAll+1-iTheta)
          Btheta = B_DX(3,iR,iPhi,nThetaAll+1-iTheta)
          Bphi   = B_DX(2,iR,iPhi,nThetaAll+1-iTheta)
          r = RadiusNode_I(iR)
          SinTheta = SinTheta_I(iTheta)
          CosTheta = cos(Theta_I(iTheta))
          SinPhi   = sin(Phi_I(iPhi))
          CosPhi   = cos(Phi_I(iPhi))

          write (UnitTmp_,fmt="(6(E14.6))") &
               r*SinTheta*CosPhi, &
               r*SinTheta*SinPhi, &
               r*CosTheta, &
               Br*SinTheta*CosPhi + Btheta*CosTheta*CosPhi - Bphi*SinPhi, &
               Br*SinTheta*SinPhi + Btheta*CosTheta*SinPhi + Bphi*CosPhi, &
               Br*CosTheta        - Btheta*SinTheta
       end do; end do; end do

       close(UnitTmp_)
    end if

    deallocate(B_DX)

  end subroutine save_potential_field

  !============================================================================

  subroutine set_boundary(x_C, x_G)

    real, intent(in):: x_C(nR,nTheta,nPhi)
    real, intent(inout):: x_G(0:nR+1,0:nTheta+1,0:nPhi+1)

    integer:: iPhi, jPhi, shift
    integer:: status(mpi_status_size)
    integer:: SendRequest010, SendRequest020, SendRequest360, SendRequest180, &
         SendRequest12, SendRequest21, SendRequest34, SendRequest43, &
         RecvRequest010, RecvRequest020, RecvRequest360, RecvRequest180, &
         RecvRequest12, RecvRequest21, RecvRequest34, RecvRequest43
    integer:: jProc
    integer:: iError

    if (.not. allocated(tmpXPhi0_II)) allocate( &
         tmpXPhi0_II(0:nR+1,nPhiAll),              &
         tmpXPhipi_II(0:nR+1,nPhiAll),             &
         sendBC010_II(0:nR+1,nPhi),        &
         sendBC180_II(0:nR+1,nPhi),        &
         sendBC12_II(0:nR+1,0:nTheta+1),   &
         sendBC21_II(0:nR+1,0:nTheta+1),   &
         sendBC34_II(0:nR+1,nPhi),         &
         sendBC43_II(0:nR+1,nPhi),         &
         sendBC020_II(0:nR+1,0:nTheta+1),  &
         sendBC360_II(0:nR+1,0:nTheta+1),  &
         recvBCLgr010_II(0:nR+1,nPhiLgr),       &
         recvBCLgr180_II(0:nR+1,nPhiLgr),       &
         recvBCSml010_II(0:nR+1,nPhiSml),       &
         recvBCSml180_II(0:nR+1,nPhiSml),       &
         recvBC12_II(0:nR+1,0:nTheta+1),   &
         recvBC21_II(0:nR+1,0:nTheta+1),   &
         recvBC34_II(0:nR+1,nPhi),         &
         recvBC43_II(0:nR+1,nPhi),         &
         recvBC020_II(0:nR+1,0:nTheta+1),  &
         recvBC360_II(0:nR+1,0:nTheta+1))

    !--------------------------------------------------------------------------
    ! Current solution inside
    x_G(1:nR,1:nTheta,1:nPhi) = x_C

    ! The slope is forced to be Br at the inner boundary
    if(UseBr)then
       x_G(0,1:nTheta,1:nPhi) = x_C(1,:,:) - dRadiusNode_I(1)*BrLocal_II
    else
       x_G(0,1:nTheta,1:nPhi) = x_C(1,:,:)
    end if

    ! The potential is zero at the outer boundary
    x_G(nR+1,:,:) = -x_G(nR,:,:)

    ! Set tmpXPhi0_II and tmpXPhipi_II which used to store the global 
    ! boundary close to the pole to the root
    if (iProcTheta ==0) then
       sendBC010_II = x_G(:,1,1:nPhi)
       call MPI_ISEND(sendBC010_II, (nR+2)*nPhi,MPI_REAL, 0, 010, &
                      iComm, SendRequest010, iError)
    end if
    if (iProcTheta == nProcTheta-1) then
       sendBC180_II = x_G(:,nTheta,1:nPhi)
       call MPI_ISEND(sendBC180_II, (nR+2)*nPhi,MPI_REAL, 0, 180, &
                      iComm, SendRequest180, iError)
    end if

    ! The root update tmpXPhi0_II/tmpXPhipi_II from all processors
    if (iProc == 0) then
       do jProc=0, nProcPhi-1
          if (jProc < nProcPhiLgr) then
             shift = jProc * nPhiLgr
             
             call MPI_IRECV(recvBCLgr010_II, (nR+2)*nPhiLgr, MPI_REAL, &
                            jProc, &
                            010, iComm, RecvRequest010, iError)
             call mpi_wait(RecvRequest010, status, iError)
             tmpXPhi0_II(:, shift + 1: shift + nPhiLgr) = recvBCLgr010_II
                
             call MPI_IRECV(recvBCLgr180_II, (nR+2)*nPhiLgr, MPI_REAL, &
                            (nProcTheta-1)*nProcPhi+jProc, &
                            180, iComm, RecvRequest180, iError)
             call mpi_wait(RecvRequest180, status, iError)
             tmpXPhipi_II(:, shift + 1: shift + nPhiLgr) = recvBCLgr180_II
          else
             shift = nProcPhiLgr*nPhiLgr + (jProc - nProcPhiLgr)*nPhiSml

             call MPI_IRECV(recvBCSml010_II, (nR+2)*nPhiSml,  MPI_REAL, &
                            jProc, &
                            010, iComm, RecvRequest010, iError)
             call mpi_wait(RecvRequest010, status, iError)
             tmpXPhi0_II(:, shift + 1: shift + nPhiSml) = recvBCSml010_II

             call MPI_IRECV(recvBCSml180_II , (nR+2)*nPhiSml,  MPI_REAL, &
                            (nProcTheta-1)*nProcPhi+jProc, &
                            180, iComm, RecvRequest180, iError)
             call mpi_wait(RecvRequest180, status, iError)
             tmpXPhipi_II(:, shift + 1: shift + nPhiSml) = recvBCSml180_II
          end if
       end do
    end if

    call  MPI_bcast(tmpXPhi0_II,  (nR+2)*nPhiAll, MPI_REAL, 0,  iComm, iError)
    call  MPI_bcast(tmpXPhipi_II, (nR+2)*nPhiAll, MPI_REAL, 0,  iComm, iError)

    ! Symmetric in Theta but shifted by nPhiAll/2, be careful about the shift!
    if (iProcTheta == 0) then
       do iPhi = 1, nPhi
          jPhi = modulo(iPhi + iPhi0 - 1 + nPhiAll/2, nPhiAll) + 1
          x_G(:,0,iPhi)        = tmpXPhi0_II(:,jPhi)
       end do
    end if
    if (iProcTheta == nProcTheta-1) then
       do iPhi = 1, nPhi
          jPhi = modulo(iPhi + iPhi0 - 1 + nPhiAll/2, nPhiAll) + 1
          x_G(:,nTheta+1,iPhi) = tmpXPhipi_II(:,jPhi)
       end do
    end if

    !Update the local theta boundary
    if (iProcTheta /= nProcTheta-1) then
       sendBC34_II = x_G(:,nTheta,1:nPhi)
       call MPI_ISEND(sendBC34_II, (nR+2)*nPhi, MPI_REAL, &
                     (iProcTheta+1)*nProcPhi+iProcPhi, &
                     34, iComm, SendRequest34,  iError)
    end if
    if (iProcTheta /= 0) then
       sendBC43_II = x_G(:,1,1:nPhi)
       call MPI_ISEND(sendBC43_II, (nR+2)*nPhi, MPI_REAL, &
                     (iProcTheta-1)*nProcPhi+iProcPhi, &
                     43, iComm,  SendRequest43,  iError)
    end if
    if (iProcTheta /= nProcTheta-1) then
       call MPI_IRECV(recvBC43_II, (nR+2)*nPhi, MPI_REAL, &
                     (iProcTheta+1)*nProcPhi+iProcPhi, &
                     43, iComm,  RecvRequest43, iError)
       call mpi_wait(RecvRequest43, status, iError)
       x_G(:,nTheta+1,1:nPhi) = recvBC43_II
    end if
    if (iProcTheta /= 0) then
       call MPI_IRECV(recvBC34_II, (nR+2)*nPhi, MPI_REAL, &
                     (iProcTheta-1)*nProcPhi+iProcPhi, &
                     34, iComm,  RecvRequest34, iError)
       call mpi_wait(RecvRequest34, status, iError)
       x_G(:,0,1:nPhi) = recvBC34_II
    end if

    ! Periodic around phi
    ! Send boundary info
    if (iProcPhi == 0) then
       sendBC020_II = x_G(:,:,1)
       call MPI_ISEND(sendBC020_II, (nR+2)*(nTheta+2), MPI_REAL, &
                      iProcTheta*nProcPhi + nProcPhi-1, &
                      020, iComm, SendRequest020, iError)
    end if
    if (iProcPhi == nProcPhi-1) then
       sendBC360_II = x_G(:,:,nPhi)
       call MPI_ISEND(sendBC360_II, (nR+2)*(nTheta+2), MPI_REAL, &
                      iProcTheta*nProcPhi ,  &
                      360, iComm, SendRequest360, iError)
    end if

    ! Update boundary info
    if (iProcPhi == 0) then
       call MPI_IRECV(recvBC360_II, (nR+2)*(nTheta+2), MPI_REAL, &
                      iProcTheta*nProcPhi + nProcPhi-1, &
                      360, iComm, RecvRequest360, iError)
       call mpi_wait(RecvRequest360, status, iError)
       x_G(:,:,0) = recvBC360_II
    end if
    if (iProcPhi == nProcPhi-1) then
       call MPI_IRECV(recvBC020_II, (nR+2)*(nTheta+2), MPI_REAL, &
                      iProcTheta*nProcPhi , &
                      020, iComm, RecvRequest020, iError)
       call mpi_wait(RecvRequest020, status, iError)
       x_G(:,:,nPhi+1) = recvBC020_II
    end if

    ! Start to send and update the local boundary
    if (iProcPhi /= nProcPhi-1) then
       sendBC12_II = x_G(:,:,nPhi)
       call MPI_ISEND(sendBC12_II, (nR+2)*(nTheta+2), MPI_REAL, &
                      iProcTheta*nProcPhi + iProcPhi+1, &
                      12,  iComm, SendRequest12, iError)
    end if
    if (iProcPhi /= 0) then
       sendBC21_II = x_G(:,:,1)
       call MPI_ISEND(sendBC21_II, (nR+2)*(nTheta+2), MPI_REAL, &
                      iProcTheta*nProcPhi + iProcPhi-1, &
                      21,  iComm, SendRequest21, iError)
    end if
    if (iProcPhi /= nProcPhi-1) then
       call MPI_IRECV(recvBC21_II, (nR+2)*(nTheta+2), MPI_REAL, &
                      iProcTheta*nProcPhi + iProcPhi+1, &
                      21,  iComm, RecvRequest21, iError)
       call mpi_wait(RecvRequest21, status, iError)
       x_G(:,:,nPhi+1) = recvBC21_II
    end if
    if (iProcPhi /= 0) then
       call MPI_IRECV(recvBC12_II, (nR+2)*(nTheta+2), MPI_REAL, &
                      iProcTheta*nProcPhi + iProcPhi-1, &
                      12,  iComm, RecvRequest12, iError)
       call mpi_wait(RecvRequest12, status, iError)
       x_G(:,:,0) = recvBC12_II
    end if

  end subroutine set_boundary

end module ModPotentialField

!==============================================================================

module ModB0Matvec

  use ModPotentialField, ONLY:  iComm, iProc, nProc, iProcTheta, iProcPhi, &
       nTheta, nPhi, &
       tmpXPhi0_II,tmpXPhipi_II, &
       nThetaLgr,nThetaSml,nPhiLgr,nPhiSml, &
       nProcThetaLgr,nProcThetaSml,nProcPhiLgr,nProcPhiSml, &
       nR, nTheta, nPhi, Radius_I, SinTheta_I, &
       dRadiusNode_I, dTheta_I, dCosTheta_I, dThetaNode_I, dPhiNode_I, &
       Br_II, set_boundary, &
       UseCosTheta, RadiusNode_I, Theta_I, SinThetaNode_I, dCosThetaNode_I, &
       iRTest, iThetaTest, iPhiTest, ThetaNode_I, Phi_I, PhiNode_I

  use ModMPI
  use ModUtilities, ONLY: flush_unit
  use ModIoUnit, ONLY: STDOUT_

  implicit none

contains

  !============================================================================

  subroutine matvec(x_C, y_C, n)

    use ModPotentialField, ONLY: B0_DF, &
         UsePreconditioner, nR, nTheta, d_I, e_I, e1_I, e2_I, f_I, f1_I, f2_I
    use ModLinearSolver, ONLY: Lhepta, Uhepta

    integer, intent(in) :: n
    real, intent(in)    :: x_C(n)
    real, intent(out)   :: y_C(n)
    !--------------------------------------------------------------------------

    ! Calculate y = laplace x in two steps
    call get_gradient(x_C, B0_DF)
    call get_divergence(B0_DF, y_C)

    ! Preconditioning: y'= U^{-1}.L^{-1}.y
    if(UsePreconditioner)then
       call Lhepta(        n, 1, nR, nR*nTheta, y_C, d_I, e_I, e1_I, e2_I)
       call Uhepta(.true., n, 1, nR, nR*nTheta, y_C,      f_I, f1_I, f2_I)
    end if

  end subroutine matvec

  !============================================================================

  subroutine get_gradient(x_C, Grad_DG)

    real, intent(in):: x_C(nR,nTheta,nPhi)
    real, intent(out):: Grad_DG(3,nR+1,nTheta+1,nPhi+1)

    real, allocatable, save :: x_G(:,:,:)

    integer:: iR, iTheta, iPhi

    ! real:: r, GradExact_D(3)

    !--------------------------------------------------------------------------
    if(.not.allocated(x_G))then
       allocate(x_G(0:nR+1,0:nTheta+1,0:nPhi+1))
       ! Initialize so that corners are all set
       x_G = 0.0
    end if

    call set_boundary(x_C, x_G)
    
    ! This initialization is only for the corners
    Grad_DG = 0.0

    do iPhi = 1, nPhi
       do iTheta = 1, nTheta
          do iR = 1, nR+1
             Grad_DG(1,iR,iTheta,iPhi) = &
                  (x_G(iR,iTheta,iPhi) - x_G(iR-1,iTheta,iPhi)) &
                  / dRadiusNode_I(iR)
          end do
       end do
    end do

    if(UseCosTheta)then
       do iPhi = 1, nPhi
          do iTheta = 1, nTheta+1
             do iR = 1, nR
                Grad_DG(2,iR,iTheta,iPhi) = &
                     (x_G(iR,iTheta,iPhi) - x_G(iR,iTheta-1,iPhi)) &
                     *SinThetaNode_I(iTheta) &
                     / (Radius_I(iR)*dCosThetaNode_I(iTheta))
             end do
          end do
       end do
    else
       do iPhi = 1, nPhi
          do iTheta = 1, nTheta+1
             do iR = 1, nR
                Grad_DG(2,iR,iTheta,iPhi) = &
                     (x_G(iR,iTheta,iPhi) - x_G(iR,iTheta-1,iPhi)) &
                     / (Radius_I(iR)*dThetaNode_I(iTheta))
             end do
          end do
       end do
    end if

    do iPhi = 1, nPhi+1
       do iTheta = 1, nTheta
          do iR = 1, nR
             Grad_DG(3,iR,iTheta,iPhi) = &
                  (x_G(iR,iTheta,iPhi) - x_G(iR,iTheta,iPhi-1)) &
                  / (Radius_I(iR) &
                  *max(1e-10,SinTheta_I(iTheta))*dPhiNode_I(iPhi))
          end do
       end do
    end do

    ! Calculate discretization error for the l=m=1 harmonics
    !iR = iRTest; iPhi = iPhiTest; iTheta = iThetaTest
    !
    !r = Radius_I(iR)
    !GradExact_D  = (/ &
    !     (1+2*rMax**3/RadiusNode_I(iR)**3)/(1+2*rMax**3) &
    !     *sin(Theta_I(iTheta))*cos(Phi_I(iPhi)), &
    !     (r-rMax**3/r**2)/(1+2*rMax**3)/r &
    !     *cos(ThetaNode_I(iTheta))*cos(Phi_I(iPhi)), &
    !     -(r-rMax**3/r**2)/(1+2*rMax**3)/r*sin(PhiNode_I(iPhi)) /)
    !
    !write(*,*) 'magnetogram at test cell=', Br_II(iTheta,iPhi)
    !do iDim = 1, 3
    !   write(*,*) 'Grad, Exact, Error=', &
    !        Grad_DG(iDim,iR,iTheta,iPhi), GradExact_D(iDim), &
    !        Grad_DG(iDim,iR,iTheta,iPhi) - GradExact_D(iDim)
    !end do

  end subroutine get_gradient

  !============================================================================

  subroutine get_divergence(b_DG, DivB_C)

    use ModPotentialField, ONLY: nR, Radius_I, dRadius_I, &
         dPhi_I, SinTheta_I, dCosTheta_I, RadiusNode_I, &
         SinThetaNode_I, nTheta, nPhi

    real, intent(in) :: b_DG(3,nR+1,nTheta+1,nPhi+1)
    real, intent(out):: DivB_C(nR,nTheta,nPhi)

    ! real:: r, DivExact_D(3), Div_D(3)
    integer:: iR, iTheta, iPhi
    !--------------------------------------------------------------------------
    do iPhi = 1, nPhi
       do iTheta = 1, nTheta
          do iR = 1, nR
             DivB_C(iR,iTheta,iPhi) = &
                  ( RadiusNode_I(iR+1)**2*b_DG(1,iR+1,iTheta,iPhi)   &
                  - RadiusNode_I(iR)**2  *b_DG(1,iR  ,iTheta,iPhi) ) &
                  / (Radius_I(iR)**2 *dRadius_I(iR)) &
                  + &
                  ( SinThetaNode_I(iTheta+1)*b_DG(2,iR,iTheta+1,iPhi)   &
                  - SinThetaNode_I(iTheta)  *b_DG(2,iR,iTheta  ,iPhi) ) &
                  / (Radius_I(iR)*dCosTheta_I(iTheta)) &
                  + &
                  ( b_DG(3,iR,iTheta,iPhi+1) &
                  - b_DG(3,iR,iTheta,iPhi) ) &
                  / (Radius_I(iR)*max(1e-10,SinTheta_I(iTheta))*dPhi_I(iPhi))
          end do
       end do
    end do

    ! Calculate discretization error for the l=m=1 harmonics
    !iR = iRTest; iPhi = iPhiTest; iTheta = iThetaTest
    !r = Radius_I(iR)
    !
    !Div_D(1) = ( RadiusNode_I(iR+1)**2*b_DG(1,iR+1,iTheta,iPhi)   &
    !     - RadiusNode_I(iR)**2  *b_DG(1,iR  ,iTheta,iPhi) ) &
    !     / (Radius_I(iR)**2 *dRadius_I(iR))
    !
    !Div_D(2) = ( SinThetaNode_I(iTheta+1)*b_DG(2,iR,iTheta+1,iPhi)   &
    !     - SinThetaNode_I(iTheta)  *b_DG(2,iR,iTheta  ,iPhi) ) &
    !     / (Radius_I(iR)*dCosTheta_I(iTheta))
    !
    !Div_D(3) = ( b_DG(3,iR,iTheta,iPhi+1) - b_DG(3,iR,iTheta,iPhi) ) &
    !     / (Radius_I(iR)*SinTheta_I(iTheta)*dPhi_I(iPhi))
    !
    !DivExact_D = &
    !     (/ 2*SinTheta_I(iTheta), &
    !     (1-SinTheta_I(iTheta)**2)/SinTheta_I(iTheta), &
    !     - 1/SinTheta_I(iTheta) /)
    !
    !DivExact_D = DivExact_D &
    !     *(r-rMax**3/r**2)/(1+2*rMax**3)/r**2*cos(Phi_I(iPhi))
    !
    !do iDim = 1, 3
    !   write(*,*) 'Div_D, Exact, Error=', Div_D(iDim), DivExact_D(iDim), &
    !        Div_D(iDim) - DivExact_D(iDim)
    !end do
    !   
    !write(*,*)'testlaplace=', DivB_C(iR,iTheta,iPhi)
    !write(*,*)'location   =', maxloc(abs(DivB_C))
    !write(*,*)'max laplace=', maxval(abs(DivB_C))
    !write(*,*)'avg laplace=', sum(abs(DivB_C))/(nR*nThetaAll*nPhiAll)
    !
    !stop

  end subroutine get_divergence

end module ModB0Matvec

!==============================================================================

program potential_field

  ! Solve 3D potential field with given Br at inner boundary,
  ! radial field at outer boundary.

  use ModPotentialField
  use ModB0Matvec, ONLY: get_gradient, get_divergence, matvec
  use ModLinearSolver, ONLY: bicgstab, prehepta
  use ModPlotFile, ONLY: save_plot_file
  use ModUtilities, ONLY: flush_unit
  use ModIoUnit, ONLY: STDOUT_
  use ModMpi

  implicit none

  integer :: nIter=10000
  real    :: r, DivBMax, DivBMaxAll
  integer :: n, i, iError, iR, iPhi, iTheta
  !--------------------------------------------------------------------------

  call MPI_init(iError)
  call MPI_comm_rank(iComm,iProc,iError)
  call MPI_comm_size(iComm,nProc,iError)

  call read_fdips_param

  if(DoReadMagnetogram .and. iProc == 0) call read_magnetogram

  call MPI_bcast(nThetaAll, 1, MPI_INTEGER, 0, iComm, iError)
  call MPI_bcast(nPhiAll,   1, MPI_INTEGER, 0, iComm, iError)
  if (.not. allocated(Br_II)) allocate(Br_II(nThetaAll,nPhiAll))

  call MPI_bcast(Br_II, nThetaAll*nPhiAll, MPI_REAL, 0,  iComm, iError)

  if(UseTiming) TimeStart = mpi_wtime()

  call init_potential_field

  if(.not.DoReadMagnetogram)then
     allocate(Br_II(nThetaAll,nPhiAll))
     do iPhi = 1, nPhiAll; do iTheta = 1, nThetaAll; 
        ! magnetogram proportional to the l=m=n harmonics
        n = 1 ! or 2
       Br_II(iTheta,iPhi) = sin(Theta_I(iTheta))**n *cos(n*Phi_I(iPhi))

        ! Exact solution
        do iR = 1, nR
           r = Radius_I(iR)
           Potential_C(iR,iTheta,iPhi) = Br_II(iTheta,iPhi) &
                * (r**n - rMax**(2*n+1)/r**(n+1)) &
                / (n    + (n+1)*rMax**(2*n+1))
        end do
     end do; end do

     write(*,*)'rTest    =',Radius_I(iRTest)
     write(*,*)'PhiTest  =',Phi_I(iPhiTest)
     write(*,*)'ThetaTest=',Theta_I(iThetaTest)
     write(*,*)'BrTest   =',Br_II(iThetaTest,iPhiTest)
     write(*,*)'PotTest  =',Potential_C(iRTest,iThetaTest,iPhiTest)

  end if

  n = nR*nTheta*nPhi
  ! write(*,*) 'iProc, n', iProc, n

  if(UsePreconditioner)then

     allocate(d_I(n), e_I(n), f_I(n), e1_I(n), f1_I(n), e2_I(n), f2_I(n))

     i = 0
     do iPhi = 1, nPhi; do iTheta = 1, nTheta; do iR = 1, nR
        i = i + 1
        e_I(i)  = RadiusNode_I(iR)**2 &
             /(Radius_I(iR)**2 * dRadiusNode_I(iR) * dRadius_I(iR))

        f_I(i)  = RadiusNode_I(iR+1)**2 &
             /(Radius_I(iR)**2 * dRadiusNode_I(iR+1) * dRadius_I(iR))

        e1_I(i) = SinThetaNode_I(iTheta)**2 / &
             (Radius_I(iR)**2 * dCosThetaNode_I(iTheta)  *dCosTheta_I(iTheta))

        !e1_I(i) = 0.0

        f1_I(i) = SinThetaNode_I(iTheta+1)**2 /&
             (Radius_I(iR)**2 * dCosThetaNode_I(iTheta+1)*dCosTheta_I(iTheta))

        !f1_I(i) = 0.0

        e2_I(i) = 1/(Radius_I(iR)**2 * SinTheta_I(iTheta)**2 &
             * dPhiNode_I(iPhi) * dPhi_I(iPhi))

        !e2_I(i) = 0.0

        f2_I(i) = 1/(Radius_I(iR)**2 * SinTheta_I(iTheta)**2 &
             * dPhiNode_I(iPhi+1) * dPhi_I(iPhi))

        !f2_I(i) = 0.0

        d_I(i)  = -(e_I(i) + f_I(i) + e1_I(i) + f1_I(i) + e2_I(i) + f2_I(i))

        if(iR     == 1)      d_I(i)  = d_I(i) + e_I(i) ! inner BC
        if(iR     == 1)      e_I(i)  = 0.0
        if(iR     == nR)     d_I(i)  = d_I(i) - f_I(i) ! outer BC
        if(iR     == nR)     f_I(i)  = 0.0

        if (iProcTheta ==0) then
           if(iTheta == 1)      e1_I(i) = 0.0
          end if
        if (iProcTheta == nProcTheta-1) then
           if(iTheta == nTheta) f1_I(i) = 0.0
        end if
        
        if (iProcPhi == 0) then
           if(iPhi   == 1)      e2_I(i) = 0.0
        end if
        if (iProcPhi == nProcPhi-1) then
           if(iPhi   == nPhi)   f2_I(i) = 0.0
        end if

     end do; end do; end do

     ! A -> LU
     call prehepta(n, 1, nR, nR*nTheta, PrecondParam, &
          d_I, e_I, f_I, e1_I, f1_I, e2_I, f2_I)

  end if

  UseBr = .true.
  call matvec(Potential_C, Rhs_C, n)
  Rhs_C = -Rhs_C
  
  UseBr = .false.

  call flush_unit(STDOUT_)
  call mpi_barrier(iComm, iError)
  
  call bicgstab(matvec, Rhs_C, Potential_C, .false., n, &
       Tolerance, 'rel', nIter, iError, DoTest=DoTestMe, iCommIn=iComm)

  UseBr = .true.

  if(UseTiming) TimeEnd = mpi_wtime()
  if(iProc == 0) &
       write(*,*)'nIter, Tolerance, iError=', nIter, Tolerance, iError

  ! report maximum divb
  call get_gradient(Potential_C, B0_DF)
  call get_divergence(B0_DF, DivB_C)
  DivbMax = maxval(abs(DivB_C))
  if(nProc > 1)then
     call MPI_reduce(DivBMax, DivBMaxAll, 1, MPI_REAL, MPI_MAX, 0, &
          iComm, iError)
     if(iProc==0) DivBMax = DivBMaxAll
  end if
  if(iProc ==0)then
     write(*,*) 'max(abs(divb)) = ', DivBMax
     write(*,*) 'nProcTheta, nProcPhi=', nProcTheta, nProcPhi
  end if
  if(UseTiming) write(*,*) 'running time=', TimeEnd - TimeStart

  if(DoSavePotential)then
     allocate(PlotVar_VC(6,nR,nTheta,nPhi))
     PlotVar_VC = 0.0
     PlotVar_VC(1,:,:,:) = Potential_C
     PlotVar_VC(2,:,:,:) = &
          0.5*(B0_DF(1,1:nR,1:nTheta,1:nPhi) + &
          B0_DF(1,2:nR+1,1:nTheta,1:nPhi))
     PlotVar_VC(3,:,:,:) = &
          0.5*(B0_DF(2,1:nR,1:nTheta,1:nPhi) + &
          B0_DF(2,1:nR,2:nTheta+1,1:nPhi))
     PlotVar_VC(4,:,:,:) = &
          0.5*(B0_DF(3,1:nR,1:nTheta,1:nPhi) + &
          B0_DF(3,1:nR,1:nTheta,2:nPhi+1))
     PlotVar_VC(5,:,:,:) = DivB_C
     PlotVar_VC(6,:,:,:) = Rhs_C

     ! Note the fake processor index to be used by redistribute.pl
     write(NameFile,'(a,2i2.2,a,i3.3,a)') &
          trim(NameFilePotential)//'_np01', nProcTheta, nProcPhi,'_', &
          iProcTheta + iProcPhi*nProcTheta, '.out'

     ! Save divb, potential and RHS for testing purposes
     call save_plot_file(NameFile, TypeFileIn=TypeFilePotential, &
          StringHeaderIn='potential field', &
          NameVarIn='r theta phi pot br btheta bphi divb rhs', &
          Coord1In_I=Radius_I(1:nR), &
          Coord2In_I=Theta_I(1:nTheta), &
          Coord3In_I=Phi_I(1:nPhi), &
          VarIn_VIII=PlotVar_VC)

     deallocate(PlotVar_VC)
  end if

  call save_potential_field

  call MPI_FINALIZE(iError)

end program potential_field
!==============================================================================
subroutine CON_stop(String)

  character(len=*), intent(in):: String

  write(*,*) 'ERROR:', String
  stop

end subroutine CON_stop
!==============================================================================
