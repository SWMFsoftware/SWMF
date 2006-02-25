module SC_ModExpansionFactors
  use SC_ModMagnetogram
  implicit none
  save
 !Distribution of the solar wind model parameters: 

  real,allocatable,dimension(:,:,:)::FiskFactor_N
  ! The value of this factor at a given grid point
  ! is equal to the value of |B_R(r=r_Sun)|/|B(r=R_Sun)|}, 
  ! where B_R is taken at the "photospheric footpoint" of the 
  ! magnetic field line, passing through this point:
  !\               Grid point
  ! \R_Sun       iR,iPhi,iTheta   !R_surface
  !  -------------+---------------!
  ! /                             !
  !/  
  ! Field line predicted by the source surface model 
  ! (this is a real curved magnetic field, not a ray Theta=const
  ! The value of the Fisk factor at the considered grid point is 
  ! defined as
  ! 
  ! FiskFactor_N(iR,iPhi,iTheta)=|B_R(R=R_Sun)|/|B|, 
  !
  ! if the magnetic field line is open, 
  ! zero otherwise.

  real,allocatable,dimension(:,:,:)::ExpansionFactorInv_N
  ! The expansion factor. !!!INVERTED!!!!
  ! The value of this factor at a given grid point
  ! is equal to the value of  
  ! B(R=R_{SourceSurface}/B(R=R_{Sun})*(R_SourceSurface/R_Sun)^2
  ! where the ratio of the magnetic field intensities is taken at the 
  ! two marginal point of the magnetic field line, passing through 
  ! the considered grid point:
  !
  !\               Grid point
  ! \R_Sun       iR,iPhi,iTheta   !R_SourceSurface
  !  -------------+---------------!
  ! /                             !
  !/ 
  ! Field line predicted by the source surface model 
  ! (this is a real curved magnetic field, not a ray Theta=const.
  ! The definition of the expansion factor is taken from:
  ! Wang & Shheley, ApJ, 1990, 355:726 and from 
  ! Arge & Pizzo, JGR, 2000, 105(A5):10,465    
  ! We use the referenced definition to define the expansion factor 
  ! for open field lines. If the field line is close, the expansion
  ! factor is set to be zero.
  real,allocatable,dimension(:,:,:)::ThetaB_N

  ! Speed distribution extracted from Wang-Sheeley-Arge 
  ! model (Arge et al. 2006):
  real,allocatable,dimension(:,:,:)::WSAspeed_N
contains
  subroutine SC_write_plot_tec

    integer :: iError,iPhi,iTheta,iUnit
    iUnit=io_unit_new()
    call write_prefix;write(iUnitOut,*)'Writing PFSSM factors  output file'
    open(unit = iUnit, file = 'SC/IO2/PFSSM_Factors.dat', form = 'formatted', &
         access = 'sequential', status = 'replace', iostat = iError )

    if ( iError /= 0 ) then
       call write_prefix;write(iUnitOut, '(a)' ) ' '
       call write_prefix;write(iUnitOut, '(a)' )&
            'TECPLOT_WRITE_OPEN - Fatal error!'
       call write_prefix;write(iUnitOut, '(a)' )&
            '  Could not open the output file.'
       call SC_stop_mpi('')
    end if

    write ( iUnit, '(a)' ) 'Title = "'     // trim ('PFSSM_Br') // '"'
    write ( iUnit, '(a)' ) &
         'Variables = ' // trim (&
         '"Longitude [Deg]", "Latitude [Deg]",  "f_s",&
         &"Fisk_factor","Theta_b","WSAspeed [Km/s]"')
    write ( iUnit, '(a)' ) ' '
    write ( iUnit, '(a,i6,a,i6,a)' )&
         'Zone I = ', N_PFSSM+1, ', J=', N_PFSSM+1,&
         ', F=point' 

    do iTheta=0,N_PFSSM
       do iPhi=0,N_PFSSM
          write ( iUnit, '(6f10.3)' )real(iPhi)*dPhi/cDegToRad,&
               (cPi*cHalf-real(iTheta)*dTheta)/cDegToRad,&
               ExpansionFactorInv_N(N_PFSSM,iPhi,iTheta),&
               FiskFactor_N(N_PFSSM,iPhi,iTheta),&
               ThetaB_N(N_PFSSM,iPhi,iTheta),&
               WSAspeed_N(N_PFSSM,iPhi,iTheta)/cE3 !Output in km/s'
       end do
    end do
    close(iUnit)

  end subroutine SC_write_plot_tec
   !==========================================================================
  subroutine set_expansion_factors
    real :: dS,dSMax
    real,dimension(nDim) :: R_D !The vector r,phi,theta
    real,dimension(nDim) :: BTemp_D,BSun_D,BSS_D
    real,dimension(nDim) :: RSS_D,RSun_D,RPlusEnd_D,RMinusEnd_D
    integer :: iR,iPhi,iTheta
    integer :: iBcast, iStart, nSize, iError,iIteration
    integer :: iNorth, iSouth
    real,allocatable,dimension(:,:) :: Phi_IJ,Theta_IJ
    real :: WSAPowerIndex=2.0/7.0

    ! Allocte factors arrays
    if(allocated(ExpansionFactorInv_N))deallocate(ExpansionFactorInv_N)
    allocate(ExpansionFactorInv_N(-nRExt:N_PFSSM,0:N_PFSSM,0:N_PFSSM))
    if(allocated(FiskFactor_N))deallocate(FiskFactor_N)
    allocate(FiskFactor_N(-nRExt:N_PFSSM,0:N_PFSSM,0:N_PFSSM))
    if(allocated(ThetaB_N))deallocate(ThetaB_N)
    allocate(ThetaB_N(-nRExt:N_PFSSM,0:N_PFSSM,0:N_PFSSM))

    !Initalize arrays:
    ExpansionFactorInv_N=cZero
    FiskFactor_N=cZero
    ThetaB_N=cZero

    if(allocated(Phi_IJ))deallocate(Phi_IJ)
    allocate(Phi_IJ(0:N_PFSSM,0:N_PFSSM))
    if(allocated(Theta_IJ))deallocate(Theta_IJ)
    allocate(Theta_IJ(0:N_PFSSM,0:N_PFSSM))

    Phi_IJ=cZero
    Theta_IJ=cZero
    do iTheta=0,n_PFSSM
       do iPhi=0,n_PFSSM
          Phi_IJ(iPhi,iTheta)=iPhi*dPhi
          Theta_IJ(iPhi,iTheta)=iTheta*dTheta
       end do
    end do

    ! Set Maximum value for dS (r=Rs_PFSSM, sin(Theta)=1)
    dSMax=sqrt(dR**2+(dPhi**2+dTheta**2)*Rs_PFSSM**2)

    ! Above +- 70 deg, the values of f_s and Fisk factors are set to be 
    ! constant with the same value as of +- 70 deg (The magnetograms 
    ! is not calculated for these latitudes).
    iNorth=max(int(20.0 *cDegToRad/dTheta),1)
    iSouth=N_PFSSM-iNorth

    !Loop by theta, each processor treats a separate part of the grid
    do iTheta=iProc*nThetaPerProc,(iProc+1)*nThetaPerProc-1
       if(iTheta>N_PFSSM)EXIT !Some processors may have less amount of work
       if(iTheta<iNorth.or.iTheta>iSouth)CYCLE
       do iPhi=0,N_PFSSM
          do iR=-nRExt,N_PFSSM
             ! Define the location of the grid point
             call start_at_grid_point(iR,iPhi,iTheta)
             iIteration=0
             !Integrate in the positive direction of the magnetic field
             do while (R_D(R_) <= Rs_PFSSM .and. R_D(R_) >= Ro_PFSSM)
                !Integrate the solution per local dS
                iIteration=iIteration+1
                call advance_line_point(R_D,+1.0)
             end do
        
             ! Save the positive end of the field line
             RPlusEnd_D = R_D

             ! Define the location of the grid point
             call start_at_grid_point(iR,iPhi,iTheta)
             iIteration=0
             !Integrate in the negative direction of the magnetic field
             do while (R_D(R_) <= Rs_PFSSM .and. R_D(R_) >= Ro_PFSSM)
                !Integrate the solution per local dS
                iIteration=iIteration+1
                call advance_line_point(R_D,-1.0)             
             end do
             ! Save the negative end of the field line
             RMinusEnd_D = R_D


             ! Check if the field line end points are at the same 
             ! radius (closed or hanging field line). If the field is 
             ! closed, the inv_expansion factor and Fisk factor are 
             ! set to zero. 
             if(abs(RPlusEnd_D(R_)-RMinusEnd_D(R_)) <= dSMax)then
                ExpansionFactorInv_N(iR,iPhi,iTheta) = cZero
                FiskFactor_N(iR,iPhi,iTheta) = cZero
             else
                ! Check which end of the field line is at the photosphere
                ! and which one is at the source surface. Then get the field
                ! components for both ends and calculate the value of 
                ! the factors.
                if(RPlusEnd_D(R_) > RMinusEnd_D(R_))then
                   RSS_D  = RPlusEnd_D
                   RSun_D = RMinusEnd_D
                else
                   RSS_D  = RMinusEnd_D
                   RSun_D = RPlusEnd_D
                end if
                call interpolate_field(RSS_D,BSS_D)
                call interpolate_field(RSun_D,BSun_D)
                ! Get factors for the grid point
                ExpansionFactorInv_N(iR,iPhi,iTheta)=&
                     (sqrt(dot_product(BSS_D,BSS_D))&
                     /sqrt(dot_product(BSun_D,BSun_D)))*&
                     (Rs_PFSSM/Ro_PFSSM)**2
                FiskFactor_N(iR,iPhi,iTheta) = abs(BSun_D(R_))/&
                     sqrt(dot_product(BSun_D,BSun_D))
             end if
          end do

       end do
    end do
    if(nProc>1)then
       do iBcast=0,nProc-1
          iStart=iBcast*nThetaPerProc
          if(iStart>N_PFSSM)EXIT
          nSize=min(nThetaPerProc,N_PFSSM+1-iStart)*(N_PFSSM+1)*&
               (N_PFSSM+1+nRExt)
          call MPI_bcast(ExpansionFactorInv_N(-nRExt,0,iStart)&
               ,nSize,MPI_REAL,iBcast,iComm,iError)
          call MPI_bcast(FiskFactor_N(-nRExt,0,iStart)&
               ,nSize,MPI_REAL,iBcast,iComm,iError)
       end do
    end if
    
    do iTheta=0,iNorth-1
       ExpansionFactorInv_N(:,:,iTheta)= ExpansionFactorInv_N(:,:,iNorth)
       ExpansionFactorInv_N(:,:,N_PFSSM-iTheta)= &
            ExpansionFactorInv_N(:,:,iSouth)
       FiskFactor_N(:,:,iTheta)= FiskFactor_N(:,:,iNorth)
       FiskFactor_N(:,:,N_PFSSM-iTheta)= &
            FiskFactor_N(:,:,iSouth)
    end do

    ! Calculate ThetaB_N
    do iTheta=iProc*nThetaPerProc,(iProc+1)*nThetaPerProc-1
       if(iTheta>N_PFSSM)EXIT !Some processors may have less amount of work
       if(iTheta<iNorth.or.iTheta>iSouth)CYCLE
       do iPhi=0,N_PFSSM
          do iR=-nRExt,N_PFSSM
             ! Define the location of the grid point
             call start_at_grid_point(iR,iPhi,iTheta)
             iIteration=0
             !Integrate in the positive direction of the magnetic field
             do while (R_D(R_) <= Rs_PFSSM .and. R_D(R_) >= Ro_PFSSM)
                !Integrate the solution per local dS
                iIteration=iIteration+1
                call advance_line_point(R_D,+1.0)
             end do
             ! Save the positive end of the field line
             RPlusEnd_D = R_D
             
             call start_at_grid_point(iR,iPhi,iTheta)
             iIteration=0
             !Integrate in the negative direction of the magnetic field
             do while (R_D(R_) <= Rs_PFSSM .and. R_D(R_) >= Ro_PFSSM)
                !Integrate the solution per local dS
                iIteration=iIteration+1
                call advance_line_point(R_D,-1.0)
             end do
             ! Save the negative end of the field line
             RMinusEnd_D = R_D


             ! Check if the field line end points are at the same 
             ! radius (closed or hanging field line). If the field is 
             ! closed, the inv_expansion factor and Fisk factor are 
             ! set to zero.                
             if(abs(RPlusEnd_D(R_)-RMinusEnd_D(R_)) <= dSMax)then
                ThetaB_N(iR,iPhi,iTheta) = cZero
             else
                ! Check which end of the field line is at the photosphere
                ! and which one is at the source surface. Then get the field
                ! components for both ends and calculate the value of 
                ! the factors.
                if(RPlusEnd_D(R_) > RMinusEnd_D(R_))then
                   RSS_D  = RPlusEnd_D
                   RSun_D = RMinusEnd_D
                else
                   RSS_D  = RMinusEnd_D
                   RSun_D = RPlusEnd_D
                end if
                ! Get ThetaB_N for the grid point
                ThetaB_N(iR,iPhi,iTheta)=&
                     theta_b(RSun_D(Phi_),RSun_D(Theta_))
             end if
          end do
       end do
    end do
    if(nProc>1)then
       do iBcast=0,nProc-1
          iStart=iBcast*nThetaPerProc
          if(iStart>N_PFSSM)EXIT
          nSize=min(nThetaPerProc,N_PFSSM+1-iStart)*(N_PFSSM+1)*&
               (N_PFSSM+1+nRExt)
          call MPI_bcast(ThetaB_N(-nRExt,0,iStart)&
               ,nSize,MPI_REAL,iBcast,iComm,iError)
       end do
    end if

    do iTheta=0,iNorth-1
       ThetaB_N(:,:,iTheta)= ThetaB_N(:,:,iNorth)
       ThetaB_N(:,:,N_PFSSM-iTheta)= &
            ThetaB_N(:,:,iSouth)
    end do
    !Transform to Deg
    ThetaB_N=ThetaB_N*cRadToDeg

    ! Get WSA speed
    if(allocated(WSAspeed_N))deallocate(WSAspeed_N)
    allocate(WSAspeed_N(-nRExt:N_PFSSM,0:N_PFSSM,0:N_PFSSM)) 
    WSAspeed_N=cZero

    ! Calculate WSA speed distribution using eq. 2 in Arge et al. 2006:
    WSAspeed_N(:,:,:)=(265.0+25.0*&
         (ExpansionFactorInv_N(:,:,:)+cTiny)**WSAPowerIndex*&
         (5.0-1.1*exp(1.0-(ThetaB_N(:,:,:)/4.0)**2))**2)& !km/s so far
         *cE3         !To get the result in SI

    if(iProc==0)call SC_write_plot_tec

  contains
    !--------------------------------------------------------------------------
    subroutine advance_line_point(RInOut_D,Dir)
      real,intent(inout),dimension(nDim)::RInOut_D
      real,intent(in)::Dir
      dS=0.5*min(dR,cTen)*cTwo**(iIteration/(2*N_PFSSM))
      !To avoid the line bouncing near null points
      RInOut_D=RInOut_D+Dir*dS*f_d(&
           RInOut_D+Dir*dS*cHalf*f_d(RInOut_D))
      call correct_angles(RInOut_D)
    end subroutine advance_line_point
    !------------------------------------------------------------------
    subroutine start_at_grid_point(iR,iPhi,iTheta)
      integer,intent(in)::iR,iPhi,iTheta
      R_D(R_)=Ro_PFSSM+real(iR)*dR
      R_D(Phi_)=real(iPhi)*dPhi
      R_D(Theta_)=real(iTheta)*dTheta

      if(iTheta==0)R_D(Theta_)=cTiny
      if(iTheta==N_PFSSM)R_D(Theta_)=cPi-cTiny
      R_D(R_)=min(max(R_D(R_),Ro_PFSSM+cQuarter*dR),Rs_PFSSM-cQuarter*dR)
    end subroutine start_at_grid_point
    ! This fucnction calculates the value of 
    ! F(i)= B(i)/|B|/(1,r*sin(theta),r)
    function f_d(RIn_D)
      real,dimension(nDim) :: f_d
      real,dimension(nDim),intent(in) :: RIn_D
      real,parameter::cTol=cOne/(cE9*cE1)

      !Get the vector (B_r,B_phi,B_theta)
      call interpolate_field(RIn_D,f_d)
     
      !Divide by the metric coefficients, to obtain
      !the vector ||B|| d (r,phi,theta)/dS along the field line

      f_d=f_d/(/cOne,RIn_D(R_)*max(sin(RIn_D(Theta_)),cTol),&
           RIn_D(R_)/)

      !Divide by some scale, to limit the displacement within the integration 
      !step
      f_d=f_d/sqrt(sum(f_d**2))
    end function f_d

    function theta_b(Phi,Theta)
      real,intent(in) :: Phi, Theta
      real :: theta_b

      theta_b=sqrt(minval((Phi-Phi_IJ(:,:))**2+&
           (Theta-Theta_IJ(:,:))**2,&
           mask=ExpansionFactorInv_N(0,:,:)<0.001))

    end function theta_b

  end subroutine set_expansion_factors
end module SC_ModExpansionFactors
!=================================INTERFACE==================================
subroutine SC_set_expansion_factors
  use SC_ModExpansionFactors
  implicit none
  call set_expansion_factors
end subroutine SC_set_expansion_factors

subroutine SC_get_bernoulli_integral(xInput,yInput,zInput,Output)
  use SC_ModExpansionFactors
  implicit none
  real, intent(in):: xInput,yInput,zInput
  real, intent(out):: Output
  real:: Rin_PFSSM,Theta_PFSSM,Phi_PFSSM

  integer::Node_D(nDim)
  real::Res_D(nDim)

  real::Weight_III(0:1,0:1,0:1)
  real:: R_PFSSM

  !--------------------------------------------------------------------------
  !\
  ! Calculate cell-centered spherical coordinates::
  !/
  Rin_PFSSM   = sqrt(xInput**2+yInput**2+zInput**2)
  !\
  ! Avoid calculating inside a critical radius = 0.5*Rsun
  !/
  if (Rin_PFSSM <max(Ro_PFSSM-dR*cE1,0.90*Ro_PFSSM)) then
     Output= cZero
     RETURN
  end if
  Theta_PFSSM = acos(zInput/Rin_PFSSM)
  Phi_PFSSM   = atan2(yInput,xInput)

  !\
  ! Set the source surface radius::
  ! The inner boundary in the simulations starts at a height
  ! H_PFSSM above that of the magnetic field measurements!
  !/

  R_PFSSM =min(Rin_PFSSM+H_PFSSM, Rs_PFSSM)


  !\
  ! Transform Phi_PFSSM from the component's frame to the magnetogram's frame.
  !/

  Phi_PFSSM = Phi_PFSSM - Phi_Shift*cDegToRad

  !\
  ! Take a residual for the bi-linear interpolation
  !/
  Res_D=(/R_PFSSM,Phi_PFSSM,Theta_PFSSM/)

  !Limit a value of R:
  Res_D(R_)=max(min(Res_D(R_),Rs_PFSSM-cTiny),Ro_PFSSM-nRExt*dR+cTiny)

  Res_D(R_)=Res_D(R_)-Ro_PFSSM

  call correct_angles(Res_D)
  Res_D=Res_D*dInv_D
  Node_D=floor(Res_D)
  if(Node_D(R_)==N_PFSSM)Node_D(R_)=Node_D(R_)-1
  if(Node_D(Theta_)==N_PFSSM)Node_D(Theta_)=Node_D(Theta_)-1
  Res_D=Res_D-real(Node_D)
  if(Node_D(Phi_)==N_PFSSM)Node_D(Phi_)=0

  Weight_III(0,:,:)=cOne-Res_D(R_)
  Weight_III(1,:,:)=Res_D(R_)
  Weight_III(:,0,:)=Weight_III(:,0,:)*(cOne-Res_D(Phi_))
  Weight_III(:,1,:)=Weight_III(:,1,:)*Res_D(Phi_)
  Weight_III(:,:,0)=Weight_III(:,:,0)*(cOne-Res_D(Theta_))
  Weight_III(:,:,1)=Weight_III(:,:,1)*Res_D(Theta_)

  Output=&
       sum(Weight_III*WSASpeed_N(&
       Node_D(R_):Node_D(R_)+1,&
       Node_D(Phi_):Node_D(Phi_)+1,&
       Node_D(Theta_):Node_D(Theta_)+1))
end subroutine SC_get_bernoulli_integral
