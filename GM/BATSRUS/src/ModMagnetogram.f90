module ModMagnetogram
  use ModNumConst
  use ModReadParam
  use ModMpi
  use ModIoUnit,   ONLY: io_unit_new
  use CON_axes,    ONLY: dLongitudeHgrDeg

  !Dependecies to be removed
  use ModProcMH,   ONLY: iProc,nProc,iComm
  use ModIO,       ONLY: iUnitOut, write_prefix


  implicit none

  save

  private !Except

  public::read_magnetogram_file
  public::get_magnetogram_field

 !Dependecies to be removed
  public::Rs_PFSSM
  !\
  ! PFSSM related variables::
  !/
  ! Maximum order of spherical harmonics
  integer:: N_PFSSM=90

  ! Number of header lines in the file
  integer:: iHead_PFSSM


  ! Name of the input file
  character (LEN=32):: File_PFSSM='mf.dat'

  ! Header string in the input file
  character (LEN=80):: Head_PFSSM=''


  ! Rs - radius of outer source surface where field is taken to be 0.
  ! Ro - radius of inner boundary for the potential
  ! H  - height of ??
  !  
  
  real:: Rs_PFSSM=2.50,Ro_PFSSM=1.0,H_PFSSM=0.0

  ! Units of the magnetic field in the file including corrections
  ! relative to the magnetic units in the SC (Gauss)
  real :: UnitB=3.0

  ! Rotation angle around z-axis, in degrees,
  ! from the coordinate system of the component
  ! towards the coordinate system of the magnetic 
  ! map in the positive direction - hence, is positive. 
  ! Is set automatically to be equal to the 
  ! H(eliographic) L(ongitude) of C(entral) M(eridian)
  ! of the M(ap) minus 180 Deg, if in the input file 
  ! Phi_Shift is negative
  !
  real:: Phi_Shift=-1.0

  !Global arrays at the magnetogram grid
  integer,parameter::nDim=3,x_=1,y_=2,z_=3,R_=1,Phi_=2,Theta_=3

  real,allocatable,private,dimension(:,:,:,:)::B_DN
  real::dR=cOne,dPhi=cOne,dTheta=cOne,dInv_D(nDim)=cOne
contains
  !===========================================================================
  subroutine read_magnetogram_file
    call read_var('Ro_PFSSM'   ,Ro_PFSSM)
    call read_var('Rs_PFSSM'   ,Rs_PFSSM)
    call read_var('H_PFSSM'    ,H_PFSSM)
    call read_var('File_PFSSM' ,File_PFSSM)
    call read_var('iHead_PFSSM',iHead_PFSSM)
    call read_var('Phi_Shift'  ,Phi_Shift)
    call read_var('UnitB'      ,UnitB)

    if (iProc==0) then
       call write_prefix
       write(iUnitOut,*) 'Norder = ',N_PFSSM
       call write_prefix
       write(iUnitOut,*) 'Entered coefficient file name :: ',File_PFSSM
       call write_prefix
       write(iUnitOut,*) 'Entered number of header lines:: ',iHead_PFSSM
    endif


    call set_magnetogram

  end subroutine read_magnetogram_file
  !===========================================================================
  subroutine set_magnetogram
    !
    !---------------------------------------------------------------------------
    ! This subroutine computes PFSS (Potential Field Source Surface)
    ! field model components in spherical coordinates at user-specified
    ! r (solar radii units, r>1) and theta, phi (both radians).
    ! The subroutine requires a file containing the spherical
    ! harmonic coefficients g(n,m) & h(n.m) obtained from a separate analysis 
    ! of a photospheric synoptic map, in a standard ("radial")
    ! format and normalization used by Stanford.
    !
    ! The PFSS field model assumes no currents in the corona and
    ! a pure radial field at the source surface, here R=2.5 Rsun
    !
    ! Get solar coefficients from Todd Hoeksema's files:
    !    1. Go to http://solar.stanford.edu/~wso/forms/prgs.html
    !    2. Fill in name and email as required
    !    3. Chose Carrington rotation (leave default 180 center longitude)
    ! For most requests of integer CRs with order < 20, result will come back
    ! immediately on the web.
    !    4. Count header lines before 1st (0,0) coefficient -this will be asked!
    !---------------------------------------------------------------------------
    ! Notes:
    !
    ! In the calling routine you must initialize one variable: istart=0 (it is a 
    ! flag used to tell the subroutine to read the coefficient file the first 
    ! time only). The first time around (DoFirst=0), the subroutine will ask for
    ! the coefficient file name, the order of the expansion to use (N_PFSSM=40 or 
    ! less*, but the coeff file can contain more orders than you use), and the 
    ! number of lines in the coefficient file header. (*note computation time 
    ! increases greatly with order used).
    !
    ! The source surface surface radius has been set at Rs=2.5*Ro in the 
    ! subroutine. PFSS fields at R>Rs are radial.(br,bthet,bphi) are the resulting
    ! components. Note the units of the B fields will differ with observatory used
    ! for the coefficients. Here we assume use of the wso coefficients so units are
    ! microT. The computation of the B fields is taken mainly from Altschuler, 
    ! Levine, Stix, and Harvey, "High Resolutin Mapping of the Magnetic Field of
    ! the Solar Corona," Solar Physics 51 (1977) pp. 345-375. It uses Schmidt
    ! normalized Legendre polynomials and the normalization is explained in the 
    ! paper. The field expansion in terms of the Schmidt normalized Pnm and dPnm's
    ! is best taken from Todd Hoeksema's notes which can be downloaded from the Web
    ! http://quake.stanford.edu/~wso/Description.ps
    ! The expansions  used to get include radial factors to make the field become
    ! purely radial at the source surface. The expans. in Altschuler et al assumed
    ! that the the coefficient g(n,m) and h(n,m) were the LOS coefficients -- but 
    ! the g(n,m) and h(n,m) now available are radial (according to Janet Luhman). 
    ! Therefore, she performs an initial correction to the g(n,m) and h(n,m) to 
    ! make them consistent with the the expansion. There is no reference for this
    ! correction.
    !---------------------------------------------------------------------------

    integer:: iError
    integer:: i,n,m,iTheta,iPhi,iR
    real:: gtemp,htemp
    real:: c_n
    real:: SinPhi,CosPhi
    real:: SinPhi_I(0:N_PFSSM),CosPhi_I(0:N_PFSSM)
    real:: CosTheta,SinTheta
    real:: stuff1,stuff2,stuff3
    real:: SumR,SumT,SumP,SumPsi
    real:: Theta,Phi,R_PFSSM
    real:: Psi_PFSSM
    integer::iBcast, iStart, nSize

    ! Weights of the spherical harmonics
    real, dimension(N_PFSSM+1,N_PFSSM+1):: g_nm,h_nm
    ! Temporary variable
    real, dimension(N_PFSSM+1):: FactRatio1
    real, dimension(N_PFSSM+1,N_PFSSM+1):: p_nm,dp_nm
    integer :: iUnit

    !\
    ! Optimization by G. Toth::
    !/
    integer, parameter:: MaxInt=10000
    real, save:: Sqrt_I(MaxInt)
    !\
    !
    !/


    real:: Xy
    real:: SinThetaM, SinThetaM1
    integer:: delta_m0
    integer::nThetaPerProc
    real, dimension(-1:N_PFSSM+2):: RoRsPower_I, RoRPower_I, rRsPower_I

    !\
    ! Initialize once g(n+1,m+1) & h(n+1,m+1) by reading a file
    ! created from Web data::
    !/ 

    !\
    ! Formats adjusted for wso CR rad coeffs::
    !/
    iUnit = io_unit_new()
    open(iUnit,file=File_PFSSM,status='old',iostat=iError)
    if (iHead_PFSSM /= 0) then
       do i=1,iHead_PFSSM
          read(iUnit,'(a)') Head_PFSSM
          if(Phi_Shift<-cTiny)call get_hlcmm	
       enddo
    endif

    !\
    ! Initialize all coefficient arrays::
    !/
    g_nm(:,:) = cZero; h_nm(:,:)  = cZero
    p_nm(:,:) = cZero; dp_nm(:,:) = cZero
    !\
    ! Read file with coefficients, g_nm and h_nm::
    !/
    do
       read(iUnit,*,iostat=iError) n,m,gtemp,htemp
       if (iError /= 0) EXIT
       if (n > N_PFSSM .or. m > N_PFSSM) CYCLE
       g_nm(n+1,m+1) = gtemp
       h_nm(n+1,m+1) = htemp
    enddo
    close(iUnit)

    !\
    ! Add correction factor for radial, not LOS, coefficients::
    ! Note old "coefficients" file are LOS, all new coeffs and 
    ! files are radial)
    !/
    do n=0,N_PFSSM
       stuff1 = cOne/real(n+1+(n/(Rs_PFSSM**(2*n+1))))
       do m=0,n
          g_nm(n+1,m+1) = g_nm(n+1,m+1)*stuff1
          h_nm(n+1,m+1) = h_nm(n+1,m+1)*stuff1
       enddo
    enddo




    if(.not.allocated(B_DN))&
         allocate(B_DN(R_:Theta_,-10:N_PFSSM,0:N_PFSSM,0:N_PFSSM))

    B_DN=cZero

    dR=(Rs_PFSSM-Ro_PFSSM)/real(N_PFSSM)
    dPhi=cTwoPi/real(N_PFSSM)
    dTheta=cPi/real(N_PFSSM)


    !\
    ! Calculate sqrt(integer) from 1 to 10000::
    !/
    do m=1,MaxInt
       Sqrt_I(m) = sqrt(real(m))
    end do
    !\
    ! Calculate the ratio sqrt(2m!)/(2^m*m!)::
    !/
    factRatio1(:) = cZero; factRatio1(1) = cOne
    do m=1,N_PFSSM
       factRatio1(m+1) = factRatio1(m)*Sqrt_I(2*m-1)/Sqrt_I(2*m)
    enddo





    nThetaPerProc=N_PFSSM/nProc+1

    do iTheta=iProc*nThetaPerProc,(iProc+1)*nThetaPerProc-1

       if(iTheta>N_PFSSM)EXIT
       Theta=iTheta*dTheta
       CosTheta=cos(Theta)
       SinTheta=max(sin(Theta), cOne/(cE9*cE1))

       !\
       ! Calculate polynomials with appropriate normalization
       ! for Theta_PFSSMa::
       !/
       SinThetaM  = cOne
       SinThetaM1 = cOne

       do m=0,N_PFSSM
          if (m == 0) then
             delta_m0 = 1
          else
             delta_m0 = 0
          endif
          !\
          ! Eq.(27) from Altschuler et al. 1976::
          !/
          p_nm(m+1,m+1) = factRatio1(m+1)*Sqrt_I((2-delta_m0)*(2*m+1))* &
               SinThetaM
          !\
          ! Eq.(28) from Altschuler et al. 1976::
          !/
          if (m < N_PFSSM) p_nm(m+2,m+1) = p_nm(m+1,m+1)*Sqrt_I(2*m+3)* &
               CosTheta
          !\
          ! Eq.(30) from Altschuler et al. 1976::
          !/
          dp_nm(m+1,m+1) = factRatio1(m+1)*Sqrt_I((2-delta_m0)*(2*m+1))*&
               m*CosTheta*SinThetaM1
          !\
          ! Eq.(31) from Altschuler et al. 1976::
          !/
          if (m < N_PFSSM) &
               dp_nm(m+2,m+1) = Sqrt_I(2*m+3)*(CosTheta*&
               dp_nm(m+1,m+1)-SinTheta*p_nm(m+1,m+1))

          SinThetaM1 = SinThetaM
          SinThetaM  = SinThetaM*SinTheta

       enddo
       do m=0,N_PFSSM-2; do n=m+2,N_PFSSM
          !\
          ! Eq.(29) from Altschuler et al. 1976::
          !/
          stuff1         = Sqrt_I(2*n+1)/Sqrt_I(n**2-m**2)
          stuff2         = Sqrt_I(2*n-1)
          stuff3         = Sqrt_I((n-1)**2-m**2)/Sqrt_I(2*n-3)
          p_nm(n+1,m+1)  = stuff1*(stuff2*CosTheta*p_nm(n,m+1)-  &
               stuff3*p_nm(n-1,m+1))
          !\
          ! Eq.(32) from Altschuler et al. 1976::
          !/
          dp_nm(n+1,m+1) = stuff1*(stuff2*(CosTheta*dp_nm(n,m+1)-&
               SinTheta*p_nm(n,m+1))-stuff3*dp_nm(n-1,m+1))
       enddo; enddo
       !\
       ! Apply Schmidt normalization::
       !/
       do m=0,N_PFSSM; do n=m,N_PFSSM
          !\
          ! Eq.(33) from Altschuler et al. 1976::
          !/
          stuff1 = cOne/Sqrt_I(2*n+1)
          !\
          ! Eq.(34) from Altschuler et al. 1976::
          !/
          p_nm(n+1,m+1)  = p_nm(n+1,m+1)*stuff1
          dp_nm(n+1,m+1) = dp_nm(n+1,m+1)*stuff1
       enddo; enddo



       do iPhi=0,N_PFSSM
          Phi=real(iPhi)*dPhi
          do m=0,N_PFSSM
             CosPhi_I(m)=cos(m*Phi)
             SinPhi_I(m)=sin(m*phi)
          end do


          do iR=-10,N_PFSSM
             !\
             ! Calculate powers of the ratios of radii
             !/
             R_PFSSM=Ro_PFSSM+dR*iR
             rRsPower_I(-1) = Rs_PFSSM/R_PFSSM 
             ! This one can have negative power.
             rRsPower_I(0)  = cOne
             RoRsPower_I(0) = cOne
             RoRPower_I(0)  = cOne
             do m=1,N_PFSSM+2
                RoRsPower_I(m) = RoRsPower_I(m-1) * (Ro_PFSSM/Rs_PFSSM)
                RoRPower_I(m)  = RoRPower_I(m-1)  * (Ro_PFSSM/R_PFSSM)
                rRsPower_I(m)  = rRsPower_I(m-1)  * (R_PFSSM /Rs_PFSSM)
             end do

             !\
             ! Initialize the values of SumR,SumT,SumP, and SumPsi::
             !/
             SumR = cZero; SumT   = cZero
             SumP = cZero; SumPsi = cZero
             !\
             ! Leave out monopole (n=0) term::
             !/
             g_nm(1,1) = cZero
             !\
             ! Calculate B for (R_PFSSM,Phi_PFSSM)::
             ! Also calculate magnetic potential Psi_PFSSM
             !/
             do m=0,N_PFSSM
                do n=m,N_PFSSM
                   !\
                   ! c_n corresponds to Todd's c_l::
                   !/
                   c_n    = -RoRsPower_I(n+2)
                   !\
                   ! Br_PFSSM = -d(Psi_PFSSM)/dR_PFSSM::
                   !/
                   stuff1 = (n+1)*RoRPower_I(n+2)-c_n*n*rRsPower_I(n-1)
                   stuff2 = g_nm(n+1,m+1)*CosPhi_I(m)+h_nm(n+1,m+1)*SinPhi_I(m)
                   SumR   = SumR + p_nm(n+1,m+1)*stuff1*stuff2
                   !\
                   ! Bt_PFSSM = -(1/R_PFSSM)*d(Psi_PFSSM)/dTheta_PFSSM::
                   !/
                   stuff1 = RoRPower_I(n+2)+c_n*rRsPower_I(n-1)
                   SumT   = SumT-dp_nm(n+1,m+1)*stuff1*stuff2
                   !\
                   ! Psi_PFSSM::
                   !/
                   SumPsi = SumPsi+R_PFSSM*p_nm(n+1,m+1)*stuff1*stuff2
                   !\
                   ! Bp_PFSSM = -(1/R_PFSSM)*d(Psi_PFSSM)/dPhi_PFSSM::
                   !/
                   stuff2 = g_nm(n+1,m+1)*SinPhi_I(m)-h_nm(n+1,m+1)*CosPhi_I(m)
                   SumP   = SumP + p_nm(n+1,m+1)*m/SinTheta*stuff1*stuff2
                enddo
             enddo
             !\
             ! Compute (Br_PFSSM,Btheta_PFSSM,Bphi_PFSSM) and Psi_PFSSM::
             !/
             Psi_PFSSM    = SumPsi
             B_DN(R_,iR,iPhi,iTheta)     = SumR
             B_DN(Phi_,iR,iPhi,iTheta) = SumP
             B_DN(Theta_,iR,iPhi,iTheta)= SumT
          end do
       end do
    end do
    if(nProc>1)then
       do iBcast=0,nProc-1
          iStart=iBcast*nThetaPerProc
          if(iStart>N_PFSSM)EXIT
          nSize=min(nThetaPerProc,N_PFSSM+1-iStart)*(N_PFSSM+1)*(N_PFSSM+11)*3
          call MPI_bcast(B_DN(1,-10,0,iStart),nSize,MPI_REAL,iBcast,iComm,iError)
       end do
    end if
  end subroutine set_magnetogram
  !=================================================================
  ! SUBROUTINE get_hlcmm
  ! Read H(eliographic) L(ongitude) of the C(entral) M(eridian) of 
  ! the M(ap). Assign Phi_Shift=HLCMM-180
  subroutine get_hlcmm
    real::HLCMM     !Heliographic Longitudee of the central meridian of map
    integer::iHLCMM !The same, but HLCMM is integer at WSO magnetograms
    integer::iErrorRead,iPosition
    iPosition=index(Head_PFSSM,'Centered')	
    if (iPosition>0.and.Phi_Shift<cZero)then	
       Head_PFSSM(1:len(Head_PFSSM)-iPosition)=&
            Head_PFSSM(iPosition+1:len(Head_PFSSM))
       iPosition=index(Head_PFSSM,':')
       Head_PFSSM(1:len(Head_PFSSM)-iPosition)=&
            Head_PFSSM(iPosition+1:len(Head_PFSSM))
       iPosition=index(Head_PFSSM,':')
       Head_PFSSM(1:len(Head_PFSSM)-iPosition)=&
            Head_PFSSM(iPosition+1:len(Head_PFSSM))
       read(Head_PFSSM,'(i3)',iostat=iErrorRead)iHLCMM
       if(iErrorRead>0)call stop_mpi(&
            'Can nod find HLCMM, '//File_PFSSM//&
            ' is not a true WSO magnetogram')
       Phi_Shift=modulo(iHLCMM-180-dLongitudeHgrDeg, 360.0) 
       if(iProc==0)then
          call write_prefix;write(iUnitOut,*)'Phi_Shift=',Phi_Shift
       end if
       return
    end if
    iPosition=index(Head_PFSSM,'Central')
    if(iPosition>0.and.Phi_Shift<cZero)then
       Head_PFSSM(1:len(Head_PFSSM)-iPosition)=&
            Head_PFSSM(iPosition+1:len(Head_PFSSM))
       iPosition=index(Head_PFSSM,':')
       Head_PFSSM(1:len(Head_PFSSM)-iPosition)=&
            Head_PFSSM(iPosition+1:len(Head_PFSSM))
       read(Head_PFSSM,*,iostat=iErrorRead)HLCMM
       if(iErrorRead>0)call stop_mpi(&
            'Can nod find HLCMM, '//File_PFSSM//&
            ' is not a true MDI magnetogram')
       Phi_Shift=modulo(HLCMM-180-dLongitudeHgrDeg, 360.0) 
       if(iProc==0)then
          call write_prefix;write(iUnitOut,*)'Phi_Shift=',Phi_Shift
       end if
    end if
  end subroutine get_hlcmm
  !==========================================================================
  subroutine interpolate_field(R_D,BMap_D)
    real,intent(in)::R_D(nDim)
    real,intent(out)::BMap_D(nDim)
    integer::Node_D(nDim)
    real::Res_D(nDim)
    integer::iDim,i,j,k  
    real::Weight_III(0:1,0:1,0:1)
    Res_D=R_D
    Res_D(R_)=Res_D(R_)-Ro_PFSSM
    Res_D=Res_D*dInv_D
    Node_D=floor(Res_D*dInv_D)
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
    do iDim=1,nDim
       BMap_D(iDim)=&
            sum(Weight_III*B_DN(iDim,&
            Node_D(R_):Node_D(R_)+1,&
            Node_D(Phi_):Node_D(Phi_)+1,&
            Node_D(Theta_):Node_D(Theta_)+1))
    end do  
  end subroutine interpolate_field
  !==========================================================================
  subroutine get_magnetogram_field(xInput,yInput,zInput,B0_D)
    real, intent(in):: xInput,yInput,zInput
    real, intent(out), dimension(3):: B0_D
  
    real:: Rin_PFSSM,Theta_PFSSM,Phi_PFSSM
    real:: BMap_D(nDim)
    real:: SinPhi,CosPhi,XY,R_PFSSM
    real:: CosTheta,SinTheta
    !--------------------------------------------------------------------------
    !\
    ! Calculate cell-centered spherical coordinates::
    !/
    Rin_PFSSM   = sqrt(xInput**2+yInput**2+zInput**2)
    !\
    ! Avoid calculating B0 inside a critical radius = 0.5*Rsun
    !/
    if (Rin_PFSSM <cOne-dR*cE1) then
       B0_D= cZero
       RETURN
    end if
    Theta_PFSSM = acos(zInput/Rin_PFSSM)
    Xy          = sqrt(xInput**2+yInput**2+cTiny**2)
    Phi_PFSSM   = atan2(yInput,xInput)
    SinTheta    = Xy/Rin_PFSSM
    CosTheta    = zInput/Rin_PFSSM
    SinPhi      = yInput/Xy
    CosPhi      = xInput/Xy
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

    call interpolate_field((/R_PFSSM,Phi_PFSSM,Theta_PFSSM/),BMap_D)
    !\
    ! Magnetic field components in global Cartesian coordinates::
    ! Set B0x::
    !/
    B0_D(1) = BMap_D(R_)*SinTheta*CosPhi+    &
         BMap_D(Theta_)*CosTheta*CosPhi-&
         BMap_D(Phi_)*SinPhi
    !\
    ! Set B0y::
    !/
    B0_D(2) =  BMap_D(R_)*SinTheta*SinPhi+    &
         BMap_D(Theta_)*CosTheta*SinPhi+&
         BMap_D(Phi_)*CosPhi
    !\
    ! Set B0z::
    !/
    B0_D(3) = BMap_D(R_)*CosTheta-           &
         BMap_D(Theta_)*SinTheta
    !\
    ! Apply field strength normalization::
    ! UnitB contains the units of the CR file relative to 1 Gauss
    ! and a possible correction factor (e.g. 1.8 or 1.7).
    !/
    B0_D=B0_D*UnitB

    if (Rin_PFSSM > Rs_PFSSM) &
         B0_D  =  B0_D*(Rs_PFSSM/Rin_PFSSM)**2
  end subroutine get_magnetogram_field
end module ModMagnetogram
