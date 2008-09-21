module ModMagHarmonics
  implicit none  
    
  include 'mpif.h'

  integer :: iProc, nProc, iComm

  ! Name of input file
  character (len=*), parameter:: NameFileIn='fitsfile.dat'

  ! Name of output file
  character (len=*), parameter:: NameFileOut='harmonics.dat'

  ! This Module reads a raw (RADIAL, not LOS!!!) magnetogram data file and
  ! generates a magnetogram file in the form of spherical
  ! harmonics to be use by SWMF. 
  
  ! ************************ Data Links ********************************
  ! * MDI:   http://soi.stanford.edu/magnetic/index6.html              *
  ! * WSO:   http://wso.stanford.edu/forms/prgs.html                   *  
  ! * GONG:  http://gong.nso.edu/data/magmap/QR/mqs/                   *
  ! * SOLIS: ftp://solis.nso.edu/synoptic/level3/vsm/merged/carr-rot   *
  ! * MWO:   ftp://howard.astro.ucla.edu/pub/obs/synoptic_charts       *
  ! ********************************************************************
  real, parameter:: &
       cPi        =  3.1415926535897932384626433832795,       &
       cTwoPi     =  2*cPi,                                   &
       cHalfPi    =  0.5*cPi,                                 &
       cRadToDeg  =  180./cPi,                                &
       cDegToRad  =  cPi/180.

  real,allocatable,dimension(:,:)::g_nm,h_nm, Br_DD
  real::dR=1.0,dPhi=1.0,dTheta,dSinTheta=1.0
  integer:: nThetaPerProc
  integer:: nPhi=72, nTheta=29
  integer:: nHarmonics=180
  !----------------------------------------------------------------

  integer:: i,n,m,iTheta,iPhi,iR,mm,nn,CarringtonRotation
  real:: c_n
  real:: SinPhi,CosPhi
  real:: CosTheta,SinTheta
  real:: stuff1,stuff2,stuff3
  real:: Theta,Phi,R_PFSSM
  !----------------------------------------------------------------
  real :: SumArea,da
  real :: NormalizationFactor!,R_n
  integer :: iUnit2,NmPerProc,iBcast, iStart, nSize,iNM,ll
  real,allocatable,dimension(:) :: gArray, hArray
  integer,allocatable,dimension(:) :: nArray, mArray
  integer :: SizeOfnm, ArrPerProc,EndProc,SizeLastProc
  real, allocatable, dimension(:,:) :: p_nm, dp_nm
  !----------------------------------------------------------------
  real :: Factorial1,Factorial2,delta
  real:: SinThetaM, SinThetaM1 
  integer:: delta_m0
  real, allocatable:: FactRatio1(:)
  integer, parameter:: MaxInt=100000
  real:: Sqrt_I(MaxInt)

contains

  !=================================================================
  real function sin_latitude(iTheta)
    integer,intent(in)::iTheta
    sin_latitude=(real(iTheta)+0.5)*dSinTheta-1.0
  end function sin_latitude
  !=================================================================
  real function r_latitude(iTheta)
    integer,intent(in)::iTheta
    r_latitude=asin(sin_latitude(iTheta))
  end function r_latitude
  !=================================================================
  real function colatitude(iTheta)
    integer,intent(in)::iTheta
    colatitude=cPi*0.5-r_latitude(iTheta)
  end function colatitude
  !=================================================================
  subroutine read_raw_magnetogram

    ! Read the raw magnetogram file into a 2d array
    
    integer :: iRM,jRM,iUnit,iError,nHarmonicsIn
    character (len=100) :: line
    real, allocatable:: tempBr(:)
    !----------------------------------------------------------
    
    iUnit = 1
    open(iUnit,file=NameFileIn,status='old',iostat=iError)
    
    do 
       read(iUnit,'(a)', iostat = iError ) line
       if(index(line,'#CR')>0)then
          read(iUnit,*) CarringtonRotation
       endif
       if(index(line,'#nMax')>0)then
          read(iUnit,*) nHarmonicsIn
       endif
       if(index(line,'#ARRAYSIZE')>0)then
          read(iUnit,*) nPhi
          read(iUnit,*) nTheta
       endif
       if(index(line,'#START')>0) EXIT
    end do
    
    if(iProc==0) write(*,*)'Magnetogram size - Theta,Phi: ',nTheta,nPhi
    
    ! Setting the order on harmonics to be equal to the 
    ! latitudinal resolution.
    if(nHarmonicsIn > 0 .and. nHarmonicsIn <180)then
       nHarmonics=nHarmonicsIn
    else
       nHarmonics=min(nTheta,180)
    endif
    if(iProc==0) write(*,*)'Order of harmonics: ',nHarmonics
    
    dPhi=cTwoPi/nPhi
    dTheta=cPi/nTheta
    dSinTheta=2.0/(nTheta+1)
    
    ! Allocate the harmonic coefficients arrays
    allocate( &
         p_nm(nHarmonics+1,nHarmonics+1), dp_nm(nHarmonics+1,nHarmonics+1),&
         g_nm(nHarmonics+1,nHarmonics+1), h_nm(nHarmonics+1,nHarmonics+1), &
         FactRatio1(nHarmonics+1))

    !Allocate the magnetic field array, at the spherical grid.
    allocate( tempBr(0:nPhi*nTheta-1), Br_DD(0:nPhi-1,0:nTheta-1))
    
    p_nm = 0.0
    dp_nm = 0.0
    g_nm = 0.0
    h_nm = 0.0
    
    tempBr = 0.0 
    Br_DD  = 0.0
    
    iRM=0
    do
       read(iUnit,*,iostat=iError) tempBr(iRM)
       if (iError /= 0) EXIT
       iRM=iRM+1
    end do
    
    close(iUnit)
    
    do iRM=0,nTheta-1
       do jRM=0,nPhi-1
          ! The MDI magnetogram is saturated for magnetic field larger than 
          ! 1900 gauss.
          if (abs(tempBr(iRM*nPhi+jRM)) > 1900.0) &
               tempBr(iRM*nPhi+jRM)=1900.0*sign(1.,tempBr(iRM*nPhi+jRM))
          Br_DD(jRM,iRM) = tempBr(iRM*nPhi+jRM)
       end do
    end do
    
    deallocate(tempBr)
  end subroutine read_raw_magnetogram

  !=========================================================================
  subroutine calc_harmonics

    ! This suroutine calculates the spherical harmonics from the raw 
    ! magnetogram data

    integer :: iUnit, iError, nError
    real, parameter :: Rs_PFSSM=2.5
    !-------------------------------------------------------------------------

    if(iProc==0)write(*,*)'Calculating harmonic coefficients'

    ! Creating an array with the size of SizeOfnm (the total g_nm and h_nm) for
    ! parallel calculation of the coefficients
    SizeOfnm=1
    do iNM=1,nHarmonics
       SizeOfnm=SizeOfnm+(iNM+1)
    enddo

    ArrPerProc=int(SizeOfnm/nProc)

    ! Allocating and initializing arrays
    if(allocated(gArray))deallocate(gArray)
    allocate(gArray(0:SizeOfnm-1))
    if(allocated(hArray))deallocate(hArray)
    allocate(hArray(0:SizeOfnm-1))
    if(allocated(nArray))deallocate(nArray)
    allocate(nArray(0:SizeOfnm-1))
    if(allocated(mArray))deallocate(mArray)
    allocate(mArray(0:SizeOfnm-1))

    gArray=0.0; hArray=0.0
    nArray=0.0; mArray=0.0

    ! Create arrays for n and m for parallelization
    iNM=0
    do nn=0,nHarmonics
       iNM=iNM+nn
       do mm=0,nn
          nArray(iNM+mm)=nn
          mArray(iNM+mm)=mm
       enddo
    end do

    SizeLastProc = SizeOfnm-ArrPerProc*(nProc-1)
    if(iProc==nProc-1)then 
       EndProc=SizeOfnm-1
    else
       EndProc=(iProc+1)*ArrPerProc-1
    end if

    ! Each processor gets part of the array 
    do iNM=iProc*ArrPerProc,EndProc

       ! The proper normalization factor is (2n+1)/R_n, where R_n=n+1+n(1/Rs)**(2n+1).
       ! However, in this code the coefficients are normalized only with 2n+1 to reproduce 
       ! the coefficients provided by Stanford. The division by R_n is done after
       ! the coefficients are been read in ModMagnetogram.f90.
       ! R_n=(nArray(iNM)+1.0)+nArray(iNM)*(1.0/Rs_PFSSM)**(2*nArray(iNM)+1)
       NormalizationFactor=(2*nArray(iNM)+1)

       SumArea=0.0

       do iTheta=0,nTheta-1

          Theta=colatitude(iTheta) 
          CosTheta=cos(Theta)
          SinTheta=max(sin(Theta), 1E-9)
          da=SinTheta*dTheta*dPhi
          ! Calculate the set of Legandre polynoms (with Schmidt normalization), 
          ! for a given CosTheta,SinTheta
          ! For non-radial magnetogram (LOS), a division in SinTheta is needed.
          call calc_Legandre_polynoms

          do iPhi=0,nPhi-1
             Phi=(iPhi)*dPhi
             gArray(iNM) = gArray(iNM)+&
                  Br_DD(iPhi,iTheta)*da*p_nm(nArray(iNM)+1,mArray(iNM)+1)*&
                  cos((mArray(iNM))*Phi)!/SinTheta
             hArray(iNM) = hArray(iNM)+&
                  Br_DD(iPhi,iTheta)*da*p_nm(nArray(iNM)+1,mArray(iNM)+1)*&
                  sin((mArray(iNM))*Phi)!/SinTheta
             SumArea = SumArea+da
          end do
       end do
       gArray(iNM) = &
            NormalizationFactor*gArray(iNM)/SumArea
       hArray(iNM) = &
            NormalizationFactor*hArray(iNM)/SumArea
    end do

    ! Each processor broadcasts his part of the g_nm, h_nm arrays
    if(nProc>1)then
       do iBcast=0,nProc-1
          iStart=iBcast*ArrPerProc
          if(iStart>SizeOfnm)EXIT
          if(iBcast==nProc-1)then
             nSize = SizeLastProc
          else
             nSize =  ArrPerProc
          end if
          call MPI_bcast(gArray(iStart),nSize,MPI_REAL,iBcast,iComm,iError)
          call MPI_bcast(hArray(iStart),nSize,MPI_REAL,iBcast,iComm,iError)
       end do
    end if

    if(iProc==0)then

       iNM=0
       do nn=0,nHarmonics
          iNM=iNM+nn
          do mm=0,nn
             g_nm(nn+1,mm+1) = gArray(iNM+mm)
             h_nm(nn+1,mm+1) = hArray(iNM+mm) 
          enddo
       end do

       !\
       ! Leave out monopole (n=0) term::
       !/
       g_nm(1,1) = 0.0

       write(*,*)'Done Calculating harmonic coefficients' 

       !\
       ! Writing spherical harmonics file
       !/

       write(*,*)'Writing harmonic coefficients file, named ',NameFileOut
       iUnit = 2
       open ( unit = iUnit, &
            file = NameFileOut, &
            form = 'formatted', &
            access = 'sequential', &
            status = 'replace', iostat = iError )

       if ( iError /= 0 ) then
          write (*,*)' Could not output file ', NameFileOut
          call MPI_abort(iComm, nError, iError)
       end if

       write ( iUnit, '(a19,I3,a10,I4,a4)' ) 'Coefficients order=',nHarmonics,' center=CT',CarringtonRotation,':180'
       write ( iUnit, '(a)' ) 'Observation time'
       write ( iUnit, '(a45,I3)' ) 'B0 angle & Nmax:        0          ',nHarmonics
       write ( iUnit, * )
       write ( iUnit, * )
       write ( iUnit, * )
       write ( iUnit, * )
       write ( iUnit, * )
       write ( iUnit, '(a19,I3,a26)' ) 'Max Harmonic Order:',nHarmonics,' Units: Gauss'
       write ( iUnit, '(a)' ) ' '
       write ( iUnit, '(a)' ) '  l   m      g(uT)      h(uT)'
       write ( iUnit, '(a)' ) ' '


       do nn=0,nHarmonics
          do mm=0,nn
             write(iUnit, '(2I5,2f15.3)') nn,mm,g_nm(nn+1,mm+1),h_nm(nn+1,mm+1)
          enddo
       end do

       close(iUnit)
    end if

  end subroutine calc_harmonics
 
  !\
  ! This subroutine calculates the Legendre polynoms for a particular latitude
  !/
  subroutine calc_Legandre_polynoms

    
    ! Calculate sqrt(integer) from 1 to 10000::
     
    do m=1,MaxInt
       Sqrt_I(m) = sqrt(real(m))
    end do

    ! Calculate the ratio sqrt(2m!)/(2^m*m!)::
    
    factRatio1(:) = 0.0; factRatio1(1) = 1.0
    do m=1,nHarmonics
       factRatio1(m+1) = factRatio1(m)*Sqrt_I(2*m-1)/Sqrt_I(2*m)
    enddo

    ! Calculate polynomials with appropriate normalization
    ! for Theta_PFSSMa::

    SinThetaM  = 1.0
    SinThetaM1 = 1.0
    p_nm(:,:)  = 0.0
    dp_nm(:,:) = 0.0
    
    do m=0,nHarmonics
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
       if (m < nHarmonics) p_nm(m+2,m+1) = p_nm(m+1,m+1)*Sqrt_I(2*m+3)* &
            CosTheta
       !\
       ! Eq.(30) from Altschuler et al. 1976::
       !/
       dp_nm(m+1,m+1) = factRatio1(m+1)*Sqrt_I((2-delta_m0)*(2*m+1))*&
            m*CosTheta*SinThetaM1
       !\
       ! Eq.(31) from Altschuler et al. 1976::
       !/
       if (m < nHarmonics) &
            dp_nm(m+2,m+1) = Sqrt_I(2*m+3)*(CosTheta*&
            dp_nm(m+1,m+1)-SinTheta*p_nm(m+1,m+1))
       
       SinThetaM1 = SinThetaM
       SinThetaM  = SinThetaM*SinTheta
       
    enddo
    do m=0,nHarmonics-2; do n=m+2,nHarmonics
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
    do m=0,nHarmonics; do n=m,nHarmonics
       !\
       ! Eq.(33) from Altschuler et al. 1976::
       !/
       stuff1 = 1.0/Sqrt_I(2*n+1)
       !\
       ! Eq.(34) from Altschuler et al. 1976::
       !/
       p_nm(n+1,m+1)  = p_nm(n+1,m+1)*stuff1
       dp_nm(n+1,m+1) = dp_nm(n+1,m+1)*stuff1
    enddo; enddo
  end subroutine calc_Legandre_polynoms
     
end module ModMagHarmonics
