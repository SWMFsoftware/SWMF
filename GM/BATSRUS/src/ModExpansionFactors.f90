module ModExpansionFactors
  use ModMagnetogram
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
contains
  subroutine write_plot_tec
      
    integer :: ios,iPhi,iTheta
    
    open ( unit = 1, file = 'SC/Test_WSA.dat', form = 'formatted', &
         access = 'sequential', status = 'replace', iostat = ios )
    
    if ( ios /= 0 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'TECPLOT_WRITE_OPEN - Fatal error!'
       write ( *, '(a)' ) '  Could not open the output file.'
       stop
    end if
    
    write ( 1, '(a)' ) 'Title = "'     // trim ('Test_WSA_Factors') // '"'
    write ( 1, '(a)' ) &
         'Variables = ' // trim (&
         '"Phi", "Theta",  "f_s", "Fisk_factor", "Theta_b"')
    write ( 1, '(a)' ) ' '
    write ( 1, '(a,i6,a,i6,a)' ) 'Zone I = ', N_PFSSM, ', J=', N_PFSSM, &
         ', F=point' 
    
    do iTheta=0,N_PFSSM
       do iPhi=0,N_PFSSM
          write ( 1, '(2f10.3,3f5.3)' )real(iPhi)*dPhi,&
               real(iTheta)*dTheta,ExpansionFactorInv_N(0,iPhi,iTheta),&
               FiskFactor_N(0,iPhi,iTheta),ThetaB_N(0,iPhi,iTheta)
       end do
    end do
    close(1)
    
  end subroutine Write_plot_tec
   !==========================================================================
  subroutine set_expansion_factors
    real :: dS,dSMax
    real,dimension(nDim) :: R_D !The vector r,phi,theta
    real,dimension(nDim) :: BTemp_D,BSun_D,BSS_D
    real,dimension(nDim) :: RSS_D,RSun_D,RPlusEnd_D,RMinusEnd_D
    integer :: iR,iPhi,iTheta
    integer :: iBcast, iStart, nSize, iError

    ! Allocte factors arrays
    if(allocated(ExpansionFactorInv_N))deallocate(ExpansionFactorInv_N)
    allocate(ExpansionFactorInv_N(-10:N_PFSSM,0:N_PFSSM,0:N_PFSSM))
    if(allocated(FiskFactor_N))deallocate(FiskFactor_N)
    allocate(FiskFactor_N(-10:N_PFSSM,0:N_PFSSM,0:N_PFSSM))
    if(allocated(ThetaB_N))deallocate(ThetaB_N)
    allocate(ThetaB_N(-10:N_PFSSM,0:N_PFSSM,0:N_PFSSM))

    !Initalize arrays:
    ExpansionFactorInv_N=cZero
    FiskFactor_N=cZero
    ThetaB_N=cZero

    ! Set Maximum value for dS (r=Rs_PFSSM, sin(Theta)=1)
    dSMax=sqrt(dR**2+(dPhi**2+dTheta**2)*Rs_PFSSM**2)

    !Loop by theta, each processor treats a separate part of the grid
    do iTheta=iProc*nThetaPerProc,(iProc+1)*nThetaPerProc-1
       if(iTheta>N_PFSSM)EXIT !Some processors may have less amount of work
       if(iTheta==0.or.iTheta==N_PFSSM)CYCLE
       do iPhi=0,N_PFSSM
          do iR=-10,N_PFSSM
             ! Define the location of the grid point
             R_D(R_)=Ro_PFSSM+real(iR)*dR
             R_D(Phi_)=real(iPhi)*dPhi
             R_D(Theta_)=real(iTheta)*dTheta
             ! Correct the radius if inside R_0
             if(R_D(R_) <= Ro_PFSSM)then
                call interpolate_field(R_D,BTemp_D)
                if(BTemp_D(R_) == cZero)then 
                   ExpansionFactorInv_N(iR,iPhi,iTheta)=cZero
                   FiskFactor_N(iR,iPhi,iTheta)=cZero
                else
                   R_D(R_)=Ro_PFSSM+cHalf*dR!*&
                   !(/cOne,BTemp_D(Phi_)/BTemp_D(R_),&
                   ! BTemp_D(Theta_)/BTemp_D(R_)/)
                end if
             end if
             call interpolate_field(R_D,BTemp_D)
             if(BTemp_D(R_) == cZero)then 
                ExpansionFactorInv_N(iR,iPhi,iTheta)=cZero
                FiskFactor_N(iR,iPhi,iTheta)=cZero
             else
                ! Correct the radius if inside R_0
                if(R_D(R_) <= Ro_PFSSM)then
                   R_D=R_D+(Ro_PFSSM+cHalf*dR-R_D(R_))*&
                        (/cOne,BTemp_D(Phi_)/BTemp_D(R_),&
                        BTemp_D(Theta_)/BTemp_D(R_)/)
                end if

                !Integrate in the positive direction of the magnetic field
                do while (R_D(R_) <= Rs_PFSSM .and. R_D(R_) >= Ro_PFSSM)
                   !Integrate the solution per local dS
                   dS=sqrt(dR**2+&
                        (R_D(R_)*sin(R_D(Theta_))*dPhi)**2&
                        +(R_D(R_)*dTheta)**2)
                   R_D=R_D+dS*f_d(R_D+dS*cHalf*f_d(R_D))
                   call correct_angles(R_D)
                end do
                ! Save the positive end of the field line
                RPlusEnd_D = R_D
                !Integrate in the negative direction of the magnetic field
                do while (R_D(R_) <= Rs_PFSSM .and. R_D(R_) >= Ro_PFSSM)
                   !Integrate the solution per local dS
                   dS=sqrt(dR**2+&
                        (R_D(R_)*sin(R_D(Theta_))*dPhi)**2&
                        +(R_D(R_)*dTheta)**2)
                   R_D=R_D-dS*f_d(R_D-dS*cHalf*f_d(R_D))
                   call correct_angles(R_D)
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
             end if
          end do
       end do
    end do


    ! Calculate ThetaB_N
    do iTheta=iProc*nThetaPerProc,(iProc+1)*nThetaPerProc-1
       if(iTheta>N_PFSSM)EXIT !Some processors may have less amount of work
       do iPhi=0,N_PFSSM
          do iR=-10,N_PFSSM
             ! Define the location of the grid point
             R_D(R_)=Ro_PFSSM+real(iR)*dR
             R_D(Phi_)=real(iPhi)*dPhi
             R_D(Theta_)=real(iTheta)*dTheta

             call interpolate_field(R_D,BTemp_D)
             if(BTemp_D(R_) == cZero)then 
                ThetaB_N(iR,iPhi,iTheta)=cZero
             else
                ! Correct the radius if inside R_0
                if(R_D(R_) <= Ro_PFSSM)then
                   R_D=R_D+(Ro_PFSSM+cHalf*dR-R_D(R_))*&
                        (/cOne,BTemp_D(Phi_)/BTemp_D(R_),&
                        BTemp_D(Theta_)/BTemp_D(R_)/)
                end if

                !Integrate in the positive direction of the magnetic field
                do while (R_D(R_) <= Rs_PFSSM .and. R_D(R_) >= Ro_PFSSM)
                   !Integrate the solution per local dS
                   dS=sqrt(dR**2+&
                        (R_D(R_)*sin(R_D(Theta_))*dPhi)**2&
                        +(R_D(R_)*dTheta)**2)
                   R_D=R_D+dS*f_d(R_D+dS*cHalf*f_d(R_D))
                   call correct_angles(R_D)
                end do
                ! Save the positive end of the field line
                RPlusEnd_D = R_D
                !Integrate in the negative direction of the magnetic field
                do while (R_D(R_) <= Rs_PFSSM .and. R_D(R_) >= Ro_PFSSM)
                   !Integrate the solution per local dS
                   dS=sqrt(dR**2+&
                        (R_D(R_)*sin(R_D(Theta_))*dPhi)**2&
                        +(R_D(R_)*dTheta)**2)
                   R_D=R_D-dS*f_d(R_D-dS*cHalf*f_d(R_D))
                   call correct_angles(R_D)
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
             end if
          end do
       end do
    end do
    if(nProc>1)then
       do iBcast=0,nProc-1
          iStart=iBcast*nThetaPerProc
          if(iStart>N_PFSSM)EXIT
          nSize=min(nThetaPerProc,N_PFSSM+1-iStart)*(N_PFSSM+1)*(N_PFSSM+11)
          call MPI_bcast(ExpansionFactorInv_N(-10,0,iStart)&
               ,nSize,MPI_REAL,iBcast,iComm,iError)
          call MPI_bcast(FiskFactor_N(-10,0,iStart)&
               ,nSize,MPI_REAL,iBcast,iComm,iError)
          call MPI_bcast(ThetaB_N(-10,0,iStart)&
               ,nSize,MPI_REAL,iBcast,iComm,iError)
       end do
    end if

  contains

    ! This fucnction calculates the value of 
    ! F(i)= B(i)/|B|/(1,r*sin(theta),r)
    function f_d(RIn_D)
      real,dimension(nDim) :: f_d
      real,dimension(nDim),intent(in) :: RIn_D
      real,parameter::cTol=cOne/(cE9*cE1)

      !Get the vector (B_r,B_phi,B_theta)
      call interpolate_field(RIn_D,f_d)
      !      write(*,*)'f_d from function',f_d
      !Divide by the metric coefficients, to obtain
      !the vector ||B|| d (r,phi,theta)/dS along the field line

      f_d=f_d/(/cOne,RIn_D(R_)*max(sin(RIn_D(Theta_)),cTol),&
           RIn_D(R_)/)

      !      write(*,*)'f_d/metric coeff',f_d

      !Divide by some scale, to limit the displacement within the integration 
      !step
      f_d=f_d/sqrt(sum(f_d**2))
      !      write(*,*)'f_d normalized',f_d
      !           max(sqrt(f_d(R_)**2+f_d(Phi_)**2+f_d(Theta_)**2),&
      !          cOne/(cE9*cE1))
      !    if(RIn_D(R_) == cZero)then
      !       f_d=f_d/(/cOne,max(sin(RIn_D(Theta_)),cOne/(cE9*cE1))*&
      !            cOne/(cE9*cE1),cOne/(cE9*cE1)/) 
      !    else

      !    end if  
    end function f_d

    function theta_b(Phi,Theta)
      real,intent(in) :: Phi, Theta
      real,allocatable,dimension(:,:) :: Phi_IJ,Theta_IJ
      real :: theta_b

      if(allocated(Phi_IJ))deallocate(Phi_IJ)
      allocate(Phi_IJ(0:N_PFSSM,0:N_PFSSM))
      if(allocated(Theta_IJ))deallocate(Theta_IJ)
      allocate(Theta_IJ(0:N_PFSSM,0:N_PFSSM))

      Phi_IJ=cZero
      Theta_IJ=cZero

      theta_b=sqrt(minval((Phi-Phi_IJ(:,:))**2+&
           (Theta-Theta_IJ(:,:))**2,&
           mask=ExpansionFactorInv_N(0,:,:)<0.001))

    end function theta_b

  end subroutine set_expansion_factors
end module ModExpansionFactors
