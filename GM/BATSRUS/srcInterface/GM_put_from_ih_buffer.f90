module ModIhBuffer

  implicit none
  save

  character(len=3) :: NameCoord
  integer          :: nY, nZ
  real             :: yMin, yMax, zMin, zMax
  real, allocatable:: State_VII(:,:,:)

end module ModIhBuffer
!=============================================================================
subroutine read_ih_buffer(yIn,zIn,State_V)

  use ModVarIndexes, ONLY: nVar, RhoUy_, RhoUz_, By_, Bz_
  use ModIhBuffer
  use ModMain, ONLY: TypeCoordSystem, Time_Simulation
  use CON_axes, ONLY: transform_matrix

  implicit none

  real, intent(in) :: yIn, zIn
  real, intent(out):: State_V(nVar)

  real :: Xyz_D(3), GmIh_DD(3,3), Yz_D(2), y, z
  real, save :: GmIh_II(2,2), TimeSimLast = -1.0
  integer :: I1, I2, J1, J2
  !---------------------------------------------------------------------------

  if(nY == 1 .or. nZ == 1)then
     ! Cannot interpolate
     State_V = State_VII(:,1,1)
  else
     if(TypeCoordSystem == NameCoord)then
        y = yIn
        z = zIn
     else
        if(Time_Simulation /= TimeSimLast)then
           ! Update GSM-GSE transformation if necessary
           GmIh_DD = transform_matrix(Time_Simulation, &
                NameCoord, TypeCoordSystem)
           GmIh_II = GmIh_DD(2:3,2:3)
           TimeSimLast = Time_Simulation
        end if
        ! Convert from GM coordinates to buffer coordinates
        Yz_D    = matmul( (/yIn, zIn/), GmIh_II )
        y = Yz_D(1)
        z = Yz_D(2)
     end if

     ! Normalize coordinates to indexes: yMin --> 1, yMax --> nY
     y = (nY-1)*(y - yMin)/(yMax - yMin) + 1
     z = (nZ-1)*(z - zMin)/(zMax - zMin) + 1

     ! Convert to indexes
     I1 = floor(y)
     J1 = floor(z)

     ! Make sure that all indexes are inside
     I1 = max(1, min(I1, nY-1))
     J1 = max(1, min(J1, nZ-1))
     I2 = I1 + 1
     J2 = J1 + 1
  
     ! Interpolate (or possibly extrapolate?)
     State_V = &
          (J2-z)*((I2-y)*State_VII(:,I1,J1) + (y-I1)*State_VII(:,I2,J1)) + &
          (z-J1)*((I2-y)*State_VII(:,I1,J2) + (y-I1)*State_VII(:,I2,J2))

     ! Debugging
     !if(abs(yIn-yTest_mod)+abs(zIn-zTest_mod)<0.001)then
     !   write(*,*)'!!! y,zTestMod = ',yTest_mod,zTest_mod
     !   write(*,*)'!!! y,zIn=',yIn,zIn
     !   write(*,*)'!!! y,z  =',y,z
     !   write(*,*)'!!! I1,J1=',I1,J1
     !   write(*,*)'!!! I2,J2=',I2,J2
     !   write(*,*)'!!! Array=',State_VII(1,I1:I2,J1:J2)
     !   write(*,*)'!!! State=',State_V(1)
     !end if
  end if

  ! Convert vector variables from buffer coordinates to GM coordinates
  ! Note: ignore the GSM velocity for now, because the GSM inertial forces
  !       are neglected in BATSRUS !!!
  if(TypeCoordSystem /= NameCoord)then
     State_V(RhoUy_:RhoUz_) = matmul(GmIh_II, State_V(RhoUy_:RhoUz_))
     State_V(By_:Bz_)       = matmul(GmIh_II, State_V(By_:Bz_))
  end if

end subroutine read_ih_buffer

!=============================================================================

subroutine GM_put_from_ih_buffer( &
     NameCoordIn, nYIn, nZIn, yMinIn, yMaxIn, zMinIn, zMaxIn, Buffer_VII)

  use ModVarIndexes
  use ModPhysics, ONLY: Si2No_V, No2Io_V, UnitX_,UnitRho_,UnitU_,UnitB_,UnitP_
  use ModMain, ONLY: TypeBc_I, west_
  use ModKind, ONLY: Real8_
  use ModIhBuffer

  implicit none

  character(len=*), intent(in) :: NameCoordIn
  integer,          intent(in) :: nYIn, nZIn
  real,             intent(in) :: yMinIn, yMaxIn, zMinIn, zMaxIn
  real,             intent(in) :: Buffer_VII(nVar, nYIn, nZIn)

  integer                      :: j, k
  character(len=*), parameter  :: NameSub = 'GM_put_from_ih_buffer.f90'
  !---------------------------------------------------------------------------
  if(.not.allocated(State_VII)) then
     ! Check coordinate system. Only GSM and GSE make sense.
     if(NameCoordIn /= 'GSM' .and. NameCoord /= 'GSE') &
          call CON_stop(NameSub//': cannot handle coord system='//NameCoordIn)
     ! Store grid information
     NameCoord = NameCoordIn
     yMin = yMinIn * Si2No_V(UnitX_)
     yMax = yMaxIn * Si2No_V(UnitX_)
     zMin = zMinIn * Si2No_V(UnitX_)
     zMax = zMaxIn * Si2No_V(UnitX_)
     nY   = nYIn
     nZ   = nZIn
     ! Allocate space
     allocate(State_VII(nVar,nY,nZ))

     ! Make sure that GM uses the IH buffer
     TypeBc_I(west_) = 'ihbuffer'

     ! Debugging
     !write(*,*)'!!! NameCoord, nY, nZ=',NameCoord,nY,nZ
     !write(*,*)'!!! yMin, yMax, zMin, zMax=',yMin, yMax, zMin, zMax
  end if

  ! Store input data
  State_VII = Buffer_VII

  ! Convert units and velocity to momentum
  do k=1,nZ; do j=1,nY
     State_VII(Rho_,j,k)          = State_VII(Rho_,j,k)    * Si2No_V(UnitRho_)
     State_VII(RhoUx_:RhoUz_,j,k) = &
          State_VII(Rho_,j,k)*State_VII(Rhoux_:RhoUz_,j,k) * Si2No_V(UnitU_)
     State_VII(Bx_:Bz_,j,k)       = State_VII(Bx_:Bz_,j,k) * Si2No_V(UnitB_)
     State_VII(P_,j,k)            = State_VII(P_,j,k)      * Si2No_V(UnitP_)
  end do; end do

  !write(*,*)'GM_put_from_ih_buffer finished'
  !write(*,*)'Rho=',State_VII(Rho_,1,1)*No2Io_V(UnitRho_)
  !write(*,*)'U=',State_VII(RhoUx_:RhoUz_,1,1)/State_VII(Rho_,1,1)*No2Io_V(UnitU_)
  !write(*,*)'B=',State_VII(Bx_:Bz_,1,1)*No2Io_V(UnitB_)
  !write(*,*)'P=',State_VII(p_,1,1)*No2Io_V(UnitP_)

end subroutine GM_put_from_ih_buffer
