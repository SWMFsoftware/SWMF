!==============================================================================
subroutine read_ih_buffer(yIn,zIn,State_V)

  use ModVarIndexes, ONLY: nVar, RhoUy_, RhoUz_, By_, Bz_
  use GM_couple_ih, ONLY: NameCoord, nY, nZ, yMin, yMax, zMin, zMax, State_VII
  use ModMain, ONLY: TypeCoordSystem, Time_Simulation
  use CON_axes, ONLY: transform_matrix

  implicit none

  real, intent(in) :: yIn, zIn
  real, intent(out):: State_V(nVar)

  real :: GmIh_DD(3,3), Yz_D(2), y, z
  real, save :: GmIh_II(2,2), TimeSimLast = -1.0
  integer :: I1, I2, J1, J2
  !--------------------------------------------------------------------------

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
