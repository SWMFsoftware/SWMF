module ModSpice
  use ModKind, ONLY: Real8_

  ! Interface for the SPICE library that provides coordinate system
  ! transformation for various things. SPICE stands for:
  ! Spacecraft, Planets, Instruments, C-matrix (pointing), Events
  ! http://naif.jpl.nasa.gov/naif/aboutspice.html

  implicit none
  save
  private

  public:: spice_init            ! read in SPICE "kernels", define start time
  public:: spice_rot_matrix      ! return 3x3 rotation matrix
  public:: spice_rot_vel_matrix  ! return 6x6 matrix for position and 
  !                              ! velocity transform
  public:: spice_get_two_body_distance ! return the distance between two bodies
  !                                    ! in the solar system

  ! Number of seconds between SWMF and SPICE base times:
  ! SWMF: 1965-01-01T00:00:00
  ! SPICE: 2000-01-01T11:58:56
  ! SpiceTime = SwmfTime + DtSpiceSwmf
  real, parameter,  public::   DtSpiceSwmf = -1104494359.0 

  ! Local variables
  logical     :: DoInitialize = .true.
  real(Real8_):: tStartSpice  = 0.0

contains
  !============================================================================
  subroutine spice_init(tStart, NameDirIn)

    real(Real8_),     intent(in):: tStart
    character(len=*), intent(in), optional:: NameDirIn
    character(len=100):: NameDir

    character(len=*), parameter:: NameSub = 'spice_init'
    !--------------------------------------------------------------------------
    if(.not.DoInitialize)then
       write(*,*) NameSub,' WARNING: has been initialized'
       RETURN
    end if

    tStartSpice = tStart

    NameDir = 'Param/Spice/'
    if(present(NameDirIn))then
       if(NameDirIn /= '') NameDir = NameDirIn // '/'
    end if

    CALL FURNSH(NameDir//'naif0010.tls')
    CALL FURNSH(NameDir//'pck00010.tpc')
    CALL FURNSH(NameDir//'msgr_de405_de423s.bsp')
    CALL FURNSH(NameDir//'MSO.tf')

    DoInitialize = .false.

  end subroutine spice_init
  !============================================================================
  subroutine spice_rot_matrix(tSimulation, NameCoord1, NameCoord2, Rot_DD)
    ! Return rotation matrix from coordinate system NamCoord1 to NameCoord2 at 
    ! simulation time tSimulation

    real,             intent(in) :: tSimulation
    character(len=*), intent(in) :: NameCoord1
    character(len=*), intent(in) :: NameCoord2
    real,             intent(out):: Rot_DD(3,3)

    real(Real8_):: tSpice

    character(len=*), parameter:: NameSub = 'spice_rot_matrix'
    !--------------------------------------------------------------------------
    if(DoInitialize)call CON_stop(NameSub//': ModSpice was not initialized')
    tSpice = tStartSpice + tSimulation
    call PXFORM(NameCoord1, NameCoord2, tSpice, Rot_DD)

  end subroutine spice_rot_matrix
  !============================================================================
  subroutine spice_rot_vel_matrix(tSimulation, NameCoord1, NameCoord2, Rot_II)
    ! Return rotation for position AND velocity matrix 
    ! from coordinate system NamCoord1 to NameCoord2
    ! at simulation time tSimulation

    real,             intent(in) :: tSimulation
    character(len=*), intent(in) :: NameCoord1
    character(len=*), intent(in) :: NameCoord2
    real,             intent(out):: Rot_II(6,6)

    real(Real8_):: tSpice

    character(len=*), parameter:: NameSub = 'spice_rot_vel_matrix'
    !--------------------------------------------------------------------------
    if(DoInitialize)call CON_stop(NameSub//': ModSpice was not initialized')
    tSpice = tStartSpice + tSimulation
    call SXFORM(NameCoord1, NameCoord2, tSpice, Rot_II)

  end subroutine spice_rot_vel_matrix
  !=============================================================================
  subroutine spice_get_two_body_distance(tSimulation, NameBodyTarget, NameBodyObserver, Dist_I)
    ! Return the distance, Dist_I, between two solar system bodies based on 
    ! (NameBodyTarget and NameBodyObserver) at simulation time tSimulation

    real,             intent(in) :: tSimulation
    character(len=*), intent(in) :: NameBodyTarget
    character(len=*), intent(in) :: NameBodyObserver
    real,             intent(out):: Dist_I !Unit in km

    real(Real8_):: tSpice
    real(Real8_):: posTarget(3)   ! Position of target.
    real(Real8_):: lightTime      ! One way light time between observer and target.

    character(len=*), parameter:: NameSub = 'spice_get_two_body_distance'
    !--------------------------------------------------------------------------
    if(DoInitialize)call CON_stop(NameSub//': ModSpice was not initialized')
    tSpice = tStartSpice + tSimulation
    Call SPKPOS(NameBodyTarget,tSpice,'J2000','NONE',NameBodyObserver,posTarget,lightTime)
    Dist_I=Sqrt(Sum(posTarget(:)**2))

  end subroutine spice_get_two_body_distance
end module ModSpice
