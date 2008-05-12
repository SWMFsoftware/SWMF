module ModPhysics

  !\
  ! Physical constants and parameters.
  !/

  ! Basic units for space and time
  real :: unit_x, unit_t

  ! thermal/total energy ratio limits for Pcorrect
  real    :: Pratio_min=1.e-6,Pratio_lo=0.01, Pratio_hi=0.1 

  real :: g,gm1,gm2,gp1,inv_g,inv_gm1              !gamma and gamma derived values
  real :: cLIGHT,c2LIGHT,inv_c2LIGHT,boris_cLIGHT_factor
  
  real, parameter :: pi     = 3.141592653589793    ! Pi
  real, parameter :: CON_k  = 1.3807E-23           ! Boltzmann constant [J/K]
  real, parameter :: CON_mp = 1.6726E-27           ! Proton Mass [kg]
  real, parameter :: CON_mu = pi*4.0E-7            ! Permeability of free space [H/m]
  real, parameter :: CON_e  = 1.6022E-19           ! fundamental charge [C (coulomb)]

  !\
  ! Dipole and multipole expansion terms.
  !/
  real :: Bdp,Bdpx,Bdpy,Bdpz,Bdp_dim      ! constants for the dipole moment of B0
  real, dimension(1:3,1:3) :: Qqp         ! parameters for the quadrupole moment of B0
  real, dimension(1:3,1:3,1:3) :: Oop     ! parameters for the octupole moment of B0
  real :: THETAtilt, &                    ! tilt angle of magnetic axis
          sinTHETAtilt,cosTHETAtilt  
  real :: Bdp2, Bdp2_dim

  !\
  ! The following are some notes on how to pick the Q's.  I have used the
  ! cartesian version of the quadrupole magnetic potential because it was
  ! the easiest to differentiate in our coordinate system.  The two ways 
  ! to write the potential are as follows (in SI - see Jackson pp.136-8):
  !
  !    spherical:       phi = (1/5)*mu * sum(m=-2..2) qlm * Ylm/r^3
  !    cartesian:       phi = 1/(8*pi)*mu * sum(i,j=1..3) Qqpij * xi*xj/r^5
  !
  !  the coefficients are related to each other as follows (note that I am 
  !  using the relations in Jackson.  He uses Gaussian units.  I would guess
  !  that these relations are the same for both systems but I have not checked
  !  them.  They are most usefull in getting the theta and phi dependance that 
  !  you want and are not really used to do any type of converting):
  !
  !               q22 = (1/12)*sqrt(15/2/pi)*(Qqp11-i*Qqp12 - Qqp22)
  !               q21 = -(1/3)*sqrt(15/8/pi)*(Qqp13-i*Qqp23)
  !               q20 = (1/2)*sqrt(5/4/pi)*Qqp33
  !
  !  Note that Qqp is TRACELESS.  Also note that it is symmetric.  The q's have
  !  the following property:
  !
  !                  ql(-m) = (-1)^m  *  complex_conjugate(qlm)   
  !
  !/

  !\
  ! Heliosphere terms.
  !/
  real ::  RGuni, Avogadro, Cgrav, MWproton, Mp, RGp
  real ::  Tsun, Presun, NDsun, RHOsun, SSPsun, Tsunrot
  real ::  Msun, Rsun, Gsun, Velsun, OMEGAsun
  real ::  Qsun, Theat, Rheat, SIGMAheat
  real ::  Xplanet, Yplanet, Zplanet, Rplanet, RADIUSplanet, &
           Tplanet, OMEGAplanet, THETAplanet, PHIplanet, &
           PHIplanet_offset_orbit, PHIplanet_offset_rotframe, &
           NXplanet, NYplanet, NZplanet
  real ::  Xspacecraft, Yspacecraft, Zspacecraft, Rspacecraft
  integer :: iPEplanet, iPEspacecraft
 
  !\
  ! Far field solar wind solution variables.
  !/
  real ::      SW_T_dim  , &
               SW_a_dim  , &
       SW_rho, SW_rho_dim, &
       SW_p  , SW_p_dim  , &
       SW_Ux , SW_Ux_dim , &
       SW_Uy , SW_Uy_dim , &
       SW_Uz , SW_Uz_dim , &
       SW_Bx , SW_Bx_dim , &
       SW_By , SW_By_dim , &
       SW_Bz , SW_Bz_dim , &
       SW_B_factor

  real, dimension(0:1) :: &
       SW_rho_t,  &
       SW_p_t  ,  &
       SW_Ux_t ,  &
       SW_Uy_t ,  &
       SW_Uz_t ,  &
       SW_Bx_t ,  &
       SW_By_t ,  &
       SW_Bz_t ,  &
       SW_time_t
  
  !\
  ! Earth magnetosphere terms + generic magnetosphere terms.
  !/
  real ::  SSPearth
  real ::  Rearth, Rbody, Rcurrents, Mearth
  real ::  Body_rho_dim, Body_T_dim, Body_rho, Body_p
  real ::  Tearthrot, OMEGAearth, OMEGAbody, Gearth, Gbody

  !\
  ! Saturn/Jupiter magnetosphere and neutral torus terms.
  !/
  real ::  SSPsaturn, SSPjupiter
  real ::  Rsaturn, Rjupiter, Msaturn, Mjupiter
  real ::  Tsaturnrot, Tjupiterrot, OMEGAsaturn, OMEGAjupiter, &
           Gsaturn, Gjupiter
  real ::  Rio, Rio_torus, Qio, OMEGAio
  real ::  Rtitan, Rtitan_orbit, OMEGAtitan_orbit

  !\
  ! Magnetic Diffusion test case.
  !/
  real ::  SSPmagdiff, Rmagdiff

  !\
  ! General variables for the second body
  !/
  real :: Rbody2, Rcurrents2, xBody2, yBody2, zBody2

  !\
  ! Venus ionosphere terms.
  !/
  real ::  SSPvenus
  real ::  Rvenus, Mvenus
  real ::  Tvenusrot, OMEGAvenus, Gvenus
  
end module ModPhysics
