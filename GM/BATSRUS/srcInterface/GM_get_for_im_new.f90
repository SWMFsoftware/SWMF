!===========================================================================

subroutine GM_get_for_im_new

  use ModProcMH,   ONLY: iProc
  use CON_coupler, ONLY: Grid_C, IM_, GM_
  use ModPhysics,  ONLY: rBody
  use ModNumConst, ONLY: cDegToRad
  use CON_physics, ONLY: map_planet_field, transform_matrix
  use CON_time,    ONLY: tSimulation
  use ModCoordTransform, ONLY: sph_to_xyz

  implicit none
  integer :: nLat, nLon, iLat, iLon, iProcFound, iBlockFound, i, j, k
  real    :: r,Theta, Phi, XyzIono_D(3), Xyz_D(3), GmIm_DD(3,3)
  integer :: iHemisphere
  !---------------------------------------------------------------------------

  ! Transformation matrix from IM to GM
  GmIm_DD = transform_matrix(tSimulation,&
       Grid_C(IM_) % TypeCoord, Grid_C(GM_) % TypeCoord)

  ! Grid size of IM grid
  nLat = Grid_C(IM_) % nCoord_D(1)
  nLon = Grid_C(IM_) % nCoord_D(2)

  ! Integrate rays starting from each IM grid point
  r = 1.0
  do iLat = 1,nLat
     Theta =  Grid_C(IM_) % Coord1_I(iLat)
     do iLon = 1,nLon
        Phi = Grid_C(IM_) % Coord2_I(iLon)

        ! Convert to SMG Cartesian coordinates
        call sph_to_xyz(r,Theta,Phi,XyzIono_D)

        ! Map from the ionosphere to rBody
        call map_planet_field(tSimulation, XyzIono_D, 'SMG NORM', &
             rBody, Xyz_D, iHemisphere)

        ! Check if the mapping is on the north hemisphere
        if(iHemisphere /= 1) CYCLE

        ! Convert to the GM coordinate system
        Xyz_D = matmul(GmIm_DD,Xyz_D)

        ! Check if the point is on this processor
        call xyz_to_peblk(Xyz_D(1),Xyz_D(2),Xyz_D(3),iProcFound,iBlockFound,&
             .false.,i,j,k)

        if(iProcFound /= iProc) CYCLE

        ! Put the initial volumes into the arrays
        !! !!! to be implemented

        ! Follow ray towards south, which is opposite of the field line
        ! iRay = 2
        call follow_ray(2, iLat, iLon, 0, 0)
     end do
  end do

end subroutine GM_get_for_im_new

