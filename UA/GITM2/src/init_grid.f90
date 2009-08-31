
subroutine init_grid

  use ModGITM
  use ModInputs
  use ModConstants
  use ModPlanet
  use ModSphereInterface
  use ModTime
  use ModEUV,     only: init_mod_euv
  use ModSources, only: init_mod_sources
  implicit none

  type (UAM_ITER) :: r_iter

  integer :: iBlock

  logical :: IsOk, IsDone, DoTouchSouth, DoTouchNorth

  call report("init_grid",1)

  if (.not. Is1D) then

     if (IsFullSphere) then
        call UAM_module_setup(iCommGITM, &
             nLons, nLats, nAlts, &
             nBlocksLon, nBlocksLat, &
             -pi/2.0, pi/2.0, &
             .true., .true., &
             0.0, 0.0, 0.0, &
             RBody+AltMin, 5000.0, &
             ok=IsOk)
     else

        DoTouchNorth = .false.
        DoTouchSouth = .false.

        if (LatEnd >= pi/2) then
           LatEnd = pi/2
           DoTouchNorth = .true.
        endif
        if (LatStart <= -pi/2) then
           LatStart = -pi/2
           DoTouchSouth = .true.
        endif

        call UAM_module_setup(iCommGITM, &
             nLons, nLats, nAlts, &
             nBlocksLon, nBlocksLat, &
             LatStart, LatEnd, &
             DoTouchSouth, DoTouchNorth, &
             0.0, 0.0, 0.0, &
             RBody+AltMin, 5000.0, &
             ok=IsOk)
     endif

     if (.not.IsOk) call stop_gitm("Error in trying to create grid.")

     call UAM_XFER_create(ok=IsOk)
     if (.not. IsOk) then
        call UAM_write_error()
        call stop_gitm("Error with UAM_XFER_create")
     endif

     call UAM_ITER_create(r_iter)
     call UAM_ITER_reset(r_iter,iBlock,IsDone)

     nBlocks = 0
     do while (.not. IsDone)
        nBlocks = nBlocks + 1
        call UAM_ITER_next(r_iter,iBlock,IsDone)
     enddo

  else
     nBlocks = 1
     Latitude = LatStart
     Longitude = LonStart
  endif

  call init_mod_gitm
  call init_mod_euv
  call init_mod_sources

end subroutine init_grid
