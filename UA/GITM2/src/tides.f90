
subroutine read_tides

  use ModGitm
  use ModTime
  use ModInputs
  use ModTides

  implicit none

  integer :: iComp
  logical :: IsFirstTime = .true.
  real                           :: version 
  integer                        :: nAlts_gswm, nVars, iVar, iHour
  character(len=40), allocatable :: varName(:)
  integer                        :: time(7)
  real, allocatable              :: param(:,:,:,:,:)
  integer :: iError

  call report("read_tides",1)
  call start_timing("read_tides")

  if (IsFirstTime) then
     GSWM_file_name(1) ="UA/DataIn/diur_mig_99km_12.bin" 
     GSWM_file_name(2) ="UA/DataIn/diur_nonmig_99km_12.bin" 
     GSWM_file_name(3) ="UA/DataIn/semidiur_mig_99km_12.bin" 
     GSWM_file_name(4) ="UA/DataIn/semidiur_nonmig_99km_12.bin" 
     IsFirstTime = .false.
  endif

  iError = 0

  do iComp = 1,4

     if (UseGswmComp(iComp)) then

        write(*,*) "Reading File : ",gswm_file_name(iComp)
        open(iInputUnit_,file=gswm_file_name(iComp),&
             status="old",form="unformatted")
          
        read(iInputUnit_) version
        read(iInputUnit_) nLonsGswm, nLatsGswm, nAltsGswm, nvars
          
        if(.not.allocated(varname)) allocate(varname(nvars))
  
        do ivar = 1, nvars 
           read(iInputUnit_) varname(ivar) 
        enddo
          
        read(iInputUnit_) time

        if(.not.allocated(param)) &
             allocate(param(nLonsGswm, nLatsGswm, nAltsGswm, nvars, 24))
        if(.not.allocated(u_gswm)) then
           allocate(u_gswm(nLonsGswm, nLatsGswm, nAltsGswm, 24))
           u_gswm(:,:,:,:) = 0.
        endif
        
        if(.not.allocated(v_gswm)) then
           allocate(v_gswm(nLonsGswm, nLatsGswm, nAltsGswm, 24))
           v_gswm(:,:,:,:) = 0.
        endif

        if(.not.allocated(t_gswm)) then
           allocate(t_gswm(nLonsGswm, nLatsGswm, nAltsGswm, 24))
           t_gswm(:,:,:,:) = 0. 
        endif

        if(.not.allocated(lon_gswm)) allocate(lon_gswm(nLonsGswm))
        if(.not.allocated(lat_gswm)) allocate(lat_gswm(nLatsGswm))
        
        do ihour = 1, 24
           do ivar = 1, nvars 
              read(iInputUnit_) param(:,:,:,ivar,ihour) 
           enddo
        enddo

        lon_gswm(:)     = param(:,1,1,1,1)
        lat_gswm(:)     = param(1,:,1,2,1)
        ! The GSWM grid is uniform, so let's just assume this!
        dLatGswm = lat_gswm(2) - lat_gswm(1)
        dLonGswm = lon_gswm(2) - lon_gswm(1)
        t_gswm(:,:,:,:) = t_gswm(:,:,:,:) + param(:,:,:,4,:) 
        u_gswm(:,:,:,:) = u_gswm(:,:,:,:) + param(:,:,:,5,:) 
        v_gswm(:,:,:,:) = v_gswm(:,:,:,:) + param(:,:,:,6,:) 
        
        close(iInputUnit_)

        deallocate(varname)
        deallocate(param)

     endif

  enddo

  call end_timing("read_tides")

end subroutine read_tides

subroutine update_tides

  use ModGitm
  use ModTime, only: iTimeArray
  use ModInputs, only: iDebugLevel
  use ModTides

  implicit none

  integer :: iLonHalf, iBlock, iLon, iLat, iAlt
  integer :: iLonGSWM, iLatGSWM, iFac1, iFac2
  real :: rfac

  iLonHalf = pi / dLonGswm

  rfac = float(iTimeArray(5))/60.0 + float(iTimeArray(6))/3600.0

  iFac1 = iTimeArray(4)+1            ! this will go from 1-24
  iFac2 = mod(iTimeArray(4)+1,24)+1  ! this will go from 2-25 -> 2-24 back to 1

  do iBlock = 1, nBlocks
     do iAlt = 1,2
        do iLat=-1,nLats+2
           do iLon=-1,nLons+2
           
              iLonGSWM = floor(Longitude(iLon, iBlock) / dLonGSWM)
              iLonGSWM = mod(iLonGSWM + iLonHalf + nLonsGSWM, nLonsGswm) + 1
              iLatGSWM = floor((Latitude(iLat, iBlock) + Pi/2) / dLatGswm) + 1 

              iLatGSWM = max(1,iLatGSWM)
              iLatGSWM = min(nLatsGSWM,iLatGSWM)

              if (iDebugLevel > 5) then
                 write(*,*) "iLon, iLat : ", iLonGSWM, iLatGSWM
                 write(*,*) "lat/lon : ",Latitude(iLat, iBlock)*180/pi, &
                      Longitude(iLon, iBlock)*180/pi
                 write(*,*) "gswm lat/lon : ", &
                      lat_gswm(iLatGSWM)*180/pi, lon_gswm(iLonGSWM)*180/pi
              endif

              TidesTemp(iLon,iLat,iAlt,iBlock) = &
                   (1-rfac) * t_gswm(iLonGSWM, iLatGSWM, 1, iFac1) + &
                   (  rfac) * t_gswm(iLonGSWM, iLatGSWM, 1, iFac2)
              TidesEast(iLon,iLat,iAlt,iBlock) = &
                   (1-rfac) * u_gswm(iLonGSWM, iLatGSWM, 1, iFac1) + &
                   (  rfac) * u_gswm(iLonGSWM, iLatGSWM, 1, iFac2)
              TidesNorth(iLon,iLat,iAlt,iBlock) = &
                   (1-rfac) * v_gswm(iLonGSWM, iLatGSWM, 1, iFac1) + &
                   (  rfac) * v_gswm(iLonGSWM, iLatGSWM, 1, iFac2)

           enddo
        enddo
     enddo
  enddo

end subroutine update_tides
