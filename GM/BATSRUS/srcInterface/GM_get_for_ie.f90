!^CFG COPYRIGHT UM
!^CMP FILE IE
!============================================
!                                           |
!     Magnetosphere/Ionosphere Coupling     |
!     Subroutines here are used with MPI    |
!============================================
Module ModFieldAlignedCurrent

  Use ModMappingParam

  implicit none
  save
  integer :: IONO_NORTH_nMagBndPts=-1, IONO_SOUTH_nMagBndPts=-1
  real, dimension(:), allocatable ::                          &
       MAG_NORTH_X,MAG_NORTH_Y,MAG_NORTH_Z,MAG_NORTH_R,       & !Magnetospheric coordinates
       MAG_NORTH_Theta,MAG_NORTH_Psi,                         & !
       MAG_SOUTH_X,MAG_SOUTH_Y,MAG_SOUTH_Z,MAG_SOUTH_R,       & !
       MAG_SOUTH_Theta,MAG_SOUTH_Psi,                         & !
       MAG_NORTH_JR,MAG_NORTH_Jx,MAG_NORTH_Jy,MAG_NORTH_Jz,   & !Magnetospheric current
       MAG_SOUTH_JR,MAG_SOUTH_Jx,MAG_SOUTH_Jy,MAG_SOUTH_Jz

  real, dimension(:,:), allocatable ::   &
       MAG_NORTH_MagField, MAG_SOUTH_MagField

  real, dimension(:,:), allocatable ::   &
       MAG_NORTH_IONO_LOC, MAG_SOUTH_IONO_LOC

  !\
  ! Magnetosphere inner boundary current solution variable definitions.
  !/
  integer :: nMagBndPts_North, nMagBndPts_South
  integer, dimension(:), allocatable :: nMagBndPts_North_PE, nMagBndPts_South_PE
  real, dimension(:), allocatable ::             &
       Xmag_North,Ymag_North,Zmag_North,         & ! Magnetospheric coordinates
       Xmag_South,Ymag_South,Zmag_South,         & !
       JXmag_North,JYmag_North,JZmag_North,      & ! Magnetospheric current
       JXmag_South,JYmag_South,JZmag_South

end Module ModFieldAlignedCurrent
!==========================================================================
subroutine GM_get_for_ie(Buffer_IV,nPoint,nVar,NameVar)

  use ModFieldAlignedCurrent,ONLY:&
       MAG_NORTH_IONO_LOC, MAG_SOUTH_IONO_LOC,&
       MAG_NORTH_Jx,MAG_NORTH_Jy,MAG_NORTH_Jz,MAG_NORTH_MagField,&
       MAG_SOUTH_Jx,MAG_SOUTH_Jy,MAG_SOUTH_Jz,MAG_SOUTH_MagField
  implicit none

  character (len=*), parameter :: NameSub='GM_get_for_ie'

  integer, intent(in) :: nPoint,nVar
  real, intent(out), dimension(nPoint,nVar) :: Buffer_IV
  character (len=*), intent(in) :: NameVar

  logical :: DoTest, DoTestMe
  !--------------------------------------------------------------------------
  call CON_set_do_test(NameSub,DoTest, DoTestMe)
  if(DoTest)write(*,*)NameSub,': starting with NameVar=',NameVar
  select case(NameVar)
  case('LocNorth3')
     Buffer_IV = MAG_NORTH_IONO_LOC
  case('LocSouth3')
     Buffer_IV = MAG_SOUTH_IONO_LOC
  case('jNorth3:BinfoNorth5')
     Buffer_IV(:,1)  =MAG_NORTH_Jx
     Buffer_IV(:,2)  =MAG_NORTH_Jy
     Buffer_IV(:,3)  =MAG_NORTH_Jz
     Buffer_IV(:,4:8)=MAG_NORTH_MagField
  case('jSouth3:BinfoSouth5')
     Buffer_IV(:,1)  =MAG_SOUTH_Jx
     Buffer_IV(:,2)  =MAG_SOUTH_Jy
     Buffer_IV(:,3)  =MAG_SOUTH_Jz
     Buffer_IV(:,4:8)=MAG_SOUTH_MagField
  case default
     call CON_stop(NameSub//' invalid NameVar='//NameVar)
  end select
  if(DoTest)write(*,*)NameSub,': finished with NameVar=',NameVar

end subroutine GM_get_for_ie


subroutine magnetosphere_allocate

  use ModProcMH
  use ModMain
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK,R_BLK,Rmin_BLK
  use ModPhysics,ONLY:RCurrents
  use ModFieldAlignedCurrent
  use ModMpi

  implicit none

  integer :: i, j, k, nn, ns, iBLK, iPE, iError, &
       itag = 99, status(MPI_STATUS_SIZE)
  real, dimension(3) :: Mag_Loc, Iono_Loc
  real :: iono_r

  call set_mapping_param

  allocate(nMagBndPts_North_PE(nProc), stat = iError)
  if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
       " allocation error for nMagBndPts_North_PE", &
       & " nProc = ",nProc
  nMagBndPts_North_PE = 0

  allocate(nMagBndPts_South_PE(nProc), stat = iError)
  if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
       " allocation error for nMagBndPts_South_PE", &
       & " nProc = ",nProc
  nMagBndPts_South_PE = 0

  nMagBndPts_North = 0
  nMagBndPts_South = 0
  do iBLK = 1,nBlockMax
     if (.not.unusedBLK(iBLK) .and. Rmin_BLK(iBLK) < Rcurrents) then
        do k = 1, nK
           do j = 1, nJ
              do i = 1, nI

                 if (R_BLK(i,j,k,iBLK) > Rcurrents .and. &
                      ((R_BLK(i+1,j,k,iBLK) <= Rcurrents) .or. &
                      (R_BLK(i-1,j,k,iBLK) <= Rcurrents) .or. &
                      (R_BLK(i,j+1,k,iBLK) <= Rcurrents) .or. &
                      (R_BLK(i,j-1,k,iBLK) <= Rcurrents) .or. &
                      (R_BLK(i,j,k+1,iBLK) <= Rcurrents) .or. &
                      (R_BLK(i,j,k-1,iBLK) <= Rcurrents)) ) then

                    ! Figure out whether the point maps to the
                    !  northern hemisphere or the southern hemisphere
                    Mag_Loc(1) = x_BLK(i,j,k,iBLK)
                    Mag_Loc(2) = y_BLK(i,j,k,iBLK)
                    Mag_Loc(3) = z_BLK(i,j,k,iBLK)
                    call Get_Mapping_Point(Mag_Loc, 1.0, Iono_Loc)
                    iono_r = sum(Iono_Loc**2)

                    if (iono_r >= 0.9) then
                       nMagBndPts_North = nMagBndPts_North + 1 
                    else
                       nMagBndPts_South = nMagBndPts_South + 1 
                    end if
                 end if

              end do !end i loop
           end do !end j loop
        end do !end k loop

     end if
  end do

  if (nMagBndPts_North > 0) then
     allocate(Xmag_North(nMagBndPts_North), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for Xmag_North", &
          & " nMagBndPts_North = ",nMagBndPts_North 
     Xmag_North = 0.00
     allocate(Ymag_North(nMagBndPts_North), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for Ymag_North", &
          & " nMagBndPts_North = ",nMagBndPts_North 
     Ymag_North = 0.00
     allocate(Zmag_North(nMagBndPts_North), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for Zmag_North", &
          & " nMagBndPts_North = ",nMagBndPts_North
     Zmag_North = 0.00
     allocate(JXmag_North(nMagBndPts_North), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for JXmag_North", &
          & " nMagBndPts_North = ",nMagBndPts_North
     JXmag_North = 0.00
     allocate(JYmag_North(nMagBndPts_North), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for JYmag_North", &
          & " nMagBndPts_North = ",nMagBndPts_North
     JYmag_North = 0.00
     allocate(JZmag_North(nMagBndPts_North), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for JZmag_North", &
          & " nMagBndPts_North = ",nMagBndPts_North
     JZmag_North = 0.00

  end if

  if (nMagBndPts_South > 0) then
     allocate(Xmag_South(nMagBndPts_South), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for Xmag_South", &
          & " nMagBndPts_South = ",nMagBndPts_South
     Xmag_South = 0.00
     allocate(Ymag_South(nMagBndPts_South), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for Ymag_South", &
          & " nMagBndPts_South = ",nMagBndPts_South
     Ymag_South = 0.00
     allocate(Zmag_South(nMagBndPts_South), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for Zmag_South", &
          & " nMagBndPts_South = ",nMagBndPts_South
     Zmag_South = 0.00
     allocate(JXmag_South(nMagBndPts_South), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for JXmag_South", &
          & " nMagBndPts_South = ",nMagBndPts_South
     JXmag_South = 0.00
     allocate(JYmag_South(nMagBndPts_South), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for JYmag_South", &
          & " nMagBndPts_South = ",nMagBndPts_South
     JYmag_South = 0.00
     allocate(JZmag_South(nMagBndPts_South), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for JZmag_South", &
          & " nMagBndPts_South = ",nMagBndPts_South
     JZmag_South = 0.00

  end if

  nn = 0
  ns = 0
  do iBLK = 1,nBlockMax
     if (.not.unusedBLK(iBLK) .and. Rmin_BLK(iBLK) < Rcurrents) then
        do k = 1, nK
           do j = 1, nJ
              do i = 1, nI

                 if (R_BLK(i,j,k,iBLK) > Rcurrents .and. &
                      ((R_BLK(i+1,j,k,iBLK) <= Rcurrents) .or. &
                      (R_BLK(i-1,j,k,iBLK) <= Rcurrents) .or. &
                      (R_BLK(i,j+1,k,iBLK) <= Rcurrents) .or. &
                      (R_BLK(i,j-1,k,iBLK) <= Rcurrents) .or. &
                      (R_BLK(i,j,k+1,iBLK) <= Rcurrents) .or. &
                      (R_BLK(i,j,k-1,iBLK) <= Rcurrents)) ) then

                    ! Figure out whether the point maps to the
                    !  northern hemisphere or the southern hemisphere
                    Mag_Loc(1) = x_BLK(i,j,k,iBLK)
                    Mag_Loc(2) = y_BLK(i,j,k,iBLK)
                    Mag_Loc(3) = z_BLK(i,j,k,iBLK)
                    call Get_Mapping_Point(Mag_Loc, 1.0, Iono_Loc)
                    iono_r = sum(Iono_Loc**2)

                    if (iono_r >= 0.9) then

                       nn = nn + 1
                       Xmag_North(nn) = x_BLK(i,j,k,iBLK)
                       Ymag_North(nn) = y_BLK(i,j,k,iBLK)
                       Zmag_North(nn) = z_BLK(i,j,k,iBLK)
                    else
                       ns = ns + 1
                       Xmag_South(ns) = x_BLK(i,j,k,iBLK)
                       Ymag_South(ns) = y_BLK(i,j,k,iBLK)
                       Zmag_South(ns) = z_BLK(i,j,k,iBLK)
                    end if
                 end if

              end do !end i loop
           end do !end j loop
        end do !end k loop

     end if
  end do

  call MPI_ALLGATHER(nMagBndPts_North, 1, MPI_INTEGER, &
       nMagBndPts_North_PE(1), 1, MPI_INTEGER, iComm, iError)

  call MPI_ALLGATHER(nMagBndPts_South, 1, MPI_INTEGER, &
       nMagBndPts_South_PE(1), 1, MPI_INTEGER, iComm, iError)

  if (iProc == 0) then
     IONO_NORTH_nMagBndPts = 0
     do iPE = 1, nProc 
        IONO_NORTH_nMagBndPts = IONO_NORTH_nMagBndPts + nMagBndPts_North_PE(iPE)
     end do

     IONO_South_nMagBndPts = 0
     do iPE = 1, nProc 
        IONO_South_nMagBndPts = IONO_South_nMagBndPts + nMagBndPts_South_PE(iPE)
     end do

     !
     ! North
     !

     allocate(MAG_NORTH_X(IONO_NORTH_nMagBndPts), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for MAG_NORTH_X"
     MAG_NORTH_X = 0.00
     allocate(MAG_NORTH_Y(IONO_NORTH_nMagBndPts), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for MAG_NORTH_Y"
     MAG_NORTH_Y = 0.00
     allocate(MAG_NORTH_Z(IONO_NORTH_nMagBndPts), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for MAG_NORTH_Z"
     MAG_NORTH_Z = 0.00
     allocate(MAG_NORTH_Theta(IONO_NORTH_nMagBndPts), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for MAG_NORTH_Theta"
     MAG_NORTH_Theta = 0.00
     allocate(MAG_NORTH_Psi(IONO_NORTH_nMagBndPts), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for MAG_NORTH_Psi"
     MAG_NORTH_Psi = 0.00
     allocate(MAG_NORTH_Jx(IONO_NORTH_nMagBndPts), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for MAG_NORTH_Jx"
     MAG_NORTH_Jx = 0.00
     allocate(MAG_NORTH_Jy(IONO_NORTH_nMagBndPts), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for MAG_NORTH_Jy"
     MAG_NORTH_Jy = 0.00
     allocate(MAG_NORTH_Jz(IONO_NORTH_nMagBndPts), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for MAG_NORTH_Jz"
     MAG_NORTH_Jz = 0.00

     allocate(MAG_NORTH_JR(IONO_NORTH_nMagBndPts), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for MAG_NORTH_JR"
     MAG_NORTH_JR = 0.00

     !
     ! The 5 points in this array are:
     !   bx
     !   by
     !   bz
     !   |Biono|/|Bmagneto|
     !   cos(btilt in ionosphere)
     !

     allocate(MAG_NORTH_MagField(IONO_NORTH_nMagBndPts,5), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for MAG_NORTH_MagField"
     MAG_NORTH_MagField = 0.00

     !
     ! The 3 points in this array are:
     !   X in ionosphere
     !   Y in ionosphere
     !   Z in ionosphere
     !

     allocate(MAG_NORTH_IONO_LOC(IONO_NORTH_nMagBndPts,3), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for MAG_NORTH_IONO_LOC"
     MAG_NORTH_IONO_LOC = 0.00


     !
     ! South
     !

     allocate(MAG_SOUTH_X(IONO_SOUTH_nMagBndPts), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for MAG_SOUTH_X"
     MAG_SOUTH_X = 0.00
     allocate(MAG_SOUTH_Y(IONO_SOUTH_nMagBndPts), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for MAG_SOUTH_Y"
     MAG_SOUTH_Y = 0.00
     allocate(MAG_SOUTH_Z(IONO_SOUTH_nMagBndPts), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for MAG_SOUTH_Z"
     MAG_SOUTH_Z = 0.00
     allocate(MAG_SOUTH_Theta(IONO_SOUTH_nMagBndPts), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for MAG_SOUTH_Theta"
     MAG_SOUTH_Theta = 0.00
     allocate(MAG_SOUTH_Psi(IONO_SOUTH_nMagBndPts), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for MAG_SOUTH_Psi"
     MAG_SOUTH_Psi = 0.00
     allocate(MAG_SOUTH_Jx(IONO_SOUTH_nMagBndPts), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for MAG_SOUTH_Jx"
     MAG_SOUTH_Jx = 0.00
     allocate(MAG_SOUTH_Jy(IONO_SOUTH_nMagBndPts), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for MAG_SOUTH_Jy"
     MAG_SOUTH_Jy = 0.00
     allocate(MAG_SOUTH_Jz(IONO_SOUTH_nMagBndPts), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for MAG_SOUTH_Jz"
     MAG_SOUTH_Jz = 0.00
     allocate(MAG_SOUTH_JR(IONO_SOUTH_nMagBndPts), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for MAG_SOUTH_JR"
     MAG_SOUTH_JR = 0.00

     allocate(MAG_SOUTH_MagField(IONO_SOUTH_nMagBndPts,5), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for MAG_SOUTH_MagField"
     MAG_SOUTH_MagField = 0.00

     allocate(MAG_SOUTH_IONO_LOC(IONO_SOUTH_nMagBndPts,3), stat = iError)
     if (iError > 0) write(*,*) "magnetosphere_fac: PE = ",iProc, &
          " allocation error for MAG_SOUTH_IONO_LOC"
     MAG_SOUTH_IONO_LOC = 0.00

     nn = 0
     if (nMagBndPts_North_PE(1) > 0) then
        MAG_NORTH_X(nn+1:nn+nMagBndPts_North_PE(1)) = Xmag_North
        MAG_NORTH_Y(nn+1:nn+nMagBndPts_North_PE(1)) = Ymag_North
        MAG_NORTH_Z(nn+1:nn+nMagBndPts_North_PE(1)) = Zmag_North
        nn = nn + nMagBndPts_North_PE(1)
     end if
     do iPE = 2, nProc 
        if (nMagBndPts_North_PE(iPE) > 0) then
           call MPI_RECV(MAG_NORTH_X(nn+1), nMagBndPts_North_PE(iPE), &
                MPI_REAL, iPE-1, itag, iComm, status, iError)
           call MPI_RECV(MAG_NORTH_Y(nn+1), nMagBndPts_North_PE(iPE), &
                MPI_REAL, iPE-1, itag, iComm, status, iError)
           call MPI_RECV(MAG_NORTH_Z(nn+1), nMagBndPts_North_PE(iPE), &
                MPI_REAL, iPE-1, itag, iComm, status, iError)
           nn = nn + nMagBndPts_North_PE(iPE)
        end if
     end do

     ns = 0
     if (nMagBndPts_South_PE(1) > 0) then
        MAG_SOUTH_X(ns+1:ns+nMagBndPts_South_PE(1)) = Xmag_South
        MAG_SOUTH_Y(ns+1:ns+nMagBndPts_South_PE(1)) = Ymag_South
        MAG_SOUTH_Z(ns+1:ns+nMagBndPts_South_PE(1)) = Zmag_South
        ns = ns + nMagBndPts_South_PE(1)
     end if
     do iPE = 2, nProc 
        if (nMagBndPts_South_PE(iPE) > 0) then
           call MPI_RECV(MAG_SOUTH_X(ns+1), nMagBndPts_South_PE(iPE), &
                MPI_REAL, iPE-1, itag, iComm, status, iError)
           call MPI_RECV(MAG_SOUTH_Y(ns+1), nMagBndPts_South_PE(iPE), &
                MPI_REAL, iPE-1, itag, iComm, status, iError)
           call MPI_RECV(MAG_SOUTH_Z(ns+1), nMagBndPts_South_PE(iPE), &
                MPI_REAL, iPE-1, itag, iComm, status, iError)
           ns = ns + nMagBndPts_South_PE(iPE)
        end if
     end do
  else
     if (nMagBndPts_North > 0) then
        call MPI_SEND(Xmag_North(1), nMagBndPts_North, &
             MPI_REAL, 0, itag, iComm, iError)
        call MPI_SEND(Ymag_North(1), nMagBndPts_North, &
             MPI_REAL, 0, itag, iComm, iError)
        call MPI_SEND(Zmag_North(1), nMagBndPts_North, &
             MPI_REAL, 0, itag, iComm, iError)
     end if

     if (nMagBndPts_South > 0) then
        call MPI_SEND(Xmag_South(1), nMagBndPts_South, &
             MPI_REAL, 0, itag, iComm, iError)
        call MPI_SEND(Ymag_South(1), nMagBndPts_South, &
             MPI_REAL, 0, itag, iComm, iError)
        call MPI_SEND(Zmag_South(1), nMagBndPts_South, &
             MPI_REAL, 0, itag, iComm, iError)
     end if
  end if

end subroutine magnetosphere_allocate


subroutine magnetosphere_deallocate
  use ModProcMH
  use ModFieldAlignedCurrent
  implicit none

  if (.not.allocated(nMagBndPts_North_PE)) return

  deallocate(nMagBndPts_North_PE)
  deallocate(nMagBndPts_South_PE)

  if (nMagBndPts_North > 0) then
     deallocate(Xmag_North)
     deallocate(Ymag_North)
     deallocate(Zmag_North)
     deallocate(JXmag_North)
     deallocate(JYmag_North)
     deallocate(JZmag_North)
  end if

  if (nMagBndPts_South > 0) then
     deallocate(Xmag_South)
     deallocate(Ymag_South)
     deallocate(Zmag_South)
     deallocate(JXmag_South)
     deallocate(JYmag_South)
     deallocate(JZmag_South)
  end if

  if (iProc == 0) then

     deallocate(MAG_NORTH_X)
     deallocate(MAG_NORTH_Y)
     deallocate(MAG_NORTH_Z)
     deallocate(MAG_NORTH_Theta)
     deallocate(MAG_NORTH_Psi)
     deallocate(MAG_NORTH_Jx)
     deallocate(MAG_NORTH_Jy)
     deallocate(MAG_NORTH_Jz)
     deallocate(MAG_NORTH_JR)
     deallocate(MAG_NORTH_MagField)
     deallocate(MAG_NORTH_IONO_LOC)

     deallocate(MAG_SOUTH_X)
     deallocate(MAG_SOUTH_Y)
     deallocate(MAG_SOUTH_Z)
     deallocate(MAG_SOUTH_Theta)
     deallocate(MAG_SOUTH_Psi)
     deallocate(MAG_SOUTH_Jx)
     deallocate(MAG_SOUTH_Jy)
     deallocate(MAG_SOUTH_Jz)
     deallocate(MAG_SOUTH_JR)
     deallocate(MAG_SOUTH_MagField)
     deallocate(MAG_SOUTH_IONO_LOC)

  end if

end subroutine magnetosphere_deallocate


!==============================================================================
subroutine GM_calc_fac(IsNewFacPoint, nPoint_I, ColatLim_B)

  !\
  ! Calculate the incoming and outgoing field-aligned currents
  ! (FAC) to the ionosphere at the mapping points in the magnetosphere
  ! defined at a radius of Rcurrents and send this solution 
  ! information to the PE that is performing the ionosphere
  ! potential solution calculation.
  !/

  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,nBlockMax,unusedBLK,Time_Simulation,lVerbose, &
       iNewGrid,iNewDecomposition
  use ModAdvance, ONLY : State_VGB, Bx_, By_, Bz_
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK,R_BLK,Rmin_BLK,dx_BLK,dy_BLK,dz_BLK
  use ModPhysics, ONLY : Rbody,Rcurrents, unitSI_J
  use ModFieldAlignedCurrent
  use ModMpi
  use ModIo, ONLY: write_prefix, iUnitOut
  use CON_physics, ONLY: get_axes

  implicit none

  ! Logical for new positions and 
  ! the number of mapped FAC points for North and South hemispheres
  logical, intent(out) :: IsNewFacPoint
  integer, intent(out) :: nPoint_I(2)
  real, intent(out)    :: ColatLim_B(2)

  integer :: i, j, k, nn, ns, iBLK, iPE, iError, &
       itag = 99, status(MPI_STATUS_SIZE)
  real, dimension(3) :: Mag_Loc, Iono_Loc, CurrentDensity_D
  real :: iono_r

  integer :: iLastGrid = -1, iLastDecomposition = -1
  real    :: LastMagAxisTilt = -999.0
  real    :: MagAxisTilt

  character(len=*),parameter :: NameSub='GM_calc_fac'
  logical :: DoTest, DoTestMe
  !--------------------------------------------------------------------------

  call CON_set_do_test(NameSub,DoTest,DoTestMe)
  if(DoTestMe)write(*,*)NameSub,': starting'

  call get_axes(Time_Simulation, MagAxisTiltGsmOut = MagAxisTilt)

  IsNewFacPoint = iNewGrid/=iLastGrid .or.         &
       iNewDecomposition/=iLastDecomposition .or.  &
       abs(MagAxisTilt-LastMagAxisTilt) > 0.00001

  if(IsNewFacPoint)then
     call magnetosphere_deallocate
     call magnetosphere_allocate
     if (iProc == 0) call GM_map_fac(ColatLim_B)
     if (iProc == 0 .and. (DoTest.or.lVerbose>0) &
          .and. (iNewGrid/=iLastGrid .or. MagAxisTilt/=LastMagAxisTilt)) then
        call write_prefix; write(iUnitOut,'(a55,i7)') &
             'Northern Magnetospheric FAC Mapping Points = ', &
             IONO_NORTH_nMagBndPts
        call write_prefix; write(iUnitOut,'(a55,i7)') &
             'Southern Magnetospheric FAC Mapping Points = ', &
             IONO_SOUTH_nMagBndPts
        call write_prefix; write(iUnitOut,*)
     end if
     iLastGrid = iNewGrid
     iLastDecomposition = iNewDecomposition
     LastMagAxisTilt = MagAxisTilt
  endif

!!! Added for SWMF coupling
  nPoint_I(1)=IONO_NORTH_nMagBndPts
  nPoint_I(2)=IONO_SOUTH_nMagBndPts

  IONO_Radius_Mag_Boundary = Rbody*IONO_Radius

  nn = 0
  ns = 0
  do iBLK = 1,nBlockMax
     if (.not.unusedBLK(iBLK) .and. Rmin_BLK(iBLK) < Rcurrents) then
        do k = 1, nK
           do j = 1, nJ
              do i = 1, nI

                 if (R_BLK(i,j,k,iBLK) > Rcurrents .and. &
                      ((R_BLK(i+1,j,k,iBLK) <= Rcurrents) .or. &
                      (R_BLK(i-1,j,k,iBLK) <= Rcurrents) .or. &
                      (R_BLK(i,j+1,k,iBLK) <= Rcurrents) .or. &
                      (R_BLK(i,j-1,k,iBLK) <= Rcurrents) .or. &
                      (R_BLK(i,j,k+1,iBLK) <= Rcurrents) .or. &
                      (R_BLK(i,j,k-1,iBLK) <= Rcurrents)) ) then

                    !^CFG IF CARTESIAN BEGIN
                    
                    CurrentDensity_D(1) = 0.50* &              
                            ((State_VGB(Bz_,i,j+1,k,iBLK) - &
                            State_VGB(Bz_,i,j-1,k,iBLK))/dy_BLK(iBLK) - &
                            (State_VGB(By_,i,j,k+1,iBLK) - &
                            State_VGB(By_,i,j,k-1,iBLK))/dz_BLK(iBLK))

                    CurrentDensity_D(2) =  0.50* &
                            ((State_VGB(Bx_,i,j,k+1,iBLK) - &
                            State_VGB(Bx_,i,j,k-1,iBLK))/dz_BLK(iBLK) - &
                            (State_VGB(Bz_,i+1,j,k,iBLK) - &
                            State_VGB(Bz_,i-1,j,k,iBLK))/dx_BLK(iBLK))

                    CurrentDensity_D(3) = 0.50* &
                            ((State_VGB(By_,i+1,j,k,iBLK) - &
                            State_VGB(By_,i-1,j,k,iBLK))/dx_BLK(iBLK) - &
                            (State_VGB(Bx_,i,j+1,k,iBLK) - &
                            State_VGB(Bx_,i,j-1,k,iBLK))/dy_BLK(iBLK))

                    !^CFG END CARTESIAN
!              call covariant_curlb(i,j,k,iBLK,CurrentDensity_D)!^CFG IF NOT CARTESIAN
                    
                    ! Figure out whether the point maps to the
                    !  northern hemisphere or the southern hemisphere
                    Mag_Loc(1) = x_BLK(i,j,k,iBLK)
                    Mag_Loc(2) = y_BLK(i,j,k,iBLK)
                    Mag_Loc(3) = z_BLK(i,j,k,iBLK)
                    call Get_Mapping_Point(Mag_Loc, 1.0, Iono_Loc)

                    iono_r = sum(Iono_Loc**2)

                    if (iono_r >= 0.9) then
                       nn = nn + 1
                       JXmag_North(nn) = CurrentDensity_D(1)
                       JYmag_North(nn) = CurrentDensity_D(2)
                       JZmag_North(nn) = CurrentDensity_D(3) 
                    else
                       ns = ns + 1
                       JXmag_South(ns) = CurrentDensity_D(1)
                       JYmag_South(ns) = CurrentDensity_D(2) 
                       JZmag_South(ns) = CurrentDensity_D(3)
                    end if
                    
                 end if

              end do !end i loop
           end do !end j loop
        end do !end k loop

     end if
  end do

  if (iProc == 0) then
     nn = 0
     if (nMagBndPts_North_PE(1) > 0) then
        MAG_NORTH_Jx(nn+1:nn+nMagBndPts_North_PE(1)) = JXmag_North
        MAG_NORTH_Jy(nn+1:nn+nMagBndPts_North_PE(1)) = JYmag_North
        MAG_NORTH_Jz(nn+1:nn+nMagBndPts_North_PE(1)) = JZmag_North
        nn = nn + nMagBndPts_North_PE(1)
     end if
     do iPE = 2, nProc 
        if (nMagBndPts_North_PE(iPE) > 0) then
           call MPI_RECV(MAG_NORTH_Jx(nn+1), nMagBndPts_North_PE(iPE), &
                MPI_REAL, iPE-1, itag, iComm, status, iError)
           call MPI_RECV(MAG_NORTH_Jy(nn+1), nMagBndPts_North_PE(iPE), &
                MPI_REAL, iPE-1, itag, iComm, status, iError)
           call MPI_RECV(MAG_NORTH_Jz(nn+1), nMagBndPts_North_PE(iPE), &
                MPI_REAL, iPE-1, itag, iComm, status, iError)
           nn = nn + nMagBndPts_North_PE(iPE)
        end if
     end do

     if(DoTestMe)write(*,*)NameSub,': sum(abs(North Jx,y,z))=',&
          sum(abs(MAG_NORTH_Jx)),sum(abs(MAG_NORTH_Jy)),sum(abs(MAG_NORTH_Jz))

!!!FRAMEWORK convert to SI units before IE would use this
     MAG_NORTH_Jx = MAG_NORTH_Jx * unitSI_J
     MAG_NORTH_Jy = MAG_NORTH_Jy * unitSI_J
     MAG_NORTH_Jz = MAG_NORTH_Jz * unitSI_J
     

     ns = 0
     if (nMagBndPts_South_PE(1) > 0) then
        MAG_SOUTH_Jx(ns+1:ns+nMagBndPts_South_PE(1)) = JXmag_South
        MAG_SOUTH_Jy(ns+1:ns+nMagBndPts_South_PE(1)) = JYmag_South
        MAG_SOUTH_Jz(ns+1:ns+nMagBndPts_South_PE(1)) = JZmag_South
        ns = ns + nMagBndPts_South_PE(1)
     end if
     do iPE = 2, nProc 
        if (nMagBndPts_South_PE(iPE) > 0) then
           call MPI_RECV(MAG_SOUTH_Jx(ns+1), nMagBndPts_South_PE(iPE), &
                MPI_REAL, iPE-1, itag, iComm, status, iError)
           call MPI_RECV(MAG_SOUTH_Jy(ns+1), nMagBndPts_South_PE(iPE), &
                MPI_REAL, iPE-1, itag, iComm, status, iError)
           call MPI_RECV(MAG_SOUTH_Jz(ns+1), nMagBndPts_South_PE(iPE), &
                MPI_REAL, iPE-1, itag, iComm, status, iError)
           ns = ns + nMagBndPts_South_PE(iPE)
        end if
     end do

     if(DoTestMe)write(*,*)NameSub,': sum(abs(South Jx,y,z))=',&
          sum(abs(MAG_SOUTH_Jx)),sum(abs(MAG_SOUTH_Jy)),sum(abs(MAG_SOUTH_Jz))

!!!FRAMEWORK convert to SI units before IE would use this
     MAG_SOUTH_Jx = MAG_SOUTH_Jx * unitSI_J
     MAG_SOUTH_Jy = MAG_SOUTH_Jy * unitSI_J
     MAG_SOUTH_Jz = MAG_SOUTH_Jz * unitSI_J
  else
     if (nMagBndPts_North > 0) then
        call MPI_SEND(JXmag_North(1), nMagBndPts_North, &
             MPI_REAL, 0, itag, iComm, iError)
        call MPI_SEND(JYmag_North(1), nMagBndPts_North, &
             MPI_REAL, 0, itag, iComm, iError)
        call MPI_SEND(JZmag_North(1), nMagBndPts_North, &
             MPI_REAL, 0, itag, iComm, iError)
     end if

     if (nMagBndPts_South > 0) then
        call MPI_SEND(JXmag_South(1), nMagBndPts_South, &
             MPI_REAL, 0, itag, iComm, iError)
        call MPI_SEND(JYmag_South(1), nMagBndPts_South, &
             MPI_REAL, 0, itag, iComm, iError)
        call MPI_SEND(JZmag_South(1), nMagBndPts_South, &
             MPI_REAL, 0, itag, iComm, iError)
     end if
  end if

end subroutine GM_calc_fac

!-------------------------------------------------------------------------

subroutine GM_map_fac(ColatLim_B)

  use ModFieldAlignedCurrent
  use ModPhysics
  use CON_physics, ONLY : get_planet_field
  use ModMain,     ONLY : Time_Simulation, TypeCoordSystem

  implicit none

!!! Added for SWMF
  real, intent(out) :: ColatLim_B(2)

  integer              :: i, n_save, n, ierror
  real, dimension(1:3) :: Mag_Loc, Iono_Loc, MagField_Orientation, &
                          MagField_Magneto, MagField_Iono

  real :: Mag_Strength_Iono, Mag_Strength_Magneto
  real :: Theta
  real :: dir

  logical :: done

  logical :: DoTestMe, DoTest
  !--------------------------------------------------------------------------
  call set_oktest('map_fac',DoTest,DoTestMe)

  if(DoTestMe)then
     write(*,*)'map_fac starting'
     write(*,*)'map_fac nFacNorth,nFacSouth=',&
          IONO_NORTH_nMagBndPts,IONO_SOUTH_nMagBndPts
     write(*,*)'map_fac sum(abs(MAG_NORTH_X,Y,Z)=',&
          sum(abs(MAG_NORTH_X)),sum(abs(MAG_NORTH_Y)),sum(abs(MAG_NORTH_Z))
     write(*,*)'map_fac sum(abs(MAG_SOUTH_X,Y,Z)=',&
          sum(abs(MAG_SOUTH_X)),sum(abs(MAG_SOUTH_Y)),sum(abs(MAG_SOUTH_Z))
  endif
  ! set the direction according to the direction of the
  ! intrisic dipole

  if (DipoleSign > 0.0) then
     dir = -1.0
  else
     dir =  1.0
  endif

  ! -------------------------------------------------------------------
  ! First figure out the magnetosphere down to ionosphere information
  ! -------------------------------------------------------------------

  ColatLim_B = 0.0


  do i = 1, IONO_NORTH_nMagBndPts

     Mag_Loc(1) = MAG_NORTH_X(i)
     Mag_Loc(2) = MAG_NORTH_Y(i)
     Mag_Loc(3) = MAG_NORTH_Z(i)

     call Get_Mapping_Point(Mag_Loc, 1.0, Iono_Loc)

     MAG_NORTH_IONO_LOC(i,1:3) = Iono_Loc

     call Get_MagneticField_Orient(Mag_Loc, MagField_Orientation)

     MAG_NORTH_MagField(i,1:3) = MagField_Orientation

     !
     ! We need to get the magnetic field strength and orientation
     ! in the magnetosphere and ionosphere.
     !
     
     ! Here is the magnetosphere:

     call get_planet_field(Time_Simulation,Mag_Loc,TypeCoordSystem//' NORM',&
          MagField_Magneto)
     Mag_Strength_Magneto = sqrt(sum(MagField_Magneto**2))

     ! Here is the ionosphere:

     call get_planet_field(Time_Simulation,Iono_Loc,'SMG NORM',&
          MagField_Iono)

     Mag_Strength_Iono = sqrt(sum(MagField_Iono**2))
     MagField_Orientation = MagField_Iono / Mag_Strength_Iono

     ! Get the ratio between the magnitude of the ionosphere B 
     ! and the magnetosphere B

     MAG_NORTH_MagField(i,4) = -dir * Mag_Strength_Iono / Mag_Strength_Magneto

     ! Get the absolute value of the cos of the angle between the magnetic 
     ! field and radial direction
     MAG_NORTH_MagField(i,5) = abs(sum(Iono_Loc*MagField_Orientation))/rMap

     ! Get the the colatitude limit
     Theta = acos(Iono_Loc(3)/rMap)
     if (Theta > ColatLim_B(1)) ColatLim_B(1) = Theta

  enddo

  !
  ! South
  !

  ColatLim_B(2) = IONO_PI

  do i = 1, IONO_SOUTH_nMagBndPts

     Mag_Loc(1) = MAG_SOUTH_X(i)
     Mag_Loc(2) = MAG_SOUTH_Y(i)
     Mag_Loc(3) = MAG_SOUTH_Z(i)

     call Get_Mapping_Point(Mag_Loc, -1.0, Iono_Loc)

     MAG_SOUTH_IONO_LOC(i,1:3) = Iono_Loc

     call Get_MagneticField_Orient(Mag_Loc, MagField_Orientation)

     MAG_SOUTH_MagField(i,1:3) = MagField_Orientation

     !
     ! We need to get the magnetic field strength and orientation
     ! in the magnetosphere and ionosphere.
     !
     
     ! Here is the magnetosphere:

     call get_planet_field(Time_Simulation,Mag_Loc,TypeCoordSystem//' NORM',&
          MagField_Magneto)
     Mag_Strength_Magneto = sqrt(sum(MagField_Magneto**2))

     ! Here is the ionosphere:

     call get_planet_field(Time_Simulation,Iono_Loc,'SMG NORM',&
          MagField_Iono)
     Mag_Strength_Iono = sqrt(sum(MagField_Iono**2))
     MagField_Orientation = MagField_Iono / Mag_Strength_Iono

     ! Get the ratio between the magnitude of the ionosphere B 
     ! and the magnetosphere B

     MAG_SOUTH_MagField(i,4) = dir*Mag_Strength_Iono / Mag_Strength_Magneto

     ! Get the cos of the angle between the magnetic field and radial

     MAG_SOUTH_MagField(i,5) = abs(sum(Iono_Loc*MagField_Orientation))/rMap

     Theta = acos(Iono_Loc(3)/rMap)
     if (Theta < ColatLim_B(2)) ColatLim_B(2) = Theta

  enddo

  if(DoTestMe)then
     write(*,*)'map_fac sum(abs(MAG_NORTH(:,1..5)))',&
          sum(abs(MAG_NORTH_MagField(:,1))),&
          sum(abs(MAG_NORTH_MagField(:,2))),&
          sum(abs(MAG_NORTH_MagField(:,3))),&
          sum(abs(MAG_NORTH_MagField(:,4))),&
          sum(abs(MAG_NORTH_MagField(:,5)))
     write(*,*)'map_fac sum(abs(MAG_SOUTH(:,1..5)))',&
          sum(abs(MAG_SOUTH_MagField(:,1))),&
          sum(abs(MAG_SOUTH_MagField(:,2))),&
          sum(abs(MAG_SOUTH_MagField(:,3))),&
          sum(abs(MAG_SOUTH_MagField(:,4))),&
          sum(abs(MAG_SOUTH_MagField(:,5)))

     write(*,*)'map_fac maxloc(NORTH_MagField(:,1)-SOUTH_MagField(:,1))=',&
          maxloc(abs(MAG_NORTH_MagField(:,1))-abs(MAG_SOUTH_MagField(:,1))),&
          maxval(abs(MAG_NORTH_MagField(:,1))-abs(MAG_SOUTH_MagField(:,1)))

     write(*,*)'map_fac MAG_NORTH(1,:)=', MAG_NORTH_MagField(1,:)
     write(*,*)'map_fac MAG_SOUTH(1,:)=', MAG_SOUTH_MagField(1,:)
     write(*,*)'map_fac finished'
  endif

end subroutine GM_map_fac
