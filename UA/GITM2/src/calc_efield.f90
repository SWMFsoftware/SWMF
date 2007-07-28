
subroutine calc_efield(iBlock)

  use ModGITM
  use ModInputs

  implicit none

  integer, intent(in) :: iBlock

  integer :: i, j, k, imax, jmax, kmax
  real :: maxi
  real :: p1, p2

  call report("Electric Field",2)
  call start_timing("calc_efield")

  EField = 0.0
  ExB    = 0.0

  maxi = 0.0

  !!! This is only first order accurate for stretched grids ???
  do k=1,nAlts
     do i=1,nLats
        do j=1,nLons
           EField(j,i,k,iEast_) = &
               -(Potential(j+1,i,k,iBlock)-Potential(j-1,i,k,iBlock)) &
                *0.5*InvdLonDist_GB(j,i,k,iBlock)

           EField(j,i,k,iNorth_)= &
                -(Potential(j,i+1,k,iBlock)-Potential(j,i-1,k,iBlock)) &
                *0.5*InvdLatDist_GB(j,i,k,iBlock)

           if (k < nAlts) then
              p1 = (potential(j,i,k+2,iBlock) + potential(j,i,k+1,iBlock) + &
                   potential(j,i,k,iBlock))/3.0
           else
              p1 = (potential(j,i,k+1,iBlock) + potential(j,i,k,iBlock))/2.0
           endif

           if (k > 1) then
              p2 = (potential(j,i,k-2,iBlock) + potential(j,i,k-1,iBlock) + &
                   potential(j,i,k,iBlock))/3.0
           else
              p2 = (potential(j,i,k-1,iBlock) + potential(j,i,k,iBlock))/2.0
           endif

           p1 = potential(j,i,k+1,iBlock)
           p2 = potential(j,i,k-1,iBlock)

           EField(j,i,k,iUp_) = -(P1-P2) / (2*dAlt_GB(j,i,k,iBlock))

        enddo
     enddo
  enddo

  ExB(:,:,:,iEast_)  =    EField(:,:,:,iNorth_) * B0(:,:,:,iUp_,iBlock)    - &
                          EField(:,:,:,iUp_)    * B0(:,:,:,iNorth_,iBlock)

  ExB(:,:,:,iNorth_) = - (EField(:,:,:,iEast_)  * B0(:,:,:,iUp_,iBlock)    - &
                        EField(:,:,:,iUp_)    * B0(:,:,:,iEast_,iBlock))

  ExB(:,:,:,iUp_)    =    EField(:,:,:,iEast_)  * B0(:,:,:,iNorth_,iBlock) - &
                        EField(:,:,:,iNorth_) * B0(:,:,:,iEast_,iBlock)

  ExB(:,:,:,iEast_)  = ExB(:,:,:,iEast_)  / (B0(:,:,:,iMag_,iBlock)**2)
  ExB(:,:,:,iNorth_) = ExB(:,:,:,iNorth_) / (B0(:,:,:,iMag_,iBlock)**2)
  ExB(:,:,:,iUp_)    = ExB(:,:,:,iUp_)    / (B0(:,:,:,iMag_,iBlock)**2)

  if (iDebugLevel > 5) then

     kmax = 1
     jmax = 1
     imax = 1
     maxi = 0.0

     do k=0,nAlts+1
        do i=0,nLats+1
           do j=0,nLons+1

              if (abs(ExB(j,i,k,iNorth_)) > maxi) then
                 write(*,*) "======> efield : ", i,j,k, &
                      EField(j,i,k,iEast_) * 1000.0, &
                      Potential(j+1,i,k,iBlock)/1000.0, &
                      Potential(j-1,i,k,iBlock)/1000.0, &
                      dLonDist_GB(j, i, k, iBlock), &
                      ExB(j,i,k,iNorth_), B0(j,i,k,iMag_,iBlock)
                 maxi = abs(ExB(j,i,k,iNorth_))
              endif

           enddo
        enddo
     enddo

  endif

  if (iDebugLevel > 3) then

     write(*,*) "====> Min ExB : ", minval(ExB(:,:,:,iEast_)),  &
          minval(ExB(:,:,:,iNorth_)),  minval(ExB(:,:,:,iUp_))
     write(*,*) "====> Max ExB : ", maxval(ExB(:,:,:,iEast_)),  &
          maxval(ExB(:,:,:,iNorth_)),  maxval(ExB(:,:,:,iUp_))

  endif

  call end_timing("calc_efield")

end subroutine calc_efield
