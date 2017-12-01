!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

!\
! This routine finds a point on in the Spherical system, given
! a Theta, Phi: 
! LocIn(1) = Phi
! LocIn(2) = Theta
! It returns a 4 element array:
! LocOut(1) = Index of Block
! LocOut(2) = Index of Longitude
! LocOut(3) = Index of Latitude
! LocOut(4) = Multiplication factor for Longitude
! LocOut(5) = Multiplication factor for Latitude
!/

subroutine IE_FindPoint(LocIn, LocOut, IsNormalGrid, iError)

  use ModErrors
  use ModIE_Interface
  use ModConst
  use ModIonosphere

  implicit none

  real, dimension(2), intent(in)  :: LocIn
  real, dimension(5), intent(out) :: LocOut
  logical, intent(in) :: IsNormalGrid
  integer, intent(out) :: iError
  real :: MLTIn, LatIn
  integer :: j,i, iBLK

  integer :: jTemp_I(1)


  LocOut = -1.0

  iError = 0

  !\
  ! Check to see if the point is even on the grid.
  !/

  MLTIn = mod(LocIn(1),24.0)

  LatIn = LocIn(2)
  if (LatIn > 90.0-IONO_TOLER*cRadToDeg) then
     LatIn = 180.0 - 2*IONO_TOLER*cRadToDeg - LatIn
     MLTIn = mod(MLTIn+12.0,24.0)
  endif
  if (LatIn < -90.0+IONO_TOLER*cRadToDeg) then
     LatIn = -180.0 + 2*IONO_TOLER*cRadToDeg - LatIn
     MLTIn = mod(MLTIn+12.0,24.0)
  endif

  if (MLTIn > 24.0 .or. MLTIn < 0 .or.&
       LatIn > 90.0 -IONO_TOLER*cRadToDeg.or.&
       LatIn < -90.0+IONO_TOLER*cRadToDeg) then
     iError = ecPointOutofRange_
     return
  endif
  if (IsNormalGrid) then

     iBLK = 1
     i    = IEi_HavenLats
     j    = IEi_HavenMLTs

     do while (iBLK <= IEi_HavenBLKs .and. LocOut(1) < 0.0)

        if ((LatIn <  maxval(IEr3_HaveLats(1,:,iBLK)) .and. &
             LatIn >= minval(IEr3_HaveLats(1,:,iBLK))) .and. &
             ((MLTIn <  maxval(IEr3_HaveMLTs(:,1,iBLK)) .and. &
             MLTIn >= minval(IEr3_HaveMLTs(:,1,iBLK))).or.&
             (MLTIn>=maxval(IEr3_HaveMLTs(:,1,iBLK)).and.MLTIn<=24.0))) then
           LocOut(1) = iBLK
        else
           iBLK = iBLK + 1
        endif

     enddo
     if(LocOut(1)<=cZero)then
        write(*,*)'Wrong LocOut(1)',LocOut(1)
        write(*,*)'LatIn,MLTIn',LatIn,MLTIn
        write(*,*)'Lat limits:',minval(IEr3_HaveLats(1,:,:)),&
              maxval(IEr3_HaveLats(1,:,:))
        write(*,*)'MLT limits:',minval(IEr3_HaveMLTs(:,1,:)),&
              maxval(IEr3_HaveMLTs(:,1,:))
        call CON_stop('IE_FindPoint: Point is out of range.')
     end if
     j    = 1
     iBLK = nint(LocOut(1))
     i    = IEi_HavenLats

     do while (j < IEi_HavenMLTs .and. LocOut(2) <0.0)

        jTemp_I = maxloc(IEr3_HaveMLTs(:,1,iBLK))

        if ((MLTIn <  IEr3_HaveMLTs(j+1,1,iBLK) .and. &
            MLTIn >= IEr3_HaveMLTs(j,1,iBLK)).or.&
            MLTIn>=maxval(IEr3_HaveMLTs(:,1,iBLK)).and.&
            j==jTemp_I(1))then
           LocOut(2) = j 
        else
           j = j + 1
        endif

     enddo

     i = 1
     iBLK = nint(LocOut(1))
     j    = IEi_HavenMLTs

     do while (i < IEi_HavenLats .and. LocOut(3) < 0.0)

        if ((LatIn <  IEr3_HaveLats(j,i+1,iBLK) .and. &
             LatIn >= IEr3_HaveLats(j,i,iBLK)) .or. &
            (LatIn >  IEr3_HaveLats(j,i+1,iBLK) .and. &
             LatIn <= IEr3_HaveLats(j,i,iBLK))) then
           LocOut(3) = i
        else
           i = i + 1 
        endif

     enddo

     iBLK = nint(LocOut(1))
     j    = nint(LocOut(2))
     i    = nint(LocOut(3))

     if (LocOut(1) > -cHalf .and. LocOut(2) > -cHalf .and. LocOut(3) > -cHalf) then
        
        if(MLTIn<maxval(IEr3_HaveMLTs(:,1,iBLK)))then
           LocOut(4) = (MLTIn                    -IEr3_HaveMLTs(j,i,iBLK))/&
                       (IEr3_HaveMLTs(j+1,i,iBLK)-IEr3_HaveMLTs(j,i,iBLK))
        else
           LocOut(4) = (MLTIn                    -IEr3_HaveMLTs(j,i,iBLK))/&
                       (24.0-IEr3_HaveMLTs(j,i,iBLK))
        end if
        LocOut(5) = (LatIn                    -IEr3_HaveLats(j,i,iBLK))/&
                    (IEr3_HaveLats(j,i+1,iBLK)-IEr3_HaveLats(j,i,iBLK))

     endif

  else

     iBLK = 1
     do while (iBLK <= IEi_HavenBLKs .and. LocOut(1) == -1.0)
        j = 1
        do while (j < IEi_HavenMLTs .and. LocOut(2) == -1.0)
           i = 1
           do while (i < IEi_HavenLats .and. LocOut(3) == -1.0)

              !\
              ! Check to see if the point is within the current cell
              !/

              if (((LatIn <  IEr3_HaveLats(j,i+1,iBLK) .and. &
                   LatIn >= IEr3_HaveLats(j,i,iBLK)) .or. &
                   (LatIn >  IEr3_HaveLats(j,i+1,iBLK) .and. &
                   LatIn <= IEr3_HaveLats(j,i,iBLK))) .and. &
                   MLTIn <  IEr3_HaveMLTs(j+1,i,iBLK) .and. &
                   MLTIn >= IEr3_HaveMLTs(j,i,iBLK)) then

                 !\
                 ! If it is, then store the cell number and calculate
                 ! the interpolation coefficients.
                 !/

                 LocOut(1) = iBLK
                 LocOut(2) = j
                 LocOut(3) = i

                 LocOut(4) = (MLTIn                    -IEr3_HaveMLTs(j,i,iBLK))/&
                             (IEr3_HaveMLTs(j+1,i,iBLK)-IEr3_HaveMLTs(j,i,iBLK))
                 LocOut(5) = (LatIn                    -IEr3_HaveLats(j,i,iBLK))/&
                             (IEr3_HaveLats(j,i+1,iBLK)-IEr3_HaveLats(j,i,iBLK))

                 iBLK = IEi_HavenBLKs
                 j = IEi_HavenMLTs
                 i = IEi_HavenLats

              endif

              i = i + 1

           enddo

           j = j + 1

        enddo

        iBLK = iBLK + 1

     enddo

  endif

end subroutine IE_FindPoint
