!==========================================================================
subroutine getMinMaxLatLon(minLat,maxLat, minLon,maxLon, rayface_2d, shifted)
  use ModMain, ONLY : nI,nJ
  implicit none

  real, intent(inout) :: minLat,maxLat, minLon,maxLon
  real, intent(inout), dimension(3,2,1:nI+1,1:nJ+1):: rayface_2d
  logical, intent(inout) :: shifted

  shifted=.false.
  minLat=minval(rayface_2d(1,1,1:nI+1,1:nJ+1), &
       MASK=   (rayface_2d(3,1,1:nI+1,1:nJ+1)==3.))
  maxLat=maxval(rayface_2d(1,1,1:nI+1,1:nJ+1), &
       MASK=   (rayface_2d(3,1,1:nI+1,1:nJ+1)==3.))
  minLon=minval(rayface_2d(2,1,1:nI+1,1:nJ+1), &
       MASK=   (rayface_2d(3,1,1:nI+1,1:nJ+1)==3.))
  maxLon=maxval(rayface_2d(2,1,1:nI+1,1:nJ+1), &
       MASK=   (rayface_2d(3,1,1:nI+1,1:nJ+1)==3.))
  if(maxLon-minLon > 180.)then
     shifted=.true.
     where(rayface_2d(2,1,1:nI+1,1:nJ+1)<180.) &
          rayface_2d(2,1,1:nI+1,1:nJ+1)=rayface_2d(2,1,1:nI+1,1:nJ+1)+360.
     minLon=minval(rayface_2d(2,1,1:nI+1,1:nJ+1), &
          MASK=   (rayface_2d(3,1,1:nI+1,1:nJ+1)==3.))
     maxLon=maxval(rayface_2d(2,1,1:nI+1,1:nJ+1), &
          MASK=   (rayface_2d(3,1,1:nI+1,1:nJ+1)==3.))
  end if
end subroutine getMinMaxLatLon

!==========================================================================
!==========================================================================
!==========================================================================
! Given rayface values on a face of a block of data, determine if and where
!  a Lat/Lon pair intersect the face.
subroutine find_rayface_location(dbg,iBLK,rLat,rLon,xL,yL,rayface_2d)
  use ModProcMH
  use ModMain, ONLY : nI,nJ
  use ModRaytrace
  implicit none

  integer, intent(in) :: iBLK
  real, intent(in) :: rLat,rLon
  real, intent(inout) :: xL,yL
  real, intent(in), dimension(3,2,1:nI+1,1:nJ+1):: rayface_2d
  logical, intent(in) :: dbg

  integer, parameter :: k=1
  integer :: i,j, i2,j2, iL,jL, nPts,nPts0
  real :: xp(3),yp(3),Lon(3), gradLat

  logical, dimension(nI+1,nJ+1,2) :: Lnear
  logical, dimension(nI+1,nJ+1) :: Tnear

  Lnear = .false.
  Tnear = .false.

  xL=0.; yL=0.

  do i=1,nI+1; do j=1,nJ+1
     if(rayface_2d(3,1,i,j)==3.) then
        !Check Lat
        if(abs(rayface_2d(1,1,i,j)-rLat)<1.E-6)then
           !Set point and all 8 around it to true
           do iL=-1,1; do jL=-1,1
              i2=max(1,min(nI+1,i+iL)) ; j2=max(1,min(nJ+1,j+jL))
              if(rayface_2d(3,1,i2,j2)==3.) Lnear(i2,j2,1)=.true.
           end do; end do
        else
           !Loop over all 8 around to see if Lat crosses
           do iL=-1,1; do jL=-1,1
              i2=max(1,min(nI+1,i+iL)) ; j2=max(1,min(nJ+1,j+jL))
              if(rayface_2d(3,1,i2,j2)==3.)then
                 if(  (rayface_2d(1,1,i2,j2)-rLat)/ &
                      (rayface_2d(1,1,i ,j )-rLat) <= 0. )then
                    Lnear(i,j,1)=.true. ; Lnear(i2,j2,1)=.true.
                 end if
              end if
           end do; end do
        end if

        !Check Lon
        if(abs(rayface_2d(2,1,i,j)-rLon)<1.E-6)then
           !Set point and all 8 around it to true
           do iL=-1,1; do jL=-1,1
              i2=max(1,min(nI+1,i+iL)) ; j2=max(1,min(nJ+1,j+jL))
              if(rayface_2d(3,1,i2,j2)==3.) Lnear(i2,j2,2)=.true.
           end do; end do
        else
           !Loop over all 8 around to see if Lon crosses
           do iL=-1,1; do jL=-1,1
              i2=max(1,min(nI+1,i+iL)) ; j2=max(1,min(nJ+1,j+jL))
              if(rayface_2d(3,1,i2,j2)==3.)then
                 if(  (rayface_2d(2,1,i2,j2)-rLon)/ &
                      (rayface_2d(2,1,i ,j )-rLon) <= 0. )then
                    Lnear(i,j,2)=.true. ; Lnear(i2,j2,2)=.true.
                 end if
              end if
           end do; end do
        end if
     end if
  end do; end do

  Tnear = (Lnear(:,:,1) .and. Lnear(:,:,2))

  if(dbg)then
     write(*,*)' '
     write(*,*)'PE=',iProc,' BLK=',iBLK
     write(*,*)'  Searching for value at Lat=',rLat,' Lon=',rLon
     write(*,*)' '
     write(*,*)' Lnear-1'
     do i=1,nI+1
        write(*,*) 'I=',i,' J=1- K=1  Lnear(i,j,1)=',Lnear(i,:,1)
     end do
     write(*,*)' Lnear-2'
     do i=1,nI+1
        write(*,*) 'I=',i,' J=1- K=1  Lnear(i,j,2)=',Lnear(i,:,2)
     end do
     write(*,*)' Tnear'
     do i=1,nI+1
        write(*,*) 'I=',i,' J=1- K=1  Tnear(i,j)=  ',Tnear(i,:)
     end do
     write(*,*)' X'
     do i=1,nI+1
        write(*,*) 'I=',i,' J=1- K=1  rayface_2d(1,2,i,j)=',rayface_2d(1,2,i,:)
     end do
     write(*,*)' Y'
     do i=1,nI+1
        write(*,*) 'I=',i,' J=1- K=1  rayface_2d(2,2,i,j)=',rayface_2d(2,2,i,:)
     end do
     write(*,*)' Lat'
     do i=1,nI+1
        write(*,*) 'I=',i,' J=1- K=1  rayface_2d(1,1,i,j)=',rayface_2d(1,1,i,:)
     end do
     write(*,*)' Lon'
     do i=1,nI+1
        write(*,*) 'I=',i,' J=1- K=1  rayface_2d(2,1,i,j)=',rayface_2d(2,1,i,:)
     end do
     write(*,*)' Status'
     do i=1,nI+1
        write(*,*) 'I=',i,' J=1- K=1  rayface_2d(3,1,i,j)=',rayface_2d(3,1,i,:)
     end do
     write(*,*)' '
  end if

  !Return if no points found
  if(.not.any(Tnear)) return

  !Loop over all Tnear points trying triangles
  do i=1,nI+1; do j=1,nJ+1
     if(Tnear(i,j))then
        !Try all four triangles surrounding this point
        if(dbg)write(*,*)'Trying point: ',i,j
        do iL=-1,1,2
           i2=i+iL
           if(i2<1 .or. i2>nI+1) CYCLE
           if(.not.Tnear(i2,j)) CYCLE

           nPts=0

           !Found good point, try it.
           if(dbg)write(*,*)'  1st point:',i2,j
           if(abs(rayface_2d(1,1,i2,j)-rayface_2d(1,1,i,j))>1.E-6)then
              gradLat=(rLat-rayface_2d(1,1,i,j))/ &
                   (rayface_2d(1,1,i2,j)-rayface_2d(1,1,i,j))
              if(gradLat>=0. .and. gradLat<=1.)then
                 nPts=nPts+1
                 xp(nPts) =    rayface_2d(1,2,i ,j ) + &
                      gradLat*(rayface_2d(1,2,i2,j )-rayface_2d(1,2,i,j))
                 yp(nPts) =    rayface_2d(2,2,i ,j ) + &
                      gradLat*(rayface_2d(2,2,i2,j )-rayface_2d(2,2,i,j))
                 Lon(nPts)=    rayface_2d(2,1,i ,j ) + &
                      gradLat*(rayface_2d(2,1,i2,j )-rayface_2d(2,1,i,j))
                 if(dbg)write(*,*)'  Using point:',nPts,xp(nPts),yp(nPts),Lon(nPts)
              end if
           end if
           nPts0=nPts

           do jL=-1,1,2
              nPts=nPts0
              j2=j+jL
              if(j2<1 .or. j2>nJ+1) CYCLE
              if(.not.Tnear(i,j2)) CYCLE

              !Found another good point, try it.
              if(dbg)write(*,*)'  2nd point:',i,j2
              if(abs(rayface_2d(1,1,i,j2)-rayface_2d(1,1,i,j))>1.E-6)then
                 gradLat=(rLat-rayface_2d(1,1,i,j))/ &
                      (rayface_2d(1,1,i,j2)-rayface_2d(1,1,i,j))
                 if(gradLat>=0. .and. gradLat<=1.)then
                    nPts=nPts+1
                    xp(nPts) =    rayface_2d(1,2,i ,j ) + &
                         gradLat*(rayface_2d(1,2,i ,j2)-rayface_2d(1,2,i,j))
                    yp(nPts) =    rayface_2d(2,2,i ,j ) + &
                         gradLat*(rayface_2d(2,2,i ,j2)-rayface_2d(2,2,i,j))
                    Lon(nPts)=    rayface_2d(2,1,i ,j ) + &
                         gradLat*(rayface_2d(2,1,i ,j2)-rayface_2d(2,1,i,j))
                    if(dbg)write(*,*)'  Using point:',nPts,xp(nPts),yp(nPts),Lon(nPts)
                 end if
              end if

              !Now try diagonal
              if(dbg)write(*,*)'  diagonal:'
              if(abs(rayface_2d(1,1,i,j2)-rayface_2d(1,1,i2,j))>1.E-6)then
                 gradLat=(rLat-rayface_2d(1,1,i2,j))/ &
                      (rayface_2d(1,1,i,j2)-rayface_2d(1,1,i2,j))
                 if(gradLat>=0. .and. gradLat<=1.)then
                    nPts=nPts+1
                    xp(nPts) =    rayface_2d(1,2,i2,j ) + &
                         gradLat*(rayface_2d(1,2,i ,j2)-rayface_2d(1,2,i2,j))
                    yp(nPts) =    rayface_2d(2,2,i2,j ) + &
                         gradLat*(rayface_2d(2,2,i ,j2)-rayface_2d(2,2,i2,j))
                    Lon(nPts)=    rayface_2d(2,1,i2,j ) + &
                         gradLat*(rayface_2d(2,1,i ,j2)-rayface_2d(2,1,i2,j))
                    if(dbg)write(*,*)'  Using point:',nPts,xp(nPts),yp(nPts),Lon(nPts)
                 end if
              end if

              if(dbg)write(*,*)' # Intersections:',nPts
              if(nPts<2) CYCLE

              !It hit right on it, return with answer
              if(Lon(2)==rLon)then
                 xL=xp(2); yL=yp(2)
                 if(dbg)then
                    write(*,*)' '
                    write(*,*)' Final Location w/o interpolation: ',xL,yL
                    write(*,*)' '
                 end if
                 return
              end if

              !Interpolate Lon to get final point
              if( (Lon(1)-rLon)/(Lon(2)-rLon) <= 0. )then
                 xL=xp(1)+((rLon-Lon(1))/(Lon(2)-Lon(1)))*(xp(2)-xp(1))
                 yL=yp(1)+((rLon-Lon(1))/(Lon(2)-Lon(1)))*(yp(2)-yp(1))
                 if(dbg)then
                    write(*,*)' '
                    write(*,*)' Final Location w/  interpolation: ',xL,yL
                    write(*,*)' '
                 end if
                 return
              end if
              if(dbg)write(*,*)' missed triangle:',rLon,Lon(1:2)
           end do
        end do
     end if
  end do; end do

end subroutine find_rayface_location

!==========================================================================
!==========================================================================
!==========================================================================
subroutine find_eqB(dbg,iBLK,xL,yL,zL,eqB)
  use ModMain, ONLY : nI,nJ,nK
  use ModAdvance, ONLY : Bx_,By_,Bz_,State_VGB
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK
  implicit none

  integer, intent(in) :: iBLK
  real, intent(in) :: xL,yL,zL
  real, intent(out) :: eqB
  logical, intent(in) :: dbg

  integer :: i,j,k, i1,i2, j1,j2, k1,k2
  real :: qb(3), dx1,dx2, dy1,dy2, dz1,dz2

  ! Get B0 values for location
  call get_b0(xL,yL,zL,qb)
  if(dbg)write(*,*)'B0 qb:',qb

  ! Determine cell indices corresponding to location qx
  i1=-9
  do i=0,nI; do j=0,nJ; do k=0,nK
     if(  x_BLK(i,j,k,iBLK) <= xL .and. x_BLK(i+1,j,k,iBLK) > xL .and. &
          y_BLK(i,j,k,iBLK) <= yL .and. y_BLK(i,j+1,k,iBLK) > yL .and. &
          z_BLK(i,j,k,iBLK) <= zL .and. z_BLK(i,j,k+1,iBLK) > zL )then
        i1=i; i2=i+1
        j1=j; j2=j+1
        k1=k; k2=k+1
     end if
  end do; end do; end do
  if(dbg)write(*,*)'IN:',xL,yL,zL, &
       '  X:',x_BLK(0,0,0,iBLK),x_BLK(nI+1,nJ+1,nK+1,iBLK), &
       '  Y:',y_BLK(0,0,0,iBLK),y_BLK(nI+1,nJ+1,nK+1,iBLK), &
       '  Z:',z_BLK(0,0,0,iBLK),z_BLK(nI+1,nJ+1,nK+1,iBLK), &
       '  ijk:',i1,i2,j1,j2,k1,k2
  if(i1==-9)then
     write(*,*)'ERROR in find_eqB'
     write(*,*)'IN:',xL,yL,zL, &
          '  X:',x_BLK(0,0,0,iBLK),x_BLK(nI+1,nJ+1,nK+1,iBLK), &
          '  Y:',y_BLK(0,0,0,iBLK),y_BLK(nI+1,nJ+1,nK+1,iBLK), &
          '  Z:',z_BLK(0,0,0,iBLK),z_BLK(nI+1,nJ+1,nK+1,iBLK), &
          '  ijk:',i1,i2,j1,j2,k1,k2
     stop
  end if

  ! Distance relative to the cell centers
  dx1=(xL-x_BLK(i1,j1,k1,iBLK))/(x_BLK(i2,j2,k2,iBLK)-x_BLK(i1,j1,k1,iBLK)); dx2=1.-dx1
  dy1=(yL-y_BLK(i1,j1,k1,iBLK))/(y_BLK(i2,j2,k2,iBLK)-y_BLK(i1,j1,k1,iBLK)); dy2=1.-dy1
  dz1=(zL-z_BLK(i1,j1,k1,iBLK))/(z_BLK(i2,j2,k2,iBLK)-z_BLK(i1,j1,k1,iBLK)); dz2=1.-dz1
  if(dbg)write(*,*)'dels:',dx1,dx2,dy1,dy2,dz1,dz2

  ! Add in interpolated B1 values
  qb(1)=qb(1)+&
       dx1*( dy1*( dz1*State_VGB(Bx_,i2,j2,k2,iBLK) + dz2*State_VGB(Bx_,i2,j2,k1,iBLK))+&
             dy2*( dz1*State_VGB(Bx_,i2,j1,k2,iBLK) + dz2*State_VGB(Bx_,i2,j1,k1,iBLK)))+&
       dx2*( dy1*( dz1*State_VGB(Bx_,i1,j2,k2,iBLK) + dz2*State_VGB(Bx_,i1,j2,k1,iBLK))+&
             dy2*( dz1*State_VGB(Bx_,i1,j1,k2,iBLK) + dz2*State_VGB(Bx_,i1,j1,k1,iBLK)))
  qb(2)=qb(2)+&
       dx1*( dy1*( dz1*State_VGB(By_,i2,j2,k2,iBLK) + dz2*State_VGB(By_,i2,j2,k1,iBLK))+&
             dy2*( dz1*State_VGB(By_,i2,j1,k2,iBLK) + dz2*State_VGB(By_,i2,j1,k1,iBLK)))+&
       dx2*( dy1*( dz1*State_VGB(By_,i1,j2,k2,iBLK) + dz2*State_VGB(By_,i1,j2,k1,iBLK))+&
             dy2*( dz1*State_VGB(By_,i1,j1,k2,iBLK) + dz2*State_VGB(By_,i1,j1,k1,iBLK)))
  qb(3)=qb(3)+&
       dx1*( dy1*( dz1*State_VGB(Bz_,i2,j2,k2,iBLK) + dz2*State_VGB(Bz_,i2,j2,k1,iBLK))+&
             dy2*( dz1*State_VGB(Bz_,i2,j1,k2,iBLK) + dz2*State_VGB(Bz_,i2,j1,k1,iBLK)))+&
       dx2*( dy1*( dz1*State_VGB(Bz_,i1,j2,k2,iBLK) + dz2*State_VGB(Bz_,i1,j2,k1,iBLK))+&
             dy2*( dz1*State_VGB(Bz_,i1,j1,k2,iBLK) + dz2*State_VGB(Bz_,i1,j1,k1,iBLK)))
  if(dbg)write(*,*)'qb:',qb
  eqB=sqrt(qb(1)**2+qb(2)**2+qb(3)**2)
  if(dbg)write(*,*)'eqB:',eqB

end subroutine find_eqB
!==========================================================================
