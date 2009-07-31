Module ModWriteTec
  implicit none
  logical :: DoWriteTec = .false., IsFirst
  
contains
  subroutine write_tec(Time, Flux_C, eBound_I)
    use rbe_cgrid, ONLY: Sina_C => gridy, Energy_I => gride, &
         nLat => ir, nLon => ip, nAngle => ig, nEnergy => je
    use rbe_cfield,ONLY: Xeq_C => xo, Yeq_C => yo
    use ModIoUnit, ONLY: UnitTmp_
    implicit none
    real, intent(in) :: Time
    real, intent(in) :: Flux_C(nLat,nLon,nEnergy,nAngle), eBound_I(nEnergy+1)
    
    integer :: iLat, iLon, iEnergy, iAngle, iE ! loop variables 
    Character(len=100) :: NameFile             !output file name
    Character(len=600) :: NameVar              !output Var name
    Character(len=100)  :: NameVarTemp
    real, allocatable  :: Cosa_C(:), dE_I(:), AvFlux_C(:,:,:), &
         Anisotropy_C(:,:,:), dAngle_I(:)
    real :: RightSina, LeftSina, FluxPara, FluxPerp, ParaPlusPerp
    integer :: TimeOut
    !--------------------------------------------------------------------------

    ! Allocate needed variables
    if (.not.allocated(Cosa_C)) allocate(Cosa_C(nAngle))
    if (.not.allocated(dE_I)) allocate(dE_I(nEnergy))
    if (.not.allocated(AvFlux_C)) allocate(AvFlux_C(nLat,nLon,nEnergy))
    if (.not.allocated(Anisotropy_C)) allocate(Anisotropy_C(nLat,nLon,nEnergy))
    if (.not.allocated(dAngle_I)) allocate(dAngle_I(nAngle))

    ! Set plot file name and open file
    TimeOut=int(Time)
    write(NameFile,"(a,i8.8,a)") &
         'RB/EquatorialFlux_',TimeOut,'.dat'
    open(UnitTmp_,FILE=NameFile)
    

    do iEnergy=1,nEnergy
       if (iEnergy == 1) then
          write(NameVar,"(a,f7.2,a,f7.2,a)") &
               '"X", "Y", "Z", "',eBound_I(iEnergy),'-',eBound_I(iEnergy+1),'/cm2/s", ' 
       else
          write(NameVarTemp,"(a,f7.2,a,f7.2,a)") &
               ', "',eBound_I(iEnergy),'-',eBound_I(iEnergy+1),'/cm2/s", '
          NameVar = trim(NameVar)//trim(NameVarTemp)
       endif
       write(NameVarTemp,"(a,f7.2,a,f7.2,a)") &
            ' "',eBound_I(iEnergy),'-',eBound_I(iEnergy+1),'PA"'
       NameVar = trim(NameVar)//trim(NameVarTemp)
    enddo
    
    ! Write file header
    write(UnitTmp_,'(a)') 'VARIABLES = '//NameVar
    
    write(UnitTmp_,'(a,i8,a,i8,a)') &
         'Zone I=', nLat, ', J=', nLon,', DATAPACKING=POINT'
    
    ! Get Energy spacing
    do iEnergy = 1,nEnergy
       dE_I(iEnergy) = eBound_I(iEnergy+1) - eBound_I(iEnergy)
    enddo
    
    ! Get Angle Space
    Cosa_C(:)=cos(asin(Sina_C(:)))
    do iAngle = 1,nAngle
       if(iAngle == 1) LeftSina=0
       if(iAngle >  1) LeftSina=0.5*(Sina_C(iAngle)+Sina_C(iAngle-1))
       if(iAngle == nAngle) RightSina=1.0
       if(iAngle <  nAngle) RightSina=0.5*(Sina_C(iAngle)+Sina_C(iAngle+1))
       dAngle_I(iAngle) = sqrt(1.0-LeftSina**2.0)-sqrt(1.0-RightSina**2.0)
    enddo

    AvFlux_C(:,:,:) = 0.0
    do iLat = 1,nLat ; do iLon =1,nLon ; do iE=1,nEnergy
       FluxPara = 0.0
       FluxPerp = 0.0
       do iAngle = 1,nAngle
          ! PA averaged flux as a function of lat, lon, energy bin
          AvFlux_C(iLat,iLon,iE) =    &
               AvFlux_C(iLat,iLon,iE) &
               + Flux_C(iLat,iLon,iE,iAngle)*dAngle_I(iAngle)
          FluxPara =    &
              FluxPara &
              + Flux_C(iLat,iLon,iE,iAngle)*dAngle_I(iAngle)*Cosa_C(iAngle)**2.0
          FluxPerp =    &
              FluxPerp &
              + 0.5*Flux_C(iLat,iLon,iE,iAngle)&
              * dAngle_I(iAngle)*Sina_C(iAngle)**2.0
       enddo
       ParaPlusPerp = FluxPara+FluxPerp
       if (ParaPlusPerp > 0.0) Anisotropy_C(iLat,iLon,iE) = &
            (FluxPerp-FluxPara)/ParaPlusPerp
       
       ! Get average flux over PA and Energy 
       AvFlux_C(iLat,iLon,iE) = AvFlux_C(iLat,iLon,iE)*dE_I(iE)
    enddo; enddo; enddo
    
    ! Write data to plot file
    do iLat = 1,nLat 
       do iLon =1,nLon
          write(UnitTmp_,fmt="(30(E14.6))") &
            Xeq_C(iLat,iLon), Yeq_C(iLat,iLon), 0.0,&
            (AvFlux_C(iLat,iLon,iE), Anisotropy_C(iLat, iLon, iE), iE=1,nEnergy)
       enddo
    enddo
    close(UnitTmp_)

    ! deallocate variables to save memory
    deallocate(Cosa_C,dE_I,AvFlux_C,Anisotropy_C,dAngle_I)

    
  end subroutine write_tec
  
end Module ModWriteTec
