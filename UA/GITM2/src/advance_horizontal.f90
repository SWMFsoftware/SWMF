subroutine advance_horizontal(iBlock)

  use ModConstants, only : Gamma, pi
  use ModSizeGitm
  use ModPlanet, only : nSpecies, nIonsAdvect, OmegaBody
  use ModGITM
  use ModInputs, only : UseIonAdvection, iDebugLevel
  
  implicit none

  integer, intent(in) :: iBlock

  integer :: iAlt, iIon, iLon, iLat, iSpecies
  logical :: IsFound

  real :: MaxDiff

  real :: cp_C(1:nLons,1:nLats)
  real :: Rho_C(-1:nLons+2,-1:nLats+2)
  real :: Temp_C(-1:nLons+2,-1:nLats+2)
  real :: Vel_CD(-1:nLons+2,-1:nLats+2,3)
  real :: IVel_CD(-1:nLons+2,-1:nLats+2,3)
  real :: Num_CV(-1:nLons+2,-1:nLats+2,nSpecies)
  real :: INum_CV(-1:nLons+2,-1:nLats+2,nIonsAdvect)

  real :: NewRho_C(nLons,nLats)
  real :: NewTemp_C(nLons,nLats)
  real :: NewVel_CD(nLons,nLats, 3)
  real :: NewNum_CV(nLons,nLats, nSpecies)
  real :: NewINum_CV(nLons,nLats, nIonsAdvect)

  MaxDiff = 0.0

  call report("advance_horizontal",2)

  do iAlt=1,nAlts

     cp_c   = cp(:,:,iAlt,iBlock)
     Rho_C  = Rho(:,:,iAlt,iBlock)
     Vel_CD = Velocity(:,:,iAlt,:,iBlock)
     Num_CV = nDensityS(:,:,iAlt,1:nSpecies,iBlock)
     Temp_C = Temperature(:,:,iAlt,iBlock)

     IVel_CD = IVelocity(:,:,iAlt,:,iBlock)
     INum_CV = IDensityS(:,:,iAlt,1:nIonsAdvect,iBlock)
     
     NewRho_C  = Rho_C(1:nLons,1:nLats)
     NewVel_CD = Vel_CD(1:nLons,1:nLats,:)
     NewNum_CV = Num_CV(1:nLons,1:nLats,:)
     NewTemp_C = Temp_C(1:nLons,1:nLats)
     NewINum_CV = INum_CV(1:nLons,1:nLats,:)

     call horizontal_solver

     Rho(1:nLons,1:nLats,iAlt,iBlock)                     = NewRho_C
     Velocity(1:nLons,1:nLats,iAlt,:,iBlock)              = NewVel_CD
     Temperature(1:nLons,1:nLats,iAlt,iBlock)             = NewTemp_C

     if (minval(NewNum_CV) < 0.0) then
        write(*,*) "Negative Density after horizontal advection!!"
        write(*,*) "Correcting...."
        do iLon = 1, nLons
           do iLat = 1, nLats
              IsFound = .false.
              do iSpecies = 1, nSpecies
                 if (NewNum_CV(iLon, iLat, iSpecies) < 0.0) then
                 write(*,*) "Species : ", iSpecies, iLon, iLat, iBlock
stop
                    NewNum_CV(iLon, iLat, iSpecies) = 1.0
                    IsFound=.true.
                 endif
              enddo
              If (IsFound) then
                 Rho(iLon,iLat,iAlt,iBlock) = 0.0
                 do iSpecies = 1, nSpecies
                    Rho(iLon,iLat,iAlt,iBlock) = &
                         Rho(iLon,iLat,iAlt,iBlock) + &
                         NewNum_CV(iLon, iLat, iSpecies)*Mass(iSpecies)
                 enddo
              endif
           enddo
        enddo
     endif
 
     nDensityS(1:nLons,1:nLats,iAlt,1:nSpecies,iBlock)    = NewNum_CV

    do iLon = 1, nLons
        do iLat = 1, nLats
           if (abs(IDensityS(iLon,iLat,iAlt,1,iBlock) - &
                NewINum_CV(iLon,iLat,1))/NewINum_CV(iLon,iLat,1) &
                > MaxDiff) then
!              write(*,*) "MaxDiff: ",iLon, iLat, iAlt, &
!                   MaxDiff, &
!                   IDensityS(iLon,iLat,iAlt,1,iBlock), &
!                   NewINum_CV(iLon,iLat,1),&
!                   IVel_CD(iLon,iLat,1),IVel_CD(iLon,iLat,2)
              MaxDiff = abs(IDensityS(iLon,iLat,iAlt,1,iBlock) - &
                NewINum_CV(iLon,iLat,1))/NewINum_CV(iLon,iLat,1)
           endif
        enddo
     enddo

!     MaxDiff = max(MaxDiff, &
!          maxval(abs((IDensityS(1:nLons,1:nLats,iAlt,1:nIonsAdvect,iBlock) - &
!          NewINum_CV)/NewINum_CV)))

     if (UseIonAdvection) then

        if (minval(NewINum_CV) < 1.0e2) then
!           write(*,*) "Negative Ion Density after horizontal advection!!"
!           write(*,*) "Correcting...."
           do iLon = 1, nLons
              do iLat = 1, nLats
                 do iIon = 1, nIonsAdvect
                    if (NewINum_CV(iLon, iLat, iIon) < 1.0e2) then
!                       write(*,*) "Location : ", &
!                            Longitude(iLon,iBlock)*180/pi, &
!                            Latitude(iLon,iBlock)*180/pi, &
!                            Altitude(iAlt)/1000.0, iIon
                       NewINum_CV(iLon, iLat, iIon) = 1.0e2
                    endif
                 enddo
              enddo
           enddo
        endif

        IDensityS(1:nLons,1:nLats,iAlt,1:nIonsAdvect,iBlock) = NewINum_CV
        !\
        ! New Electron Density
        !/
        IDensityS(1:nLons,1:nLats,iAlt,ie_,iBlock) = 0.0
        do iIon = 1, nIons-1
           IDensityS(1:nLons,1:nLats,iAlt,ie_,iBlock) = &
                IDensityS(1:nLons,1:nLats,iAlt,ie_,iBlock) + &
                IDensityS(1:nLons,1:nLats,iAlt,iIon,iBlock)
        enddo
     endif

  end do

  if (iDebugLevel > 2) &
       write(*,*) "===> MaxDiff : ", MaxDiff, maxval(abs(IVel_CD)), &
       maxval(IDensityS(1:nLons,1:nLats,1:nAlts,1:nIonsAdvect,iBlock)), dt

contains

  subroutine horizontal_solver

    use ModInputs, only: UseCoriolis

    ! Solve horizontal equations for a single block

    integer :: iDim, iSpc

    real, dimension(nLons,nLats)          :: &
         GradLonRho_C, GradLatRho_C,   DiffRho_C
    real, dimension(nLons,nLats)          :: &
         GradLonTemp_C, GradLatTemp_C, DiffTemp_C, DivVel_C
    real, dimension(nLons,nLats,3)        :: &
         GradLonVel_CD, GradLatVel_CD, DiffVel_CD
    real, dimension(nLons,nLats,nSpecies) :: &
         GradLonNum_CV, GradLatNum_CV, DiffNum_CV
    real, dimension(nLons,nLats,nIonsAdvect) :: &
         GradLonINum_CV, GradLatINum_CV, DiffINum_CV

    real :: CoriolisSin, CoriolisCos, CentrifugalParameter

    real :: RhoTest

    integer :: iLon, iLat

    call calc_rusanov_horizontal( &
         Rho_C,  GradLonRho_C,  GradLatRho_C,  DiffRho_C)
    call calc_rusanov_horizontal( &
         Temp_C, GradLonTemp_C, GradLatTemp_C, DiffTemp_C)

    do iDim = 1,3
       call calc_rusanov_horizontal( Vel_CD(:,:,iDim), &
            GradLonVel_CD(:,:,iDim), GradLatVel_CD(:,:,iDim), &
            DiffVel_CD(:,:,iDim))
    end do

    do iLat = 1, nLats
       DivVel_C(:,iLat) = &
            GradLatVel_CD(:,iLat,iNorth_) + GradLonVel_CD(:,iLat,iEast_) &
            + TanLatitude(iLat,iBlock) * Vel_CD(1:nLons,iLat,iNorth_) / &
            RadialDistance(iAlt)
    end do

    do iSpc = 1,nSpecies
       call calc_rusanov_horizontal(Num_CV(:,:,iSpc), &
            GradLonNum_CV(:,:,iSpc), GradLatNum_CV(:,:,iSpc), &
            DiffNum_CV(:,:,iSpc))
    end do

    do iSpc = 1,nIonsAdvect
       call calc_rusanov_horizontal(INum_CV(:,:,iSpc), &
            GradLonINum_CV(:,:,iSpc), GradLatINum_CV(:,:,iSpc), &
            DiffINum_CV(:,:,iSpc))
    end do

    do iLat=1,nLats

       CoriolisSin = sin(Latitude(iLat,iBlock)) * 2 * OmegaBody
       CoriolisCos = cos(Latitude(iLat,iBlock)) * 2 * OmegaBody

       CentrifugalParameter = OmegaBody**2 * cos(Latitude(iLat,iBlock)) * &
            sin(Latitude(iLat,iBlock))

       do iLon=1,nLons

          ! drho/dt = -(rho*div V + grad rho * V)
          !         = -[rho*(dv/dLat - v*tan(Lat) + (du/dLon)/cos(Lat))
          !             drho/dLat * vLat + drho/dLon * vLon/cos(Lat)]
          !            / r

          NewRho_C(iLon,iLat) = NewRho_C(iLon,iLat) - Dt * ( &
               Rho_C(iLon,iLat) * DivVel_C(iLon,iLat) &
               + GradLatRho_C(iLon,iLat)*Vel_CD(iLon,iLat,iNorth_) & 
               + GradLonRho_C(iLon,iLat)*Vel_CD(iLon,iLat,iEast_)) & 
               + Dt * DiffRho_C(iLon,iLat)

          do iSpc = 1, nSpecies
             NewNum_CV(iLon,iLat,iSpc) = NewNum_CV(iLon,iLat,iSpc) - Dt * ( &
                  Num_CV(iLon,iLat,iSpc) * DivVel_C(iLon,iLat) &
                  + GradLatNum_CV(iLon,iLat,iSpc)*Vel_CD(iLon,iLat,iNorth_) & 
                  + GradLonNum_CV(iLon,iLat,iSpc)*Vel_CD(iLon,iLat,iEast_)) & 
                  + Dt * DiffNum_CV(iLon,iLat,iSpc)
          enddo

          do iSpc = 1, nIonsAdvect
             NewINum_CV(iLon,iLat,iSpc) = NewINum_CV(iLon,iLat,iSpc) - Dt * ( &
!                  INum_CV(iLon,iLat,iSpc) * DivVel_C(iLon,iLat) &
                  + GradLatINum_CV(iLon,iLat,iSpc)*IVel_CD(iLon,iLat,iNorth_) &
                  + GradLonINum_CV(iLon,iLat,iSpc)*IVel_CD(iLon,iLat,iEast_)) &
                  + Dt * DiffINum_CV(iLon,iLat,iSpc)
          enddo

          RhoTest = sum(Mass(1:nSpecies) * NewNum_CV(iLon,iLat,1:nSpecies))

          if (abs(RhoTest - NewRho_C(iLon,iLat))/RhoTest > 0.1) then
             write(*,*) "Problem!! ", RhoTest, NewRho_C(iLon,iLat)
!             call stop_gitm("Have to stop")
          endif

          ! dv_phi/dt = -(V grad V + (1/rho) grad P)_phi
          ! (1/rho) grad p = grad T + T/rho grad rho 

          NewVel_CD(iLon,iLat,iEast_) = NewVel_CD(iLon,iLat,iEast_) - Dt * ( &
               Vel_CD(iLon,iLat,iNorth_)*GradLatVel_CD(iLon,iLat,iEast_) + &
               Vel_CD(iLon,iLat,iEast_ )*GradLonVel_CD(iLon,iLat,iEast_) + &
               Vel_CD(iLon,iLat,iEast_)*(Vel_CD(iLon,iLat,iUp_) &
               - TanLatitude(iLat,iBlock)*Vel_CD(iLon,iLat,iNorth_)) &
               / RadialDistance(iAlt) + &
               GradLonTemp_C(iLon,iLat) + & 
               GradLonRho_C(iLon,iLat)*Temp_C(iLon,iLat)/Rho_C(iLon,iLat)) &
               + Dt * DiffVel_CD(iLon,iLat,iEast_)


          ! dv_theta/dt = -(V grad V + (1/rho) grad P)_theta
          ! (1/rho) grad p = grad T + T/rho grad rho 

!          if (iLon == 1) &
!          write(*,*) "vel before adv : ",NewVel_CD(iLon,iLat,iNorth_) 

          NewVel_CD(iLon,iLat,iNorth_) = NewVel_CD(iLon,iLat,iNorth_) &
               - Dt * ( &
               Vel_CD(iLon,iLat,iNorth_)*GradLatVel_CD(iLon,iLat,iNorth_) + &
               Vel_CD(iLon,iLat,iEast_ )*GradLonVel_CD(iLon,iLat,iNorth_) + &
               (Vel_CD(iLon,iLat,iNorth_)*Vel_CD(iLon,iLat,iUp_) &
               + TanLatitude(iLat,iBlock)*Vel_CD(iLon,iLat,iEast_)**2 &
               ) / RadialDistance(iAlt) + &
               GradLatTemp_C(iLon,iLat) + & 
               GradLatRho_C(iLon,iLat)*Temp_C(iLon,iLat)/Rho_C(iLon,iLat)) &
               + Dt * DiffVel_CD(iLon,iLat,iNorth_)

!          if (iLon == 1) then
!             write(*,*) "vel before cor : ",NewVel_CD(iLon,iLat,iNorth_) 
!          endif

          ! dv_r/dt = -(V grad V)_r

          ! Need to add the Vertical Velocity here!!!!!
          NewVel_CD(iLon,iLat,iUp_) = NewVel_CD(iLon,iLat,iUp_) &
               - Dt * ( &
               Vel_CD(iLon,iLat,iNorth_)*GradLatVel_CD(iLon,iLat,iUp_) + &
               Vel_CD(iLon,iLat,iEast_ )*GradLonVel_CD(iLon,iLat,iUp_)) &
               + Dt * DiffVel_CD(iLon,iLat,iUp_)

          if (UseCoriolis) then

             NewVel_CD(iLon,iLat,iEast_) = NewVel_CD(iLon,iLat,iEast_) + Dt*( &
                  + CoriolisSin * Vel_CD(iLon,iLat,iNorth_) &
                  - CoriolisCos * Vel_CD(iLon,iLat,iUp_))

             NewVel_CD(iLon,iLat,iNorth_) = NewVel_CD(iLon,iLat,iNorth_) - &
                  Dt*( &
                  CentrifugalParameter * RadialDistance(iAlt) &
                  + CoriolisSin * Vel_CD(iLon,iLat,iEast_))
          endif

!          if (iLon == 1) then
!             write(*,*) "dt : ",dt, iAlt
!             write(*,*) "gradlon: ",GradLonVel_CD(iLon,iLat,iNorth_)
!             write(*,*) "tanlat : ",&
!                  TanLatitude(iLat,iBlock)*Vel_CD(iLon,iLat,iEast_)**2 &
!                  / RadialDistance(iAlt)
!             write(*,*) "gradt  : ",GradLatTemp_C(iLon,iLat)
!             write(*,*) "gradrho: ", &
!                  GradLatRho_C(iLon,iLat)*Temp_C(iLon,iLat)/Rho_C(iLon,iLat)
!             write(*,*) "diff   : ",DiffVel_CD(iLon,iLat,iNorth_)
!            write(*,*) "centr  : ",CentrifugalParameter * RadialDistance(iAlt)
!            write(*,*) "cori   : ",CoriolisParameter *Vel_CD(iLon,iLat,iEast_)
!             write(*,*) "vel after cor : ",NewVel_CD(iLon,iLat,iNorth_)
!          endif


          ! dT/dt = -(V.grad T + (gamma - 1) T div V

!          NewTemp_C(iLon,iLat) = NewTemp_C(iLon,iLat)

!          NewTemp_C(iLon,iLat) = NewTemp_C(iLon,iLat) - Dt * ( &
!               + GradLatTemp_C(iLon,iLat)*Vel_CD(iLon,iLat,iNorth_) & 
!               + GradLonTemp_C(iLon,iLat)*Vel_CD(iLon,iLat,iEast_)) & 
!               + Dt * DiffTemp_C(iLon,iLat)
          
!          NewTemp_C(iLon,iLat) = NewTemp_C(iLon,iLat) - Dt * ( &
!               Temp_C(iLon,iLat) * DivVel_C(iLon,iLat)/cp_c(iLon,iLat) &
!               + GradLatTemp_C(iLon,iLat)*Vel_CD(iLon,iLat,iNorth_) & 
!               + GradLonTemp_C(iLon,iLat)*Vel_CD(iLon,iLat,iEast_)) & 
!               + Dt * DiffTemp_C(iLon,iLat)

          NewTemp_C(iLon,iLat) = NewTemp_C(iLon,iLat) - Dt * ( &
               (gamma-1) * Temp_C(iLon,iLat) * DivVel_C(iLon,iLat) &
               + GradLatTemp_C(iLon,iLat)*Vel_CD(iLon,iLat,iNorth_) & 
               + GradLonTemp_C(iLon,iLat)*Vel_CD(iLon,iLat,iEast_)) & 
               + Dt * DiffTemp_C(iLon,iLat)

       end do
    end do

  end subroutine horizontal_solver


  !=====================================================================
  subroutine calc_rusanov_horizontal(Var, GradLonVar, GradLatVar, DiffVar)

    use ModSizeGitm
    use ModGITM

    implicit none

    real, intent(in)  :: Var(-1:nLons+2, -1:nLats+2)
    real, intent(out) :: GradLonVar(nLons, nLats), GradLatVar(nLons, nLats)
    real, intent(out) :: DiffVar(nLons, nLats)

    real, dimension(1:nLats+1) :: VarNorth, VarSouth, DiffLatFlux
    real, dimension(1:nLons+1) :: VarEast,  VarWest,  DiffLonFlux

    integer :: iLon, iLat

    ! Calculate gradient and diffusive flux with respect to latitude
    do iLon = 1, nLons
       call calc_GITM_facevalues(nLats, latitude(:,iBlock), Var(iLon,:), &
            VarSouth, VarNorth)
       ! Gradient based on averaged Left/Right values

       GradLatVar(iLon,:) = 0.5 * &
            (VarSouth(2:nLats+1)+VarNorth(2:nLats+1) - &
            VarSouth(1:nLats)-VarNorth(1:nLats)) &
            / dLatDist_GB(1:nLats,iAlt,iBlock)

       ! Rusanov/Lax-Friedrichs diffusive term
       DiffLatFlux = 0.5 * max(&
            cMax_GDB(iLon,0:nLats  ,iAlt,iNorth_,iBlock),&
            cMax_GDB(iLon,1:nLats+1,iAlt,iNorth_,iBlock)) &
            * (VarNorth - VarSouth)

       DiffVar(iLon,:) = (DiffLatFlux(2:nLats+1) - DiffLatFlux(1:nLats)) &
            / dLatDist_GB(1:nLats,iAlt,iBlock)
    end do

    ! Calculate gradient and diffusive flux with respect to longitude
    do iLat = 1, nLats
       call calc_GITM_facevalues(nLons, Longitude(:,iBlock), Var(:,iLat), &
            VarWest, VarEast)

       ! Gradient based on averaged West/East values

       GradLonVar(:,iLat) = 0.5 * &
            (VarWest(2:nLons+1) + VarEast(2:nLons+1) - &
            VarWest(1:nLons)    - VarEast(1:nLons)) &
            / dLonDist_GB(iLat,iAlt,iBlock)

       ! Rusanov/Lax-Friedrichs diffusive term
       DiffLonFlux = 0.5 * max(&
            cMax_GDB(0:nLons,   iLat, iAlt, iEast_, iBlock), &
            cMax_GDB(1:nLons+1, iLat, iAlt, iEast_, iBlock)) &
            * (VarEast - VarWest)

       DiffVar(:,iLat) = DiffVar(:,iLat) + &
            (DiffLonFlux(2:nLons+1) - DiffLonFlux(1:nLons)) &
            / dLonDist_GB(iLat, iAlt, iBlock)
    end do

  end subroutine calc_rusanov_horizontal

end subroutine advance_horizontal
