!^CFG COPYRIGHT UM

module ModUser

  use ModSize,     ONLY: nI,nJ,nK,gcn,nBLK

  use ModUserEmpty, &
       IMPLEMENTED1 => user_calc_sources

  include 'user_module.h' !list of public methods

  real, parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = &
       'Non-Gyrotropic Reconnection Model, Masha Kuznetsova'


  real, dimension(1:nI,1:nJ,1:nK):: &
       Srho,SrhoUx,SrhoUy,SrhoUz,SBx,SBy,SBz,Sp,SEr,SE
  integer, parameter :: nii=81, njj=170
  !			       njj=145
  !                              njj=120
  !                              njj=192
  !
  ! taken from small scale runs
  real, parameter :: alpha = 0.156250
  !
  !!         real, parameter :: alpha = 0.036
  integer, dimension(1:nii,1:njj) :: iBLK_rec, iPE_rec
  integer, dimension(1:njj) :: iBLK_rec1, iBLK_rec2, iPE_rec1, iPE_rec2

  real, dimension(1:nii) :: XXR
  real, dimension(1:njj) :: YYR, XX0, DWZ_0, DWX_0, EF0_0, &
       BZPR_0, VXPR_0, Cos_Thet_0, Sin_Thet_0
  real, dimension(1:njj) :: XX1, DWZ_1, DWX_1, EF0_1, &
       BZPR_1, VXPR_1, Cos_Thet_1, Sin_Thet_1


contains

  subroutine user_calc_sources
    use ModMain
    use ModVarIndexes
    use ModAdvance
    use ModGeometry
    use ModConst
    use ModNumConst
    use ModPhysics
    use ModProcMH
    use ModPhysics
    use ModParallel

    ! User declared local variables go here
    integer :: i,j,k, jj
    real :: xt, yt, zt
    real :: xxx, zzz, yyy, xxx1
    real :: DWZ0
    real :: delta_tau
    ! DWZ for rec2 test     DWZ = 1
    !       DWZ = 0.25
    DWZ0 = alpha*sqrt(2.)*2.

    do k = 1, nK
       do j = 1, nJ
          do i = 1, nI
             ! User coding of the source terms should be inside this loop
             ! Note that user source terms should be added to the the
             ! values already loaded in the arrays.  For example
             !     Srho(i,j,k) = Srho(i,j,k) + user source term for rho
             xt = x_BLK(i,j,k,globalBLK)
             yt = y_BLK(i,j,k,globalBLK)
             zt = z_BLK(i,j,k,globalBLK)
             zzz = 2.* zt/DWZ0
             yyy = yt/25.
             !       jj = 1 + floor((abs(yt) - 0.125)/0.25)
             !!        jj = 1 + floor((abs(yt) - 0.5*dy_BLK(globalBLK))/dy_BLK(globalBLK))
             jj = 1 + floor((abs(yt) - 0.0625)/0.125)
             jj = max(1,jj)
             !!        if (jj <= floor(18./dy_BLK(globalBLK))-1) then
             if (jj <= 143) then
                if(xt > -99. .and. xt < 0. .and. XX0(jj) > -86. .and. XX0(jj) < -3. &
                     .and. DWX_0(jj) > 0. .and. EF0_0(jj) > 0.) then
                   xxx = 2.*(xt - XX0(jj))/DWX_0(jj)
                   !!        zzz = 2.*zt/DWZ(jj)


                   !!	Source_VC(Bx_,i,j,k) = Source_VC(Bx_,i,j,k) &
                   !!        - EF0_0(jj)*Cos_Thet_0(jj)/DWZ0*tanh(zzz)/cosh(zzz)/cosh(xxx)/cosh(yyy)
                   !
                   Source_VC(Bx_,i,j,k) = Source_VC(Bx_,i,j,k) &
                        - 2.*EF0_0(jj)/DWZ0*tanh(zzz)/cosh(zzz)/cosh(xxx)/cosh(yyy)


                   !! 
                   !!	Source_VC(Bx_,i,j,k) = Source_VC(Bx_,i,j,k) &
                   !!       + EF0_0(jj)*Sin_Thet_0(jj)/DWZ0*tanh(zzz)/cosh(zzz)/cosh(xxx)/cosh(yyy)
                   !
!!!        Source_VC(Bx_,i,j,k) = Source_VC(Bx_,i,j,k) &
!!!        - EF0_0(jj)/DWZ(jj)*tanh(zzz)/cosh(zzz)/cosh(xxx)/cosh(yyy)
                   !
                   !
                   Source_VC(Bz_,i,j,k) = Source_VC(Bz_,i,j,k) + &
                        2.* EF0_0(jj)/DWX_0(jj)*tanh(xxx)/cosh(zzz)/cosh(xxx)/cosh(yyy)
                   !
                endif
                !
                !! Add term from second reconnection site
                if(xt > -99. .and. xt < 0. .and. xt < XX0(jj) .and. XX1(jj) > -86. &
                     .and. DWX_1(jj) > 0.  .and. XX1(jj) < -5. .and. DWX_0(jj) > 0. &
                     .and. XX1(jj) < XX0(jj)-DWX_0(jj)*0.5  .and. EF0_1(jj) > 0.) then
                   xxx1 = 2.*(xt - XX1(jj))/DWX_1(jj)
                   !!        zzz = 2.*zt/DWZ(jj)


                   !!	Source_VC(Bx_,i,j,k) = Source_VC(Bx_,i,j,k) &
                   !!        - EF0_1(jj)*Cos_Thet_1(jj)/DWZ0*tanh(zzz)/cosh(zzz)/cosh(xxx1)/cosh(yyy)
                   !
                   Source_VC(Bx_,i,j,k) = Source_VC(Bx_,i,j,k) &
                        - 2.*EF0_1(jj)/DWZ0*tanh(zzz)/cosh(zzz)/cosh(xxx1)/cosh(yyy)


                   !! 
                   !!	Source_VC(Bx_,i,j,k) = Source_VC(Bx_,i,j,k) &
                   !!       + EF0_1(jj)*Sin_Thet_1(jj)/DWZ0*tanh(zzz)/cosh(zzz)/cosh(xxx1)/cosh(yyy)
                   !
!!!        Source_VC(Bx_,i,j,k) = Source_VC(Bx_,i,j,k) &
!!!        - EF0_1(jj)/DWZ(jj)*tanh(zzz)/cosh(zzz)/cosh(xxx1)/cosh(yyy)
                   !
                   !
                   Source_VC(Bz_,i,j,k) = Source_VC(Bz_,i,j,k) + &
                        2.* EF0_1(jj)/DWX_1(jj)*tanh(xxx1)/cosh(zzz)/cosh(xxx1)/cosh(yyy)
                   !
                endif


             endif
          end do     ! end the i loop
       end do     ! end the j loop
    end do     ! end the k loop

  end subroutine user_calc_sources

  !========================================================================
  subroutine set_rec_line
    use ModMain
    use ModVarIndexes
    use ModSize
    use ModIO
    use ModPhysics
    use ModGeometry
    use ModAdvance
    use ModParallel
    use ModProcMH
    use ModConst
    use ModMpi

    implicit none

    integer :: isat, iPE,iPEtmp, iPEmax, iBLK, iBLKtemp, iBLKmax
    integer :: iPE1, iPE2, iPEtmp1, iPEtmp2, iBLKtemp1, iBLKtemp2

    integer :: i,j,k
    logical :: oktest, oktest_me
    integer :: ii,jj
    integer :: j1, i1, iX, jY
    integer :: iError
    real :: XX0tmp, XX0max, XX0max_all,XX0tmp_all, yj1
    real :: Rho1, Rho2, Rho3, Rho4, Rho0, P0, BXPR, DWX0, RhoU1, RhoU2
    real :: Jxx, Jyy, Jxy, CosThet, SinThet
    real :: VXPRX, VXPRY, VYPRX, VYPRY
    real, dimension(1:3) :: Xrec1, Xrec2
    real,  dimension(1:7) :: qval
    real, external :: point_value1
    real,dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn):: &
         BX,BY,BZ,B2
    !



    !--------------------------------------------------------------------------


    do jj=1, njj
       !   YYR(jj)=0.03125+0.0625*(jj -1)
       YYR(jj)=0.0625+0.125*(jj-1)
       XX0(jj) = -100.
       XX0tmp = -200.
       XX0max = -100.
       XX0tmp_all = -100.
       XX0max_all = -100.
       !       BXPR = -1.
       BZPR_0(jj) = -1.
       VXPR_0(jj) = 0.
       DWX_0(jj) = -1.
       DWZ_0(jj) = -1.
       EF0_0(jj) = 0.
       Cos_Thet_0(jj) = 1.
       Sin_Thet_0(jj) = 0.
       qval(1) = -1.
       qval(2) = 0.
       qval(3) = -1.
       qval(4) = 0.
       qval(5) = -1.
       qval(6) = 1.
       qval(7) = 0.
       !
       !
       !
       iPE = -1
       iPE1 = -1
       iPE2 = -1
       iBLKtemp = -1
       iBLKtemp1 = -1
       iBLKtemp2 = -1
       iPEmax = -1
       iBLKmax = -1

       do iBLK = 1,nBlockMax



          if (unusedBLK(iBLK)) cycle


          if (x_BLK(nI,1,1,iBLK) > -87. .and. &
               x_BLK(nI,1,1,iBLK) < -3.  .and.&
               YYR(jj) >= y_BLK(1,1,1,iBLK)       - 0.50*dy_BLK(iBLK) .and. &
               YYR(jj) < y_BLK(1,nJ,1,iBLK) + 0.50*dy_BLK(iBLK) .and. &
               0.25*dz_BLK(iBLK) >=  z_BLK(1,1,1,iBLK)       - 0.50*dz_BLK(iBLK) .and. &
               0.25*dz_BLK(iBLK) < z_BLK(1,1,nK,iBLK) + 0.50*dz_BLK(iBLK) ) then

             yj1=(YYR(jj)-y_BLK(1,1,1,iBLK))/dy_BLK(iBLK)+1.
             j1=floor(yj1)


             !       write(*,*)'me_world=',iproc,'iBLK=',iBLK,'YYR=',yyR(jj),'j1=',j1,'yj1=',yj1
             do i1=1,nI
                i = nI+1-i1
                !
                if (B0zCell_BLK(i+1,j1,1,iBLK)+State_VGB(Bz_,i+1,j1,1,iBLK) > 0. .and. &
                     B0zCell_BLK(i-1,j1,1,iBLK)+State_VGB(Bz_,i-1,j1,1,iBLK) <=0.) then


                   XX0tmp = x_BLK(i,j1,1,iBLK)
                   if (XX0tmp > XX0max) then
                      XX0max=XX0tmp
                      iX=i
                      iBLKmax = iBLK
                      jY = j1
                   endif
                   !       write(*,*)'me_world=',iproc,'iBLK=',iBLK,'XX0tmp=',XX0tmp
                   goto 20
                endif
             enddo

          end if
20     enddo

       call MPI_ALLREDUCE(XX0max, XX0max_all, 1, &
            MPI_REAL, MPI_MAX, iComm,iError)
       XX0(jj) = XX0max_all
       !
       if (XX0max == XX0max_all .and. XX0max_all > -99.) then
          iPE = iproc
          !       write(*,*),jj,'Y=',YYR(jj),'XX0=',XX0(jj),'iBLKmax=',iBLKmax,'iX=',iX,'jY=',jY
          !       endif
          !
          BX(:,:,:) = State_VGB(Bx_,:,:,:,iBLKmax)+B0xCell_BLK(:,:,:,iBLKmax)
          BY(:,:,:) = State_VGB(By_,:,:,:,iBLKmax)+B0yCell_BLK(:,:,:,iBLKmax)
          BZ(:,:,:) = State_VGB(Bz_,:,:,:,iBLKmax)+B0zCell_BLK(:,:,:,iBLKmax)
          !        B2(:,:,:) = BX**2+BY**2+BZ**2
          Jyy=0.5*((BX(iX,jY,2)-BX(iX,jY,0))/dz_BLK(iBLKmax) - &
               (BZ(iX+1,jY,1)-BZ(iX-1,jY,1))/dx_BLK(iBLKmax))
          Jxx=0.5*((BZ(iX,jY+1,1)-BZ(iX,jY-1,1))/dy_BLK(iBLKmax)     - &
               (BY(iX,jY,2)-BY(iX,jY,0))/dz_BLK(iBLKmax))
          Jxy=sqrt(Jyy**2+Jxx**2)
          CosThet=Jyy/Jxy
          SinThet=Jxx/Jxy
          !!        CosThet=1.
          !!        SinThet=0.
          qval(6)=CosThet
          qval(7)=SinThet
          !!        qval(1) =(BZ(iX+1,jY,1)-BZ(iX-1,jY,1))/(2.*dx_BLK(iBLKmax))*CosThet- &
          !!        (BZ(iX,jY+1,1)-BZ(iX,jY-1,1))/(2.*dy_BLK(iBLKmax))*SinThet
          !!        write(*,*)'CosThet=',CosThet,'SinThet=',SinThet, &
          !!         'BZPR=',(BZ(iX,jY+1,1)-BZ(iX,jY-1,1))/(2.*dy_BLK(iBLKmax))
          !!       write(*,*)'Jyy=',Jyy,'Jxx=',Jxx
          qval(1) =(BZ(iX+1,jY,1)-BZ(iX-1,jY,1))/(2.*dx_BLK(iBLKmax))
          !     
          !
          !        
          Rho1 = State_VGB(rho_,iX-1,jY,1,iBLKmax)
          Rho2 = State_VGB(rho_,iX+1,jY,1,iBLKmax)
          Rho0 = State_VGB(rho_,iX,jY,1,iBLKmax)
          Rho3 = State_VGB(rho_,iX,jY-1,1,iBLKmax)
          Rho4 = State_VGB(rho_,iX,jY+1,1,iBLKmax)
          P0 = State_VGB(P_,iX,jY,1,iBLKmax)
          VXPRX=(State_VGB(rhoUx_,iX+1,jY,1,iBLKmax)/Rho2 - &
               State_VGB(rhoUx_,iX-1,jY,1,iBLKmax)/Rho1)/(2.*dx_BLK(iBLKmax)) 
          VYPRX=(State_VGB(rhoUy_,iX+1,jY,1,iBLKmax)/Rho2 - &
               State_VGB(rhoUy_,iX-1,jY,1,iBLKmax)/Rho1)/(2.*dx_BLK(iBLKmax)) 
          VXPRY=(State_VGB(rhoUx_,iX,jY+1,1,iBLKmax)/Rho4 - &
               State_VGB(rhoUx_,iX,jY-1,1,iBLKmax)/Rho3)/(2.*dy_BLK(iBLKmax))
          VYPRY=(State_VGB(rhoUy_,iX,jY+1,1,iBLKmax)/Rho4 - &
               State_VGB(rhoUy_,iX,jY-1,1,iBLKmax)/Rho3)/(2.*dy_BLK(iBLKmax))
          !!        qval(2)=VXPRX*CosThet**2-CosThet*SinThet*(VYPRX+VXPRY)+ &
          !!               VYPRY*SinThet**2
          qval(2) = VXPRX
          !!        write(*,*)'VXPRX=',VXPRX,'VYPRX=',VYPRY,'VXPRY=',VXPRY,'VYPRY=',VYPRY


          BXPR =((B0xCell_BLK(iX,jY,2,iBLKmax)+State_VGB(Bx_,iX,jY,2,iBLKmax)) &
               -  (B0xCell_BLK(iX,jY,0,iBLKmax)+State_VGB(Bx_,iX,jY,0,iBLKmax)))/ &
               (2.* dz_BLK(iBLKmax))
          !
          if (qval(1) > 0) qval(3) = sqrt(alpha/qval(1)*sqrt(2.*P0/Rho0))
          if (BXPR > 0) then
             qval(5) = max(dz_BLK(iBLKmax)*4.,sqrt(alpha/BXPR*sqrt(2.*P0/Rho0)))
             !       qval(5) = min(qval(5),alpha*sqrt(2.)*2.)
          else
             qval(5) = dz_BLK(iBLKmax)*4.
          endif
          qval(4) = alpha*sqrt(2.*P0/Rho0)
          !       write(*,*)jj,' Y=',YYR(jj),'iPemax=',iproc,'XX0=', XX0max_all, 'BZPR =',qval(1), &
          !       'VXPR=',qval(2), 'DWX = ', qval(3), 'EF0=',qval(4)
       else
          iPE = -1
       endif
       ! Tell everyone which processor contains iBLKmax
       call MPI_ALLREDUCE(iPE, iPEmax, 1, &
            MPI_INTEGER, MPI_MAX, iComm,iError)
       !
       DWX0 = qval(3)
       !!        if (iPEmax >= 0) then
       !!        call MPI_Bcast(DWX0,1,MPI_REAL,iPEmax,iComm,iError)
       !!        DWX_0(jj) = DWX0
       !!           endif
       !
       !
       !
       if (iPEmax >= 0) then
          call MPI_Bcast(qval,7,MPI_REAL,iPEmax,iComm,iError)
          BZPR_0(jj) = qval(1)
          VXPR_0(jj) = qval(2)
          DWX_0(jj)= qval(3)
          !   EF0 was increased on factor 2 for run hr4
          EF0_0(jj) = qval(4)*qval(2)
          DWZ_0(jj) = qval(5)
          Cos_Thet_0(jj) = qval(6)
          Sin_Thet_0(jj) = qval(7)
       endif


       !
!!!     if (iproc == 0) then
!!!     write(*,*)jj,' Y=',YYR(jj),'iPEmax=',iPEmax,'XX0=', XX0(jj), 'BZPR =',BZPR(jj), &
!!!       'VXPR=',BZPR(jj), 'DWX = ', DWX(jj), 'EF0=', EF0(jj)
!!!     endif
    enddo
    if (iproc == 0) then
       write(*,*)'XX0=', XX0(1), 'BZPR =',BZPR_0(1), &
            'VXPR=',VXPR_0(1), 'DWX = ', DWX_0(1), 'DWZ =',DWZ_0(1), 'EF0=', EF0_0(1), &
            'CosThet=',Cos_Thet_0(1), 'SinThet=', Sin_Thet_0(1)
       !!
       write(*,*)'XX0_48=', XX0(48), 'BZPR =',BZPR_0(48), &
            'VXPR=',VXPR_0(48), 'DWX = ', DWX_0(48), 'DWZ =',DWZ_0(48), 'EF0=', EF0_0(48), &
            'CosThet=',Cos_Thet_0(48), 'SinThet=', Sin_Thet_0(48)

       !!
       write(*,*)'XX0_96=', XX0(96), 'BZPR =',BZPR_0(96), &
            'VXPR=',VXPR_0(96), 'DWX = ', DWX_0(96), 'DWZ =',DWZ_0(96), 'EF0=', EF0_0(96), &
            'CosThet=',Cos_Thet_0(96), 'SinThet=', Sin_Thet_0(96)

    endif

!!!!!!!!


  end subroutine set_rec_line

  !======================================================================

  subroutine set_rec_line1
    use ModMain
    use ModVarIndexes
    use ModSize
    use ModIO
    use ModPhysics
    use ModGeometry
    use ModAdvance
    use ModParallel
    use ModProcMH
    use ModConst
    use ModMpi

    implicit none

    integer :: isat, iPE,iPEtmp, iPEmax, iBLK, iBLKtemp, iBLKmax
    integer :: iPE1, iPE2, iPEtmp1, iPEtmp2, iBLKtemp1, iBLKtemp2

    integer :: i,j,k
    logical :: oktest, oktest_me
    integer :: ii,jj
    integer :: j1, i1, iX, jY
    integer :: iError
    real :: XX0tmp, XX0max, XX0max_all,XX0tmp_all, yj1
    real :: Rho1, Rho2, Rho3, Rho4, Rho0, P0, BXPR, DWX0, RhoU1, RhoU2
    real :: Jxx, Jyy, Jxy, CosThet, SinThet
    real :: VXPRX, VXPRY, VYPRX, VYPRY
    real, dimension(1:3) :: Xrec1, Xrec2
    real,  dimension(1:7) :: qval
    real, external :: point_value1
    real,dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn):: &
         BX,BY,BZ,B2
    !



    !--------------------------------------------------------------------------


    do jj=1, njj
       !   YYR(jj)=0.03125+0.0625*(jj -1)
       YYR(jj)=0.0625+0.125*(jj-1)
       XX1(jj) = -100.
       XX0tmp = -200.
       XX0max = -100.
       XX0tmp_all = -100.
       XX0max_all = -100.
       !       BXPR = -1.
       BZPR_1(jj) = -1.
       VXPR_1(jj) = 0.
       DWX_1(jj) = -1.
       DWZ_1(jj) = -1.
       EF0_1(jj) = 0.
       Cos_Thet_1(jj) = 1.
       Sin_Thet_1(jj) = 0.
       qval(1) = -1.
       qval(2) = 0.
       qval(3) = -1.
       qval(4) = 0.
       qval(5) = -1.
       qval(6) = 1.
       qval(7) = 0.
       !
       !
       !
       iPE = -1
       iPE1 = -1
       iPE2 = -1
       iBLKtemp = -1
       iBLKtemp1 = -1
       iBLKtemp2 = -1
       iPEmax = -1
       iBLKmax = -1

       do iBLK = 1,nBlockMax



          if (unusedBLK(iBLK)) cycle


          if (x_BLK(nI,1,1,iBLK) > -87. .and. &
               x_BLK(nI,1,1,iBLK) < -5.  .and.  x_BLK(nI,1,1,iBLK) < XX0(jj) - 2. .and. &
               YYR(jj) >= y_BLK(1,1,1,iBLK)       - 0.50*dy_BLK(iBLK) .and. &
               YYR(jj) < y_BLK(1,nJ,1,iBLK) + 0.50*dy_BLK(iBLK) .and. &
               0.25*dz_BLK(iBLK) >=  z_BLK(1,1,1,iBLK)       - 0.50*dz_BLK(iBLK) .and. &
               0.25*dz_BLK(iBLK) < z_BLK(1,1,nK,iBLK) + 0.50*dz_BLK(iBLK) ) then

             yj1=(YYR(jj)-y_BLK(1,1,1,iBLK))/dy_BLK(iBLK)+1.
             j1=floor(yj1)


             !       write(*,*)'me_world=',iproc,'iBLK=',iBLK,'YYR=',yyR(jj),'j1=',j1,'yj1=',yj1
             do i1=1,nI
                i = nI+1-i1
                !
                if (B0zCell_BLK(i+1,j1,1,iBLK)+State_VGB(Bz_,i+1,j1,1,iBLK) > 0. .and. &
                     B0zCell_BLK(i-1,j1,1,iBLK)+State_VGB(Bz_,i-1,j1,1,iBLK) <=0.) then


                   XX0tmp = x_BLK(i,j1,1,iBLK)
                   if (XX0tmp > XX0max) then
                      XX0max=XX0tmp
                      iX=i
                      iBLKmax = iBLK
                      jY = j1
                   endif
                   !       write(*,*)'me_world=',iproc,'iBLK=',iBLK,'XX0tmp=',XX0tmp
                   goto 20
                endif
             enddo

          end if
20     enddo

       call MPI_ALLREDUCE(XX0max, XX0max_all, 1, &
            MPI_REAL, MPI_MAX, iComm,iError)
       XX1(jj) = XX0max_all
       !
       if (XX0max == XX0max_all .and. XX0max_all > -99.) then
          iPE = iproc
          !       write(*,*),jj,'Y=',YYR(jj),'XX0=',XX0(jj),'iBLKmax=',iBLKmax,'iX=',iX,'jY=',jY
          !       endif
          !
          BX(:,:,:) = State_VGB(Bx_,:,:,:,iBLKmax)+B0xCell_BLK(:,:,:,iBLKmax)
          BY(:,:,:) = State_VGB(By_,:,:,:,iBLKmax)+B0yCell_BLK(:,:,:,iBLKmax)
          BZ(:,:,:) = State_VGB(Bz_,:,:,:,iBLKmax)+B0zCell_BLK(:,:,:,iBLKmax)
          !        B2(:,:,:) = BX**2+BY**2+BZ**2
          Jyy=0.5*((BX(iX,jY,2)-BX(iX,jY,0))/dz_BLK(iBLKmax) - &
               (BZ(iX+1,jY,1)-BZ(iX-1,jY,1))/dx_BLK(iBLKmax))
          Jxx=0.5*((BZ(iX,jY+1,1)-BZ(iX,jY-1,1))/dy_BLK(iBLKmax)     - &
               (BY(iX,jY,2)-BY(iX,jY,0))/dz_BLK(iBLKmax))
          Jxy=sqrt(Jyy**2+Jxx**2)
          CosThet=Jyy/Jxy
          SinThet=Jxx/Jxy
          !!        CosThet=1.
          !!        SinThet=0.
          qval(6)=CosThet
          qval(7)=SinThet
          !!        qval(1) =(BZ(iX+1,jY,1)-BZ(iX-1,jY,1))/(2.*dx_BLK(iBLKmax))*CosThet- &
          !!        (BZ(iX,jY+1,1)-BZ(iX,jY-1,1))/(2.*dy_BLK(iBLKmax))*SinThet
          !!        write(*,*)'CosThet=',CosThet,'SinThet=',SinThet, &
          !!         'BZPR=',(BZ(iX,jY+1,1)-BZ(iX,jY-1,1))/(2.*dy_BLK(iBLKmax))
          !!       write(*,*)'Jyy=',Jyy,'Jxx=',Jxx
          qval(1) =(BZ(iX+1,jY,1)-BZ(iX-1,jY,1))/(2.*dx_BLK(iBLKmax))
          !     
          !
          !        
          Rho1 = State_VGB(rho_,iX-1,jY,1,iBLKmax)
          Rho2 = State_VGB(rho_,iX+1,jY,1,iBLKmax)
          Rho0 = State_VGB(rho_,iX,jY,1,iBLKmax)
          Rho3 = State_VGB(rho_,iX,jY-1,1,iBLKmax)
          Rho4 = State_VGB(rho_,iX,jY+1,1,iBLKmax)
          P0 = State_VGB(P_,iX,jY,1,iBLKmax)
          VXPRX=(State_VGB(rhoUx_,iX+1,jY,1,iBLKmax)/Rho2 - &
               State_VGB(rhoUx_,iX-1,jY,1,iBLKmax)/Rho1)/(2.*dx_BLK(iBLKmax)) 
          VYPRX=(State_VGB(rhoUy_,iX+1,jY,1,iBLKmax)/Rho2 - &
               State_VGB(rhoUy_,iX-1,jY,1,iBLKmax)/Rho1)/(2.*dx_BLK(iBLKmax)) 
          VXPRY=(State_VGB(rhoUx_,iX,jY+1,1,iBLKmax)/Rho4 - &
               State_VGB(rhoUx_,iX,jY-1,1,iBLKmax)/Rho3)/(2.*dy_BLK(iBLKmax))
          VYPRY=(State_VGB(rhoUy_,iX,jY+1,1,iBLKmax)/Rho4 - &
               State_VGB(rhoUy_,iX,jY-1,1,iBLKmax)/Rho3)/(2.*dy_BLK(iBLKmax))
          !!        qval(2)=VXPRX*CosThet**2-CosThet*SinThet*(VYPRX+VXPRY)+ &
          !!               VYPRY*SinThet**2
          qval(2) = VXPRX
          !!        write(*,*)'VXPRX=',VXPRX,'VYPRX=',VYPRY,'VXPRY=',VXPRY,'VYPRY=',VYPRY


          BXPR =((B0xCell_BLK(iX,jY,2,iBLKmax)+State_VGB(Bx_,iX,jY,2,iBLKmax)) &
               -  (B0xCell_BLK(iX,jY,0,iBLKmax)+State_VGB(Bx_,iX,jY,0,iBLKmax)))/ &
               (2.* dz_BLK(iBLKmax))
          !
          if (qval(1) > 0) qval(3) = sqrt(alpha/qval(1)*sqrt(2.*P0/Rho0))
          if (BXPR > 0) then
             qval(5) = max(dz_BLK(iBLKmax)*4.,sqrt(alpha/BXPR*sqrt(2.*P0/Rho0)))
             !       qval(5) = min(qval(5),alpha*sqrt(2.)*2.)
          else
             qval(5) = dz_BLK(iBLKmax)*4.
          endif
          qval(4) = alpha*sqrt(2.*P0/Rho0)
          !       write(*,*)jj,' Y=',YYR(jj),'iPemax=',iproc,'XX0=', XX0max_all, 'BZPR =',qval(1), &
          !       'VXPR=',qval(2), 'DWX = ', qval(3), 'EF0=',qval(4)
       else
          iPE = -1
       endif
       ! Tell everyone which processor contains iBLKmax
       call MPI_ALLREDUCE(iPE, iPEmax, 1, &
            MPI_INTEGER, MPI_MAX, iComm,iError)
       !
       DWX0 = qval(3)
       !!        if (iPEmax >= 0) then
       !!        call MPI_Bcast(DWX0,1,MPI_REAL,iPEmax,iComm,iError)
       !!        DWX_0(jj) = DWX0
       !!           endif
       !
       !
       !
       if (iPEmax >= 0) then
          call MPI_Bcast(qval,7,MPI_REAL,iPEmax,iComm,iError)
          BZPR_1(jj) = qval(1)
          VXPR_1(jj) = qval(2)
          DWX_1(jj)= qval(3)
          !   EF0 was increased on factor 2 for run hr4
          EF0_1(jj) = qval(4)*qval(2)
          DWZ_1(jj) = qval(5)
          Cos_Thet_1(jj) = qval(6)
          Sin_Thet_1(jj) = qval(7)
       endif


       !
!!!     if (iproc == 0) then
!!!     write(*,*)jj,' Y=',YYR(jj),'iPEmax=',iPEmax,'XX0=', XX0(jj), 'BZPR =',BZPR(jj), &
!!!       'VXPR=',BZPR(jj), 'DWX = ', DWX(jj), 'EF0=', EF0(jj)
!!!     endif
    enddo
    if (iproc == 0) then
       write(*,*)'XX1=', XX1(1), 'BZPR =',BZPR_1(1), &
            'VXPR=',VXPR_1(1), 'DWX = ', DWX_1(1), 'DWZ =',DWZ_1(1), 'EF0=', EF0_1(1), &
            'CosThet=',Cos_Thet_1(1), 'SinThet=', Sin_Thet_1(1)
       !
       write(*,*)'XX1_48=', XX1(48), 'BZPR =',BZPR_1(48), &
            'VXPR=',VXPR_1(48), 'DWX = ', DWX_1(48), 'DWZ =',DWZ_1(48), 'EF0=', EF0_1(48), &
            'CosThet=',Cos_Thet_1(48), 'SinThet=', Sin_Thet_1(48)
       !
       write(*,*)'XX1_96=', XX1(96), 'BZPR =',BZPR_1(96), &
            'VXPR=',VXPR_1(96), 'DWX = ', DWX_1(96), 'DWZ =',DWZ_1(96), 'EF0=', EF0_1(96), &
            'CosThet=',Cos_Thet_1(96), 'SinThet=', Sin_Thet_1(96)

    endif

    !


  end subroutine set_rec_line1


end module ModUser

