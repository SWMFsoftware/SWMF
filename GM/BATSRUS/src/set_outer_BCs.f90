!^CFG COPYRIGHT UM

subroutine set_outer_BCs(iBlock, time_now, DoSetEnergy)

  ! Set ghost cells values rho, U, B, and P for iBLK. 
  ! Set E if DoSetEnergy is true.
  use ModSetOuterBC
  use ModProcMH
  use ModMain
  use ModVarIndexes
  use ModAdvance, ONLY: State_VGB, Erad_
  use ModParallel, ONLY: NOBLK, NeiLev
  use ModGeometry, ONLY: far_field_BCs_BLK, MaxBoundary, TypeGeometry, XyzMin_D
  use ModPhysics
  use ModUser, ONLY: user_set_outerBCs
  use ModMultiFluid, ONLY: iFluid, nFluid, iRhoUx_I, iRhoUy_I, iRhoUz_I
  use ModEnergy, ONLY: calc_energy
  use ModGrayDiffusion, ONLY: set_gray_outflow_bc !^CFG IF IMPLICIT
  use ModMultiGroupDiffusion, ONLY: set_rad_outflow_bc !^CFG IF IMPLICIT
  implicit none

  integer, intent(in) :: iBlock
  real,    intent(in) :: time_now
  logical, intent(in) :: DoSetEnergy

  integer :: iStart, iVar, iSide

  integer :: iGhost, jGhost, kGhost, iGhost2, jGhost2, kGhost2

  logical :: DoTest, DoTestMe, IsFound

  character(len=*), parameter :: NameSub = 'set_outer_bcs'
  !--------------------------------------------------------------------------
  iBLK=iBlock
  if(iBLK==BLKtest.and.iProc==PROCtest)then
     call set_oktest(NameSub, DoTest, DoTestMe)
  else
     DoTest=.false.; DoTestMe=.false.
  endif

  if(.not.far_field_BCs_BLK(iBLK))then
     write(*,*) NameSub,' warning: iBLK=',iBLK,' is not far_field block'
     RETURN
  end if

  if(DoTestMe)write(*,*)NameSub,': iBLK, set_E, neiLEV=',&
       iBLK,DoSetEnergy,neiLEV(:,iBLK)

  if(DoTestMe)then
     Ighost =Itest; Jghost=Jtest;  Kghost=Ktest
     Ighost2=Itest; Jghost2=Jtest; Kghost2=Ktest
     select case(DimTest)
     case(x_)
        if(iTest== 1)then; iGhost=0;    iGhost2 = -1;   endif
        if(iTest==nI)then; iGhost=nI+1; iGhost2 = nI+2; endif
     case(y_)
        if(jTest== 1)then; jGhost=0;    jGhost2 = -1;   endif
        if(jTest==nJ)then; jGhost=nJ+1; jGhost2 = nJ+2; endif
     case(z_)
        if(kTest== 1)then; kGhost=0;    kGhost2 = -1;   endif
        if(kTest==nK)then; kGhost=nK+1; kGhost2 = nK+2; endif
     end select

     write(*,*)'iTest,  jTest,  kTest   =',iTest,  jTest,  kTest
     write(*,*)'iGhost, jGhost, kGhost  =',iGhost, jGhost, kGhost
     write(*,*)'iGhost2,jGhost2,kGhost2 =',iGhost2,jGhost2,kGhost2
     do iVar=1,nVar
        write(*,*)'initial',NameVar_V(iVar),   'cell,ghost,ghost2=',&
             State_VGB(iVar,Itest,Jtest,Ktest,iBLK),&
             State_VGB(iVar,Ighost,Jghost,Kghost,iBLK), &
             State_VGB(iVar,Ighost2,Jghost2,Kghost2,iBLK)
     end do
  end if

  iStart=max(MaxBoundary+1,1)
  do iSide = iStart, Top_

     ! Check if this side of the block is indeed an outer boundary
     if(neiLEV(iSide,iBLK)/=NOBLK) CYCLE

     ! Do not apply cell boundary conditions at the pole 
     ! This is either handled by message passing or supercell
     if(TypeGeometry(1:9) == 'spherical' .and. iSide >= Bot_) CYCLE
     if(TypeGeometry      == 'cylindrical' .and. &
          XyzMin_D(1) == 0.0 .and. iSide == 1) CYCLE

     ! Set index limits
     imin1g=-1; imax1g=nI+2; imin2g=-1; imax2g=nI+2
     jmin1g=-1; jmax1g=nJ+2; jmin2g=-1; jmax2g=nJ+2
     kmin1g=-1; kmax1g=nK+2; kmin2g=-1; kmax2g=nK+2

     imin1p=-1; imax1p=nI+2; imin2p=-1; imax2p=nI+2
     jmin1p=-1; jmax1p=nJ+2; jmin2p=-1; jmax2p=nJ+2
     kmin1p=-1; kmax1p=nK+2; kmin2p=-1; kmax2p=nK+2

     select case(iSide)
     case(east_)
        imin1g=0; imax1g=0; imin2g=-1; imax2g=-1
        imin1p=1; imax1p=1; imin2p= 2; imax2p= 2
     case(west_)
        imin1g=nI+1; imax1g=nI+1; imin2g=nI+2; imax2g=nI+2
        imin1p=nI  ; imax1p=nI  ; imin2p=nI-1; imax2p=nI-1
     case(south_)
        jmin1g=0; jmax1g=0; jmin2g=-1; jmax2g=-1
        jmin1p=1; jmax1p=1; jmin2p= 2; jmax2p= 2
     case(north_)
        jmin1g=nJ+1; jmax1g=nJ+1; jmin2g=nJ+2; jmax2g=nJ+2
        jmin1p=nJ  ; jmax1p=nJ  ; jmin2p=nJ-1; jmax2p=nJ-1
     case(bot_)
        kmin1g=0; kmax1g=0; kmin2g=-1; kmax2g=-1
        kmin1p=1; kmax1p=1; kmin2p= 2; kmax2p= 2
     case(top_)
        kmin1g=nK+1; kmax1g=nK+1; kmin2g=nK+2; kmax2g=nK+2
        kmin1p=nK  ; kmax1p=nK  ; kmin2p=nK-1; kmax2p=nK-1
     end select

     select case(TypeBc_I(iSide))
     case('coupled')
        ! For SC-IH coupling the extra wave energy variable needs a BC
        if(NameThisComp == 'SC') call BC_cont(ScalarFirst_, ScalarLast_)
     case('periodic')
        call stop_mpi('The neighbors are not deifned at the periodic boundary')
     case('float','outflow')
        call BC_cont(1,nVar)
        !^CFG IF IMPLICIT BEGIN
        if(UseGrayDiffusion.or.UseRadDiffusion) &
             call set_radiation_outflow_bc(WaveFirst_, WaveLast_, iSide)
        !^CFG END IMPLICIT
     case('raeder')
        call BC_cont(1,nVar)
        if(iSide==north_.or.iSide==south_)then
           call BC_fixed(By_,By_,DefaultState_V)
        elseif(iSide==bot_.or.iSide==top_)then
           call BC_fixed(Bz_,Bz_,DefaultState_V)
        end if
     case('reflect')
        ! Scalars are symmetric
        call BC_symm(1,nVar)
        ! Normal vector components are mirror symmetric
        if(iSide==east_.or.iSide==west_)then
           do iFluid = 1, nFluid
              call BC_asymm(iRhoUx_I(iFluid), iRhoUx_I(iFluid))
           end do
           if(UseB)call BC_asymm(Bx_,Bx_)
        endif
        if(iSide==south_.or.iSide==north_)then
           do iFluid = 1, nFluid
              call BC_asymm(iRhoUy_I(iFluid), iRhoUy_I(iFluid))
           end do
           if(UseB)call BC_asymm(By_,By_)
        endif
        if(iSide==bot_.or.iSide==top_)then
           do iFluid = 1, nFluid
              call BC_asymm(iRhoUz_I(iFluid), iRhoUz_I(iFluid))
           end do
           if(UseB)call BC_asymm(Bz_,Bz_)
        endif
     case('linetied')
        call BC_symm(rho_,rho_)
        call BC_asymm(rhoUx_,rhoUz_)
        call BC_cont(rhoUz_+1,nVar)
     case('fixed','inflow','vary','ihbuffer')
        if(time_accurate &
             .and. (TypeBc_I(iSide)=='vary'.or.TypeBc_I(iSide)=='inflow'))then
           call BC_solar_wind(time_now)
        else if(TypeBc_I(iSide)=='ihbuffer'.and.time_loop)then
           call BC_solar_wind_buffer
        else
           call BC_fixed(1,nVar,CellState_VI(:,iSide))
           if(UseB0)call BC_fixed_B
        end if
     case('fixedB1','fixedb1')
        call BC_fixed(1,nVar,CellState_VI(:,iSide))
     case('shear')
        call BC_shear(1,nVar,iSide)
     case('none')
     case default
        IsFound=.false.
        if(UseUserOuterBcs .or. TypeBc_I(iSide) == 'user')&
           call user_set_outerBCs(iBLK,iSide,TypeBc_I(iSide),IsFound)

        if(.not. IsFound) call stop_mpi( &
             NameSub // ': unknown TypeBc_I=' //TypeBc_I(iSide))
     end select

     if(DoSetEnergy)then
        call calc_energy(imin1g,imax1g,jmin1g,jmax1g,kmin1g,kmax1g,iBLK,&
             1,nFluid)
        call calc_energy(imin2g,imax2g,jmin2g,jmax2g,kmin2g,kmax2g,iBLK,&
             1,nFluid)
     end if
  end do

  if(DoTestMe)then
     do iVar=1,nVar
        write(*,*)'final',NameVar_V(iVar),'   cell,ghost,ghost2=',&
             State_VGB(iVar,Itest,Jtest,Ktest,iBLK),&
             State_VGB(iVar,Ighost,Jghost,Kghost,iBLK),&
             State_VGB(iVar,Ighost2,Jghost2,Kghost2,iBLK)
     end do
  end if

end subroutine set_outer_BCs

!==========================================================================  
subroutine BC_cont(iVarStart,iVarLast)
  use ModVarIndexes
  use ModAdvance, ONLY : State_VGB
  use ModSetOuterBC

  ! Continuous: q_BLK(ghost)= q_BLK(phys1)

  integer, intent(in) :: iVarStart,iVarLast

  State_VGB(iVarStart:iVarLast,imin1g:imax1g,jmin1g:jmax1g,kmin1g:kmax1g,iBLK)=&
       State_VGB(iVarStart:iVarLast,imin1p:imax1p,jmin1p:jmax1p,kmin1p:kmax1p,iBLK)
  State_VGB(iVarStart:iVarLast,imin2g:imax2g,jmin2g:jmax2g,kmin2g:kmax2g,iBLK)=&
       State_VGB(iVarStart:iVarLast,imin1p:imax1p,jmin1p:jmax1p,kmin1p:kmax1p,iBLK)

end subroutine BC_cont

!==========================================================================
subroutine BC_shear(iVarStart, iVarLast ,iSide)
  use ModVarIndexes
  use ModAdvance, ONLY : State_VGB
  use ModSetOuterBC
  use ModSize
  use ModPhysics
  implicit none
  ! Shear: q_BLK(ghost)= q_BLK(phys1+shear)

  integer, intent(in) :: iVarStart, iVarLast, iSide
  integer :: iVar, Dn
  !------------------------------------------------------------------------
  ! For the corners or bot_ and top_ fill with unsheared data first
  call BC_cont(iVarStart,iVarLast)

  ! If the shock is not tilted, there is nothing to do
  if(abs(ShockSlope)<cTiny) RETURN

  do iVar = iVarStart, iVarLast
     ! Shear according to ShockSlope
     if(ShockSlope < -cTiny)then
        call stop_mpi('ShockSlope must be positive!')
     elseif(ShockSlope >= cOne)then
        Dn = nint(ShockSlope)
        if(abs(Dn-ShockSlope)>cTiny)&
             call stop_mpi('ShockSlope > 1 should be a round number!')
        select case(iSide)
           ! Shift parallel to Y by 1 but copy from distance Dn in X
        case(east_)
           State_VGB(iVar,     imin1g,    jmin1g+1:jmax1g,   kmin1g:kmax1g,iBLK)=&
                State_VGB(iVar,imin1g+Dn, jmin1p  :jmax1p-1, kmin1p:kmax1p,iBLK)

           State_VGB(iVar,     imin2g,    jmin2g+1:jmax2g,   kmin2g:kmax2g,iBLK)=&
                State_VGB(iVar,imin2g+Dn, jmin1p  :jmax1p-1, kmin1p:kmax1p,iBLK)
        case(west_)
           State_VGB(iVar,     imin1g,    jmin1g  :jmax1g-1, kmin1g:kmax1g,iBLK)=&
                State_VGB(iVar,imin1g-Dn, jmin1p+1:jmax1p,   kmin1p:kmax1p,iBLK)

           State_VGB(iVar,     imin2g,    jmin2g  :jmax2g-1, kmin2g:kmax2g,iBLK)=&
                State_VGB(iVar,imin2g-Dn, jmin1p+1:jmax1p,   kmin1p:kmax1p,iBLK)

           ! Shift parallel to X by Dn and 2*Dn
        case(south_)
           State_VGB(iVar,     imin1g+Dn:imax1g,   jmin1g, kmin1g:kmax1g,iBLK)=&
                State_VGB(iVar,imin1p:imax1p-Dn,   jmin1p, kmin1p:kmax1p,iBLK)
           State_VGB(iVar,     imin2g+2*Dn:imax2g, jmin2g, kmin2g:kmax2g,iBLK)=&
                State_VGB(iVar,imin1p:imax1p-2*Dn, jmin1p, kmin1p:kmax1p,iBLK)
        case(north_)
           State_VGB(iVar,     imin1g:imax1g-Dn,   jmin1g, kmin1g:kmax1g,iBLK)=&
                State_VGB(iVar,imin1p+Dn:imax1g,   jmin1p, kmin1p:kmax1p,iBLK)
           State_VGB(iVar,     imin2g:imax2g-2*Dn, jmin2g, kmin2g:kmax2g,iBLK)=&
                State_VGB(iVar,imin1p+2*Dn:imax1p, jmin1p, kmin1p:kmax1p,iBLK)
        end select
     else
        ! ShockSlope < 1
        Dn = nint(cOne/ShockSlope)
        if(abs(Dn-cOne/ShockSlope)>cTiny)call stop_mpi( &
             'ShockSlope < 1 should be the inverse of a round number!')
        select case(iSide)
           ! Shift parallel to Y by Dn
        case(east_)
           State_VGB(iVar,     imin1g, jmin1g+Dn:jmax1g,   kmin1g:kmax1g,iBLK)=&
                State_VGB(iVar,imin1p, jmin1p:jmax1p-Dn,   kmin1p:kmax1p,iBLK)

           State_VGB(iVar,     imin2g, jmin2g+2*Dn:jmax2g, kmin2g:kmax2g,iBLK)=&
                State_VGB(iVar,imin1p, jmin1p:jmax1p-2*Dn, kmin1p:kmax1p,iBLK)
        case(west_)
           State_VGB(iVar,    imin1g,  jmin1g:jmax1g-Dn,   kmin1g:kmax1g,iBLK)=&
                State_VGB(iVar,imin1p, jmin1p+Dn:jmax1p,   kmin1p:kmax1p,iBLK)
           State_VGB(iVar,     imin2g, jmin2g:jmax2g-2*Dn, kmin2g:kmax2g,iBLK)=&
                State_VGB(iVar,imin1p, jmin1p+2*Dn:jmax1p, kmin1p:kmax1p,iBLK)

           ! Shift parallel to X by 1, but copy from distance Dn in Y
        case(south_)
           State_VGB(iVar,     imin1g+1:imax1g,   jmin1g,    kmin1g:kmax1g,iBLK)=&
                State_VGB(iVar,imin1p  :imax1p-1, jmin1g+Dn, kmin1p:kmax1p,iBLK)
           State_VGB(iVar,     imin2g+1:imax2g,   jmin2g,    kmin2g:kmax2g,iBLK)=&
                State_VGB(iVar,imin1p  :imax1p-1, jmin2g+Dn, kmin1p:kmax1p,iBLK)
        case(north_)
           State_VGB(iVar,     imin1g  :imax1g-1, jmin1g,    kmin1g:kmax1g,iBLK)=&
                State_VGB(iVar,imin1p+1:imax1p,   jmin1g-Dn, kmin1p:kmax1p,iBLK)
           State_VGB(iVar,     imin2g  :imax2g-1, jmin2g,    kmin2g:kmax2g,iBLK)=&
                State_VGB(iVar,imin1p+1:imax1p,   jmin2g-Dn, kmin1p:kmax1p,iBLK)
        end select
     end if
  end do
end subroutine BC_shear

subroutine BC_symm(iVarStart,iVarLast)
  use ModSetOuterBC
  use ModVarIndexes
  use ModAdvance, ONLY : State_VGB
  implicit none
  ! Mirror symmetry: q_BLK(ghost)= -q_BLK(phys)

  integer, intent(in) :: iVarStart,iVarLast

  State_VGB(iVarStart:iVarLast, imin1g:imax1g,jmin1g:jmax1g,kmin1g:kmax1g,iBLK)=&
       State_VGB(iVarStart:iVarLast,imin1p:imax1p,jmin1p:jmax1p,kmin1p:kmax1p,iBLK)
  State_VGB(iVarStart:iVarLast,imin2g:imax2g,jmin2g:jmax2g,kmin2g:kmax2g,iBLK)=&
       State_VGB(iVarStart:iVarLast,imin2p:imax2p,jmin2p:jmax2p,kmin2p:kmax2p,iBLK)

end subroutine BC_symm

!==========================================================================  
subroutine BC_asymm(iVarStart,iVarLast)
  use ModSetOuterBC
  use ModVarIndexes
  use ModAdvance, ONLY : State_VGB
  ! Mirror symmetry with sign change: q_BLK(ghost)= -q_BLK(phys)

  integer, intent(in) :: iVarStart,iVarLast

  State_VGB(iVarStart:iVarLast, imin1g:imax1g,jmin1g:jmax1g,kmin1g:kmax1g,iBLK)=&
       - State_VGB(iVarStart:iVarLast, imin1p:imax1p,jmin1p:jmax1p,kmin1p:kmax1p,iBLK)
  State_VGB(iVarStart:iVarLast,imin2g:imax2g,jmin2g:jmax2g,kmin2g:kmax2g,iBLK)=&
       -  State_VGB(iVarStart:iVarLast,imin2p:imax2p,jmin2p:jmax2p,kmin2p:kmax2p,iBLK)

end subroutine BC_asymm

!==========================================================================  
subroutine BC_fixed(iVarStart,iVarLast,q)
  use ModVarIndexes
  use ModAdvance, ONLY : State_VGB
  use ModSetOuterBC
  ! Set q_B=q in ghost cells

  integer, intent(in) :: iVarStart,iVarLast
  real,dimension(nVar), intent(in)    :: q
  integer::i,j,k
  do k=kmin1g,kmax1g; do j=jmin1g,jmax1g; do i=imin1g,imax1g
     State_VGB(iVarStart:iVarLast,i,j,k,iBLK)=&
          q(iVarStart:iVarLast)
  end do;end do;end do
  do k=kmin2g,kmax2g; do j=jmin2g,jmax2g; do i=imin2g,imax2g
     State_VGB(iVarStart:iVarLast,i,j,k,iBLK)=&
          q(iVarStart:iVarLast)
  end do;end do;end do

end subroutine BC_fixed

!==========================================================================  
subroutine BC_fixed_B
  use ModSetOuterBC
  use ModVarIndexes
  use ModAdvance, ONLY : State_VGB,B0_DGB
  ! Set q_B=q-q_B0 in ghost cells

  State_VGB(Bx_:Bz_,imin1g:imax1g,jmin1g:jmax1g,kmin1g:kmax1g,iBLK)= &
       State_VGB(Bx_:Bz_,imin1g:imax1g,jmin1g:jmax1g,kmin1g:kmax1g,iBLK)&
       - B0_DGB(:,imin1g:imax1g,jmin1g:jmax1g,kmin1g:kmax1g,iBLK)
 
  State_VGB(Bx_:Bz_,imin2g:imax2g,jmin2g:jmax2g,kmin2g:kmax2g,iBLK)= &
       State_VGB(Bx_:Bz_,imin2g:imax2g,jmin2g:jmax2g,kmin2g:kmax2g,iBLK)&
       - B0_DGB(:,imin2g:imax2g,jmin2g:jmax2g,kmin2g:kmax2g,iBLK)

end subroutine BC_fixed_B

!==========================================================================  

subroutine BC_solar_wind(time_now)

  use ModGeometry,ONLY:x_BLK, z_BLK,y_BLK, x2
  use ModVarIndexes
  use ModAdvance, ONLY : State_VGB, B0_DGB
  use ModSetOuterBC
  use ModMultiFluid, ONLY: &
       iRho_I, iUx_I, iUy_I, iUz_I, iRhoUx_I, iRhoUy_I, iRhoUz_I
  use ModPhysics, ONLY: LowDensityRatio
  use ModNumConst, ONLY: cTiny
  use ModSolarwind, ONLY: get_solar_wind_point
  implicit none

  ! Current simulation time in seconds
  real, intent(in) :: time_now 

  ! index and location of a single point
  integer :: i,j,k
  real :: x, y, z
  ! Varying solar wind parameters
  real :: SolarWind_V(nVar)
  !-----------------------------------------------------------------------
  
  do k = kmin1g, kmax1g 
     z = z_BLK(1,1,k,iBLK)
     do j = jmin1g, jmax2g
        y = y_BLK(1,j,1,iBLK)
        do i = imin1g, imax2g, sign(1,imax2g-imin1g)

           ! x= x_BLK(i,j,k,iBLK) ! for cell based BC this would be best
           x=x2 ! for face based and for west side as the inflow
           call get_solar_wind_point(time_now, (/x, y, z/), SolarWind_V)

           State_VGB(:,i,j,k,iBLK) = SolarWind_V

           ! Convert velocities to momenta
           State_VGB(iRhoUx_I, i,j,k,iBLK) = &
                State_VGB(iUx_I, i,j,k,iBLK)*State_VGB(iRho_I,i,j,k,iBLK)
           State_VGB(iRhoUy_I, i,j,k,iBLK) = &
                State_VGB(iUy_I, i,j,k,iBLK)*State_VGB(iRho_I,i,j,k,iBLK)
           State_VGB(iRhoUz_I, i,j,k,iBLK) = &
                State_VGB(iUz_I, i,j,k,iBLK)*State_VGB(iRho_I,i,j,k,iBLK)

           ! Subtract B0:   B1 = B - B0
           State_VGB(Bx_:Bz_,i,j,k,iBLK)    = &
                State_VGB(Bx_:Bz_,i,j,k,iBLK) - B0_DGB(:,i,j,k,iBLK)
        end do
     end do
  end do
  
end subroutine BC_solar_wind

!==========================================================================  

subroutine BC_solar_wind_buffer

  use ModGeometry, ONLY: z_BLK, y_BLK
  use ModVarIndexes, ONLY: Bx_, By_, Bz_
  use ModAdvance, ONLY : State_VGB, B0_DGB
  use ModSetOuterBC

  implicit none

  ! index and location of a single point
  integer :: i, j, k
  real    :: y, z
  !-----------------------------------------------------------------------
  do k=kmin1g,kmax1g 
     z = z_BLK(1,1,k,iBLK)
     do j=jmin1g,jmax2g
        y = y_BLK(1,j,1,iBLK)
        do i=imin1g,imax2g,sign(1,imax2g-imin1g)
           call read_ih_buffer(y,z,State_VGB(:,i,j,k,iBlk))
           ! Subtract B0
           State_VGB(Bx_:Bz_,i,j,k,iBLK) = State_VGB(Bx_:Bz_,i,j,k,iBLK) &
                - B0_DGB(:,i,j,k,iBLK)
        end do
     end do
  end do

end subroutine BC_solar_wind_buffer
!==============================================================================

subroutine set_radiation_outflow_bc(iVarFirst, iVarLast, iSide)

  use ModSetOuterBC
  use ModVarIndexes
  use ModAdvance,  ONLY: State_VGB, nI, nJ, nK, nOpacity
  use ModGeometry, ONLY: dx_BLK, dy_BLK, dz_BLK
  use ModPhysics,  ONLY: Si2No_V, UnitX_
  use ModUser,     ONLY: user_material_properties
  implicit none

  integer, intent(in) :: iVarFirst, iVarLast, iSide

  integer :: iVar, i, j, k, iOpacity
  real :: OpacitySi_I(nOpacity), Coef
  !----------------------------------------------------------------------------

  select case(iSide)
  case(1,2)
     do k = 1, nK; do j = 1, nJ
        call user_material_properties(State_VGB(:,imin1p,j,k,iBLK), &
             imin1p, j, k, iBLK, DiffusionOpacitySiOut_I=OpacitySi_I)

        do iVar = iVarFirst, iVarLast
           iOpacity = iVar - iVarFirst + 1
           Coef = 2/sqrt( &
                (3*OpacitySi_I(iOpacity)/Si2No_V(UnitX_) * dx_BLK(iBLK))**2 &
                + ((State_VGB(iVar,imin2p,j,k,iBLK) &
                -   State_VGB(iVar,imin1p,j,k,iBLK)) &
                /  State_VGB(iVar,imin1p,j,k,iBLK))**2)
           State_VGB(iVar,imin1g,j,k,iBLK)=State_VGB(iVar,imin1p,j,k,iBLK) &
                *(Coef - 0.5)/(Coef + 0.5)
           State_VGB(iVar,imin2g,j,k,iBLK) &
                = 2*State_VGB(iVar,imin1g,j,k,iBLK) &
                -   State_VGB(iVar,imin1p,j,k,iBLK)
        end do
     end do; end do

  case(3,4)
     do k = 1, nK; do i = 1, nI
        call user_material_properties(State_VGB(:,i,jmin1p,k,iBLK), &
             i, jmin1p, k, iBLK, DiffusionOpacitySiOut_I=OpacitySi_I)

        do iVar = iVarFirst, iVarLast
           iOpacity = iVar - iVarFirst + 1
           Coef = 2/sqrt( &
                (3*OpacitySi_I(iOpacity)/Si2No_V(UnitX_) * dy_BLK(iBLK))**2 &
                + ((State_VGB(iVar,i,jmin2p,k,iBLK) &
                -   State_VGB(iVar,i,jmin1p,k,iBLK)) &
                /  State_VGB(iVar,i,jmin1p,k,iBLK))**2)
           State_VGB(iVar,i,jmin1g,k,iBLK)=State_VGB(iVar,i,jmin1p,k,iBLK) &
                *(Coef - 0.5)/(Coef + 0.5)
           State_VGB(iVar,i,jmin2g,k,iBLK) &
                = 2*State_VGB(iVar,i,jmin1g,k,iBLK) &
                -   State_VGB(iVar,i,jmin1p,k,iBLK)
        end do
     end do; end do

  case(5,6)
     do j = 1, nJ; do i = 1, nI
        call user_material_properties(State_VGB(:,i,j,kmin1p,iBLK), &
             i, j, kmin1p, iBLK, DiffusionOpacitySiOut_I=OpacitySi_I)

        do iVar = iVarFirst, iVarLast
           iOpacity = iVar - iVarFirst + 1
           Coef = 2/sqrt( &
                (3*OpacitySi_I(iOpacity)/Si2No_V(UnitX_) * dz_BLK(iBLK))**2 &
                + ((State_VGB(iVar,i,j,kmin2p,iBLK) &
                -   State_VGB(iVar,i,j,kmin1p,iBLK)) &
                /  State_VGB(iVar,i,j,kmin1p,iBLK))**2)
           State_VGB(iVar,i,j,kmin1g,iBLK)=State_VGB(iVar,i,j,kmin1p,iBLK) &
                *(Coef - 0.5)/(Coef + 0.5)
           State_VGB(iVar,i,j,kmin2g,iBLK) &
                = 2*State_VGB(iVar,i,j,kmin1g,iBLK) &
                -   State_VGB(iVar,i,j,kmin1p,iBLK)
        end do
     end do; end do
  end select

end subroutine set_radiation_outflow_bc
!==============================================================================

