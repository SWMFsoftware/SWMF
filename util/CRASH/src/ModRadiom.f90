!------- radiom.f90
! \
module M_RADIOM
  ! /
  ! The module solves the main equations of the non-LTE
  ! model
  !
  use CRASH_M_projE
  !Projects the radiation energy densities in user-defined groups 
  !onto a refined grid  
  implicit none
  PRIVATE !Except

  !\
  ! Parameters for writing the version
  !/
  character(LEN=*),parameter :: author='M.Busquet' &
       ,wher='ARTEP,inc' &
       ,version='CRASH-2011.005.g' &
       ,library='RADIOM' &
       ,dateModified='2011-08-09' 
  character(LEN=*),parameter :: &
       compiled='2011.09.02 13:38:10' &
       ,collab='M.Klapisch & I.Sokolov' 
  character(LEN=20),parameter :: &
       credits='D.Fyfe & J.H.Gardner'
  !  Flags :
  logical,save :: notPrepUbar = .false.,dbgRadiom=.false.
  !  Options :
  logical,save :: jhg=.false.,sumMode=.false.,projMode=.true. &
       ,doNew=.true.
  real,save :: Auger=1,Auger_1oB=20d0,Auger_B=0,Gaunt=0.2 &
       ,Auger_A=0,Auger_1oA=0
  ! work arrays for "corrUbar" :
  integer,parameter :: mxN=1+20*5 , mxu=1+10*5 ,mxGr=300
  real :: QJ_1,QJ_2,dQJ ,corrUs(mxN),QJs(mxN)
  logical,save :: setXubar=.true.

  real,parameter :: aSaha=6.02e21, b_CovR=1.34e13*0.2 ,one=1.0
  public:: caltz0,prep_projE,PrepCorrUBar

  !-------
contains
  !-------
  subroutine printVersion()
    write(*,*)'...................................'
    write(*,*)'library=',library
    write(*,*)'version=',version
    write(*,*)'modified on ',dateModified
    write(*,*)'compiled on ',compiled
    write(*,*)'author=',author
    write(*,*)' at ',wher
    write(*,*)' in collaboration with : ',collab
    write(*,*)' credits also to : ',credits
    write(*,*)'...................................'
  end subroutine printVersion
  !-------
  subroutine setDbgRadiom(flag)
    logical,intent(IN) :: flag
    dbgRadiom=flag
    write(*,*)'- - flag  dbgRadiom  set to ',dbgRadiom
  end subroutine setDbgRadiom
  !-------
  subroutine xubar(Te_in,Ne_in ,ngr,grBound,EoB,ubar)
    !  find the center of C.S.distribution in normalized energy,  U=hnu/Te
    !  it uses a "non-overflow, non-sensitive for tails, algorithm" derived by G.Schurtz
    !  Un=
    !  Rn=
    !  ubar=lim(n=Max;UN/Rn)
    !
    use CRASH_M_projE	! ,only : mxOut
    use CRASH_M_expTab
    real,intent(Out)   :: ubar
    integer,intent(IN) :: ngr
    real,intent(IN)    :: Te_in,Ne_in
    real,intent(IN)    :: grBound(0:ngr),EoB(ngr)

    real,parameter :: A_LTE_limit=7.38905609	! =exp(2)
    real,parameter :: smallR=1d-200,two=2

    integer :: nfgp,nOut,ig
    real :: ts,at32s,QJ,betapm ,bu3
    real :: du,duu,s,rm,frad
    logical :: gotoLTE
    real :: u,ey		! for 'exp_tab'	 : u -> ey=exp(-u)

    logical,save :: dbg=.false.,lastWasPrep=.false.

    nfgp=ngr+1
    ts=sqrt(Te_in)
    at32s= aSaha*Te_in*ts/Ne_in
    if(at32s.le.A_LTE_limit) then
       ubar=log(at32s)	! thus Tz/Te will be  1
       return
    end if
    QJ=log(at32s)
    betapm= (b_CovR/Gaunt/aSaha)*at32s*Te_in**2        


    !	     --------
    if(notPrepUbar)then
       if(Efirst.ne.grBound(0) .or. Elast.ne.grBound(ngr) &
            .or. nbIn.ne.ngr) then
          write(*,783) nbIn,ngr,Efirst,Elast,grBound(0),grBound(ngr)
783       format(//,'-W- definition of groups changed : ngr :',i5,' -> ',i5 &
               ,' range:',1p,2e13.4,' -> ',2e13.4,//)
          call CON_stop('-E- xubar: groups changed')
          call prep_projE(grBound,ngr)	
       end if
       lastWasPrep=.false.
    else 
       lastWasPrep=.true.
    end if
    call projSP_E(Te_in,EoB,nOut,gotoLTE)		! ,Uout,SPout
    !	     --------

    if(gotoLTE) then
       ubar=QJ ! =log(at32s)		! thus Tz/Te will be  1
       return
    end if

    ! perform the "non overflow" alogrithm, u=including  jhg correction
    ubar=0
    s=one
    if(nbOut.lt.2) call CON_stop('-R- xubar : nbOut<2')
    duu=Uout(2)-Uout(1)
    do ig=2,nbOut+1
       if(Uout(ig).le.0) then
          exit		! 111003
       end if
       u  = (Uout(ig)+Uout(ig-1))/two
       du = (Uout(ig)-Uout(ig-1))
       bu3=betapm*u**3
       include 'exp_tab.h'
       frad=at32s*ey 					&
            *(Auger+Auger_1oB/betapm 			&
            +Auger_B*betapm 				&
            +bu3*SPout(ig-1)) 				&
            /(Auger +Auger_1oB/betapm 			&
            +Auger_B*betapm				&
            +bu3					)
       rm=frad/s*(du/duu)
       if(rm.le.smallR) exit		! 080712 ! ? save CPU
       rm=max(rm,smallR)
       duu=du
       ubar=(ubar+u*rm)/(one+rm)
       if(dbg)then
          write(*,'(1p,a,i4,a,e13.4,a,5e13.4,a,2e13.4)') &
               'group#',ig,' ubar=',ubar,' u,sp,frad,rm,s=',u,SPout(ig-1),frad,rm,s 
       end if
       s=one+one/rm		! sumMode
    end do
    ! 
    call CorrUbar(ubar,QJ)

  end subroutine xubar

  !===================

  subroutine xubar0(Te_in,Ne_in ,EoB,ubar)
    ! 
    ! same as "xubar" w/o check of Eground
    !
    use CRASH_M_projE	! ,only : mxOut
    use CRASH_M_expTab
    real :: ubar
    integer :: ngr
    real,intent(IN) :: Te_in,Ne_in
    real,intent(IN) :: EoB(nbIn)

    real,parameter :: A_LTE_limit=7.38905609	! =exp(2)
    real,parameter :: smallR=1d-200,two=2

    integer :: nfgp,nOut,ig
    real :: ts,at32s,QJ,betapm ,bu3
    real :: du,duu,s,rm,frad
    logical :: gotoLTE
    real :: u,ey		! for 'exp_tab'  u -> ey=exp(-u), ex_u=1-ey

    logical,save :: dbg=.false.,lastWasPrep=.false.

    ngr=nbIn
    nfgp=ngr+1
    ts=sqrt(Te_in)
    at32s= aSaha*Te_in*ts/Ne_in
    if(at32s.le.A_LTE_limit) then
       ubar=log(at32s)	! thus Tz/Te will be  1
       return
    end if
    QJ=log(at32s)
    betapm= (b_CovR/Gaunt/aSaha)*at32s*Te_in**2        
    !	     --------
    call projSP_E(Te_in,EoB,nOut,gotoLTE)	
    !	     --------
    if(gotoLTE) then
       ubar=QJ ! =log(at32s)		! thus Tz/Te will be  1
       return
    end if
    ! perform the "non overflow" alogrithm, u=including  jhg correction
    ubar=0
    s=one
    if(nbOut.lt.2) call CON_stop('-E- xubar0 : nbOut<2')
    duu=Uout(2)-Uout(1)
    do ig=2,nbOut+1
       u  = (Uout(ig)+Uout(ig-1))/two
       du = (Uout(ig)-Uout(ig-1))
       bu3=betapm*u**3
       include 'exp_tab.h'
       frad=at32s*ey 					&
            *(Auger+Auger_1oB/betapm 			&
            +Auger_B*betapm 				&
            +bu3*SPout(ig-1)) 				&
            /(Auger +Auger_1oB/betapm 			&
            +Auger_B*betapm				&
            +bu3					)
       rm=frad/s*(du/duu)
       if(rm.le.smallR) exit		! 080712 ! ? save CPU
       rm=max(rm,smallR)
       duu=du
       ubar=(ubar+u*rm)/(one+rm)
       s=one+one/rm		! sumMode
    end do
    call CorrUbar(ubar,QJ)
  end subroutine xubar0
  !-------
  !subroutine calTz(Te,Ne, Tz, hnug,eg,bg,ng)
    !-
    !-     THIS ROUTINE CALCULATES THE VALUE OF Tz IN BUSQUET THEORY FROM Te, IT
    !-     TAKES CALCULATED VALUES OF F^RAD(U_M) FOR EACH VALUE OF GROUP U_M;
    !-     USES THESE TO COMPUTE IONIZATION POPULATIONS AND THE FIRST MOMENT OF
    !-     U_M, UBAR. THIS UBAR IS THEN USED TO COMPUTE THE RATIO OF ELECTON TO 
    !-     IONIZATION TEMPERATURE SUCH THAT THE SAME DISTRIBUTION FUNCTION IS 
    !-     FOUND FOR THE LTE AT TZ AS THE NONLTE AT TE.
    !-
    !-     VARIABLES:
    !-         TEOTZ (DIMLESS) = "T_e/T_z" = RATIO OF ELECTRON AND "REDUCED" TEMPS
    !-         UM (DIMLESS) = "U_M" = REDUCED MEAN FREQUENCIES (H NU / K T_E)
    !-         SPINT (eV**3/h**2/c**3) = "SPECTRAL INTENSITY" OF RADIATION IN
    !-                                    REDUCED MEAN GROUP U_M; HYDRO CODE
    !-                                    USES U_\nu = {\int I_\nu d\Omega}
    !-         TE (eV) = ELECTRON TEMPERATURE
    !-         TZ        (RESULT)
    !-         DEN(1/CCM) = ELECTRON NO. DENSITY 
    !-         UBAR (DIMLESS) = "U-BAR" = ION POT / K T_E = REDUCED ION POT
    !-         FRADLN (DIMLESS) = "LN(F^RAD)" = NATURAL LOG OF BUSQUET'S F^RAD
    !-                                          FUNCTION (EQ. 53)'
    !
    !  Te is the electronic temperature in eV
    !  Ne is the estimated electronic density ( = Zbar * ro/A * avoadro),  in  cm-3
    !   Tz  will be the "ionization temperature",  in eV
    !   Zbar will then be obtained from the LTE E.O.S. package, but at Tz
    !
    !     it has to be iterated so Zbar(Tz) = estimated Zbar
    !
    !  ng  is the number of frequency groups
    !
    !   eg(1:ng)  is to be proportional to the "radiative energy density", integrated on the frequ. bin
    !   be(1:ng)  is to be proportional to the Planckian function (computed at Te)
    !   so   eg(:)/bg(:)  =  1  if the radiation field is Planckian
    !   eg(:)=0  define an "optically thin approximation"
    !
    !   hnug(0:ng)  are  the group boundaries, in same units than Te (usually eV)
    !
    !
   ! use CRASH_M_projE,only : nbIn

    !Inputs
    !real,intent(IN) :: Te,ne
    !integer,intent(IN) :: ng
    !real,dimension(ng),intent(IN) :: eg,bg
    !real,dimension(0:ng),intent(IN) :: hnug
    !real,intent(OUT) :: Tz

    !real,dimension(mxgr) :: EoB
    !real,parameter :: smallB=1d-30
    !integer :: ig

    !real :: at32s,ubar,QJ

    !if(nbIn.le.0) then
    !   write(*,*)'-P- prep_projE not done before "calTE"'
    !   call CON_stop('-E- calTZ: prep_projE not done')
    !end if
    !at32s=aSaha*te*sqrt(te)/Ne
    !tz=te
    !if(at32s <= one) return
    !if(ng.gt.0) then
    !   EOB(1:ng)=0.
    !   do ig=1,ng
    !      if(bg(ig) >smallB) EoB(ig)=eg(ig)/bg(ig)
    !   end do
    !   !	elseif(ng.lt.0) then
    !   !  may use empirical fit ...
    !   !	 return
    !else
    !   write(*,*)'-E- CALTZ: ng=0, should use an artifical gridding ...'
    !   call CON_stop('-P- no Eg(:) in CALTZ')
    !end if
    !call xubar(te,ne,ng,hnug,EoB ,ubar)	! no more a function
    !QJ=log(at32s)
    !tz=te*ubar/QJ
  !end subroutine calTz
  !-------
  subroutine calTz0(Te,Ne, Tz, EoB)
    !-
    !-     same as "calTZ" but w/o hnug(0:ng) & ubar1
    !
    !
    use CRASH_M_projE,only : nbIn
    real,intent(IN) :: Te,ne
    real,dimension(:),intent(IN) :: EoB
    real,intent(OUT) :: Tz

    real,parameter :: smallB=1d-30
    integer :: ig

    real :: at32s,ubar,QJ
    if(nbIn.le.0) then
       write(*,*)'-P- prep_projE not done before "caltz0"'
       call CON_stop('-P- calTZ: prep_projE not done')
    end if
    at32s=aSaha*te*sqrt(te)/Ne
    tz=te
    if(at32s.le.one) return
    ! this is already done in "setErad", should pass by arg ?
    call xubar0(te,ne,EoB ,ubar)	! no more a function
    QJ=log(at32s)
    tz=te*ubar/QJ
  end subroutine calTz0
  !-------
  subroutine calte(TeOld,Ne,Tz, Tenew, hnug,eg,bg,ng)
    use CRASH_M_projE,only : nbIn
    implicit none
    real,intent(IN) :: TeOld,Tz,Ne
    real,intent(OUT) :: TeNew
    integer :: ng
    real,dimension(0:ng),intent(IN) :: hnug
    real,dimension(ng),intent(IN) :: eg,bg

    real,dimension(mxgr) :: EoB
    real,parameter :: smallB=1d-30
    integer :: ig

    real :: at32s,ubar,QJ

    if(nbIn.le.0) then
       write(*,*)'-P- prep_projE not done before "calte"'
       call CON_stop('-P- calTE: prep_projE not done')
    end if
    TeNew=TeOld
    if(Ne.le.one) return
    at32s=aSaha*TeOld*sqrt(TeOld)/Ne
    if(at32s.le.one) return
    QJ=log(at32s)
    EOB(1:ng)=0.
    do ig=1,ng
       if(bg(ig).gt.smallB) EoB(ig)=eg(ig)/bg(ig)
    end do
    call xubar(TeOld,ne,ng,hnug,EoB,ubar)
    TeNew=Tz*QJ/ubar
  end subroutine calte
  !-------
  subroutine calte0(Te,Ne,Tz, EoB)
    use CRASH_M_projE,only : nbIn
    implicit none
    real,intent(IN) :: Tz,Ne
    real,intent(INOUT) :: Te
    real,dimension(:),intent(IN) :: EoB
    real,parameter :: zero=0

    real,parameter :: smallB=1d-30
    integer :: ig

    real :: at32s,ubar,QJ

    if(nbIn.le.0) then
       write(*,*)'-P- prep_projE not done before "calte0"'
       call CON_stop('-P- calTE0: prep_projE not done')
    end if

    if(Ne.le.one) return
    if(Te.le.0) then
       at32s=aSaha*Tz*sqrt(Tz)/Ne
    else
       at32s=aSaha*Te*sqrt(Te)/Ne
    end if
    Te=Tz
    if(at32s.le.one) return
    QJ=log(at32s)
    call xubar0(Te,ne,EoB,ubar)
    Te=Tz*QJ/ubar
    return
  end subroutine calte0
  !-------
  subroutine CorrUbar(ubar,QJ)

    real, intent(IN) :: QJ
    real :: ubar,d
    real, parameter :: QJ_limit=2.1
    integer :: iq
    !	
    if(setXubar) return
    if(.not.notPrepUbar) call prepCorrUbar()
    ! 
    iq=int((QJ-QJ_1)/dQJ) +1
    if(iq <= 0 .or. QJ <= QJ_limit) then
       ubar=QJ
    elseif(iq.ge.mxN) then
       ubar=ubar*corrUs(mxN)
    else
       d=(QJ-QJs(iq))/dQJ
       d=(one-d)*corrUs(iq)+d*corrUs(iq+1)
       ubar=ubar*d
    end if
    ! 
  end subroutine CorrUbar
  !-------
  subroutine prepCorrUbar()
    !
    !  tabulate correction to UBAR (see corrUbar)
    !
    use CRASH_M_projE,only : Umin,Umax	,nbIn,Efirst,Elast

    real :: ubar,ne,te,at32s,qj,r,eg1
    real,parameter :: zero=0,two=2 
    real,dimension(0:mxN) :: Ugs
    real,dimension(mxN+1) :: EovB
    integer :: iq
    ! 
    te=0.1
    r=(Umax/Umin)**(one/mxu)
    eg1=Umin*te
    do iq=1,mxN
       EovB(iq)=1	! 0 & 1 give same result
       Ugs(iq-1)=eg1
       eg1=eg1*r
    end do
    Ugs(mxN)=eg1
    notPrepUbar=.false.
    setXubar=.true.
    QJ_1=2.
    QJ_2=log(6.d4)
    dQJ=(QJ_2-QJ_1)/(mxN-1)
    qj=QJ_1
    do iq=1,mxN
       at32s=exp(qj)
       ne=te*sqrt(te)*(aSaha/at32s)
       call xubar(Te,Ne,mxu,Ugs,EovB,ubar)
       corrUs(iq)=qj/ubar
       QJs(iq)=qj
       qj=qj+dQJ
    end do
    ! 
    setXubar=.false.
    notPrepUbar=.true.	

  end subroutine prepCorrUbar
  !-------
  ! \
end module M_RADIOM
 ! /
