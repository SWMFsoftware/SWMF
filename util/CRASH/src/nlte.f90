!
! ***  nlte.f90  ***
!
!
!  code will use :
!
! initialisation:
!	call exp_tab8()		! prepare tabulated exponentials
!		call prepCorrUbar()	! or it will called automatically
!	call prep_projE(hnuGr(0:nbGr),nbGR)	! using group. definition of the run
!
!  runtime:
!       call set_kbr
!	call nLTE_EOS(...
!	call correctKnuJnu(...	! not yet implemented
!
! ToDo :
!	- check Auger
!	- use Sum(E/B...) for FB thresholds ?
!	- grey processing 
!                          AT32S * exp(-UM/gama)
!                  1 =  -----------------------  == AT32S * exp(-Um*Te/Tz/gama)
!                           1 + BT72S * UM**3
!	- analytical fit in transparent approximation
!
!-----------
! \ 
module CRASH_M_NLTE
  ! /
  !This module solves the equation for the effective energy.
  !E_{eff} = E_{tot}^{NLTE} - (3/2) k_B x rho x Z^*(T_Z) * (T_e-T_z)
  !where E_{tot} is the input parameter for the "inverse" equation of state
  !(from which the temperatures should be solved) and E_{eff} is the tabulated
  !LTE equation of state with (rho, T_Z) being the input parameters for the 
  !latter.
  use M_RADIOM,only : caltz0!,calte0,caltz,calte
  use CRASH_M_projE,only : mxOut
  use CRASH_M_localProperties
  implicit none

  logical :: useLTE=.false.
  logical :: useEElog=.false.
  logical :: useZbrent=.true.
  logical :: dbg=.false.

  !  NOTE :  use either    (useEElog=.F.,useZbrent=.F.)  or  (useEElog=.T.,useZbrent=.F.)
  !
  !  external call :   
  !		call LTE_EOS_inv
  !		call LTE_EOS_dir
  !
  !  user should fill these before calling  LTE_EOS, 
  ! Te 	: eV
  ! Ne	: el/cm3
  ! Natom : atom/cm3
  ! Etot	: internal energy
  ! Ptot	: total pressure
  !
  ! KbR_E, kBr_P :   C_subV, C_subP , including ro if required,
  !  i.e.    Etot=1.5 *(Zbar+1) * kbR_E * T
  !	   Ptot=     (Zbar+1) * kbR_P * T
  ! Efloor : lower bound of total Energy in the table.
  !	   when  Etot.LE.Efloor, LTE will be turned on
  !
  ! \
  !  CAUTION !
  ! /
  ! Ecold,Tcold : lower bound of total Energy in the table for the actual density
  !	during the process of NLTE< it may come that the input value EEeff turns less than Ecold
  !      Tcold will be returned
  !
  !  useLTE : .T. to save CPU (ex. for H at high temperature): forces LTE
  !
  integer,save :: ng_rad=-1
  ! real,save,dimension(mxout) :: Erad,Brad
  real,save,dimension(mxout) :: EoB
  real,save,dimension(0:mxout) :: hnu_rad
  ! parameter for LTE_EOS_dir  convergence on 'Zbar'
  real,save :: epsD=1d-5
  integer,save :: niterMax=20

  !================Interface to the equation of state========
  !n	interface
  !n	  subroutine  LTE_EOS_inv(te,Etot,Ptot,Zbar,Cv)
  !n	    implicit none
  !n	    real,intent(IN) :: Etot
  !n	    real,intent(OUT) :: te,Ptot,Zbar,Cv
  !n	  end subroutine LTE_EOS_inv
  !n	end interface
  !-------
contains
  !-------

  function TZdif(Tzold , Ne)
    implicit none
    real,intent(IN) :: Tzold ,Ne
    real :: Te,Tz,TZdif,Ez,TzEOS,Pe,Cv,Zbar,Etot
    !
    !- Ne    : electronic density [cm-3]
    !- Natom : atomic density [cm-3]
    !- Te    : electronic temperature in eV
    !
    TZeos=TZold
    call caltz0(Te,Ne, Tz, EoB)	! Te,Ne, Erad,Brad from M_localProperties
    Ez=max( Etot-kBr_E*(Te-Tz) , Efloor)
    !- waiting for using the real direct EOS
    call LTE_EOS_inv(tz,Ez,Pe,Zbar,Cv)	! ro : in module M_localProperties
    TZdif=Tz-TZeos
    stop 'TZdif not implemented'

    return
  end function TZdif
  !-------
  subroutine NLTE_EOS(Natom,RO_in, Te_in,Ee_in,Et_in,Pe_in,Pt_in &
       , estim_Zbar,estim_Tz,estim_Te  &
       ,Zbar_out,Tz_out  ,Te_out,Ee_out,Et_out,Pe_out,Pt_out,Cv_out)
    !
    !  to account for  Et=Ee+Ei w/ Te=Ti : use  Ee_* or Et_*
    !
    !  direct/inverse EOS is used depending of presence of  E*_in or Te_in
    ! 
    !   CALL SETERAD(Erad,Brad,hnuG,nbG)  has been set before.
    !
    use CRASH_M_localProperties,only : ro,Ni
    implicit none
    real,optional,intent(IN) :: Natom,ro_in 	! atomic density in   atom/cm3 
    real,optional,intent(IN) :: Te_in,Ee_in,Et_in,Pe_in,Pt_in 	! one and only one
    real,optional,intent(IN) :: estim_Zbar,estim_Tz,estim_Te
    real,optional,intent(OUT) :: Zbar_out,Tz_out
    real,optional,intent(OUT) :: Te_out,Ee_out,Et_out,Pe_out,Pt_out,Cv_out
    real :: Ee,Pe,Ne,d
    real :: Te,Tz,zbar,Cv,Eeff ,Elow,Ehigh,difLo,difHi
    integer :: niter
    real,parameter :: x_0=0,x_3o2=1.5d0,x_1=1d0,x_1o2=0.5d0,two=2d0,three=3d0

    ! variables for NLTE EOS_inv
    real,external :: zbrentEE
    !  smallF  : root declared found when |diff| < smallF
    !	real,parameter :: smallF=1e-4	! OK for "smooth EOS", like  ZTF_EOS
    real,parameter :: smallF=1e-5	! with "tlunv-Au", little defect at ro=300, Te.GE.6e3, w/ useEElog=.F.
    real,parameter :: tol_0=0
    real,parameter :: tol_EE=1d-6	! convergence param. for zbrentEE (linear function)
    real,external :: rootEE  ,EEdiff,EEdiff_ln
    logical :: inverse,direct,oneT,twoT
    real :: Ecold,Tcold

    if(present(Natom)) then
       !	 if(NIgiven .or. ROgiven) goto 104
       Ni=Natom
       if(present(RO_in)) goto 102
       ro=(Ni/avogadro)*atoMass
    else
       if(.not.present(RO_in)) then
       !   if(.not.ROgiven .and. .not.NIgiven) goto 102
       end if
       !	 if(NIgiven .or. ROgiven) goto 104
       ro=RO_in
       Ni=avogadro*(ro/atoMass)
    end if
    !
    ! check argument list consistency
    ! and check  1T (Te=Ti)  .or. 2T ?
    !
140 twoT=present(Ee_in).or. present(Ee_out) .or. present(Pe_in) .or. present(Pe_out)
    oneT=present(Et_in) .or. present(Pt_in) .or. present(Et_out) .or. present(Pt_out)
    if(twoT)then
       if(oneT) goto 103
       zion=0
       direct=present(Te_in) .or. present(Ee_out)
       inverse=present(Ee_in) .or. present(Pe_in) .or. present(Te_out) 
    else
       if(.not.oneT) goto 103
       zion=1 
       direct=present(Te_in) .or. present(Et_out)
       inverse=present(Et_in) .or. present(Pt_in) .or. present(Te_out)
    end if
    if(inverse.and.direct) goto 101
    if(.not.inverse.and..not.direct) goto 101

    if(ng_rad.lt.0) then
       write(*,*)'-P- setErad has to be called before,'
       write(*,*)'    whenever Erad or Brad are modified'
       stop '-P- setErad has to be called before NLTE_OPAC'
    end if

    !--
    !- direct EOS -
    if(direct) then
       !--
       if(dbg) write(*,*)'direct NLTE_EOS(te_in=',te_in,'  useLTE=',useLTE
       te=Te_in
       if(useLTE) then
          !- waiting for using the real direct EOS
          call LTE_EOS_dir(te,Ee,Pe,Zbar,Cv)	! ro : in module M_localProperties, can be put in arg.
          tz=te
       else
          if(present(estim_Zbar)) then
             zbar=estim_Zbar
          else
             zbar=atoNum*x_1o2
             if(zbar.eq.x_0) zbar=x_1
          end if

          d=100.
	  niter=0
          DIR:	  do while ( abs(d).gt.epsD .and. niter.lt.niterMax)
             Ne=Zbar*Ni
             call calTz0(Te,Ne, Tz, EoB)
             if(dbg) write(*,*)'niter#',niter,' calTz0(Te=',te,' Ne=',ne,' -> Tz=',tz
             call LTE_EOS_dir(tz,Ee,Pe,Zbar,Cv)	! ro : in module M_localProperties
             if(zBar<=0.0)EXIT DIR	
             d=(Ne-Ni*zbar)/(Ne+Ni*zbar)
             if(dbg) write(*,*)'oldNe,newNe=',ne,zbar*ni,'  d=',d
             niter=niter+1
	  end do DIR
          if(dbg) write(*,*)'before correct : ee,te,tz,zbar=',ee,te,tz,zbar
          call correctEOS(zbar , Te, Tz ,Ee,Pe,Cv)
       end if
       !--
    else	! if(present(TE_in))
       !--
       !- inverse EOS -
       if(present(Ee_in) ) then
          ee=Ee_in
       else
          ee=Et_in
       end if

       !Caclulate the cold equation of state.
       call getEcold(ro,Ecold,Tcold)		! 110827

       if(dbg) write(*,*)'NLTE_EOS,inv : ee,Ecold=',ee,Ecold
       if(useLTE .or. (ee.le.Ecold)) then		! 110827
          !- 
          call LTE_EOS_inv(te,ee,pe,Zbar,Cv)	! ro, {Erad,Brad} : in module M_localProperties
          Ne=Zbar*Natom
          tz=te
       else
          !- inverse EOS, non LTE -
          !  should use bracket(Tzdiff,tz1,tz2) + zbrent(tzDif,tz1,tz2)  <<<<<<<
          !   	ro= ; Erad= ; Brad= 
          !   T	e= ; Ne=  : update at each iteration ? use Cv_loc ?

          !- waiting for using the real inverse EOS
          call LTE_EOS_inv(tz,ee,pe,Zbar,Cv)	! zion flag must be dealt with !!
          if(dbg) write(*,*)'see what LTE_inv would be : ee=',ee,' -> tz,zbar=',tz,zbar
!!!
!!!		 TAKE CARE OF   Te / Te+Ti ... if not in "LTE_EOS_inv"
!!!

          Ne=Zbar*Ni
          if(Zbar.lt.0.1) then
             ! for low Z, LTE is a sensible approximation
             te=tz
             goto 200
          end if
          !
          te=tz
          call calTZ0(Te,Ne, Tz, EoB)
          if(dbg) write(*,*)'calTZ : te,ne=',te,ne,'  -> tz=',tz
          Eeff=ee-x_3o2*kbR_E*(zbar+zion)*(Te-Tz)
          if(dbg) write(*,*)'ee=',ee,' -> Eeff=',Eeff,' Efloor=',Efloor
          if(Eeff.lt.Efloor) then
             ! for low EE, LTE is a sensible approximation
             if(dbg) write(*,*)'Eeff < Efloor  -> te=tz=',tz
             te=tz
             goto 200
          end if
          !
          if(dbg) write(*,*)'calls bracket_EE(ee=',ee
          call bracket_EE(EEdiff,ee,Elow,Ehigh,Efloor,difLo,difHi)
          if(dbg) write(*,*)'bracket_ee (ee=',ee,') -> Elo,Ehi=', Elow,Ehigh,' d=',difLo,difHi
          if(Elow.eq.Ehigh) then
             call LTE_EOS_inv(tz,Elow,pe,Zbar,Cv)
          else
             if(useZbrent) then
                if(useEElog) then
                   Eeff=zbrentEE(EEdiff_ln,log(Elow),log(Ehigh),log(Efloor),tol_EE,smallF,difLo,difHi,ee ,te,tz)
                else
                   Eeff=zbrentEE(EEdiff,Elow,Ehigh,Efloor,tol_EE,tol_0,difLo,difHi,ee ,te,tz)
                end if
             else
                if(dbg) write(*,*)'259.calls rootEE ...  useEElog=',useEElog
                if(useEElog) then
                   Eeff=rootEE(EEdiff_ln,log(Elow),log(Ehigh),log(Efloor),tol_EE,smallF,difLo,difHi,ee ,te,tz)
                else
                   if(dbg) write(*,*)'calls rootEE w/ Elo,Ehi,Emi=',Elow,Ehigh,Efloor &
                        ,'  dLo,dHI=',difLo,difHi
                   Eeff=rootEE(EEdiff   ,Elow,Ehigh,Efloor               ,tol_EE,smallF,difLo,difHi,ee ,te,tz)
                   if(dbg) write(*,*)'rootEE( Eeff=',Eeff,') ->  ee,te,tz=',ee,te,tz
                end if
             end if
          end if
          goto 200
          !

       end if	! useLTE/useEEdiff
       !--
    end if 	! if(present(TE_in))
    !--
    !x		if(dbg) write(*,*)'zbar,tz,te,ee,pe,Cv=',zbar,tz,te,ee,pe,Cv
200 if(present(zbar_out)) zbar_out=zbar
    if(present(Tz_out)) Tz_out=Tz
    if(present(Te_out)) Te_out=te
    if(present(Ee_out)) Ee_out=ee
    if(present(Et_out)) Et_out=ee
    if(present(Pe_out)) Pe_out=pe
    if(present(Pt_out)) Pt_out=pe
    if(present(Cv_out)) Cv_out=Cv
    ! call set_RO_Ni()	! reset ROgiven and NIgiven to .false.
    ! to make sure that "set_RO_NI" is called before NLTE_EOS if
    ! option arguments  "ro" and "Natom" are not used 
    return
    !
103 write(*,*)'-P- LTE_EOS require  using either  Ee/Pe .xor. Et/Pt'
    stop '-P- LTE_EOS.103'
102 write(*,*)'-P- LTE_EOS require  Natom .xor. RO_in, or prior use to "set_RO_NI"'
    stop '-P- LTE_EOS.102'
101 write(*,*)'-P- LTE_EOS require  TE_in .xor. EE_in/ET_in (direct EOS .xor. inverse EOS)'
    stop '-P- LTE_EOS.101'
104 write(*,*)'-W- LTE_EOS : use either optional argument (Ni or Ro) or prior call to Set_NI_RO !! '
    !	stop '-P- LTE_EOS.104'
    goto 140
    !
  end subroutine NLTE_EOS
  !-------
  subroutine correctEOS(zbar , Te,Tz &
       ,Ee,Pe,Cv)
    ! Te,Tz, kBro_E,kBro_P,Etot,Ptot : thorugh module variables
    use CRASH_M_localProperties,only : Zsmall
    implicit none
    real,intent(IN) :: zbar,Te,Tz
    real,intent(INOUT) :: Ee
    real,optional,intent(INOUT) :: Pe
    real,optional,intent(INOUT) :: Cv
    real :: zdt,zp
    real,parameter :: x_1=1d0,x_3o2=1.5d0
    ! kbRO is the proportionnality factor in the code units (generally = Boltzmann cst * density)
    !   so  3/2* kbRo *Te  is the kinetic (translational) energy at Te

    if(Te.le.0) return
    !  small Zbar would yield diverging the "EEeff equ." (see correctEOS + EEdiff)
    !no	zp=(zbar+zion)
    zp=(max(zbar,Zsmall) + zion)
    zdt=zp * ( Te-Tz)
    Ee=Ee + x_3o2 * kBr_E * zdt
    if(present(Pe)) Pe=Pe +         kBr_P * zdt
    if(present(Cv)) then
       Cv=Cv+ x_3o2 * kBr_E * zp
    end if
    return
  end subroutine correctEOS

  !----------      ***  eediff.f90 ***

  subroutine bracket_EE(diff_func,E_in,Elow,Ehigh,Emini,difLo,difHi)
    !
    ! at this time, bracket_EE cannot be use w/ EEdiff_ln, as well as internal handling, that as results forwarded to rootEE
    !
    use CRASH_M_localProperties,only : &
         ro	&	  ! bulk density
         ,Ni		  ! atom density ( = avo*ro/atoMass)
    implicit none
    !
    ! find pair of Eeff values that brackets the solution: EEdiff will have opposite signs at each ends
    ! 
    real,intent(IN) :: E_in,Emini
    real,intent(OUT) :: Elow,Ehigh
    real,parameter :: x_1o2=0.5d0,facLo=0.5, facHi=1.2,facHi2=1.5
    real :: difLo,difHi,Ebig ,te,tz,e,d
    real,external :: diff_func
    integer :: niter
    integer,parameter :: niterMax=100

    Ehigh=facHi*E_in		!! Ehigh can be > E_in if  Erad > Brad
    Ebig=10*E_in
    if(dbg) write(*,*)'bracket_EE.432 : E_in,ro,ni=',E_in,ro,ni,' calls EEdiff(',Ehigh,')'
    difHi=diff_func(Ehigh,E_in ,te,tz)
    if(dbg) write(*,*)'bracket_EE for E_in=',E_in,' eehi=',Ehigh,' -> difHI=',difHi ,' ro,ni=',ro,ni
    if(difHi.eq.0)then
       Elow=E_in
       if(dbg) write(*,*)'401 calls EEdiff(elow=',elow,' calls EEdiff(',Elow,')'
       difLo=diff_func(Elow,E_in ,te,tz)
       if(dbg) write(*,*)'bracket_EE for E_in=',E_in,' eeLO=',Ehigh,' -> difLO=',difLO ,' ro,ni=',ro,ni
       if(difLo.eq.0) then
          Ehigh=Elow
          goto 300
       else		! gives a chance to track the real root
          Ehigh=facHi*Ehigh
          if(dbg) write(*,*)'408 calls EEdiff(Ehi=',Ehigh,' calls EEdiff'
          difHi=diff_func(Ehigh,E_in ,te,tz)
          if((difHi*difLo).gt.0) then
             write(*,*) '-R- found null extremum of EEeff'
             stop '-R- found null extremum of EEeff'
          end if
          Elow=E_in
          goto 300
       end if
    else
       e=E_in
       d=diff_func(e,E_in ,te,tz)
       if((d*difHi).gt.0) then
          Ehigh=e
          difHi=d
       end if
    end if
    Elow=E_in*facLo
    if(dbg) write(*,*)'481: Ehigh,difHi=',Ehigh,difHi
    if(dbg) write(*,*)'bracket_EE.481 : Elow=',Elow,' calls EEdiff.'
    difLo=diff_func(Elow,E_in ,te,tz)
    if(dbg) write(*,*)'bracket_EE for E_in=',E_in,' eeLO=',Elow,' -> difLO=',difLO &
         ,' ro,ni=',ro,ni,' te,tz=',te,tz
    if(difLo.eq.0) then
       Ehigh=Elow
       goto 300
    end if
    if(dbg) write(*,*)'534: difLo,difHi=',difLo,difHi,' e=',Elow,Ehigh
    if((difLo*difHi).le.0) goto 300
    !
    ! Newton-Raphson like pushing the bracket limits
    niter=0
301 niter=niter+1
    if(dbg) write(*,*)'301 : niter=',niter
    if(niter.eq.niterMax) goto 302
    e= (Elow*difHi-Ehigh*difLo)/(difHi-difLo)
    e=Max(e,Emini)
    if(dbg) write(*,*)'545.difLo,difHi=',difLo,difHi &
         ,' E:',Elow,Ehigh,' Emini=',Emini,' -> e=',e
    if(dbg) write(*,*)'bracket_EE.496 : eLo,Ehi=',eLow,Ehigh,' e=',e &
         ,' d=',difLo,difHi,' calls EEdiff'
    d=diff_func(e,E_in ,te,tz)
    if(dbg) write(*,*)'bracket_EE.301 : Elo,Ehi,e=',Elow,Ehigh,e,' dif=',difLo,difHi,d
    if(abs(difLo).le.abs(difHi))then
       Ehigh=e
       if((d*difLo).le.0) then
          ! We found a bracket
          difHi=d
          ! if lower bound is Emini, try to rise it :
          goto 300
       end if
       if(abs(d).gt.max(abs(difLo),abs(difHi))) goto 302
       difHi=d
    else	! if(abs(difLo).gt.abs(difHi))
       Elow=e
       if((d*difHi).le.0) then
          difLo=d
          goto 300
       end if
       if(abs(d).gt.max(abs(difLo),abs(difHi))) goto 302
       difLo=d
    end if
    if(Elow.eq.Ehigh) goto 302
    goto 301

300 continue
    if(dbg) write(*,*)'--> Elo,Ehi=',Elow,Ehigh,'  d=',difLo,difHi,d
    return
302 niter=0
    do while ((difLo*difHi).gt.0 .and. niter.ne.niterMax)
       Elow=max(Emini,Elow*facLo)
       Ehigh=Ehigh*facHi
       difLo=diff_func(Elow ,E_in ,te,tz)
       difHI=diff_func(Ehigh,E_in ,te,tz)
       niter=niter+1
    end do
    if(niter.ge.niterMax) goto 200
    goto 300
    ! no bracket found, abort
200 write(*,201) E_in,ro
201 format(//,'-R- bracket_EE cannot find bracket for Etot,ro=' &
         ,1p,2e13.4,//)
    stop '-R- cannot find bracket'
  end subroutine bracket_EE

  !-------

  subroutine setErad(eg_o_bg,hnug,ng)
    !Subroutine re-assigns 
    use CRASH_M_projE,only : prep_projE
    implicit none
    integer,intent(IN) :: ng                 !The dimension of array EOverB
    real,dimension(0:ng),intent(in) :: hnug  !The group boundaries, in eV
    real,dimension(ng),intent(in) :: eg_o_bg !Input array of EOverB
    integer :: i

    if(ng.ne.ng_rad) then
       call prep_projE(hnug,ng)
    end if
    ng_rad=ng
    hnu_rad(0:ng_rad)=hnug(0:ng)
    EoB(1:ng_Rad)=eg_o_bg(1:ng_Rad)

  end subroutine setErad
  !-------

  !subroutine set_RO_Ni(atoMass,RO_in,Natom)
    !
    !  This routine can be used in place of passing argument RO or NI  to NLTE_EOS
    ! with no arguments, this routine will reset  NIgiven , ROgiven
    ! implicit none
    ! real,optional, intent(IN) :: RO_in,Natom,atoMass

   ! NIgiven=.false.
   ! ROgiven=.false.

    !if(present(Natom)) then
    !   Ni=Natom
    !   if(present(RO_in)) goto 102
    !   ro=(Ni/avogadro)*atoMass
    !   ROgiven=.true.
    !else
    !   if(.not.present(RO_in)) goto 102
    !   ro=RO_in
    !   Ni=avogadro*(ro/atoMass)
    !   NIgiven=.true.
    !end if

    !return
!102 if(.not.present(atoMass))return
 !   write(*,*)'-P- set_RO_Ni require  Natom .xor. RO_in'
  !  stop '-P- set_RO_Ni.102'
  !end subroutine set_RO_Ni

  !-------
  ! \
end module CRASH_M_NLTE
! /
!-------



!------
function zbrentEE(diff_func,x1,x2,xmin,tol,smallF,dLo,dHi ,EE_in ,te,tz)
  !
  !- This routine was taken from Numerical Recipes. it uses the Brent
  !- method to find a zero of a function without given derivative.
  !- we use here this function to iterate on X till 
  !- "diff_func"  variation < tol  
  !- or |value of function| < smallF		! not in the original Numerical Recipes subroutine

  !
  use CRASH_M_localProperties,only : ro,Ni
  implicit none
  !
  real,intent(IN) :: x1,x2,xmin,tol,smallF,EE_in ,dLo,dHi
  !
  real :: zbrentEE
  real,external :: diff_func
  real ::  te,tz
  !
  integer,parameter ::  itmax=100
  integer ::iter
  real,parameter :: epsZbrent=3.D-8 &
       , tenth=0.1d0,zero=0,half=0.4d0,one=1d0,two=2d0,three=3d0,x1o2=0.5d0
  !
  real :: a,b,c,d,e,fa,fb,fc,tenn
  real :: p,q,r,s,xm,tol1
  real :: TEold,TZold ! ,xt
  integer :: nbZbrent=0
  !
  zbrentEE=(x1+x2)*half	! 080715
  !
  a=x1
  b=x2
  if(a.eq.b) then	
     a=min(a,b)-0.1	
     b=a+0.1	
  end if
  fa=dLo			! avoiding computing two times this point
21 if(fa.eq.zero) then	! 100629
     zbrentEE=a		! 100629
     return			! 100629
  end if			! 100629
  fb=dHi
22 if(fb.eq.zero) then	! 100629
     zbrentEE=b		! 100629
     return			! 100629
  end if			! 100629
1010 iter=0			! 100629
2 if((fb*fa).gt.zero)then	! 100629
     write(*,*)'-R- zbrentEE: root must be bracketed'
     write(*,*)'    x1,x2=',x1,x2,' f1,f2=',fa,fb
     stop 'stop  zbrentEE .not. bracketed'
  end if			! 100629
  !
100 format(5x,'n = ',i4,/,5x,'x1 = ',e12.5,/,5x,'x2 = ',e12.5 &
       ,/,5x,'f1 = ',e12.5,/,5x,'f2 = ',e12.5)
  !
  iter=0	! 051025
  !-  extend range if required
10 tenn=max(tenth*abs(a),tenth)
  if((fa*fb).gt.0) then
     if(iter.ge.itmax)then	! 051025
        goto 25
     end if	! 051025
     iter=iter+1
     if(a.lt.b)then
        a=max(xmin,a-tenn)
        b=b+tenn
     else
        b=max(xmin,b-tenn)
        a=a+tenn
     end if
     fb=diff_func(b,EE_in ,te,tz)
     fa=diff_func(a,EE_in ,te,tz)
     goto 10
  end if	! (fb*fa).gt.zero
  !
  !-  use  Brent algorithm
20 c=a
  fc=fa
  TEold=-1
  TZold=-1
  e = b - a
  d = e
  LOOP11 : do iter=1,itmax
     if((fb*fc).gt.zero) then	! at iter#1 fb=fc
        c=a
        fc=fa
        d=b-a
        e=d
     end if
     if(abs(fc).lt.abs(fb)) then	! at iter#1 fb=fc
        a=b
        b=c
        c=a
        fa=fb
        fb=fc
        fc=fa
     end if
     !	 tol1=two*epsZbrent*abs(b)+half*tol	! original writing
     tol1=(abs(c)+abs(b))*epsZbrent	! proposed change
     xm=half*(c-b)	! note as tol1 is prop. to |b| ....

     TEold=te
     if(abs(xm).le.tol1)then
        zbrentEE=b
        return
        
     elseif(abs(fb).le.smallF) then		! mb
        zbrentEE=b			! mb
        return
    
     elseif(abs(fa).le.smallF) then		! mb
        zbrentEE=a			! mb
        return
     end if

     if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
        s=fb/fa
        if(a.eq.c) then
           p=2*xm*s
           q=one-s
        else
           q=fa/fc
           r=fb/fc
           p=s*(2*xm*q*(q-r)-(b-a)*(r-one))
           q=(q-one)*(r-one)*(s-one)
        end if
        if(p.gt.zero) then
           q=-q
        elseif(p.le.zero) then
           p=-p
        end if
        if(2*p .lt. min(3*xm*q-abs(tol1*q),abs(e*q))) then
           e=d

           d=p/q
        else
           d=xm
           e=d
        end if
     else	! (abs(e).ge.tol1 ...
        d=xm
        e=d
     end if	! (abs(e).ge.tol1 ...
     a=b
     fa=fb
     if(abs(d) .gt. tol1) then
        b=b+d
     else
        b=b+sign(tol1,xm)
     end if
     b=max(b,xmin)		! mb.080715
     fb=diff_func(b,EE_in ,te,tz)
  end do LOOP11
  !
  write(* ,111) itmax
111 format('-R- zbrentEE exceeding maximum iterations' &
       ,',(=',i5,')')
  stop '-R- no convergence in zbrentEE'
  zbrentEE=b
  return
  !- out of range
25 write(* ,112) nbZbrent,x1,x2,tol,a,b,fa,fb,c,fc
112 format('-R- from ',i8,'# zbrentEE' &
       ,' function not bracketed ' &
       ,/,'x1,x2,tol,lfix=',1p,3e12.3,l2 &
       ,/,' a,b,fa,fb=',1p,4e12.3,' c,fc=',2e12.3)
  stop '-R- no bracketing in zbrentEE'

end function zbrentEE

!----------      ***  eediff.f90 ***
function rootEE(diff_func,x1,x2,xmin,tol,smallF,dLo,dHi ,EE_in ,te,tz)
  !
  use CRASH_M_localProperties,only : Efloor  , ro,Ni
  use CRASH_M_NLTE,only : dbg,useEElog
  implicit none
  !
  real,intent(IN) :: EE_in
  real,intent(IN) :: x1,x2,xmin,tol,smallF
  real,intent(OUT) ::  te,tz
  real ::  dLo,dHI,d
  !
  real :: rootEE
  real,external :: diff_func
  real :: eLo,eHi,e
  !
  integer,parameter ::  itmax=100
  integer :: iter
  real,parameter :: eps=1d-4
  integer :: i=0,maxiter=100

  if(dbg) write(*,*)' rooEE : x1,x2=',x1,x2,'   dLo,dHI=',dLo,dHi

  eLo=x1
  !	dLo=diff_func(eLo ,EE_in ,te,tz)	! done outside
  rootEE=eLo
  if(abs(dLo).le.eps) return
  eHi=x2
  !	dHi=diff_func(eHi ,EE_in ,te,tz)	! done outside
  if(dhi.eq.0) stop '#741.  rootEE found dHi=0'
  rootEE=eHi
  if(abs(dHi).le.eps) return

  if((dlo*dHi).gt.0) goto 100
  if(eHi.le.Efloor .or. abs(dhi).gt.0.99) then
     e=Ehi
     e=sqrt(e*Elo)
     PUSHhi:	 do while (e.gt.Efloor) 
        d=diff_func(e ,EE_in ,te,tz)		! PUSHhi : 110827-22:31
        if(abs(d).le.smallF) then
	   rootEE=e
	   return
        end if
        if((d*dLo).le.0) then
	   eHi=e
	   dHi=d
           exit PUSHhi
        end if
        e=sqrt(e*Elo)
     end do PUSHhi
  end if
888 do iter=1,maxiter
     ! overshooting secant
     if(abs(dHi).gt.abs(dLo))then
        e=(eLo*dHi -eHi*2*dLo)/(dHi-2*dLo)
     else
        e=(eLo*2*dHi -eHi*dLo)/(2*dHi-dLo)
     end if
     if(dbg) write(*,*)'iter#',iter,' eLo,eHi=',eLo,eHi,'  dLo,dHi=',dLo,dHi,' -> e=',e
     d=diff_func(e ,EE_in ,te,tz)
     if(dbg) write(*,*)'iter#',iter,' eLo,eHi,e=',eLo,eHi,e,'  dLo,dHi,d=',dLo,dHi,d,'  smallF=',smallF
     if(abs(d).le.smallF) then
        rootEE=e
        return
     elseif((dLo*d).lt.0) then
        eHi=e
        dHi=d
        if(dbg) write(*,*)'755.rootEE:EEdiff(',ehi,')=',dhi,' te,tz=',te,tz
        rootEE=eHi
        if(abs(dHi).le.smallF .or.iter.eq.maxiter) return
     else
        eLo=e
        dLo=d
        rootEE=eLo
        if(abs(dLo).le.smallF .or.iter.eq.maxiter) return
     end if
  end do
  return
  !
100 write(*,101) eLo,eHi,ro,EE_in,dLo,dHi
101 format(/,'-P- rootEE not bracketed: Eeff_Lo,Eeff_eHi=',1p,2e13.4, &
       ' for ro,Etot=',2e13.4,' diff=',2e13.3,//)
  stop '-P- rootEE not bracketed:'
end function rootEE

!----------      ***  eediff.f90 ***

function EEdiff(Eeff ,EE_in ,te,Tz_new)

  ! given  ro,Etot,Eeff, derive Te,Tz,zbar, then Eeff_new  and  compute the difference  EEdiff=Eeff_new-Eeff_old
  ! The function is used to find the "effective energy", i.e. internal energy minus the non-LTE correction
  !  3/2 * kBr_E * (Te-Tz).  kBr_E is the conversion factor from temperature to kinetic energy
  ! It may or may not include the bulk density, depending on the EOS in use.

  use CRASH_M_localProperties,only : Ni &
       ,kBr_E,kBr_P &	! conversion factor T -> E [erg/g], or P [dyne], includes "ro"
                                ! they have to be set for the current ro !
       ,zion	&	! 0 for elec. only,  1 for elec.+ion.
       ,Zsmall &
                                ! it has to be set for the actual mode (2T or 1T) of the code
       ,Efloor,TEfloor	! lower limits of EOS  or of NLTE domain
  ! E < Efloor or T < TEfloor, LTE will be assumed
  use CRASH_M_NLTE , only : EoB,ng_rad  ,dbg	! radiative enegy density (for all groups) and Planckian function
  use M_radiom, only : caltz0
  ! the calling sequence have to fill them
  implicit none
  real :: EEdiff
  real,intent(IN) :: EE_in	! total energy
  real,intent(IN) :: Eeff	! Effective energy = Etot -3/2*kb*ro*(Te-Tz)
  real,parameter :: x_3o2=1.5d0
  real :: zbar,te,tz,ne,pe,cv,zp,Tz_new,EE_new
  !
  call LTE_EOS_inv(tz,Eeff,Pe,Zbar,Cv)	! ro_in

  !  small Zbar would yield diverging the "EEeff equ." (see correctEOS + EEdiff)
  zp=max(zbar,Zsmall)+zion
  !
  if(dbg) write(*,*)'EEdiff.1102: EE_in,Eeff, kBr_E,zp,tz=',EE_in,Eeff, kBr_E,zp,tz
  te=(EE_in-Eeff) / (x_3o2 * kBr_E * zp) +tz
  if(dbg) write(*,*)'EEdiff: tz=',tz,' -> te=',te,'  zbar,zion=',zbar,zion
  if(te.gt.1e9) write(*,*)'EEdiff found  Te =',te,' > 1e9'
  te=min(te,1000*Tz)
  te=max(Te,TEfloor)
  !	Ni=avogadro*ro/atomass
  ne=zbar*Ni
  if(dbg) write(*,*)'ne,zbar,ni,te=',ne,zbar,ni,te
  if(ne.eq.0) stop '-P- EEdif: found  Ne=0'
  call calTZ0(Te,Ne, Tz_new, EoB)
  if(dbg) write(*,*)'EEdiff calls calTz, te,ne=',te,ne,' -> TZ_new=',tz_new
  zp=max(zbar,Zsmall)+zion
  EE_new=EE_in - x_3o2 * kBr_E * ( Te-Tz_new)*zp
  EEdiff=(EE_new-Eeff)/(abs(EE_new)+abs(Eeff))
  if(dbg) write(*,*)'EEdiff:  -> EE_new=',EE_new,'   (Eeff=',Eeff,')  diff=',EEdiff &
       ,'  te,tz=',Te,Tz_new
  !
  return
end function EEdiff


!----------      ***  eediff.f90 ***

function EEdiff_ln(lnEeff ,EE_in ,te,Tz_new)

  ! same as EEdiff, but the input arg. lnEeff is log(Eeff)

  use CRASH_M_expTab
  use CRASH_M_localProperties,only : Ni &
       ,kBr_E,kBr_P &	! conversion factor T -> E [erg/g], or P [dyne], includes "ro"
                                ! they have to be set for the current ro !
       ,zion	&	! 0 for elec. only,  1 for elec.+ion.
       ,Zsmall &
                                ! it has to be set for the actual mode (2T or 1T) of the code
       ,Efloor,TEfloor	! lower limits of EOS  or of NLTE domain
  ! E < Efloor or T < TEfloor, LTE will be assumed
  use CRASH_M_NLTE , only : EoB,ng_rad  	! radiative enegy density (for all groups) and Planckian function
  use M_radiom, only : caltz0
  ! the calling sequence have to fill them
  implicit none
  real :: EEdiff_ln
  real,intent(IN) :: EE_in	! total energy
  real,intent(IN) :: lnEeff	! Effective energy = Etot -3/2*kb*ro*(Te-Tz)
  real:: Eeff
  real,parameter :: x_3o2=1.5d0
  real :: zbar,te,tz,ne,pe,cv,zp,Tz_new,EE_new
  real :: y,ey    !for exp_tabm.h
  !
  y=-lnEeff
  include 'exp_tabm.h'
  Eeff=ey
  !
  call LTE_EOS_inv(tz,Eeff,Pe,Zbar,Cv)	! ro_in
  !no	if(Zbar .le. Zsmall) then
  !no	 EEdiff_ln=0	! this would lead to artificial convergence
  !no	 return
  !no	end if

  !  small Zbar would yield diverging the "EEeff equ." (see correctEOS + EEdiff)
  zp=max(zbar,Zsmall)+zion
  !
  te=(EE_in-Eeff) / (x_3o2 * kBr_E * zp) +tz
  te=max(Te,TEfloor)
  !	Ni=avogadro*ro/atomass
  ne=zbar*Ni
  !xx		if(dbg) write(*,*)'ne,zbar,ni=',ne,zbar,ni
  if(ne.eq.0) stop 'EEdif_ln: NE=0'
  call calTZ0(Te,Ne, Tz_new, EoB)
  EE_new=EE_in - x_3o2 * kBr_E * ( Te-Tz_new)*zp
  EEdiff_ln=(EE_new-Eeff)/(abs(EE_new)+abs(Eeff))
  !
  return
end function EEdiff_ln

