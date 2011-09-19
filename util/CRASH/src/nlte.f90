 !
 ! nlte.f90
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
 !	call nLTE_EOS(...
 !	call correctEOS(...
 !	call correctKnuJnu(...
 !
 !-----------
 !------
module M_localProperties
  !This module is used to
  !1. Keep as global variables the input and output
  !parameters for the model equations to be solved as 
  !well as some general physical constants and the conversion coefficients.
  !
  !Real flag to distiguish between the cases of two temperatures or a single 
  !temperature. With this flag the multipliers like (Zbar) or (ZBar+1),
  !for a single temperature or for two temperatures correspondingly,
  !may be coded in a genereic way like Zbar+Zion.
  real,save :: zion=1 	! 0 if code use 2 T (Te + Ti)
  !\\\\\\\\
  !Attention!!!!!!!!!!!!
  !////////
  !All possible options for the equation to be solved within
  !the framework of the module M_NLTE, which allow for different units
  !for tenperature and energy, as well as the choice for the energy to be related
  !to the unit of volume, of mass, or per an atomic cell, are all included into KBr_E
  !In case the equation is solved for pressure (not for the energy density the same 
  !consideration is applicable to KBr_P.
  !////////   
  real ::   &
	 atoNum,    &   !Atomic number
         Atomass,   &   !Atomic mass
         roSolid,   &   !Solid density (the input parameter in the Thomas-Fermi EOS
         Efloor=1d-5, & !The minimal internal energy density, in the same units as E_tot
         Pcold=0    &	!The conventional value for the pressure in the unshocked solid.
         ,kBr_E,kBr_P  &	! conversion factor T -> E [erg/g], or P [dyne], includes "ro"
         ,ro &		! bulk density [g/cm3] shared by all routines, and unchanged
         ,Te,Ne &	! electronic temperature [eV] and density [cm-3]
         ,Etot,Ptot &	! internal energy & pressure (e+i), if (2T) this should be modified
         ,zt		! work value for "ZTF_dif", == zbar*te
  real,parameter ::  ErgPEReV=1.8591817e-8,DYNEperEV=ErgPEReV &
       ,avogadro=6.02e23
  !------
end module M_localProperties

!\ 
module M_NLTE
  ! /
  !This module solves the equation for the effective energy.
  !E_{eff} = E_{tot}^{NLTE} - (3/2) k_B x rho x Z^*(T_Z) * (T_e-T_z)
  !where E_{tot} is the input parameter for the "inverse" equation of state
  !(from which the temperatures should be solved) and E_{eff} is the tabulated
  !LTE equation of state with (rho, T_Z) being the input parameters for the 
  !latter. 
  use M_RADIOM,only : caltz0,calte0,caltz,calte
  use M_projE,only : mxOut
  use M_localProperties
  implicit none
  
  logical :: useLTE=.false.
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
  !  useLTE : .T. to save CPU (ex. for H at high temperature)
  !
  integer,save :: ng
  real,dimension(mxout) :: Erad,Brad
  ! parameter for LTE_EOS_dir  convergence on 'Zbar'
  real,save :: epsD=1d-3
  integer,save :: niterMax=10
  !
  !================Interface to the equation of state========
  interface
     subroutine  LTE_EOS_inv(te,Etot,Ptot,Zbar,Cv)
       implicit none
       real,intent(IN) :: Etot
       real,intent(OUT) :: te,Ptot,Zbar,Cv
     end subroutine LTE_EOS_inv
  end interface
  !-------
contains
  !-------
  real function TZdif(Tzold )
    real,intent(IN) :: Tzold
    real :: Tz,Ez,TzEOS,Pe,Cv,Zbar
    !
    !- Ne    : electronic density [cm-3]
    !- Natom : atomic density [cm-3]
    !- Te    : electronic temperature in eV
    !
    TZeos=TZold
    call caltz0(Te,Ne, Tz, Erad,Brad)	! Te,Ne, Erad,Brad from M_localProperties
    Ez=max( Etot-kBr_E*(Te-Tz) , Efloor)
    !- waiting for using the real direct EOS
    call LTE_EOS_inv(tz,Ez,Pe,Zbar,Cv)	! ro : in module M_localProperties
    TZdif=Tz-TZeos
    
    return
  end function TZdif
  !-------
  subroutine NLTE_EOS(Natom,RO_in, Te_in,Ee_in,Et_in,Pe_in,Pt_in, estim_Zbar,estim_Tz &
       ,estim_Te,Zbar_out,Tz_out  ,Te_out,Ee_out,Et_out,Pe_out,Pt_out,Cv_out)
    !
    !  to account for  Et=Ee+Ei w/ Te=Ti : use  Ee_* or Et_*
    !
    !  direct/inverse EOS is used depending of presence of  E*_in or Te_in
    ! 
    use M_localProperties,only : ro
    implicit none
    real,optional,intent(IN) :: Natom,ro_in !density [atom/cm3] or [g/cm3] 
    !Herewith, e stands for electrons, t for total
    real,optional,intent(IN) :: Te_in,Ee_in,Et_in,Pe_in,Pt_in ! one and only one
    !Initial guess
    real,optional,intent(IN) :: estim_Zbar,estim_Tz,estim_Te
    real,optional,intent(OUT) :: Zbar_out,Tz_out
    real,optional,intent(OUT) :: Te_out,Ee_out,Et_out,Pe_out,Pt_out,Cv_out
   
    real :: Ni,Ez,Ee,Pe,Ne,zdt,z_old,tz_old,te_old,d
    real :: Te,Tz,zbar,zp,Cv
    integer :: niter
    real,parameter :: x_0=0,x_3o2=1.5d0,x_1=1d0,x_1o2=0.5d0,two=2d0,three=3d0
    !------------------
    !Initialize variables
    zbar=0
    tz=0
    te=0
    ee=0
    pe=0
    Cv=0
    ne=0
    if(present(Natom)) then
       Ni=Natom
       if(present(RO_in)) goto 102
       ro=(Ni/avogadro)*atoMass
    else
       if(.not.present(RO_in)) goto 102
       ro=RO_in
       ni=avogadro*(ro/atoMass)
    endif
    !--
    !- direct EOS -
    if(present(Te_in)) then
       !--
       if(present(Ee_in)) goto 101
       ! Te=Ti in force or not ?
       zion=1					! 110928
       if(present(Ee_out).or.present(Pe_Out) ) zion=0	! 110928
       te=Te_in
       if(useLTE) then
          !- waiting for using the real direct EOS
	  call LTE_EOS_dir(te,Ee,Pe,Zbar,Cv)	! ro : in module M_localProperties
          tz=te
       else
          if(present(estim_Zbar)) then
             zbar=estim_Zbar
          else
             zbar=atoNum*0.5
          endif
          ! iterates over 'Zbar'  >>
          niter=0
111       niter=niter+1
          Ne=Zbar*Ni
          call calTz0(Te,Ne, Tz, Erad,Brad)
          !- waiting for using the real direct EOS
          call LTE_EOS_dir(tz,Ee,Pe,Zbar,Cv)	! ro : in module M_localProperties
          d=(Zbar*Ni-Ne)/(Zbar*Ni+Ne)
          if(abs(d).gt.epsD.and. niter.lt.niterMax) goto 111
          ! iterates over 'Zbar' <<
          call correctEOS(zbar , Te, Tz ,Ee,Pe,Cv)
       endif
       !--
    else	! if(present(TE_in))
       !--
       !- inverse EOS -
       if(.not.(present(EE_in).or.present(Et_in))) goto 101
       ! Te=Ti in force or not ?
       if(present(Ee_in) ) then
	  zion=0		! 110928
	  ee=Ee_in
       else
	  zion=1
	  ee=Et_in
       endif

       if(useLTE) then
          !- waiting for using the real inverse EOS
	  call LTE_EOS_inv(te,ee,pe,Zbar,Cv)	! ro, {Erad,Brad} : in module M_localProperties
	  Ne=Zbar*Natom
	  tz=te
       else

          !- inverse EOS, non LTE -
          !  should use bracket(Tzdiff,tz1,tz2) + zbrent(tzDif,tz1,tz2)  <<<<<<<
          !   	ro= ; Erad= ; Brad= 
          !   T	e= ; Ne=  : update at each iteration ? use Cv_loc ?

          ! estimates  Tz,zbar,Ne to start iterations (not complete) ....

          call LTE_EOS_inv(tz,Ee_in,pe,Zbar,Cv)

          write(*,*)'check vs flowChart'
	  if(present(estim_Tz)) then
             tz=estim_Tz
             if(present(estim_zbar)) then
                zbar=estim_Zbar
             else
                call LTE_EOS_dir(tz,ee,pe,Zbar,cv)
             endif
             if(present(estim_Te)) then
                te=estim_Te
             else
                Te=2*Tz
             endif
	  elseif(present(estim_Te)) then
             tz=estim_Te*0.5
             if(present(estim_zbar)) then
                zbar=estim_Zbar
             else
                call LTE_EOS_dir(tz,ee,pe,Zbar,cv)
             endif
	  else
             call LTE_EOS_inv(tz,ee,pe,Zbar,Cv)	! ro, {Erad,Brad} : in module M_localProperties
             Te=1.5*Tz
             Tz=0.75*Tz
	  endif
	  ne=zbar*Ni
112	  z_old=zbar
	  tz_old=tz
	  te_old=te
	  zp=zbar+zion
	  zdt=zp*(Te-Tz)
	  Ee=Ee_in
	  Ez=Ee- x_3o2 * kBr_E * zdt
	  call LTE_EOS_inv(tz,Ez,pe,Zbar,Cv)	! ro, {Erad,Brad} : in module M_localProperties
       endif
       !--
    endif 	! if(present(TE_in))
    !--
    !x		write(*,*)'zbar,tz,te,ee,pe,Cv=',zbar,tz,te,ee,pe,Cv
    if(present(zbar_out)) zbar_out=zbar
    !x		write(*,*)'tz=',tz
    if(present(Tz_out)) Tz_out=Tz
    !x		write(*,*)'te=',te
    if(present(Te_out)) Te_out=te
    !x		write(*,*)'ee=',ee
    if(present(Ee_out)) Ee_out=ee
    if(present(Et_out)) Et_out=ee
    !x		write(*,*)'pe=',pe
    if(present(Pe_out)) Pe_out=pe
    if(present(Pt_out)) Pt_out=pe
    !x		write(*,*)'Cv=',Cv
    if(present(Cv_out)) Cv_out=Cv
    !x		write(*,*)'LTE_EOS  returns'
    return
102 print *,'-P- LTE_EOS require  Natom .xor. RO_in'
    return
101 print *,'-P- LTE_EOS require  TE_in .xor. EE_in'
    return
    !
  end subroutine NLTE_EOS
 !-------
	subroutine correctEOS(zbar , Te,Tz &
		,Ee,Pe,Cv)
 ! Te,Tz, kBro_E,kBro_P,Etot,Ptot : thorugh module variables
	implicit none
	real,intent(IN) :: zbar,Te,Tz
	real,intent(INOUT) :: Ee,Pe
	real,optional,intent(INOUT) :: Cv
	real :: zdt,zp
	real,parameter :: x_1=1d0,x_3o2=3d0/2d0
 ! kbRO is the proportionnality factor in the code units (generally = Boltzmann cst * density)
 !   so  3/2* kbRo *Te  is the kinetic (translational) energy at Te

	if(Te.le.0) return
	zp=(zbar+zion)
	zdt=zp * ( Te-Tz)
	Ee=Ee + x_3o2 * kBr_E * zdt
	Pe=Pe +         kBr_P * zdt
	if(present(Cv)) then
	 Cv=Cv+ x_3o2 * kBr_E * zp
	endif
	return
	end subroutine correctEOS
 !-------
 ! \
	end module M_NLTE
 ! /
 !-------
