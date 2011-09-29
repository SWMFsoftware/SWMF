 ! *** mod_local.f90 ***
 !------
module M_localProperties
	real,save :: zion=1 	! 0 if code use 2 T (Te + Ti)
	real :: &
	 atoNum,Atomass,roSolid,Efloor=1d-5,TEfloor=1./40.,Pcold=0	&
 ! density dependant , can be passed through MODULE:
	,kBr_E,kBr_P  &	! conversion factor T -> E [erg/g], or P [dyne], includes "ro"
	,ro &		! bulk density [g/cm3] shared by all routines, and unchanged
	,Ni    		! ionic density [cm-3]
	real,parameter :: kB_ztf_E=1.8591817e-8,kB_ztf_P=kb_ztf_E
	real,parameter :: ErgPEReV=1.60219e-12,DYNEperEV=ErgPEReV 
	real,parameter :: avogadro=6.02e23,EVperK = 11604.1
	
 ! flags set by CALL setRoNi
	logical,save :: roGiven=.false.,niGiven=.false.

 !  Eflor, TEfloor may vary with the EOS in use. Has to be set by the calling program
 !  kBr_E, kbR_P  is dependant on the units used.  kinetic contribution to Energy and pressure
 !  they are generally proportional to the actual bulk density  ro
 !  are   3/2*kBr_E * Zp * T  and kBr_P * Zp * T  ,  
 !  zp =Zbar if only elec. contrib,  zp=Zbar+1   with  elec+ion contib. (w/ Te=Ti) 

	real,parameter :: Zsmall=0.05d0

 !  Zsmall is a lower bound to be used in the  Eef eq. (see correctEOS and EEdiff)

 !------
end module M_localProperties

subroutine set_ZA(Z,A)
  use M_localProperties,only : atoNum,atoMass
  implicit none
  real :: Z,A
  atoNum=Z
  atoMass=A
  return
end subroutine set_ZA
