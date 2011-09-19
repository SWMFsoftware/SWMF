 !-  testNLTE.f90 -  in progress
 !	testing NLTE_EOS, using zBrent for  te,tz finding given ro,Etot & EOS=ZTF_Perrot
 !      compile with :  radiom.f90, eos_material.f90
 
 !- testZTF.f90 - checked
 !      testing  eos_material (ZTF_EOS_dir, ZTF_EOS_inv & ZTF*s
 !	compule with : eos_material.f90, testZTF.f90
 
 !- testradiom.f90 - checked
 !	testing RADIOM routines : CALTZ, CALTE, 
 !  f90 -q -O -z 2 -en -Z899 -Z124 -Z938 -Z1553 -o tn -g testNLTE.f90 radiom.o
 !	

program testNLTE
  use M_expTab,only: exp_tab8
  use M_localProperties
  use M_RADIOM, only : prep_projE,prepCorrUbar &
       ,printVersion
  use M_NLTE,only : nlte_eos ,useLTE ,ng,Erad,Brad
  implicit none
  real :: Ni_l,Te_l,Ee_l,Tz_l &
       ,zbar_l, pe_l, Cv_l
  !
  real,dimension(0:45),parameter :: gr45=(/ &
	    1., 10., 48.75, 50.5, 58.5, 67.5 ,68.5, 79.25, 91., 98.5, 105. &
	  , 127.5 ,150., 167.33, 184.67, 202., 219.33, 236.67 ,254., 267., 287., 300. &
	  , 308., 309.5 ,320., 335., 345., 361., 364., 369. ,376., 378., 403.5, 435.5 &
	  , 436.5, 446.2 ,474., 480., 540., 600., 757., 950. ,1250., 2500. &
	  , 12000., 100000. &
	  /)
  real,dimension(301) :: ones =1.d0
  integer :: ir,it
  
	
	
 
  call printVersion()
  !\
  ! Tabulates the exponential function.
  !/
  call exp_tab8()

  !\
  ! Coefficients for tranforming from the user defined grid to
  ! the refined logrithmic-unifrom internal fixed grid
  !/ 
  call prep_projE(gr45,45)
  !
  ng=45

  !\
  ! Initialize and calculate some internal arrays
  !/
  call prepCorrUbar()
  !===================Set material====================
  atoNum=79.
  atoMass=179.
  roSolid=19.3
  Efloor=1e-5
  Pcold=0
  !===================================================
  !Check the consitency of the input parameters
  call verify()
  !==================================================
  !
  !Set the limit for atomic concentration   
  LOOPro: do ir= -6,2	    ! -4,-4	!
     ! Atomic particles per cm3
     Ni_l=(3*10d0**ir)*avogadro/atomass  ![g/cm3]*[particles/mole]/[g/mole]!
     LOOPte:	do  it= 1,5 ! 2,2	!
        !For each density initialize the temperature
        Te_l=10d0**it
        print*
        !\\\\\\\\\\\\
        !============ Direct EOS, LTE
        !////////////
        write(*,*)'- direct EOS , LTE -'
        useLTE=.true.
        Ee_l=1e12
        print*
        write(*,*)'LTE,ni+te=',Ni_l,Te_l
        Ee_l=-1;tz_l=-1;zbar_l=-1;pe_l=-1;cv_l=-1	! to check
        call NLTE_EOS(Natom=Ni_l, Te_in=Te_l, Ee_out=Ee_l, Tz_out=tz_l &
             ,Zbar_out=zbar_l, Pe_out=pe_l, Cv_out=Cv_l)
        write(*,*)' -> Ee,Tz=',Ee_l,Tz_l
        write(*,*)' -> Zbar=',zbar_l,' Pe=',Pe_l,' Cv=',Cv_l
        
        print*
        !\\\\\\\\\\\\
        !============ Inverse EOS, LTE
        !////////////
        write(*,*)'- inverse EOS , LTE-'
        write(*,*)'LTE,ni+Ee=',Ni_l,Ee_l
        te_l=-1;tz_l=-1;zbar_l=-1;pe_l=-1;cv_l=-1	! to check
        call NLTE_EOS(Natom=Ni_l, EE_IN=Ee_l, Tz_out=tz_l &
             ,TE_out=te_l, Zbar_out=zbar_l, Pe_out=pe_l, Cv_out=Cv_l)
        write(*,*)' -> Te,Tz=',Te_l,Tz_l
        write(*,*)' -> Zbar=',zbar_l,' Pe=',pe_l,' Cv=',Cv_l
        !\\\\\\\\\\\\\
        !============= Direct EOS, Non LTE 
        !/////////////
        !As the first step, define, Erad, Brad arrays
        !Where Erad are the group-integrated radiation energy, in arbitrary units,
        !and Brad are the integrals of the Plankian spectrum, in the same uinits,
        !with the temperature of radiation being equal to the electron temperature,
        !at the beginning of the time step: even if the EOS is applied not in the
        !starting time for this timestep.
        Erad = 0.0; Brad = 1.0
        print*
        write(*,*)'- ---------- -'
        print*
        write(*,*)'- direct EOS , nonLTE -'
        write(*,*)'nonLTE,ni+te=',Ni_l,Te_l
        Ni_l=3*10d0**ir*avogadro/atoMass
        Te_l=10d0**it
        useLTE=.false.
        Ee_l=-1;tz_l=-1;zbar_l=-1;pe_l=-1;cv_l=-1	! to check
        call NLTE_EOS(Natom=Ni_l,Te_in=Te_l,Ee_out=Ee_l,Tz_out=tz_l &
             ,Zbar_out=zbar_l, Pe_out=pe_l, Cv_out=Cv_l)
        write(*,*)' -> Ee,Tz=',Ee_l,Tz_l
        write(*,*)' -> Zbar=',zbar_l,' Pe=',pe_l,' Cv=',Cv_l
        !
        print*
        !\\\\\\\\\\\\\
        !============= Inverse EOS, Non-LTE 
        !/////////////
        write(*,*)'- inverse EOS , nonLTE (not completed) -'
        print*
        write(*,*)'nonLTE,ni+Ee=',Ni_l,Ee_l
        te_l=-1; tz_l=-1; zbar_l=-1; pe_l=-1; cv_l=-1		! to check
        call NLTE_EOS(Natom=Ni_l, Ee_in=Ee_l, Te_out=Te_l, Tz_out=tz_l &
             ,Zbar_out=zbar_l, Pe_out=pe_l, Cv_out=Cv_l)
        write(*,*)' -> Ee,Tz=',Ee_l,Tz_l
        write(*,*)' -> Zbar=',zbar_l,' Pe=',pe_l,' Cv=',Cv_l
        
     enddo LOOPte
  enddo LOOPro

 !
  stop
end program testNLTE

 !------
	subroutine verify()	! check transmission of Z,A,...
	use M_localProperties
	implicit none
	write(*,*)'verify: Z,A=',atonum,atomass,' roS=',roSolid
	if(atonum.gt.0  .and. atomass.gt.0) return
	print *,'-P-  atoNum or atoMass not given'
	stop
	end subroutine verify

!============================================================================
! The following subroutines are here so that we can use SWMF library routines
! Also some features available in SWMF mode only require empty subroutines
! for compilation of the stand alone code.
!============================================================================
subroutine CON_stop(StringError)
  implicit none
  character (len=*), intent(in) :: StringError
  write(*,*)StringError
  stop
end subroutine CON_stop
!============================================================================
subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe
end subroutine CON_set_do_test
