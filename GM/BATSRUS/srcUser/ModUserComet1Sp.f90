!^CFG COPYRIGHT UM
!========================================================================
module ModUser

  use ModSize
  use ModUserEmpty,               &
       IMPLEMENTED1 => user_set_ics,                    &
       IMPLEMENTED2 => user_read_inputs,                &
       IMPLEMENTED3 => user_calc_sources,               &
       IMPLEMENTED4 => user_update_states,              &
       IMPLEMENTED5 => user_init_point_implicit

  include 'user_module.h' !list of public methods

  !\
  ! Here you must define a user routine Version number and a 
  ! descriptive string.
  !/
  real,              parameter :: VersionUserModule = 2.01
  character (len=*), parameter :: NameUserModule = &
       'Yingdong single species cometary MHD module, M. Rubin & K.C. Hansen, Feb 2008'

  real, public, dimension(0:nI+1, 0:nJ+1, 0:nK+1,nBLK, 4) :: Neutral_BLK

  integer :: NR_neutral=121, NTheta_neutral=121, NPhi_neutral=121
  integer, parameter :: NR_ = 4, NTh_ = 5, NPhi_ = 6
  real*8, dimension(61, 61, 61, 10) :: NeutralN
  logical :: ReadNeutral = .False. , DoInitialize = .True.
  character*16 :: NeutralFile='test_3d.dat'
  character (len=*), parameter :: Name='user_heat_source'

  real ::  jet_ln2, Qprod_day, Qprod_nit &
       , Qprod_jet, Qprod_jeta
  logical :: DoTest, DoTestMe

  logical ::  UseMultiSpecies=.false.
  integer ::  nSpecies=1
  integer, parameter :: MaxSpecies=2, MaxNuSpecies=1,  &
      MaxReactions=3
  integer ::  nNuSpecies=1

  integer, parameter :: &       !reaction number
    H2O_hv__H2Op_em_=1 ,&       !H2O+hv-->H2Op+em: photoinozation
    H2Op_em__neutral_=2   ,&    !H2Op+em-->neutral: recombination
    Hp_H2O__H2Op_H_=3           !Hp+H2O-->H2Op+H: change exchange

  real, dimension(MaxReactions) :: Rate_I
  real, dimension(MaxReactions) :: &
       Ratedim_I=(/ 5.42e-7, 7.0e-7, 8.2e-9 /)  !cm^3 s^(-1)

  integer, parameter :: &	! order of ion species
       Hp_   =1, &
       H2Op_ =2

  character (len=10), dimension(MaxSpecies):: &
     ion_name_I=(/'Hp   ', 'H2Op '/)

  real, dimension(MaxSpecies)::  &
     MassSpecies_I=(/1,18 /)  !atm

  integer, parameter :: & ! order of Neutral species
       H2O_=1

  real, dimension(MaxNuSpecies):: CrossSection_dim_I=(/8.2e-9/),&
       CrossSection_I
  real:: Productrate0,Optdep
  real, dimension(MaxNuSpecies)::  NuMassSpecies_I=(/18/), &
       HNuSpecies_I, BodynDenNuSpecies_I, &  
       BodynDenNuSpdim_I=(/1.e8/)

  integer, parameter :: & ! other numbers
    em_=-1 ,&
    hv_=-2

  ! Define some variables and set defaults
  real :: &  
    kin=1.7E-9,&
    kin_in=1.,&
    mbar=17.,&
    Unr=1.,&
    Unr_in=1.,&
    Qprod=7.E29,&
    ionization_rate=1.E-6

  real :: &
    Qprodd=0.5,&
    Qprodn=0.5,&
    Qprodj=.0,&
    Qprodja=.0,&
    jet_width=25.0,&
    Ajet=.25,&
    Ajet1=.25,&
    Tion=180.0

  integer :: jpattern=0


contains

  !========================================================================
  !  SUBROUTINE user_read_inputs
  !========================================================================
  subroutine user_read_inputs

    use ModMain
    use ModProcMH,    ONLY: iProc
    use ModReadParam
    use ModIO,        ONLY: write_prefix, write_myname, iUnitOut

    integer:: i, j, k
    character (len=100) :: NameCommand, line


    !-------------------------------------------------------------------------

    if(iProc==0.and.lVerbose > 0)then
       call write_prefix; write(iUnitOut,*)'User read_input COMET starts'
    endif

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)

         case("#COMET")
            call read_var('Qprod' ,Qprod)
            call read_var('Unr_in' ,Unr_in)
            call read_var('mbar',mbar)
            call read_var('ionization_rate',ionization_rate)
            call read_var('kin_in',kin_in)
            !\
            ! Convert comet input parameters to SI
            !/
            kin=kin_in*1E-6
            Unr=Unr_in*1E3

         case("#COMETPARAM")
              call read_var('Qprodd' ,Qprod_day)
              call read_var('Qprodn' ,Qprod_nit)
              call read_var('Qprodj' ,Qprod_jet)
              call read_var('Qprodja' ,Qprod_jeta)
              call read_var('jet_width' ,jet_width)
              call read_var('Ajet' ,Ajet)
              call read_var('Ajet1' ,Ajet1)
              call read_var('jpattern' ,jpattern)
              call read_var('Tion', Tion)
              DoInitialize = .True.

       case("#NeutralComa")
         call read_var('readneutral',ReadNeutral)
         if(ReadNeutral)then
           call read_var('neutralfile', NeutralFile)
           open(125, file=NeutralFile, status="old")
           do k=1, 9; read(125,*) line ; enddo
           read(125, '(A3, I6, A4, I6, A4, I6, A40)' ) line, NR_neutral, line, NTheta_neutral, &
                       line, NPhi_neutral, line
           read(125, *) line
           read(125, *) line
           do k = 1, NR_neutral; do j=1, NTheta_neutral;  do i=1, NPhi_neutral
              read(125,*) NeutralN(i,j,k,:)
           end do;  end do;  end do
           close(125)
         end if
       case('#USERINPUTEND')
         if(iProc==0.and.lVerbose > 0)then
             call write_prefix;
             write(iUnitOut,*)'User read_input COMET ends'
         endif
         EXIT
       case default
         if(iProc==0) then
             call write_myname; write(*,*) &
                  'ERROR: Invalid user defined #COMMAND in user_read_inputs. '
             write(*,*) '--Check user_read_inputs for errors'
             write(*,*) '--Check to make sure a #USERINPUTEND command was used'
             write(*,*) '  *Unrecognized command was: '//NameCommand
             call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
         end if
       end select
    end do

!    call stop_user(Name)
  end subroutine user_read_inputs




  !========================================================================
  subroutine user_set_ICs
  !========================================================================

    use ModMain, ONLY: nBlock, UnusedBlk, ProcTest, body1_, globalBLK
    use ModProcMH, ONLY: iProc
    use ModPhysics
    use ModNumConst
    use ModVarIndexes, ONLY: Bx_, By_, Bz_
    use ModGeometry,ONLY: x_BLK,y_BLK,z_BLK,R_BLK
    use ModIO,ONLY:restart

    integer :: iBlock
    logical :: oktest, oktest_me

  real :: Theta, Phi, xR, xTheta,xPhi, unr_o, unr_i, unr0
  integer, parameter :: jTh_axis=4
  integer :: i,j,k
  integer:: iR, jTheta, kPhi, iRp1,jThetap1,kPhip1, jTh_axr, &
	kPhi_conj, kPhip1_conj


    !--------------------------------------------------------------------------

  if(restart)return!jet and Bx parameters has to be set in 1st session.
    if( DoInitialize ) then
      if(iProc==PROCtest)then
        write(*,*)'Initializing Comet Jet Data'
        write(*,*)'Parameters:'
        write(*,*)'jet_width=',jet_width,'Qprod_day =',Qprod_day
        write(*,*)'Qprod_nit=',Qprod_nit,' Qprod_jet =',Qprod_jet
        write(*,*)'Qprod_jeta =',Qprod_jeta, ' SW_By=',SW_By

        call set_oktest('user_initial_perturbation',oktest,oktest_me)
      else
        oktest=.false.; oktest_me=.false.
      end if

!       jet_width = jet_width*cPi/180.
!       Qprod_day = Qprod_day*2.0
!       Qprod_nit = Qprod_nit*2.0
!       if ( jpattern == 1 ) Qprod_jet = Qprod_jet*4.0*cPi
!       if ( jpattern == 2 ) Qprod_jeta = Qprod_jeta*4.0*cPi
!       if ( Ajet < cTiny ) then
!         Qprod_jet = cZero
!       else
!         Qprod_jet = Qprod_jet/Ajet
!       endif
!       if ( Ajet1 < cTiny ) then
!         Qprod_jeta = cZero
!       else
!         Qprod_jeta = Qprod_jeta/Ajet1
!       endif


!this is set after set_IC_comet, but this part is not used there.
      FaceState_VI(Bx_,  body1_) = SW_Bx
      FaceState_VI(By_,  body1_) = SW_By
      FaceState_VI(Bz_,  body1_) = SW_Bz

      DoInitialize = .False.
    endif

    if(ReadNeutral) then
  jTh_axr = NTheta_neutral-jTh_axis

  do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1

     if( R_BLK(i,j,k,globalBLK) < 7.0E8/NO2SI_V(UnitX_) .and.  &
		R_BLK(i,j,k,globalBLK) > 2000.0/NO2SI_V(UnitX_) )then     !interpolate

        iR = int( (NR_neutral-1)/6.*log10(R_BLK(i,j,k,globalBLK)*NO2SI_V(UnitX_)/1E3) ) + 1
        iR = max(iR,1)
        iR = min(iR, NR_neutral-1)
	iRp1 = iR+1
        Theta = acos( x_BLK(i,j,k,globalBLK)/R_BLK(i,j,k,globalBLK) )/cPi
        jTheta = int( Theta*(NTheta_neutral-1) ) + 1
	jThetap1 = jTheta+1
        Phi = acos( y_BLK(i,j,k,globalBLK)/sqrt( &
		y_BLK(i,j,k,globalBLK)*y_BLK(i,j,k,globalBLK)+ &
		z_BLK(i,j,k,globalBLK)*z_BLK(i,j,k,globalBLK) )+cTiny8/1E3 )
	Phi = Phi*sign( cOne,z_BLK(i,j,k,globalBLK) ) + &
		cPi*( cOne-sign(cOne,z_BLK(i,j,k,globalBLK)) )
        Phi = Phi / cTwoPi
        kPhi = int( Phi*(NPhi_neutral-1) ) + 1

        kPhi = max(kPhi,1)
	kPhip1 = mod( kPhi, NPhi_neutral ) + 1
        kPhi_conj   = mod( kPhi  +(NPhi_neutral-1)/2,NPhi_neutral )+1
        kPhip1_conj = mod( kPhip1+(NPhi_neutral-1)/2,NPhi_neutral )+1

        xR = ( R_BLK(i,j,k,globalBLK)*NO2SI_V(UnitX_)-NeutralN(iR,jTheta,kPhi,NR_) )/ &
		( NeutralN(iRp1,jTheta,kPhi,NR_)-NeutralN(iR,jTheta,kPhi,NR_) )
        xPhi = ( Phi*360.-NeutralN(iR,jTheta,kPhi,NPhi_) )/ &
		( NeutralN(iR,jTheta,kPhip1,NPhi_)-NeutralN(iR,jTheta,kPhi,NPhi_) )

	if(jTheta>=jTh_axis .and. jtheta<=jTh_axr) then
!	if( .true. ) then
           xTheta = ( Theta*180.-NeutralN(iR,jTheta,kPhi,NTh_) )/ &
                ( NeutralN(iR,jThetap1,kPhi,NTh_)-NeutralN(iR,jTheta,kPhi,NTh_) )
           Neutral_BLK(i,j,k,globalBLK,:)= &
	     ( (NeutralN(iR  ,jTheta  ,kPhi  ,7:10)*(cOne-xPhi)+ 	        &
		NeutralN(iR  ,jTheta  ,kPhip1,7:10)*      xPhi)*(cOne-xTheta)+  &
               (NeutralN(iR  ,jThetap1,kPhi  ,7:10)*(cOne-xPhi)+ 		&
		NeutralN(iR  ,jThetap1,kPhip1,7:10)* xPhi)*xTheta )*(cOne-xR) + &
             ( (NeutralN(iRp1,jTheta  ,kPhi  ,7:10)*(cOne-xPhi)+ 		&
		NeutralN(iRp1,jTheta  ,kPhip1,7:10)*      xPhi)*(cOne-xTheta)+  &
               (NeutralN(iRp1,jThetap1,kPhi  ,7:10)*(cOne-xPhi)+ 		&
		NeutralN(iRp1,jThetap1,kPhip1,7:10)*      xPhi)*xTheta )*xR
	elseif( jtheta<jTh_axis ) then
           unr_i = sqrt( &
		NeutralN(iR  ,jTh_axis,1,7)*NeutralN(iR  ,jTh_axis,1,7) + &
                NeutralN(iR  ,jTh_axis,1,8)*NeutralN(iR  ,jTh_axis,1,8) + &
                NeutralN(iR  ,jTh_axis,1,9)*NeutralN(iR  ,jTh_axis,1,9) )
           unr_o = sqrt( &
		NeutralN(iRp1,jTh_axis,1,7)*NeutralN(iRp1,jTh_axis,1,7) + &
                NeutralN(iRp1,jTh_axis,1,8)*NeutralN(iRp1,jTh_axis,1,8) + &
                NeutralN(iRp1,jTh_axis,1,9)*NeutralN(iRP1,jTh_axis,1,9) )
	   unr0 = unr_i*(cOne-xR)+unr_o*xR
           Neutral_BLK(i,j,k,globalBLK,1) = unr0*cos(Theta*cPi)
           Neutral_BLK(i,j,k,globalBLK,2) = unr0*sin(Theta*cPi)*cos(Phi*cTwoPi)
           Neutral_BLK(i,j,k,globalBLK,3) = unr0*sin(Theta*cPi)*sin(Phi*cTwoPi)
           xTheta =  (NeutralN(iR,jTh_axis,kPhi,NTh_)-Theta*180.)* &
		cHalf/NeutralN(iR,jTh_axis,kPhi,NTh_)
           Neutral_BLK(i,j,k,globalBLK,4) =   &
             ( (NeutralN(iR  ,jTh_axis,kPhi       ,10)*(cOne-xPhi)+ 		   &
                NeutralN(iR  ,jTh_axis,kPhip1     ,10)*   xPhi)*(cOne-xTheta)+     &
               (NeutralN(iR  ,jTh_axis,kPhi_conj  ,10)*(cOne-xPhi)+ 		   &
		NeutralN(iR  ,jTh_axis,kPhip1_conj,10)* xPhi)*xTheta )*(cOne-xR) + &
             ( (NeutralN(iRp1,jTh_axis,kPhi       ,10)*(cOne-xPhi)+ 		   &
                NeutralN(iRp1,jTh_axis,kPhip1     ,10)*      xPhi)*(cOne-xTheta)+  &
               (NeutralN(iRp1,jTh_axis,kPhi_conj  ,10)*(cOne-xPhi)+ 		   &
		NeutralN(iRp1,jTh_axis,kPhip1_conj,10)*      xPhi)*xTheta )*xR
	   Neutral_BLK(i,j,k,globalBLK,4) = Neutral_BLK(i,j,k,globalBLK,4) * &
		( cOne+1.5*(NeutralN(iR,jTh_axis,kPhi,NTh_)/180.-Theta) )
	else
           unr_i = sqrt( &
		NeutralN(iR  ,jTh_axr,1,7)*NeutralN(iR  ,jTh_axr,1,7) + &
                NeutralN(iR  ,jTh_axr,1,8)*NeutralN(iR  ,jTh_axr,1,8) + &
                NeutralN(iR  ,jTh_axr,1,9)*NeutralN(iR  ,jTh_axr,1,9) )
           unr_o = sqrt( &
		NeutralN(iRp1,jTh_axr,1,7)*NeutralN(iRp1,jTh_axr,1,7) + &
                NeutralN(iRp1,jTh_axr,1,8)*NeutralN(iRp1,jTh_axr,1,8) + &
                NeutralN(iRp1,jTh_axr,1,9)*NeutralN(iRp1,jTh_axr,1,9) )
           unr0 = unr_i*(cOne-xR)+unr_o*xR
           Neutral_BLK(i,j,k,globalBLK,1) = unr0*cos(Theta*cPi)
           Neutral_BLK(i,j,k,globalBLK,2) = unr0*sin(Theta*cPi)*cos(Phi*cTwoPi)
           Neutral_BLK(i,j,k,globalBLK,3) = unr0*sin(Theta*cPi)*sin(Phi*cTwoPi)
           xTheta = (Theta*180.-NeutralN(iR,jTh_axr,kPhi,NTh_))*cHalf/ &
		NeutralN(iR,jTh_axis+1,kPhi,NTh_)
           Neutral_BLK(i,j,k,globalBLK,4) =  &
             ( (NeutralN(iR  ,jTh_axr,kPhi       ,10)*(cOne-xPhi)+ 		  &
                NeutralN(iR  ,jTh_axr,kPhip1     ,10)*      xPhi)*(cOne-xTheta)+  &
               (NeutralN(iR  ,jTh_axr,kPhi_conj  ,10)*(cOne-xPhi)+ 		  &
		NeutralN(iR  ,jTh_axr,kPhip1_conj,10)* xPhi)*xTheta )*(cOne-xR) + &
             ( (NeutralN(iRp1,jTh_axr,kPhi	 ,10)*(cOne-xPhi)+		  &
                NeutralN(iRp1,jTh_axr,kPhip1     ,10)*      xPhi)*(cOne-xTheta)+  &
               (NeutralN(iRp1,jTh_axr,kPhi_conj  ,10)*(cOne-xPhi)+ 		  &
		NeutralN(iRp1,jTh_axr,kPhip1_conj,10)*      xPhi)*      xTheta )*xR
          Neutral_BLK(i,j,k,globalBLK,4) = Neutral_BLK(i,j,k,globalBLK,4) * &
                ( cOne+1.2*(NeutralN(iR,jTh_axr,kPhi,NTh_)/180.-Theta) )
	endif

     else if( R_BLK(i,j,k,globalBLK) > 7.0E8/NO2SI_V(UnitX_)) then	!Extrapolate beyond 10^6km

        Theta = acos( x_BLK(i,j,k,globalBLK)/R_BLK(i,j,k,globalBLK) )/cPi
        jTheta = int( Theta*(NTheta_neutral-1) ) + 1
        jTheta = max(jTheta,jTh_axis)          ! temperary boundary solution
        jTheta = min(jTheta, NTheta_neutral-jTh_axis)
        Phi = acos( y_BLK(i,j,k,globalBLK)/sqrt( &
		y_BLK(i,j,k,globalBLK)*y_BLK(i,j,k,globalBLK)+ &
		z_BLK(i,j,k,globalBLK)*z_BLK(i,j,k,globalBLK))+cTiny8/1E3 )
        Phi = Phi*sign( cOne,z_BLK(i,j,k,globalBLK) )
        Phi = ( Phi + cPi*(cOne-sign(cOne,z_BLK(i,j,k,globalBLK))) ) / cTwoPi
        kPhi = int( Phi*(NPhi_neutral-1) ) + 1
        kPhi = max(kPhi,1)

        Neutral_BLK(i,j,k,globalBLK,1:3)= NeutralN(NR_Neutral-2,jTheta,kPhi,7:9)
        Neutral_BLK(i,j,k,globalBLK,4)= .8*NeutralN(NR_Neutral-2,jTheta,kPhi,10)/ &
		R_BLK(i,j,k,globalBLK)/R_BLK(i,j,k,globalBLK)* &
		exp( -R_BLK(i,j,k,globalBLK)*NO2SI_V(UnitX_)*ionization_rate/ &
		sqrt(Neutral_BLK(i,j,k,globalBLK,1)*Neutral_BLK(i,j,k,globalBLK,1) + &
		     Neutral_BLK(i,j,k,globalBLK,2)*Neutral_BLK(i,j,k,globalBLK,2) + &
		     Neutral_BLK(i,j,k,globalBLK,3)*Neutral_BLK(i,j,k,globalBLK,3)) )

     else	!Extrapolate within 3km
if(.false.) then
        Theta = acos( x_BLK(i,j,k,globalBLK)/R_BLK(i,j,k,globalBLK) )/cPi
        jTheta = int( Theta*(NTheta_neutral-1) ) + 1
        jThetap1 = jTheta+1
        Phi = acos( y_BLK(i,j,k,globalBLK)/sqrt( &
		y_BLK(i,j,k,globalBLK)*y_BLK(i,j,k,globalBLK)+ &
		z_BLK(i,j,k,globalBLK)*z_BLK(i,j,k,globalBLK))+cTiny8/1E3 )
        Phi = Phi*sign( cOne,z_BLK(i,j,k,globalBLK) ) + &
		cPi*( cOne-sign(cOne,z_BLK(i,j,k,globalBLK)) )
        Phi = Phi / cTwoPi
        kPhi = int( Phi*(NPhi_neutral-1) ) + 1
        kPhi = max(kPhi,1)
        kPhi = min(kPhi, NPhi_neutral-1)
        kPhip1 = mod( kPhi, NPhi_neutral ) + 1

        xPhi = ( Phi*360.-NeutralN(2,jTheta,kPhi,NPhi_) )/ &
                ( NeutralN(2,jTheta,kPhip1,NPhi_)-NeutralN(2,jTheta,kPhi,NPhi_) )

        if(jTheta>=jTh_axis .and. jtheta<=(NTheta_neutral-jTh_axis)) then
           xTheta = ( Theta*180.-NeutralN(2,jTheta,kPhi,NTh_) )/ &
                ( NeutralN(2,jThetap1,kPhi,NTh_)-NeutralN(2,jTheta,kPhi,NTh_) )
           Neutral_BLK(i,j,k,globalBLK,:) = &
             ( (NeutralN(2,jTheta  ,kPhi  ,7:10)*(cOne-xPhi)+ 		&
                NeutralN(2,jTheta  ,kPhip1,7:10)*  xPhi)*(cOne-xTheta)+ &
               (NeutralN(2,jThetap1,kPhi  ,7:10)*(cOne-xPhi)+ 		&
                NeutralN(2,jThetap1,kPhip1,7:10)*      xPhi)*  xTheta )
        elseif( jtheta<jTh_axis ) then
           xTheta = (NeutralN(2,jTh_axis,kPhi,NTh_)-Theta*180.)*cHalf/ &
		     NeutralN(2,jTh_axis,kPhi,NTh_)
           Neutral_BLK(i,j,k,globalBLK,:) = &
             ( (NeutralN(2,jTh_axis,kPhi       ,7:10)*(cOne-xPhi)+ 	    &
                NeutralN(2,jTh_axis,kPhip1     ,7:10)* xPhi)*(cOne-xTheta)+ &
               (NeutralN(2,jTh_axis,kPhi_conj  ,7:10)*(cOne-xPhi)+ 	    &
		NeutralN(2,jTh_axis,kPhip1_conj,7:10)*      xPhi)* xTheta )
        else
           xTheta = (Theta*180.-NeutralN(2,jTh_axr,kPhi,NTh_))* &
		cHalf/NeutralN(2,jTh_axis+1,kPhi,NTh_)
           Neutral_BLK(i,j,k,globalBLK,:) = &
             ( (NeutralN(2,jTh_axr,kPhi       ,7:10)*(cOne-xPhi)+ 	    &
                NeutralN(2,jTh_axr,kPhip1     ,7:10)*  xPhi)*(cOne-xTheta)+ &
               (NeutralN(2,jTh_axr,kPhi_conj  ,7:10)*(cOne-xPhi)+ 	    &
		NeutralN(2,jTh_axr,kPhip1_conj,7:10)*      xPhi)*  xTheta )
        endif

endif

!write(*,*) '3rd choice', R_BLK(i,j,k,globalBLK)
        Neutral_BLK(i,j,k,globalBLK,1:3)= cOne
        Neutral_BLK(i,j,k,globalBLK,4) = 40.*1E9*mbar
     end if

     Neutral_BLK(i,j,k,globalBLK,4) = max( Neutral_BLK(i,j,k,globalBLK,4), 50000.0 )

  end do;  end do;  end do

  Neutral_BLK(:,:,:,globalBLK,4) = Neutral_BLK(:,:,:,globalBLK,4)*17
!  Neutral_BLK(:,:,:,globalBLK,4) = Neutral_BLK(:,:,:,globalBLK,4)/22.


!!!now Neutral_BLK(v) is dimensionless, while n is still dimensional!!!
    Neutral_BLK(:,:,:,globalBLK,1:3) = &
	Neutral_BLK(:,:,:,globalBLK,1:3)/NO2SI_V(UnitU_)

!write(*,*) 'set_ICs ', Neutral_BLK(2,2,2,globalBLK, :)

     endif
  end subroutine user_set_ICs

  subroutine user_init_point_implicit

    use ModVarIndexes
    use ModPointImplicit, ONLY: iVarPointImpl_I, IsPointImplMatrixSet

    ! Allocate and set iVarPointImpl_I
    allocate(iVarPointImpl_I(4+nSpecies))

    iVarPointImpl_I = (/Rho_, RhoUx_, RhoUy_, RhoUz_, p_/)

    ! Note that energy is not an independent variable for the 
    ! point implicit scheme. The pressure is an independent variable,
    ! and in this example there is no implicit pressure source term.

    ! Tell the point implicit scheme if dS/dU will be set analytically
    ! If this is set to true the DsDu_VVC matrix has to be set below.
    ! Initialization for comet implicit sources: false => numerical ptimplicit.
    IsPointImplMatrixSet = .false.

  end subroutine user_init_point_implicit


  !========================================================================
  !  SUBROUTINE user_calc_sources
  !========================================================================
  subroutine user_calc_sources

    ! Evaluate the explicit or implicit or both source terms.
    ! If there is no explicit source term, the subroutine user_expl_source 
    ! and the corresponding calls can be removed.

    use ModPointImplicit, ONLY: UsePointImplicit, IsPointImplSource

    if(.not.UsePointImplicit)then
       ! Add all source terms if we do not use the point implicit scheme
       call user_expl_source
       call user_impl_source
    elseif(IsPointImplSource)then
       ! Add implicit sources only
       call user_impl_source
    else
       ! Add explicit sources only
       call user_expl_source
    end if

  end subroutine user_calc_sources
  !==========================================================================
  subroutine user_expl_source

    ! Here come the explicit source terms

  end subroutine user_expl_source
  !==========================================================================
  subroutine user_impl_source
  
    ! This is a test and example for using point implicit source terms
    ! Apply friction relative to some medium at rest
    ! The friction force is proportional to the velocity and the density.

    use ModPointImplicit, ONLY: &
         UsePointImplicit, iVarPointImpl_I, IsPointImplMatrixSet, DsDu_VVC

    use ModMain,    ONLY: GlobalBlk,nI,nJ,nK,n_step,iTest,jTest,kTest,BlkTest
    use ModAdvance, ONLY: State_VGB, Source_VC, &
         Rho_, RhoUx_, RhoUy_, RhoUz_, Bx_,By_,Bz_, p_, Energy_
    use ModGeometry,ONLY: x_BLK,y_BLK,z_BLK,R_BLK

    use ModPhysics
    use ModProcMH

    integer :: iBlock, i, j, k
    real    :: Coef

    integer, save :: step = 0
    logical, save :: FirstCall = .TRUE.


!****************    yingdong defined for comet    **************
!    real, parameter ::  x_jet = 0.855363194, y_jet = -0.49384417, &
!       z_jet = 0.156434465 	!borrelly jet pattern
    real, parameter ::  x_jet = 1., y_jet = 0., z_jet = 0.
	!encke sunward jet pattern
    real :: lambda, sMassdn, sMassj, &
		jtheta, rcyl
    real, dimension(1:nI,1:nJ,1:nK) :: &
       usqr,sMass,Te,alphaTe,term1,term2,ne, eta, Rkm, logR, fi &
	,sMasseta, Losse, chargexchg, Unx, Uny, Unz, ux, uy, uz

    !-------------------------------------------------------------------------

    iBlock = GlobalBlk

    if (FirstCall) then
       FirstCall = .False.
       step = n_step

       jet_width = jet_width*cPi/180.
       Qprod_day = Qprod_day*2.0
       Qprod_nit = Qprod_nit*2.0
       if ( jpattern == 1 ) Qprod_jet = Qprod_jet*4.0*cPi
       if ( jpattern == 2 ) Qprod_jeta = Qprod_jeta*4.0*cPi
       if ( Ajet < cTiny ) then
          Qprod_jet = cZero
       else
          Qprod_jet = Qprod_jet/Ajet
       endif
       if ( Ajet1 < cTiny ) then
          Qprod_jeta = cZero
       else
          Qprod_jeta = Qprod_jeta/Ajet1
       endif
       
!        if (iProc == 0) then
!          open(unit=321,file='Source_VCs.log',status='unknown', action ='write', position='rewind')
!          write(321,*) "iter rho_ rhoUx_ rhoUy_ rhoUz_ p_ Energy_"
!          close(321)
!        end if

    end if

!*********************************************************
    ! Add implicit source here
    ! In this example a simple friction term is added to the momentum and 
    ! energy equtaions.

    ux=State_VGB(rhoUx_,1:nI,1:nJ,1:nK,globalBLK) / &
          State_VGB(rho_,1:nI,1:nJ,1:nK,globalBLK)
    uy=State_VGB(rhoUy_,1:nI,1:nJ,1:nK,globalBLK) / &
          State_VGB(rho_,1:nI,1:nJ,1:nK,globalBLK)
    uz=State_VGB(rhoUz_,1:nI,1:nJ,1:nK,globalBLK) / &
          State_VGB(rho_,1:nI,1:nJ,1:nK,globalBLK)
    usqr = ux*ux+uy*uy+uz*uz
    Rkm = R_BLK(1:nI,1:nJ,1:nK,globalBLK)*NO2SI_V(UnitX_)/1E3	!km unit
    logR = log10(Rkm)

    if (jpattern == 50 ) then
       sMass = cZero
       lambda= 1.80000*1E6
       where(R_BLK(1:nI,1:nJ,1:nK,globalBLK)<6.) &
	  sMass = Qprod*(R_BLK(1:nI,1:nJ,1:nK,globalBLK)*NO2SI_V(UnitX_)/lambda)**(-3.5)

    else
       lambda = Unr/ionization_rate
      !! fi multiplicator value for the ionization frequency (including enhanced electron impact in ion pile up reagion)
      !       fi = cOne	!       set f_i=1 for all r
       do k=1,nK ;   do j=1,nJ ;   do i=1,nI       
	  if (rkm(i,j,k) >= 5000. .and. rkm(i,j,k) < 10000.) then
             fi(i,j,k) = 1.0+0.77*log(rkm(i,j,k)/5000.)
	  elseif (rkm(i,j,k) >= 10000. .and. rkm(i,j,k) < 50000.) then
             fi(i,j,k) = 1.5-0.31067*log(rkm(i,j,k)/10000.)
	  else
	     fi(i,j,k) = cOne
	  endif
       end do ;  end do ; end do

       sMass = Qprod * mbar * fi * &
          exp(-R_BLK(1:nI,1:nJ,1:nK,globalBLK)*NO2SI_V(UnitX_)/lambda) / &
          (4.0*cPi*lambda*R_BLK(1:nI,1:nJ,1:nK,globalBLK)**2*NO2SI_V(UnitX_)**2)

    !yingdong 060705 neutral num density
       if(ReadNeutral) sMass = fi * ionization_rate * &
	  Neutral_BLK(1:nI,1:nJ,1:nK,globalBLK,4)

     ! sMass(amu/s/m^3)/NO2SI_V(UnitN_)*NO2SI_V(UnitT_) 
       sMass = sMass/NO2SI_V(UnitN_)*NO2SI_V(UnitT_)

     ! eta is already non-dimensional
       eta = kin/mbar/fi/ionization_rate*NO2SI_V(UnitN_)
! added yingdong Aug22,03 for increased hydrogen CX.
!       chargexchg = eta* ( 1.+exp(11.*R_BLK(1:nI,1:nJ,1:nK,globalBLK)* &
!                       NO2SI_V(UnitX_)/12./lambda)/6. )
        chargexchg = eta        !original term
!********   modification Apr 03 yingdong for Barrelley /end ***********
!***** 3 etas r changed into chargexchg *****
!**************************************

    endif

     ! ne is the dimensionless electron density
    ne = State_VGB(rho_,1:nI,1:nJ,1:nK,globalBLK)*NO2SI_V(UnitN_)/mbar

     ! Calculate electron temperature from pressure
     ! assuming Te=Ti
     ! for comet borrelly
!    Te=State_VGB(p_,1:nI,1:nJ,1:nK,globalBLK) * NO2SI_V(UnitP_) * mbar * cProtonMass / &
!          ( 2.0 * unitSI_rho * cBoltzmann * State_VGB(rho_,1:nI,1:nJ,1:nK,globalBLK) )

    ! Standard Profile
    where  (rkm <= 1584.893) 
       Te = 1.E+2
    end where
    where (rkm > 1584.893 .and. rkm <= 6918.310) 
       Te = 10.**( 1.143  * logR -  1.667 )
    end where
    where (rkm > 6918.310 .and. rkm <= 1.E+4)
       Te = 10.**(10.965  * logR - 39.3735)
    end where
    where (rkm > 1.E+4 .and. rkm <= 1.E+5)
       Te = 10.**( 0.5135 * logR +  2.43)
    end where
    where (rkm > 1.E+5)
       Te = 1.E+5
    end where


!********   modification Apr 03 yingdong for Borelley //start  *********
    if ( jet_width < cTiny ) then
	jet_ln2 = cZero
    else
	jet_ln2 = -dlog(2.0)/jet_width/jet_width
    endif

!     do k=1,nK ;   do j=1,nJ ;    do i=1,nI	!version with jet
    do k=1,0 ;   do j=1,0 ;    do i=1,0	!version without jet

     	!** the jet is located at 30 degree towards the +y in the xy plane   ***
     	!** and 9 degree towards the +z in the xz plane                      ***
     	!** the vector components of this vector are hard coded here         ***

	jtheta = ( X_BLK(i,j,k,globalBLK)*x_jet+Y_BLK(i,j,k,globalBLK)*y_jet+ &
			Z_BLK(i,j,k,globalBLK)*z_jet ) / R_BLK(i,j,k,globalBLK)
	if (jpattern /= 6 ) jtheta = dmax1(cZero, jtheta)
	jtheta = dmin1(cOne,  jtheta)		!jtheta = cos(theta) now

	if ( X_BLK(i,j,k,globalBLK) >= cZero ) then
	  sMassdn = sMass(i,j,k)*Qprod_day
	else
	  sMassdn = sMass(i,j,k)*Qprod_nit
	endif

	if( jpattern == 1 ) then			!exp jet
              jtheta = dacos(jtheta)
              sMass(i,j,k) = sMassdn+sMass(i,j,k)*Qprod_jet*exp(jet_ln2*jtheta*jtheta)
	elseif ( jpattern == 0 .or. jpattern == 7 ) then		!dayside cos jet
              sMass(i,j,k) = sMassdn+sMass(i,j,k)*Qprod_jet*jtheta
	elseif ( jpattern == 6 ) then		!all cos jet
              sMass(i,j,k) = sMassdn+sMass(i,j,k)*Qprod_jet*(cOne+jtheta)
	elseif ( jpattern == 4 ) then		!dayside liner jet
              jtheta = dacos(jtheta)
	    sMass(i,j,k) = sMassdn+sMass(i,j,k)*Qprod_jet*(cOne-2.0*jtheta/cPi)
	elseif ( jpattern == 5 ) then		!dayside cos2 jet
              sMass(i,j,k) = sMassdn+sMass(i,j,k)*Qprod_jet*jtheta*jtheta
	elseif ( jpattern == 2 ) then
              sMassj = sMassdn+sMass(i,j,k)*Qprod_jet*jtheta
              jtheta = dacos(jtheta)
              sMass(i,j,k) = sMassj+sMass(i,j,k)*Qprod_jeta* &
		exp(jet_ln2*jtheta*jtheta)*.5*( 1.-tanh((.15-R_BLK(i,j,k,globalBLK))*64.) )
	endif

    end do ;  end do ; end do 

    if (jpattern == 50 ) then
	Losse = cZero
	sMasseta = sMass*ionization_rate*state_VGB(rho_,1:nI,1:nJ,1:nK,globalBLK)
!	sMasseta = sMass*ionization_rate
!        sMasseta = 0.
    else
     ! Define alpha.
       where  (Te < 200.) 
	        alphaTe = 7.E-7*sqrt(300./Te)
       elsewhere
	        alphaTe = 2.342*7.E-7*Te**(0.2553-0.1633*log10(Te))
       end where
     !normalize alphaTe
     !alpha Te has units [cm^3/s]
       alphaTe=alphaTe/1E6

     ! Compute source terms.
       Losse    = alphaTe*ne*NO2SI_V(UnitT_)	!added yingdong Jan 02 and Modified Mar.03.
       sMasseta = sMass*chargexchg	!modified yingdong Oct. 04 to seperate ionisation/friction
  
    endif

     if( jpattern == 3 ) then		!shade by body
	do k=1,nK ;   do j=1,nJ ;    do i=1,nI
          rcyl=sqrt( Z_BLK(i,j,k,globalBLK)**2 + Y_BLK(i,j,k,globalBLK)**2 )
          if ( X_BLK(i,j,k,globalBLK) <= cZero .and. rcyl <= Rbody )  sMass(i,j,k) = 0.
	end do; end do; end do
     end if

     term1 = sMass+sMasseta*State_VGB(rho_,1:nI,1:nJ,1:nK,globalBLK)
     term2 = sMasseta+Losse

     if( jpattern == 7 ) then           !debug pattern
        do k=1,nK ;   do j=1,nJ ;    do i=1,nI
          if ( R_BLK(i,j,k,globalBLK) <= 3.e-4 .and. y_BLK(i,j,k,globalBLK) >= -1.e-5 .and. &
            R_BLK(i,j,k,globalBLK) >= 8.e-5 .and. abs(z_BLK(i,j,k,globalBLK)) <= 2.e-5 &
                .and. z_BLK(i,j,k,globalBLK) >= -1.e-6 .and. n_step == 4001 ) &
              write(*,*) 'xyzr, smass, Losse, Losse*rho, sMasseta, term1, term2,', &
                ' term2*rhou, p', &
                X_BLK(i,j,k,globalBLK), Y_BLK(i,j,k,globalBLK), Z_BLK(i,j,k,globalBLK), &
                R_BLK(i,j,k,globalBLK), sMass(i,j,k), Losse(i,j,k), &
                Losse(i,j,k)*State_VGB(rho_,i,j,k,globalBLK), sMasseta(i,j,k), term1(i,j,k), &
                term2(i,j,k), term2(i,j,k)*State_VGB(rhoUx_,i,j,k,globalBLK), usqr(i,j,k), &
                State_VGB(p_,i,j,k,globalBLK)
        end do; end do; end do
     end if

     Source_VC(rho_,:,:,:) = Source_VC(rho_,:,:,:) + &
	(sMass - Losse*State_VGB(rho_,1:nI,1:nJ,1:nK,globalBLK))

     if(ReadNeutral) then       ! yingdong 060605 neutral profile
        Unx = Neutral_BLK(1:nI,1:nJ,1:nK,globalBLK,1)
        Uny = Neutral_BLK(1:nI,1:nJ,1:nK,globalBLK,2)
        Unz = Neutral_BLK(1:nI,1:nJ,1:nK,globalBLK,3)

!write(*,*) 'user_src ',  Neutral_BLK(2,2,2,globalBLK,:)

        do k=1,nK ;   do j=1,nJ ;    do i=1,nI
          unr = Unx(i,j,k)*Unx(i,j,k) + Uny(i,j,k)*Uny(i,j,k) +  &
		Unz(i,j,k)*Unz(i,j,k)
!note: this unr is unr*unr*dimensionless

          Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) + &
                ( term1(i,j,k)* (0.5*unr) -  &
                term2(i,j,k)*(0.5*State_VGB(rho_,i,j,k,globalBLK)* &
                usqr(i,j,k) + 1.5*State_VGB(p_,i,j,k,globalBLK)) )	!Pi->P
!           usqr(i,j,k) + cHalf*1.5*State_VGB(p_,i,j,k,globalBLK)) )
        end do; enddo; enddo

     else
        Unx=Unr/NO2SI_V(UnitU_)*X_BLK(1:nI,1:nJ,1:nK,globalBLK)/ &
          R_BLK(1:nI,1:nJ,1:nK,globalBLK)
        Uny=Unr/NO2SI_V(UnitU_)*Y_BLK(1:nI,1:nJ,1:nK,globalBLK)/ &
          R_BLK(1:nI,1:nJ,1:nK,globalBLK)
        Unz=Unr/NO2SI_V(UnitU_)*Z_BLK(1:nI,1:nJ,1:nK,globalBLK)/ &
          R_BLK(1:nI,1:nJ,1:nK,globalBLK)

     	Source_VC(Energy_,:,:,:) = Source_VC(Energy_,:,:,:) + &
	  ( term1*(0.5*(Unr/NO2SI_V(UnitU_))**2) - &
          term2*(0.5*State_VGB(rho_,1:nI,1:nJ,1:nK,globalBLK)* &
          usqr + 1.5*State_VGB(p_,1:nI,1:nJ,1:nK,globalBLK)) )	!Pi->P
!          usqr + cHalf*1.5*State_VGB(p_,1:nI,1:nJ,1:nK,globalBLK)) )

    endif
    Source_VC(rhoUx_,:,:,:) = Source_VC(rhoUx_,:,:,:) + &
	( term1*Unx - term2*State_VGB(rhoUx_,1:nI,1:nJ,1:nK,globalBLK) )
    Source_VC(rhoUy_,:,:,:) = Source_VC(rhoUy_,:,:,:) + &
	( term1*Uny - term2*State_VGB(rhoUy_,1:nI,1:nJ,1:nK,globalBLK) )
    Source_VC(rhoUz_,:,:,:) = Source_VC(rhoUz_,:,:,:) + &
	( term1*Unz - term2*State_VGB(rhoUz_,1:nI,1:nJ,1:nK,globalBLK) )
    Source_VC(p_,:,:,:) = Source_VC(p_,:,:,:) + term1* &
        1.0/3.0*( (Unx-ux)**2+(Uny-uy)**2+(Unz-uz)**2 ) - &
        term2*State_VGB(p_,1:nI,1:nJ,1:nK,globalBLK)	!Pi->P
!          term2*State_VGB(p_,1:nI,1:nJ,1:nK,globalBLK)*cHalf

!     open(unit=321,file='Source_VCs.log',status='old', action ='write', position='append')
!     if(globalBLK==BLKtest) then
!        do k=1,nK ;  
!           do j=1,nJ ;    
!              do i=1,nI ;
!                 if(itest==i.and.jtest==j.and.ktest==k) then
!                    if (step == n_step) then
!                       write(321,123) n_step, Source_VC(rho_,i,j,k), Source_VC(rhoUx_,i,j,k),&
!                            Source_VC(rhoUy_,i,j,k),Source_VC(rhoUz_,i,j,k),&
!                            Source_VC(p_,i,j,k), Source_VC(Energy_,i,j,k)  
!                       123 format (i7,6(1x,E16.10))
!                       step = n_step + 1
!                    end if
!                 endif
!              enddo
!           enddo
!        enddo
!     endif
!    close(321)


    if(IsPointImplMatrixSet)then
       ! Set the non-zero dS/dU matrix elements here
!      term3    = (5.-3.*g)*(sMasseta+Losse)
!      term4    = 1.5*(sMasseta+Losse)	!for energy source
!      term4    = term2		!for pressure source

      DsDu_VVC = cZero

      DsDu_VVC(1,1,:,:,:) = - Losse
      DsDu_VVC(2,1,:,:,:) = sMasseta*Unx
      DsDu_VVC(2,2,:,:,:) = - term2
      DsDu_VVC(3,1,:,:,:) = sMasseta*Uny
      DsDu_VVC(3,3,:,:,:) = - term2
      DsDu_VVC(4,1,:,:,:) = sMasseta*Unz
      DsDu_VVC(4,4,:,:,:) = - term2
      DsDu_VVC(5,2,:,:,:) = -term1*(Unx-ux)*2.0*1.0/3.0
      DsDu_VVC(5,3,:,:,:) = -term1*(Uny-uy)*2.0*1.0/3.0
      DsDu_VVC(5,4,:,:,:) = -term1*(Unz-uz)*2.0*1.0/3.0
      DsDu_VVC(5,5,:,:,:) = - term2

      if(ReadNeutral) then       ! yingdong 060605 neutral profile
        do k=1,nK ;   do j=1,nJ ;    do i=1,nI
          unr = Unx(i,j,k)*Unx(i,j,k) + Uny(i,j,k)*Uny(i,j,k) +  &
		Unz(i,j,k)*Unz(i,j,k)
!note: this unr is unr*unr*dimensionless
          DsDu_VVC(5,1,i,j,k) = sMasseta(i,j,k)*1.0/3.0*(unr-usqr(i,j,k)) + &
             sMass(i,j,k)*2.0*1.0/3.0/State_VGB(rho_,i,j,k,globalBLK)* &
             ( Unx(i,j,k)*ux(i,j,k) + Uny(i,j,k)*uy(i,j,k) + &
	     Unz(i,j,k)*uz(i,j,k) - usqr(i,j,k) )
	enddo; enddo; enddo

      else
        DsDu_VVC(5,1,:,:,:) = sMasseta*1.0/3.0*(unr*unr/NO2SI_V(UnitU_)/NO2SI_V(UnitU_)-usqr) + &
             sMass*2.0*1.0/3.0/State_VGB(rho_,1:nI,1:nJ,1:nK,globalBLK)*( &
             Unx*ux+Uny*uy+Unz*uz-usqr)

      endif

    end if

  end subroutine user_impl_source

  !========================================================================
  !  SUBROUTINE user_update_states(iStage,iBlock)
  !========================================================================
  subroutine user_update_states(iStage,iBlock)
    use ModVarIndexes
    use ModSize
    use ModAdvance, ONLY: State_VGB
    use ModPhysics
    use ModEnergy
    integer,intent(in):: iStage,iBlock
    integer:: i,j,k

    call update_states_MHD(iStage,iBlock)

    !\
    ! Begin update check of temperature::
    ! now check to see if the temperature is less than some
    ! prescribed minimum. If it is set it to the minimum value
    !/

    where( State_VGB(p_,1:nI,1:nJ,1:nK,iBlock)*NO2SI_V(UnitP_) < &
           (State_VGB(rho_,1:nI,1:nJ,1:nK,iBlock)*NO2SI_V(UnitN_)/mbar)*cBoltzmann*Tion )
       State_VGB(p_,1:nI,1:nJ,1:nK,iBlock) = &
           (State_VGB(rho_,1:nI,1:nJ,1:nK,iBlock)*NO2SI_V(UnitN_)/mbar)*cBoltzmann*Tion/NO2SI_V(UnitP_)
    end where

    call calc_energy_cell(iBlock)

  end subroutine user_update_states


end module ModUser


