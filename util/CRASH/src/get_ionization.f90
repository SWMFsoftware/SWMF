!^CFG COPYRIGHT UM
! In the initial CRASH treatment of materials,
!
! 1. We can treat any given material as having a single average Z.
!
! 2. Mixtures can be treated either as an average material or as mixtures.
!
! 3. If mixtures are treated as mixtures, then collisional rates should be 
! calculated using the "effective Z", which is the average of Z squared 
! divided by the average Z.
!
! 4. The average Z can be determined from the Saha equation, but must not 
! exceed the nuclear charge of the material in question.
!
! 5. We do not need to account for electron degeneracy in the initial model.
!
! 6. In our regime of interest, the electrons behave as an ideal gas in an 
! ion-sphere environment within which Coulomb interactions do affect the 
! electron pressure and internal energy. The electron pressure and internal 
! energy are best calculated using equations 3.47 through 3.50 in 
! R. P. Drake, High Energy Density Physics
!
! 7. The ion pressure is the ideal gas pressure. The ion internal energy 
! includes the particle energy of random motion and the energy of ionization. 
! The model of eqs. 3.74 through 3.76 in the mentioned book is acceptable. 
! Alternatively, a more complex model using actual ionization energies would 
! be acceptable.
!
! 8. The materials that matter are
! Beryllium
! Xenon
! Polyimide (C_22 H_10 N_2 O_5)
!
!***********************************************************************
!    calculation of ionization equilibrium
!    for given - concentration of  atomic particles Co+Ciz [1/cm^3]
!              - electron temperature Te[eV]
!***********************************************************************

module ModSaha
  use ModConst
  implicit NONE
  include 'CRASH_definitions.h'
  real,  parameter ::  &
       cToleranceHere =  cTiny**2 ,     &  !012, 009 givs up to 4-5 digits
       cTwoPiInv = cOne /cTwoPi   ,     &
       Hmax  = cPlankH  /cErg     ,     &  !   erg*s Planck
       HmaxT = cPlankHBar                  !   the same/2Pi


  real:: &  
       vNe        , &   ! cm^{-3}
       vTe        , &   ! eV
       vNatomII   , &   ! a sum over concentrations of atoms+ions, cm^{-3}
       Z          , &   ! vNe/sum(C_I) <= only ions  
       C_I(0:99)  , &   ! C_I[0]-neutrals resided after ionization
       U_I(1:99)  , &   ! U_I(iZ) - the energy to make an ion iZ+ from (iZ-1)+    
       CUITotal         ! sum(Uiz*Ciz)/Ce  <~> ionization energy per 1 electron
  integer :: nCi = 54   ! for Xenon, must be a var for other elements
end module ModSaha
!=============================================================================

subroutine    ConcNafter  ! for given summ(No+Niz) on boundary II, calc. C_I[z]
  use ModSaha
  implicit NONE

  integer   :: iter 
  real      :: &
       vNTotal  ,   vNTotalInv,   & !  
       x3, x1,x2,f1,f2, df            

  vNTotal   =     vNatomII     
  vNTotalInv= 1.0/vNatomII   ! Inverted value, often used       
  C_I(0)    = vNTotal

  ! debuT:    !write (*,*) 'ConcN: NatomII=', NatomII, ' Te=', vTe

  x1= C_I(0)*1.001;  f1=( Conc(x1) -vNTotal )*vNTotalInv; ! loop initiation 
  x2= C_I(0)*0.011;  f2=( Conc(x2) -vNTotal )*vNTotalInv; ! to find  

  iter = 0

  neutral:  do while(  (f2*f2 > cToleranceHere).AND.( iter < 55 )   )
     df = (f2 -f1)
     if (0.0 == df) then
        df = cToleranceHere    ! if (x1!=x2) must  be additional analysis
        EXIT                   ! befor  ABORT
     end if
              x3 = (-x2*f1 +  x1*f2 )/df !   (f2 - f1);
     if(      x3 < 1.00d-22 )  then
              x3 = 1.00d-22  ; ! bottom limit for trial concentration
     else if( x3 > 1.00d+33 )  then  
              x3 = 1.00d+33  ; ! upper limit, , to avoid precision lost or overflow  
     end if
     if(f1*f1 > f2*f2 ) then 
        x1=x2;  f1=f2; 
     end if
     x2 = x3;
     f2 = ( Conc(x3) - vNTotal )*vNTotalInv;
     
     if( 0.0 == f2  )  exit
     iter =iter +1

  ! debuT:  !write (*,*) 'ConcN: Iter=',iter

  end do  neutral
  return
  !..........
contains
 ! 
 !  To solve an equation: sum(C)=Na, where Na is a given consentration
 !  of heavy component (residual neatrals + all ion states) we
 !  take a trial value of the residual neutral, N0After, and calculate
 !  sum(C), with this trial value. At high temperatures, instead of 
 !  the residual nuetrals we take the concentration of ions in such
 !  charge state that the states with lower Z can be treated as non-populated.  
 !
  function   conc( vN0after) result(CTotal)
    real, intent(in) :: vN0after ! conc [cm^{-3}] of neautrals(*) after ionization
    integer,parameter:: nW=5    ! range of most populated Zi cncentrations above iZmin
    integer:: &        
         iW      ,  &  ! index within nW
         iZMin   ,  &  ! start of slider   iZmin{|+0,|+1,|+2,|+3,|+4...!nW-1}
         iter, iter2    
    real ::         &
         vNORM   ,  &  !  == vN0after, 
         Te32    ,  &  !  Te32  = Te*sqrt(Te) 
         CTotal,                        &
         a4,                            &
         NeIterated   ,     NeInv   ,   &
         DeltaNe   ,DeltaC   ,          &
         PopulationRatio_I(nW)      ,   &
         PopulationRatioNonid_I(nW) ,   &
         StatWeight_I(nW)           ,   & 
         CZTotal    , CZ2Total ,        &
         UW_I(1:5), UITotal          

    real,parameter :: cElectronStatWeight=6.050e21 !To be clarified,<~> CGS units!


    CZTotal  = 0.0
    C_I(0) = vN0after

    Te32   = vTe*SQRT(vTe)      
    vNORM  = vN0after
    a4     = Te32/vNORM *  cElectronStatWeight


    !debuT:       
    write (*,*) 'Conc: Norm=', vNORM

    iZMin =1  ! for the start of SLIDing
    !....................................................

    SLIDER: do while ( iZmin <= (11-(nW-1))  ) 

       !debuT   
       write(*,*)'SLIDer==== iZMin=',iZmin

       C_I(1 : nCi)=0.0   ! Initialization
       CTotal=1.0         ! To start do loop

       !Get nw values of the ionization potentials from the table	
       UW_I = U_I(iZmin : iZmin+nW-1)  

       !Population ratio, according to the Boltzmann principle
       PopulationRatio_I = a4*exp(- UW_I/vTe)


       !Trial value for the electron concentration
       NeIterated = 1.3*sqrt( PopulationRatio_I(1));

       iter2   = 1
       DeltaC  = 1.      ! To start do loop

       WINDOW: do while ( ( DeltaC>cToleranceHere ).AND.( iter2 < 30 )  )

          PopulationRatioNonid_I=PopulationRatio_I !To be mofified for non-ideal plasma

          iter =0 
          DeltaNe    =1.   ! To start do loop

          NEWNe : do while(  (iter <30  ).AND.( DeltaNe*DeltaNe > cToleranceHere)) 

             NeInv = 1. / NeIterated			
             StatWeight_I(1) = PopulationRatioNonid_I(1)*NeInv
             do iW=2,nW
                StatWeight_I(iW)=StatWeight_I(iW-1)*&
                     PopulationRatioNonid_I(iW)*NeInv
             end do

             CZTotal=0.0; CZ2Total=0.0
             do iW=1,nW
                CZTotal  = CZTotal  + StatWeight_I(iW)*(iZMin+iW-1)       
                CZ2Total = CZ2Total + StatWeight_I(iW)*(iZMin+iW-1)**2  
             end do

             DeltaNe    =  NeIterated - CZTotal 
             NeIterated =( NeIterated *(CZTotal    +CZ2Total)/&
                  (NeIterated + CZ2Total) )
             iter= iter +1

          end do NEWNe

          Z  = CZTotal /sum(StatWeight_I)  

          DeltaC = (NeIterated - CTotal)/(NeIterated + CTotal)
          DeltaC = DeltaC*DeltaC;
          CTotal = NeIterated; 
          iter2  = iter2+1  
       end do WINDOW


       C_I(iZmin : iZmin+4 ) = StatWeight_I* vNORM ![cm^{-3}]
       vNe    = NeIterated * vNORM ;
       CTotal =   sum(C_I( 0 : iZMin+nW-1 )) 

       if( C_I(iZMin+nW-1) < C_I(iZMin+0) ) then
       !debuT:             write (*,*) 'Conc: Norm=', vNORM
              EXIT
       end if 

       iZMin=iZMin+1

    end do SLIDER

    !Finalize: get sum_i(C_i*sum(U_I(1:i)))

    UITotal  = sum( U_I(1:iZMin) )  ![eV]

    CUITotal = UITotal *C_I(iZMin)  ![eV]*[cm^{-3}]

    do   iw = iZMin+1, iZMin+nW-1
       UITotal  = UITotal  +U_I(iw)
       CUITotal = CUITotal +UITotal*C_I(iw)
    end do

    CUITotal = CUITotal / vNe   ! [eV] per electron    
  end function conc
  !........................
end subroutine ConcNafter


!=======================================================================
subroutine      get_ionization_equilibrium

  !****   All   input PARAMETERS 
  !****   MUST   BE transformed 
  !****   before  CONCnAFTER call

end subroutine  get_ionization_equilibrium

!=======================================================================
