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
       cToleranceHere =  cTiny**2  ,    &  !012, 009 givs up to 4-5 poin
       cTwoPiInv = cOne/cTwo       ,    &
       Hmax  = cPlankH/cErg,            &  !   erg*s Planck
       HmaxT = cPlankHBar                  !   the same/2Pi


  real:: &  
       vNe        , &   ! cm^{-3}
       vTe        , &   ! eV
       NatomII    , &   ! a sum over concentrations of atoms+ions, cm^{-3}
       Z          , &   ! vNe/NatomII
       C_I(0:99)  , &   ! C_I[0]-neutrals resided after ionization
       U_I(99)     , &  ! U_I(iZ) - the energy to make an ion iZ+ from (iZ-1)+    
       SUi              ! sum(Uiz*Ciz)
  integer :: nCi = 54   ! for Xenon, must be a var for other elements
end module ModSaha
!=============================================================================

subroutine    ConcNafter  ! for given summ(No+Niz) on boundary II, calc. C_I[z]
  use ModSaha
  implicit NONE

  integer   :: iter 
  real      :: &
       vNTotal  ,   vNTotalInv,        &  
       x3, x1,x2,f1,f2, df 

  vNTotal   =     NatomII     
  vNTotalInv= 1.0/NatomII 
  C_I(0)    = vNTotal

  ! debuT:       
  write (*,*) 'ConcN: NatomII=', NatomII, ' Te=', vTe

  x1= C_I(0)*1.001;    f1=( Conc(x1) -vNTotal )*vNTotalInv;
  x2= C_I(0)*0.011;    f2=( Conc(x2) -vNTotal )*vNTotalInv;

  iter = 0

  neutral:   do while(  (f2*f2 > cToleranceHere).AND.( iter < 55 )   )
     df = (f2 -f1)
     if (0.0 == df) then
        df = cToleranceHere
     end if
              x3 = (-x2*f1 +  x1*f2 )/df !   (f2 - f1);
     if(      x3 < 1.00d-22 ) then
              x3 = 1.00d-22  ;
     else if( x3 > 1.00d+33 ) then  
              x3 = 1.00d+33  ;
     end if
     if(f1*f1 > f2*f2 ) then 
        x1=x2;  f1=f2; 
     end if
     x2 = x3;
     f2 = ( Conc(x3) - vNTotal )*vNTotalInv;
     

     if( 0.0 == f2  )  exit
     iter =iter +1

  end do  neutral
  return
  !..........
contains

  real function   Conc( N0after )
    real, intent(in) :: N0after  !  conc [1/cm3] of neautrals after ionization
    integer,parameter::nW=5
    integer:: &  
         k=1     ,  &
         iW      ,  &
         iter    ,  &
         iter2   ,  &
         iZer    ,  & !  start of slider   zb{|+0,|+1,|+2,|+3,|+4}
         zb
    real ::         &
         NORM    ,  & != N0after, 
         Te32    ,  & != vTe ,     Te32  = te*sqrt(te) 
         CTotal,                        &
         a4,                            &
         c1   ,  c2,                    & 
         bc   , d   ,dc   ,             &
         PopulationRatio_I(nW)      ,   &
         PopulationRatioNonid_I(nW) ,   &
         StatWeight_I(nW)           ,   & 
         s    , s1  , s2  , st  ,       &
         UW_I(1:5), Uis          ;


    s     = 0.0
    C_I(0)= N0after

    Te32 = vTe*SQRT(vTe)      
    NORM = N0after
    a4   =  6.050e21*Te32/NORM 


    !debuT:       
    write (*,*) 'Conc: Norm=', NORM

    iZer =1  ! for the start of SLIDing
    !....................................................

!    C_I(iZer+4)=10.  ; C_I(iZer4)=1.;  !only 11 Ui in use now, 5 is the wide of sliding window 
    SLIDER: do while ( izer <= (11-4)  ) 

       !debuT   
       write(*,*)'SLIDer==== iZer=',izer

       CTotal=1.0

       c2  =1.0

       C_I(1 : nCi)=0.0
 

       UW_I = U_I(izer : izer+nW-1)	
       PopulationRatio_I = a4*exp(- UW_I/vTe)
     
       c1 = 1.3*sqrt(c2*PopulationRatio_I(1));

       iter2= 1
       dc   = 1.
       zb   = iZer


       WINDOW: do while ( ( dc>cToleranceHere ).AND.( iter2 < 30 )  )

          PopulationRatioNonid_I=PopulationRatio_I !To be mofified for non-ideal plasma
          
          iter =0 
          d    =1.

          NEWTON : do while(  (iter <30  ).AND.( d*d > cToleranceHere)) 
        
             bc = 1. / c1			
             StatWeight_I(1) = PopulationRatioNonid_I(1)*bc
             do iW=2,nW
                StatWeight_I(iW)=StatWeight_I(iW-1)*&
                         PopulationRatioNonid_I(iW)*bc
             end do
             s=0.0; s1=0
             do iW=1,nW
                s  = s  + StatWeight_I(iW)*(zb+iW-1)       
                s1 = s1 + StatWeight_I(iW)*(zb+iW-1)**2  
             end do
             d  =  c1 -c2 * s 
             c1 =( c1 *c2 *(s +s1)/(c1 +c2*s1) )
             iter= iter +1
          end do NEWTON

          s2 = sum(StatWeight_I)         
          Z  = s /s2;

          dc  = (c1 -CTotal)/(c1+CTotal)
          dc  = dc*dc;
          CTotal = c1; 
          iter2  = iter2+1  
       end do WINDOW


       C_I(izer : izer+4 ) = StatWeight_I* NORM
       vNe    = c1* NORM ;

       CTotal =   sum(C_I( 0 : iZer+4 )) 
     

       if( C_I(izer+4) < C_I(izer+0) ) EXIT
       iZer=iZer+1

    end do SLIDER

    Conc = CTotal

    SUi  = 0.0  ![eV]

    Uis  = sum( U_I(1:iZer) )  ![eV]

    SUi = Uis *C_I(iZer)   

    do   k = iZer+1, iZer+nW-1
       Uis = Uis +U_I(k)
       SUi = SUi +Uis*C_I(k)
    end do

    SUi = SUi/vNe   ! [eV]*[cm^3]     

    return 
  end function Conc
  !........................
end subroutine ConcNafter


!=======================================================================
subroutine      get_ionization_equilibrium


  !****   All   input PARAMETERS 
  !****   MUST   BE transformed 
  !****   before  CONCnAFTER call

end subroutine  get_ionization_equilibrium

!=======================================================================
