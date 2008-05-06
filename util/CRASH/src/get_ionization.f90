
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

!//================================================
!//    calculation of ionization equilibrium
!//    for given  concentration of  atomic particles Ca+Ciz [1/cm^3]
!//                  and electron temperature Te[eV]
!//=================================================

  module ModSaha
    implicit NONE
    include 'CRASH_definitions.h'
    real,  parameter ::  &
        Err =  1.00e-012  ,         &  !*f9 9 //09        //009 work up to 4-5 poin
        Pi  = (3.1415926) ,        &
        Pi2 = (0.5/Pi)    ,       &
         Tev     = (1./11610.),    &   !   ( 1 K in eV)
         Hmax    = (6.6200e-027),   &  !   erg*cek Planck
         HmaxT   = (6.6200e-027*Pi2)   !   the same/2Pi


     real:: &  
           Neplasma  , &
           Telectrons, &  
	   NatomII,    &  ! that is sum conc (atoms+ions) at the ioniz.reg.
	   Z     ,     &
	   Ci(0:99  ),  & !  Ci[0]-neutrals AFTER ionization
           Ui(99)
     integer :: nCi = 67 !Xenon, must be a var for other elements
   end module ModSaha
!........................................................................
!==============================================/
!===========================================/
!========================================/
!=====================================/
!==================================/
 subroutine      get_ionization_equilibrium
 
 
 !****   All   input PARAMETERS 
 !****   MUST   BE transformed 
 !****   before  CONCAFTER call

 end subroutine  get_ionization_equilibrium

!=================================\
!===================================\
!=====================================\
!=======================================\
!=========================================\
!===========================================\
!=============================================\
!===============================================\
 real function   Conc( N0after )
    use ModSaha
    implicit NONE
    real, intent(in) :: N0after !// conc neautrals

   integer::  k=1     ,  &   
              iter    ,  &
              iter2=0 ,  &
              iZer =1 ,  &   !  start of slider   zb+{|+0,|+1,|+2,|+3,|+4}
              zb

    real ::       &
        NORM  ,      &!= N0after,         &
	te, t32 ,    &!= Telectrons ,     	t32  = te*sqrt(te) 
	a4,          &!
	c1   , c11, c2, & 
	bc   , d   ,dc   ,            &
        a11  ,a21  ,a31  ,a41  ,a51, &
        a12  ,a22  ,a32  ,a42  ,a52, &
	a7   , b7  , c7  , d7  , e7, & 
	s    , s1  , s2  , st  ,     &
        U1,U2,U3,U4,U5      ;
	
	s     = 0.0
	Ci(0) = N0after



     te   = Telectrons
     t32  = te*SQRT(te)      
     NORM = N0after
     a4   =  6.050e21*t32/NORM 


  !debuT:       
           write (*,*) 'Conc: Norm=', NORM


     izer =1  ! for the start of SLIDing
!......................................... make SLIDEr  below ............



        Ci(5)=10.  ; Ci(4)=1.;             !only 11 Ui now
 SLIDER: do while ( izer <= (11-4)    ) 
   
  write(*,*)'SLIDer====== iZer=',izer


	 c11=1.0
         c2 =1.0

     
         do    k = 1,  nCi 
            Ci(k)= 0.
         end do

         U1 = Ui(izer )	
         U2 = Ui(izer+1)	
         U3 = Ui(izer+2)	
         U4 = Ui(izer+3)	
         U5 = Ui(izer+4)	

	st=U1/te;  a11= a4*exp(-st ) !*1./2. ;
	st=U2/te;  a21= a4*exp(-st ) !*6./1. ;
	st=U3/te;  a31= a4*exp(-st ) !*9./6. ;
	st=U4/te;  a41= a4*exp(-st ) !*6./1. ;
 	st=U5/te;  a51= a4*exp(-st ) !*6./1. ;

	c1=1.3*sqrt(c2*a11);

            iter2= 1
            dc   = 1.
            zb   = iZer
      OKNO: do while ( (dc>err ).AND.(iter2 < 30)  )!zzzzzzzzzzzzzzzzzzzzzzzz 

            a12 = a11   !// *exp(dui   ); 
	    a22 = a21   !// *exp(dui*2.);
	    a32 = a31   !// *exp(dui*3.); 
	    a42 = a41   !// *exp(dui*4.);
            a52 = a51   !<================ plasma NonIDeality would be  later

	      iter =0; 
              d    =1.;
              
	     NEWTON : do while(  (iter <30  ).AND.( d*d > Err)) !{ xxxxxxxxxxxx
			bc = 1. / c1;			
			a7 =        a12*bc;
                        b7 = (a7) * a22*bc; 
                        c7 = (b7) * a32*bc;
		        d7 = (c7) * a42*bc;
		        e7 = (d7) * a52*bc;

                        s  = a7* zb            +b7*(zb+1.)           &  
                            +c7*(zb+2.)        +d7*(zb+3.)           &
                            +e7*(zb+4.)                              
		        s1 = a7* zb*zb          +b7*(zb+1.)*(zb+1.)  &
                            +c7*(zb+2.)*(zb+2.) +d7*(zb+3.)*(zb+3.)  &
                            +e7*(zb+4.)*(zb+4.) 

			d  =  c1 -c2 * s  ;
			c1 =( c1 *c2 *(s +s1)/(c1 +c2*s1) );
                     iter=iter+1
		 ! }xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
              end do NEWTON

!	    Ce =c1 * NORM

	    s2 = a7+b7+c7+d7+e7;         
	    Z  = s /s2;
        
	   dc  = (c1 -c11)/(c1+c11);   
	   dc  = dc*dc;
	   c11 = c1; 
           iter2 = iter2+1  
	!}zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
        end do OKNO


	  Ci(izer  ) = a7* NORM
	  Ci(izer+1) = b7* NORM
	  Ci(izer+2) = c7* NORM
	  Ci(izer+3) = d7* NORM
	  Ci(izer+4) = e7* NORM
          Neplasma   = c1* NORM ;


	    c11=0.0;  
       BCE: do  k  = 0, nCi  !.................................................   
            c11= c11+Ci(k);
!debuT:  
            if( 1.00<Ci(k))then
               write(*,*)'#### ',k,'  summa=',c11 ! ,' 0>', Ci(0), ' <| ', NORM
            end if 
        end do BCE !............................................................


!debuT:
        write(*,'(a,10(1pg12.6,2x) )') 'Conc: Ci=',    &
                    Ci(0), Ci(izer), Ci(izer+1), Ci(izer+2), Ci(izer+3), Ci(izer+4)

     if( Ci(izer+4) < Ci(izer+0) ) then
!debuT:
                  write(*,*)'####### @ exit  Izer=',izer , &
                             'Ci(',izer+0,')=',Ci(izer+0), &
                             'Ci(',izer+1,')=',Ci(izer+1), &
                             'Ci(',izer+2,')=',Ci(izer+2), &
                             'Ci(',izer+3,')=',Ci(izer+3), &
                             'Ci(',izer+4,')=',Ci(izer+4), &
                             '######'                     
                  exit
     else 
                  izer=izer+1
     end if


 end do SLIDER 

      Conc =c11

  return 
end function Conc
!.......................................












!*//======================================
!// exp( 2.303)=10. <~> exp(57.565)=1.e25;




 subroutine    ConcNafter  ! // for given summ(No+Niz) on boundary II, calc. Ci[z]
    use ModSaha
    implicit NONE

    integer          ::    iter 
    real ::               &
	         nsum  ,   bs,        &  
                 x3, x1,x2,f1,f2, df, &
                 Conc
  
        nsum =     NatomII     
        bs   = 1.0/NatomII 
        Ci(0)= nsum

! debuT:       
        write (*,*) 'ConcN: NatomII=', NatomII, ' Te=', Telectrons

               



	x1= Ci(0)*1.001; 
	    f1=( Conc(x1) -nsum)*bs;

	x2=Ci(0)*0.011; 
	    f2=( Conc(x2) -nsum)*bs;

         iter = 0

   neutral:   do while(  (f2*f2 > Err).AND.( iter < 55 )   )
          !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
                    df = (f2 -f1)
            if (0.==df) then
                    df = Err
            end if
		     x3 = (-x2*f1 +  x1*f2 )/df !   (f2 - f1);
            if(      x3 < 1.00d-22 ) then
                     x3 = 1.00d-22  ;
	    else if( x3 > 1.00d+33 ) then  
                     x3 = 1.00d+33  ;
            end if
	    if( f1*f1 > f2*f2 ) then 
                x1=x2;  f1=f2; 
            end if
	    x2 = x3;
	    f2 = ( Conc(x3)-nsum )*bs;
	  if( 0.0 == f2  )  exit
          iter =iter +1

!debuT:
       write (*,*) '...f2=', DABS (f2) ,' <2-1> ', DABS(f1) 

      end do  neutral   
  return
end subroutine ConcNafter





