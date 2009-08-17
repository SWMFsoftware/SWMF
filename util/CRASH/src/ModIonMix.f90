module CRASH_ModIonMix
  !This is an interface to some methods of the code iomix.f
  ! IONMIX, A CODE FOR COMPUTING THE EQUATION OF STATE AND   
  ! RADIATIVE PROPERTIES OF LTE AND NON-LTE PLASMAS.  J.J. MACFARLANE.  
  ! REF. IN COMP. PHYS. COMMUN. 56 (1989) 259  
  implicit none
  SAVE
  real,parameter:: con(9)=(/&
       0.0,&
       1.0e-10,& !min. species concentration to compute bb and bf transitions
       1.0e-10,& !min. ionization concentration to compute bb and bf transitions
       1.0e-10,& !min. concentration of an atomic state 4 bb and bf transitions
       1.0e10, & !range, in # of line widths (fwhm)to compute contribution from bb
       10.0,   & !width, in # of line widths (fwhm),of line core           
       10.0,   & ! abs. cfs. are weighted by Planck fn. when 1/con7 < hv/kT < con7
       0.0,    &
       1.0     /)! multiplier for mesh pt spacing about lines
contains
 
  !===============================================
  real function gint ( igint,xStart,xEnd )
    real,intent(in) :: xStart, xEnd
    integer,intent(in) :: iGInt
    !                                                  
    ! ... this function returns the integral <from xmin to xmax> of  
    !     GINTn (where n = 1 thru 6)                                 

    logical::DoTest                                                   

    real,parameter :: xMin =  .005 , xMax =  30. , DxLog = .17754 , &   
         xmnlog =  -5.2983   
    !Misc
    real:: xLog, XII, EX
    integer :: II 
    !----------------                                               



    select case(iGInt)                                
    case(1)                                                                  
       gint = gint1(xStart) - gint1(xEnd)                                  
    case(2)                                                         
       gint = gint2(xStart) - gint2(xEnd)                                  
    case(3)                                                         
       gint = gint3(xStart) - gint3(xEnd)                                  
    case(4)                                                        
       gint = gint4(xStart) - gint4(xEnd)                              
    case(5)
       gint = gint5(xStart) - gint5(xEnd)                               
    case(6)                                                         
       gint = gint6(xStart) - gint6(xEnd)                                 
    case default
       call CON_stop('Error in g_int')                                                        
    end select
  contains

    real function gint1 ( x )                                              
      real,intent(in) :: x                                                                           
      ! ... this routine returns the integral (from 0 to x) of the            
      !     function:                                                         
      !                                                                       
      !                     x_                                                
      !           gint1  = _/  dy exp(-y)                                     
      !                   0                                                   
      !                                                                       
      !---------------------                                                  

      gint1 = -exp( -x )                                                


    end function gint1
    !==================
    real function gint2 ( x )                                         
      real,intent(in) :: x   

      ! ... this routine returns the integral (from 0 to x) of the            
      !     function:                                                         
      !                                                                       
      !                     x_                                                
      !           gint2  = _/  dy y**3 exp(-y)                                
      !                   0                                                   
      !                                                                       
      gint2 = -exp( -x ) * ( x**3 + 3.*x**2 + 6.*x + 6. )               
    end function gint2
    !==================
    real function gint3 ( x )                                              
      real,intent(in) :: x   

      ! ... this routine returns the integral (from 0 to x) of the            
      !     function:                                                         
      !                                                                       
      !                     x_                                                
      !           gint3  = _/  dy y**4 exp(-y) / ( 1-exp(-y) )**3             
      !                   0                                                   
      !                                                                       
      ! ... the table values are the logarithms of the integrals for          
      !     equally spaced values of  log(x).                                 
      !  
      real,parameter:: gtable(50) = (/&                                              
           -1.1288E+01,-1.0933E+01,-1.0577E+01,-1.0222E+01,-9.8661E+00,& 
           -9.5103E+00,-9.1545E+00,-8.7984E+00,-8.4422E+00,-8.0858E+00,&
           -7.7292E+00,-7.3722E+00,-7.0149E+00,-6.6571E+00,-6.2988E+00,&
           -5.9399E+00,-5.5803E+00,-5.2199E+00,-4.8584E+00,-4.4958E+00,&
           -4.1318E+00,-3.7661E+00,-3.3986E+00,-3.0290E+00,-2.6568E+00,&
           -2.2819E+00,-1.9038E+00,-1.5224E+00,-1.1374E+00,-7.4873E-01,&
           -3.5677E-01, 3.7811E-02, 4.3361E-01, 8.2825E-01, 1.2180E+00,&
           1.5975E+00, 1.9592E+00, 2.2940E+00, 2.5913E+00, 2.8406E+00,&
           3.0339E+00, 3.1686E+00, 3.2499E+00, 3.2903E+00, 3.3058E+00,&
           3.3101E+00, 3.3109E+00, 3.3110E+00, 3.3110E+00, 3.3110E+00/)   
      !---------------------------------
      if ( x .lt. xmin ) then
         ! ...    analytic form using expansion of exp(x)    
         gint3 = x**2 / 2.                                              
      else if ( x .gt. xmax*0.9999 ) then                               

         ! ...    the integral as x -> infinity                                  
         gint3 = exp( 3.3110 )                                         

      else                                                           

         ! ...    otherwise, interpolate using table                             
         xlog = log( x )                                                
         ii = 1 + ( xlog-xmnlog ) / dxlog                               
         xii = xmnlog + ( ii-1 )*dxlog                                  
         ex = gtable(ii) + ( xlog-xii ) *  &                             
              ( gtable(ii+1) - gtable(ii) ) / dxlog              
         gint3 = exp( ex )                                             

      end if
    end function gint3
    !=================
    real function gint4 ( x )
      real,intent(in) :: x                           
      ! ... this routine returns the integral (from 0 to x) of the    
      !     function:                                                  
      !                                                                
      !                     x_                                         
      !           gint4  = _/  dy y**7 exp(-y) / ( 1-exp(-y) )**3      
      !                   0                                                   
      !                                                                
      ! ... the table values are the logarithms of the integrals for  
      !     equally spaced values of  log(x).                          
      !                                                                 

      real,parameter:: gtable(50)= (/&                                              
           -2.8099E+01,-2.7211E+01,-2.6323E+01,-2.5434E+01,-2.4546E+01,&   
           -2.3657E+01,-2.2769E+01,-2.1880E+01,-2.0991E+01,-2.0101E+01,&     
           -1.9212E+01,-1.8322E+01,-1.7431E+01,-1.6540E+01,-1.5648E+01,&     
           -1.4756E+01,-1.3863E+01,-1.2968E+01,-1.2073E+01,-1.1176E+01,&     
           -1.0277E+01,-9.3763E+00,-8.4734E+00,-7.5678E+00,-6.6594E+00,&     
           -5.7477E+00,-4.8324E+00,-3.9134E+00,-2.9904E+00,-2.0638E+00,&     
           -1.1342E+00,-2.0272E-01, 7.2834E-01, 1.6554E+00, 2.5731E+00,&     
           3.4736E+00, 4.3462E+00, 5.1772E+00, 5.9497E+00, 6.6448E+00,&     
           7.2434E+00, 7.7287E+00, 8.0902E+00, 8.3287E+00, 8.4606E+00,&     
           8.5172E+00, 8.5343E+00, 8.5375E+00, 8.5379E+00, 8.5379E+00/)
      !------------------------ 
      if ( x .lt. xmin ) then                                           

         ! ...    analytic form using expansion of exp(x)                       
         gint4 = x**5 / 5.                                              

      else if ( x .gt. xmax*0.9999 ) then                               

         ! ...    the integral as x -> infinity                                  
         gint4 = exp( 8.5379 )                                          

      else                                                              

         ! ...    otherwise, interpolate using table                            
         xlog = log( x )                                                
         ii = 1 + ( xlog-xmnlog ) / dxlog                               
         xii = xmnlog + ( ii-1 )*dxlog                                  
         ex = gtable(ii) + ( xlog-xii ) *  &                             
              ( gtable(ii+1) - gtable(ii) ) / dxlog              
         gint4 = exp( ex )
      end if
    end function gint4
    !==================
    real  function gint5 ( x )                                          
      real, intent(in):: x                                                 
      ! ... this routine returns the integral (from 0 to x) of the           
      !     function:                                                         
      !                                                                       
      !                     x_                                                
      !           gint5  = _/  dy y**3 / ( exp(y)-1 )                         
      !                   0                                                   
      !                                                                       
      ! ... the table values are the logarithms of the integrals for          
      !     equally spaced values of  log(x).                                 
      !                                                                                                                  
      real,parameter:: gtable(50) = (/&                                                     
           -1.6995E+01,-1.6463E+01,-1.5931E+01,-1.5399E+01,-1.4867E+01,&    
           -1.4335E+01,-1.3803E+01,-1.3272E+01,-1.2740E+01,-1.2209E+01,&     
           -1.1678E+01,-1.1148E+01,-1.0618E+01,-1.0088E+01,-9.5594E+00,&     
           -9.0312E+00,-8.5039E+00,-7.9775E+00,-7.4524E+00,-6.9289E+00,&     
           -6.4070E+00,-5.8874E+00,-5.3703E+00,-4.8563E+00,-4.3460E+00,&     
           -3.8402E+00,-3.3399E+00,-2.8462E+00,-2.3605E+00,-1.8846E+00,&     
           -1.4204E+00,-9.7072E-01,-5.3858E-01,-1.2780E-01, 2.5710E-01,&     
           6.1092E-01, 9.2788E-01, 1.2021E+00, 1.4283E+00, 1.6031E+00,&     
           1.7268E+00, 1.8043E+00, 1.8456E+00, 1.8634E+00, 1.8693E+00,&     
           1.8706E+00, 1.8708E+00, 1.8709E+00, 1.8709E+00, 1.8709E+00 /)
      !========================                                                                  
      if ( x .lt. xmin ) then                                           

         ! ...    analytic form using expansion of exp(x)                        
         gint5 = x**3 / 3.                                             

      else if ( x .gt. xmax*0.9999 ) then                               

         ! ...    the integral as x -> infinity                                  
         gint5 = exp( 1.8709 )                                         

      else                                                              

         ! ...    otherwise, interpolate using table                             
         xlog = log( x )                                                
         ii = 1 + ( xlog-xmnlog ) / dxlog                              
         xii = xmnlog + ( ii-1 )*dxlog                                 
         ex = gtable(ii) + ( xlog-xii ) * &                             
              ( gtable(ii+1) - gtable(ii) ) / dxlog              
         gint5 = exp( ex )                                             

      endif

    end function gint5
    !-----------------------------------------------------------------
    real function gint6 ( x )
      real,intent(in) :: x                                             

      ! ... this routine returns the integral (from 0 to x) of the            
      !     function:                                                         
      !                                                                       
      !                     x_                                                
      !           gint6  = _/  dy y**4 exp(-y) / ( 1-exp(-y) )**2             
      !                   0                                                   
      !                                                                       
      ! ... the table values are the logarithms of the integrals for          
      !     equally spaced values of  log(x).                                 
      !                                                                       


      real,parameter:: gtable(50)= (/   &                                                 
           -1.6994E+01,-1.6461E+01,-1.5928E+01,-1.5396E+01,-1.4863E+01,&
           -1.4330E+01,-1.3798E+01,-1.3265E+01,-1.2733E+01,-1.2200E+01,&
           -1.1667E+01,-1.1135E+01,-1.0602E+01,-1.0070E+01,-9.5370E+00,&
           -9.0045E+00,-8.4720E+00,-7.9395E+00,-7.4071E+00,-6.8748E+00,&
           -6.3426E+00,-5.8106E+00,-5.2789E+00,-4.7476E+00,-4.2169E+00,&
           -3.6869E+00,-3.1581E+00,-2.6309E+00,-2.1060E+00,-1.5843E+00,&
           -1.0671E+00,-5.5646E-01,-5.4768E-02, 4.3443E-01, 9.0650E-01,&
           1.3554E+00, 1.7736E+00, 2.1524E+00, 2.4822E+00, 2.7542E+00,&
           2.9625E+00, 3.1063E+00, 3.1926E+00, 3.2353E+00, 3.2516E+00,&
           3.2562E+00, 3.2571E+00, 3.2572E+00, 3.2572E+00, 3.2572E+00 /)    
      !---------------------------

      if ( x .lt. xmin ) then                     
         ! ...    analytic form using expansion of exp(x)               
         gint6 = x**3 / 3.                                              

      else if ( x .gt. xmax*0.9999 ) then                               

         ! ...    the integral as x -> infinity                                  
         gint6 = exp( 3.2572 )                                          

      else                                                              

         ! ...    otherwise, interpolate using table                             
         xlog = log( x )                                                
         ii = 1 + ( xlog-xmnlog ) / dxlog                               
         xii = xmnlog + ( ii-1 )*dxlog                                  
         ex = gtable(ii) + ( xlog-xii ) * &                             
              ( gtable(ii+1) - gtable(ii) ) / dxlog               
         gint6 = exp( ex )
      end if
    end function gint6
  end function gint
end module CRASH_ModIonMix
