C IRIFUN.FOR
C
C**************************************************************  
C********** INTERNATIONAL REFERENCE IONOSPHERE ****************  
C**************************************************************  
C****************  IRIFU9, DECEMBER 1988  *********************  
C****************  FUNCTIONS,SUBROUTINES  *********************
C**************************************************************
C** NE: 	XE1, DXE1N, XE2, XE3, XE4, XE5, XE6, XE
C** TE: 	TEBA, SPHARM, ELTE, TEDE
C** TN/TI:	TUNCAL, TN, DTNDH, TI, TEDER
C** NI:		RPID, RDHHE, RDNO, KOEFP1, KOEFP2, KOEFP3, SUFE
C** PEAKS:	F2OUT, HMF2ED, FOF1ED, FOEEDI, XMDED, GAMMA1
C** MAG. FIELD: GGM, FIELDG
C** FUNCTIONS: 	REGFA1, TAL
C** TIME:	SOCO,HPOL, MODA
C**************************************************************  
C  
C**************************************************************  
C***  -------------------ADDRESSES------------------------  ***
C***  I  PROF. K. RAWER             DR. D. BILITZA       I  ***
C***  I  HERRENSTR. 43              GSFC CODE 633        I  ***
C***  I  7801 MARCH 1               GREENBELT MD 20771   I  ***
C***  I  F.R.G.                     USA                  I  ***
C***  I                             TEL. (301)286-9536   I  ***  
C***  ----------------------------------------------------  ***
C**************************************************************  
C**************************************************************  
C        
C*************************************************************   
C*************** ELECTRON DENSITY ****************************   
C*************************************************************   
C
C
      FUNCTION XE1(H)    
C REPRESENTING ELECTRON DENSITY(M-3) IN THE TOPSIDE IONOSPHERE   
C (H=HMF2....1000 KM) BY HARMONIZED BENT-MODEL ADMITTING 
C VARIABILITY OFGLOBAL PARAMETER ETA,ZETA,BETA,DELTA WITH        
C GEOM. LATITUDE, SMOOTHED SOLAR FLUX AND CRITICAL FREQUENCY     
C (SEE MAIN PROGRAM).    
C [REF.:K.RAWER,S.RAMAKRISHNAN,1978]     
      REAL NMF2       
      COMMON/BLOCK1/HMF2,NMF2,HMF1/BLO10/BETA,ETA,DELTA,ZETA                    
      DXDH=(1000.-HMF2)/700.
      X=(H-HMF2)/DXDH+300.0-DELTA
      E1=EXP((X-394.5)/BETA) 
      E2=EXP((X-300.0)/100.)        
      Y=DXDH*(BETA*ETA*ALOG((1.+E1)/(1.+EXP((-94.5-DELTA)/BETA)))               
     &+ZETA*(100.*ALOG((1.+E2)/(1.+EXP(-DELTA/100.)))-X+300.-DELTA))            
      XE1=NMF2*EXP(-Y)                             
      RETURN          
      END             
C 
C
      FUNCTION DXE1N(H)                            
C LOGARITHMIC DERIVATIVE OF FUNCTION XE1 (KM-1).   
      REAL NMF2       
      COMMON/BLOCK1/HMF2,NMF2,HMF1/BLO10/BETA,ETA,DELTA,ZETA                    
      X=(H-HMF2)/(1000.0-HMF2)*700.0+300.0-DELTA   
      DXE1N=(ZETA*(1.-1./(1.+EXP((300.-X)/100.)))- 
     &ETA*1./(1.+EXP((394.5-X)/BETA)))             
      RETURN          
      END             
C
C
      REAL FUNCTION XE2(H)                         
C ELECTRON DENSITY FOR THE BOTTOMSIDE F-REGION (HMF1...HMF2).                   
      REAL NMF2       
      COMMON/BLOCK1/HMF2,NMF2,HMF1/BLOCK2/B0,B1,C1,HZ,T,HST,STR                 
      X=(HMF2-H)/B0
      XE2=NMF2*EXP(-X**B1)/COSH(X)                 
      RETURN          
      END             
C
C
      REAL FUNCTION XE3(H)                         
C ELECTRON DENSITY FOR THE F1-LAYER (HZ.....HMF1). 
      REAL NMF2       
      COMMON/BLOCK1/HMF2,NMF2,HMF1/BLOCK2/B0,B1,C1,HZ,T,HST,STR
      XE3=XE2(H)+NMF2*C1*SQRT(ABS(HMF1-H)/B0)      
      RETURN          
      END             
C
C
      REAL FUNCTION XE4(H)                         
C ELECTRON DENSITY FOR THE INDERMEDIUM REGION (HEF..HZ).                        
      COMMON/BLOCK2/B0,B1,C1,HZ,T,HST,STR/BLOCK3/HDX,HME,  
     &  XME,HMD,XMD,HEF,D1,XKK,FP30,FP3U,FP1,FP2 
      IF(HST.LT.0.) GOTO 100                       
      XE4=XE3(HZ+T/2.0-SIGN(1.0,T)*SQRT(T*(HZ-H+T/4.0)))                        
      RETURN          
100   XE4=XME+T*(H-HEF)                            
      RETURN          
      END             
C
C
      REAL FUNCTION XE5(H)                         
C ELECTRON DENSITY FOR THE E AND VALLEY REGION (HME..HEF).   
      LOGICAL NIGHT   
      COMMON/BLOCK3/HDX,HME,XNME,HMD,XNMD,HEF,D1,XKK,FP30,
     &  FP3U,FP1,FP2/BLOCK6/NIGHT,E(4)                        
      T3=H-HME        
      T1=T3*T3*(E(1)+T3*(E(2)+T3*(E(3)+T3*E(4))))  
      IF(NIGHT) GOTO 100                           
      XE5=XNME*(1+T1)  
      RETURN          
100   XE5=XNME*EXP(T1)                              
      RETURN          
      END             
C
C
      REAL FUNCTION XE6(H)                         
C ELECTRON DENSITY FOR THE D REGION (HA...HME).    
      COMMON/BLOCK3/HDX,HME,XNME,HMD,XNMD,HEF,D1,XKK,FP30,
     &  FP3U,FP1,FP2
      IF(H.GT.HDX) GOTO 100                        
      Z=H-HMD         
      FP3=FP3U        
      IF(Z.GT.0.0) FP3=FP30                        
      XE6=XNMD*EXP(Z*(FP1+Z*(FP2+Z*FP3)))           
      RETURN          
100   Z=HME-H         
      XE6=XNME*EXP(-D1*Z**XKK)
      RETURN          
      END             
C
C
      REAL FUNCTION XE(H)                          
C ELECTRON DENSITY BEETWEEN HA(KM) AND 1000 KM     
C SUMMARIZING PROCEDURES  NE1....6;                
      COMMON/BLOCK1/HMF2,XNMF2,HMF1/BLOCK2/B0,B1,C1,HZ,T,
     &  HST,STR/BLOCK3/HDX,HME,XNME,HMD,XNMD,HEF,D1,XKK,
     &  FP30,FP3U,FP1,FP2
      IF(H.LT.HMF2) GOTO 100                       
      XE=XE1(H)     
      RETURN          
100   IF(H.LT.HMF1) GOTO 300                       
      XE=XE2(H)       
      RETURN          
300   IF(H.LT.HZ) GOTO 400                         
      XE=XE3(H)       
      RETURN          
400   IF(H.LT.HEF) GOTO 500                        
      XE=XE4(H)       
      RETURN          
500   IF(H.LT.HME) GOTO 600                        
      XE=XE5(H)       
      RETURN          
600   XE=XE6(H)       
      RETURN          
      END             
C                     
C**********************************************************                     
C***************** ELECTRON TEMPERATURE ********************                    
C**********************************************************                     
C
      SUBROUTINE TEBA(DIPL,SLT,IS,TE) 
C CALCULATES ELECTRON TEMPERATURES TE(1) TO TE(4) AT ALTITUDES
C 300, 400, 1400 AND 3000 KM FOR DIP-LATITUDE DIPL/DEG AND 
C LOCAL SOLAR TIME SLT/H USING THE BRACE-THEIS-MODELS (J. ATMOS.
C TERR. PHYS. 43, 1317, 1981); IS=1 EQUINOX, IS=2 SUMMER.
C ALSO CALCULATED ARE THE TEMPERATURES AT 400 KM ALTITUDE FOR
C MIDNIGHT (TE(5)) AND NOON (TE(6)).   
      DIMENSION C(4,2,81),A(82),TE(6)
      COMMON/CONST/UMR
      DATA (C(1,1,J),J=1,81)/                      
     &.3100E1,-.3215E-2,.2440E+0,-.4613E-3,-.1711E-1,.2605E-1,                  
     &-.9546E-1,.1794E-1,.1270E-1,.2791E-1,.1536E-1,-.6629E-2,                  
     &-.3616E-2,.1229E-1,.4147E-3,.1447E-2,-.4453E-3,-.1853,                    
     &-.1245E-1,-.3675E-1,.4965E-2,.5460E-2,.8117E-2,-.1002E-1,                 
     &.5466E-3,-.3087E-1,-.3435E-2,-.1107E-3,.2199E-2,.4115E-3,                 
     &.6061E-3,.2916E-3,-.6584E-1,.4729E-2,-.1523E-2,.6689E-3,                  
     &.1031E-2,.5398E-3,-.1924E-2,-.4565E-1,.7244E-2,-.8543E-4,                 
     &.1052E-2,-.6696E-3,-.7492E-3,.4405E-1,.3047E-2,.2858E-2,                  
     &-.1465E-3,.1195E-2,-.1024E-3,.4582E-1,.8749E-3,.3011E-3,                  
     &.4473E-3,-.2782E-3,.4911E-1,-.1016E-1,.27E-2,-.9304E-3,                   
     &-.1202E-2,.2210E-1,.2566E-2,-.122E-3,.3987E-3,-.5744E-1,                  
     &.4408E-2,-.3497E-2,.83E-3,-.3536E-1,-.8813E-2,.2423E-2,                   
     &-.2994E-1,-.1929E-2,-.5268E-3,-.2228E-1,.3385E-2,                         
     &.413E-1,.4876E-2,.2692E-1,.1684E-2/          
      DATA (C(1,2,J),J=1,81)/.313654E1,.6796E-2,.181413,.8564E-1,               
     &-.32856E-1,-.3508E-2,-.1438E-1,-.2454E-1,.2745E-2,.5284E-1,               
     &.1136E-1,-.1956E-1,-.5805E-2,.2801E-2,-.1211E-2,.4127E-2,                 
     &.2909E-2,-.25751,-.37915E-2,-.136E-1,-.13225E-1,.1202E-1,                 
     &.1256E-1,-.12165E-1,.1326E-1,-.7123E-1,.5793E-3,.1537E-2,                 
     &.6914E-2,-.4173E-2,.1052E-3,-.5765E-3,-.4041E-1,-.1752E-2,                
     &-.542E-2,-.684E-2,.8921E-3,-.2228E-2,.1428E-2,.6635E-2,-.48045E-2,        
     &-.1659E-2,-.9341E-3,.223E-3,-.9995E-3,.4285E-1,-.5211E-3,                 
     &-.3293E-2,.179E-2,.6435E-3,-.1891E-3,.3844E-1,.359E-2,-.8139E-3,          
     &-.1996E-2,.2398E-3,.2938E-1,.761E-2,.347655E-2,.1707E-2,.2769E-3,         
     &-.157E-1,.983E-3,-.6532E-3,.929E-4,-.2506E-1,.4681E-2,.1461E-2,           
     &-.3757E-5,-.9728E-2,.2315E-2,.6377E-3,-.1705E-1,.2767E-2,                 
     &-.6992E-3,-.115E-1,-.1644E-2,.3355E-2,-.4326E-2,.2035E-1,.2985E-1/        
      DATA (C(2,1,J),J=1,81)/.3136E1,.6498E-2,.2289,.1859E-1,-.3328E-1,         
     &-.4889E-2,-.3054E-1,-.1773E-1,-.1728E-1,.6555E-1,.1775E-1,                
     &-.2488E-1,-.9498E-2,.1493E-1,.281E-2,.2406E-2,.5436E-2,-.2115,            
     &.7007E-2,-.5129E-1,-.7327E-2,.2402E-1,.4772E-2,-.7374E-2,                 
     &-.3835E-3,-.5013E-1,.2866E-2,.2216E-2,.2412E-3,.2094E-2,.122E-2           
     &,-.1703E-3,-.1082,-.4992E-2,-.4065E-2,.3615E-2,-.2738E-2,                 
     &-.7177E-3,.2173E-3,-.4373E-1,-.375E-2,.5507E-2,-.1567E-2,                 
     &-.1458E-2,-.7397E-3,.7903E-1,.4131E-2,.3714E-2,.1073E-2,                  
     &-.8991E-3,.2976E-3,.2623E-1,.2344E-2,.5608E-3,.4124E-3,.1509E-3,          
     &.5103E-1,.345E-2,.1283E-2,.7238E-3,-.3464E-4,.1663E-1,-.1644E-2,          
     &-.71E-3,.5281E-3,-.2729E-1,.3556E-2,-.3391E-2,-.1787E-3,.2154E-2,         
     &.6476E-2,-.8282E-3,-.2361E-1,.9557E-3,.3205E-3,-.2301E-1,                 
     &-.854E-3,-.1126E-1,-.2323E-2,-.8582E-2,.2683E-1/                          
      DATA (C(2,2,J),J=1,81)/.3144E1,.8571E-2,.2539,.6937E-1,-.1667E-1,         
     &.2249E-1,-.4162E-1,.1201E-1,.2435E-1,.5232E-1,.2521E-1,-.199E-1,          
     &-.7671E-2,.1264E-1,-.1551E-2,-.1928E-2,.3652E-2,-.2019,.5697E-2,          
     &-.3159E-1,-.1451E-1,.2868E-1,.1377E-1,-.4383E-2,.1172E-1,                 
     &-.5683E-1,.3593E-2,.3571E-2,.3282E-2,.1732E-2,-.4921E-3,-.1165E-2         
     &,-.1066,-.1892E-1,.357E-2,-.8631E-3,-.1876E-2,-.8414E-4,.2356E-2,         
     &-.4259E-1,-.322E-2,.4641E-2,.6223E-3,-.168E-2,-.1243E-3,.7393E-1,         
     &-.3143E-2,-.2362E-2,.1235E-2,-.1551E-2,.2099E-3,.2299E-1,.5301E-2         
     &,-.4306E-2,-.1303E-2,.7687E-5,.5305E-1,.6642E-2,-.1686E-2,                
     &.1048E-2,.5958E-3,.4341E-1,-.8819E-4,-.333E-3,-.2158E-3,-.4106E-1         
     &,.4191E-2,.2045E-2,-.1437E-3,-.1803E-1,-.8072E-3,-.424E-3,                
     &-.26E-1,-.2329E-2,.5949E-3,-.1371E-1,-.2188E-2,.1788E-1,                  
     &.6405E-3,.5977E-2,.1333E-1/                  
      DATA (C(3,1,J),J=1,81)/.3372E1,.1006E-1,.1436,.2023E-2,-.5166E-1,         
     &.9606E-2,-.5596E-1,.4914E-3,-.3124E-2,-.4713E-1,-.7371E-2,                
     &-.4823E-2,-.2213E-2,.6569E-2,-.1962E-3,.3309E-3,-.3908E-3,                
     &-.2836,.7829E-2,.1175E-1,.9919E-3,.6589E-2,.2045E-2,-.7346E-2             
     &,-.89E-3,-.347E-1,-.4977E-2,.147E-2,-.2823E-5,.6465E-3,                   
     &-.1448E-3,.1401E-2,-.8988E-1,-.3293E-4,-.1848E-2,.4439E-3,                
     &-.1263E-2,.317E-3,-.6227E-3,.1721E-1,-.199E-2,-.4627E-3,                  
     &.2897E-5,-.5454E-3,.3385E-3,.8432E-1,-.1951E-2,.1487E-2,                  
     &.1042E-2,-.4788E-3,-.1276E-3,.2373E-1,.2409E-2,.5263E-3,                  
     &.1301E-2,-.4177E-3,.3974E-1,.1418E-3,-.1048E-2,-.2982E-3,                 
     &-.3396E-4,.131E-1,.1413E-2,-.1373E-3,.2638E-3,-.4171E-1,                  
     &-.5932E-3,-.7523E-3,-.6883E-3,-.2355E-1,.5695E-3,-.2219E-4,               
     &-.2301E-1,-.9962E-4,-.6761E-3,.204E-2,-.5479E-3,.2591E-1,                 
     &-.2425E-2,.1583E-1,.9577E-2/                 
      DATA (C(3,2,J),J=1,81)/.3367E1,.1038E-1,.1407,.3622E-1,-.3144E-1,         
     &.112E-1,-.5674E-1,.3219E-1,.1288E-2,-.5799E-1,-.4609E-2,                  
     &.3252E-2,-.2859E-3,.1226E-1,-.4539E-2,.1310E-2,-.5603E-3,                 
     &-.311,-.1268E-2,.1539E-1,.3146E-2,.7787E-2,-.143E-2,-.482E-2              
     &,.2924E-2,-.9981E-1,-.7838E-2,-.1663E-3,.4769E-3,.4148E-2,                
     &-.1008E-2,-.979E-3,-.9049E-1,-.2994E-2,-.6748E-2,-.9889E-3,               
     &.1488E-2,-.1154E-2,-.8412E-4,-.1302E-1,-.4859E-2,-.7172E-3,               
     &-.9401E-3,.9101E-3,-.1735E-3,.7055E-1,.6398E-2,-.3103E-2,                 
     &-.938E-3,-.4E-3,-.1165E-2,.2713E-1,-.1654E-2,.2781E-2,                    
     &-.5215E-5,.2258E-3,.5022E-1,.95E-2,.4147E-3,.3499E-3,                     
     &-.6097E-3,.4118E-1,.6556E-2,.3793E-2,-.1226E-3,-.2517E-1,                 
     &.1491E-3,.1075E-2,.4531E-3,-.9012E-2,.3343E-2,.3431E-2,                   
     &-.2519E-1,.3793E-4,.5973E-3,-.1423E-1,-.132E-2,-.6048E-2,                 
     &-.5005E-2,-.115E-1,.2574E-1/                 
      DATA (C(4,1,J),J=1,81)/.3574E1,.0,.7537E-1,.0,-.8459E-1,                  
     &0.,-.294E-1,0.,.4547E-1,-.5321E-1,0.,.4328E-2,0.,.6022E-2,                
     &.0,-.9168E-3,.0,-.1768,.0,.294E-1,.0,.5902E-3,.0,-.9047E-2,               
     &.0,-.6555E-1,.0,-.1033E-2,.0,.1674E-2,.0,.2802E-3,-.6786E-1               
     &,.0,.4193E-2,.0,-.6448E-3,.0,.9277E-3,-.1634E-1,.0,-.2531E-2              
     &,.0,.193E-4,.0,.528E-1,.0,.2438E-2,.0,-.5292E-3,.0,.1555E-1               
     &,.0,-.3259E-2,.0,-.5998E-3,.3168E-1,.0,.2382E-2,.0,-.4078E-3              
     &,.2312E-1,.0,.1481E-3,.0,-.1885E-1,.0,.1144E-2,.0,-.9952E-2               
     &,.0,-.551E-3,-.202E-1,.0,-.7283E-4,-.1272E-1,.0,.2224E-2,                 
     &.0,-.251E-2,.2434E-1/                        
      DATA (C(4,2,J),J=1,81)/.3574E1,-.5639E-2,.7094E-1,                        
     &-.3347E-1,-.861E-1,-.2877E-1,-.3154E-1,-.2847E-2,.1235E-1,                
     &-.5966E-1,-.3236E-2,.3795E-3,-.8634E-3,.3377E-2,-.1071E-3,                
     &-.2151E-2,-.4057E-3,-.1783,.126E-1,.2835E-1,-.242E-2,                     
     &.3002E-2,-.4684E-2,-.6756E-2,-.7493E-3,-.6147E-1,-.5636E-2                
     &,-.1234E-2,-.1613E-2,-.6353E-4,-.2503E-3,-.1729E-3,-.7148E-1              
     &,.5326E-2,.4006E-2,.6484E-3,-.1046E-3,-.6034E-3,-.9435E-3,                
     &-.2385E-2,.6853E-2,.151E-2,.1319E-2,.9049E-4,-.1999E-3,                   
     &.3976E-1,.2802E-2,-.103E-2,.5599E-3,-.4791E-3,-.846E-4,                   
     &.2683E-1,.427E-2,.5911E-3,.2987E-3,-.208E-3,.1396E-1,                     
     &-.1922E-2,-.1063E-2,.3803E-3,.1343E-3,.1771E-1,-.1038E-2,                 
     &-.4645E-3,-.2481E-3,-.2251E-1,-.29E-2,-.3977E-3,-.516E-3,                 
     &-.8079E-2,-.1528E-2,.306E-3,-.1582E-1,-.8536E-3,.1565E-3,                 
     &-.1252E-1,.2319E-3,.4311E-2,.1024E-2,.1296E-5,.179E-1/                    
      COLAT=UMR*(90.-DIPL)                    
      AZ=.2618*SLT    
      CALL SPHARM(A,8,8,COLAT,AZ)                  
      DO 2 K=1,4      
	STE=0.          
 	DO 1 I=1,81     
1     		STE=STE+A(I)*C(K,IS,I)                       
2     	TE(K)=10.**STE  

C---------- TEMPERATURE AT 400KM AT MIDNIGHT AND NOON
      DO 4 J=1,2      
      	STE=0.          
      	AZ=.2618*(J-1)*12.                           
      	CALL SPHARM(A,8,8,COLAT,AZ)                  
      	DO 3 I=1,81     
3     		STE=STE+A(I)*C(2,IS,I)                       
4     	TE(J+4)=10.**STE                             
      RETURN          
      END             
C
      SUBROUTINE SPHARM(C,L,M,COLAT,AZ)            
C CALCULATES THE COEFFICIENTS OF THE SPHERICAL HARMONIC                         
C EXPANSION THAT WAS USED FOR THE BRACE-THEIS-MODELS.                           
      DIMENSION C(82)                              
      C(1)=1.         
      K=2             
      X=COS(COLAT)    
      C(K)=X          
      K=K+1           
      DO 10 I=2,L     
      C(K)=((2*I-1)*X*C(K-1)-(I-1)*C(K-2))/I       
10    K=K+1           
      Y=SIN(COLAT)    
      DO 20 MT=1,M    
      CAZ=COS(MT*AZ)  
      SAZ=SIN(MT*AZ)  
      C(K)=Y**MT      
      K=K+1           
      IF(MT.EQ.L) GOTO 16                          
      C(K)=C(K-1)*X*(2*MT+1)                       
      K=K+1           
      IF((MT+1).EQ.L) GOTO 16                      
      DO 15 I=2+MT,L  
      C(K)=((2*I-1)*X*C(K-1)-(I+MT-1)*C(K-2))/(I-MT)                            
15    K=K+1           
16    N=L-MT+1        
      DO 18 I=1,N     
      C(K)=C(K-N)*CAZ                              
      C(K-N)=C(K-N)*SAZ                            
18    K=K+1           
20    CONTINUE        
      RETURN          
      END             
C
C
      REAL FUNCTION ELTE(H)                        
C ELECTRON TEMPERATURE PROFILE BASED ON THE TEMPERATURES AT 120                 
C HMAX,300,400,600,1400,3000 KM ALTITUDE. INBETWEEN CONSTANT                    
C GRADIENT IS ASSUMED. ARGMAX IS MAXIMUM ARGUMENT ALLOWED FOR
C EXP-FUNCTION.
C
      COMMON /BLOTE/AH(7),ATE1,ST(6),D(5)   /ARGEXP/ARGMAX          
C
      SUM=ATE1+ST(1)*(H-AH(1))                     
      DO 1 I=1,5      
      A=(H-AH(I+1))/D(I)
	AA=0.
	IF(A.GE.ARGMAX) AA=A     
      IF(ABS(A).LT.ARGMAX) AA=ALOG(1.+EXP(A))
      B=(AH(1)-AH(I+1))/D(I)
	BB=0.
	IF(A.GE.ARGMAX) BB=B                              
      IF(ABS(B).LT.ARGMAX) BB=ALOG(1.+EXP(B))
1     SUM=SUM+(ST(I+1)-ST(I))*(AA-BB)*D(I)                
      ELTE=SUM        
      RETURN          
      END             
C
C
      FUNCTION TEDE(H,DEN,COV)                     
C ELECTRON TEMEPERATURE MODEL AFTER BRACE,THEIS .  
C FOR NEG. COV THE MEAN COV-INDEX (3 SOLAR ROT.) IS EXPECTED.                   
C DEN IS THE ELECTRON DENSITY IN M-3.              
      Y=1051.+(17.01*H-2746.)*                     
     &EXP(-5.122E-4*H+(6.094E-12-3.353E-14*H)*DEN) 
      ACOV=ABS(COV)   
      YC=1.+(.117+2.02E-3*ACOV)/(1.+EXP(-(ACOV-102.5)/5.))                      
      IF(COV.LT.0.)   
     &YC=1.+(.123+1.69E-3*ACOV)/(1.+EXP(-(ACOV-115.)/10.))                      
      TEDE=Y*YC       
      RETURN          
      END             
C
C                     
C*************************************************************                  
C**************** NEUTRAL AND ION TEMPERATURE ******************                
C*************************************************************                  
C
C
      FUNCTION TUNCAL(COV,XLATI,SD,SLSTA)          
C CALCULATES THE EXOSPHERIC TEMPERATURE FOR COVINGTON-INDEX                     
C COV,LATITUDE XLATI,SOLAR DEKLINATION SD AND LOCAL SOLAR TIME
C ANGLE SLSTA USING CIRA72-MODEL                   
      COMMON/CONST/UMR
      TC=379.0+3.24*COV                            
      ETA=ABS(XLATI-SD)/2.0                        
      THETA=ABS(XLATI+SD)/2.0                      
      H1=COS(ETA*UMR)**2.2                         
      H2=SIN(THETA*UMR)**2.2                       
      TD=1.0+0.3*H1   
      TN=1.0+0.3*H2   
      A=(TD-TN)/TN    
      TAU=SLSTA-37.0+6.0*SIN((SLSTA+43.0)*UMR)     
      X=COS(TAU/2.0*UMR)                           
      TUNCAL=TC*TN*(1.0+A*SIGN(1.0,X)*ABS(X)**3.0) 
      RETURN          
      END             
C
C
      REAL FUNCTION TN(X,TTNX,ATN,CCTN)            
C NEUTRAL TEMPERATURE                              
C D. BILITZA, 1978, NEUTRAL TEMPERATURE PROFILE    
C APPROXIMATING ROUGHLY CIRA72-MODEL;              
      DIMENSION CCTN(3)                            
      Z=X-125.0       
      IF(Z.LE.0.0) GOTO 100                        
      Y=Z**2.5        
      Y=(1.0+(4.5E-6)*Y)*Z                         
      YY=CCTN(1)/ATN  
      Y=YY*Y          
      Y=ATAN(Y)*ATN   
      TN=TTNX+Y       
      RETURN          
100   TN=TTNX+Z*(CCTN(1)+Z*Z*(CCTN(2)+Z*CCTN(3)))  
200   RETURN          
      END             
C
C
      FUNCTION DTNDH(H,ATN,CCTN)                   
C DERIVATIVE OF NEUTRAL TEMPERATURE                
C D.BILITZA, 1978. APPROXIMATE HEIGHT DERIVATIVE OF NEUTRAL                     
C CIRA-72 MODEL TEMPERATURE PROFILE                
      DIMENSION CCTN(3)                            
      Z=H-125.0       
      IF(Z.GT.0.0) GOTO 100                        
      DTNDH=CCTN(1)+Z*Z*(3.0*CCTN(2)+Z*4.0*CCTN(3))                             
      RETURN          
100   H1=CCTN(1)/ATN  
      H2=Z**2.5       
      H3=H1*Z*(1.0+(4.5E-6)*H2)                    
      DTNDH=ATN/(1.0+H3*H3)*H1*(1.0+(15.75E-6)*H2) 
      RETURN          
      END             
C
C
      REAL FUNCTION TI(H)                          
C ION TEMPERATURE FOR HEIGHTS NOT GREATER 1000 KM AND NOT LESS HS               
C EXPLANATION SEE FUNCTION RPID.                   
      REAL 		MM
      COMMON  /BLOCK8/HS,TNHS,XSM(4),MM(5),G(4),M  /ARGEXP/ARGMAX       
      SUM=MM(1)*(H-HS)+TNHS                        
      DO 100 I=1,M-1  
      	A=(H-XSM(I))/G(I)
      	AA=0.0
      	IF(A.GE.ARGMAX) AA=A
      	IF(ABS(A).LT.ARGMAX) AA=ALOG(1.0+EXP(A))
      	B=(HS-XSM(I))/G(I)
      	BB=0.0
      	IF(B.GE.ARGMAX) BB=B
      	IF(ABS(B).LT.ARGMAX) BB=ALOG(1.0+EXP(B))
100   	SUM=SUM+(MM(I+1)-MM(I))*(AA-BB)*G(I)                
      TI=SUM          
      RETURN          
      END             
C
C
      REAL FUNCTION TEDER(H)                       
C THIS FUNCTION ALONG WITH PROCEDURE REGFA1 ALLOWS TO FIND                      
C THE  HEIGHT ABOVE WHICH TN BEGINS TO BE DIFFERENT FROM TI                     
      COMMON/BLOCK5/ZX,TNX,ATN,CTN(3)/BLOCK8/HS,TNHS,XSM(4),
     &  XMM(5),G(4),M               
      TNH=TN(H,TNX,ATN,CTN)                        
      DTDX=DTNDH(H,ATN,CTN)                        
      TEDER=DTDX*(XSM(1)-H)+TNH                    
      RETURN          
      END             
C
C                     
C*************************************************************                  
C************* ION RELATIVE PRECENTAGE DENSITY *****************                
C*************************************************************                  
C
C
      REAL FUNCTION RPID (H, H0, N0, M, ST, ID, XS)                             
C D.BILITZA,1977,THIS ANALYTIC FUNCTION IS USED TO REPRESENT THE                
C RELATIVE PRECENTAGE DENSITY OF ATOMAR AND MOLECULAR OXYGEN IONS.              
C THE M+1 HEIGHT GRADIENTS ST(M+1) ARE CONNECTED WITH EPSTEIN-                  
C STEP-FUNCTIONS AT THE STEP HEIGHTS XS(M) WITH TRANSITION                      
C THICKNESSES ID(M). RPID(H0,H0,N0,....)=N0.       
C INSTEAD OF 88.0 YOU CAN USE THE HIGHEST ALLOWED ARGUMENT                      
C FOR EXP AT YOUR COMPUTER.                        
C
      REAL 		N0         
      DIMENSION 	ID(4), ST(5), XS(4)                
      COMMON  /ARGEXP/ARGMAX
      SUM=(H-H0)*ST(1)                             
      DO 100  I=1,M   
	      XI=ID(I)                              
 	      A=(H-XS(I))/XI
              AA=0
	      IF(A.GE.ARGMAX) AA=A
	      IF(ABS(A).LT.ARGMAX) AA=ALOG(1.0+EXP(A))   
	      B=(H0-XS(I))/XI
	      BB=0
	      IF(B.GE.ARGMAX) BB=B
	      IF(ABS(B).LT.ARGMAX) BB=ALOG(1.0+EXP(B)) 
100	      SUM=SUM+(ST(I+1)-ST(I))*(AA-BB)*XI 
      SM=0.0
      IF(SUM.GE.ARGMAX) SM=N0*EXP(88.)
      IF(ABS(SUM).LT.ARGMAX) SM=N0*EXP(SUM)
      RPID=SM        
      RETURN          
      END             
C
C
      SUBROUTINE RDHHE (H,HB,RDOH,RDO2H,RNO,PEHE,RDH,RDHE)                      
C BILITZA,FEB.82,H+ AND HE+ RELATIVE PERECENTAGE DENSITY BELOW                  
C 1000 KM. THE O+ AND O2+ REL. PER. DENSITIES SHOULD BE GIVEN                   
C (RDOH,RDO2H). HB IS THE ALTITUDE OF MAXIMAL O+ DENSITY. PEHE                  
C IS THE PRECENTAGE OF HE+ IONS COMPARED TO ALL LIGHT IONS.                     
C RNO IS THE RATIO OF NO+ TO O2+DENSITY AT H=HB.   
      RDHE=0.0        
      RDH=0.0         
      IF(H.LE.HB) GOTO 100                         
      REST=100.0-RDOH-RDO2H-RNO*RDO2H              
      RDH=REST*(1.-PEHE/100.)                      
      RDHE=REST*PEHE/100.                          
100   RETURN          
      END             
C
C
      REAL FUNCTION RDNO(H,HB,RDO2H,RDOH,RNO)      
C D.BILITZA, 1978. NO+ RELATIVE PERCENTAGE DENSITY ABOVE 100KM.                 
C FOR MORE INFORMATION SEE SUBROUTINE RDHHE.       
      IF (H.GT.HB) GOTO 200                        
      RDNO=100.0-RDO2H-RDOH                        
      RETURN          
200   RDNO=RNO*RDO2H  
      RETURN          
      END
C
C
      SUBROUTINE  KOEFP1(PG1O)                     
C THIEMANN,1979,COEFFICIENTS PG1O FOR CALCULATING  O+ PROFILES                  
C BELOW THE F2-MAXIMUM. CHOSEN TO APPROACH DANILOV-                             
C SEMENOV'S COMPILATION.                           
      DIMENSION PG1O(80)                           
      REAL FELD (80)  
      DATA FELD/-11.0,-11.0,4.0,-11.0,0.08018,     
     &0.13027,0.04216,0.25  ,-0.00686,0.00999,     
     &5.113,0.1 ,170.0,180.0,0.1175,0.15,-11.0,    
     &1.0 ,2.0,-11.0,0.069,0.161,0.254,0.18,0.0161,                             
     &0.0216,0.03014,0.1,152.0,167.0,0.04916,      
     &0.17,-11.0,2.0,2.0,-11.0,0.072,0.092,0.014,0.21,                          
     &0.01389,0.03863,0.05762,0.12,165.0,168.0,0.008,                           
     &0.258,-11.0,1.0,3.0,-11.0,0.091,0.088,       
     &0.008,0.34,0.0067,0.0195,0.04,0.1,158.0,172.0,                            
     &0.01,0.24,-11.0,2.0,3.0, -11.0,0.083,0.102,  
     &0.045,0.03,0.00127,0.01,0.05,0.09,167.0,185.0,                            
     &0.015,0.18/     
      K=0             
      DO 10 I=1,80    
      K=K+1           
10    PG1O(K)=FELD(I)                              
      RETURN          
      END             
C
C
      SUBROUTINE KOEFP2(PG2O)                      
C THIEMANN,1979,COEFFICIENTS FOR CALCULATION OF O+ PROFILES                     
C ABOVE THE F2-MAXIMUM (DUMBS,SPENNER:AEROS-COMPILATION)                        
      DIMENSION PG2O(32)                           
      REAL FELD(32)   
      DATA FELD/1.0,-11.0,-11.0,1.0,695.0,-.000781,                             
     &-.00264,2177.0,1.0,-11.0,-11.0,2.0,570.0,    
     &-.002,-.0052,1040.0,2.0,-11.0,-11.0,1.0,695.0,                            
     &-.000786,-.00165,3367.0,2.0,-11.0,-11.0,2.0, 
     &575.0,-.00126,-.00524,1380.0/                
      K=0             
      DO 10 I=1,32    
      K=K+1           
10    PG2O(K)=FELD(I)                              
      RETURN          
      END             
C
C
      SUBROUTINE  KOEFP3(PG3O)                     
C THIEMANN,1979,COEFFICIENTS FOR CALCULATING O2+ PROFILES.                      
C CHOSEN AS TO APPROACH DANILOV-SEMENOV'S COMPILATION.                          
      DIMENSION PG3O(80)                           
      REAL FELD(80)   
      DATA FELD/-11.0,1.0,2.0,-11.0,160.0,31.0,130.0,                           
     &-10.0,198.0,0.0,0.05922,-0.07983,            
     &-0.00397,0.00085,-0.00313,0.0,-11.0,2.0,2.0,-11.0,                        
     &140.0,30.0,130.0,-10.0,                      
     &190.0,0.0,0.05107,-0.07964,0.00097,-0.01118,-0.02614,                     
     &-0.09537,       
     &-11.0,1.0,3.0,-11.0,140.0,37.0,125.0,0.0,182.0,                           
     &0.0,0.0307,-0.04968,-0.00248,                
     &-0.02451,-0.00313,0.0,-11.0,2.0,3.0,-11.0,   
     &140.0,37.0,125.0,0.0,170.0,0.0,              
     &0.02806,-0.04716,0.00066,-0.02763,-0.02247,-0.01919,                      
     &-11.0,-11.0,4.0,-11.0,140.0,45.0,136.0,-9.0, 
     &181.0,-26.0,0.02994,-0.04879,                
     &-0.01396,0.00089,-0.09929,0.05589/           
      K=0             
      DO 10 I=1,80    
      K=K+1           
10    PG3O(K)=FELD(I)                              
      RETURN          
      END             
C
C
      SUBROUTINE SUFE (FIELD,RFE,M,FE)             
C SELECTS THE REQUIRED ION DENSITY PARAMETER SET.
C THE INPUT FIELD INCLUDES DIFFERENT SETS OF DIMENSION M EACH                
C CARACTERISED BY 4 HEADER NUMBERS. RFE(4) SHOULD CONTAIN THE                   
C CHOSEN HEADER NUMBERS.FE(M) IS THE CORRESPONDING SET.                         
      DIMENSION RFE(4),FE(12),FIELD(80),EFE(4)     
      K=0             
100   DO 101 I=1,4    
      K=K+1           
101   EFE(I)=FIELD(K)                              
      DO 111 I=1,M    
      K=K+1           
111   FE(I)=FIELD(K)  
      DO 120 I=1,4    
      IF((EFE(I).GT.-10.0).AND.(RFE(I).NE.EFE(I))) GOTO 100                     
120   CONTINUE        
      RETURN          
      END             
C
C                     
C*************************************************************                  
C************* PEAK VALUES ELECTRON DENSITY ******************                  
C*************************************************************                  
C
C
      SUBROUTINE F2OUT(XMODIP,XLATI,XLONGI,FF0,XM0,UT,
     &                              FOF2,XM3000)
C CALCULATES FOF2/MHZ AND M3000 USING THE CCIR-MAPS.                            
C INPUT: MODIFIED DIP LATITUDE XMODIP, GEOG. LATITUDE XLATI,                    
C LONGITUDE XLONGI (ALL IN DEG.), SMOOTHED SUNSPOT NUMBER R,
C MONTH AND UNIVERSAL TIME UT (DEC. HOURS).                    
C D.BILITZA,JULY 85.  
      DIMENSION FF0(988),XM0(441)                       
      INTEGER QM(7),QF(9)                          
      DATA QF/11,11,8,4,1,0,0,0,0/,QM/6,7,5,2,1,0,0/
      FOF2=GAMMA1(XMODIP,XLATI,XLONGI,UT,6,QF,9,76,13,988,FF0)                  
      XM3000=GAMMA1(XMODIP,XLATI,XLONGI,UT,4,QM,7,49,9,441,XM0)                 
      RETURN          
      END             
C
C
      REAL FUNCTION HMF2ED(XMAGBR,R,X,XM3)         
C CALCULATES THE PEAK HEIGHT HMF2/KM FOR THE MAGNETIC                           
C LATITUDE XMAGBR/DEG. AND THE SMOOTHED ZUERICH SUNSPOT                         
C NUMBER R USING CCIR-M3000 XM3 AND THE RATIO X=FOF2/FOE.                       
C [REF. D.BILITZA ET AL., TELECOMM.J., 46, 549-553, 1979]                       
C D.BILITZA,1980.     
      F1=(2.32E-3)*R+0.222                         
      F2=1.2-(1.16E-2)*EXP((2.39E-2)*R)            
      F3=0.096*(R-25.0)/150.0                      
      DELM=F1*(1.0-R/150.0*EXP(-XMAGBR*XMAGBR/1600.0))/(X-F2)+F3                
      HMF2ED=1490.0/(XM3+DELM)-176.0               
      RETURN          
      END             
C
C
      REAL FUNCTION FOF1ED(YLATI,R,CHI)
C CALCULATES THE F1 PEAK PLASMA FREQUENCY (FOF1/MHZ)
C FOR DIP-LATITUDE (YLATI/DEG.), SMOOTHED 
C ZURICH SUNSPOT NUMBER (R), AND SOLAR ZENITH ANGLE (CHI/DEG.).
C [REF: E.D.DUCHARME ET AL., RADIO SCI., 6,369-378, 1971                        
C AND 8, 837-839, 1973; HOWEVER WITH MAGNETIC INSTEAD OF GEO-
C MAGNETIC LATITUDE,R.EYFRIG,1979]                    
C D. BILITZA, 1988.   
	COMMON/CONST/UMR
	FOF1 = 0.0
	DLA = ABS(YLATI)
		CHI0 = 49.84733 + 0.349504 * DLA
		CHI100 = 38.96113 + 0.509932 * DLA
		CHIM = ( CHI0 + ( CHI100 - CHI0 ) * R / 100. )
		IF(CHI.GT.CHIM) GOTO 1 
  	F0 = 4.35 + DLA * ( 0.0058 - 1.2E-4 * DLA ) 
	F100 = 5.348 + DLA * ( 0.011 - 2.3E-4 * DLA )
	FS = F0 + ( F100 - F0 ) * R / 100.0
	XMUE = 0.093 + DLA * ( 0.0046 - 5.4E-5 * DLA ) + 3.0E-4 * R
	FOF1 = FS * COS( CHI * UMR ) ** XMUE
1	FOF1ED = FOF1     
	RETURN
	END             
C
C
      REAL FUNCTION FOEEDI(COV,XHI,XHIM,SLATI)
C-------------------------------------------------------
C CALCULATES FOE/MHZ BY THE EDINBURGH-METHOD.      
C INPUT: MEAN 10.7CM SOLAR RADIO FLUX (COV), GEOGRAPHIC
C LATITUDE (SLATI/DEG), SOLAR ZENITH ANGLE (XHI/DEG AND 
C XHIM/DEG AT NOON).
C REFERENCE: 
C 	KOURIS-MUGGELETON, CCIR DOC. 6/3/07, 1973
C 	TROST, J. GEOPHYS. RES. 84, 2736, 1979 (was used
C		to improve the nighttime varition)
C D.BILITZA--------------------------------- AUGUST 1986.    
      COMMON/CONST/UMR
C variation with solar activity (factor A) ...............
      A=1.0+0.0094*(COV-66.0)                      
C variation with noon solar zenith angle (B) and with latitude (C)
      XLATI=ABS(SLATI)                             
      SL=COS(XLATI*UMR)
	IF(XLATI.LT.32.0) THEN
		SM=-1.93+1.92*SL                             
		C=23.0+116.0*SL                              
 	ELSE
	 	SM=0.11-0.49*SL                              
	 	C=92.0+35.0*SL  
  	ENDIF
	bb=COS(XHIM*UMR)
	if(bb.le.0.0) then
		if(bb.ne.0.0) then
			bb=-bb
		else
			bb=0.0001
		endif
	endif
	B=bb**SM
C variation with solar zenith angle (D) ..........................        
 	IF(XLATI.GT.12.0) THEN
		SP=1.2
	ELSE
		SP=1.31         
	ENDIF
C adjusted solar zenith angle during nighttime (XHIC) .............
      XHIC=XHI-3.*ALOG(1.+EXP((XHI-89.98)/3.))   
      D=COS(XHIC*UMR)**SP       
C determine foE**4 ................................................
      R4FOE=A*B*C*D     
C minimum allowable foE (sqrt[SMIN])...............................
      SMIN=0.121+0.0015*(COV-60.)
      SMIN=SMIN*SMIN
      IF(R4FOE.LT.SMIN) R4FOE=SMIN                     
      FOEEDI=R4FOE**0.25                           
      RETURN          
      END   
C
C
      REAL FUNCTION XMDED(XHI,R,YW)                
C D. BILITZA, 1978, CALCULATES ELECTRON DENSITY OF D MAXIMUM.                   
C XHI/DEG. IS SOLAR ZENITH ANGLE, R SMOOTHED ZURICH SUNSPOT NUMBER              
C AND YW/M-3 THE ASSUMED CONSTANT NIGHT VALUE.     
C [REF.: D.BILITZA, WORLD DATA CENTER A REPORT UAG-82,7,                        
C       BOULDER,1981]                              
      COMMON/CONST/UMR
      Y=6.05E8+0.088E8*R                           
      Z=(-0.1/(ALOG(YW/Y)))**0.3704                
      SUXHI=ACOS(Z)  
      IF (SUXHI.LT.1.0472) SUXHI=1.0472            
      XXHI=XHI*UMR    
      IF (XXHI.GT.SUXHI) GOTO 100                  
      X=COS(XXHI)     
      XMDED=Y*EXP(-0.1/X**2.7)                     
      RETURN          
100   XMDED=YW        
      RETURN          
      END
C
C
      REAL FUNCTION GAMMA1(SMODIP,SLAT,SLONG,HOUR,IHARM,NQ,
     &  			K1,M,MM,M3,SFE)      
C CALCULATES GAMMA1=FOF2 OR M3000 USING CCIR NUMERICAL MAP                      
C COEFFICIENTS SFE(M3) FOR MODIFIED DIP LATITUDE (SMODIP/DEG)
C GEOGRAPHIC LATITUDE (SLAT/DEG) AND LONGITUDE (SLONG/DEG)  
C AND UNIVERSIAL TIME (HOUR/DECIMAL HOURS).
C NQ(K1) IS AN INTEGER ARRAY GIVING THE HIGHEST DEGREES IN 
C LATITUDE FOR EACH LONGITUDE HARMONIC.                  
C M=1+NQ1+2(NQ2+1)+2(NQ3+1)+... .                  
C SHEIKH,4.3.77.      
      REAL*8 C(12),S(12),COEF(100),SUM             
      DIMENSION NQ(K1),XSINX(13),SFE(M3)           
      COMMON/CONST/UMR
      HOU=(15.0*HOUR-180.0)*UMR                    
      S(1)=SIN(HOU)   
      C(1)=COS(HOU)   
      DO 250 I=2,IHARM                             
      C(I)=C(1)*C(I-1)-S(1)*S(I-1)                 
      S(I)=C(1)*S(I-1)+S(1)*C(I-1)                 
250   CONTINUE        
      DO 300 I=1,M    
      MI=(I-1)*MM     
      COEF(I)=SFE(MI+1)                            
      DO 300 J=1,IHARM                             
      COEF(I)=COEF(I)+SFE(MI+2*J)*S(J)+SFE(MI+2*J+1)*C(J)                       
300   CONTINUE        
      SUM=COEF(1)     
      SS=SIN(SMODIP*UMR)                           
      S3=SS           
      XSINX(1)=1.0    
      INDEX=NQ(1)     
      DO 350 J=1,INDEX                             
      SUM=SUM+COEF(1+J)*SS                         
      XSINX(J+1)=SS   
      SS=SS*S3        
350   CONTINUE        
      XSINX(NQ(1)+2)=SS                            
      NP=NQ(1)+1      
      SS=COS(SLAT*UMR)                             
      S3=SS           
      DO 400 J=2,K1   
      S0=SLONG*(J-1.)*UMR                          
      S1=COS(S0)      
      S2=SIN(S0)      
      INDEX=NQ(J)+1   
      DO 450 L=1,INDEX                             
      NP=NP+1         
      SUM=SUM+COEF(NP)*XSINX(L)*SS*S1              
      NP=NP+1         
      SUM=SUM+COEF(NP)*XSINX(L)*SS*S2              
450   CONTINUE        
      SS=SS*S3        
400   CONTINUE        
      GAMMA1=SUM      
      RETURN          
      END             
C
C                     
C************************************************************                   
C*************** EARTH MAGNETIC FIELD ***********************                   
C**************************************************************                 
C
C
      SUBROUTINE GGM(ART,LONG,LATI,MLONG,MLAT)            
C CALCULATES GEOMAGNETIC LONGITUDE (MLONG) AND LATITUDE (MLAT) 
C FROM GEOGRAFIC LONGITUDE (LONG) AND LATITUDE (LATI) FOR ART=0
C AND REVERSE FOR ART=1. ALL ANGLES IN DEGREE.
C LATITUDE:-90 TO 90. LONGITUDE:0 TO 360 EAST.         
      INTEGER ART     
      REAL MLONG,MLAT,LONG,LATI
      COMMON/CONST/FAKTOR
      ZPI=FAKTOR*360.                              
      CBG=11.4*FAKTOR                              
      CI=COS(CBG)     
      SI=SIN(CBG)
      IF(ART.EQ.0) GOTO 10                         
      CBM=COS(MLAT*FAKTOR)                           
      SBM=SIN(MLAT*FAKTOR)                           
      CLM=COS(MLONG*FAKTOR)                          
      SLM=SIN(MLONG*FAKTOR)
      SBG=SBM*CI-CBM*CLM*SI                     
      LATI=ASIN(SBG)   
      CBG=COS(LATI)     
      SLG=(CBM*SLM)/CBG  
      CLG=(SBM*SI+CBM*CLM*CI)/CBG
        IF(ABS(CLG).GT.1.) CLG=SIGN(1.,CLG)                  
      LONG=ACOS(CLG)  
      IF(SLG.LT.0.0) LONG=ZPI-ACOS(CLG)            
      LATI=LATI/FAKTOR    
      LONG=LONG/FAKTOR  
      LONG=LONG-69.8    
      IF(LONG.LT.0.0) LONG=LONG+360.0                 
      RETURN          
10    YLG=LONG+69.8    
      CBG=COS(LATI*FAKTOR)                           
      SBG=SIN(LATI*FAKTOR)                           
      CLG=COS(YLG*FAKTOR)                          
      SLG=SIN(YLG*FAKTOR)                          
      SBM=SBG*CI+CBG*CLG*SI                        
      MLAT=ASIN(SBM)   
      CBM=COS(MLAT)     
      SLM=(CBG*SLG)/CBM                            
      CLM=(-SBG*SI+CBG*CLG*CI)/CBM
        IF(CLM.GT.1.) CLM=1.                 
      MLONG=ACOS(CLM)  
      IF(SLM.LT..0) MLONG=ZPI-ACOS(CLM)             
      MLAT=MLAT/FAKTOR    
      MLONG=MLONG/FAKTOR  
      RETURN          
      END             
C
C
      SUBROUTINE FIELDG(DLAT,DLONG,ALT,X,Y,Z,F,DIP,DEC,SMODIP)                  
C THIS IS A SPECIAL VERSION OF THE POGO 68/10 MAGNETIC FIELD                    
C LEGENDRE MODEL. TRANSFORMATION COEFF. G(144) VALID FOR 1973.                  
C INPUT: DLAT, DLONG=GEOGRAPHIC COORDINATES/DEG.(-90/90,0/360),                 
C        ALT=ALTITUDE/KM.                          
C OUTPUT: F TOTAL FIELD (GAUSS), Z DOWNWARD VERTICAL COMPONENT                  
C        X,Y COMPONENTS IN THE EQUATORIAL PLANE (X TO ZERO LONGITUDE).          
C        DIP INCLINATION ANGLE(DEGREE). SMODIP RAWER'S MODFIED DIP.             
C SHEIK,1977.         
      DIMENSION H(144),XI(3),G(144),FEL1(72),FEL2(72)
      COMMON/CONST/UMR                           
      DATA FEL1/0.0, 0.1506723,0.0101742, -0.0286519, 0.0092606,                
     & -0.0130846, 0.0089594, -0.0136808,-0.0001508, -0.0093977,                
     & 0.0130650, 0.0020520, -0.0121956, -0.0023451, -0.0208555,                
     & 0.0068416,-0.0142659, -0.0093322, -0.0021364, -0.0078910,                
     & 0.0045586,  0.0128904, -0.0002951, -0.0237245,0.0289493,                 
     & 0.0074605, -0.0105741, -0.0005116, -0.0105732, -0.0058542,               
     &0.0033268, 0.0078164,0.0211234, 0.0099309, 0.0362792,                     
     &-0.0201070,-0.0046350,-0.0058722,0.0011147,-0.0013949,                    
     & -0.0108838,  0.0322263, -0.0147390,  0.0031247, 0.0111986,               
     & -0.0109394,0.0058112,  0.2739046, -0.0155682, -0.0253272,                
     &  0.0163782, 0.0205730,  0.0022081, 0.0112749,-0.0098427,                 
     & 0.0072705, 0.0195189, -0.0081132, -0.0071889, -0.0579970,                
     & -0.0856642, 0.1884260,-0.7391512, 0.1210288, -0.0241888,                 
     & -0.0052464, -0.0096312, -0.0044834, 0.0201764,  0.0258343,               
     &0.0083033,  0.0077187/                       
      DATA FEL2/0.0586055,0.0102236,-0.0396107,    
     & -0.0167860, -0.2019911, -0.5810815,0.0379916,  3.7508268,                
     & 1.8133030, -0.0564250, -0.0557352, 0.1335347, -0.0142641,                
     & -0.1024618,0.0970994, -0.0751830,-0.1274948, 0.0402073,                  
     &  0.0386290, 0.1883088,  0.1838960, -0.7848989,0.7591817,                 
     & -0.9302389,-0.8560960, 0.6633250, -4.6363869, -13.2599277,               
     & 0.1002136,  0.0855714,-0.0991981, -0.0765378,-0.0455264,                 
     &  0.1169326, -0.2604067, 0.1800076, -0.2223685, -0.6347679,               
     &0.5334222, -0.3459502,-0.1573697,  0.8589464, 1.7815990,                  
     &-6.3347645, -3.1513653, -9.9927750,13.3327637, -35.4897308,               
     &37.3466339, -0.5257398,  0.0571474, -0.5421217,  0.2404770,               
     & -0.1747774,-0.3433644, 0.4829708,0.3935944, 0.4885033,                   
     &  0.8488121, -0.7640999, -1.8884945, 3.2930784,-7.3497229,                
     & 0.1672821,-0.2306652, 10.5782146, 12.6031065, 8.6579742,                 
     & 215.5209961, -27.1419220,22.3405762,1108.6394043/                        
      K=0             
      DO 10 I=1,72    
      K=K+1           
      G(K)=FEL1(I)    
10    G(72+K)=FEL2(I)                              
      RLAT=DLAT*UMR   
      CT=SIN(RLAT)    
      ST=COS(RLAT)    
      NMAX=11         
      D=SQRT(40680925.0-272336.0*CT*CT)            
      RLONG=DLONG*UMR                              
      CP=COS(RLONG)   
      SP=SIN(RLONG)   
      ZZZ=(ALT+40408589.0/D)*CT/6371.2             
      RHO=(ALT+40680925.0/D)*ST/6371.2             
      XXX=RHO*CP      
      YYY=RHO*SP      
      RQ=1.0/(XXX*XXX+YYY*YYY+ZZZ*ZZZ)             
      XI(1)=XXX*RQ    
      XI(2)=YYY*RQ    
      XI(3)=ZZZ*RQ    
      IHMAX=NMAX*NMAX+1                            
      LAST=IHMAX+NMAX+NMAX                         
      IMAX=NMAX+NMAX-1                             
      DO 100 I=IHMAX,LAST                          
100   H(I)=G(I)       
      DO 200 K=1,3,2  
      I=IMAX          
      IH=IHMAX        
300   IL=IH-I         
      F1=2./(I-K+2.)  
      X1=XI(1)*F1     
      Y1=XI(2)*F1     
      Z1=XI(3)*(F1+F1)                             
      I=I-2           
      IF((I-1).LT.0) GOTO 400                      
      IF((I-1).EQ.0) GOTO 500                      
      DO 600 M=3,I,2  
      H(IL+M+1)=G(IL+M+1)+Z1*H(IH+M+1)+X1*(H(IH+M+3)-H(IH+M-1))-                
     &Y1*(H(IH+M+2)+H(IH+M-2))                     
      H(IL+M)=G(IL+M)+Z1*H(IH+M)+X1*(H(IH+M+2)-H(IH+M-2))+                      
     &Y1*(H(IH+M+3)+H(IH+M-1))                     
600   CONTINUE        
500   H(IL+2)=G(IL+2)+Z1*H(IH+2)+X1*H(IH+4)-Y1*(H(IH+3)+H(IH))                  
      H(IL+1)=G(IL+1)+Z1*H(IH+1)+Y1*H(IH+4)+X1*(H(IH+3)-H(IH))                  
400   H(IL)=G(IL)+Z1*H(IH)+2.0*(X1*H(IH+1)+Y1*H(IH+2))                          
700   IH=IL           
      IF(I.GE.K) GOTO 300                          
200   CONTINUE        
      S=0.5*H(1)+2.0*(H(2)*XI(3)+H(3)*XI(1)+H(4)*XI(2))                         
      XT=(RQ+RQ)*SQRT(RQ)                          
      X=XT*(H(3)-S*XXX)                            
      Y=XT*(H(4)-S*YYY)                            
      Z=XT*(H(2)-S*ZZZ)                            
      F=SQRT(X*X+Y*Y+Z*Z)                          
      BRH0=Y*SP+X*CP  
      Y=Y*CP-X*SP     
      X=Z*ST-BRH0*CT  
      Z=-Z*CT-BRH0*ST                              
      DIP=ASIN(Z/F)  
      DEC=ASIN(Y/SQRT(X*X+Y*Y))                   
      SMODIP=ASIN(DIP/SQRT(DIP*DIP+ST))           
      DIP=DIP/UMR     
      DEC=DEC/UMR     
      SMODIP=SMODIP/UMR                            
      RETURN          
      END             
C
C
C************************************************************                   
C*********** INTERPOLATION AND REST ***************************                 
C**************************************************************                 
C
C
      SUBROUTINE REGFA2(X11,X22,EPS,FW,F,SCHALT,X) 
C REGULA-FALSI-PROCEDURE TO FIND X WITH F(X)-FW=0. X1,X2 ARE THE                
C STARTING VALUES. THE COMUTATION ENDS WHEN THE X-INTERVAL                      
C HAS BECOME LESS THAN EPS . IF SIGN(F(X1)-FW)= SIGN(F(X2)-FW)                  
C THEN SCHALT=.TRUE.  
      LOGICAL LINKS,SCHALT !,L1                    
      SCHALT=.FALSE.  
      X1=X11          
      X2=X22          
      F1=F(X1)-FW     
      F2=F(X2)-FW     
      IF(F1*F2.LE.0.0) GOTO 200
	X=0.0           
	SCHALT=.TRUE.
	RETURN   
200   X=(X1*F2-X2*F1)/(F2-F1)
      FX=F(X)-FW 
      LINKS=(F1*FX.GT.0.0)  
      IF(LINKS) THEN
	X1=X            
 	F1=FX           
      ELSE
	X2=X            
 	F2=FX
      ENDIF  
700   IF(ABS(X2-X1).GT.EPS) GOTO 200 
      RETURN          
      END             
C
C
      SUBROUTINE REGFA1(X11,X22,EPS,FW,F,SCHALT,X) 
C REGULA-FALSI-PROCEDURE TO FIND X WITH F(X)-FW=0. X1,X2 ARE THE                
C STARTING VALUES. THE COMUTATION ENDS WHEN THE X-INTERVAL                      
C HAS BECOME LESS THAN EPS . IF SIGN(F(X1)-FW)= SIGN(F(X2)-FW)                  
C THEN SCHALT=.TRUE.  
      LOGICAL L1,LINKS,K,SCHALT                    
      SCHALT=.FALSE.
      EP=EPS  
      X1=X11          
      X2=X22          
      F1=F(X1)-FW     
      F2=F(X2)-FW     
      K=.FALSE.       
      NG=2       
      LFD=0     
      IF(F1*F2.LE.0.0) GOTO 200
   	X=0.0           
    	SCHALT=.TRUE.   
      	RETURN
200   X=(X1*F2-X2*F1)/(F2-F1)                      
      GOTO 400        
300   	L1=LINKS        
	DX=(X2-X1)/NG
     	IF(.NOT.LINKS) DX=DX*(NG-1)
     	X=X1+DX
400   FX=F(X)-FW
      LFD=LFD+1
      IF(LFD.GT.20) THEN
	EP=EP*10.
	LFD=0
      ENDIF 
      LINKS=(F1*FX.GT.0.0)
      K=.NOT.K        
      IF(LINKS) THEN
	X1=X            
 	F1=FX           
      ELSE
	X2=X 
	F2=FX 
      ENDIF   
      IF(ABS(X2-X1).LE.EP) GOTO 800               
      IF(K) GOTO 300  
      IF((LINKS.AND.(.NOT.L1)).OR.(.NOT.LINKS.AND.L1)) NG=2*NG                  
      GOTO 200        
800   RETURN          
      END             
C
C
      SUBROUTINE TAL(SHABR,SDELTA,SHBR,SDTDH0,AUS6,SPT)                         
C CALCULATES THE COEFFICIENTS SPT FOR THE POLYNOMIAL
C Y(X)=1+SPT(1)*X**2+SPT(2)*X**3+SPT(3)*X**4+SPT(4)*X**5               
C TO FIT THE VALLEY IN Y, REPRESENTED BY:                
C Y(X=0)=1, THE X VALUE OF THE DEEPEST VALLEY POINT (SHABR),                    
C THE PRECENTAGE DEPTH (SDELTA), THE WIDTH (SHBR) AND THE                       
C DERIVATIVE DY/DX AT THE UPPER VALLEY BOUNDRY (SDTDH0).                        
C IF THERE IS AN UNWANTED ADDITIONAL EXTREMUM IN THE VALLEY                     
C REGION, THEN AUS6=.TRUE., ELSE AUS6=.FALSE..     
C FOR -SDELTA THE COEFF. ARE CALCULATED FOR THE FUNCTION                        
C Y(X)=EXP(SPT(1)*X**2+...+SPT(4)*X**5).           
      DIMENSION SPT(4)                             
      LOGICAL AUS6    
      Z1=-SDELTA/(100.0*SHABR*SHABR)               
      IF(SDELTA.GT.0.) GOTO 500                    
      SDELTA=-SDELTA  
      Z1=ALOG(1.-SDELTA/100.)/(SHABR*SHABR)        
500   Z3=SDTDH0/(2.*SHBR)                          
      Z4=SHABR-SHBR   
      SPT(4)=2.0*(Z1*(SHBR-2.0*SHABR)*SHBR+Z3*Z4*SHABR)/                        
     &  (SHABR*SHBR*Z4*Z4*Z4)                        
      SPT(3)=Z1*(2.0*SHBR-3.0*SHABR)/(SHABR*Z4*Z4)-
     &  (2.*SHABR+SHBR)*SPT(4)          
      SPT(2)=-2.0*Z1/SHABR-2.0*SHABR*SPT(3)-3.0*SHABR*SHABR*SPT(4)              
      SPT(1)=Z1-SHABR*(SPT(2)+SHABR*(SPT(3)+SHABR*SPT(4)))                      
      AUS6=.FALSE.    
      B=4.*SPT(3)/(5.*SPT(4))+SHABR                
      C=-2.*SPT(1)/(5*SPT(4)*SHABR)                
      Z2=B*B/4.-C     
      IF(Z2.LT.0.0) GOTO 300                       
      Z3=SQRT(Z2)     
      Z1=B/2.         
      Z2=-Z1+Z3       
      IF(Z2.GT.0.0.AND.Z2.LT.SHBR) AUS6=.TRUE.     
      IF (ABS(Z3).GT.1.E-15) GOTO 400              
      Z2=C/Z2         
      IF(Z2.GT.0.0.AND.Z2.LT.SHBR) AUS6=.TRUE.     
      RETURN          
400   Z2=-Z1-Z3       
      IF(Z2.GT.0.0.AND.Z2.LT.SHBR) AUS6=.TRUE.     
300   RETURN          
      END             
C
C
C******************************************************************
C********** ZENITH ANGLE, DAY OF YEAR, TIME ***********************
C******************************************************************
C
C
	subroutine soco (lly,ld,t, flat,Elon,reflon,
     &		DECLINATION, ZENITH, SUNRISE, SUNSET)
c--------------------------------------------------------------------
c	s/r to calculate the solar declination, zenith angle, and
c	sunrise & sunset times  - based on Newbern Smith's algorithm
c	[leo mcnamara, 1-sep-86, last modified 16-jun-87]
c
c in:	lly	local time year [last to digits only]
c	ld	local day of year
c	t	local hour (decimal)
c	flat	northern latitude in degrees
c	elon	east longitude in degrees
c	reflon	east longitude to which the hour & minute are referred
c
c out:	declination declination of the sun in degrees
c	zenith	    zenith angle of the sun in degrees
c	sunrise	    local time of sunrise in hours 
c	sunset	    local time of sunset in hours 
c-------------------------------------------------------------------
c
	common/const/	dtr
c amplitudes of Fourier coefficients  --  1955 epoch.................
	data  	p1,p2,p3,p4,p6 /
     & 	0.017203534,0.034407068,0.051610602,0.068814136,0.103221204 /
c
	lh=int(t)
	lm=int((lh-t)*60)
	flm = lm
	if(lly.gt.100) ly=lly-1900
	ly=lly-60
c
c s/r is formulated in terms of WEST longitude.......................
	Eelon = elon
	if(elon.le.0.) Eelon = 360. + elon
	wlon = 360. - Eelon
	Ereflon = reflon
	if(reflon.le.0.) Ereflon = 360. + reflon
	Wreflon = 360. - Ereflon
c
c time of equinox for that year......................................
	if (ly) 7, 7, 9
7	il = int((3.1 - ly)/4.)
	go to 11
9	il = - int((ly + 0.1)/4.)
11	dt = -0.7819 + 0.24225*ly + il
	td = ld + (lh + flm/60. + Wreflon/15.) / 24.
	te = td - dt
c
c declination of the sun..............................................
	dcl = 23.256 * sin(p1*(te-82.242)) + 0.381 * sin(p2*(te-44.855))
     &      + 0.167 * sin(p3*(te-23.355)) - 0.013 * sin(p4*(te+11.97))
     &      + 0.011 * sin(p6*(te-10.41)) + 0.339137
	declination = dcl
	decrad = dcl * dtr
c
c the equation of time................................................
	tf = te - 0.5
	eqt = -7.38*sin(p1*(tf-4.)) - 9.87*sin(p2*(tf+9.))
     &      + 0.27*sin(p3*(tf-53.)) - 0.2*cos(p4*(tf-17.))
	et = eqt * dtr / 4.
c
	fa = flat * dtr
	fb = wlon * dtr
	fc = Wreflon * dtr
	pho = fc - fb + et
	phi = 0.26179939 * ( t - 12.) + pho
c
	dc = decrad
	a = sin(fa) * sin(dc)
	b = cos(fa) * cos(dc)
	cosx = a + b * cos(phi)
	zenith = acos(cosx) / dtr
c
c calculate sunrise and sunset times --  at the ground...........
c see Explanatory Supplement to the Ephemeris (1961) pg 401......
c sunrise at height h metres is at...............................
c	chi(h) = 90.83 + 0.0347 * sqrt(h)........................
c this includes corrections for horizontal refraction and........
c semi-diameter of the solar disk................................
	ch = cos(90.83 * dtr)
	cosphi = (ch -a ) / b
c if abs(secphi) > 1., sun does not rise/set.....................
c allow for sun never setting - high latitude summer.............
	secphi = 999999.
	if(cosphi.ne.0.) secphi = 1./cosphi
	sunset = 99.
	sunrise = 99.
	if(secphi.gt.-1.0.and.secphi.le.0.) return
c allow for sun never rising - high latitude winter..............
	sunset = -99.
	sunrise = -99.
	if(secphi.gt.0.0.and.secphi.lt.1.) return
c
	phi = acos(cosphi)
	et = et / 0.26179939
	srlt = 12. - phi / 0.26179939
	term = (wreflon - wlon)/15. + et
	sunrise = srlt - term
	sunset = 24. - srlt - term
	if(sunrise.lt.0.) sunrise = sunrise + 24.
	if(sunrise.ge.24.) sunrise = sunrise - 24.
	if(sunset.lt.0.) sunset = sunset + 24.
	if(sunset.ge.24.) sunset = sunset - 24.
c
	return
	end
c
C
      FUNCTION HPOL(HOUR,TW,XNW,SA,SU,DSA,DSU)            
C-------------------------------------------------------
C PROCEDURE FOR SMOOTH TIME-INTERPOLATION USING EPSTEIN  
C STEP FUNCTION AT SUNRISE (SA) AND SUNSET (SU). THE 
C STEP-WIDTH FOR SUNRISE IS DSA AND FOR SUNSET DSU.
C TW,NW ARE THE DAY AND NIGHT VALUE OF THE PARAMETER TO 
C BE INTERPOLATED. SA AND SU ARE TIME OF SUNRISE AND 
C SUNSET IN DECIMAL HOURS.
C BILITZA----------------------------------------- 1979.
	IF(ABS(SU).GT.25.) THEN
		IF(SU.GT.0.0) THEN
			HPOL=TW
		ELSE
			HPOL=XNW
		ENDIF
		RETURN
	ENDIF
 	IF(ABS(SU-SA).GE.24.0) THEN 
		IF(SA.LT.1.0) THEN
			HPOL=TW
		ELSE
			HPOL=XNW        
             	ENDIF           
		RETURN  
	ENDIF
      HPOL=XNW+(TW-XNW)*EPST(HOUR,DSA,SA)+
     &	(XNW-TW)*EPST(HOUR,DSU,SU) 
      RETURN          
      END       
C      
C
	SUBROUTINE MODA(IN,MONTH,IDAY,IYEAR,IDOY)
C-------------------------------------------------------------------
C CALCULATES DAY OF YEAR (IDOY) FROM MONTH (MONTH), DAY (IDAY) AND
C YEAR (IYEAR) IF IN=0, OR MONTH (MONTH) AND DAY (IDAY) FROM DAY OF 
C YEAR (IDOY) AND YEAR (IYEAR), IF IN=1. 
C-------------------------------------------------------------------
	DIMENSION	MO(12),ML(12)
	DATA		MO/0,31,59,90,120,151,181,212,243,273,304,334/
	DATA		ML/0,31,60,91,121,152,182,213,244,274,305,335/
	IMO=0
	MOBE=0
	IF(IYEAR.GT.100) IYEAR=IYEAR-1900
	IF(IN.GT.0) GOTO 5
		IF(MOD(IYEAR,4).EQ.0) IDOY=ML(MONTH)+IDAY
		IF(MOD(IYEAR,4).NE.0) IDOY=MO(MONTH)+IDAY
	RETURN
5		IMO=IMO+1
		MOOLD=MOBE
		IF(IMO.GT.12) GOTO 55
		IF(MOD(IYEAR,4).EQ.0) MOBE=ML(IMO)
		IF(MOD(IYEAR,4).NE.0) MOBE=MO(IMO)
		IF(MOBE.LT.IDOY) GOTO 5
55		MONTH=IMO-1
		IDAY=IDOY-MOOLD
	RETURN
	END		
C
	REAL FUNCTION EPST ( X, B1, X1 )
C -------------------------------------------------------------- STEP
	COMMON/ARGEXP/ARGMAX
	D1 = ( X - X1 ) / B1
	IF (ABS(D1).LT.ARGMAX) GOTO 1
	EPST = 0.
	IF (D1.GT.0.0) EPST = 1.
	RETURN
1	EPST = 1. / ( 1. + EXP( -D1 ))
	RETURN
	END
C
C
	REAL FUNCTION EPSTEP ( C, A, B1, X1, X)
C------------------------------------------------ STEP FROM A TO C	
	EPSTEP = A + ( C - A ) * EPST ( X, B1, X1)
	RETURN
	END
C*****************************************************************
C********* INTERNATIONAL REFERENCE IONOSPHERE (IRI). *************
C*****************************************************************
C**************** IRIMP9, MARCH  1987 ****************************
C****************   DRIVER PROGRAM    ****************************
C*****************************************************************
C*** THIS PROGRAM PRODUCES PROFILES OF                         ***
C***      ELECTRON DENSITY                                     ***
C***      NORMALISED ELECTRON DENSITY (TO F-MAXIMUM)           ***
C***      NEUTRAL TEMPERATURE (CIRA 72)                        ***
C***      ELECTRON TEMPERATURE                                 ***
C***      ION TEMPERATURE                                      ***
C***      ELECTRON TO ION TEMPERATURE RATIO                    ***
C***      RELATIVE PERCENTAGE DENSITIES OF THE IONS            ***
C***           ATOMIC OXYGEN, HYDROGEN, HELIUM,                ***
C***           MOLECULAR OXYGEN AND NITROGEN OXYD (NO+)        ***
C*****************************************************************
C*** FOR SPECIFIED                                             ***
C***      GEOG. OR GEOM. LATITUDE AND LONGITUDE                ***
C***      ALTITUDE                                             ***
C***      ZURICH SUNSPOT NUMBER (12-MONTH-RUNNING MEAN)        ***
C***      MONTH-NUMBER                                         ***
C***      LOCAL OR UNIVERSAL TIME                              ***
C***      YOUR CHOICE OF VARIABLE                              ***
C*****************************************************************
C*** THE ALTITUDE LIMITS ARE:     LOWER           UPPER        ***
C***     ELECTRON DENSITY      60(DAY)/80 KM     1000 KM       ***
C***     TEMPERATURES             120 KM         3000 KM       ***
C***     ION DENSITIES            100 KM         1000 KM       ***
C*****************************************************************
C*     --------------------ADDRESSES------------------------     *
C*     I  PROF. K. RAWER              DR. D. BILITZA       I     *
C*     I  HERRENSTR. 43               GSFC/NSSDC CODE 633  I     *
C*     I  D-7801 MARCH                GREENBELT MD 20771   I     *
C*     I  F.R.G.                      USA                  I     *
C*     -----------------------------------------------------     *
C*****************************************************************
C*****************************************************************
C*****************************************************************
C*********       ALL ANGLES ARE IN DEGREE           **************
C*********       ALL DENSITIES ARE IN M-3           **************
C*********       ALL ALTITUDES ARE IN KM            **************
C*********     ALL TEMPERATURES ARE IN KELVIN       **************
C*********     ALL TIMES ARE IN DECIMAL HOURS       **************
C*****************************************************************
C********************  OPTIONS  **********************************
C*****************************************************************
C* FOR HMF2=0 OR FOF2=0 THE F2 PEAK VALUES ARE CALCULATED WITH   *
C* THE CCIR MODELS. THE CCIR COEFFICIENT SET FOR THE MONTH "mm"  *
C* IS EXPECTED IN THE BINARY FILE "CCIRmm.BIN". IF YOU USE THE   *
C* ASCII CODED FILES "CCIRmm.ASC", YOU HAVE TO INCORPORATE THE   *
C* CHANGES INDICTED IN THE PART: "READ CCIR COEFFICIENT SET FOR  *
C* CHOSEN MONTH."                                                * 
C*****************************************************************
C* HOUR       IS   LOCAL TIME   OR   UNIVERSAL TIME + 25         *
C* JF(11)     ARE FLAGS TO SELECT THE IRI PARAMETERS TO BE CALC.:*
C*            NE,NE/NMF2,TN,TE,TI,TE/TI,O+,H+,HE+,O2+,NO+        *
C* JMAG       =0 GEOGRAFIC LATITUDE AND LONGITUDE EXPECTED       *   
C*            =1 GEOMAGNETIC LATITUDE AND LONGITUDE EXPECTED     *
C*****************************************************************
C* NE(...)    =0 AE-C AND AEROS MODELS ARE USED FOR ELECTRON     *
C*            TEMPERATURE OTHERWISE THE ANTI-CORRELATION-MODEL   *
C*            IS USED WITH THE GIVEN ELECTRON DENSITY VALUES     *
C*****************************************************************
C*****************************************************************
C***   THIS VERSION DIFFERS FROM IRIMP8 BY CONTAINING:         ***
C***    - RESTRICTION OF RZ12/COVM TO VALUES BELOW 145/185     ***
C***    - CORRECTED HARMONIZED BENT                            ***
C***    - NEW ELECTRON TEMPERATURE MODEL                       ***
C***    - NEWER FOE NIGHT-TIME                                 ***
C*****************************************************************
C*****************************************************************

*
*  this is hacked to create a subroutine to give O^+ at a given altitude
*
      subroutine iri(altIn, jMagIn, LatIn, LongIn, rIn, monthIn
     &	  , hourIn, Oplus, TiOut, TeOut, magbr)
	real altIn, LatIn, LongIn, rIn, hourIn, Oplus, TiOut, TeOut, Ne
	integer jMagIn, month
	logical first

      INTEGER 		EGNR,AGNR,DAYNR,DDO,DO2,SEASON
      REAL 		LATI,LONGI,MO2,MO,MODIP,NMF2,MAGBR,
     &  		NMF1,NME,NMD,NEI,MM,MLAT,MLONG,NOBO2
      CHARACTER*48	FILENAME
      CHARACTER*5	ITEXT(6)
      CHARACTER*4	LTEX,IMZ(6)
      DIMENSION F(3),B0F(32),OUTF(11),RIF(4),XSM(4),HOA(3),JF(11),
     & DDO(4),DO2(2),MM(5),SIPL(2),XNAR(3),ATE(7),XM0(441),
     & PF1O(12),PF2O(4),PF3O(12),HO(4),MO2(3),MO(5),E(4),SIPH(2),
     & TEA(6),HO2(2),PG1O(80),PG2O(32),PG3O(80),CTN(3),
     & XDELS(4),DNDS(4),DTI(4),DTE(5),AHH(7),STTE(6),PARMAX(6),
     & FF0(988),F2(13,76,2),FM3(9,49,2),XVAR(6),JFM(11),PARMIN(6) !CTNN(3)
      LOGICAL		EXT,SCHALT,NIGHT,F1REG,UNTI,NOTBEG,CCIRNO
      COMMON	/BLOCK1/HMF2,NMF2,HMF1
     &		/BLOCK2/B0,B1,C1,HZ,T,HST,STR
     &  /BLOCK3/HDX,HME,NME,HMD,NMD,HEF,D1,XKK,FP30,FP3U,FP1,FP2
     &  /BLOCK5/ZX,TNX,ATN,CTN  	/BLOCK6/NIGHT,E
     &  /BLOTE/AHH,ATE1,STTE,DTE	/CONST/UMR
     &  /BLOCK8/HS,TNHS,XSM,MM,DTI,MXSM
     &	/BLO10/BETA,ETA,DELTA,ZETA	/ARGEXP/ARGMAX
      EXTERNAL 		XE1,XE2,XE3,XE4,XE5,XE6,TEDER
      DATA B0F/114.,64.,134.,77.,128.,66.,75.,73.,113.,115.,
     &  150.,116.,138.,123.,94.,132.,72.,84.,83.,89.,75.,85.,
     &  57.,76.,102.,100.,120.,110.,107.,103.,76.,86./,
     &  HOA/300.,400.,600./,XDELS/3*5.,10./,IMZ/'GEOD','GEOD',
     &  2*'    ','L.T.','    '/,DNDS/.016,.01,2*.016/,DDO/9,5,5,50/,
     &  DO2/5,5/,ITEXT/' LATI',' LONG',' RZ12','MONTH',' HOUR',
     &  ' H/KM'/,LATI,LONGI,HEIGHT,R,MONTH,HOUR,IVAR,BVAR,EVAR,
     &  SVAR/45.1,293.1,100,100,10,12.5,6,100,1000,100/,
     &  JF/11*1/,JFM/2,10*1/,XNAR/3*0.0/,JAGNR,JAGNRM,ISELEC/0,1,0/
      DATA PARMIN/-90.,-360.,0.,1.,0.,60./,PARMAX/90.,360.,250.,
     &  12.,49.,3000./
C
C PROGAM CONSTANTS
C

        first = .true.

	ARGMAX=88.0
     	UMR=ATAN(1.0)*4./180.
      	ALOG2=ALOG(2.)
      	MIN0=0
     	MAX1=1
	XNMAX=1.E15
	XMIN0=0.
C
C Insert new code to aleviate block data problem.
C
        HMF2=0
        FOF2=0
        AHH(1)=120.
        AHH(2)=0.
        AHH(3)=300.
        AHH(4)=400.
        AHH(5)=600.
        AHH(6)=1400.
        AHH(7)=3000.
        DTE(1)=5.
        DTE(2)=5.
        DTE(3)=10.
        DTE(4)=20.
        DTE(5)=20.
        DTI(1)=10.
        DTI(2)=10.
        DTI(3)=20.
        DTI(4)=20.
C
C FIRST SPECIFY YOUR COMPUTERS CHANNEL NUMBERS ....................
C EGNR=INPUT, MONITO=MONITOR, KONSOL=MESSAGES......................
C AGNR=OUTPUT (OUTPUT IS DISPLAYED OR STORED IN FILE OUTPUT.IRI)...
C IUCCIR=UNIT NUMBER FOR CCIR COEFFICIENTS ........................
C IUOUT=UNIT NUMBER FOR OUTPUT FILE OUTPUT.IRI ....................
C IUMESS=UNIT NUMBER FOR MESSAGE OUTPUT FILE KONSOL.IRI ...........
C
      LTEX='(LT)'
      ISTART=1
      EGNR=5
      MONITO=6
      IUCCIR=10
      IUOUT=16
      IUMESS=15
*--------------------------------------------------------------------
      NOTBEG=.FALSE.
      MONTHO=0
      RGO=0
      GOTO 5508
3293  DLATI=LATI
      DLONG=LONGI
      IF(JMAG.EQ.1) THEN
	DLATI=MLAT
	DLONG=MLONG
      ENDIF
      DHOU=HOUR
      LTEX='(LT)'
      IF(UNTI) THEN
	DHOU=UT
	LTEX='(UT)'
      ENDIF
      WRITE(MONITO,5602) DLATI,JMAG,DLONG,ITEXT(IVAR),HEIGHT,BVAR,
     &  EVAR,SVAR,R,MONTH,LTEX,DHOU
      MAXI=13
      ISTART=ISTART+1
5602  FORMAT(1X//' **** WHICH PARAMETER DO YOU WANT TO ',
     & 'CHANGE?'/1X,60('-')/' 0  NO FURTHER CHANGES, CALCULATE',
     & ' PROFILE'/' 1  LATITUDE  #',F5.1,'#',8X,'7  GEOD.(0) OR',
     & ' GEOM.(1)    #',I1,'#'/' 2  LONGITUDE #',F5.1,'#',8X,
     & '8  SELECTION OF VARIABLE   #',A5,'#'/' 3  ALTITUDE  #',
     & F6.1,'#',7X,'9  VARIABLE RANGE  #',3F8.1,'#'/' 4  RZ12',
     & 6X,'#',F5.1,'#',7X,'10  SELECTION OF IRI PARAMETERS'/
     & ' 5  MONTH     #',I2,'#',10X,'11  F PEAK VALUES'/' 6  HOUR ',
     & A4,' #',F5.1,'#',7X,'12  TE(NE) MODEL'/28X,'13  DISPLAY',
     & ' OF PROFILE AND MESSAGES'/1X,70('-')/' ENTER NUMBER')   
8102  READ(EGNR,*,ERR=8100,END=3330) ISWIT
      IF((ISWIT.LE.MAXI).AND.(ISWIT.GE.MIN0)) GOTO 8104
8100  	WRITE(MONITO,8101) MIN0,MAXI
8101  	FORMAT(1X,'Your input value is outside the value range:',
     &		I2,' to ',I2,', try again')
      	GOTO 8102
8104  GOTO (5300,3329,3339,3331,3500,5501,5502,3332,5503,5504,
     &    5505,5506,5507,5508) ISWIT+1
*--------------------------------------------------------------------
5508  continue
      MAXI=2
8105  jagnr = 0
      IF((JAGNR.LE.MAXI).AND.(JAGNR.GE.MIN0)) GOTO 8107
8106  	WRITE(MONITO,*) ' error condition'
      	GOTO 8105
8107  IVARNR=0
      AGNR=MONITO
      IF(JAGNR.GT.0) THEN
        OPEN(UNIT=IUOUT,FILE='OUTPUT.IRI',STATUS='NEW',FORM='FORMATTED')
        IF(JAGNR.EQ.1) AGNR=IUOUT
      ENDIF
*--------------------------------------------------------------------
      MAXI=1
      IF(JAGNR.GT.0) MAXI=2
8111  jagnrm = 0
      IF((JAGNRM.GE.MIN0).AND.(JAGNRM.LE.MAXI)) GOTO 8109
8110  	WRITE(MONITO,*) ' error condition'
      	GOTO 8111
8109  KONSOL=MONITO
      IF(JAGNRM.EQ.1) THEN
        OPEN(UNIT=IUMESS,FILE='KONSOL.IRI',STATUS='NEW',
     &		FORM='FORMATTED')
        KONSOL=IUMESS
      	ENDIF
      IF(JAGNRM.EQ.2) KONSOL=IUOUT
      IF(NOTBEG) GOTO 3293
*--------------------------------------------------------------------
5503  continue
      MAXI=6
      MINI=1
8113  ivar = 6
      IF((IVAR.LE.MAXI).AND.(IVAR.GE.MINI)) GOTO 5504
8112  	WRITE(MONITO,*) ' error condition'
      	GOTO 8113
*--------------------------------------------------------------------
5504  continue
      XMAX=PARMAX(IVAR)
      XMIN=PARMIN(IVAR)
8114  bvar = altIn
      evar = altIn
      svar = 100.0
      IF((EVAR.LE.XMAX).AND.(BVAR.GE.XMIN)) GOTO 8116
8115  	WRITE(MONITO,*) ' error condition'
     	GOTO 8114
8116  IF(IVAR.NE.5) GOTO 2929
	IMZ(5)='L.T.'
	UNTI=.FALSE.
      IF(BVAR.GT.24.1) THEN
	UNTI=.TRUE.
	IMZ(5)='U.T.'
	BVAR=BVAR-25.
      ENDIF
2929  LANZ=INT((EVAR-BVAR)/SVAR)+1
      IF(NOTBEG) GOTO 3293
*--------------------------------------------------------------------
3332  continue
8118  jmag = jMagIn
      IF((JMAG.EQ.MAX1).OR.(JMAG.EQ.MIN0)) GOTO 8120
8119  	WRITE(MONITO,*) ' error condition'
      	GOTO 8118
8120  IMZ(1)='GEOD'
      IMZ(2)='GEOD'
      IF(JMAG.EQ.1) THEN
	IMZ(1)='GEOM'
 	IMZ(2)='GEOM'
      	ENDIF
      IF(NOTBEG) GOTO 3293
      IVARNR=IVARNR+1
      IF(IVARNR.EQ.IVAR) GOTO 7339
3329  DLATI=LATI
      IF(JMAG.EQ.1) DLATI=MLAT
*--------------------------------------------------------------------
      XMAX=PARMAX(1)
      XMIN=PARMIN(1)
8121  dlati = LatIn
      IF((DLATI.LE.XMAX).AND.(DLATI.GE.XMIN)) GOTO 8123
8122  	WRITE(MONITO,*) ' error condition'
      	GOTO 8121
8123  IF(JMAG.EQ.0) THEN
	LATI=DLATI
      ELSE
	MLAT=DLATI
      ENDIF
      IF(NOTBEG) GOTO 3293
7339  IVARNR=IVARNR+1
      IF(IVARNR.EQ.IVAR) GOTO 7500
3339  DLONG=LONGI
      IF(JMAG.EQ.1) DLONG=MLONG
*--------------------------------------------------------------------
      XMAX=PARMAX(2)
      XMIN=PARMIN(2)
8124  dlongi = longIn
      IF((DLONGI.LE.XMAX).AND.(DLONGI.GE.XMIN)) GOTO 8126
8125  	WRITE(MONITO,*) ' error condition'
      	GOTO 8124
8126  IF(JMAG.EQ.0) THEN
	LONGI=DLONGI
      ELSE
	MLONG=DLONGI
      ENDIF
      IF(NOTBEG) GOTO 3293
7500  IVARNR=IVARNR+1
      IF(IVARNR.EQ.IVAR) GOTO 7701
*--------------------------------------------------------------------
3500  continue
      XMAX=PARMAX(3)
      XMIN=PARMIN(3)
8127  r = rIn
      IF((R.LE.XMAX).AND.(R.GE.XMIN)) GOTO 8129
8128  	WRITE(MONITO,*) ' error condition'
      	GOTO 8127
8129  IF(NOTBEG) GOTO 3293
7701  IVARNR=IVARNR+1
      IF(IVARNR.EQ.IVAR) GOTO 7502
*--------------------------------------------------------------------
5501  continue
      MAXI=12
      MINI=1
8130  month = monthIn
      IF((MONTH.LT.13).AND.(MONTH.GT.0)) GOTO 8132
8131  	WRITE(MONITO,*) ' error condition'
      	GOTO 8130
8132  IF(NOTBEG) GOTO 3293
7502  IVARNR=IVARNR+1
      IF(IVARNR.EQ.IVAR) GOTO 7331
5502  DHOU=HOUR
      IF(UNTI) DHOU=UT
*--------------------------------------------------------------------
      XMAX=PARMAX(5)
      XMIN=PARMIN(5)
8133  dhour = hourIn				! WARNING, if UT add 25.0
      IF((DHOUR.LE.XMAX).AND.(DHOUR.GE.XMIN)) GOTO 8135
8134  	WRITE(MONITO,*) ' error condition'
      	GOTO 8133
8135  UNTI=.FALSE.
      IF(DHOUR.GT.24.1) THEN
	UNTI=.TRUE.
	UT=DHOUR-25.
      ELSE
	HOUR=DHOUR
      ENDIF
      IF(NOTBEG) GOTO 3293
7331  IVARNR=IVARNR+1
      IF(IVARNR.EQ.IVAR) GOTO 5505
3331  WRITE(MONITO,6002) HEIGHT
6002  FORMAT(1X//1X,'ALTITUDE ?    [KM]',34X,'#',F6.1,'#'//
     & ' !! ALTITUDE=0   profiles of F2,F1,E, and D peak altitudes'/
     & ' !!',14X,'and densities are calculated')
      WRITE(MONITO,8640)
8640  FORMAT(1X,60('-')/' Enter /, to continue with current value(s);',
     &  ' Ctrl Z, to exit')
      XMAX=PARMAX(6)
8136  READ(EGNR,*,ERR=8137,END=3330) HEIGHT
      IF((HEIGHT.LE.XMAX).AND.(HEIGHT.GE.XMIN0)) GOTO 8138
8137  	WRITE(MONITO,8117) XMIN0,XMAX
8117  	FORMAT(1X,'Your input value is outside the value range:',
     &		F7.1,' to ',F7.1,', try again')
      	GOTO 8136
8138  IF(NOTBEG) GOTO 3293
*--------------------------------------------------------------------
5505  continue
8139  iselec = 1
      IF((ISELEC.EQ.MAX1).OR.(ISELEC.EQ.MIN0)) GOTO 8141
8140  	WRITE(MONITO,*) ' error condition'
      	GOTO 8139
8141  IF(ISELEC.EQ.0) GOTO 5533
*--------------------------------------------------------------------
      do jj = 1, 11
	  jf(jj) = 0
      end do
      jf(1) = 1
      jf(4) = 1
      jf(5) = 1
      jf(7) = 1
5533  INECH=JF(1)+JF(2)
      KOMB=JF(3)+JF(4)+JF(5)+JF(6)
      IOND=JF(7)+JF(8)+JF(9)+JF(10)+JF(11)
      IF(NOTBEG) GOTO 3293
5506  IF(INECH.LT.1) GOTO 5507
*--------------------------------------------------------------------
8145  fof2 = 0.0
      IF((FOF2.LE.XNMAX).AND.(FOF2.GE.XMIN0)) GOTO 8147
8146  	WRITE(MONITO,*) ' error condition'
8148  	FORMAT(1X,'Your input value is outside the value range:',
     &		F3.1,' to ',E9.1,', try again')
      	GOTO 8145
*--------------------------------------------------------------------
8147  continue
      XMAX=700.
8149  hmf2 = 0.0
      IF((HMF2.LE.XMAX).AND.(HMF2.GE.XMIN0)) GOTO 8151
8150  	WRITE(MONITO,*) ' error condition'
      	GOTO 8149
8151  CCIRNO=.TRUE.
      IF((HMF2.LT.1.0).OR.(FOF2.LT.1.0)) CCIRNO=.FALSE.
      IF(NOTBEG) GOTO 3293
5507  IF(KOMB.LT.1) GOTO 5300
*--------------------------------------------------------------------
8152  continue
      xnar(1) = 0.
      xnar(2) = 0.
      xnar(3) = 0.
      DO 8154 JXNAR=1,3
  	DLA=XNAR(JXNAR)
8154 	IF((DLA.GT.XNMAX).OR.(DLA.LT.XMIN0)) GOTO 8153
      GOTO 8155
8153  WRITE(MONITO,*) ' error condition'
      GOTO 8152
8155  IF(NOTBEG) GOTO 3293
5300  NOTBEG=.TRUE.
*      IF(HEIGHT.LT.1.0) THEN
*        WRITE(AGNR,8192) ITEXT(IVAR),IMZ(IVAR)
*      ELSE
*        WRITE(AGNR,8193) ITEXT(IVAR),IMZ(IVAR)
*      ENDIF
      IF(JAGNR.LT.2) GOTO 8444
	IF(AGNR.EQ.MONITO) THEN
		AGNR=IUOUT
		GOTO 5300
	ELSE
		AGNR=MONITO
	ENDIF
8192  FORMAT(//////////3X,A5,'  PEAK ALTITUDES IN KM',6X,'PEAK DENSI',
     &  'TIES IN CM-3'/4X,A4,'   HMF2 HMF1  HME  HMD     NMF2   ',
     &  'NMF1    NME     NMD')
8193  FORMAT(//////////3X,A5,'  ELECTRON DENSITY      TEMPERATURES',
     &  7X,'ION PERCENTAGE DENSITIES'/4X,A4,'  NE/CM-3 NE/NMF2  ',
     &  'TN/K  TI/K  TE/K  TE/TI   O+   H+  He+  O2+  NO+')
8444    XVAR(1)=LATI
        XVAR(2)=LONGI
        IF(JMAG.EQ.1) THEN
	  XVAR(1)=MLAT
	  XVAR(2)=MLONG
        ENDIF
        XVAR(3)=R
        XVAR(4)=MONTH
        XVAR(5)=HOUR
        IF(UNTI) XVAR(5)=UT
        XVAR(6)=HEIGHT
        LFD=0
        XVAR(IVAR)=BVAR-SVAR
2123    XVAR(IVAR)=XVAR(IVAR)+SVAR
        LFD=LFD+1
        IF(JMAG.EQ.1) THEN
           MLAT=XVAR(1)
           MLONG=XVAR(2)
        ELSE
           LATI=XVAR(1)
           LONGI=XVAR(2)
        ENDIF
        R=XVAR(3)
        IF(XVAR(4).LT.13.) MONTH=XVAR(4)
        IF(UNTI) THEN
	  UT=XVAR(5)
	ELSE
          HOUR=XVAR(5)
        ENDIF
        HEIGHT=XVAR(6)
C
C CALCULATION OF MEAN F10.7CM SOLAR RADIO FLUX (COV)................
C CALCULATION OF RESTRICTED SOLAR ACTIVITIES (RG,COVG)..............
C
      COV=63.75+R*(0.728+R*0.00089)
      RG=R
      COVG=COV
      IF(R.GT.150.) RG=150.
      IF(COV.GT.193.) COVG=193.
C
C CALCULATION OF GEOG. OR GEOM. COORDINATES IN DEG....................
C CALCULATION OF MAGNETIC INCLINATION (DIP), DECLINATION (DEC)........
C   DIP LATITUDE (MAGBR) AND MODIFIED DIP (MODIP). ALL IN DEGREE......
C
      IF((LFD.EQ.1).OR.(IVAR.LT.3)) THEN
        CALL GGM(JMAG,LONGI,LATI,MLONG,MLAT)
        ABSLAT=ABS(LATI)
        CALL FIELDG(LATI,LONGI,300.0,XMA,YMA,ZMA,BET,DIP,DEC,MODIP)
        MAGBR=ATAN(0.5*TAN(DIP*UMR))/UMR
      ENDIF     
C
C CALCULATION OF SEASON (SUMMER=2, WINTER=4)..........................
C CALCULATION OF DAY OF YEAR AND SUN DECLINATION......................
C
      DAYNR=INT(MONTH*30.42-15.21)
      SUNDEC=-0.40915*COS(.0172142063*(DAYNR+8.))
      SEASON=INT((DAYNR+45.0)/92.0)
      IF(SEASON.LT.1) SEASON=4
      NSEASON=SEASON
      IF(LATI.GT.0.0) GOTO 5592
   	SEASON=SEASON-2
    	IF(SEASON.LT.1) SEASON=SEASON+4
C
C CALCULATION OF SOLAR ZENITH ANGLE (XHI/DEG).........................
C NOON VALUE (XHINON).................................................
C
5592  IF(UNTI) THEN
	HOUR=UT+LONGI/15.
	IF(HOUR.GT.24.) HOUR=HOUR-24.
      ELSE
        UT=HOUR-LONGI/15.
        IF(UT.LT.0.) UT=UT+24.
      ENDIF
2918  Z1=SIN(SUNDEC)*SIN(LATI*UMR)
      Z2=COS(SUNDEC)*COS(LATI*UMR)
      IF(ABS(Z2).GT.0.0) GOTO 120
      SAX=24.0
      IF(Z1.GE.0.0) SAX=0.0
      GOTO 140
120   IF(ABS(Z1/Z2).LE.1.0) GOTO 510
      SAX=24.0
      IF((Z1+Z2).GE.0.0) SAX=0.0
      GOTO 140
510   SAX=12.0-ACOS(-Z1/Z2)/(UMR*15.0)
140   SUX=24.0-SAX
      XLSTA=15.0*(HOUR-12.0)
      COSXHI=Z1+Z2*COS(XLSTA*UMR)
      XHI=ACOS(COSXHI)/UMR
      XHINON=ACOS(Z1+Z2)/UMR
      IF(XHINON.LE.89.0) GOTO 503
      WRITE(KONSOL,5569) XHINON
5569  FORMAT(1X,'!!!!! YOUR NOON ZENITH ANGLE IS',F7.1,
     &  '; RE-ADJUSTED TO 89.0 !!!!!')
      XHINON=89.0
503   NIGHT=.FALSE.
      IF((HOUR.GT.SUX).OR.(HOUR.LT.SAX)) NIGHT=.TRUE.
      SUNDEC=SUNDEC/UMR
C
C CALCULATION OF ELECTRON DENSITY PARAMETERS................
C
      HNEA=65.
      IF(NIGHT) HNEA=80.
      HNEE=2000.
      IF(INECH.LT.1) GOTO 4933
      DELA=4.32
      IF(ABS(MODIP).GE.18.) DELA=1.0+EXP(-(ABS(MODIP)-30.0)/10.0)
      DELL=1+EXP(-(ABSLAT-20.)/10.)
C!!!!!!! F-REGION PARAMETERS AND E-PEAK !!!!!!!!!!!!!!!!!!!!!!!!!!
      FOE=FOEEDI(COV,XHI,XHINON,LATI)
      IF(CCIRNO) GOTO 501
      IF((MONTH.EQ.MONTHO).AND.(RG.EQ.RGO)) GOTO 4292
      IF(MONTH.EQ.MONTHO) GOTO 4291
C
C READ CCIR COEFFICIENT SET FOR CHOSEN MONTH....................
C
*
* WARNING, this will have to edited for a different computer or directory
*
        WRITE(FILENAME,104) MONTH+10
104     FORMAT('[craig.empirical.iri]CCIR',I2,'.BIN')
        OPEN(IUCCIR,FILE=FILENAME,STATUS='OLD',ERR=8448,
     &		FORM='UNFORMATTED')
        READ(IUCCIR) F2,FM3
C !!!!!!!!!!!!!!! FOR ASCII CODED CCIR COEFFIECENTS FILES
C !!!!!!!!!!!!!!! SUBSTITUTE THE LAST 3 STATEMENTS BY:
C104     FORMAT('CCIR',I2,'.ASC')
C        OPEN(IUCCIR,FILE=FILENAME,STATUS='OLD',ERR=8448,
C     &		FORM='FORMATTED')
C        READ(IUCCIR,4689) F2,FM3
C4689    FORMAT(1X,4E15.8)
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CLOSE(IUCCIR)
	GOTO 4291
8448	WRITE(MONITO,8449) FILENAME
8449	FORMAT(1X////,
     &	  ' !!!! The file ',A10,' is not in your directory,'//
     &	  ' !!!!   try a different MONTH (enter: 5),'/
     &    ' !!!!   or exit and load file (enter: Ctrl Z)')
	GOTO 3293
C
C LINEAR INTERPOLATION IN SOLAR ACTIVITY
C
4291    RR2=RG/100.
        RR1=1.-RR2
        DO 20 I=1,76
        DO 20 J=1,13
        K=J+13*(I-1)
20      FF0(K)=F2(J,I,1)*RR1+F2(J,I,2)*RR2
        DO 30 I=1,49
        DO 30 J=1,9
        K=J+9*(I-1)
30      XM0(K)=FM3(J,I,1)*RR1+FM3(J,I,2)*RR2
4292  CALL F2OUT(MODIP,LATI,LONGI,FF0,XM0,UT,FOF2,XM3000)
      HMF2=HMF2ED(MAGBR,RG,FOF2/FOE,XM3000)
501   IF(FOF2.GT.100.0) FOF2=SQRT(FOF2/1.24E10)
      NMF2=1.24E10*FOF2*FOF2
      COS2=COS(MLAT*UMR)
      COS2=COS2*COS2
      FLU=(COVG-40.0)/30.0
      ETA1=-0.0070305*COS2
      IF(JF(1).EQ.2) GOTO 5566
        EX=EXP(-MLAT/15.)
        EX1=EX+1
        EPIN=4.*EX/(EX1*EX1)
        ETA1=-0.02*EPIN
5566  ETA=0.058798+ETA1+FLU*(-0.014065+0.0069724*COS2)+
     &(0.0024287+0.0042810*COS2-0.00015280*FOF2)*FOF2
      ZETA=0.078922-0.0046702*COS2+FLU*(-0.019132+0.0076545*COS2)+
     &(0.0032513+0.0060290*COS2-0.00020872*FOF2)*FOF2
      BETA=-128.03+20.253*COS2+FLU*(-8.0755-0.65896*COS2)+(0.44041
     &+0.71458*COS2-0.042966*FOF2)*FOF2
      Z=EXP(94.45/BETA)
      Z1=Z+1
      Z2=Z/(BETA*Z1*Z1)
      DELTA=(ETA/Z1-ZETA/2.0)/(ETA*Z2+ZETA/400.0)
1501  B1=3.0
C!!!!!!! INTERPOLATION FOR B0 OUT OF ARRAY B0F !!!!!!!!!!!!!!!!!!!!!
      DO 7033 ISR=1,2
      DO 7034 ISL=1,2
      I=(ISR-1)*8+(ISL-1)*16+SEASON*2
7034  SIPH(ISL)=HPOL(HOUR,B0F(I-1),B0F(I),SAX,SUX,1.,1.)
7033  SIPL(ISR)=SIPH(1)+(SIPH(2)-SIPH(1))/DELA
      B0=SIPL(1)+(SIPL(2)-SIPL(1))/90.*(R-10.)
C!!!!!!! F1-REGION PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      F1REG=.FALSE.
      HMF1=0.
      NMF1=0.
      C1=0.
      XHI0=49.84733+0.349504*ABSLAT
      XHI100=38.96113+0.509932*ABSLAT
      XHIM=(XHI0+(XHI100-XHI0)*R/100.0)
      IF((XHI.GT.XHIM).OR.NIGHT.OR.(SEASON.EQ.4)) GOTO 150
      F1REG=.TRUE.
      C1=.09+.11/DELA
      FOF1=FOF1ED(MAGBR,R,COSXHI)
      NMF1=1.24E10*FOF1*FOF1
C!!!!!!! PARAMETER FOR E AND VALLEY-REGION !!!!!!!!!!!!!!!!!!!!!
150   NME=1.24E10*FOE*FOE
      HME=105.0
      XDEL=XDELS(SEASON)/DELA
      DNDHBR=DNDS(SEASON)/DELA
      HDEEP=HPOL(HOUR,10.5/DELA,28.,SAX,SUX,1.,1.)
      WIDTH=HPOL(HOUR,17.8/DELA,45.+22./DELA,SAX,SUX,1.,1.)
      DEPTH=HPOL(HOUR,XDEL,81.,SAX,SUX,1.,1.)
      DLNDH=HPOL(HOUR,DNDHBR,.06,SAX,SUX,1.,1.)
      IF(DEPTH.LT.1.0) GOTO 600
      IF(NIGHT) DEPTH=-DEPTH
      CALL TAL(HDEEP,DEPTH,WIDTH,DLNDH,EXT,E)
      IF(.NOT.EXT) GOTO 667
      WRITE(KONSOL,650)
650   FORMAT(1X,'*NE* E-REGION VALLEY CAN NOT BE MODELLED')
600   WIDTH=.0
667   HEF=HME+WIDTH
C!!!!!!!D-REGION PARAMETER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      NMD=XMDED(XHI,R,4.0E8)
      HMD=HPOL(HOUR,81.0,88.0,SAX,SUX,1.,1.)
      F(1)=HPOL(HOUR,0.02+0.03/DELA,0.05,SAX,SUX,1.,1.)
      F(2)=HPOL(HOUR,4.6,4.5,SAX,SUX,1.,1.)
      F(3)=HPOL(HOUR,-11.5,-4.0,SAX,SUX,1.,1.)
      FP1=F(1)
      FP2=-FP1*FP1/2.0
      FP30=(-F(2)*FP2-FP1+1.0/F(2))/(F(2)*F(2))
      FP3U=(-F(3)*FP2-FP1-1.0/F(3))/(F(3)*F(3))
      HDX=HMD+F(2)
      X=HDX-HMD
      XDX=NMD*EXP(X*(FP1+X*(FP2+X*FP30)))
      DXDX=XDX*(FP1+X*(2.0*FP2+X*3.0*FP30))
      X=HME-HDX
      XKK=-DXDX*X/(XDX*ALOG(XDX/NME))
      D1=DXDX/(XDX*XKK*X**(XKK-1.0))
C
C SEARCH FOR HMF1 AND HST.....................................
C
3801  IF(F1REG) GOTO 924
      HMF1=HMF2
      GOTO 380
924   H=10.0
133   H=H+10.0
      IF(H.GT.(HMF2-HEF)) GOTO 135
      CALL REGFA1(HMF2-H,HMF2,0.001,NMF1,XE2,SCHALT,HMF1)
      IF(SCHALT) GOTO 133
      GOTO 137
135   WRITE(KONSOL,11)
11    FORMAT(1X,'*NE* HMF1 IS NOT EVALUATED BY THE FUNCTION XE2')
3985  IF (HMF2-HEF.GE.B0) GOTO 922
      B0X=(HMF2-HEF)/1.1
      IF(B0X.LT.10.0) GOTO 922
      WRITE(KONSOL,923) B0,B0X
923   FORMAT(6X,'CORR.: B0(OLD)=',F5.1,' B0(NEW)=',F5.1)
      B0=B0X
      GOTO 924
922   B1X=B1+.5
      IF(B1X.GT.4) GOTO 7398
      WRITE(KONSOL,902) B1,B1X
902   FORMAT(6X,'CORR.: B1(OLD)=',F4.1,' B1(NEW)=',F4.1)
      B1=B1X
      GOTO 924
7398  WRITE(KONSOL,9269)
9269  FORMAT(1X,'CORR.: NO F1 REGION, B1=3, C1=0.0')
      HMF1=0.
      NMF1=0.
      C1=0.0
      B1=3.
      F1REG=.FALSE.
      GOTO 380
137   IF(HMF1.GT.HEF) GOTO 380
      WRITE(KONSOL,9)
9     FORMAT(1X,'*NE* UPPER E-VALLEY LIMIT(HEF) ABOVE HMF1')
      IF((HMF2-HEF.LT.B0).OR.(B1.LT.4)) GOTO 3985
      WRITE(KONSOL,9853)
9853  FORMAT(1X,'CORR.: NO F1 REGION, B1=3, C1=0.0')
      HMF1=0.
      NMF1=0.
      C1=0.0
      B1=3.
      F1REG=.FALSE.
C!!!!!!! NE3(HST)=NME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
380   H=-3.
125   H=H+3.0
      IF(H.GT.(HMF1-HEF)) GOTO 900
      CALL REGFA1(HEF+H,HMF1,0.001,NME,XE3,SCHALT,HST)
      STR=HST
      IF(SCHALT) GOTO 125
      GOTO 360
900   WRITE(KONSOL,100)
100   FORMAT(1X,'*NE* HST IS NOT EVALUATED BY THE FUNCTION XE3')
      IF (HMF2-HEF.GE.B0) GOTO 9221
      B0X=(HMF2-HEF)/1.1
      IF(B0X.LT.10.0) GOTO 9221
      WRITE(KONSOL,9231) B0,B0X
9231  FORMAT(6X,'CORR.: B0(OLD)=',F5.1,' B0(NEW)=',F5.1)
      B0=B0X
      GOTO 3801
9221  B1X=B1+.5
      IF(B1X.GT.5) GOTO 7391
      WRITE(KONSOL,9021) B1,B1X
9021  FORMAT(6X,'CORR.: B1(OLD)=',F4.1,' B1(NEW)=',F4.1)
      B1=B1X
      GOTO 3801
7391  CONTINUE
      HZ=HEF+.75*(HMF1-HEF)
      WRITE(KONSOL,901) HZ,HEF
901   FORMAT(6X,'CORR.: LIN. APP. BETWEEN HZ=',F5.1,
     &  ' AND HEF=',F5.1)
      XNEHZ=XE(HZ)
      T=(XNEHZ-NME)/(HZ-HEF)
      HST=-333.
      GOTO 6153
360   HZ=(HST+HMF1)/2.0
      D=HZ-HST
      T=D*D/(HZ-HEF-D)
6153  IF(.NOT.F1REG) HMF1=0.
C
C CALCULATION OF NEUTRAL TEMPERATURE PARAMETER...........
C
4933  HTA=120.0
      HTE=3000.0
      HTNA=90.0
      TNA=183.0
      ZX=125.0
      Z1=XLSTA
      HDEL=ZX-HTNA
      HD2=HDEL*HDEL
      TUN=TUNCAL(COV,LATI,SUNDEC,Z1)
      TNX=371.6678+0.0518806*TUN-294.3505*EXP(-0.00216222*TUN)
      ATN=0.63662*(TUN-TNX)
      TDEL=TNX-TNA
      CTN(1)=1.9*TDEL/HDEL
      CTN(3)=3.0*TDEL/(HD2*HD2)
      CTN(2)=CTN(3)*1.333333*HDEL-CTN(1)/HD2
C
C CALCULATION OF ELECTRON TEMPERATURE PARAMETER................
C
881   CONTINUE
C !!!!!!!!!! TE(120KM)=TN(120KM) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ATE(1)=TN(AHH(1),TNX,ATN,CTN)
C !!!!!!!!!! TE-MAXIMUM (JICAMARCA,ARECIBO) !!!!!!!!!!!!!!!!!!!!
      HM=60.*EXP(-(MLAT/22.41)**2)+210.
      AHH(2)=HPOL(HOUR,HM,150.,SAX,SUX,1.,1.)
      TM=800.*EXP(-(MLAT/33.)**2)+1500.
C !!!!!!!!!! TE(300,400KM)=TE-AE-C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!! TE(1400,3000KM)=TE-ISIS !!!!!!!!!!!!!!!!!!!!!!!!!!!
      DIPLAT=MAGBR
      IF(NSEASON.LT.3) THEN
	IS=NSEASON
      ELSE IF(NSEASON.GT.3) THEN
	IS=2
	DIPLAT=-MAGBR
      ELSE
	IS=1
      ENDIF
C      IS=SEASON-SEASON/3*2
      CALL TEBA(DIPLAT,HOUR,IS,TEA)
      ATE(2)=HPOL(HOUR,TM,TEA(5),SAX,SUX,1.,1.)
      ATE(3)=TEA(1)
      ATE(4)=TEA(2)
      ATE(6)=TEA(3)
      ATE(7)=TEA(4)
C !!!!!!!!!! TE(600KM)=TE-AEROS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ETT=EXP(-MLAT/11.35)
      TET=2900.-5600.*ETT/((ETT+1)**2.)
      TEN=839.+1161./(1.+EXP(-(ABS(MLAT)-45.)/5.))
      ATE(5)=HPOL(HOUR,TET,TEN,SAX,SUX,1.5,1.5)
C !!!!!!!!!! OPTION TO USE TE-NE-RELATION !!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!! AT 300, 400 OR 600 KM  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO 3395 I=1,3
3395  IF(XNAR(I).GT.0.) ATE(I+2)=TEDE(HOA(I),XNAR(I),-COV)
C !!!!!!!!!! TE'S ARE CORRECTED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!! ALSO TE > TN ENFORCED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      TNAHH2=TN(AHH(2),TNX,ATN,CTN)
      IF(ATE(2).LT.TNAHH2) ATE(2)=TNAHH2
      STTE1=(ATE(2)-ATE(1))/(AHH(2)-AHH(1))
      DO 1901 I=2,6
       TNAHHI=TN(AHH(I+1),TNX,ATN,CTN)
       IF(ATE(I+1).LT.TNAHHI) ATE(I+1)=TNAHHI
       STTE2=(ATE(I+1)-ATE(I))/(AHH(I+1)-AHH(I))
       ATE(I)=ATE(I)-(STTE2-STTE1)*DTE(I-1)*ALOG2
1901  STTE1=STTE2
C !!!!!!!!!! GRADIENTS ARE CALCULATED WITH !!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!! CORRECTED REGION BOUNDARIES !!!!!!!!!!!!!!!!!!!!!!
      DO 1902 I=1,6
1902  STTE(I)=(ATE(I+1)-ATE(I))/(AHH(I+1)-AHH(I))
      ATE1=ATE(1)
C
C CALCULATION OF ION TEMPERATURE PARAMETERS....................
C
887   CONTINUE
C !!!!!!!!!! TI(430KM,DAY)=TI-AEROS !!!!!!!!!!!!!!!!!!!!!!!!!!!
      XSM(1)=430.0
      EZ1=EXP(-0.09*MLAT)
      EZ2=EZ1+1.
      TID1 = 1240. - 1400. * EZ1 / (EZ2*EZ2)
      MM(2)=HPOL(HOUR,3.0,0.0,SAX,SUX,1.,1.)
C ----------- TN < TI < TE ------------------------------------
        TED1=TEA(6)+30       ! 90-60; 90 because 400; 60 for Epst.
        TND1=TUNCAL(COV,LATI,SUNDEC,0.0)
        IF(TED1.LT.TND1) TED1=TND1
        IF(TID1.GT.TED1) TID1=TED1
        IF(TID1.LT.TND1) TID1=TND1
C !!!!!!!!!! TI(430KM,NIGHT)=TI-AEROS !!!!!!!!!!!!!!!!!!!!!!!!!
      Z1=ABS(MLAT)
      Z2=Z1*(0.47+Z1*0.024)*UMR
      Z3=COS(Z2)
      TIN1=1200.0-300.0*SIGN(1.0,Z3)*SQRT(ABS(Z3))
C ----------- TN < TI < TE ------------------------------------
        TEN1=TEA(5)-30
        TNN1=TUNCAL(COV,LATI,SUNDEC,180.0)
        IF(TEN1.LT.TNN1) TEN1=TNN1
        IF(TIN1.GT.TEN1) TIN1=TEN1
        IF(TIN1.LT.TNN1) TIN1=TNN1
C !!!!!!!!!! TI(430KM,LT) FROM STEP FUNCTION !!!!!!!!!!!!!!!!!!
      TI1=TIN1  
      IF(TID1.GT.TIN1) TI1=HPOL(HOUR,TID1,TIN1,SAX,SUX,1.,1.)
C !!!!!!!!!! TANGENT ON TN DETERMINES HS !!!!!!!!!!!!!!!!!!!!!!
      CALL REGFA1(130.0,500.0,0.1,TI1,TEDER,SCHALT,HS)
      IF(SCHALT) HS=200.
      TNHS=TN(HS,TNX,ATN,CTN)
      MM(1)=DTNDH(HS,ATN,CTN)
      IF(SCHALT) MM(1)=(TI1-TNHS)/(XSM(1)-HS)
      MXSM=2
C !!!!!!!!!! XTETI ALTITTUDE WHERE TE=TI !!!!!!!!!!!!!!!!!!!!!!
2391    XTTS=500.
        X=500.
2390    X=X+XTTS
        IF(X.GE.AHH(7)) GOTO 240
        TEX=ELTE(X)
        TIX=TI(X)
        IF(TIX.LT.TEX) GOTO 2390
        X=X-XTTS
        XTTS=XTTS/10.
        IF(XTTS.GT.0.1) GOTO 2390
        XTETI=X+XTTS*5.
C !!!!!!!!!! TI=TE ABOVE XTETI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        MXSM=3
        MM(3)=STTE(6)
        XSM(2)=XTETI
        IF(XTETI.GT.AHH(6)) GOTO 240
        MXSM=4
        MM(3)=STTE(5)
        MM(4)=STTE(6)
        XSM(3)=AHH(6)
        IF(XTETI.GT.AHH(5)) GOTO 240
        MXSM=5
        DTI(1)=5.
        DTI(2)=5.
        MM(3)=STTE(4)
        MM(4)=STTE(5)
        MM(5)=STTE(6)
        XSM(3)=AHH(5)
        XSM(4)=AHH(6)
C
C INPUT OF THE ION DENSITY PARAMETER ARRAYS PF1O,PF2O AND PF3O......
C
240   IF(IOND.LT.1) GOTO 141
      RIF(1)=2.
      IF(ABSLAT.LT.30.0) RIF(1)=1.
      RIF(2)=2.
      IF(COV.LT.100.0) RIF(2)=1.
      RIF(3)=SEASON
      IF(SEASON.EQ.1) RIF(3)=3.
      RIF(4)=1.
      IF(NIGHT) RIF(4)=2.
      CALL KOEFP1(PG1O)
      CALL KOEFP2 (PG2O)
      CALL KOEFP3 (PG3O)
      CALL SUFE (PG1O,RIF,12,PF1O)
      CALL SUFE(PG2O,RIF,4,PF2O)
      CALL SUFE(PG3O,RIF,12,PF3O)
C
C CALCULATION OF ION DENSITY PARAMETER..................
C
      HNIA=100.
      HNIE=2000.
      ZZZ1=0.0
      IF(XHI.LE.90.0) ZZZ1=COSXHI
      HFIXO=300.0
      IF((RIF(2).EQ.2.).AND.(RIF(3).EQ.2.)) HFIXO=249.0
      MO(1)=EPSTEP(PF1O(1),PF1O(2),PF1O(3),PF1O(4),ZZZ1)
      MO(2)=EPSTEP(PF1O(5),PF1O(6),PF1O(7),PF1O(8),ZZZ1)
      MO(3)=0.0
      HO(1)=EPSTEP(PF1O(9),PF1O(10),PF1O(11),PF1O(12),ZZZ1)
      HO(4)=PF2O(1)
      MO(4)=PF2O(2)
      MO(5)=PF2O(3)
7100  HO(2)=290.0
      IF((RIF(2).EQ.2.).AND.(RIF(3).EQ.2.)) HO(2)=237.0
      HO(3)=(4.60517-MO(5)*(HO(4)-PF2O(4)))/MO(4)+HO(4)
      IF(HO(2).LT.HO(3)) GOTO 7101
      MO(4)=MO(4)-0.001
      GOTO 7100
7101  Z2=5.0
      X=HO(2)
      Z1=0.0
7102  X=X+Z2
      Y=RPID(X,HFIXO,98.0,4,MO,DDO,HO)
      IF(Y.LE.Z1) GOTO 7103
      Z1=Y
      GOTO 7102 
7103  IF(Z2.LE.1.0) GOTO 7104
      X=X-Z2
      Z2=1.0
      GOTO 7102
7104  H0O=X-0.5
      DO 7105 I=2,4,2
      L=I/2
      HO2(L)=PF3O(1+I)+PF3O(2+I)*ZZZ1
7105  MO2(L+1)=PF3O(7+I)+PF3O(8+I)*ZZZ1
      MO2(1)=PF3O(7)+PF3O(8)*ZZZ1
7106  Y=RPID(H0O,PF3O(1),PF3O(2),2,MO2,DO2,HO2)
      IF(Y.LE.0.1) GOTO 1899
      MO2(3)=MO2(3)-0.02
      GOTO 7106
C!!!!!!! RATIO OF NO+ TO O2+ DENSITY AT O+ MAXIMUM IS USED !!!!!!
C!!!!!!! FOR NO+ CALCULATION ABOVE THE O+ MAXIMUM (H0O) !!!!!!!!!
1899  Z1B= RPID (H0O,HFIXO,98.0,4,MO,DDO,HO)
      Z2B= RPID (H0O,PF3O(1),PF3O(2),2,MO2,DO2,HO2)
      NOBO2= (100.0-Z1B-Z2B)/Z2B
      IF(Z2B.LT.1.) NOBO2=0.0
C
C CALCULATION FOR THE REQUIRED HEIGHT RANGE.......................
C
141   IF(.NOT.F1REG) HMF1=HZ
      XCOR=XVAR(IVAR)
      IF(HEIGHT.LT.1.0) THEN 
       WRITE(AGNR,3910) XCOR,INT(HMF2+.5),
     &  INT(HMF1+.5),INT(HME+.5),INT(HMD+.5),INT(NMF2/1.E6+.5),
     &  INT(NMF1/1.E6+.5),INT(NME/1.E6+.5),INT(NMD/1.E6+.5)
       IF(JAGNR.EQ.2) WRITE(IUOUT,3910) XCOR,INT(HMF2+.5),
     &  INT(HMF1+.5),INT(HME+.5),INT(HMD+.5),INT(NMF2/1.E6+.5),
     &  INT(NMF1/1.E6+.5),INT(NME/1.E6+.5),INT(NMD/1.E6+.5)
3910   FORMAT(1X,F7.1,2X,4I5,3X,4I7)
       GOTO 2135
      ENDIF
      X=HEIGHT
300   DO 7397 KI=1,11
7397  OUTF(KI)=-1.6
      OUTF(1)=-1.
      OUTF(2)=-1.
      OUTF(6)=-1.
      IF((INECH.LT.1).OR.(X.LT.HNEA)) GOTO 330
      NEI=XE(X)
      IF(JF(1).GT.0) OUTF(1)=NEI/1.E6
      IF(JF(2).GT.0) OUTF(2)=NEI/NMF2
330   IF(KOMB.LT.1) GOTO 7108
      IF((X.GT.HTE).OR.(X.LT.HTA)) GOTO 7108
      TNH=TN(X,TNX,ATN,CTN)
      TIH=TNH
      IF(X.GE.HS) TIH=TI(X)
      TEH=ELTE(X)
      IF(JF(3).GT.0) OUTF(3)=TNH
      IF(JF(4).GT.0) OUTF(4)=TIH
      IF(JF(5).GT.0) OUTF(5)=TEH
      IF(JF(6).GT.0) OUTF(6)=TEH/TIH
7108  IF(IOND.LT.1) GOTO 7118
      IF((X.GT.HNIE).OR.(X.LT.HNIA)) GOTO 7118
      Z1=RPID(X,HFIXO,98.0,4,MO,DDO,HO)
      Z2=RPID(X,PF3O(1),PF3O(2),2,MO2,DO2,HO2)
      IF(JF(7).GT.0) OUTF(7)=Z1
      IF((JF(8).GT.0).OR.(JF(9).GT.0)) 
     & CALL RDHHE(X,H0O,Z1,Z2,NOBO2,10.,OUTF(8),OUTF(9))
      IF(JF(10).GT.0) OUTF(10)=Z2
      IF(JF(11).GT.0) OUTF(11)=RDNO(X,H0O,Z2,Z1,NOBO2)
7118  IF(IVAR.EQ.6) XCOR=X
*-------------------------- main output -----------------------------

      jj = 1
      if (jj .eq. 1) then
	Ne = int(outf(1)+0.5)
	Oplus = Ne*int(outf(7)+0.5)/100.0
	TiOut = int(outf(4) + .5)
	TeOut = int(outf(5) + .5)
	return						! this is some hack
      end if

      WRITE(AGNR,7117) XCOR,INT(OUTF(1)+.5),OUTF(2),INT(OUTF(3)+.5),
     &  INT(OUTF(4)+.5),INT(OUTF(5)+.5),OUTF(6),INT(OUTF(7)+.5),
     &  INT(OUTF(8)+.5),INT(OUTF(9)+.5),INT(OUTF(10)+.5),
     &  INT(OUTF(11)+.5)
      IF(JAGNR.EQ.2) WRITE(IUOUT,7117)
     &  XCOR,INT(OUTF(1)+.5),OUTF(2),INT(OUTF(3)+.5),
     &  INT(OUTF(4)+.5),INT(OUTF(5)+.5),OUTF(6),INT(OUTF(7)+.5),
     &  INT(OUTF(8)+.5),INT(OUTF(9)+.5),INT(OUTF(10)+.5),
     &  INT(OUTF(11)+.5)
7117  FORMAT(1X,F7.1,I9,F8.4,3I6,F7.2,5I5)
      IF(IVAR.EQ.6) THEN
        X=X+SVAR
        IF(X.LE.EVAR) GOTO 300
        GOTO 2289
      ENDIF
2135  HMF2=0.0
      MONTHO=MONTH
      RGO=RG
      IF(XCOR.LT.EVAR) GOTO 2123
2289  WRITE(AGNR,2193) LATI,LONGI,HEIGHT,R,MONTH,HOUR,MLAT,
     &  MLONG,DIP,COV,UT
      IF(JAGNR.EQ.2) WRITE(IUOUT,2193) LATI,LONGI,HEIGHT,R,MONTH,HOUR,
     &  MLAT,MLONG,DIP,COV,UT
2193  FORMAT(1X,74('-')/' LATI/LONG=',F5.1,'/',F5.1,'    H=',F6.1,
     &  '   RZ12=',F5.1,'  MONTH:',I2,'  L.T.:',F4.1/' MLAT/MLON=',
     &  F5.1,'/',F5.1,'  DIP=',F5.1,'   F10.7=',F5.1,12X,'U.T.:',
     &  F4.1/1X,74('-'))
      WRITE(MONITO,5600)
5600  FORMAT(1X/' **** DO YOU WANT TO CONTINUE?'/1X,60('-')/
     &  ' "0"  QUIT AND EXIT        "1"  NEW PARAMETERS')
	MAXI=1
      IF(IVAR.EQ.6) THEN
	WRITE(MONITO,5607)
5607    FORMAT(1X,'"2"  DIFFERENT ALTITUDE RANGE')
	MAXI=2
	ENDIF
      WRITE(MONITO,5666)
5666  FORMAT(1X,60('-'))
8672  READ(EGNR,*,ERR=8670,END=3330) IALL
      IF((IALL.GE.MIN0).AND.(IALL.LE.MAXI)) GOTO 8671
8670  	WRITE(MONITO,8101) MIN0,MAXI
      	GOTO 8672
8671  IF(IALL.EQ.0) GOTO 3330
      IF((IALL.EQ.2).AND.(IVAR.EQ.6)) THEN
        WRITE(MONITO,4618) BVAR,EVAR,SVAR
4618    FORMAT(1X,'BEGIN, END, STEPWIDTH   [KM]',10X,'#',F6.1,','
     &   ,F6.1,',',F6.1,'#')
      	WRITE(MONITO,8640)
	XMAX=PARMAX(6)
	XMIN=PARMIN(6)
8691 	READ(EGNR,*,ERR=8692,END=3330) BVAR,EVAR,SVAR
    	IF((EVAR.LE.XMAX).AND.(BVAR.GE.XMIN)) GOTO 8693
8692  	WRITE(MONITO,8117) XMIN,XMAX
     	GOTO 8691
8693    X=BVAR
        WRITE(AGNR,1793) ITEXT(IVAR),IMZ(IVAR)
        IF(JAGNR.EQ.2) WRITE(IUOUT,1793) ITEXT(IVAR),IMZ(IVAR)
1793  FORMAT(//////////3X,A5,'  ELECTRON DENSITY      TEMPERATURES',
     &  7X,'ION PERCENTAGE DENSITIES'/4X,A4,'  NE/CM-3 NE/NMF2  ',
     &  'TN/K  TI/K  TE/K  TE/TI   O+   H+  He+  O2+  NO+')
        GOTO 300
      ENDIF
      GOTO 3293
3330  CONTINUE
      STOP
      END

