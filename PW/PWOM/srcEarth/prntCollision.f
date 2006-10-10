
C ============================================================================

      SUBROUTINE prntCollision

      use ModCommonVariables

      write (iUnitCollision,"(a79)") 'Collision Frequencies'
      write (iUnitCollision,"(i7,1pe13.5,3i3)") nint(time/dt),time,1,1,12
      write (iUnitCollision,"(3i4)") nDim
      write (iUnitCollision,"(100(1pe13.5))") Gamma
      write (iUnitCollision,"(a79)") 
     &'alt OXN2 OXO2 OXO OXHE OXH OXHD OXHL OXEL HN2 HO2 HHE HO gamma'

C CFHH CFHOX CFHHL CFHEL CFHEN2 CFHEO2 CFHEHE CFHEO CFHEH CFHEO CFHEHD CFHEEL 
C CFELN2 CFELO2 CFELHE CFELO CFELH CFELOX CFELHL CFELHD gamma'
      
      DO K=1,NDIM
         QS1=CFOXN2(K) 
         QS2=CFOXO2(K) 
         QS3=CFOXO(K) 
         QS4=CFOXHE(K) 
         QS5=CFOXH(K) 
         QS6=CFOXHD(K) 
         QS7=CFOXHL(K) 
         QS8=CFOXEL(K) 
         QS9=CFHN2(K) 
         QS10=CFHO2(K) 
         QS11=CFHHE(K) 
         QS12=CFHO(K) 
C         QS13=CFHH(K) 
C         QS14=CFHOX(K) 
C         QS15=CFHHL(K) 
C         QS16=CFHEL(K) 
C         QS17=CFHEN2(K) 
C         QS18=CFHEO2(K) 
C         QS19=CFHEHE(K) 
C         QS20=CFHEO(K) 
C         QS21=CFHEH(K) 
C         QS22=CFHEOX(K) 
C         QS23=CFHEHD(K) 
C         QS24=CFHEEL(K) 
C         QS25=CFELN2(K) 
C         QS26=CFELO2(K) 
C         QS27=CFELHE(K) 
C         QS28=CFELO(K) 
C         QS29=CFELH(K) 
C         QS30=CFELOX(K) 
C         QS31=CFELHL(K) 
C         QS32=CFELHD(K)


         WRITE (iUnitCollision,"(100(1pe18.10))") 
     &        ALTD(K),QS1,QS2,QS3,QS4,QS5,QS6,QS7,QS8,QS9,QS10,QS11,QS12, gamma

C     &        QS13,QS14,QS15,QS16,QS17,QS18,QS19,QS20,QS21,QS22,QS23,QS24,
C     &        QS25,QS26,QS27,QS28,QS29,QS30,QS31,QS32,gamma
      enddo

      RETURN
      END
