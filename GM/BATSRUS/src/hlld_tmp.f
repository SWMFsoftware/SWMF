      SUBROUTINE HLLD_TMP(PrimLeft_V,PrimRight_V,Flux_V)

      use ModFaceFlux, ONLY: DoTestCell
      use ModVarIndexes
      use ModPhysics, ONLY: gamma => g

      real, intent(in):: PrimLeft_V(nVar), PrimRight_V(nVar)
      real, intent(out):: Flux_V(nVar+1)

      DIMENSION PL(nVar),PR(nVar),FM(nVar),F(nVar+1)
      DIMENSION FFL(nVar),FFR(nVar)

C     Change index order
      PL(1)      =PrimLeft_V(p_)
      PL(Ux_:Bz_)=PrimLeft_V(Ux_:Bz_)
      PL(8)      =PrimLeft_V(Rho_)

      PR(1)      =PrimRight_V(p_)
      PR(Ux_:Bz_)=PrimRight_V(Ux_:Bz_)
      PR(8)      =PrimRight_V(Rho_)

      PRL=PL(1)
      VXL=PL(Ux_)
      VYL=PL(Uy_)
      VZL=PL(Uz_)
      BXL=PL(Bx_)
      BYL=PL(By_)
      BZL=PL(Bz_)
      DEL=PL(8)

      PRR=PR(1)
      VXR=PR(Ux_)
      VYR=PR(Uy_)
      VZR=PR(Uz_)
      BXR=PR(Bx_)
      BYR=PR(By_)
      BZR=PR(Bz_)
      DER=PR(8)

      if(DoTestCell)then
         write(*,*)'hlld: PL=',PL
         write(*,*)'hlld: PR=',PR
      endif

      CALL SOUND(PRL,DEL,BXL,BYL,BZL,CL,CSL,CFL,CAL)
      CALL SOUND(PRR,DER,BXR,BYR,BZR,CR,CSR,CFR,CAR)

C      SL =MIN(VXL,VXR)-MAX(CFL,CFR)
C      SR =MAX(VXL,VXR)+MAX(CFL,CFR)
      SL =MIN(VXL-CFL,VXR-CFR)
      SR =MAX(VXL+CFL,VXR+CFR)

      if(DoTestCell)write(*,*)'hlld SL,SR=',SL,SR

      IF(SL.GT.0.) THEN
        CALL PRIMF(PL,FM)
        GOTO 100
      END IF

      IF(SR.LT.0.) THEN
        CALL PRIMF(PR,FM)
        GOTO 100
      END IF

      PBL=0.5*(BXL**2+BYL**2+BZL**2)
      PBR=0.5*(BXR**2+BYR**2+BZR**2)
      PKL=0.5*DEL*(VXL**2+VYL**2+VZL**2)
      PKR=0.5*DER*(VXR**2+VYR**2+VZR**2)

      PTL=PRL+PBL
      PTR=PRR+PBR
      BX =(BXL+BXR)/2.
      S=SIGN(1.,BX)

      EL =PRL/(GAMMA-1.)+PBL+PKL
      ER =PRR/(GAMMA-1.)+PBR+PKR
      VBL=VXL*BXL+VYL*BYL+VZL*BZL
      VBR=VXR*BXR+VYR*BYR+VZR*BZR

      DXR=SR-VXR
      DXL=SL-VXL
      SM =(DXR*DER*VXR-DXL*DEL*VXL-PTR+PTL)/(DXR*DER-DXL*DEL)

      DLA =DEL*DXL/(SL-SM)
      DRA =DER*DXR/(SR-SM)
      SLA=SM-ABS(BXL)/SQRT(DLA)
      SRA=SM+ABS(BXR)/SQRT(DRA)

      DENL=DEL*DXL*(SL-SM)-BXL**2
      DENR=DER*DXR*(SR-SM)-BXR**2
      if(denl.eq.0.) then
        VLA =VYL
        WLA =VZL
        BYLA=BYL
        BZLA=BZL
      else
        VLA =VYL-BXL*BYL*(SM-VXL)/DENL
        WLA =VZL-BXL*BZL*(SM-VXL)/DENL
        BYLA=BYL*(DEL*DXL**2-BXL**2)/DENL
        BZLA=BZL*(DEL*DXL**2-BXL**2)/DENL
      end if
      if(denr.eq.0.) then
        VRA =VYR
        WRA =VZR
        BYRA=BYR
        BZRA=BZR
      else
        VRA =VYR-BXR*BYR*(SM-VXR)/DENR
        WRA =VZR-BXR*BZR*(SM-VXR)/DENR
        BYRA=BYR*(DER*DXR**2-BXR**2)/DENR
        BZRA=BZR*(DER*DXR**2-BXR**2)/DENR
      end if
      PTA =(DXR*DER*PTL-DXL*DEL*PTR+DEL*DER*DXR*DXL*(VXR-VXL))
     +        /(DXR*DER-DXL*DEL)

      IF(SL.LE.0..AND.SLA.GE.0.) THEN
        VBLA=SM*BXL+VLA*BYLA+WLA*BZLA
        ELA =(SM*(EL+PTA)-VXL*(EL+PTL)+BXL*(VBL-VBLA))/(SL-SM)
        CALL PRIMF(PL,FFL)
        FM(1)=FFL(1)+SL*ELA
        FM(2)=FFL(2)+SL*(SM*DLA-VXL*DEL)
        FM(3)=FFL(3)+SL*(VLA*DLA-VYL*DEL)
        FM(4)=FFL(4)+SL*(WLA*DLA-VZL*DEL)
        FM(5)=0.
        FM(6)=FFL(6)+SL*(BYLA-BYL)
        FM(7)=FFL(7)+SL*(BZLA-BZL)
        FM(8)=FFL(8)+SL*(DLA-DEL)
        GOTO 100
      END IF
      IF(SRA.LE.0..AND.SR.GE.0.) THEN
        VBRA=SM*BXR+VRA*BYRA+WRA*BZRA
        ERA =(SM*(ER+PTA)-VXR*(ER+PTR)+BXR*(VBR-VBRA))/(SR-SM)
        CALL PRIMF(PR,FFR)
        FM(1)=FFR(1)+SR*ERA
        FM(2)=FFR(2)+SR*(SM*DRA-VXR*DER)
        FM(3)=FFR(3)+SR*(VRA*DRA-VYR*DER)
        FM(4)=FFR(4)+SR*(WRA*DRA-VZR*DER)
        FM(5)=0.
        FM(6)=FFR(6)+SR*(BYRA-BYR)
        FM(7)=FFR(7)+SR*(BZRA-BZR)
        FM(8)=FFR(8)+SR*(DRA-DER)
        GOTO 100
      END IF

      SDLA=SQRT(DLA)
      SDRA=SQRT(DRA)
      DENA=SDLA+SDRA
      SSD=SDLA*SDRA

      VAA =(SDLA*VLA +SDRA*VRA +(BYRA-BYLA)*S)/DENA
      WAA =(SDLA*WLA +SDRA*WRA +(BZRA-BZLA)*S)/DENA
      BYAA=(SDLA*BYRA+SDRA*BYLA+(VRA -VLA )*S*SSD)/DENA
      BZAA=(SDLA*BZRA+SDRA*BZLA+(WRA -WLA )*S*SSD)/DENA
      VBAA=SM*BX+VAA*BYAA+WAA*BZAA

      IF(SLA.LE.0..AND.SM.GE.0.) THEN
        CALL PRIMF(PL,FFL)
        VBLA=SM*BXL+VLA*BYLA+WLA*BZLA
        ELA =(SM*(EL+PTA)-VXL*(EL+PTL)+BXL*(VBL-VBLA))/(SL-SM)
        ELAA =-SQRT(DLA)*(VBLA-VBAA)*S
        FM(1)=FFL(1)+SLA*ELAA         +SL*ELA
        FM(2)=FFL(2)                  +SL*(SM*DLA-VXL*DEL)
        FM(3)=FFL(3)+SLA*DLA*(VAA-VLA)+SL*(VLA*DLA-VYL*DEL)
        FM(4)=FFL(4)+SLA*DLA*(WAA-WLA)+SL*(WLA*DLA-VZL*DEL)
        FM(5)=0.
        FM(6)=FFL(6)+SLA*(BYAA-BYLA)  +SL*(BYLA-BYL)
        FM(7)=FFL(7)+SLA*(BZAA-BZLA)  +SL*(BZLA-BZL)
        FM(8)=FFL(8)                  +SL*(DLA-DEL)
        GOTO 100
      END IF
      IF(SM.LE.0..AND.SRA.GE.0.) THEN
        CALL PRIMF(PR,FFR)
        VBRA=SM*BXR+VRA*BYRA+WRA*BZRA
        ERA =(SM*(ER+PTA)-VXR*(ER+PTR)+BXR*(VBR-VBRA))/(SR-SM)
        ERAA =SQRT(DRA)*(VBRA-VBAA)*S
        FM(1)=FFR(1)+SRA*ERAA         +SR*ERA
        FM(2)=FFR(2)                  +SR*(SM*DRA-VXR*DER)
        FM(3)=FFR(3)+SRA*DRA*(VAA-VRA)+SR*(VRA*DRA-VYR*DER)
        FM(4)=FFR(4)+SRA*DRA*(WAA-WRA)+SR*(WRA*DRA-VZR*DER)
        FM(5)=0.
        FM(6)=FFR(6)+SRA*(BYAA-BYRA)  +SR*(BYRA-BYR)
        FM(7)=FFR(7)+SRA*(BZAA-BZRA)  +SR*(BZRA-BZR)
        FM(8)=FFR(8)                  +SR*(DRA-DER)
        GOTO 100
      END IF
      PRINT*,'ERROR IN HLLD'
      STOP
 
 100  CONTINUE

      if(DoTestCell)write(*,*)'hlld: FM=',FM

      Flux_V(Energy_) = FM(1)
      Flux_V(Ux_)     = FM(2)
      Flux_V(Uy_)     = FM(3)
      Flux_V(Uz_)     = FM(4)
      Flux_V(Bx_)     = FM(5)
      Flux_V(By_)     = FM(6)
      Flux_V(Bz_)     = FM(7)
      Flux_V(Rho_)    = FM(8)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 

      SUBROUTINE SOUND(P,D,BX,BY,BZ,C,CS,CF,CA)

      use ModPhysics, ONLY: gamma => g

C     sound velocity

      C2=GAMMA*P/D
      C=SQRT(C2)

C     alfven velocity

      CA2=BX**2/D
      CA=SQRT(CA2)

C     fast magnetosonic velocity

      BT2=BY**2+BZ**2
      CT2=BT2/D
      RAIZ=SQRT(((C-CA)**2+CT2)*((C+CA)**2+CT2))

      CF=SQRT(0.5*(C2+CT2+CA2+RAIZ))

      CF=MAX(CF,CA)
      CF=MAX(CF,C)

C     slow magnetosonic velocity

      CS=CA*C/CF

      IF(CS.EQ.CF) THEN
        PRINT*,'CS=CF!!!'
        STOP
      END IF

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 

      SUBROUTINE PRIMF(P,F)

      use ModVarIndexes, ONLY: nVar
      use ModPhysics, ONLY: gamma => g

C    calculates the fluxes as a function of the primitive
C    variables

      
      DIMENSION F(nVar),P(nVar)
 
      BTOT=P(5)**2+P(6)**2+P(7)**2
      VV=P(2)**2+P(3)**2+P(4)**2
      F(1)=(P(1)*GAMMA/(GAMMA-1.)+BTOT+0.5*P(8)*VV)*P(2)
     &       -P(5)*(P(5)*P(2)+P(6)*P(3)+P(7)*P(4))
      F(2)=P(2)*P(2)*P(8)-P(5)*P(5)+P(1)+BTOT*0.5
      F(3)=P(2)*P(3)*P(8)-P(5)*P(6)
      F(4)=P(2)*P(4)*P(8)-P(5)*P(7)
      F(5)=0.
      F(6)=P(2)*P(6)-P(5)*P(3)
      F(7)=P(2)*P(7)-P(5)*P(4)
      F(8)=P(2)*P(8)

      RETURN
      END

