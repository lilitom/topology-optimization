*USER SUBROUTINES
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
	INCLUDE 'ABA_PARAM.INC'
C     
      REAL Et, Stress_Max, EkEk

       CHARACTER*8 CMNAME
       DIMENSION STRESS(NTENS),STATEV(NSTATEV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 Et(3),Stress_Max(3,3),EkEk(1),DFGRDM1_INV(3,3),En(3),
     5 TEMP(1),DTEMP(1),DFGRDP(3),DFGRDM1(3, 3),DFGRDM0(3, 3),
     6 CBAR(3,3),DDSDDE_Max(6,6)
C
      DIMENSION BBAR(6),DISTGR(3,3)
C
      PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, FOUR=4.D0)
C
C ----------------------------------------------------------------
C    UMAT FOR COMPRESSIBLE NEO-HOOKEAN HYPERELASTICITY
C    CANNOT BE USED FOR PLANE STRESS
C ----------------------------------------------------------------
C    PROPS(1) - C10
C    PROPS(2) - D1
C    PROPS(3) - EPSILON
C    PROPS(4) - En(1)/abs(En)
C    PROPS(5) - En(2)/abs(En)
C    PROPS(6) - En(3)/abs(En)
C    PROPS(7) - DFGRDP(1,1)
C    PROPS(8) - DFGRDP(2,2)
C    PROPS(9) - DFGRDP(3,3)
C ----------------------------------------------------------------
C
C   Read Material PROPERTIES
C
      C10=PROPS(1)
      D1 =PROPS(2)
      EPSILON=PROPS(3)
      
C
C   Read The Prestretch
C      
      do J=1,3
      DFGRDP(J)=0
      end do  
           
      IF (KSTEP==1) then
      DFGRDP(1)=exp(log(PROPS(7))*time(2)) 
      DFGRDP(2)=exp(log(PROPS(8))*time(2))
      DFGRDP(3)=exp(log(PROPS(9))*time(2))
      else
      DFGRDP(1)=PROPS(7) 
      DFGRDP(2)=PROPS(8)
      DFGRDP(3)=PROPS(9)
      endif
C
C  Update the deformation gradient tensor
C     
      do I=1,3
      do J=1,3
      DFGRDM1(I, J)=DFGRD1(I, J)*DFGRDP(J)
      DFGRDM0(I, J)=DFGRD0(I, J)*DFGRDP(J)
      enddo
      enddo
      
C
C    JACOBIAN AND DISTORTION TENSOR
C
      DET=DFGRDM1(1, 1)*DFGRDM1(2, 2)*DFGRDM1(3, 3)
     1   -DFGRDM1(1, 2)*DFGRDM1(2, 1)*DFGRDM1(3, 3)
      IF(NSHR.EQ.3) THEN
        DET=DET+DFGRDM1(1, 2)*DFGRDM1(2, 3)*DFGRDM1(3, 1)
     1         +DFGRDM1(1, 3)*DFGRDM1(3, 2)*DFGRDM1(2, 1)
     2         -DFGRDM1(1, 3)*DFGRDM1(3,1)*DFGRDM1(2, 2)
     3         -DFGRDM1(2, 3)*DFGRDM1(3, 2)*DFGRDM1(1, 1)
      END IF
      SCALE=DET**(-ONE/THREE)
      DO K1=1, 3
        DO K2=1, 3
          DISTGR(K2, K1)=SCALE*DFGRDM1(K2, K1)
        END DO
      END DO
C
C    CALCULATE LEFT CAUCHY-GREEN TENSOR
C
      BBAR(1)=DISTGR(1, 1)**2+DISTGR(1, 2)**2+DISTGR(1, 3)**2
      BBAR(2)=DISTGR(2, 1)**2+DISTGR(2, 2)**2+DISTGR(2, 3)**2
      BBAR(3)=DISTGR(3, 3)**2+DISTGR(3, 1)**2+DISTGR(3, 2)**2
      BBAR(4)=DISTGR(1, 1)*DISTGR(2, 1)+DISTGR(1, 2)*DISTGR(2, 2)
     1       +DISTGR(1, 3)*DISTGR(2, 3)
      IF(NSHR.EQ.3) THEN
        BBAR(5)=DISTGR(1, 1)*DISTGR(3, 1)+DISTGR(1, 2)*DISTGR(3, 2)
     1         +DISTGR(1, 3)*DISTGR(3, 3)
        BBAR(6)=DISTGR(2, 1)*DISTGR(3, 1)+DISTGR(2, 2)*DISTGR(3, 2)
     1         +DISTGR(2, 3)*DISTGR(3, 3)
      END IF      
C
C    CALCULATE THE STRESS
C
      TRBBAR=(BBAR(1)+BBAR(2)+BBAR(3))/THREE
      EG=TWO*C10/DET
      EK=TWO/D1*(TWO*DET-ONE)
      PR=TWO/D1*(DET-ONE)
      DO K1=1,NDI
        STRESS(K1)=EG*(BBAR(K1)-TRBBAR)+PR
      END DO
      DO K1=NDI+1,NDI+NSHR
        STRESS(K1)=EG*BBAR(K1)
      END DO
      TMPT=STRESS(3)
C
C    CALCULATE THE STIFFNESS DUE TO ELASTICITY 
C
      EG23=EG*TWO/THREE
      DDSDDE(1, 1)= EG23*(BBAR(1)+TRBBAR)+EK
      DDSDDE(2, 2)= EG23*(BBAR(2)+TRBBAR)+EK
      DDSDDE(3, 3)= EG23*(BBAR(3)+TRBBAR)+EK
      DDSDDE(1, 2)=-EG23*(BBAR(1)+BBAR(2)-TRBBAR)+EK
      DDSDDE(1, 3)=-EG23*(BBAR(1)+BBAR(3)-TRBBAR)+EK
      DDSDDE(2, 3)=-EG23*(BBAR(2)+BBAR(3)-TRBBAR)+EK
      DDSDDE(1, 4)= EG23*BBAR(4)/TWO
      DDSDDE(2, 4)= EG23*BBAR(4)/TWO
      DDSDDE(3, 4)=-EG23*BBAR(4)
      DDSDDE(4, 4)= EG*(BBAR(1)+BBAR(2))/TWO
      IF(NSHR.EQ.3) THEN
      DDSDDE(1, 5)= EG23*BBAR(5)/TWO
      DDSDDE(2, 5)=-EG23*BBAR(5)
      DDSDDE(3, 5)= EG23*BBAR(5)/TWO
      DDSDDE(1, 6)=-EG23*BBAR(6)
      DDSDDE(2, 6)= EG23*BBAR(6)/TWO
      DDSDDE(3, 6)= EG23*BBAR(6)/TWO
      DDSDDE(5, 5)= EG*(BBAR(1)+BBAR(3))/TWO
      DDSDDE(6, 6)= EG*(BBAR(2)+BBAR(3))/TWO
      DDSDDE(4,5)= EG*BBAR(6)/TWO
      DDSDDE(4,6)= EG*BBAR(5)/TWO
      DDSDDE(5,6)= EG*BBAR(4)/TWO
      END IF
      DO K1=1, NTENS
        DO K2=1, K1-1
          DDSDDE(K1, K2)=DDSDDE(K2, K1)
        END DO
      END DO
C
C     Calculate the inverse of deformation gradient
C
      DFGRDM1_INV(1,1)=DFGRDM1(2,2)*DFGRDM1(3,3)
     1 -DFGRDM1(2,3)*DFGRDM1(3,2) 
      DFGRDM1_INV(1,2)=-DFGRDM1(1,2)*DFGRDM1(3,3)
     1 +DFGRDM1(1,3)*DFGRDM1(3,2)  
      DFGRDM1_INV(1,3)=DFGRDM1(1,2)*DFGRDM1(2,3)
     1 -DFGRDM1(1,3)*DFGRDM1(2,2)
 
      DFGRDM1_INV(2,1)=-DFGRDM1(2,1)*DFGRDM1(3,3)
     1 +DFGRDM1(2,3)*DFGRDM1(3,1) 
      DFGRDM1_INV(2,2)=DFGRDM1(1,1)*DFGRDM1(3,3)
     1 -DFGRDM1(1,3)*DFGRDM1(3,1) 
      DFGRDM1_INV(2,3)=-DFGRDM1(1,1)*DFGRDM1(2,3)
     1 +DFGRDM1(1,3)*DFGRDM1(2,1)

      DFGRDM1_INV(3,1)=DFGRDM1(2,1)*DFGRDM1(3,2)
     1 -DFGRDM1(2,2)*DFGRDM1(3,1) 
      DFGRDM1_INV(3,2)=-DFGRDM1(1,1)*DFGRDM1(3,2)
     1 +DFGRDM1(1,2)*DFGRDM1(3,1)  
      DFGRDM1_INV(3,3)=DFGRDM1(1,1)*DFGRDM1(2,2)
     1 -DFGRDM1(1,2)*DFGRDM1(2,1)
       
      do I=1,3
      do J=1,3
      DFGRDM1_INV(I,J)=DFGRDM1_INV(I,J)/DET
      end do
      end do
C      
C     Calculate the current true electric field
C
      En(1)=(TEMP(1)+DTEMP(1))*PROPS(4)     
      En(2)=(TEMP(1)+DTEMP(1))*PROPS(5)
      En(3)=(TEMP(1)+DTEMP(1))*PROPS(6)  
	  
      open(unit=60, file='F:\Lambda.csv')
      write(60,*) DFGRDM1(1,1),En(3)/sqrt(2*12550/4.1595e-11)
	  
      do I=1,3
      Et(I)=0
      end do 
      
      do I=1,3
      do J=1,3
      Et(J)=Et(J)+En(I)*DFGRDM1_INV(I,J)
      end do
      end do

      EkEk=0.0D0
      DO I=1, 3
         EkEk=EkEk+Et(I)*Et(I)
      ENDDO 

      DO I=1,3
        DO J=1,3
        Stress_Max(I,J)=0.0D0
        ENDDO
      ENDDO
C
C     Update the true stress due to polarization
C
      DO I=1,3
        DO J=1,3
        IF (I==J) THEN
        Stress_Max(I,J)=EPSILON*Et(I)*Et(J)-0.50D0*EPSILON*EkEk(1)
        ELSE
        Stress_Max(I,J)=EPSILON*Et(I)*Et(J)
        ENDIF
        ENDDO
      ENDDO
      
      STRESS(1)=STRESS(1)+Stress_Max(1,1)
      STRESS(2)=STRESS(2)+Stress_Max(2,2)
      STRESS(3)=STRESS(3)+Stress_Max(3,3)
      STRESS(4)=STRESS(4)+Stress_Max(1,2)
      STRESS(5)=STRESS(5)+Stress_Max(1,3)
      STRESS(6)=STRESS(6)+Stress_Max(2,3)
      
      RETURN
      END 


