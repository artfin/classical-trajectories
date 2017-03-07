C PES for CO2-Ar by Hutson (J.Chem.Phys.,1996,V.105(20),9130-9140)
C 1) call POTDAT4() once (initialization)
C 2) call EXTPOT(R,Theta) repeatedly
C R: Bohr; Theta: radian; Result: Hartree
C (in other words, everything in atomic units)
C ---
C original PES code edited first by S.V.Ivanov,
C then by S.E.Lokshtanov (Mar 2017)... Hope it works!
C
C-------------------------------------------------------
	SUBROUTINE POTDAT4
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/POTHUT/ RLAM(20),BLAM(20),OLAM(20),
	1RDISPC,RREF,BFAC,RDISPO,BOND,C8O2,ISITE,C10O,C8O,C6O,C10C2,C8C2,
     2C6C2,C10C,C8C,C6C,SITE,C10O2,SFAC2,SFAC0,C6O2,SFACE,SFACC,
     3NLEGC,NLEGO
C
      DATA ISITE,SFAC0,SFAC2,BFAC/1, 4.6454D0, 0.6288D0, 1.D0/
      DATA NLEGC,NLEGO,RREF,BOND/5, 5, 7.D0, 2.196D0/
	DATA (RLAM(I),BLAM(I),I=1,10)/
	1 32.27581D0, 2.03171D0,
	2  0.D0     ,      0.D0,
	3 11.85757D0, 1.79571D0,
	4  0.D0     ,      0.D0,
	5 -9.41362D0, 2.09529D0,
	6 32.36857D0, 2.12256D0,
	7 -3.40252D0, 1.87293D0,
	8 -1.68672D0, 2.18922D0,
	9 -1.91165D0, 1.96196D0,
	1  1.10092D0, 1.98533D0/
      DATA C6,C6AFAC/101.D0, 0.270D0/
      DATA SITE,RDISP,FRACC/1.5874D0, -1.D0, 0.000000001D0/
      DATA C8RAT,C8AFAC/29.21D0, 999.D0/
C
        DO 100 I=1,NLEGC+NLEGO
  100   RLAM(I)=RLAM(I)*1.D-6
C
      C6C=FRACC*C6
      C6O=0.5D0*(C6-C6C)
      FRACC2=FRACC
      C62=C6AFAC*C6
      C6C2=FRACC2*C62
      C6O2=0.5D0*(C62-C6C2)
      C10RAT=49.D0/40.D0
      C8RATH=19.1412D0
      C8RATC=C8RAT
      C8C=C8RATC*C6C
      C10C=C10RAT*C8RATC*C8C
      C8C2=C8RATC*C6C2
      C10C2=C10RAT*C8RATC*C8C2
      RDISPC=RDISP
      IF(RDISP.LT.0.D0) RDISPC=10.01D0*DSQRT(C8RATC/C8RATH)
      C8RATO=C8RAT
      C8O=C8RATO*C6O
      C10O=C10RAT*C8RATO*C8O
      C8O2=C8RATO*C6O2
      C8=C8C+C8O+C8O+SITE**2*(10.D0*C6O+2.D0*C6O2)
      C10O2=C10RAT*C8RATO*C8O2
      RDISPO=RDISP
      IF(RDISP.LT.0.D0) RDISPO=10.01D0*DSQRT(C8RATO/C8RATH)
      C82=C8C+C8O+C8O+SITE**2*(32.D0*C6O+(88.D0/7.D0)*C6O2)
      IF(ISITE.GT.0) THEN	 ! Split repulsion potential
        SFACC=SFAC0+0.5D0*SFAC2
        SFACE=SFAC0-0.5D0*SFAC2
      ENDIF
C
	RETURN
	END
C------------------------------------------------------------
      FUNCTION EXTPOT(R,COSTH)
C------------------------------------------------------------
	IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/POTHUT/ RLAM(20),BLAM(20),OLAM(20),
	1RDISPC,RREF,BFAC,RDISPO,BOND,C8O2,ISITE,C10O,C8O,C6O,C10C2,C8C2,
     2C6C2,C10C,C8C,C6C,SITE,C10O2,SFAC2,SFAC0,C6O2,SFACE,SFACC,
     3NLEGC,NLEGO
      DATA AUE/219474.6354D0/
      DATA AUR/0.52917706D0/
C	BC=1.38066D-2 ! Boltzmann constant 
      RR=R   ! /AUR  ! Angstrom -> a.u.
      IF(ISITE.GT.0) THEN	 ! Split repulsion potential 
        T=-BOND/RR
        TEMP=DSQRT(1.D0+T*(T+COSTH+COSTH))
        RO1=RR*TEMP
        COSO1=(COSTH+T)/TEMP
C
        T=+BOND/RR
        TEMP=DSQRT(1.D0+T*(T+COSTH+COSTH))
        RO2=RR*TEMP
        COSO2=(COSTH+T)/TEMP
C
        DO 6 I=1,NLEGC
    6   OLAM(I)=RLAM(I)*DEXP(BFAC*BLAM(I)*(RREF-RR))
        OLAPC=SUMLEG(OLAM,NLEGC,COSTH)
        DO 7 I=1,NLEGO
    7   OLAM(I)=RLAM(NLEGC+I)*DEXP(BFAC*BLAM(NLEGC+I)*(RREF-RO1))
        OLAPO1=SUMLEG(OLAM,NLEGO,COSO1)
        DO 8 I=1,NLEGO
    8   OLAM(I)=RLAM(NLEGC+I)*DEXP(BFAC*BLAM(NLEGC+I)*(RREF-RO2))
        OLAPO2=SUMLEG(OLAM,NLEGO,-COSO2)
        V=SFACC*OLAPC+SFACE*(OLAPO1+OLAPO2)
      ENDIF
C
      T=-SITE/RR
      TEMP=DSQRT(1.D0+T*(T+COSTH+COSTH))
      RO1=RR*TEMP
      COSO1=(COSTH+T)/TEMP
C
      T=+SITE/RR
      TEMP=DSQRT(1.D0+T*(T+COSTH+COSTH))
      RO2=RR*TEMP
      COSO2=(COSTH+T)/TEMP
C
      R2=1.D0/(RR*RR)
      DISP=-(C6C+(C8C+C10C*R2)*R2
     1   +(C6C2+(C8C2+C10C2*R2)*R2)*(1.5D0*COSTH*COSTH-0.5D0))*R2**3
      IF(RR.GT.RDISPC) THEN
        V=V+DISP
      ELSE
        V=V+DISP*DEXP(-0.4D0*(RDISPC/RR-1.D0)**2)
      ENDIF
C
      R2=1.D0/(RO1*RO1)
      DISP=-(C6O+(C8O+C10O*R2)*R2
     1   +(C6O2+(C8O2+C10O2*R2)*R2)*(1.5D0*COSO1*COSO1-0.5D0))*R2**3
      IF(RO1.GT.RDISPO) THEN
        V=V+DISP
      ELSE
        V=V+DISP*DEXP(-0.4D0*(RDISPO/RO1-1.D0)**2)
      ENDIF
C
      R2=1.D0/(RO2*RO2)
      DISP=-(C6O+(C8O+C10O*R2)*R2
     1   +(C6O2+(C8O2+C10O2*R2)*R2)*(1.5D0*COSO2*COSO2-0.5D0))*R2**3
      IF(RO2.GT.RDISPO) THEN
        V=V+DISP
      ELSE
        V=V+DISP*DEXP(-0.4D0*(RDISPO/RO2-1.D0)**2)
      ENDIF
C
      EXTPOT=V  ! *AUE ! to cm-1
C      EXTPOT=EXTPOT*1.4388*BC ! Convert energy to Tsukanov's units
      RETURN
	END
C--------------------------------------------------------------------
      FUNCTION SUMLEG(A,N,CX)
	IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N)
C Sums the weigted Legendre polynomials A(I)*Pn(x), x=cos(gamma)
C  n=0,1,2,3,4
        DCX=CX*CX
        DDCX=DCX*DCX
       VA0=1.D0  ! P0(x)
       VA1=CX    ! P1(x)
       VA2=1.5D0*DCX-0.5D0  ! P2(x)
       VA3=0.5D0*CX*(5.D0*DCX-3.D0)  ! P3(x)
       VA4=4.375D0*DDCX-3.75D0*DCX+0.375D0   ! P4(x)
	SUMLEG=A(1)*VA0+A(2)*VA1+A(3)*VA2+A(4)*VA3+A(5)*VA4
	RETURN
	END