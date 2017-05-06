C edited 02.2005:
C   potn2n2 -> function;
C   external data files removed (assignment blocks added);
C   special comments for MS compiler/linker added (to build a .DLL)
C ==================================================================
C
C PES AVDA ET AL JCP84(3), 1629 (1986)
C CE PGM CALCULE CETTE SURFACE V(R,THETA1,THETA2,PHI)
C COEF LUS (TABLE I ET II, JCP) EN KJ.MOL-1 ET R EN NM
C FT 21/02/03
C
C     en sortie R en ua et V en cm-1
C
c  a compiler avec -L /u/libPOWER2  -lbibli -l slatec
c
c ==================================================================
c  finally:  R - atomic units,  V - cm-1,  three angles - radians
c ==================================================================
c
      subroutine potinit
	  !MS$ATTRIBUTES DLLEXPORT :: potinit
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MAT/MATIDX(20,3), COEFOVP(20,3)
      common/CN/CNDISP(5,3)
c -----------------------------------------------------
c
c      OPEN(UNIT=11,FILE='coex_n4.dat',STATUS='OLD')
c      DO LIGN=1,20
c         READ(11,*) (MATIDX(LIGN,KOL),KOL=1,3),
c     &              (COEFOVP(LIGN,KOL),KOL=1,3)
c      ENDDO
c      CLOSE(11)
c
c      OPEN(UNIT=11,FILE='cn_n4.dat')
c      DO LIGN=1,5
c	 READ(11,*) L1, L2, L3, (CNDISP(LIGN,KOL),KOL=1,3)
c      ENDDO
c      CLOSE(11)
c
 	MATIDX(1,1)=0
	MATIDX(1,2)=0
	MATIDX(1,3)=0
	COEFOVP(1,1)=4.26386D5
	COEFOVP(1,2)=26.4733
	COEFOVP(1,3)=12.9778
c
      MATIDX(2,1)=2
	MATIDX(2,2)=0
	MATIDX(2,3)=2
	COEFOVP(2,1)=1.97571D5
	COEFOVP(2,2)=25.7967
	COEFOVP(2,3)=13.8556
c
	MATIDX(3,1)=2
	MATIDX(3,2)=2
	MATIDX(3,3)=0
	COEFOVP(3,1)=2.80720D3
	COEFOVP(3,2)=15.1933
	COEFOVP(3,3)=25.40
c
	MATIDX(4,1)=2
	MATIDX(4,2)=2
	MATIDX(4,3)=2
	COEFOVP(4,1)=-1.03825D4
	COEFOVP(4,2)=19.2333
	COEFOVP(4,3)=21.1333
c
      MATIDX(5,1)=2
	MATIDX(5,2)=2
	MATIDX(5,3)=4
	COEFOVP(5,1)=1.84021D5
	COEFOVP(5,2)=28.15
	COEFOVP(5,3)=11.8556
c
	MATIDX(6,1)=4
	MATIDX(6,2)=0
	MATIDX(6,3)=4
	COEFOVP(6,1)=4.24274D4
	COEFOVP(6,2)=27.2733
	COEFOVP(6,3)=11.3444
c
	MATIDX(7,1)=4
	MATIDX(7,2)=2
	MATIDX(7,3)=2
	COEFOVP(7,1)=2.47313D1
	COEFOVP(7,2)=2.7533
	COEFOVP(7,3)=38.9556
c
	MATIDX(8,1)=4
	MATIDX(8,2)=2
	MATIDX(8,3)=4
	COEFOVP(8,1)=-1.10124D3
	COEFOVP(8,2)=17.67
	COEFOVP(8,3)=22.30
c
	MATIDX(9,1)=4
	MATIDX(9,2)=2
	MATIDX(9,3)=6
	COEFOVP(9,1)=1.00149D5
	COEFOVP(9,2)=32.5267
	COEFOVP(9,3)=6.70
c
	MATIDX(10,1)=4
	MATIDX(10,2)=4
	MATIDX(10,3)=0
	COEFOVP(10,1)=2.80812D-6
	COEFOVP(10,2)=-60.3767
	COEFOVP(10,3)=116.0667
c
	MATIDX(11,1)=4
	MATIDX(11,2)=4
	MATIDX(11,3)=2
	COEFOVP(11,1)=-1.70878D-7
	COEFOVP(11,2)=-78.040
	COEFOVP(11,3)=140.6444
c
	MATIDX(12,1)=4
	MATIDX(12,2)=4
	MATIDX(12,3)=4
	COEFOVP(12,1)=1.54777D-1
	COEFOVP(12,2)=-11.490
	COEFOVP(12,3)=56.1667
c
	MATIDX(13,1)=4
	MATIDX(13,2)=4
	MATIDX(13,3)=6
	COEFOVP(13,1)=-1.61499D2
	COEFOVP(13,2)=16.2467
	COEFOVP(13,3)=24.3444
c
	MATIDX(14,1)=4
	MATIDX(14,2)=4
	MATIDX(14,3)=8
	COEFOVP(14,1)=1.55978D5
	COEFOVP(14,2)=40.150
	COEFOVP(14,3)=-1.3333
c
	MATIDX(15,1)=6
	MATIDX(15,2)=0
	MATIDX(15,3)=6
	COEFOVP(15,1)=1.97922D4
	COEFOVP(15,2)=36.270
	COEFOVP(15,3)=-0.2556
c
	MATIDX(16,1)=6
	MATIDX(16,2)=2
	MATIDX(16,3)=4
	COEFOVP(16,1)=7.70923D-1
	COEFOVP(16,2)=1.270
	COEFOVP(16,3)=38.8222
c
	MATIDX(17,1)=6
	MATIDX(17,2)=2
	MATIDX(17,3)=6
	COEFOVP(17,1)=-4.83277D2
	COEFOVP(17,2)=26.980
	COEFOVP(17,3)=10.000
c
	MATIDX(18,1)=6
	MATIDX(18,2)=2
	MATIDX(18,3)=8
	COEFOVP(18,1)=1.13852D5
	COEFOVP(18,2)=44.0067
	COEFOVP(18,3)=-7.0333
c
	MATIDX(19,1)=6
	MATIDX(19,2)=4
	MATIDX(19,3)=10
	COEFOVP(19,1)=4.78423D5
	COEFOVP(19,2)=54.5467
	COEFOVP(19,3)=-17.5778
c
	MATIDX(20,1)=6
	MATIDX(20,2)=6
	MATIDX(20,3)=12
	COEFOVP(20,1)=3.97250D6
	COEFOVP(20,2)=71.7233
	COEFOVP(20,3)=-35.7889
c     -1 0 0
c
c	 0 0 0
	CNDISP(1,1)=4.119D-3
	CNDISP(1,2)=3.795D-4
	CNDISP(1,3)=3.491D-5
c      2 0 2
	CNDISP(2,1)=2.165D-4
	CNDISP(2,2)=9.780D-5
	CNDISP(2,3)=1.38D-5
c      2 2 0
	CNDISP(3,1)=5.454D-6
	CNDISP(3,2)=1.458D-6
	CNDISP(3,3)=3.150D-7
c      2 2 2
	CNDISP(4,1)=6.518D-6
	CNDISP(4,2)=-2.179D-6
	CNDISP(4,3)=-8.662D-7
c
	L1=2
	L2=2
	L3=4
	CNDISP(5,1)=5.247D-5
	CNDISP(5,2)=1.133D-5
	CNDISP(5,3)=5.778D-6
c     -1 0 0
c -----------------------------------------------------
	return
	end
c =====================================================
c
      real*8 function potn2n2(rr,theta1,theta2,phi)
	  !MS$ATTRIBUTES DLLEXPORT :: potn2n2
c
      implicit real*8(a-h,o-z)
      LOGICAL FIRST
      CHARACTER*1 REP
      DIMENSION C_ELEC(6), IDXC_ELE(12)
      COMMON /PREMIER/ FIRST
      common/CN/CNDISP(5,3)
      COMMON /RPUI/ XPUI(0:10)
      COMMON/MAT/MATIDX(20,3), COEFOVP(20,3)
      DATA PI /3.141592653589793D0/, BOHR /.5291772083D0/,
     &         ECONV /1.196264D-2/, XSCALE /1.5D0/ , ZSCALE /1.0509D0/
      DATA C_ELEC /1.624D-3, 6.771D-5, 5.411D-6, 8.024D-7, 1.046D-7,
     &                2.996D-9/
      DATA IDXC_ELE /2,2,4,2,4,4,6,2,6,4,6,6/
c
      XNORM=dsqrt(4D0*PI)
      XNORM=XNORM**3
      RCUT=14D0                 ! EN UA
c
      rep='o'
c
c     THETA1=PI*THETA1/180d0
c     THETA2=PI*THETA2/180d0
c     PHI   =PI*PHI/180d0
c
         R=RR*BOHR/10D0		    ! BOHR --> ANGSTROM --> NM
         CALL PUISS(R)
c
         VPOT=0D0
         VELE=0D0
         VOVL=0D0
         VDIS=0D0
c
      DO ITERM=1,20
         FIRST=.TRUE.
         LA=MATIDX(ITERM,1)
         LB=MATIDX(ITERM,2)
         LL=MATIDX(ITERM,3)
c
         ALALBL=PR3JYY(LA,LB,LL,THETA1,THETA2,PHI)
cc
*        write(12,*) 'iterm, alalbl', iterm, alalbl
cc
         FLALBL=COEFOVP(ITERM,1)*XSCALE
         ALPHA =COEFOVP(ITERM,2)*ZSCALE
         BETA  =COEFOVP(ITERM,3)*ZSCALE*ZSCALE
         TBETA=BETA
c
C ON CHERCHE L'INDICE POUR LA CONTRIBUTION ELECT
         VELECT=0D0
	 ICHECK=0
         IF (ITERM.LT.5) GOTO 30
	 LVER=LA+LB
	 IF (LVER.EQ.LL) THEN
	    DO IC=1,6
	       LPA=IDXC_ELE(2*IC-1)
	       LPB=IDXC_ELE(2*IC)
	       IF (LPA.EQ.LA.AND.LPB.EQ.LB) ICHECK=IC
	    ENDDO
	 ENDIF
c
C CONTRIBUTION ELECTROSTATIQUE
         IF (ICHECK.NE.0) THEN
            IEXPO =-(LA+LB+1)
            VELECT=C_ELEC(ICHECK)*R**IEXPO
         ENDIF
c
C CONTRIBUTION D'ECHANGE ... (OVERLAP)
 30      VOVERL=0D0
         R2=XPUI(2)
         IF (BETA.LT.0D0.AND.RR.GT.RCUT) tbeta=0.
C ceci affecte vdis via fdamp de maniere tres tres legere
         ARGUM =-(ALPHA*R+TBETA*R2)
         FDEXP =DEXP(ARGUM)
         IF (BETA.LT.0D0.AND.RR.GT.RCUT) GOTO 40
         VOVERL=FLALBL*FDEXP

C CONTRIBUTION DISPERSION
 40      VDISP=0D0
         IF (ITERM.LE.5) THEN
            R6 =XPUI(6)
            R8 =XPUI(8)
            R10=XPUI(10)
            f6 =1d0
            f8 =1d0
            f10=1d0
            if (rep.eq.'n') goto 45
            F6 =FDAMP(6,ALPHA,TBETA,FDEXP)
            F8 =FDAMP(8,ALPHA,TBETA,FDEXP)
            F10=FDAMP(10,ALPHA,TBETA,FDEXP)
 45         VD6 =F6*CNDISP(ITERM,1)/R6
            VD8 =F8*CNDISP(ITERM,2)/R8
            VD10=F10*CNDISP(ITERM,3)/R10
            VDISP=-(VD6+VD8+VD10)
         ENDIF

C POTENTIEL TOTAL (V_LA,LB,L(R)) "RADIAL"
 50      VELECT=VELECT/ECONV
         VDISP =VDISP/ECONV
         VOVERL=VOVERL/ECONV
         DPOT=VELECT+VOVERL+VDISP                    ! EN CM-1
c
C POTENTIEL V(R,T1,T2,PHI)
         VELE=VELE+ALALBL*VELECT*XNORM
         VDIS=VDIS+ALALBL*VDISP*XNORM
         VOVL=VOVL+ALALBL*VOVERL*XNORM
         VPOT=VPOT+ALALBL*DPOT*XNORM

      ENDDO			! ITERM
	potn2n2=VPOT
      return
      END
****************************************************
C
      FUNCTION FDAMP(N,CTE1,CTE2,CTE3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NFCTMX=1001)
      LOGICAL FIRST
      DIMENSION AMAT1(0:10), AMAT2(0:10)
      COMMON /LOGFAC/FCT(NFCTMX)
      COMMON /RPUI/ XPUI(0:10)
      COMMON /PREMIER/ FIRST
      DATA NTIMES /0/, EPS /0.01/
      SAVE AMAT1, SUM, NTIMES

 5    IF (NTIMES .EQ.0) CALL FACLOG
      NTIMES = NTIMES+1

      IF (FIRST) THEN
         DO K=0,10
            AMAT1(K)=0D0
            KMIN=(K+1)/2
            DO I=KMIN,K
               AMAT1(K)=AMAT1(K)+(CTE1**(2*I-K))*(CTE2**(K-I))
     &                           /DEXP(FCT(2*I-K+1)+FCT(K-I+1))
            ENDDO               ! I
         ENDDO			! K
         FIRST=.FALSE.
      ENDIF

      IF (N.GT.6) GOTO 30

      DO K=0,10
         AMAT2(K)=AMAT1(K)*XPUI(K)
      ENDDO

      SUM=0.
      DO K=0,6
         SUM=SUM+AMAT2(K)
      ENDDO
      GOTO 40

 30   IF (N.GT.8) GOTO 35

      DO K=7,8
	 SUM=SUM+AMAT2(K)
      ENDDO
      GOTO 40

 35    DO K=9,10
	 SUM=SUM+AMAT2(K)
      ENDDO

 40   FDAMP=1D0-SUM*CTE3

      RETURN
      END
****************************************************
C
      SUBROUTINE PUISS(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /RPUI/ XPUI(0:10)

      XPUI(0)=1D0
      DO K=1,10
         XPUI(K)=XPUI(K-1)*X
      ENDDO
      RETURN
      END
****************************************************
C
      FUNCTION PR3JYY(L1,L2,L3,THTARD1,THTARD2,PHIRD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MMX=6,MDIM=MMX+1)
      LOGICAL PERMUT1, PERMUT2
      DIMENSION SPLG1(MDIM), SPLG2(MDIM), ISLG(MDIM), fact3jm(0:MMX)
      DATA ID /4/, LDIFF /0/, PI /3.141592653589793D0/
      DATA EPS /1D-5/

      PERMUT1=.FALSE.
      PERMUT2=.FALSE.

      THT1=THTARD1
      THT2=THTARD2

      DL1=DFLOAT(L1)
      DL2=DFLOAT(L2)
      DL3=DFLOAT(L3)

      PR3JYY=0.
      NCAS=2                    ! L1.NE.L2
      IF (L1.EQ.L2) NCAS=1      ! donc pas besoin de delta12
      PI_2=PI/2D0
      FNORM=DSQRT(2D0*PI)
      FACT =DSQRT(DL3+0.5D0)

      MMAX=MIN0(L1,L2)
      MMIN=0
      MDIFF=MMAX-MMIN+1 ! *** A VERIFIER ***
      IF (MDIFF.GT.MDIM) STOP '??? AUGMENTER DIMENSIONS ???'

      DO M=MMIN,MMAX
         DM=DFLOAT(M)
         FACT3JM(M)=F3J(DL1,DL2,DL3,DM,-DM,0D0)
      ENDDO
C
C     APPEL A UNE SUBROUTINE DE SLATEC
C     SUBROUTINE DXLEGF (DNU1, NUDIFF, MU1, MU2, THETA, ID, PQA, IPQA,
C    1   IERROR)
C     ID=4 ASSOCIATED NORMALIZED LEGENDRE POLYNOMIALS
C     MU==M >=0, NU==L >=0
C    THETA is DOUBLE PRECISION and in the half-open interval (0,PI/2];
C
      IF (THT1.LT.0D0.OR.THT2.LT.0D0) STOP 'STOP1 SLATEC'
      SIGN1=1.
      IF (THT1.GT.PI_2) THEN
         THT1=PI-THT1
         SIGN1=PPARITYY (L1)       ! tj >0 ici
         PERMUT1=.TRUE.
      ENDIF
      SIGN2=1.
      IF (THT2.GT.PI_2) THEN
         THT2=PI-THT2
         SIGN2=PPARITYY (L2)       ! tj >0 ici
         PERMUT2=.TRUE.
      ENDIF
      IF (THT1.LT.EPS) THT1=EPS
      IF (THT2.LT.EPS) THT2=EPS

      ISLG=0                ! ATTENTION MATRICE
      IERROR=0

      sign11=sign1
      sign22=sign2

      TH1=THT1
      TH2=THT2
      DO ICAS=1,NCAS
         IF (ICAS.EQ.2) THEN
            TH1=THT2
            TH2=THT1
         ENDIF

         SPLG1=0.            ! ATTENTION MATRICE
         SPLG2=0.            ! ATTENTION MATRICE

         CALL DXLEGF (DL1, LDIFF, MMIN, MMAX, TH1, ID, SPLG1, ISLG,
     1                IERROR)
         IF (IERROR.NE.0) STOP 'IERROR'
         DO II=1,MDIFF
            IF (ISLG(II).NE.0) STOP 'ISLG'
         ENDDO

         CALL DXLEGF (DL2, LDIFF, MMIN, MMAX, TH2, ID, SPLG2, ISLG,
     2                IERROR)
         IF (IERROR.NE.0) STOP 'IERROR'
         DO II=1,MDIFF
            IF (ISLG(II).NE.0) STOP 'ISLG'
         ENDDO

         SUM=SIGN1*SPLG1(1)*SIGN2*SPLG2(1)*FACT3JM(0)
         DO M=MMIN+1,MMAX
            sign0=PPARITYY (m)
            if (permut1) sign11=sign1*sign0
c                           		      ! car Plm(-x)=(-1)**(l+m)*Plm(x)
            if (permut2) sign22=sign2*sign0             ! idem
            DM=DFLOAT(M)
            TMP=2D0*sign0*FACT3JM(M)
     &          *SIGN11*SPLG1(M+1)*SIGN22*SPLG2(M+1)*DCOS(DM*PHIRD)
            SUM=SUM+TMP
         ENDDO
         PR3JYY=PR3JYY+SUM
      ENDDO                     ! ICAS
      PR3JYY=PR3JYY*FACT/(FNORM**3)
      RETURN
      END
****************************************************
C
      FUNCTION PPARITYY (I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PPARITYY =1.D0
      IF((I/2)*2-I.NE.0) PPARITYY =-1.D0
      RETURN
      END
****************************************************
c
      double precision function f3j (fj1,fj2,fj3, fm1,fm2,fm3)
c#
c#    calculates 3j coefficients from racah formula
c#    (messiah: t2, p 910; formula 21) .
c#    clebsch-gordan coefficients are given by (p. 908, formula 12) :
c#                         j -j +m                |j    j     j|
c#    <j j m m |j m> = (-1) 1  2   (2*j+1)**(0.5) | 1    2     |
c#      1 2 1 2                                   |m    m    -m|
c#                                                | 1    2     |
c#
c#    has been tested for j up to 200.
c#    logfac should contain the logarithms of factorials
c#    fj-fm integer not checked
c#    j.m.l. (1975)
c#
      implicit double precision (a-h,o-z)
      integer t,tmin,tmax
      parameter (nfctmx=1001)
      data tiny,zero,one /0.01d0,0.d0,1.d0/ ,ntimes /1/
      common /logfac/ fct(nfctmx)
      if (ntimes .eq. 1) call faclog
      ntimes = ntimes+1
      cc = zero
      if (fj3 . gt. (fj1+fj2+tiny))      go to 100
      if (dabs(fj1-fj2) .gt. (fj3+tiny)) go to 100
      if (dabs(fm1+fm2+fm3) .gt. tiny)   go to 100
      if (dabs(fm1) .gt. (fj1+tiny))     go to 100
      if (dabs(fm2) .gt. (fj2+tiny))     go to 100
      if (dabs(fm3) .gt. (fj3+tiny))     go to 100
      fk1 = fj3-fj2+fm1
      fk2 = fj3-fj1-fm2
      fk3 = fj1-fm1
      fk4 = fj2+fm2
      fk5 = fj1+fj2-fj3
      fk1m = fk1-tiny
      fk2m = fk2-tiny
      fk1p = fk1+tiny
      fk2p = fk2+tiny
      if (fk1m .lt. zero) k1 = fk1m
      if (fk1p .gt. zero) k1 = fk1p
      if (fk2m .lt. zero) k2 = fk2m
      if (fk2p .gt. zero) k2 = fk2p
      k3 = fk3+tiny
      k4 = fk4+tiny
      k5 = fk5+tiny
      tmin = 0
      if (k1+tmin .lt. 0) tmin = -k1
      if (k2+tmin .lt. 0) tmin = -k2
      tmax = k3
      if (k4-tmax .lt. 0) tmax = k4
      if (k5-tmax .lt. 0) tmax = k5
      n1 = fj1+fj2-fj3+one+tiny
      n2 = fj2+fj3-fj1+one+tiny
      n3 = fj3+fj1-fj2+one+tiny
      n4 = fj1+fm1+one+tiny
      n5 = fj2+fm2+one+tiny
      n6 = fj3+fm3+one+tiny
      n7 = fj1-fm1+one+tiny
      n8 = fj2-fm2+one+tiny
      n9 = fj3-fm3+one+tiny
      n10 = fj1+fj2+fj3+2.d0+tiny
      x = fct(n1)+fct(n2)+fct(n3)+fct(n4)+fct(n5)+fct(n6)
     &   +fct(n7)+fct(n8)+fct(n9)-fct(n10)
      x = 0.5d0*x
      do 10  t = tmin,tmax
         phase = one
         if (mod(t,2) .ne. 0) phase = -one
         cc = cc+phase*dexp(-fct(t+1)   -fct(k1+t+1)-fct(k2+t+1)
     &                      -fct(k3-t+1)-fct(k4-t+1)-fct(k5-t+1)+x)
 10   continue
      fsp = dabs(fj1-fj2-fm3)+tiny
      ns = fsp
      if (mod(ns,2) .gt. 0) cc = -cc
 100  f3j = cc
      return
      end
C-----------------------------------------------------------------
C
      subroutine faclog
c#    initialisation of logarithms of factorials array
      implicit double precision (a-h,o-z)
      parameter (nfctmx=1001)
      common /logfac/ fct(nfctmx)
      data ntimes /0/
c
      ntimes = ntimes+1
      if (ntimes .gt. 1) return
      fct(1) = 0.d0
      do 10 i = 1,nfctmx-1
         ai = i
         fct(i+1) = fct(i)+dlog(ai)
 10   continue
c
      return
      end
C======================================================================
C-------------------------SLATEC-PROCEDURES----------------------------
C======================================================================
C
      SUBROUTINE DXLEGF (DNU1, NUDIFF, MU1, MU2, THETA, ID, PQA, IPQA,
     1   IERROR)
      DOUBLE PRECISION PQA,DNU1,DNU2,SX,THETA,X,PI2
      DIMENSION PQA(*),IPQA(*)
C
C***FIRST EXECUTABLE STATEMENT  DXLEGF
      IERROR=0
      CALL DXSET (0, 0, 0.0D0, 0,IERROR)
      IF (IERROR.NE.0) RETURN
      PI2=2.D0*ATAN(1.D0)
C
C        ZERO OUTPUT ARRAYS
C
      L=(MU2-MU1)+NUDIFF+1
      DO 290 I=1,L
      PQA(I)=0.D0
  290 IPQA(I)=0
C
C        CHECK FOR VALID INPUT VALUES
C
      IF(NUDIFF.LT.0) GO TO 400
      IF(DNU1.LT.-.5D0) GO TO 400
      IF(MU2.LT.MU1) GO TO 400
      IF(MU1.LT.0) GO TO 400
      IF(THETA.LE.0.D0.OR.THETA.GT.PI2) GO TO 420
      IF(ID.LT.1.OR.ID.GT.4) GO TO 400
      IF((MU1.NE.MU2).AND.(NUDIFF.GT.0)) GO TO 400
C
C        IF DNU1 IS NOT AN INTEGER, NORMALIZED P(MU,DNU,X)
C        CANNOT BE CALCULATED.  IF DNU1 IS AN INTEGER AND
C        MU1.GT.DNU2 THEN ALL VALUES OF P(+MU,DNU,X) AND
C        NORMALIZED P(MU,NU,X) WILL BE ZERO.
C
      DNU2=DNU1+NUDIFF
      IF((ID.EQ.3).AND.(MOD(DNU1,1.D0).NE.0.D0)) GO TO 295
      IF((ID.EQ.4).AND.(MOD(DNU1,1.D0).NE.0.D0)) GO TO 400
      IF((ID.EQ.3.OR.ID.EQ.4).AND.MU1.GT.DNU2) RETURN
  295 CONTINUE
C
      X=COS(THETA)
      SX=1.D0/SIN(THETA)
      IF(ID.EQ.2) GO TO 300
      IF(MU2-MU1.LE.0) GO TO 360
C
C        FIXED NU, VARIABLE MU
C        CALL DXPMU TO CALCULATE P(-MU1,NU,X),....,P(-MU2,NU,X)
C
      CALL DXPMU(DNU1,DNU2,MU1,MU2,THETA,X,SX,ID,PQA,IPQA,IERROR)
      IF (IERROR.NE.0) RETURN
      GO TO 380
C
  300 IF(MU2.EQ.MU1) GO TO 320
C
C        FIXED NU, VARIABLE MU
C        CALL DXQMU TO CALCULATE Q(MU1,NU,X),....,Q(MU2,NU,X)
C
      CALL DXQMU(DNU1,DNU2,MU1,MU2,THETA,X,SX,ID,PQA,IPQA,IERROR)
      IF (IERROR.NE.0) RETURN
      GO TO 390
C
C        FIXED MU, VARIABLE NU
C        CALL DXQNU TO CALCULATE Q(MU,DNU1,X),....,Q(MU,DNU2,X)
C
  320 CALL DXQNU(DNU1,DNU2,MU1,THETA,X,SX,ID,PQA,IPQA,IERROR)
      IF (IERROR.NE.0) RETURN
      GO TO 390
C
C        FIXED MU, VARIABLE NU
C        CALL DXPQNU TO CALCULATE P(-MU,DNU1,X),....,P(-MU,DNU2,X)
C
  360 CALL DXPQNU(DNU1,DNU2,MU1,THETA,ID,PQA,IPQA,IERROR)
      IF (IERROR.NE.0) RETURN
C
C        IF ID = 3, TRANSFORM P(-MU,NU,X) VECTOR INTO
C        P(MU,NU,X) VECTOR.
C
  380 IF(ID.EQ.3) CALL DXPMUP(DNU1,DNU2,MU1,MU2,PQA,IPQA,IERROR)
      IF (IERROR.NE.0) RETURN
C
C        IF ID = 4, TRANSFORM P(-MU,NU,X) VECTOR INTO
C        NORMALIZED P(MU,NU,X) VECTOR.
C
      IF(ID.EQ.4) CALL DXPNRM(DNU1,DNU2,MU1,MU2,PQA,IPQA,IERROR)
      IF (IERROR.NE.0) RETURN
C
C        PLACE RESULTS IN REDUCED FORM IF POSSIBLE
C        AND RETURN TO MAIN PROGRAM.
C
  390 DO 395 I=1,L
      CALL DXRED(PQA(I),IPQA(I),IERROR)
      IF (IERROR.NE.0) RETURN
  395 CONTINUE
      RETURN
C
C        *****     ERROR TERMINATION     *****
C
  400 CALL XERMSG ('SLATEC', 'DXLEGF',
     +             'DNU1, NUDIFF, MU1, MU2, or ID not valid', 210, 1)
      IERROR=210
      RETURN
  420 CALL XERMSG ('SLATEC', 'DXLEGF', 'THETA out of range', 211, 1)
      IERROR=211
      RETURN
      END
C----------------------------------------------------------------------
C
      SUBROUTINE DXPMU (NU1, NU2, MU1, MU2, THETA, X, SX, ID, PQA, IPQA,
     1   IERROR)
      DOUBLE PRECISION PQA,NU1,NU2,P0,X,SX,THETA,X1,X2
      DIMENSION PQA(*),IPQA(*)
C
C        CALL DXPQNU TO OBTAIN P(-MU2,NU,X)
C
C***FIRST EXECUTABLE STATEMENT  DXPMU
      IERROR=0
      CALL DXPQNU(NU1,NU2,MU2,THETA,ID,PQA,IPQA,IERROR)
      IF (IERROR.NE.0) RETURN
      P0=PQA(1)
      IP0=IPQA(1)
      MU=MU2-1
C
C        CALL DXPQNU TO OBTAIN P(-MU2-1,NU,X)
C
      CALL DXPQNU(NU1,NU2,MU,THETA,ID,PQA,IPQA,IERROR)
      IF (IERROR.NE.0) RETURN
      N=MU2-MU1+1
      PQA(N)=P0
      IPQA(N)=IP0
      IF(N.EQ.1) GO TO 300
      PQA(N-1)=PQA(1)
      IPQA(N-1)=IPQA(1)
      IF(N.EQ.2) GO TO 300
      J=N-2
  290 CONTINUE
C
C        BACKWARD RECURRENCE IN MU TO OBTAIN
C              P(-MU2,NU1,X),P(-(MU2-1),NU1,X),....P(-MU1,NU1,X)
C              USING
C              (NU-MU)*(NU+MU+1.)*P(-(MU+1),NU,X)=
C                2.*MU*X*SQRT((1./(1.-X**2))*P(-MU,NU,X)-P(-(MU-1),NU,X)
C
      X1=2.D0*MU*X*SX*PQA(J+1)
      X2=-(NU1-MU)*(NU1+MU+1.D0)*PQA(J+2)
      CALL DXADD(X1,IPQA(J+1),X2,IPQA(J+2),PQA(J),IPQA(J),IERROR)
      IF (IERROR.NE.0) RETURN
      CALL DXADJ(PQA(J),IPQA(J),IERROR)
      IF (IERROR.NE.0) RETURN
      IF(J.EQ.1) GO TO 300
      J=J-1
      MU=MU-1
      GO TO 290
  300 RETURN
      END
C----------------------------------------------------------------------
C
      SUBROUTINE DXPMUP (NU1, NU2, MU1, MU2, PQA, IPQA, IERROR)
      DOUBLE PRECISION DMU,NU,NU1,NU2,PQA,PROD
      DIMENSION PQA(*),IPQA(*)
C***FIRST EXECUTABLE STATEMENT  DXPMUP
      IERROR=0
      NU=NU1
      MU=MU1
      DMU=MU
      N=INT(NU2-NU1+.1D0)+(MU2-MU1)+1
      J=1
      IF(MOD(REAL(NU),1.).NE.0.) GO TO 210
  200 IF(DMU.LT.NU+1.D0) GO TO 210
      PQA(J)=0.D0
      IPQA(J)=0
      J=J+1
      IF(J.GT.N) RETURN
C        INCREMENT EITHER MU OR NU AS APPROPRIATE.
      IF(NU2-NU1.GT..5D0) NU=NU+1.D0
      IF(MU2.GT.MU1) MU=MU+1
      GO TO 200
C
C        TRANSFORM P(-MU,NU,X) TO P(MU,NU,X) USING
C        P(MU,NU,X)=(NU-MU+1)*(NU-MU+2)*...*(NU+MU)*P(-MU,NU,X)*(-1)**MU
C
  210 PROD=1.D0
      IPROD=0
      K=2*MU
      IF(K.EQ.0) GO TO 222
      DO 220 L=1,K
      PROD=PROD*(DMU-NU-L)
  220 CALL DXADJ(PROD,IPROD,IERROR)
      IF (IERROR.NE.0) RETURN
  222 CONTINUE
      DO 240 I=J,N
      IF(MU.EQ.0) GO TO 225
      PQA(I)=PQA(I)*PROD*(-1)**MU
      IPQA(I)=IPQA(I)+IPROD
      CALL DXADJ(PQA(I),IPQA(I),IERROR)
      IF (IERROR.NE.0) RETURN
  225 IF(NU2-NU1.GT..5D0) GO TO 230
      PROD=(DMU-NU)*PROD*(-DMU-NU-1.D0)
      CALL DXADJ(PROD,IPROD,IERROR)
      IF (IERROR.NE.0) RETURN
      MU=MU+1
      DMU=DMU+1.D0
      GO TO 240
  230 PROD=PROD*(-DMU-NU-1.D0)/(DMU-NU-1.D0)
      CALL DXADJ(PROD,IPROD,IERROR)
      IF (IERROR.NE.0) RETURN
      NU=NU+1.D0
  240 CONTINUE
      RETURN
      END
C----------------------------------------------------------------------
C
      SUBROUTINE DXPNRM (NU1, NU2, MU1, MU2, PQA, IPQA, IERROR)
      DOUBLE PRECISION C1,DMU,NU,NU1,NU2,PQA,PROD
      DIMENSION PQA(*),IPQA(*)
C***FIRST EXECUTABLE STATEMENT  DXPNRM
      IERROR=0
      L=(MU2-MU1)+(NU2-NU1+1.5D0)
      MU=MU1
      DMU=MU1
      NU=NU1
C
C         IF MU .GT.NU, NORM P =0.
C
      J=1
  500 IF(DMU.LE.NU) GO TO 505
      PQA(J)=0.D0
      IPQA(J)=0
      J=J+1
      IF(J.GT.L) RETURN
C
C        INCREMENT EITHER MU OR NU AS APPROPRIATE.
C
      IF(MU2.GT.MU1) DMU=DMU+1.D0
      IF(NU2-NU1.GT..5D0) NU=NU+1.D0
      GO TO 500
C
C         TRANSFORM P(-MU,NU,X) INTO NORMALIZED P(MU,NU,X) USING
C              NORM P(MU,NU,X)=
C                 SQRT((NU+.5)*FACTORIAL(NU+MU)/FACTORIAL(NU-MU))
C                              *P(-MU,NU,X)
C
  505 PROD=1.D0
      IPROD=0
      K=2*MU
      IF(K.LE.0) GO TO 520
      DO 510 I=1,K
      PROD=PROD*SQRT(NU+DMU+1.D0-I)
  510 CALL DXADJ(PROD,IPROD,IERROR)
      IF (IERROR.NE.0) RETURN
  520 DO 540 I=J,L
      C1=PROD*SQRT(NU+.5D0)
      PQA(I)=PQA(I)*C1
      IPQA(I)=IPQA(I)+IPROD
      CALL DXADJ(PQA(I),IPQA(I),IERROR)
      IF (IERROR.NE.0) RETURN
      IF(NU2-NU1.GT..5D0) GO TO 530
      IF(DMU.GE.NU) GO TO 525
      PROD=SQRT(NU+DMU+1.D0)*PROD
      IF(NU.GT.DMU) PROD=PROD*SQRT(NU-DMU)
      CALL DXADJ(PROD,IPROD,IERROR)
      IF (IERROR.NE.0) RETURN
      MU=MU+1
      DMU=DMU+1.D0
      GO TO 540
  525 PROD=0.D0
      IPROD=0
      MU=MU+1
      DMU=DMU+1.D0
      GO TO 540
  530 PROD=SQRT(NU+DMU+1.D0)*PROD
      IF(NU.NE.DMU-1.D0) PROD=PROD/SQRT(NU-DMU+1.D0)
      CALL DXADJ(PROD,IPROD,IERROR)
      IF (IERROR.NE.0) RETURN
      NU=NU+1.D0
  540 CONTINUE
      RETURN
      END
C----------------------------------------------------------------------
C
      SUBROUTINE DXPQNU (NU1, NU2, MU, THETA, ID, PQA, IPQA, IERROR)
      DOUBLE PRECISION A,NU,NU1,NU2,PQ,PQA,DXPSI,R,THETA,W,X,X1,X2,XS,
     1 Y,Z
      DOUBLE PRECISION DI,DMU,PQ1,PQ2,FACTMU,FLOK
      DIMENSION PQA(*),IPQA(*)
      COMMON /DXBLK1/ NBITSF
      SAVE /DXBLK1/
C
C        J0, IPSIK, AND IPSIX ARE INITIALIZED IN THIS SUBROUTINE.
C        J0 IS THE NUMBER OF TERMS USED IN SERIES EXPANSION
C        IN SUBROUTINE DXPQNU.
C        IPSIK, IPSIX ARE VALUES OF K AND X RESPECTIVELY
C        USED IN THE CALCULATION OF THE DXPSI FUNCTION.
C
C***FIRST EXECUTABLE STATEMENT  DXPQNU
      IERROR=0
      J0=NBITSF
      IPSIK=1+(NBITSF/10)
      IPSIX=5*IPSIK
      IPQ=0
C        FIND NU IN INTERVAL [-.5,.5) IF ID=2  ( CALCULATION OF Q )
      NU=MOD(NU1,1.D0)
      IF(NU.GE..5D0) NU=NU-1.D0
C        FIND NU IN INTERVAL (-1.5,-.5] IF ID=1,3, OR 4  ( CALC. OF P )
      IF(ID.NE.2.AND.NU.GT.-.5D0) NU=NU-1.D0
C        CALCULATE MU FACTORIAL
      K=MU
      DMU=MU
      IF(MU.LE.0) GO TO 60
      FACTMU=1.D0
      IF=0
      DO 50 I=1,K
      FACTMU=FACTMU*I
   50 CALL DXADJ(FACTMU,IF,IERROR)
      IF (IERROR.NE.0) RETURN
   60 IF(K.EQ.0) FACTMU=1.D0
      IF(K.EQ.0) IF=0
C
C        X=COS(THETA)
C        Y=SIN(THETA/2)**2=(1-X)/2=.5-.5*X
C        R=TAN(THETA/2)=SQRT((1-X)/(1+X)
C
      X=COS(THETA)
      Y=SIN(THETA/2.D0)**2
      R=TAN(THETA/2.D0)
C
C        USE ASCENDING SERIES TO CALCULATE TWO VALUES OF P OR Q
C        FOR USE AS STARTING VALUES IN RECURRENCE RELATION.
C
      PQ2=0.0D0
      DO 100 J=1,2
      IPQ1=0
      IF(ID.EQ.2) GO TO 80
C
C        SERIES FOR P ( ID = 1, 3, OR 4 )
C        P(-MU,NU,X)=1./FACTORIAL(MU)*SQRT(((1.-X)/(1.+X))**MU)
C                *SUM(FROM 0 TO J0-1)A(J)*(.5-.5*X)**J
C
      IPQ=0
      PQ=1.D0
      A=1.D0
      IA=0
      DO 65 I=2,J0
      DI=I
      A=A*Y*(DI-2.D0-NU)*(DI-1.D0+NU)/((DI-1.D0+DMU)*(DI-1.D0))
      CALL DXADJ(A,IA,IERROR)
      IF (IERROR.NE.0) RETURN
      IF(A.EQ.0.D0) GO TO 66
      CALL DXADD(PQ,IPQ,A,IA,PQ,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
   65 CONTINUE
   66 CONTINUE
      IF(MU.LE.0) GO TO 90
      X2=R
      X1=PQ
      K=MU
      DO 77 I=1,K
      X1=X1*X2
   77 CALL DXADJ(X1,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
      PQ=X1/FACTMU
      IPQ=IPQ-IF
      CALL DXADJ(PQ,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
      GO TO 90
C
C        Z=-LN(R)=.5*LN((1+X)/(1-X))
C
   80 Z=-LOG(R)
      W=DXPSI(NU+1.D0,IPSIK,IPSIX)
      XS=1.D0/SIN(THETA)
C
C        SERIES SUMMATION FOR Q ( ID = 2 )
C        Q(0,NU,X)=SUM(FROM 0 TO J0-1)((.5*LN((1+X)/(1-X))
C    +DXPSI(J+1,IPSIK,IPSIX)-DXPSI(NU+1,IPSIK,IPSIX)))*A(J)*(.5-.5*X)**J
C
C        Q(1,NU,X)=-SQRT(1./(1.-X**2))+SQRT((1-X)/(1+X))
C             *SUM(FROM 0 T0 J0-1)(-NU*(NU+1)/2*LN((1+X)/(1-X))
C                 +(J-NU)*(J+NU+1)/(2*(J+1))+NU*(NU+1)*
C     (DXPSI(NU+1,IPSIK,IPSIX)-DXPSI(J+1,IPSIK,IPSIX))*A(J)*(.5-.5*X)**J
C
C        NOTE, IN THIS LOOP K=J+1
C
      PQ=0.D0
      IPQ=0
      IA=0
      A=1.D0
      DO 85 K=1,J0
      FLOK=K
      IF(K.EQ.1) GO TO 81
      A=A*Y*(FLOK-2.D0-NU)*(FLOK-1.D0+NU)/((FLOK-1.D0+DMU)*(FLOK-1.D0))
      CALL DXADJ(A,IA,IERROR)
      IF (IERROR.NE.0) RETURN
   81 CONTINUE
      IF(MU.GE.1) GO TO 83
      X1=(DXPSI(FLOK,IPSIK,IPSIX)-W+Z)*A
      IX1=IA
      CALL DXADD(PQ,IPQ,X1,IX1,PQ,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
      GO TO 85
   83 X1=(NU*(NU+1.D0)*(Z-W+DXPSI(FLOK,IPSIK,IPSIX))+(NU-FLOK+1.D0)
     1  *(NU+FLOK)/(2.D0*FLOK))*A
      IX1=IA
      CALL DXADD(PQ,IPQ,X1,IX1,PQ,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
   85 CONTINUE
      IF(MU.GE.1) PQ=-R*PQ
      IXS=0
      IF(MU.GE.1) CALL DXADD(PQ,IPQ,-XS,IXS,PQ,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
      IF(J.EQ.2) MU=-MU
      IF(J.EQ.2) DMU=-DMU
   90 IF(J.EQ.1) PQ2=PQ
      IF(J.EQ.1) IPQ2=IPQ
      NU=NU+1.D0
  100 CONTINUE
      K=0
      IF(NU-1.5D0.LT.NU1) GO TO 120
      K=K+1
      PQA(K)=PQ2
      IPQA(K)=IPQ2
      IF(NU.GT.NU2+.5D0) RETURN
  120 PQ1=PQ
      IPQ1=IPQ
      IF(NU.LT.NU1+.5D0) GO TO 130
      K=K+1
      PQA(K)=PQ
      IPQA(K)=IPQ
      IF(NU.GT.NU2+.5D0) RETURN
C
C        FORWARD NU-WISE RECURRENCE FOR F(MU,NU,X) FOR FIXED MU
C        USING
C        (NU+MU+1)*F(MU,NU,X)=(2.*NU+1)*F(MU,NU,X)-(NU-MU)*F(MU,NU-1,X)
C        WHERE F(MU,NU,X) MAY BE P(-MU,NU,X) OR IF MU IS REPLACED
C        BY -MU THEN F(MU,NU,X) MAY BE Q(MU,NU,X).
C        NOTE, IN THIS LOOP, NU=NU+1
C
  130 X1=(2.D0*NU-1.D0)/(NU+DMU)*X*PQ1
      X2=(NU-1.D0-DMU)/(NU+DMU)*PQ2
      CALL DXADD(X1,IPQ1,-X2,IPQ2,PQ,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
      CALL DXADJ(PQ,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
      NU=NU+1.D0
      PQ2=PQ1
      IPQ2=IPQ1
      GO TO 120
C
      END
C----------------------------------------------------------------------
C
      SUBROUTINE DXQMU (NU1, NU2, MU1, MU2, THETA, X, SX, ID, PQA, IPQA,
     1   IERROR)
      DIMENSION PQA(*),IPQA(*)
      DOUBLE PRECISION DMU,NU,NU1,NU2,PQ,PQA,PQ1,PQ2,SX,X,X1,X2
      DOUBLE PRECISION THETA
C***FIRST EXECUTABLE STATEMENT  DXQMU
      IERROR=0
      MU=0
C
C        CALL DXPQNU TO OBTAIN Q(0.,NU1,X)
C
      CALL DXPQNU(NU1,NU2,MU,THETA,ID,PQA,IPQA,IERROR)
      IF (IERROR.NE.0) RETURN
      PQ2=PQA(1)
      IPQ2=IPQA(1)
      MU=1
C
C        CALL DXPQNU TO OBTAIN Q(1.,NU1,X)
C
      CALL DXPQNU(NU1,NU2,MU,THETA,ID,PQA,IPQA,IERROR)
      IF (IERROR.NE.0) RETURN
      NU=NU1
      K=0
      MU=1
      DMU=1.D0
      PQ1=PQA(1)
      IPQ1=IPQA(1)
      IF(MU1.GT.0) GO TO 310
      K=K+1
      PQA(K)=PQ2
      IPQA(K)=IPQ2
      IF(MU2.LT.1) GO TO 330
  310 IF(MU1.GT.1) GO TO 320
      K=K+1
      PQA(K)=PQ1
      IPQA(K)=IPQ1
      IF(MU2.LE.1) GO TO 330
  320 CONTINUE
C
C        FORWARD RECURRENCE IN MU TO OBTAIN
C                  Q(MU1,NU,X),Q(MU1+1,NU,X),....,Q(MU2,NU,X) USING
C             Q(MU+1,NU,X)=-2.*MU*X*SQRT(1./(1.-X**2))*Q(MU,NU,X)
C                               -(NU+MU)*(NU-MU+1.)*Q(MU-1,NU,X)
C
      X1=-2.D0*DMU*X*SX*PQ1
      X2=(NU+DMU)*(NU-DMU+1.D0)*PQ2
      CALL DXADD(X1,IPQ1,-X2,IPQ2,PQ,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
      CALL DXADJ(PQ,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
      PQ2=PQ1
      IPQ2=IPQ1
      PQ1=PQ
      IPQ1=IPQ
      MU=MU+1
      DMU=DMU+1.D0
      IF(MU.LT.MU1) GO TO 320
      K=K+1
      PQA(K)=PQ
      IPQA(K)=IPQ
      IF(MU2.GT.MU) GO TO 320
  330 RETURN
      END

C----------------------------------------------------------------------
C
      SUBROUTINE DXQNU (NU1, NU2, MU1, THETA, X, SX, ID, PQA, IPQA,
     1   IERROR)
      DIMENSION PQA(*),IPQA(*)
      DOUBLE PRECISION DMU,NU,NU1,NU2,PQ,PQA,PQ1,PQ2,SX,X,X1,X2
      DOUBLE PRECISION THETA,PQL1,PQL2
C***FIRST EXECUTABLE STATEMENT  DXQNU
      IERROR=0
      K=0
      PQ2=0.0D0
      IPQ2=0
      PQL2=0.0D0
      IPQL2=0
      IF(MU1.EQ.1) GO TO 290
      MU=0
C
C        CALL DXPQNU TO OBTAIN Q(0.,NU2,X) AND Q(0.,NU2-1,X)
C
      CALL DXPQNU(NU1,NU2,MU,THETA,ID,PQA,IPQA,IERROR)
      IF (IERROR.NE.0) RETURN
      IF(MU1.EQ.0) RETURN
      K=(NU2-NU1+1.5D0)
      PQ2=PQA(K)
      IPQ2=IPQA(K)
      PQL2=PQA(K-1)
      IPQL2=IPQA(K-1)
  290 MU=1
C
C        CALL DXPQNU TO OBTAIN Q(1.,NU2,X) AND Q(1.,NU2-1,X)
C
      CALL DXPQNU(NU1,NU2,MU,THETA,ID,PQA,IPQA,IERROR)
      IF (IERROR.NE.0) RETURN
      IF(MU1.EQ.1) RETURN
      NU=NU2
      PQ1=PQA(K)
      IPQ1=IPQA(K)
      PQL1=PQA(K-1)
      IPQL1=IPQA(K-1)
  300 MU=1
      DMU=1.D0
  320 CONTINUE
C
C        FORWARD RECURRENCE IN MU TO OBTAIN Q(MU1,NU2,X) AND
C              Q(MU1,NU2-1,X) USING
C              Q(MU+1,NU,X)=-2.*MU*X*SQRT(1./(1.-X**2))*Q(MU,NU,X)
C                   -(NU+MU)*(NU-MU+1.)*Q(MU-1,NU,X)
C
C              FIRST FOR NU=NU2
C
      X1=-2.D0*DMU*X*SX*PQ1
      X2=(NU+DMU)*(NU-DMU+1.D0)*PQ2
      CALL DXADD(X1,IPQ1,-X2,IPQ2,PQ,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
      CALL DXADJ(PQ,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
      PQ2=PQ1
      IPQ2=IPQ1
      PQ1=PQ
      IPQ1=IPQ
      MU=MU+1
      DMU=DMU+1.D0
      IF(MU.LT.MU1) GO TO 320
      PQA(K)=PQ
      IPQA(K)=IPQ
      IF(K.EQ.1) RETURN
      IF(NU.LT.NU2) GO TO 340
C
C              THEN FOR NU=NU2-1
C
      NU=NU-1.D0
      PQ2=PQL2
      IPQ2=IPQL2
      PQ1=PQL1
      IPQ1=IPQL1
      K=K-1
      GO TO 300
C
C         BACKWARD RECURRENCE IN NU TO OBTAIN
C              Q(MU1,NU1,X),Q(MU1,NU1+1,X),....,Q(MU1,NU2,X)
C              USING
C              (NU-MU+1.)*Q(MU,NU+1,X)=
C                       (2.*NU+1.)*X*Q(MU,NU,X)-(NU+MU)*Q(MU,NU-1,X)
C
  340 PQ1=PQA(K)
      IPQ1=IPQA(K)
      PQ2=PQA(K+1)
      IPQ2=IPQA(K+1)
  350 IF(NU.LE.NU1) RETURN
      K=K-1
      X1=(2.D0*NU+1.D0)*X*PQ1/(NU+DMU)
      X2=-(NU-DMU+1.D0)*PQ2/(NU+DMU)
      CALL DXADD(X1,IPQ1,X2,IPQ2,PQ,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
      CALL DXADJ(PQ,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
      PQ2=PQ1
      IPQ2=IPQ1
      PQ1=PQ
      IPQ1=IPQ
      PQA(K)=PQ
      IPQA(K)=IPQ
      NU=NU-1.D0
      GO TO 350
      END
C----------------------------------------------------------------------
C
      SUBROUTINE DXADD (X, IX, Y, IY, Z, IZ, IERROR)
      DOUBLE PRECISION X, Y, Z
      INTEGER IX, IY, IZ
      DOUBLE PRECISION RADIX, RADIXL, RAD2L, DLG10R
      INTEGER L, L2, KMAX
      COMMON /DXBLK2/ RADIX, RADIXL, RAD2L, DLG10R, L, L2, KMAX
      SAVE /DXBLK2/
      DOUBLE PRECISION S, T
C
C   THE CONDITIONS IMPOSED ON L AND KMAX BY THIS SUBROUTINE
C ARE
C     (1) 1 .LT. L .LE. 0.5D0*LOGR(0.5D0*DZERO)
C
C     (2) NRADPL .LT. L .LE. KMAX/6
C
C     (3) KMAX .LE. (2**NBITS - 4*L - 1)/2
C
C THESE CONDITIONS MUST BE MET BY APPROPRIATE CODING
C IN SUBROUTINE DXSET.
C
C***FIRST EXECUTABLE STATEMENT  DXADD
      IERROR=0
      IF (X.NE.0.0D0) GO TO 10
      Z = Y
      IZ = IY
      GO TO 220
   10 IF (Y.NE.0.0D0) GO TO 20
      Z = X
      IZ = IX
      GO TO 220
   20 CONTINUE
      IF (IX.GE.0 .AND. IY.GE.0) GO TO 40
      IF (IX.LT.0 .AND. IY.LT.0) GO TO 40
      IF (ABS(IX).LE.6*L .AND. ABS(IY).LE.6*L) GO TO 40
      IF (IX.GE.0) GO TO 30
      Z = Y
      IZ = IY
      GO TO 220
   30 CONTINUE
      Z = X
      IZ = IX
      GO TO 220
   40 I = IX - IY
      IF (I) 80, 50, 90
   50 IF (ABS(X).GT.1.0D0 .AND. ABS(Y).GT.1.0D0) GO TO 60
      IF (ABS(X).LT.1.0D0 .AND. ABS(Y).LT.1.0D0) GO TO 70
      Z = X + Y
      IZ = IX
      GO TO 220
   60 S = X/RADIXL
      T = Y/RADIXL
      Z = S + T
      IZ = IX + L
      GO TO 220
   70 S = X*RADIXL
      T = Y*RADIXL
      Z = S + T
      IZ = IX - L
      GO TO 220
   80 S = Y
      IS = IY
      T = X
      GO TO 100
   90 S = X
      IS = IX
      T = Y
  100 CONTINUE
C
C  AT THIS POINT, THE ONE OF (X,IX) OR (Y,IY) THAT HAS THE
C LARGER AUXILIARY INDEX IS STORED IN (S,IS). THE PRINCIPAL
C PART OF THE OTHER INPUT IS STORED IN T.
C
      I1 = ABS(I)/L
      I2 = MOD(ABS(I),L)
      IF (ABS(T).GE.RADIXL) GO TO 130
      IF (ABS(T).GE.1.0D0) GO TO 120
      IF (RADIXL*ABS(T).GE.1.0D0) GO TO 110
      J = I1 + 1
      T = T*RADIX**(L-I2)
      GO TO 140
  110 J = I1
      T = T*RADIX**(-I2)
      GO TO 140
  120 J = I1 - 1
      IF (J.LT.0) GO TO 110
      T = T*RADIX**(-I2)/RADIXL
      GO TO 140
  130 J = I1 - 2
      IF (J.LT.0) GO TO 120
      T = T*RADIX**(-I2)/RAD2L
  140 CONTINUE
C
C  AT THIS POINT, SOME OR ALL OF THE DIFFERENCE IN THE
C AUXILIARY INDICES HAS BEEN USED TO EFFECT A LEFT SHIFT
C OF T.  THE SHIFTED VALUE OF T SATISFIES
C
C       RADIX**(-2*L) .LE. ABS(T) .LE. 1.0D0
C
C AND, IF J=0, NO FURTHER SHIFTING REMAINS TO BE DONE.
C
      IF (J.EQ.0) GO TO 190
      IF (ABS(S).GE.RADIXL .OR. J.GT.3) GO TO 150
      IF (ABS(S).GE.1.0D0) GO TO (180, 150, 150), J
      IF (RADIXL*ABS(S).GE.1.0D0) GO TO (180, 170, 150), J
      GO TO (180, 170, 160), J
  150 Z = S
      IZ = IS
      GO TO 220
  160 S = S*RADIXL
  170 S = S*RADIXL
  180 S = S*RADIXL
  190 CONTINUE
C
C   AT THIS POINT, THE REMAINING DIFFERENCE IN THE
C AUXILIARY INDICES HAS BEEN USED TO EFFECT A RIGHT SHIFT
C OF S.  IF THE SHIFTED VALUE OF S WOULD HAVE EXCEEDED
C RADIX**L, THEN (S,IS) IS RETURNED AS THE VALUE OF THE
C SUM.
C
      IF (ABS(S).GT.1.0D0 .AND. ABS(T).GT.1.0D0) GO TO 200
      IF (ABS(S).LT.1.0D0 .AND. ABS(T).LT.1.0D0) GO TO 210
      Z = S + T
      IZ = IS - J*L
      GO TO 220
  200 S = S/RADIXL
      T = T/RADIXL
      Z = S + T
      IZ = IS - J*L + L
      GO TO 220
  210 S = S*RADIXL
      T = T*RADIXL
      Z = S + T
      IZ = IS - J*L - L
  220 CALL DXADJ(Z, IZ,IERROR)
      IF (IERROR.NE.0) RETURN
      RETURN
      END
C---------------------------------------------------------------
C
      SUBROUTINE DXADJ (X, IX, IERROR)
      DOUBLE PRECISION X
      INTEGER IX
      DOUBLE PRECISION RADIX, RADIXL, RAD2L, DLG10R
      INTEGER L, L2, KMAX
      COMMON /DXBLK2/ RADIX, RADIXL, RAD2L, DLG10R, L, L2, KMAX
      SAVE /DXBLK2/
C
C   THE CONDITION IMPOSED ON L AND KMAX BY THIS SUBROUTINE
C IS
C     2*L .LE. KMAX
C
C THIS CONDITION MUST BE MET BY APPROPRIATE CODING
C IN SUBROUTINE DXSET.
C
C***FIRST EXECUTABLE STATEMENT  DXADJ
      IERROR=0
      IF (X.EQ.0.0D0) GO TO 50
      IF (ABS(X).GE.1.0D0) GO TO 20
      IF (RADIXL*ABS(X).GE.1.0D0) GO TO 60
      X = X*RAD2L
      IF (IX.LT.0) GO TO 10
      IX = IX - L2
      GO TO 70
   10 IF (IX.LT.-KMAX+L2) GO TO 40
      IX = IX - L2
      GO TO 70
   20 IF (ABS(X).LT.RADIXL) GO TO 60
      X = X/RAD2L
      IF (IX.GT.0) GO TO 30
      IX = IX + L2
      GO TO 70
   30 IF (IX.GT.KMAX-L2) GO TO 40
      IX = IX + L2
      GO TO 70
   40 CALL XERMSG ('SLATEC', 'DXADJ', 'overflow in auxiliary index',
     +             207, 1)
      IERROR=207
      RETURN
   50 IX = 0
   60 IF (ABS(IX).GT.KMAX) GO TO 40
   70 RETURN
      END
C---------------------------------------------------------------
C
      SUBROUTINE DXRED (X, IX, IERROR)
      DOUBLE PRECISION X
      INTEGER IX
      DOUBLE PRECISION RADIX, RADIXL, RAD2L, DLG10R, XA
      INTEGER L, L2, KMAX
      COMMON /DXBLK2/ RADIX, RADIXL, RAD2L, DLG10R, L, L2, KMAX
      SAVE /DXBLK2/
C
C***FIRST EXECUTABLE STATEMENT  DXRED
      IERROR=0
      IF (X.EQ.0.0D0) GO TO 90
      XA = ABS(X)
      IF (IX.EQ.0) GO TO 70
      IXA = ABS(IX)
      IXA1 = IXA/L2
      IXA2 = MOD(IXA,L2)
      IF (IX.GT.0) GO TO 40
   10 CONTINUE
      IF (XA.GT.1.0D0) GO TO 20
      XA = XA*RAD2L
      IXA1 = IXA1 + 1
      GO TO 10
   20 XA = XA/RADIX**IXA2
      IF (IXA1.EQ.0) GO TO 70
      DO 30 I=1,IXA1
        IF (XA.LT.1.0D0) GO TO 100
        XA = XA/RAD2L
   30 CONTINUE
      GO TO 70
C
   40 CONTINUE
      IF (XA.LT.1.0D0) GO TO 50
      XA = XA/RAD2L
      IXA1 = IXA1 + 1
      GO TO 40
   50 XA = XA*RADIX**IXA2
      IF (IXA1.EQ.0) GO TO 70
      DO 60 I=1,IXA1
        IF (XA.GT.1.0D0) GO TO 100
        XA = XA*RAD2L
   60 CONTINUE
   70 IF (XA.GT.RAD2L) GO TO 100
      IF (XA.GT.1.0D0) GO TO 80
      IF (RAD2L*XA.LT.1.0D0) GO TO 100
   80 X = SIGN(XA,X)
   90 IX = 0
  100 RETURN
      END
C--------------------------------------------------------------------
C
      SUBROUTINE DXSET (IRAD, NRADPL, DZERO, NBITS, IERROR)
      INTEGER IRAD, NRADPL, NBITS
      DOUBLE PRECISION DZERO, DZEROX
      COMMON /DXBLK1/ NBITSF
      SAVE /DXBLK1/
      DOUBLE PRECISION RADIX, RADIXL, RAD2L, DLG10R
      INTEGER L, L2, KMAX
      COMMON /DXBLK2/ RADIX, RADIXL, RAD2L, DLG10R, L, L2, KMAX
      SAVE /DXBLK2/
      INTEGER NLG102, MLG102, LG102
      COMMON /DXBLK3/ NLG102, MLG102, LG102(21)
      SAVE /DXBLK3/
      INTEGER IFLAG
      SAVE IFLAG
C
      DIMENSION LOG102(20), LGTEMP(20)
      SAVE LOG102
C
C   LOG102 CONTAINS THE FIRST 60 DIGITS OF LOG10(2) FOR USE IN
C CONVERSION OF EXTENDED-RANGE NUMBERS TO BASE 10 .
      DATA LOG102 /301,029,995,663,981,195,213,738,894,724,493,026,768,
     * 189,881,462,108,541,310,428/
C
C FOLLOWING CODING PREVENTS DXSET FROM BEING EXECUTED MORE THAN ONCE.
C THIS IS IMPORTANT BECAUSE SOME SUBROUTINES (SUCH AS DXNRMP AND
C DXLEGF) CALL DXSET TO MAKE SURE EXTENDED-RANGE ARITHMETIC HAS
C BEEN INITIALIZED. THE USER MAY WANT TO PRE-EMPT THIS CALL, FOR
C EXAMPLE WHEN I1MACH IS NOT AVAILABLE. SEE CODING BELOW.
      DATA IFLAG /0/
C***FIRST EXECUTABLE STATEMENT  DXSET
      IERROR=0
      IF (IFLAG .NE. 0) RETURN
      IRADX = IRAD
      NRDPLC = NRADPL
      DZEROX = DZERO
      IMINEX = 0
      IMAXEX = 0
      NBITSX = NBITS
C FOLLOWING 5 STATEMENTS SHOULD BE DELETED IF I1MACH IS
C NOT AVAILABLE OR NOT CONFIGURED TO RETURN THE CORRECT
C MACHINE-DEPENDENT VALUES.
      IF (IRADX .EQ. 0) IRADX = I1MACH (10)
      IF (NRDPLC .EQ. 0) NRDPLC = I1MACH (14)
      IF (DZEROX .EQ. 0.0D0) IMINEX = I1MACH (15)
      IF (DZEROX .EQ. 0.0D0) IMAXEX = I1MACH (16)
      IF (NBITSX .EQ. 0) NBITSX = I1MACH (8)
      IF (IRADX.EQ.2) GO TO 10
      IF (IRADX.EQ.4) GO TO 10
      IF (IRADX.EQ.8) GO TO 10
      IF (IRADX.EQ.16) GO TO 10
      CALL XERMSG ('SLATEC', 'DXSET', 'IMPROPER VALUE OF IRAD', 201, 1)
      IERROR=201
      RETURN
   10 CONTINUE
      LOG2R=0
      IF (IRADX.EQ.2) LOG2R = 1
      IF (IRADX.EQ.4) LOG2R = 2
      IF (IRADX.EQ.8) LOG2R = 3
      IF (IRADX.EQ.16) LOG2R = 4
      NBITSF=LOG2R*NRDPLC
      RADIX = IRADX
      DLG10R = LOG10(RADIX)
      IF (DZEROX .NE. 0.0D0) GO TO 14
      LX = MIN ((1-IMINEX)/2, (IMAXEX-1)/2)
      GO TO 16
   14 LX = 0.5D0*LOG10(DZEROX)/DLG10R
C RADIX**(2*L) SHOULD NOT OVERFLOW, BUT REDUCE L BY 1 FOR FURTHER
C PROTECTION.
      LX=LX-1
   16 L2 = 2*LX
      IF (LX.GE.4) GO TO 20
      CALL XERMSG ('SLATEC', 'DXSET', 'IMPROPER VALUE OF DZERO', 202, 1)
      IERROR=202
      RETURN
   20 L = LX
      RADIXL = RADIX**L
      RAD2L = RADIXL**2
C    IT IS NECESSARY TO RESTRICT NBITS (OR NBITSX) TO BE LESS THAN SOME
C UPPER LIMIT BECAUSE OF BINARY-TO-DECIMAL CONVERSION. SUCH CONVERSION
C IS DONE BY DXC210 AND REQUIRES A CONSTANT THAT IS STORED TO SOME FIXED
C PRECISION. THE STORED CONSTANT (LOG102 IN THIS ROUTINE) PROVIDES
C FOR CONVERSIONS ACCURATE TO THE LAST DECIMAL DIGIT WHEN THE INTEGER
C WORD LENGTH DOES NOT EXCEED 63. A LOWER LIMIT OF 15 BITS IS IMPOSED
C BECAUSE THE SOFTWARE IS DESIGNED TO RUN ON COMPUTERS WITH INTEGER WORD
C LENGTH OF AT LEAST 16 BITS.
      IF (15.LE.NBITSX .AND. NBITSX.LE.63) GO TO 30
      CALL XERMSG ('SLATEC', 'DXSET', 'IMPROPER VALUE OF NBITS', 203, 1)
      IERROR=203
      RETURN
   30 CONTINUE
      KMAX = 2**(NBITSX-1) - L2
      NB = (NBITSX-1)/2
      MLG102 = 2**NB
      IF (1.LE.NRDPLC*LOG2R .AND. NRDPLC*LOG2R.LE.120) GO TO 40
      CALL XERMSG ('SLATEC', 'DXSET', 'IMPROPER VALUE OF NRADPL', 204,
     +             1)
      IERROR=204
      RETURN
   40 CONTINUE
      NLG102 = NRDPLC*LOG2R/NB + 3
      NP1 = NLG102 + 1
C
C   AFTER COMPLETION OF THE FOLLOWING LOOP, IC CONTAINS
C THE INTEGER PART AND LGTEMP CONTAINS THE FRACTIONAL PART
C OF LOG10(IRADX) IN RADIX 1000.
      IC = 0
      DO 50 II=1,20
        I = 21 - II
        IT = LOG2R*LOG102(I) + IC
        IC = IT/1000
        LGTEMP(I) = MOD(IT,1000)
   50 CONTINUE
C
C   AFTER COMPLETION OF THE FOLLOWING LOOP, LG102 CONTAINS
C LOG10(IRADX) IN RADIX MLG102. THE RADIX POINT IS
C BETWEEN LG102(1) AND LG102(2).
      LG102(1) = IC
      DO 80 I=2,NP1
        LG102X = 0
        DO 70 J=1,NB
          IC = 0
          DO 60 KK=1,20
            K = 21 - KK
            IT = 2*LGTEMP(K) + IC
            IC = IT/1000
            LGTEMP(K) = MOD(IT,1000)
   60     CONTINUE
          LG102X = 2*LG102X + IC
   70   CONTINUE
        LG102(I) = LG102X
   80 CONTINUE
C
C CHECK SPECIAL CONDITIONS REQUIRED BY SUBROUTINES...
      IF (NRDPLC.LT.L) GO TO 90
      CALL XERMSG ('SLATEC', 'DXSET', 'NRADPL .GE. L', 205, 1)
      IERROR=205
      RETURN
   90 IF (6*L.LE.KMAX) GO TO 100
      CALL XERMSG ('SLATEC', 'DXSET', '6*L .GT. KMAX', 206, 1)
      IERROR=206
      RETURN
  100 CONTINUE
      IFLAG = 1
      RETURN
      END
C------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION DXPSI (A, IPSIK, IPSIX)
      DOUBLE PRECISION A,B,C,CNUM,CDENOM
      DIMENSION CNUM(12),CDENOM(12)
      SAVE CNUM, CDENOM
C
C        CNUM(I) AND CDENOM(I) ARE THE ( REDUCED ) NUMERATOR
C        AND 2*I*DENOMINATOR RESPECTIVELY OF THE 2*I TH BERNOULLI
C        NUMBER.
C
      DATA CNUM(1),CNUM(2),CNUM(3),CNUM(4),CNUM(5),CNUM(6),CNUM(7),
     1CNUM(8),CNUM(9),CNUM(10),CNUM(11),CNUM(12)
     2    / 1.D0,     -1.D0,    1.D0,     -1.D0, 1.D0,
     3   -691.D0,  1.D0,     -3617.D0, 43867.D0, -174611.D0, 77683.D0,
     4   -236364091.D0/
      DATA CDENOM(1),CDENOM(2),CDENOM(3),CDENOM(4),CDENOM(5),CDENOM(6),
     1 CDENOM(7),CDENOM(8),CDENOM(9),CDENOM(10),CDENOM(11),CDENOM(12)
     2/12.D0,120.D0,   252.D0,   240.D0,132.D0,
     3  32760.D0, 12.D0,  8160.D0, 14364.D0, 6600.D0, 276.D0, 65520.D0/
C***FIRST EXECUTABLE STATEMENT  DXPSI
      N=MAX(0,IPSIX-INT(A))
      B=N+A
      K1=IPSIK-1
C
C        SERIES EXPANSION FOR A .GT. IPSIX USING IPSIK-1 TERMS.
C
      C=0.D0
      DO 12 I=1,K1
      K=IPSIK-I
   12 C=(C+CNUM(K)/CDENOM(K))/B**2
      DXPSI=LOG(B)-(C+.5D0/B)
      IF(N.EQ.0) GO TO 20
      B=0.D0
C
C        RECURRENCE FOR A .LE. IPSIX.
C
      DO 15 M=1,N
   15 B=B+1.D0/(N-M+A)
      DXPSI=DXPSI-B
   20 RETURN
      END
C--------------------------------------------------------------------
C
      SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8 XLIBR, XSUBR
      CHARACTER*72  TEMP
      CHARACTER*20  LFIRST
C***FIRST EXECUTABLE STATEMENT  XERMSG
      LKNTRL = J4SAVE (2, 0, .FALSE.)
      MAXMES = J4SAVE (4, 0, .FALSE.)
C
C       LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL.
C       MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE
C          SHOULD BE PRINTED.
C
C       WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN
C          CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE,
C          AND THE LEVEL SHOULD BE BETWEEN 0 AND 2.
C
      IF (NERR.LT.-9999999 .OR. NERR.GT.99999999 .OR. NERR.EQ.0 .OR.
     *   LEVEL.LT.-1 .OR. LEVEL.GT.2) THEN
         CALL XERPRN (' ***', -1, 'FATAL ERROR IN...$$ ' //
     *      'XERMSG -- INVALID ERROR NUMBER OR LEVEL$$ '//
     *      'JOB ABORT DUE TO FATAL ERROR.', 72)
         CALL XERSVE (' ', ' ', ' ', 0, 0, 0, KDUMMY)
         CALL XERHLT (' ***XERMSG -- INVALID INPUT')
         RETURN
      ENDIF
C
C       RECORD THE MESSAGE.
C
      I = J4SAVE (1, NERR, .TRUE.)
      CALL XERSVE (LIBRAR, SUBROU, MESSG, 1, NERR, LEVEL, KOUNT)
C
C       HANDLE PRINT-ONCE WARNING MESSAGES.
C
      IF (LEVEL.EQ.-1 .AND. KOUNT.GT.1) RETURN
C
C       ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG.
C
      XLIBR  = LIBRAR
      XSUBR  = SUBROU
      LFIRST = MESSG
      LERR   = NERR
      LLEVEL = LEVEL
      CALL XERCNT (XLIBR, XSUBR, LFIRST, LERR, LLEVEL, LKNTRL)
C
      LKNTRL = MAX(-2, MIN(2,LKNTRL))
      MKNTRL = ABS(LKNTRL)
C
C       SKIP PRINTING IF THE CONTROL FLAG VALUE AS RESET IN XERCNT IS
C       ZERO AND THE ERROR IS NOT FATAL.
C
      IF (LEVEL.LT.2 .AND. LKNTRL.EQ.0) GO TO 30
      IF (LEVEL.EQ.0 .AND. KOUNT.GT.MAXMES) GO TO 30
      IF (LEVEL.EQ.1 .AND. KOUNT.GT.MAXMES .AND. MKNTRL.EQ.1) GO TO 30
      IF (LEVEL.EQ.2 .AND. KOUNT.GT.MAX(1,MAXMES)) GO TO 30
C
C       ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A
C       MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS)
C       AND SENDING IT OUT VIA XERPRN.  PRINT ONLY IF CONTROL FLAG
C       IS NOT ZERO.
C
      IF (LKNTRL .NE. 0) THEN
         TEMP(1:21) = 'MESSAGE FROM ROUTINE '
         I = MIN(LEN(SUBROU), 16)
         TEMP(22:21+I) = SUBROU(1:I)
         TEMP(22+I:33+I) = ' IN LIBRARY '
         LTEMP = 33 + I
         I = MIN(LEN(LIBRAR), 16)
         TEMP(LTEMP+1:LTEMP+I) = LIBRAR (1:I)
         TEMP(LTEMP+I+1:LTEMP+I+1) = '.'
         LTEMP = LTEMP + I + 1
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
C
C       IF LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE
C       PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE
C       FROM EACH OF THE FOLLOWING THREE OPTIONS.
C       1.  LEVEL OF THE MESSAGE
C              'INFORMATIVE MESSAGE'
C              'POTENTIALLY RECOVERABLE ERROR'
C              'FATAL ERROR'
C       2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE
C              'PROG CONTINUES'
C              'PROG ABORTED'
C       3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK
C           MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS
C           WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.)
C              'TRACEBACK REQUESTED'
C              'TRACEBACK NOT REQUESTED'
C       NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT
C       EXCEED 74 CHARACTERS.
C       WE SKIP THE NEXT BLOCK IF THE INTRODUCTORY LINE IS NOT NEEDED.
C
      IF (LKNTRL .GT. 0) THEN
C
C       THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL.
C
         IF (LEVEL .LE. 0) THEN
            TEMP(1:20) = 'INFORMATIVE MESSAGE,'
            LTEMP = 20
         ELSEIF (LEVEL .EQ. 1) THEN
            TEMP(1:30) = 'POTENTIALLY RECOVERABLE ERROR,'
            LTEMP = 30
         ELSE
            TEMP(1:12) = 'FATAL ERROR,'
            LTEMP = 12
         ENDIF
C
C       THEN WHETHER THE PROGRAM WILL CONTINUE.
C
         IF ((MKNTRL.EQ.2 .AND. LEVEL.GE.1) .OR.
     *       (MKNTRL.EQ.1 .AND. LEVEL.EQ.2)) THEN
            TEMP(LTEMP+1:LTEMP+14) = ' PROG ABORTED,'
            LTEMP = LTEMP + 14
         ELSE
            TEMP(LTEMP+1:LTEMP+16) = ' PROG CONTINUES,'
            LTEMP = LTEMP + 16
         ENDIF
C
C       FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK.
C
         IF (LKNTRL .GT. 0) THEN
            TEMP(LTEMP+1:LTEMP+20) = ' TRACEBACK REQUESTED'
            LTEMP = LTEMP + 20
         ELSE
            TEMP(LTEMP+1:LTEMP+24) = ' TRACEBACK NOT REQUESTED'
            LTEMP = LTEMP + 24
         ENDIF
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
C
C       NOW SEND OUT THE MESSAGE.
C
      CALL XERPRN (' *  ', -1, MESSG, 72)
C
C       IF LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A
C          TRACEBACK.
C
      IF (LKNTRL .GT. 0) THEN
         WRITE (TEMP, '(''ERROR NUMBER = '', I8)') NERR
         DO 10 I=16,22
            IF (TEMP(I:I) .NE. ' ') GO TO 20
   10    CONTINUE
C
   20    CALL XERPRN (' *  ', -1, TEMP(1:15) // TEMP(I:23), 72)
         CALL FDUMP
      ENDIF
C
C       IF LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE.
C
      IF (LKNTRL .NE. 0) THEN
         CALL XERPRN (' *  ', -1, ' ', 72)
         CALL XERPRN (' ***', -1, 'END OF MESSAGE', 72)
         CALL XERPRN ('    ',  0, ' ', 72)
      ENDIF
C
C       IF THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE
C       CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN.
C
   30 IF (LEVEL.LE.0 .OR. (LEVEL.EQ.1 .AND. MKNTRL.LE.1)) RETURN
C
C       THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A
C       FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR
C       SUMMARY IF THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT.
C
      IF (LKNTRL.GT.0 .AND. KOUNT.LT.MAX(1,MAXMES)) THEN
         IF (LEVEL .EQ. 1) THEN
            CALL XERPRN
     *         (' ***', -1, 'JOB ABORT DUE TO UNRECOVERED ERROR.', 72)
         ELSE
            CALL XERPRN(' ***', -1, 'JOB ABORT DUE TO FATAL ERROR.', 72)
         ENDIF
         CALL XERSVE (' ', ' ', ' ', -1, 0, 0, KDUMMY)
         CALL XERHLT (' ')
      ELSE
         CALL XERHLT (MESSG)
      ENDIF
      RETURN
      END
C-------------------------------------------------------------------
C
      SUBROUTINE XERCNT (LIBRAR, SUBROU, MESSG, NERR, LEVEL, KONTRL)
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
C***FIRST EXECUTABLE STATEMENT  XERCNT
      RETURN
      END
C------------------------------------------------------------------
C
      SUBROUTINE XERHLT (MESSG)
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERHLT
      STOP
      END
C------------------------------------------------------------------
C
      SUBROUTINE XERPRN (PREFIX, NPREF, MESSG, NWRAP)
      CHARACTER*(*) PREFIX, MESSG
      INTEGER NPREF, NWRAP
      CHARACTER*148 CBUFF
      INTEGER IU(5), NUNIT
      CHARACTER*2 NEWLIN
      PARAMETER (NEWLIN = '$$')
C***FIRST EXECUTABLE STATEMENT  XERPRN
      CALL XGETUA(IU,NUNIT)
C
C       A ZERO VALUE FOR A LOGICAL UNIT NUMBER MEANS TO USE THE STANDARD
C       ERROR MESSAGE UNIT INSTEAD.  I1MACH(4) RETRIEVES THE STANDARD
C       ERROR MESSAGE UNIT.
C
      N = I1MACH(4)
      DO 10 I=1,NUNIT
         IF (IU(I) .EQ. 0) IU(I) = N
   10 CONTINUE
C
C       LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE
C       BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING
C       THE REST OF THIS ROUTINE.
C
      IF ( NPREF .LT. 0 ) THEN
         LPREF = LEN(PREFIX)
      ELSE
         LPREF = NPREF
      ENDIF
      LPREF = MIN(16, LPREF)
      IF (LPREF .NE. 0) CBUFF(1:LPREF) = PREFIX
C
C       LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE
C       TIME FROM MESSG TO PRINT ON ONE LINE.
C
      LWRAP = MAX(16, MIN(132, NWRAP))
C
C       SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS.
C
      LENMSG = LEN(MESSG)
      N = LENMSG
      DO 20 I=1,N
         IF (MESSG(LENMSG:LENMSG) .NE. ' ') GO TO 30
         LENMSG = LENMSG - 1
   20 CONTINUE
   30 CONTINUE
C
C       IF THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE.
C
      IF (LENMSG .EQ. 0) THEN
         CBUFF(LPREF+1:LPREF+1) = ' '
         DO 40 I=1,NUNIT
            WRITE(IU(I), '(A)') CBUFF(1:LPREF+1)
   40    CONTINUE
         RETURN
      ENDIF
C
C       SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING
C       STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL.
C       WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT.
C       WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED.
C
C       WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE
C       INDEX INTRINSIC FUNCTION RETURNS ZERO IF THERE IS NO OCCURRENCE
C       OR IF THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH
C       OF THE SECOND ARGUMENT.
C
C       THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE
C       FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER
C       OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT
C       POSITION NEXTC.
C
C       LPIECE .EQ. 0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE
C                       REMAINDER OF THE CHARACTER STRING.  LPIECE
C                       SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC,
C                       WHICHEVER IS LESS.
C
C       LPIECE .EQ. 1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC:
C                       NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE
C                       PRINT NOTHING TO AVOID PRODUCING UNNECESSARY
C                       BLANK LINES.  THIS TAKES CARE OF THE SITUATION
C                       WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF
C                       EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE
C                       SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC
C                       SHOULD BE INCREMENTED BY 2.
C
C       LPIECE .GT. LWRAP+1  REDUCE LPIECE TO LWRAP.
C
C       ELSE            THIS LAST CASE MEANS 2 .LE. LPIECE .LE. LWRAP+1
C                       RESET LPIECE = LPIECE-1.  NOTE THAT THIS
C                       PROPERLY HANDLES THE END CASE WHERE LPIECE .EQ.
C                       LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY
C                       AT THE END OF A LINE.
C
      NEXTC = 1
   50 LPIECE = INDEX(MESSG(NEXTC:LENMSG), NEWLIN)
      IF (LPIECE .EQ. 0) THEN
C
C       THERE WAS NO NEW LINE SENTINEL FOUND.
C
         IDELTA = 0
         LPIECE = MIN(LWRAP, LENMSG+1-NEXTC)
         IF (LPIECE .LT. LENMSG+1-NEXTC) THEN
            DO 52 I=LPIECE+1,2,-1
               IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
                  LPIECE = I-1
                  IDELTA = 1
                  GOTO 54
               ENDIF
   52       CONTINUE
         ENDIF
   54    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSEIF (LPIECE .EQ. 1) THEN
C
C       WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1).
C       DON'T PRINT A BLANK LINE.
C
         NEXTC = NEXTC + 2
         GO TO 50
      ELSEIF (LPIECE .GT. LWRAP+1) THEN
C
C       LPIECE SHOULD BE SET DOWN TO LWRAP.
C
         IDELTA = 0
         LPIECE = LWRAP
         DO 56 I=LPIECE+1,2,-1
            IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
               LPIECE = I-1
               IDELTA = 1
               GOTO 58
            ENDIF
   56    CONTINUE
   58    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSE
C
C       IF WE ARRIVE HERE, IT MEANS 2 .LE. LPIECE .LE. LWRAP+1.
C       WE SHOULD DECREMENT LPIECE BY ONE.
C
         LPIECE = LPIECE - 1
         CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC  = NEXTC + LPIECE + 2
      ENDIF
C
C       PRINT
C
      DO 60 I=1,NUNIT
         WRITE(IU(I), '(A)') CBUFF(1:LPREF+LPIECE)
   60 CONTINUE
C
      IF (NEXTC .LE. LENMSG) GO TO 50
      RETURN
      END
C---------------------------------------------------------------
C
      SUBROUTINE XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL,
     +   ICOUNT)
      PARAMETER (LENTAB=10)
      INTEGER LUN(5)
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8  LIBTAB(LENTAB), SUBTAB(LENTAB), LIB, SUB
      CHARACTER*20 MESTAB(LENTAB), MES
      DIMENSION NERTAB(LENTAB), LEVTAB(LENTAB), KOUNT(LENTAB)
      SAVE LIBTAB, SUBTAB, MESTAB, NERTAB, LEVTAB, KOUNT, KOUNTX, NMSG
      DATA KOUNTX/0/, NMSG/0/
C***FIRST EXECUTABLE STATEMENT  XERSVE
C
      IF (KFLAG.LE.0) THEN
C
C        Dump the table.
C
         IF (NMSG.EQ.0) RETURN
C
C        Print to each unit.
C
         CALL XGETUA (LUN, NUNIT)
         DO 20 KUNIT = 1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
C
C           Print the table header.
C
            WRITE (IUNIT,9000)
C
C           Print body of table.
C
            DO 10 I = 1,NMSG
               WRITE (IUNIT,9010) LIBTAB(I), SUBTAB(I), MESTAB(I),
     *            NERTAB(I),LEVTAB(I),KOUNT(I)
   10       CONTINUE
C
C           Print number of other errors.
C
            IF (KOUNTX.NE.0) WRITE (IUNIT,9020) KOUNTX
            WRITE (IUNIT,9030)
   20    CONTINUE
C
C        Clear the error tables.
C
         IF (KFLAG.EQ.0) THEN
            NMSG = 0
            KOUNTX = 0
         ENDIF
      ELSE
C
C        PROCESS A MESSAGE...
C        SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
C        OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
C
         LIB = LIBRAR
         SUB = SUBROU
         MES = MESSG
         DO 30 I = 1,NMSG
            IF (LIB.EQ.LIBTAB(I) .AND. SUB.EQ.SUBTAB(I) .AND.
     *         MES.EQ.MESTAB(I) .AND. NERR.EQ.NERTAB(I) .AND.
     *         LEVEL.EQ.LEVTAB(I)) THEN
                  KOUNT(I) = KOUNT(I) + 1
                  ICOUNT = KOUNT(I)
                  RETURN
            ENDIF
   30    CONTINUE
C
         IF (NMSG.LT.LENTAB) THEN
C
C           Empty slot found for new message.
C
            NMSG = NMSG + 1
            LIBTAB(I) = LIB
            SUBTAB(I) = SUB
            MESTAB(I) = MES
            NERTAB(I) = NERR
            LEVTAB(I) = LEVEL
            KOUNT (I) = 1
            ICOUNT    = 1
         ELSE
C
C           Table is full.
C
            KOUNTX = KOUNTX+1
            ICOUNT = 0
         ENDIF
      ENDIF
      RETURN
C
C     Formats.
C
 9000 FORMAT ('0          ERROR MESSAGE SUMMARY' /
     +   ' LIBRARY    SUBROUTINE MESSAGE START             NERR',
     +   '     LEVEL     COUNT')
 9010 FORMAT (1X,A,3X,A,3X,A,3I10)
 9020 FORMAT ('0OTHER ERRORS NOT INDIVIDUALLY TABULATED = ', I10)
 9030 FORMAT (1X)
      END
C--------------------------------------------------------------------
C
      SUBROUTINE XGETUA (IUNITA, N)
      DIMENSION IUNITA(5)
C***FIRST EXECUTABLE STATEMENT  XGETUA
      N = J4SAVE(5,0,.FALSE.)
      DO 30 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
      RETURN
      END
C--------------------------------------------------------------------
C
      FUNCTION J4SAVE (IWHICH, IVALUE, ISET)
      LOGICAL ISET
      INTEGER IPARAM(9)
      SAVE IPARAM
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
C***FIRST EXECUTABLE STATEMENT  J4SAVE
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
      END
C-----------------------------------------------------------------
C
      INTEGER FUNCTION I1MACH (I)
C
      INTEGER IMACH(16),OUTPUT
      SAVE IMACH
      EQUIVALENCE (IMACH(4),OUTPUT)
C
C     MACHINE CONSTANTS FOR THE IBM PC  (32 bit ???)
C
      DATA IMACH( 1) /          5 /
      DATA IMACH( 2) /          6 /
      DATA IMACH( 3) /          0 /
      DATA IMACH( 4) /          0 /
      DATA IMACH( 5) /         32 /
      DATA IMACH( 6) /          4 /
      DATA IMACH( 7) /          2 /
      DATA IMACH( 8) /         31 /
      DATA IMACH( 9) / 2147483647 /
      DATA IMACH(10) /          2 /
      DATA IMACH(11) /         24 /
      DATA IMACH(12) /       -125 /
      DATA IMACH(13) /        127 /
      DATA IMACH(14) /         53 /
      DATA IMACH(15) /      -1021 /
      DATA IMACH(16) /       1023 /
C
C***FIRST EXECUTABLE STATEMENT  I1MACH
      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 10
C
      I1MACH = IMACH(I)
      RETURN
C
   10 CONTINUE
      WRITE (UNIT = OUTPUT, FMT = 9000)
 9000 FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')
C
C     CALL FDUMP
C
      STOP
      END
C-------------------------------------------------------------------
C
      SUBROUTINE FDUMP
C 'Symbolic Dump (should be locally written)' ???
C***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
      END
C
C  So, final procedure
C  FDUMP 'should be locally written', and it is NOT written.
C  Obviously SLATEC error handling would not work.
C  But using in fact the only SLATEC procedure DXLEGF we probably
C  will not encounter run-time errors at all.
C  Thus all error-handling procedures are left in the code
C  just to please a linker (another way would be to replace
C  upper level error-handling procedures with their void analogs).
