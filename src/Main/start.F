      SUBROUTINE START
*
*
*       Initialization of data & polynomials.
*       ------------------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'timing.h'
      DIMENSION XRELPL(3),VRELPL(3),EI(3)
      EXTERNAL SCALE,MERGE
*      COMMON/NNBCOT/ nnbave
      PARAMETER  (NS=12)
*
*
*       Initialize global scalars, counters & useful constants.
      CALL ZERO
*
*       Read input parameters.
      CALL INPUT
*
*       Open all Files.
      CALL FILE_INIT(0)
*
*       Set initial conditions: BODY(I), X(K,I), XDOT(K,I); I=1,N & K=1,3.
      CALL DATA
*
*       Scale initial conditions to new units.
      CALL SCALE
*
*       Set total mass in case routines DATA & SCALE are not used.
      ZMASS = 0.0D0
      DO 10 I = 1,N
          ZMASS = ZMASS + BODY(I)
   10 CONTINUE
*     --16/05/07 01:16-qishu-debug----------------------------*
***** Note:--------------------------------------------------**
c warning when nmass >= N
          if (nmass .ge. N) then
              print *,'ERROR: NUMBER OF MASSLESS PARTICLES=',nmass,
     & 'IS LARGER THAN N=',N
              stop
          end if
*     --16/05/07 01:16-qishu-end------------------------------*
*
*       Define mean mass in scaled units and solar mass conversion factor.
      BODYM = ZMASS/FLOAT(N-nmass)
      ZMBAR = ZMBAR/BODYM
*
*       Introduce scaling factors DAYS, YRS, SU, RAU, SMU, TSTAR & VSTAR.
      CALL UNITS
*
*       Check option for external force.
      IF (KZ(14).GT.0) THEN
#ifdef TT
*** FlorentR - initialize the tensors
          IF (KZ(14).EQ.9) CALL TTINIT
*** FRenaud
#endif
          CALL XTRNL0
      END IF
*
*       Check optional scaling to hot system.
      IF (KZ(29).GT.0) THEN
          CALL HOTSYS
      END IF
*
*       Check option for initial binaries.
      IF (KZ(8).EQ.1.OR.KZ(8).GE.3) THEN
          CALL BINPOP
      END IF
*
*       Include stable primordial triples.
      IF (KZ(18).GT.1.AND.KZ(8).GT.0) THEN
          CALL HIPOP
      END IF

*       Initial massive black hole input
      IF (KZ(24).EQ.1) CALL IMBHINIT
*
*       Set sequential name, maximum mass & primary velocity, initial single/binary total mass
      BODY1 = 0.0
      ZSMASS0 = 0.0
      ZBMASS0 = 0.0
      DO 20 I = 1,N
*     exclude massive black hole mass in initial cluster mass
          IF(KZ(24).EQ.1.AND.I.EQ.N) GO TO 19
          IF(I.LE.2*NBIN0) THEN
             ZBMASS0 = ZBMASS0 + BODY(I)
          ELSE
             ZSMASS0 = ZSMASS0 + BODY(I)
          END IF
 19       NAME(I) = I
*       Input data should mark planets with mass <= smallMass
c          print *,i,body(i)
          if(body(i).gt.0.0d0)then
              if (body(i).le.smallMass * 1.1)then
                  nomass(i) = 1
              else
                  nomass(i) = 0
              end if
          else
              print *, "particle mass is not a positive value!!!"
              nomass(i) = 1
          end if
          BODY1 = MAX(BODY1,BODY(I))
          DO 15 K = 1,3
              X0DOT(K,I) = XDOT(K,I)
   15     CONTINUE
   20 CONTINUE
*
*       Initialize fixed block steps (64 levels).
      CALL IBLOCK
*
*       Create table of inverse Stumpff coefficients.
      DO 30 I = 1,NS
          SCOEFF(I) = 1.0/((I + 1)*(I + 2))
   30 CONTINUE
*
*       Set optional stellar evolution parameters.
      IF (KZ(19).GT.2) THEN
          CALL INSTAR
      ELSE IF (KZ(14).GT.1) THEN
          DT = 1.0E-03/TSCALE
          CALL STEPK(DT,DTN)
          STEPX = DTN
      END IF
*
*       Initialize optional cloud parameters.
      IF (KZ(13).GT.0) THEN
          CALL CLOUD0
      END IF
*
*       Regularize any hard primordial binaries (assume sequential ordering).
      IF (KZ(8).GT.0.AND.NBIN0.GT.0) THEN
          SMMIN = 1.D30
          SMMAX = 0.D0
          XMMIN = 1.D30
          XMMAX = 0.D0
          TMMIN = 1.D30
          TMMAX = 0.D0
*
*     temperately set kstart = -1 for KSREG
          KSDUM = KSTART
          KSTART = -1
          DO 50 IPAIR = 1,NBIN0
              ICOMP = 2*IPAIR - 1
              JCOMP = 2*IPAIR
              RIJ2 = 0.0
*       Include standard distance criterion.
              DO 45 K = 1,3
                  RIJ2 = RIJ2 + (X(K,ICOMP) - X(K,JCOMP))**2
   45         CONTINUE
              IF(RIJ2.LT.SMMIN) SMMIN=RIJ2
              IF(RIJ2.GT.SMMAX) SMMAX=RIJ2
              XMBIN = BODY(ICOMP) + BODY(JCOMP)
              PERIOD = RIJ2**0.75/DSQRT(XMBIN)
              IF(XMBIN.LT.XMMIN) XMMIN=XMBIN
              IF(XMBIN.GT.XMMAX) XMMAX=XMBIN
              IF(PERIOD.LT.TMMIN) TMMIN=PERIOD
              IF(PERIOD.GT.TMMAX) TMMAX=PERIOD
              IF (RIJ2.LT.RMIN**2) THEN
                 CALL KSPREG
              END IF
   50     CONTINUE
*     recover kstart
          KSTART = KSDUM
      END IF
*     Diagnostic
      if(rank.eq.0) write(6,61) SQRT(SMMIN)*RAU, SQRT(SMMAX)*RAU,
     &      XMMIN*ZMBAR, XMMAX*ZMBAR, TMMIN*DAYS, TMMAX*DAYS
 61   format(/,' Binary parameters: R12(min)[AU]',F15.6,
     &     ' R12(max)[AU]',F15.6,' MCM(min)[M*]',F15.6,
     &     ' MCM(max)[M*]',F15.6,' P[min][days]',1P,E15.6,
     &     ' P[max][days]',E15.6,0P)

*       Set initial neighbour list & corresponding radius.
      RS0 = RC
      call cputim(tt1)
C#ifdef GPU
      CALL FPOLY0(RS0)
      call cputim(tt2)
      if(rank.eq.0) write(6,52) (tt2-tt1)*60.0
 52   format(/,' Fpoly0 + Neighbor List Timing: ',F15.6)
      call flush(6)
c#else      
c      NNB0 = 0
c*
c      DO 40 I = IFIRST,NTOT
c          CALL NBLIST(I,RS0)
c          NNB0 = NNB0 + LIST(1,I)
c   40 CONTINUE
c*
c*       Obtain force & first derivative.
c#ifdef PARALLEL
c      CALL FPOLY1_MPI(IFIRST,NTOT,0)
c#else
c      CALL FPOLY1(IFIRST,NTOT,0)
c#endif
c      call cputim(tt2)
c      if(rank.eq.0) write(6,53) (tt2-tt1)*60.0
c 53   format(/,' Fpoly1 + Neighbor List Timing: ',F15.6)
c      call flush(6)
c*
c#endif      

      IF (KZ(40).GT.0) THEN
*     Obtain second & third force derivatives and set time-steps.
         call cputim(tt1)
#ifdef PARALLEL
         CALL FPOLY2_MPI(IFIRST,NTOT,0)
#else
         CALL FPOLY2(IFIRST,NTOT,0)
#endif
         call cputim(tt2)
         if(rank.eq.0) write(6,54) (tt2-tt1)*60.0
 54      format(/,' Fpoly2 Timing: ',F15.6)
         call flush(6)
      ELSE
*     Obtain time step only based on xdot, F and Fdot
         CALL STEPS2(IFIRST,NTOT)
      END IF
*
*     Generate perturber list for KS
      CALL KSPINIT
      call cputim(tt3) 
      if(rank.eq.0) write(6,55) (tt3-tt2)*60.0
 55   format(/,' KS Perturbers List Initialization Timing: ',F15.6)
      call flush(6)
*
*     --11/08/14 10:12-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$      DO II = IFIRST,NTOT
c$$$         write(113,*) 'II',II,'N',NAME(II),'X',X(1:3,II),
c$$$     &        'XD',XDOT(1:3,II),'F',F(1:3,II),'FD',FDOT(1:3,II),
c$$$     &        'D2',D2(1:3,II),'STEP',STEP(II),'NNB',LIST(1,II),
c$$$     &        'LIST',LIST(2:5,II)
c$$$      END DO
c$$$      call flush(113)
c$$$      DO II = 1, NPAIRS
c$$$         write(114,*) 'II',II,'H',H(II),'R',R(II),'STEP',STEP(2*II-1),
c$$$     &        'NP',LIST(1,2*II-1),'LIST',LIST(2:5,2*II-1)
c$$$      END DO
c$$$      call flush(114)
c$$$      call abort()
*     --11/08/14 10:12-lwang-end----------------------------------------*
*       Include optional regularization of primordial triples.
      IF (KZ(18).GT.1.AND.NHI0.GT.0) THEN
          KSPAIR = 1
*       Note that each KS pair will move to the top of the queue.
   60     ICOMP = 2*KSPAIR - 1
          ICM = KSPAIR + N
          RX2 = 1.0
*       Find index of closest outer component without any assumption.
          DO 70 J = IFIRST,N
              RIJ2 = 0.0
              DO 65 K = 1,3
                  RIJ2 = RIJ2 + (X(K,ICM) - X(K,J))**2
   65         CONTINUE
              IF (RIJ2.LT.RX2) THEN
                  RX2 = RIJ2
                  JCOMP = J
              END IF
   70     CONTINUE
*       Increase KS pointer for rejected separation (include possible exit).
          IF (SQRT(RX2).GT.RMIN) THEN
              KSPAIR = KSPAIR + 1
              IF (KSPAIR.GT.NPAIRS) GO TO 80
              GO TO 60
          END IF
*       Evaluate PCRIT for R0(NPAIRS) in MERGE since IMPACT is bypassed.
          CALL HISTAB(KSPAIR,JCOMP,PMIN,RSTAB)
*       Initialize the triple (constructed to be stable in HIPOP).
          IPHASE = 6
          CALL MERGE
          IF (NMERGE.LT.NHI0) THEN
              GO TO 60
          END IF
      END IF
*
*       Adjust NNBMAX (R.Sp.)
 80   NNBMAX = MIN(N/2,LMAX - 3)
      ZNBMIN = MAX(0.01*FLOAT(NNBMAX),1.0)
      ZNBMAX = 0.9*FLOAT(NNBMAX)
*
*
* @qs
*       Regularize any close planet.
      IF (nmass.GT.0) THEN
          SMMIN = 1.D30
          SMMAX = 0.D0
          ECCMIN = 1.D30
          ECCMAX = 0.D0
          IPAIR = 0
          DO 51 I = N - nmass + 1,N
c use smallMass*1.1, because body(i) = 0.999999999999 * smallMass           
              if(body(i).gt.smallMass*1.1) go to 51
*             if(nomass(i).eq.0)go to 51
c NNB: number of neighbours
              NNB = LIST(1,I)
c              print *,I,'NNB:',NNB,LIST(1:NNB+1,I)
*       Select only planets having close neighbour
              IF (NNB.GT.0) THEN
                  IPAIR = IPAIR + 1
                  JCOMP = I
                  NNB1 = NNB + 1
                  RXMIN2 = 1.D30
*       Search for host star among neighbours
                  DO 56 L = 2,NNB1
                      J = LIST(L,I)
*     --16/05/07 22:23-qishu-debug----------------------------*
***** Note:--------------------------------------------------**
              if (nomass(J).eq.1) then 
*                  print *, "WARNING!! p in neigh" 
                  
                  go to 56
              end if
*                      if (nomass(J).eq.1)go to 56 
*     --16/05/07 22:23-qishu-end------------------------------*
*       Distance determination
                      RIJ2 = 0.0
                      DO 57 K = 1,3
                          RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
   57                 CONTINUE
*       Include standard distance criterion.
                      IF (RIJ2.LT.RXMIN2) THEN
                          ICOMP = J
                          RXMIN2 = RIJ2
                      END IF
*
   56             CONTINUE
              RIJ2 = 0.0D0
              VIJ2 = 0.0D0
              RVI = 0.0D0
*       Try regularization for all binaries.
              DO 58 K = 1,3
              XRELPL(K) = X(K,ICOMP) - X(K,JCOMP)
              VRELPL(K) = XDOT(K,ICOMP) - XDOT(K,JCOMP)
              RIJ2 = RIJ2 + XRELPL(K)**2
              VIJ2 = VIJ2 + VRELPL(K)**2
              RVI = RVI + XRELPL(K)*VRELPL(K)
   58         CONTINUE
*
          XMBIN = BODY(ICOMP) + BODY(JCOMP)
          XMPROD = BODY(ICOMP)*BODY(JCOMP)
          XMUE = XMPROD/XMBIN
*     --16/05/16 17:03-qishu-debug----------------------------*
***** Note:--------------------------------------------------**
*          EBp only for planets
          EBp = VIJ2/2.D0 - XMBIN/DSQRT(RIJ2)
          SEMI = XMBIN/(2.D0*DABS(EBp))
*     --16/05/16 17:03-qishu-end------------------------------*
          ECC2 = (1.0 - DSQRT(RIJ2)/SEMI)**2 +
     &                             RVI**2/(SEMI*XMBIN)
          ECC = DSQRT(ECC2) 
          PERIOD = SEMI**1.5D0/DSQRT(XMBIN)
*     --16/05/07 15:38-qishu-debug----------------------------*
***** Note:--------------------------------------------------**
c        if (ECC .lt. 1) then
        if(rank.eq.0)then
        if (ECC .lt. 1 .or.  DSQRT(RIJ2)*206264.81 .lt. 10000) then
            print *,'------------->>>>>>>>>>'
            print *,'i,j',ICOMP,JCOMP
            print *,'body i,j:',body(ICOMP),body(JCOMP)
            print *,'X,Y,Z i,j:',X(1:3,ICOMP),X(1:3,JCOMP)
            print *,'VX,VY,VZ i,j:',XDOT(1:3,ICOMP),XDOT(1:3,JCOMP)
c            print *,'RAU',RAU
            PRINT*,'distance/AU, Semi/AU, Ecc, Period/YRS=',
     &      DSQRT(RIJ2)*206264.81, SEMI*RAU,ECC,PERIOD*YRS
        end if
        end if
*     --16/05/07 15:38-qishu-end------------------------------*
*
*
        IF(SEMI.LT.SMMIN)SMMIN=SEMI
        IF(SEMI.GT.SMMAX)SMMAX=SEMI
        IF(ECC.LT.ECCMIN)ECCMIN=ECC
        IF(ECC.GT.ECCMAX)ECCMAX=ECC
        IF(PERIOD.LT.TMMIN)TMMIN=PERIOD
        IF(PERIOD.GT.TMMAX)TMMAX=PERIOD
        EBH = -1.D10
        RMIN2 = RMIN*RMIN
        IF(RXMIN2.LT.RMIN2)CALL KSREG
        CALL FLUSH(6)
        if(rank.eq.0.AND.SEMI.LT.RMIN2)then
                  CALL FLUSH(6)
        end if
              END IF
   51     CONTINUE
*
      if(rank.eq.0)then
          PRINT*,' Min/Max Semi in AU=',SMMIN*RAU,SMMAX*RAU
          PRINT*,' Min/Max Ecc=',ECCMIN,ECCMAX
          PRINT*,' Min/Max Period in Yrs=',TMMIN*YRS,TMMAX*YRS
      end if
*
      END IF
*
      IF (NBIN0.GT.0.OR.nmass.gt.0) THEN
*
*       Check initial neighbour lists again.
          DO 95 I = IFIRST,NTOT
              IF(I.GT.N)THEN
              ICOMP = 2*IPAIR - 1
              JCOMP = 2*IPAIR
                  RS0 = RS(ICOMP)
              ELSE
                  RS0 = RC
              END IF
                  CALL NBLIST(I,RS0)
   95     CONTINUE
*
      END IF
*       Check the average neighbour number.
      NNB0 = 0
      DO 86 I = IFIRST,NTOT
          NNB0 = NNB0 + LIST(1,I)
 86    CONTINUE
      ZNB = FLOAT(NNB0)/FLOAT(N)
*      nnbave = ZNB
      IF (ZNB.LT.0.25*ZNBMAX) THEN
          if(rank.eq.0)WRITE (6,90)  ZNB
   90     FORMAT (/,12X,'WARNING!   SMALL NEIGHBOUR NUMBERS   <NNB> =',
     &                                                             F5.1)
      END IF
*
*       Check option for writing the initial conditions on unit 10.
      IF (KZ(22).EQ.1.and.rank.eq.0) THEN
          DO 85 I = 1,N
             WRITE (10,*) BODY(I),(X(K,I),K=1,3),(XDOT(K,I),K=1,3)
   85     CONTINUE
*  100 FORMAT(1X,1P,7(1X,E26.10))
      END IF
      CALL FLUSH(10)
      CLOSE(10)
*
      RETURN
*
      END
