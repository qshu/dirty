      SUBROUTINE CALC_REG_ON_GPU(IREG,istart,iend)
*
*
*     Calculate the REG force on local GPU memory.
*     -------------------------------------------
*
      PARAMETER (maxthr=1024)
      INCLUDE 'common6.h'
      INCLUDE 'timing.h'
      INCLUDE 'omp_lib.h'
      COMMON/POTDEN/  RHO(NMAX),XNDBL(NMAX),PHIDBL(NMAX)
      
      INTEGER IREG(NMAX), istart, iend
      INTEGER LISTGP(LMAX,maxthr),NOFL,M_FLAG
#ifdef GPU
      REAL*8 GPU_RS(maxthr),GPU_DTR(maxthr)
      REAL*8 GPU_POT(maxthr)
#endif
      REAL*8 x_i(3,maxthr), v_i(3,maxthr)
      REAL*8 GPU_A(3,maxthr), GPU_JERK(3,maxthr)

      INTEGER qsi,qsj, qsflag
*      REAL*8 FI_HOST(3,maxthr),FID_HOST(3,maxthr)
*     If it's in parallel loop, flag_p = .true.
*      LOGICAL FLAG_P

*     --16/05/14 00:28-qishu-debug----------------------------*
***** Note:--------------------------------------------------**
c      initial GPU_A(:,:) = 0
      GPU_A(:,:) = 0.
*     --16/05/14 00:28-qishu-end------------------------------*
      M_FLAG = 0
      NB_FLAG = 0
*     Use mass ratio as neighbor member criterion
      IF(KZ(39).EQ.2.OR.KZ(39).EQ.3) M_FLAG = 1
      IF(KZ(39).EQ.0.OR.KZ(39).EQ.2) NB_FLAG = 1
      
      NOFL = 0
*     calculate the pot, a & jerk on the GPU
      DO L=istart,iend,maxthr

         ni = maxthr
         IF( (L+ni).GT.iend ) ni = iend - L + 1

#ifdef GPU
*     fill the GPU "working" arrays
!$omp parallel do private(ii,idi)
         DO ii=1,ni
            idi = IREG(L+ii-1)
*    Obtain irregular & regular force and determine current neighbours.
*    Take volume between inner and outer radius equal to basic sphere.
            IF(M_FLAG.EQ.1) THEN
*     Use mass ratio as neighbor member criterion
               GPU_RS(ii) = RS(idi)**2/BODYM
            ELSE
               GPU_RS(ii) = RS(idi)**2
            END IF
            GPU_DTR(ii) = STEPR(idi)
            x_i(1:3,ii) = X(1:3,idi)
            v_i(1:3,ii) = XDOT(1:3,idi)
         END DO
!$omp end parallel do

*     GPU first call with first values of GPU_POT, GPU_A & GPU_JERK
c
*     --16/05/13 18:32-qishu-debug----------------------------*
***** Note:--------------------------------------------------**
c         GPU_A maybe nan !!!
 550     call cputim(tt53)
c         print *,"ni = ",ni,'size:',size(x_i),size(v_i)
        do qsi = 1, n-1
        do qsj = qsi+1, n
        if (x(1,qsi) .eq. x(1,qsj)) stop "smae positive here!!!"
        if (x(2,qsi) .eq. x(2,qsj)) stop "smae positive here!!!"
        if (x(3,qsi) .eq. x(3,qsj)) stop "smae positive here!!!"
        end do
        end do
c            print *, x_i(qsi)
c        stop "end test"

c         print *,ni,GPU_RS,GPU_DTR,x_i,v_i,GPU_A,GPU_JERK,
c     &   GPU_POT,lmax,nnbmax,LISTGP,M_FLAG
c        end if

c        if (abs (GPU_A(1,1) - (-0.55306506536636946)) .lt. 0.0005) then
        do qsi = 1,n
        if (x_i(1,1) .eq. x(1,qsi) .and. x_i(2,1) .eq. x(2,qsi)) then
        print *,'>>>>>>>>>>>>>>',ni
           print *,'x_i(1) in x(i):',qsi
        end if
        end do
         do qsi = 1, size(GPU_RS)
         if (isnan(GPU_RS(qsi))) print *,'NAN GPU_RS:',qsi
         end do
         do qsi = 1, size(GPU_A)/3
         do qsj = 1,3
         if (isnan(GPU_A(qsj,qsi))) print *,'NAN GPU_A:',qsi,qsj
         end do
         end do
c             print *, ni, GPU_A(:,1:2)
c             print *,GPU_A
c             ,GPU_JERK,
c     &        GPU_POT,lmax,nnbmax,LISTGP,M_FLAG
c         end if
c         print *,'>>ni:',ni,'438:', GPU_A(:,438)
c         print *, abs (GPU_A(1,1) - (-0.55306506536636946))
c             print *, ni, GPU_A(:,1:2)

c        print *,'>>>>>>>>>>>>>>',ni
c         do qsi = 1, size(GPU_A)/3
c         do qsj = 1,3
c         if (isnan(GPU_A(qsj,qsi))) then
c             print *,'nan detected before gpunb_regf'
c             print *,qsi,qsj
c             print *, ni, GPU_A(:,qsi)
c         end if
c
c         end do
c         end do
         CALL gpunb_regf(ni,GPU_RS,GPU_DTR,x_i,v_i,GPU_A,GPU_JERK,
     &        GPU_POT,lmax,nnbmax,LISTGP,M_FLAG)
         do qsi = 1, size(GPU_A)/3
         do qsj = 1,3
         if (isnan(GPU_A(qsj,qsi))) then
             print *,qsi,qsj
             print *, ni, GPU_A(:,qsi)
             stop "NAN"
         end if

         end do
         end do
c         print *,'<<ni:',ni,'438:', GPU_A(:,438)

*     --16/05/13 18:32-qishu-end------------------------------*
         call cputim(tt54)
         ttgpu = ttgpu + (tt54-tt53)*60.0

*     Save new neighbour list on intermediate vector for regint
         DO II = 1,NI
            idi = IREG(L+ii-1)
            NNB = LISTGP(1,ii)

*     Warning if neighbour list overflows.
            IF (NNB.LT.0) THEN
               RI = sqrt(x_i(1,II)**2 + x_i(2,II)**2 + x_i(3,II)**2)
               WRITE (41,56) idi, NAME(idi), LIST(1,idi), -NNB, RS(idi),
     &              RI, time+toff
 56            FORMAT (' OVERFLOW! I',I8,'  NAME',I8,'  NBPRE',I6,
     &              '  NBNEW',I6,'  RNB[NB]',F10.4,
     &              '  RI[NB]',F10.4,'  Time[NB]',1P,D12.5,0P)
*     Define all particles as neighbours for overflow for full force loop in regint.
*     LIST_G(1,idi) = NNB
               NOFL = NOFL + 1
               IF(NB_FLAG.EQ.1) THEN
*     For smooth density system, use NNBOPT/NNB to estimate the neighbor radius
                  IF(-NNB.GT.NNBOPT) THEN
                     RS(idi)=(-FLOAT(NNBOPT)/FLOAT(NNB))**0.333*RS(idi)
                  ELSE
                     RS(idi)=(FLOAT(NNBOPT)/FLOAT(NNBMAX))**0.4*RS(idi)
                  END IF
               ELSE
                  RS(idi) = 0.9*RS(idi)
               END IF
               IF(M_FLAG.EQ.1) THEN
*     Use mass ratio as neighbor member criterion
                  GPU_RS(ii) = RS(idi)**2/BODYM
               ELSE
                  GPU_RS(ii) = RS(idi)**2
               END IF
            END IF
            IF (NOFL.GT.0) THEN
               NOFL = 0
               GO TO 550
            END IF
            L1 = 1
            DO LL = 2,NNB+1
*     Note GPU address starts from 0 (hence add IFIRST to neighbour list).
               ITEMP = LISTGP(LL,II) + IFIRST
               IF(ITEMP.NE.IDI) THEN
                  L1 = L1 + 1
                  LISTGP(L1,II) = ITEMP
               END IF
            END DO
            LISTGP(1,II) = L1-1
            PHIDBL(idi) = GPU_POT(ii)
C$$$#ifdef SIMD
C$$$            CALL IRR_SIMD_SET_LIST(IDI,LISTGP(1,II))
C$$$#endif
         END DO
#else
*     OpemMP version of regular force, fdot, potential and neighbor list
         call cputim(tt53)
         CALL regint(L,NI,IREG,IFIRST,NTOT,X,XDOT,BODY,STEPR,BODYM,
     &        M_FLAG,NB_FLAG,NNBOPT,RS,GPU_A,GPU_JERK,LISTGP,PHIDBL)
         call cputim(tt54)
         ttgpu = ttgpu + (tt54-tt53)*60.0
#endif
         CALL FLUSH(41)

C$$$#ifdef SIMD
*       Evaluate current irregular forces by vector procedure.
C$$$         CALL IRR_SIMD_FIRR_VEC(TIME,NI,IREG(L),FI_HOST,FID_HOST)
*     --04/26/14 11:13-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$         write(6,*) 'L L',IREG(L),'N',NAME(IREG(L)),'FI',FI_HOST(1:3,1),
c$$$     &        'FD',FID_HOST(1:3,1)
c$$$         call flush(6)
*     --04/26/14 11:13-lwang-end----------------------------------------*
C$$$#else
*     Evaluate current irregular forces
C$$$         CALL regint_gpu(ni,L,ireg,FI_HOST,FID_HOST,LISTGP)
*     --04/26/14 11:13-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$         write(6,*) 'H L',IREG(L),'N',NAME(IREG(L)),'FI',FI_HOST(1:3,1),
c$$$     &        'FD',FID_HOST(1:3,1)
c$$$         call flush(6)
c$$$         STOP
*     --04/26/14 11:13-lwang-end----------------------------------------*
C$$$#endif

*     Correction
*!$omp parallel do if(ni.GE.ithread) private(II,I)
         DO II = 1, ni
            I = IREG(L+II-1)
#ifndef GPU
            x_i(1:3,II) = x(1:3,I)
            v_i(1:3,II) = xdot(1:3,I)
#endif
*     --04/19/14 12:44-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
c$$$            IF(NAME(I).EQ.8317.AND.TIME.EQ.1.32592773437500000) THEN
c$$$               write(102+rank,*) 'I',I,'N',NAME(I),'FR',GPU_A(1,II),
c$$$     &              'FRD',GPU_JERK(1,II),'FI',FI_HOST(1,II),
c$$$     &              'FID',FID_HOST(1,II),'T',TIME
c$$$               DO III = IFIRST,NTOT
c$$$                  write(102+rank,*) 'III',III,'N',NAME(III),'X',X(1,III)
c$$$     &                 ,'XD',XDOT(1,III),'T',TIME,'F',F(1,III),
c$$$     &                 'FD',FDOT(1,III),'X0',X0(1,III),'X0D',
c$$$     ^                 X0DOT(1,III)
c$$$               END DO
c$$$               call flush(102+rank)
c$$$            END IF
*     --04/19/14 12:44-lwang-end----------------------------------------*
            CALL regcor_gpu(I,x_i(1,II),v_i(1,II),GPU_A(1,II),
     &           GPU_JERK(1,II),LISTGP(1,II))
c$$$#ifdef SIMD
*     Only when Neighbor list is modified, update it for AVX/SSE library
C$$$            IF(ICFLAG.GT.0.AND..NOT.FLAG_P) THEN
c$$$            CALL IRR_SIMD_SET_LIST(I,LIST(1,I))
C$$$            END IF
c$$$#endif
         END DO
*!$omp end parallel do
      END DO

      RETURN
      
      END
