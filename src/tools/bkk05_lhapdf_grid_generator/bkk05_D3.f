c ... bkk05_D3.f (IJS)
c ... 26.08.2005

c ... Thresholds at Heavy Quark Masses
c ... LO and NLO OPAL FFs for D^*, D^0, D^+, D_s, Lambda_c

c ... No initialization needed anymore. Thanks to Markos Maniatis

********************************************************************
*   Returns Fragmentation Functions D_i(x,Q^2) (NOT x*FFs !)       * 
*   for D^*,D^0,D^+,D_s and Lambda_c                               *
*                                                                  *
*   NOTE: D^* = D^*+ + D^*- ; D^+ + D^- ; D_s^+ + D_s^- and        *
*         Lambda_c = Lambda_c^+ + Lambda_c^-                       *
*                                                                  *
*   INPUT:                                                         *
* ... NO = 1 (LO), 2 (NLO)                                         *
* ... NH = 4 KK05 D^* OPAL, 5 KK05 D^0, 6 KK05 D^+, 7 KK05 D_s,    * 
*          8 KK05 Lambda_c                                         *
*                                                                  *
*   OUTPUT: DH(0:10) <-> D_i(x,Q^2) (i=0,...,10)                   *
* ... DH: 0    1    2    3    4    5    6    7    8    9    10     *
* ...     g    u  ubar   d  dbar   s  sbar   c  cbar   b   bbar    *
*                                                                  *
* ... Note that: D_u=D_d=D_s; D_i = D_ibar                         *
*                                                                  * 
*            x  = frag. variable  (between  1.E-4  and  1)         *
*            Q2 = scale in GeV**2 (between  2.25 and  1.E6)        *
*                                                                  *
*   last change: 26.08.2005                                        *
********************************************************************
      subroutine bkk_05_D3(NO,NH,x,Q2,DH)
      implicit none
      integer NO,NH   ! Input
      real*8 X,Q2     ! Input
      real*8 DH(0:10) ! Output

      integer NOold,NHold  
      SAVE NOold, NHold
c      
      integer m,n,i
      integer IQ,IX
      integer npart,nx,nq,narg
      parameter(NPART=4, NX=131, NQ=29, NARG=2)
      double precision XB0,XB1
      double precision PARTON(0:NPART-1,NQ,NX-1), QS(NQ), XB(NX)
      integer NA(NARG)
      double precision XT(NARG), ARRF(NX+NQ) 
      double precision XDG(NX,NQ),XDU(NX,NQ)!,XDD(NX,NQ),XDS(NX,NQ)
      double precision XDC(NX,NQ),XDB(NX,NQ)
      character*3 chno
      character*4 chnh
      integer fragini
      SAVE fragini
c      COMMON / inifrag / fragini
      SAVE XDU, XDG, XDC, XDB, NA, ARRF
      double precision FINT
      double precision XBG,XBU,XBC,XBB
      double precision xg,xu,xc,xbot
      double precision rg,ru,rc,rb
      double precision rg1,ru1,rc1,rb1
      SAVE rg, ru, rc, rb, rg1, ru1, rc1, rb1

      double precision XB1G,XB1U,XB1C,XB1B
      double precision xg1,xu1,xc1,xbot1

c ... NQ=29
*...X AND Q**2 VALUES OF THE GRID :
       DATA QS / 2.25d0, 3d0, 5d0, 7d0, 10D0, 16D0, 25D0, 40D0, 
     #           64D0, 81d0,
     #           1.0D2, 1.5D2, 2.0D2, 2.5D2, 5.0d2,
     #           1.0D3, 1.8D3, 3.2D3, 5.7D3, 8315.18d0, 
     #           1.0D4, 1.8D4, 3.2D4, 5.7D4, 1.0d5,
     #           1.8D5, 3.2D5, 5.7D5, 1.0d6/

c ... NX=131
       DATA XB / 1.0D-4,1.4D-4,2.0D-4,3.0D-4,4.5D-4,6.7D-4,8.0d-4,
     #           1.0D-3,1.4D-3,2.0D-3,3.0D-3,4.5D-3,6.7D-3,8.0d-3,
     #           1.0D-2, 1.4D-2, 2.0D-2, 3.0D-2, 4.0d-2,
     #0.046d0,0.050d0,0.054d0,0.059d0,0.064d0,0.068d0,
     #0.072d0,0.076d0,0.080d0,0.085d0,0.09d0,0.095d0,0.1d0,0.105d0,
     #0.110d0,0.117d0,0.125d0,0.130d0,0.135d0,0.140d0,0.15d0,0.16d0,
     #0.17d0,0.18d0,0.19d0,0.2d0,0.21d0,0.22d0,0.23d0,0.24d0,
     #0.246d0,0.252d0,0.260d0,0.2660,0.272d0,0.280d0,0.29d0,0.3d0,
     #0.31d0,0.32d0,0.33d0,0.34d0,0.35d0,0.36d0,0.37d0,
     #0.38d0,0.39d0,0.40d0,0.41d0,0.42d0,0.43d0,
     #0.44d0,0.45d0,0.46d0,0.47d0,0.48d0,0.49d0,0.5d0,
     #0.51d0,0.52d0,0.53d0,0.54d0,0.55d0,0.56d0,
     #0.57d0,0.58d0,0.59d0,0.6d0,0.61d0,0.62d0,0.63d0,
     #0.64d0,0.65d0,0.66d0,0.67d0,0.68d0,0.69d0,
     #0.70d0,0.71d0,0.72d0,0.73d0,0.74d0,0.75d0,0.76d0,
     #0.77d0,0.78d0,0.79d0,0.80d0,0.81d0,0.82d0,
     #0.83d0,0.84d0,0.85d0,0.86d0,0.87d0,0.888d0,0.89d0,
     #0.90d0,0.91d0,0.92d0,0.93d0,0.94d0,0.95d0,0.96d0,
     #0.965d0,0.97d0,0.975d0,0.98d0,0.985d0,0.99d0,0.995d0,0.9999d0 /


*...CHECK OF X AND Q2 VALUES : 
       IF ( (X.LT.0.99E-4) .OR. (X.GT.1.0) ) THEN
           WRITE(6,91) 
  91       FORMAT (2X,'PARTON INTERPOLATION: x OUT OF RANGE')
       ENDIF
       IF ( (Q2.LT.2.25d0) .OR. (Q2.GT.1d4) ) THEN
           WRITE(6,92) 
  92       FORMAT (2X,'PARTON INTERPOLATION: Q2 OUT OF RANGE')
       ENDIF


*
*...INITIALIZATION :
*    SELECTION AND READING OF THE GRID : 
c check whether parameters NO or NH changed
c or this is a the first call (after the firtst call fragini is set to
c the arbitrary number 12345)
c      IF (fragini.NE.0) GOTO 16
      IF (fragini.EQ.12345 .AND. NO.EQ.NOold .AND. NH.EQ.NHold) GOTO 16

      if ((NO.EQ.1) .and. (NH.EQ.1)) then ! ALEPH LO (BKK98-D^*)
         stop
               chno='no1'
               chnh='nh01'
      elseif ((NO.EQ.1) .and. (NH.EQ.2)) then ! OPAL LO (BKK98-D^*)         
         stop
               chno='no1'
               chnh='nh02'

      elseif ((NO.EQ.1) .and. (NH.EQ.3)) then ! (KK05-D^* LO ALEPH)
         stop
               chno='no1'
               chnh='nh03'

      elseif ((NO.EQ.1) .and. (NH.EQ.4)) then ! (KK05-D^* LO OPAL)
               chno='no1'
               chnh='nh04'

               rg = -1d0
               ru = -2d0
               rc = 1d0
               rb = 1d0

               rg1 = 0d0
               ru1 = 2d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.1) .and. (NH.EQ.5)) then ! (KK05-D^0 LO)
               chno='no1'
               chnh='nh05'

               rg = -2d0
               ru = -2d0
               rc = 0d0
               rb = 1d0

               rg1 = 1d0
               ru1 = 2d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.1) .and. (NH.EQ.6)) then ! (KK05-D^+ LO)
               chno='no1'
               chnh='nh06'

               rg = -1d0
               ru = -2d0
               rc = 0d0
               rb = 1d0

               rg1 = 0d0
               ru1 = 0d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.1) .and. (NH.EQ.7)) then ! (KK05-D_s LO)
               chno='no1'
               chnh='nh07'

               rg = -2d0
               ru = -2d0
               rc = 0d0
               rb = 1d0

               rg1 = 1d0
               ru1 = 2d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.1) .and. (NH.EQ.8)) then ! (KK05-Lambda_c LO)
               chno='no1'
               chnh='nh08'

               rg = -2d0
               ru = -1d0
               rc = 0d0
               rb = 1d0

               rg1 = 1d0
               ru1 = 1d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.2) .and. (NH.EQ.1)) then ! ALEPH NLO (BKK98-D^*)
         stop
               chno='no2'
               chnh='nh01'

      elseif ((NO.EQ.2) .and. (NH.EQ.2)) then ! OPAL NLO (BKK98-D^*)
         stop
               chno='no2'
               chnh='nh02'

      elseif ((NO.EQ.2) .and. (NH.EQ.3)) then ! (KK05-D^* NLO ALEPH)
         stop
               chno='no2'
               chnh='nh03'

      elseif ((NO.EQ.2) .and. (NH.EQ.4)) then ! (KK05-D^* NLO OPAL)
               chno='no2'
               chnh='nh04'

               rg = -1d0
               ru = -2d0
               rc = 1d0
               rb = 1d0

               rg1 = 0d0
               ru1 = 2d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.2) .and. (NH.EQ.5)) then ! (KK05-D^0 NLO)
               chno='no2'
               chnh='nh05'

               rg = -2d0
               ru = -2d0
               rc = 0d0
               rb = 1d0

               rg1 = 1d0
               ru1 = 2d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.2) .and. (NH.EQ.6)) then ! (KK05-D^+ NLO)
               chno='no2'
               chnh='nh06'

               rg = -1d0
               ru = -2d0
               rc = 0d0
               rb = 1d0

               rg1 = 0d0
               ru1 = 0d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.2) .and. (NH.EQ.7)) then ! (KK05-D_S NLO)
               chno='no2'
               chnh='nh07'

               rg = -2d0
               ru = -2d0
               rc = 0d0
               rb = 1d0

               rg1 = 1d0
               ru1 = 2d0
               rc1 = 2d0
               rb1 = 3d0

      elseif ((NO.EQ.2) .and. (NH.EQ.8)) then ! (KK05-Lambda_c NLO)
               chno='no2'
               chnh='nh08'

               rg = -2d0
               ru = -1d0
               rc = 0d0
               rb = 1d0

               rg1 = 1d0
               ru1 = 1d0
               rc1 = 2d0
               rb1 = 3d0

      else
         print*,'NO =',NO,' and ','NH =',NH,' not supported'
         stop
      end if    

      fragini = 12345
      NOold = NO
      NHold = NH
c
      OPEN(UNIT=8,FILE='bkk3_d_'//chno//chnh//'.grid',STATUS='old')

            DO M = 1, NX-1  ! Loop over x
               DO N = 1, NQ ! Loop over Q^2

               READ(8,90) PARTON(0,N,M), PARTON(1,N,M), PARTON(2,N,M), 
     #                    PARTON(3,N,M)

               enddo ! N-loop
            enddo ! M-loop
  90  FORMAT (6(1PE13.5,1X))

      CLOSE(8)
*
*....ARRAYS FOR THE INTERPOLATION SUBROUTINE :
      DO 10 IQ = 1, NQ
      DO 20 IX = 1, NX-1
c      DO 20 IX = 1, NX
        XB0 = XB(IX) 
        XB1 = 1.-XB0

        XBG = XB0**rg
        XBU = XB0**ru
        XBC = XB0**rc
        XBB = XB0**rb

        XB1G = XB1**rg1
        XB1U = XB1**ru1
        XB1C = XB1**rc1
        XB1B = XB1**rb1

        XDG(IX,IQ) = PARTON(0,IQ,IX) / (XBG * XB1G)
        XDU(IX,IQ) = PARTON(1,IQ,IX) / (XBU * XB1U)
        XDC(IX,IQ) = PARTON(2,IQ,IX) / (XBC * XB1C)
        XDB(IX,IQ) = PARTON(3,IQ,IX) / (XBB * XB1B)
  20  CONTINUE
        XDG(NX,IQ) = 0.E0
        XDU(NX,IQ) = 0.E0
        XDC(NX,IQ) = 0.E0
        XDB(NX,IQ) = 0.E0
  10  CONTINUE  
      NA(1) = NX
      NA(2) = NQ
      DO 30 IX = 1, NX
        ARRF(IX) = DLOG(XB(IX))
  30  CONTINUE
      DO 40 IQ = 1, NQ
        ARRF(NX+IQ) = DLOG(QS(IQ))
  40  CONTINUE
  16  CONTINUE
*
*...INTERPOLATION :
      XT(1) = DLOG(X)
      XT(2) = DLOG(Q2)

      xg = x**rg
      xu = x**ru
      xc = x**rc
      xbot = x**rb

      xg1 = (1.-x)**rg1
      xu1 = (1.-x)**ru1
      xc1 = (1.-x)**rc1
      xbot1 = (1.-x)**rb1

c ... flavour symmetric: D_u=D_d=D_s
c ... charge symmetry: D_i = D_anti-i

      DH(0) = FINT(NARG,XT,NA,ARRF,XDG) * xg * xg1      ! D_g
      DH(1) = FINT(NARG,XT,NA,ARRF,XDU) * xu * xu1      ! D_u
      DH(2) = DH(1)                                     ! D_ubar
      DH(3) = DH(1)                                     ! D_d
      DH(4) = DH(3)                                     ! D_dbar
      DH(5) = DH(1)                                     ! D_s
      DH(6) = DH(5)                                     ! D_sbar
      DH(7) = FINT(NARG,XT,NA,ARRF,XDC) * xc * xc1      ! D_c
      DH(8) = DH(7)                                     ! D_cbar
      DH(9) = FINT(NARG,XT,NA,ARRF,XDB) * xbot * xbot1  ! D_b
      DH(10)= DH(9)                                     ! D_bbar

c ... to avoid crash for scales smaller than q2 data/         
      if (Q2.lt.2.25d0) then
c         write(6,*) 'q**2 must be > 2.25 GeV^2',qp2
         do i=0,10
            DH(i)=0d0
         enddo   
      return
      endif

      if (x.eq.1.d0) then
         do i=0,10
            DH(i)=0d0
         enddo   
      return
      endif

 60   RETURN
      END
*
*
*...CERN LIBRARY ROUTINE E104 (INTERPOLATION) :
*
      FUNCTION FINT(NARG,ARG,NENT,ENT,TABLE)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION ARG(5),NENT(5),ENT(10),TABLE(10)
      DIMENSION D(5),NCOMB(5),IENT(5)

      KD=1
      M=1
      JA=1
         DO 5 I=1,NARG
      NCOMB(I)=1
      JB=JA-1+NENT(I)
         DO 2 J=JA,JB
      IF (ARG(I).LE.ENT(J)) GO TO 3
    2 CONTINUE
      J=JB
    3 IF (J.NE.JA) GO TO 4
      J=J+1
    4 JR=J-1
      D(I)=(ENT(J)-ARG(I))/(ENT(J)-ENT(JR))
      IENT(I)=J-JA
      KD=KD+IENT(I)*M
      M=M*NENT(I)
    5 JA=JB+1
      FINT=0.
   10 FAC=1.
      IADR=KD
      IFADR=1
         DO 15 I=1,NARG
      IF (NCOMB(I).EQ.0) GO TO 12
      FAC=FAC*(1.-D(I))
      GO TO 15
   12 FAC=FAC*D(I)
      IADR=IADR-IFADR
   15 IFADR=IFADR*NENT(I)
      FINT=FINT+FAC*TABLE(IADR)
      IL=NARG
   40 IF (NCOMB(IL).EQ.0) GO TO 80
      NCOMB(IL)=0
      IF (IL.EQ.NARG) GO TO 10
      IL=IL+1
         DO 50  K=IL,NARG
   50 NCOMB(K)=1
      GO TO 10
   80 IL=IL-1
      IF(IL.NE.0) GO TO 40

      RETURN
      END
