      SUBROUTINE FRPRMN(P,N,ITER,FRET)

c************************************************************
      implicit double precision (a-h,o-z)

      include 'common.inc'

      PARAMETER (NMAX=3*nummax,ITMAX=200,EPS=1.E-10)
      DIMENSION P(N),G(NMAX),HHH(NMAX),XI(NMAX)

      call ener

      FP=FUNC(P)
      CALL DFUNC(P,XI)
      DO 11 J=1,N
        G(J)=-XI(J)
        HHH(J)=G(J)
        XI(J)=HHH(J)
11    CONTINUE
      DO 14 ITS=1,ITMAX
        ITER=ITS
        CALL LINMIN(P,XI,N,FRET)

        IF(2.*ABS(FRET-FP).LE.FTOL*(ABS(FRET)+ABS(FP)+EPS))RETURN

        do 10 i=1,nummax

           xmu(i)=p(i)
           ymu(i)=p(i+nummax)
           zmu(i)=p(i+nummax+nummax)

   10   continue

      write(6,*)'frprmn1',engpetot
        call ener
      write(6,*)'frprmn2',engpetot

        FP=FUNC(P)
        CALL DFUNC(P,XI)
        GG=0.
        DGG=0.
        DO 12 J=1,N
          GG=GG+G(J)**2
C         DGG=DGG+XI(J)**2
          DGG=DGG+(XI(J)+G(J))*XI(J)
12      CONTINUE
        IF(GG.EQ.0.)RETURN
        GAM=DGG/GG
        DO 13 J=1,N
          G(J)=-XI(J)
          HHH(J)=G(J)+GAM*HHH(J)
          XI(J)=HHH(J)
13      CONTINUE
14    CONTINUE
      PAUSE 'FRPR maximum iterations exceeded'
      RETURN
      END
