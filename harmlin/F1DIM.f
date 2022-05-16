      FUNCTION F1DIM(XXX)

c************************************************************

      implicit double precision (a-h,o-z)

      include 'common.inc'

      PARAMETER (NMAX=nummax*3)
      COMMON /F1COM/ NCOM,PCOM(NMAX),XICOM(NMAX)
      DIMENSION XT(NMAX)
      DO 11 J=1,NCOM
        XT(J)=PCOM(J)+XXX*XICOM(J)
11    CONTINUE

      do 10 i=1,nummax

         xmu(i)=xt(i)
         ymu(i)=xt(i+nummax)
         zmu(i)=xt(i+nummax+nummax)

   10 continue

      call ener

      F1DIM=FUNC(XT)
      RETURN
      END
