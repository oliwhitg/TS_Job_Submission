      SUBROUTINE LINMIN(P,XI,N,FRET)

c************************************************************

      implicit double precision (a-h,o-z)

      include 'common.inc'

      PARAMETER (NMAX=nummax*3)

      EXTERNAL F1DIM,DF1DIM
      DIMENSION P(N),XI(N)
      COMMON /F1COM/ NCOM,PCOM(NMAX),XICOM(NMAX)

      NCOM=N
cc      write(6,*)'linmin start'
cc      do 111 i=1,nanion
cc         write(6,*)'xi1 in linmin 1',xi(i)
cc  111 continue
      DO 11 J=1,N
        PCOM(J)=P(J)
        XICOM(J)=XI(J)
11    CONTINUE
      AXX=0.
      XX=1.
      BXX=2.
      CALL MNBRAK(AXX,XX,BXX,FA,FX,FB)
cc      write(6,*)'ax,xx....',ax,xx,bx,fa,fx,fb
cc      FRET=BRENT(AX,XX,BX,F1DIM,TOL,XMIN)
      FRET=DBRENT(AXX,XX,BXX,F1DIM,DF1DIM,TOL,XMIN)
     
      do 20 i=1,ncation

         xi(i+nanion)=0.0d0
         xi(i+nanion+nummax)=0.0d0
         xi(i+nanion+nummax+nummax)=0.0d0

   20 continue

      DO 12 J=1,N
        XI(J)=XMIN*XI(J)
cc        write(99,*)j,'xi in linmin',xi(j)
        P(J)=P(J)+XI(J)
12    CONTINUE

cc      write(6,*)'linmin end'

      RETURN
      END
