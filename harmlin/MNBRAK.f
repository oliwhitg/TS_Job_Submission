      SUBROUTINE MNBRAK(AXX,BXX,CX,FA,FB,FC)

c************************************************************

      implicit double precision (a-h,o-z)

      include 'common.inc'

      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.E-20)
      FA=F1DIM(AXX)
      FB=F1DIM(BXX)

      IF(FB.GT.FA)THEN
        DUM=AXX
        AXX=BXX
        BXX=DUM
        DUM=FB
        FB=FA
        FA=DUM
      ENDIF
      CX=BXX+GOLD*(BXX-AXX)
      FC=F1DIM(CX)
1     IF(FB.GE.FC)THEN
        R=(BXX-AXX)*(FB-FC)
        QQ=(BXX-CX)*(FB-FA)
        U=BXX-((BXX-CX)*QQ-(BXX-AXX)*R)
     x   /(2.*SIGN(MAX(ABS(QQ-R),TINY),QQ-R))
        ULIM=BXX+GLIMIT*(CX-BXX)
        IF((BXX-U)*(U-CX).GT.0.)THEN
          FU=F1DIM(U)
          IF(FU.LT.FC)THEN
            AXX=BXX
            FA=FB
            BXX=U
            FB=FU
            GO TO 1
          ELSE IF(FU.GT.FB)THEN
            CX=U
            FC=FU
            GO TO 1
          ENDIF
          U=CX+GOLD*(CX-BXX)
          FU=F1DIM(U)
        ELSE IF((CX-U)*(U-ULIM).GT.0.)THEN
          FU=F1DIM(U)
          IF(FU.LT.FC)THEN
            BXX=CX
            CX=U
            U=CX+GOLD*(CX-BXX)
            FB=FC
            FC=FU
            FU=F1DIM(U)
          ENDIF
        ELSE IF((U-ULIM)*(ULIM-CX).GE.0.)THEN
          U=ULIM
          FU=F1DIM(U)
        ELSE
          U=CX+GOLD*(CX-BXX)
          FU=F1DIM(U)
        ENDIF
        AXX=BXX
        BXX=CX
        CX=U
        FA=FB
        FB=FC
        FC=FU
        GO TO 1
      ENDIF
      RETURN
      END
