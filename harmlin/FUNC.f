      function FUNC(p)

c************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'

      double precision p(nummax*3)

      dipsq=0.0d0
      engptot=0.0d0
      do 10 i=1,nummax*3

         dipsq=dipsq+(p(i)*p(i))

  10  continue

cc      do 11 i=1,num

cc         engptot=engptot-(p(i)*elecx(i)
cc     x                    +p(i+nummax)*elecy(i)
cc     x                    +p(i+nummax+nummax)*elecz(i))

cc   11 continue

      FUNC=engpetot+xk1(1)*dipsq
cc      FUNC=engptot+xk1*dipsq
cc      write(6,*)'bbb',engptot,FUNC
cc      write(6,*)elecx(1),elecy(1),elecz(1)

      return

      end
