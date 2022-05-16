      subroutine DFUNC(p,xi)

c************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'

      double precision p(nummax*3)
      double precision xi(nummax*3)
      
      do 10 i=1,nummax

         xi(i)=((2.0d0*xk1(1)*xmu(i))
     x        -elecx(i))
         xi(i+nummax)=((2.0d0*xk1(1)*ymu(i))
     x               -elecy(i))
         xi(i+nummax+nummax)=((2.0d0*xk1(1)*zmu(i))
     x                      -elecz(i))

   10 continue

      do 20 i=1,ncation

         xi(i+nanion)=0.0d0
         xi(i+nanion+nummax)=0.0d0
         xi(i+nanion+nummax+nummax)=0.0d0

   20 continue

      return

      end
