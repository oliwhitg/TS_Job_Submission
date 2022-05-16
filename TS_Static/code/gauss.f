      function gauss(idum)

c************************************************************
c From Allen & Tildesley but originally from Knuth - `The
c art of computer programming', Addison-Wesley, Reading.
c Produces a normal (Gaussian) distribution with zero mean
c and unit variance.
c************************************************************

      implicit double precision(a-h,o-z)

      common/random/dummy

      integer bin(100)

      a1=3.949846138
      a3=0.252408784
      a5=0.076542912
      a7=0.008355968
      a9=0.029899776

      x=0.0d0

      do 20 j=1,12
         x=x+ran1(dummy)
   20 continue

      r=(x-6.0d0)/4.0d0
      rsq=r*r
      gauss=((((a9*rsq+a7)*rsq+a5)*rsq+a3)*rsq+a1)*r

      return

      end
