
      subroutine shake
      implicit double precision (a-h,o-z)

      include 'common.inc'

      double precision numer,detmh,lambda,denom
      do 5 it=1,10

c         do 10 j=1,3
c            do 20 k=1,3
c               h(j,k)=h3(j,k)
c 20         continue
c 10      continue

         call matinv(h3,hi3,detmh)

c         do 21 j=1,3
c            do 22 k=1,3
c               hi3(j,k)=hi(j,k)
c 22         continue
c 21      continue

         do 25 j=1,3
            do 27 k=1,3
               hih(j,k)=0.0d0
               do 30 l=1,3
                  hih(j,k)=hih(j,k)+(hi3(j,l)*h2(l,k))
 30            continue
 27         continue
 25      continue

         numer=detmh-1.0d0
         denom=detmh*(hih(1,1)+hih(2,2)+hih(3,3))
         lambda=-numer/denom

         do 40 j=1,3
            do 50 k=1,3

               h3(j,k)=h3(j,k)+(lambda*h2(j,k))

 50         continue
 40      continue

 5    continue

      return

      end
