      subroutine rattle
      implicit double precision (a-h,o-z)

      include 'common.inc'

c      double precision deth
      do 10 it=1,10
      do 40 j=1,3
         do 50 k=1,3

            vg3(j,k)=vg3(j,k)-((1.0d0/3.0d0)*iden(j,k)
     x           *(vg3(1,1)+vg3(2,2)+vg3(3,3)))

 50      continue
 40   continue
 10   continue
      return

      end
