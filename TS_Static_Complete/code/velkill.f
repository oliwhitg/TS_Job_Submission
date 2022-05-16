      subroutine velkill

c******************************************************
c Kills the particle velocities.
c******************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'

      do 10 i=1,num

         vxn(i)=0.0d0
         vyn(i)=0.0d0
         vzn(i)=0.0d0
         vx(i)=0.0d0
         vy(i)=0.0d0
         vz(i)=0.0d0

   10 continue

      write(6,*) 'Velocities killed.'

      return
      end
