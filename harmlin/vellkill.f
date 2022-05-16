      subroutine vellkill

c******************************************************
c Kills the particle velocities.
c******************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'

      n=1

      do 10 i=1,nspec

         do 20 j=1,nsp(i)

            vx(n)=0.0d0
            vy(n)=0.0d0
            vz(n)=0.0d0

            n=n+1

   20    continue

   10 continue

      write(6,*) 'Velocities killed.'

      return
      end
