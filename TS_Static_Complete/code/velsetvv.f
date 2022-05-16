      subroutine velset

c************************************************************
c   Sets up the ion velocities on a Gaussian distribution.
c************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'

      n=1

      do 10 i=1,nspec
         do 20 j=1,nsp(i)

            vx(n)=gauss(dummy)*vartrans(i)
            vy(n)=gauss(dummy)*vartrans(i)
            vz(n)=gauss(dummy)*vartrans(i)

            n=n+1

   20    continue
   10 continue
c
c Rescale the velocities.
c
      call rescale
c
c Calculate the initial value of the friction
c coefficient for the Nose-Hoover thermostat.
c Note that it should be very small due to the 
c temperature rescaling.
c
c      zeta2=(tkin-trantemp)/(2.0d0*relax*tkin)
c      zeta=(tkinres-trantemp)/(2.0d0*relax*tkinres)

        zeta2=0.0d0
        zeta=0.0d0

      dtzeta2=dtime2*zeta

      write(6,*)
      write(6,*)'**** Velocities set up ****'
      write(6,*)

      return

      end
