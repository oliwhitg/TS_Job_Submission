      subroutine update_velocities

c************************************************************
c   Updates:  i) ion positions.
c            ii) ion velocities.
c           iii) dipole components.
c            iv) dipole component velocities.
c
c   Note that velocities are at the half timesteps.
c************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'

      if (moveions) then
c
c Update ion positions and velocities.
c
      do 10 i=1,num

         if (fixedarray(i)) then
            goto 10
         endif

         vx(i)=vxn(i)
         vy(i)=vyn(i)
         vz(i)=vzn(i)

   10 continue

      endif
c
c Update the friction coefficients.
c

      zeta=zetanew

      if (cimlog) deltazeta=deltazetanew
      if (daimlog) then
         deltazeta=deltazetanew
         epszeta=epszetanew
      endif
      if (quaimlog) then
         deltazeta=deltazetanew
         epszeta=epszetanew
         quaimzeta=quaimzetanew
      endif

      do 33 i=1,nspec
	
         dipzeta(i)=dipzetanew(i)

 33   continue

        eps1=eps2
        eps2=eps3
        veps1=veps2
        veps2=veps3
c       vol2=vol3

      return
      end
