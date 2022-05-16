      subroutine ener
c
c Calls correct energy routine.
c
      implicit double precision (a-h,o-z)
      include 'common.inc'

c dipolar polarizable ion only.
      if (dippimlog) then
         call recipE_dippim
         call realE_dippim
         call sr_energy
      endif

      if (rimlog) then
         call rgdrecipE
         call rgdrealE
         call sr_energy
      endif

      if (nochargelog) then
         call sr_energy
      endif

      return

      end
