      subroutine ener
c
c Calls correct energy routine.
c
      implicit double precision (a-h,o-z)
      include 'common.inc'
      include 'lj.inc'

c quadrupolar polarizable ion only.
      if (quadpimlog) then
         call recipE_quadpim
         call realE_quadpim
         call sr_energy
      endif

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

cc added Sept. 2009. 
         engpetot=0.0d0
         do 290 i=1,num
            frrx(i)=0.0d0
            frry(i)=0.0d0
            frrz(i)=0.0d0
  290    continue

         call sr_energy
      endif

c      if (ljsubsystem) then
c         call lj_ions
c      endif

      return

      end
