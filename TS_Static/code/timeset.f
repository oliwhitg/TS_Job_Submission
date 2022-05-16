      subroutine timeset(time)

c****************************************************************
c Resets time dependent variables to do the annealing correctly.
c There is no need to change thermostat dependent variables,
c as no thermostating is done during the annealing.
c****************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'

      dtime=time
      dtimerec=1.0d0/dtime
      dtime2=dtime*0.5d0

      write (6,*) 'timestep reset to', real(dtime), 'a.u.'

      return
      end
