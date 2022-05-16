      subroutine dofmoveselect
c
c Calls correct energy routine.
c
      implicit double precision (a-h,o-z)
      include 'common.inc'

      if (dippimlog) call dipolemove

      return

      end
