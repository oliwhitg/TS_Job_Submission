      subroutine dofmoveselect
c
c Calls correct energy routine.
c
      implicit double precision (a-h,o-z)
      include 'common.inc'

      if (dippimlog) call dipolemove

      if (cimlog) then
         if (cim1log) call anionmv
         if (cim2log) call anionmv_cim2
      endif

      if (daimlog) then
         call anionmv_cim2
         call daimmv
      endif

      if (quaimlog) then
         call anionmv_cim2
         call daimmv
         call quaimmv
      endif

      if (quadpimlog) call dipquadmove

      return

      end
