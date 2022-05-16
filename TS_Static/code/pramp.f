      subroutine pramp

c************************************************************
c  pressure ramping.             
c************************************************************

      implicit double precision (a-h,o-z)

      include 'common.inc'

      pext=pext+deltapramp

      write(6,*)'pext now ',pext

      return

      end
