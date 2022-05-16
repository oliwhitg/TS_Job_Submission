      subroutine lambdaramp

c************************************************************
c  pressure ramping.             
c************************************************************

      implicit double precision (a-h,o-z)

      include 'common.inc'
      include 'lj.inc'

      xlambdalj=xlambdalj+deltalambdalj

      write(6,*)'lambdalj now ',xlambdalj

      return

      end
