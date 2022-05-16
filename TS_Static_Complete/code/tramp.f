      subroutine tramp

c************************************************************
c  temperature ramping.             
c************************************************************

      implicit double precision (a-h,o-z)

      include 'common.inc'

      trantemp=trantemp+deltatramp

      trantkb=boltz*trantemp
      gtran=3.0d0*dble(num-1)*trantkb
      gtranrec=1.0d0/gtran
      gtrantim=gtran*dtime/trantemp
      gtrankin=2.0d0/(3.0d0*dble(num-1)*boltz)

      free=3.0d0*float(num-1)
      W=(free+3.0d0)*trantkb*(relaxb**2.0d0)
      Wgo=(free+3.0d0)*trantkb*(relaxb2**2.0d0)/3.0d0
      Wrec=1.0d0/W
      Wgorec=1.0d0/Wgo

      CUEb=dom*trantkb/relaxsqrec
      CUEb2=CUEb/dom
      CUEbrec=1.0d0/CUEb
      CUEb2rec=1.0d0/CUEb2

      if (.not.nth) then
         CUEbrec=0.0d0
         CUEb2rec=0.0d0
      endif
      if (.not.nib) then
         CUEb=0.0d0
         CUEb2=0.0d0
         CUEbrec=0.0d0
         CUEb2rec=0.0d0
      endif

      if (.not.nib) then
         W=0.0d0
         Wrec=0.0d0
      endif
      if (.not.nab) then
         Wgo=0.0d0
         Wgorec=0.0d0
      endif

      if (tramprescalelog) call rescalen

      write(6,*)'Temperature now ',trantemp

      return

      end
