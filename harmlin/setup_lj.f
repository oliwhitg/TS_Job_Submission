      subroutine setup_lj

      implicit double precision(a-h,o-z)

      include 'common.inc'
      include 'lj.inc'

c=====================================================================
c Atomic unit of mass (electron mass) in grams.
      emass=9.109534d-028

      xk1lj=1.0d0/(2.0d0*polarlj)

c=====================================================================
c  Thermostat parameters.
      relaxsqlj=relaxlj*relaxlj
      relaxsqreclj=1.0d0/relaxsqlj

      chrgrlxsqlj=chrgrelaxlj*chrgrelaxlj
      chrgrlxsqreclj=1.0d0/chrgrlxsqlj

      chrgtkblj=boltz*chrgtemplj

c=====================================================================
c  Convert masses from .inpt (gmol^-1) into internal (atomic) units.
      amasslj=amasslj/(avo*emass)
      recamasslj=1.0d0/amasslj
      hmasslj=amasslj/2.0d0

      dmsdlj=0.0d0

      qmasschgreclj=1.0d0/(qmasschglj)

      gtranlj=3.0d0*float(numlj-1)*trantkb
      gchglj=float(numlj)*chrgtkblj
      gtranreclj=1.0d0/gtranlj
      gchgreclj=1.0d0/gchglj
      gtrantimlj=gtranlj*dtime/trantemplj
      gchgtimlj=gchglj*dtime/chrgtemplj
      gtrankinlj=2.0d0/(3.0d0*float(numlj-1)*boltz)
      gchgkinlj=1.0d0/(float(numlj)*boltz)

      do 333 i=1,3

         dddamp2(i,3)=dddamp(i,3)*dddamp(i,3)/2.0d0
         dddamp3(i,3)=dddamp2(i,3)*dddamp(i,3)/3.0d0
         dddamp4(i,3)=dddamp3(i,3)*dddamp(i,3)/4.0d0
         dddamp5(i,3)=dddamp4(i,3)*dddamp(i,3)/5.0d0
         dddamp6(i,3)=dddamp5(i,3)*dddamp(i,3)/6.0d0
         dddamp6(3,i)=dddamp6(i,3)
         dddamp5(3,i)=dddamp5(i,3)
         dddamp4(3,i)=dddamp4(i,3)
         dddamp3(3,i)=dddamp3(i,3)
         dddamp2(3,i)=dddamp2(i,3)

         dqdamp2(i,3)=dqdamp(i,3)*dqdamp(i,3)/2.0d0
         dqdamp3(i,3)=dqdamp2(i,3)*dqdamp(i,3)/3.0d0
         dqdamp4(i,3)=dqdamp3(i,3)*dqdamp(i,3)/4.0d0
         dqdamp5(i,3)=dqdamp4(i,3)*dqdamp(i,3)/5.0d0
         dqdamp6(i,3)=dqdamp5(i,3)*dqdamp(i,3)/6.0d0
         dqdamp7(i,3)=dqdamp6(i,3)*dqdamp(i,3)/7.0d0
         dqdamp8(i,3)=dqdamp7(i,3)*dqdamp(i,3)/8.0d0
         dqdamp8(3,i)=dqdamp8(i,3)
         dqdamp7(3,i)=dqdamp7(i,3)
         dqdamp6(3,i)=dqdamp6(i,3)
         dqdamp5(3,i)=dqdamp5(i,3)
         dqdamp4(3,i)=dqdamp4(i,3)
         dqdamp3(3,i)=dqdamp3(i,3)
         dqdamp2(3,i)=dqdamp2(i,3)

  333 continue

      write(6,*)
      write(6,*)'**** Internal LJ variables set up ****'
      write(6,*)

      return

      end
