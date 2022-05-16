      subroutine rstrun_lj

c************************************************************
c************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'
      include 'lj.inc'

      logical crdlog,vellog,chgannlog,rotannlog
      character*80 rstinlj

      nrstchan=20
      open(nrstchan,file=rstinlj,status='old')

         read(nrstchan,*)crdlog
         read(nrstchan,*)vellog
         read(nrstchan,*)chgannlog
         read(nrstchan,*)rotannlog

         if (crdlog) then
            do 10 i=1,numlj
               read(nrstchan,*)xlj(i),ylj(i),zlj(i)
   10       continue
         endif

         if (vellog) then
            do 20 i=1,numlj
               read(nrstchan,*)vxlj(i),vylj(i),vzlj(i)
   20       continue
         endif

         if (chgannlog) then
            do 40 i=1,numlj
               read(nrstchan,*)xmulj(i),ymulj(i),zmulj(i)
   40       continue
         endif

         if (rotannlog) then
            do 45 i=1,numlj
               read(nrstchan,*)qvelxlj(i),qvelylj(i),qvelzlj(i)
   45       continue
         endif

      close(nrstchan)

      write(6,*)
      write(6,*)'**** LJ Run restarted ****'
      write(6,*)

      return

      end
