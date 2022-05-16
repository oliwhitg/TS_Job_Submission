      subroutine conjgrad

c************************************************************
c  Performs an energy minimilisation on the startup configuration
c via a conjugant gradient method.
c
c  This subroutine can be slotted into existing code and called
c where conjrod would have been. The only variable to add to
c the code is 'conjtol' in readin to control the required
c energy convergence tolerance.
c************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'

c================================================
      double precision p(nummax*3)
      double precision xi(nummax*3)

      logical updatelog

      updatelog=.true.

      rotkel=0.0d0

      do 10 i=1,nummax

         p(i)=xmu(i)
         p(i+nummax)=ymu(i)
         p(i+nummax+nummax)=zmu(i)

   10 continue

cc      ftol=1.0d-10
      call frprmn(p,nummax*3,iter,fret)

      do 15 i=1,nummax

         xmu(i)=p(i)
         ymu(i)=p(i+nummax)
         zmu(i)=p(i+nummax+nummax)

   15 continue

      CALL DFUNC(P,XI)

      call ener

      do 111 i=1,10
         write(17,*)xmu(i),ymu(i),zmu(i)
  111 continue

      write(6,*)
      write(6,*)'**** Anion annealing completed ****'
      write(6,*)

      return

      end
