      subroutine shortrange

c************************************************************
c
c  Adjusts the short range tensor and inverts the matrix via
c 
c
c************************************************************

      implicit double precision(a-h,o-z)
      
      include 'common.inc'
 

c=====================================================================
c Local double precision arrays.
      double precision inptarr(3,3),outarr(3,3)

      do 111 i=1,3
         do 111 j=1,3
            inptarr(i,j)=0.0d0
            outarr(i,j)=0.0d0
  111 continue

c      write(999,*)'xx',inptarr,outarr


c=====================================================================
c Set up NAG variable.
      ifail=0

      do 10 i=1,nanion

c=====================================================================
c Add the undistorted polarizability to the trace terms.
         inptarr(1,1)=srxx(i) + polarundist
         inptarr(1,2)=srxy(i)
         inptarr(1,3)=srxz(i)
         inptarr(2,2)=sryy(i) + polarundist
         inptarr(2,3)=sryz(i)
         inptarr(3,3)=srzz(i) + polarundist
         inptarr(2,1)=inptarr(1,2)
         inptarr(3,1)=inptarr(1,3)
         inptarr(3,2)=inptarr(2,3)

cc      write(98,*)inptarr(1,1)
cc     x          ,inptarr(2,2),inptarr(3,3)
c=====================================================================
c Invert the matrix inptarr putting the result in outarr.
c old NAG call left in for comparison.
cc         call f01abf(inptarr,4,3,invarr,3,work,ifail)
        call invert(inptarr,outarr,dum)

cc      write(98,*)outarr(1,1)
cc     x          ,outarr(2,2),outarr(3,3)
c=====================================================================
c Subtract the average polarizability (polarav) from the trace terms.
         srxx(i)=outarr(1,1)-xlipol(1)
         srxy(i)=outarr(1,2)
         srxz(i)=outarr(1,3)
         sryy(i)=outarr(2,2)-xlipol(1)
         sryz(i)=outarr(2,3)
         srzz(i)=outarr(3,3)-xlipol(1)

   10 continue

      return

      end
