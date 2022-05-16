      subroutine dynmat

c************************************************************
c************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'

      double precision xstore(nummax)
     x                ,ystore(nummax)
     x                ,zstore(nummax)
      double precision xmustore(nummax)
     x                ,ymustore(nummax)
     x                ,zmustore(nummax)
      double precision fstorex(nummax)
     x                ,fstorey(nummax)
     x                ,fstorez(nummax)
      double precision D(3,nummax,3,nummax)
cc      double precision D(3,10,3,10)

c=====================================================================
c Store the initial forces.
      do 30 i=1,num

cc         fstorex(i)=frrx(i)
cc         fstorey(i)=frry(i)
cc         fstorez(i)=frrz(i)
         xstore(i)=x(i)
         ystore(i)=y(i)
         zstore(i)=z(i)
         xmustore(i)=xmu(i)
         ymustore(i)=ymu(i)
         zmustore(i)=zmu(i)

   30 continue

c=====================================================================
c Calculate the dynamical matrix.
cc      do 40 i=nmat1,nmat2
      do 40 i=1,num

         write(6,*)i
         do 50 j=1,3

            x(i)=xstore(i)+amove(1,j)*amovefac
            y(i)=ystore(i)+amove(2,j)*amovefac
            z(i)=zstore(i)+amove(3,j)*amovefac

cc            call conjgrad
cc  RIM D calc. only
            call ener

c=====================================================================
c Store the forces for the first move.
            do 35 ii=1,num

               fstorex(ii)=frrx(ii)
               fstorey(ii)=frry(ii)
               fstorez(ii)=frrz(ii)

   35       continue

            x(i)=xstore(i)-amove(1,j)*amovefac
            y(i)=ystore(i)-amove(2,j)*amovefac
            z(i)=zstore(i)-amove(3,j)*amovefac

cc            call conjgrad
cc  RIM D calc. only
            call ener

            do 60 k=1,num

               fdiffx=frrx(k)-fstorex(k)
               fdiffy=frry(k)-fstorey(k)
               fdiffz=frrz(k)-fstorez(k)

               D(j,i,1,k)=0.5d0*fdiffx/amovefac
               D(j,i,2,k)=0.5d0*fdiffy/amovefac
               D(j,i,3,k)=0.5d0*fdiffz/amovefac

   60       continue

         x(i)=xstore(i)
         y(i)=ystore(i)
         z(i)=zstore(i)
         xmu(i)=xmustore(i)
         ymu(i)=ymustore(i)
         zmu(i)=zmustore(i)

cc         write(83,*)sreng(1,1),sreng(1,2)+sreng(2,1),sreng(2,2)
cc         write(84,*)ddeng(1,1),ddeng(1,2)+ddeng(2,1),ddeng(2,2)
cc         write(85,*)erfcracc,qdipacc+qmueng,srdipeng
cc     x             ,dipdipacc+xmumueng
cc         write(86,*)engpetot

   50    continue

   40 continue

      open(98,file='dynmat.dat',status='new',form='unformatted')
      write(98)D
      close(98)

      stop

      return

      end
