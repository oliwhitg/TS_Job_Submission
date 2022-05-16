      subroutine dump

c************************************************************
c Dumps ion positions, velocities and information on
c the extra degrees of freedom to disk, for later use
c as a restart file.
c************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'

      data ndumpchan/37/

      open(ndumpchan,file=fileout,status='new')

         write(ndumpchan,*) crddumplog
         write(ndumpchan,*) veldumplog
         write(ndumpchan,*) chgdumplog
         write(ndumpchan,*) chgveldulog
         write(ndumpchan,*) fulldumplog

         if (crddumplog) then
            do 10 i=1,num
               write(ndumpchan,*)x(i),y(i),z(i)
   10       continue
         endif

         if (veldumplog) then
            do 20 i=1,num
               write(ndumpchan,*)vx(i),vy(i),vz(i)
   20       continue
         endif

         if (chgdumplog) then

            do 40 i=1,num
               write(ndumpchan,*)xmu(i),ymu(i),zmu(i)
   40       continue

         endif

         if (chgveldulog) then
            do 45 i=1,num
               write(ndumpchan,*)qvelx(i),qvely(i),qvelz(i)
   45       continue
         endif

         if (fulldumplog) then
            write(ndumpchan,*)nsofar+1

            do 50 i=1,nsofar
               write(ndumpchan,*)nstpbrk(i)
   50       continue

            write(ndumpchan,*)nrun-ntotstp


c	write out barostat data

            do 60 i=1,5
               write(ndumpchan,*)vpzeta2(i)
               write(ndumpchan,*)vpzeta3(i)
               write(ndumpchan,*)pzeta3(i)
               write(ndumpchan,*)vbzeta2(i)
               write(ndumpchan,*)vbzeta3(i)
               write(ndumpchan,*)bzeta3(i)
   60       continue

            write(ndumpchan,*)eps2,eps3
            write(ndumpchan,*)veps2,veps3

            write(ndumpchan,*)vg3(1,1),vg3(2,1),vg3(3,1)
            write(ndumpchan,*)vg3(1,2),vg3(2,2),vg3(3,2)
            write(ndumpchan,*)vg3(1,3),vg3(2,3),vg3(3,3)


c  frictions etc for additional DoF's.
            if (dippimlog) then 
               do 70 i=1,nspec
                  write(ndumpchan,*)dipzeta(i),dtdip2(i)
   70          continue
               if (quadpimlog) then
                  do 71 i=1,nspec
                     write(ndumpchan,*)quadzeta(i),dtquad2(i)
   71             continue
               endif
            endif

            if (cimlog) then 
               write(ndumpchan,*)deltazeta,dtdelta2
               if (daimlog) then
                  write(ndumpchan,*)epszeta,dteps2
                  if (quaimlog) then
                     write(ndumpchan,*)quaimzeta,dtquaim2
                  endif
               endif
            endif

         endif

c	write out new cell matrix

            write(ndumpchan,*) hlab3(1,1),hlab3(1,2),hlab3(1,3)
            write(ndumpchan,*) hlab3(2,1),hlab3(2,2),hlab3(2,3)
            write(ndumpchan,*) hlab3(3,1),hlab3(3,2),hlab3(3,3)
         write (ndumpchan,*) boxlenx
         write (ndumpchan,*) boxleny
         write (ndumpchan,*) boxlenz

      close(ndumpchan)

c
c Dump out lab-frame version of the coordinates.
c Only if writing out a full restart file. (Not an annealed one.)
c
      if (fulldumplog) then
         open(50,file='finalpos.out',status='new')

         do i=1,num
            write(50,*) h(1,1)*x(i)+h(1,2)*y(i)+h(1,3)*z(i),
     x                  h(2,1)*x(i)+h(2,2)*y(i)+h(2,3)*z(i),
     x                  h(3,1)*x(i)+h(3,2)*y(i)+h(3,3)*z(i)
         enddo

         close(50)
      endif

      if (fulldumplog) then
         write(6,*)'**** Full restart file written out ****'
      else
         write(6,*)'**** Annealed restart file written out ****'
      endif

      return

      end
