      subroutine rstrun

c************************************************************
c Reads in ion positions,velocities and dipolar values and
c 'velocities' and thermostating information
c from a file (filename in the filenames.inpt file) for restart.
c************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'

      logical crdlog,vellog,chglog
     x       ,chgvellog,fullrunlog

      data nrstchan/35/

      write (6,*) 'No. of ions =',num

      open(nrstchan,file=rstin,status='old')

         read(nrstchan,*)crdlog
         read(nrstchan,*)vellog
         read(nrstchan,*)chglog
         read(nrstchan,*)chgvellog
         read(nrstchan,*)fullrunlog

c Read in positions.
         if (crdlog) then
            do 10 i=1,num
               read(nrstchan,*)x(i),y(i),z(i)
   10       continue
         else
            write(6,*)'No coordinate information in restart file!'
            stop
         endif

c Read in velocities.
c
c we need to get the velocities from somewhere!
c
         if (vellog) then
            do 20 i=1,num
               read(nrstchan,*)vx(i),vy(i),vz(i)
   20       continue
         else
c
c We must generate some velocities if we are going to rescale,
c and currently, the code assumes that rescale is OK if velset
c has been called or the run has been restarted. Having no
c velocities is OK otherwise - of course it may be desired!
c
            if (rescalelog) call velset
         endif

c Read in dipoles.
         if (chglog) then
            do 40 i=1,num
               read(nrstchan,*)xmu(i),ymu(i),zmu(i)
   40       continue
            if (quadpimlog) then
               do 41 i=1,num
                  read(nrstchan,*)quadxx(i),quadyy(i),quadzz(i)
                  read(nrstchan,*)quadxy(i),quadxz(i),quadyz(i)
   41          continue
            endif

            if (cimlog) then
               do 42 i=1,num
                  read(nrstchan,*)delta(i)
   42          continue

               if (daimlog) then
                  do 43 i=1,num
                     read(nrstchan,*)epsilonx(i)
     x                             ,epsilony(i),epsilonz(i)
   43             continue

                  if (quaimlog) then
                     do 44 i=1,num
                        read(nrstchan,*)quaimxx(i)
     x                             ,quaimyy(i),quaimzz(i)
                        read(nrstchan,*)quaimxy(i)
     x                             ,quaimxz(i),quaimyz(i)
   44                continue
                  endif
               endif
            endif
         else
            if (.not. conjrodlog .and. .not. rim) then
               write(6,*)'Warning: No dipole variables read in,'
               write(6,*)'         and no annealing is scheduled.'
            endif
         endif

c read in dipole velocities.
         if (chgvellog) then
            do 45 i=1,num
               read(nrstchan,*)qvelx(i),qvely(i),qvelz(i)
   45       continue
            if (quadpimlog) then
               do 46 i=1,num
                  read(nrstchan,*)quadvelxx(i)
     x                ,quadvelyy(i),quadvelzz(i)
                  read(nrstchan,*)quadvelxy(i)
     x                ,quadvelxz(i),quadvelyz(i)
   46          continue
            endif

            if (cimlog) then
               do 47 i=1,num
                  read(nrstchan,*)deltavel(i)
   47          continue

               if (daimlog) then
                  do 48 i=1,num
                     read(nrstchan,*)epsvelx(i)
     x                ,epsvely(i),epsvelz(i)
   48             continue

                  if (quaimlog) then
                     do 49 i=1,num
                        read(nrstchan,*)quaimvelxx(i)
     x                          ,quaimvelyy(i)
     x                          ,quaimvelzz(i)
                        read(nrstchan,*)quaimvelxy(i)
     x                          ,quaimvelxz(i)
     x                          ,quaimvelyz(i)
   49                continue
                  endif
               endif
            endif
         endif

c read in all other information. TMP
         if (fullrunlog) then
            ntotstp=0
            read(nrstchan,*)nsofar
            do 50 i=1,nsofar

               read(nrstchan,*)nstpbrk(i)
               ntotstp=ntotstp+nstpbrk(i)

   50       continue
c
c Get correct step counting for this run.
c
            nstep=nstep+ntotstp
            nrun=nrun+ntotstp
c
c read in thermostat information:
c
c	readin barostat data

            do 60 i=1,5
               read(nrstchan,*)vpzeta2(i)
               read(nrstchan,*)vpzeta3(i)
               read(nrstchan,*)pzeta3(i)
               read(nrstchan,*)vbzeta2(i)
               read(nrstchan,*)vbzeta3(i)
               read(nrstchan,*)bzeta3(i)
   60       continue

            read(nrstchan,*)eps2,eps3
            read(nrstchan,*)veps2,veps3

            read(nrstchan,*)vg3(1,1),vg3(2,1),vg3(3,1)
            read(nrstchan,*)vg3(1,2),vg3(2,2),vg3(3,2)
            read(nrstchan,*)vg3(1,3),vg3(2,3),vg3(3,3)

c	readin new cell matrix

c  frictions etc for additional DoF's.
            if (dippimlog) then
               do 70 i=1,nspec
                  read(nrstchan,*)dipzeta(i),dtdip2(i)
   70          continue
               if (quadpimlog) then
                  do 71 i=1,nspec
                     read(nrstchan,*)quadzeta(i),dtquad2(i)
   71             continue
               endif
            endif

            if (cimlog) then
               read(nrstchan,*)deltazeta,dtdelta2
               if (daimlog) then
                  read(nrstchan,*)epszeta,dteps2
                  if (quaimlog) then
                     read(nrstchan,*)quaimzeta,dtquaim2
                  endif
               endif
            endif

         endif

	 if (boxlenfromrstlog) then

            read(nrstchan,*) h(1,1),h(1,2),h(1,3)
            read(nrstchan,*) h(2,1),h(2,2),h(2,3)
            read(nrstchan,*) h(3,1),h(3,2),h(3,3)

            read(nrstchan,*)boxlenx
            read(nrstchan,*)boxleny
            read(nrstchan,*)boxlenz

	 endif

c	obtain vol by calling dcell

	 call dcell(h,b)

         vol3=boxlenx*boxleny*boxlenz*b(10)
         do 100 j=1,3
            do 110 k=1,3
               h3(j,k)=0.0d0
               hlab3(j,k)=h(j,k)
               hlab2(j,k)=h(j,k) 
 110        continue
 100     continue

         h3(1,1)=h(1,1)*boxlenx/(vol3**(1.0d0/3.0d0))
         h3(2,1)=h(2,1)*boxlenx/(vol3**(1.0d0/3.0d0))
         h3(3,1)=h(3,1)*boxlenx/(vol3**(1.0d0/3.0d0))
         h3(1,2)=h(1,2)*boxleny/(vol3**(1.0d0/3.0d0))
         h3(2,2)=h(2,2)*boxleny/(vol3**(1.0d0/3.0d0))
         h3(3,2)=h(3,2)*boxleny/(vol3**(1.0d0/3.0d0))
         h3(1,3)=h(1,3)*boxlenz/(vol3**(1.0d0/3.0d0))
         h3(2,3)=h(2,3)*boxlenz/(vol3**(1.0d0/3.0d0))
         h3(3,3)=h(3,3)*boxlenz/(vol3**(1.0d0/3.0d0))

         call matinv(h3,hi3,detmh)

         if (boxlenfromrstlog) then
            call boxreset
            call dcell(fullh,b)
            halfboxminsq=(0.5d0*dmin1(b(7),b(8),b(9)))**2
            do i=1,3
               do j=1,3
                  fullhit(i,j)=fullhi(j,i)
               enddo
            enddo

            call kset

            rksqmax_ds=float(kmax_ds)*float(kmax_ds)*pi*pi/halfboxminsq
         endif

      close(nrstchan)

      if (restart) then

         open(49,file='startuppos.out',status='new')

         do i=1,num
            write(49,*) h(1,1)*x(i)+h(1,2)*y(i)+h(1,3)*z(i),
     x                  h(2,1)*x(i)+h(2,2)*y(i)+h(2,3)*z(i),
     x                  h(3,1)*x(i)+h(3,2)*y(i)+h(3,3)*z(i)
         enddo

         close(49)

      endif

      if (displace) then

         write (6,*) 'Randomly displacing ions off lattice.'

         do 72 i=1,num

            x(i)=x(i)+0.1d0*(ran1(idum)-0.5d0)
            y(i)=y(i)+0.1d0*(ran1(idum)-0.5d0)
            z(i)=z(i)+0.1d0*(ran1(idum)-0.5d0)

   72    continue

      endif


      write(6,*)
      write(6,*)'**** Run restarted ****'
      write(6,*)

      return

      end
