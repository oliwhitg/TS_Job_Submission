      program main
c
c*****************************************************************
c M-D code to model ionic systems with anion dipoles:
c Everything is in atomic units.
c*****************************************************************

      implicit double precision (a-h,o-z)

      include 'common.inc'
      include 'lj.inc'

      endrun=.false.
      shutdown=.false.
      nstep=1
c
c Read in the run-time parameters.
c
       write(*,*) "in main file"
      call readin
c      call readin_lj

c
c Set up the auxillary variables and I/O.
c
      call setup
c      call setup_lj

      if (restart) then
         call rstrun
c         call rstrun_lj
      endif

      if (dynam) then
         call dynmat
         stop
      endif

      if (velinit) call velset

c new option. 1.4.2000 to do CG or SD minimisation.
      if (conjrodlog) then
         dipcalclog=.true.
         call timeset(dtanneal)
         write (6,*) 'conjrod called.'
         call conjrod
         call timeset(dmaintime)
      endif

      if (conjgradlog) call conjgrad

      if (rescalelog) then
         if (.not. velinit .and. .not. restart) then
            write (6,*) 'Cannot rescale from zero temperature!'
            stop
         endif
         write(6,*) 'rescale called.' 
         call rescale
      endif

      if (conjrodlog) then

         crddumplog=.true.
         veldumplog=.true.
         chgdumplog=.true.
         chgveldulog=.false.
         fulldumplog=.false.
         fileout=crstout
         call dump

      endif

cc      if (irreadlog.and.docfcalc) then
cc         write(6,*) 'correlation functions being restored'
cc         call ircfrstrestore
cc      endif

      write(6,*) 'header file being written out'
      call header

cc      if (relaxconfig) then
cc         call relax_struct
cc         stop
cc      endif

      if (nrun.eq.0) then
         stop
      endif
c
c Main loop.

      forfl=.false.

      write(6,*) 'main loop begun'

      call ener

      do while (nstep.le.nrun)
 
         write (6,*) nstep

         if (mod(float(nstep),float(nrscale)).eq.0) then
            call rescale
            write(6,*)'Rescale called'
         endif

         call trans_vv

c ****** Periodic parts of main loop ******
c
c Periodic pressure calculation:
c
cc      if (mod(float(nstep),float(nboxmove)).eq.0) then
cc         write(6,*)'called pimpressure'
cc         call pimpressure
cc      endif
c
c Periodic energies output.
c
         if (mod(float(nstep),float(npereng)).eq.0) then
            pereng=.true.
            call output
c            call output_lj
         endif
c
c Periodic velocities output.
c
         if (mod(float(nstep),float(npervel)).eq.0) then
            pervel=.true.
            call output
         endif
c
c Periodic temperature output.
c
         if (mod(float(nstep),float(nperfri)).eq.0) then
            perfric=.true.
            call output
         endif

         if (mod(float(nstep),float(npercell)).eq.0) then
            percell=.true.
            call output
         endif
c=====================================================================
c Call light scattering routines 
 
c Periodic correlation function calculation.
c Call both the normal cf and the light scattering cf routines.
cc      if (docfcalc) then
cc         if (mod(float(nstep),float(npercf)).eq.0) then

cc            call lightfields
cc            call shortrange
cc            call lightarray
cc            call lightcfcalc
c
c update correlation fn array pointers
c
cc            ncorrtime=mod(ncorrtime,ncorr)
cc            ncorrcall=ncorrcall+1
cc            ncorrtime=ncorrtime+1
cc         endif
cc      endif

c===============================================================
c periodic T/p ramping if required.
         if ((tramplog).and.
     x       (mod(float(nstep),float(nsteptramp)).eq.0)) then
            call tramp
         endif
         if ((pramplog).and.
     x       (mod(float(nstep),float(nsteppramp)).eq.0)) then
            call pramp
         endif
c        if ((fadeljpotlog).and.
c     x       (mod(float(nstep),float(nstepfadelj)).eq.0)) then
c            call lambdaramp
c         endif

c===============================================================
c Periodic radial distribution function calculation.
c
         if (mod(float(nstep),float(nperrdf)).eq.0) then
            call rdf
c            call rdf_lj
         endif
c
c Periodic structure factor calculation.
c
         if (mod(float(nstep),float(npersk)).eq.0) then
cc            call struc
            call debye_scherer
            nskcall=nskcall+1
         endif
c
c Periodic annealing.
c
         if (mod(float(nstep),float(nperanneal)).eq.0) then
            call timeset(dtanneal)
            dipcalclog=.false.
            write (6,*) 'conjrod called.'
            call conjrod
            call timeset(dmaintime)
         endif
c
c Periodic pressure calculation.
c
         if (mod(float(nstep),float(nperpres)).eq.0) then
            call prescalc
         endif
c
c Update ion variables.
c
         call update_velocities
c
c Periodic msd output:
c
         if (msdcalllog) then
            if (mod(float(nstep),float(nmsdcalltime)).eq.0) then
               call dispout
            endif
         endif

         if (nstep.eq.nrun) endrun=.true.

         if (endrun) then
c
c Implement shutdown prodedure:
c     i) Dump final positions, velocities etc.
c    ii) Dump restart (.rst) file.
c   iii) Dump correlation function data.
c
            shutdown=.true.

            call output
            veldumplog=.true.
            crddumplog=.true.
            chgdumplog=.true.
            chgveldulog=.true.
            fulldumplog=.true.
            fileout=rstout
	    call dump
 
c=====================================================================
c  Dump the light scattering correlation functions and the store arrays required
c to restart the calculation.
cc            if (docfcalc) then
cc               call lightcfdump
cc            endif

            call rdfout
cc            call strucout
            call debye_scherer_out
cc            call strucnewout

            close(21)
            close(22)
            close(23)
            close(24)
            close(25)
            close(26)
            close(27)
            close(28)
            close(29)
            close(30)
            close(48)            
            close(58)
            close(59)

            do 963 i=1,nummon

               nnn=3*(i-1)

               close(77+nnn)
               close(78+nnn)
               close(79+nnn)

  963       continue

         endif

         nstep=nstep+1

      enddo

      stop

      end
