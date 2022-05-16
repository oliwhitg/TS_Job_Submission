      subroutine readin

c************************************************************
c   Reads in the data from the .inpt file.
c************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'
c
c Local characters.
      character*80 inptfile
      character*80 FTfile
      character*80 cffile
      character*80 lcffile
      character*80 cellfile

      character*10 runtype
      character*10 runtype2
      character*10 pottype
      character*10 rstchar
c
c Input channel nos.
c
      data nfileinpt,NFTinpt,NCinpt,ncfinpt/13,10,12,14/

c
c Input I/O filenames.
c
c 29/3/2000 sorted into inputs then outputs
c
      open(nfileinpt,file='filenames.inpt',status='old')

         read(nfileinpt,'(a)')inptfile
         read(nfileinpt,'(a)')FTfile
         read(nfileinpt,'(a)')cffile
         read(nfileinpt,'(a)')lcffile
         read(nfileinpt,'(a)')cellfile
         read(nfileinpt,'(a)')rstin
         read(nfileinpt,'(a)')cfin
         read(nfileinpt,'(a)')engout
         read(nfileinpt,'(a)')engout2
         read(nfileinpt,'(a)')cellengout
         read(nfileinpt,'(a)')engtotout
         read(nfileinpt,'(a)')pzetaout
         read(nfileinpt,'(a)')bzetaout
         read(nfileinpt,'(a)')fkeout
         read(nfileinpt,'(a)')polengout
         read(nfileinpt,'(a)')tdipout
         read(nfileinpt,'(a)')diagpresout
         read(nfileinpt,'(a)')presout
         read(nfileinpt,'(a)')cellanglesout
         read(nfileinpt,'(a)')celllenout
         read(nfileinpt,'(a)')cellvolout
         read(nfileinpt,'(a)')diagstressout
         read(nfileinpt,'(a)')xxyyzzstressout
         read(nfileinpt,'(a)')xyxzyzstressout
         read(nfileinpt,'(a)')polstressout
         read(nfileinpt,'(a)')coulsrstressout
         read(nfileinpt,'(a)')poscartout
         read(nfileinpt,'(a)')cellboxout
         read(nfileinpt,'(a)')crdout
         read(nfileinpt,'(a)')velout
         read(nfileinpt,'(a)')accnout
         read(nfileinpt,'(a)')rstout
         read(nfileinpt,'(a)')crstout
         read(nfileinpt,'(a)')head
         read(nfileinpt,'(a)')qout
         read(nfileinpt,'(a)')zetout
         read(nfileinpt,'(a)')fluxout
         read(nfileinpt,'(a)')tempout
         read(nfileinpt,'(a)')eleout
         read(nfileinpt,'(a)')cfout
         read(nfileinpt,'(a)')momout
         read(nfileinpt,'(a)')fgout

      close(nfileinpt)
c
c Input run-time parameters.
c Modifed 1.4.2000 to gather together parameters related
c to similar functions.
      open(NCinpt,file=inptfile,status='old')

c run-time
         read(NCinpt,*)nrun
         read(NCinpt,*)trantemp
c....use rmm as a check for amass and use to calculate box lengths if required.
         read(NCinpt,*)rmm
         read(NCinpt,*)nionunit
         read(NCinpt,*)nspec

         do 10 i=1,nspec
            read(NCinpt,*)nsp(i)
            read(NCinpt,*)chg(i)
            read(NCinpt,*)amass(i)
            read(NCinpt,*)polarizablelog(i)
   10    continue

         read(NCinpt,*)dmaintime

c type of run

         read(NCinpt,'(a)')runtype
         if (runtype(1:3).eq.'rim') rimlog=.true.
         if (runtype(1:6).eq.'dippim') dippimlog=.true.
         if (runtype(1:7).eq.'quadpim') quadpimlog=.true.
         if (runtype(1:8).eq.'nocharge') nochargelog=.true.
         read(NCinpt,'(a)')runtype2
         if (runtype2(1:3).eq.'epp') epplog=.true.
         if (runtype2(1:3).eq.'rvp') then
            epplog=.true.
            rvplog=.true.
          endif
         if (runtype2(1:2).eq.'ts') then
            epplog=.true.
            tslog=.true.
         endif
         if (runtype2(1:3).eq.'cim') cimlog=.true.
         if (runtype2(1:4).eq.'daim') daimlog=.true.
         if (runtype2(1:5).eq.'quaim') quaimlog=.true.

c anneal options
         read (NCinpt,*) conjrodlog
         if (conjrodlog) then
            read(NCinpt,*)dtanneal
            read(NCinpt,*)nrodstep
            read(NCinpt,*)annealrad
         endif
         read (NCinpt,*) conjgradlog
         if (conjgradlog) then
c TOLERANCES HERE....
            read(NCinpt,*)tol
            read(NCinpt,*)ftol
cc            read(NCinpt,*)dummy
         endif

c perform restart? if so , where do we get the boxlengths from?
         read (NCinpt,*) restart
         if (restart) then
            read(NCinpt,'(a)')rstchar
            if (rstchar(1:4).eq.'dens') then
               boxlenfromdenslog=.true.
               boxlenfrominptlog=.false.
               boxlenfromrstlog=.false.
               read(NCinpt,*)dens
            else if (rstchar(1:4).eq.'inpt') then
               boxlenfromdenslog=.false.
               boxlenfrominptlog=.true.
               boxlenfromrstlog=.false.
            else if (rstchar(1:3).eq.'rst') then
               boxlenfromdenslog=.false.
               boxlenfrominptlog=.false.
               boxlenfromrstlog=.true.
            else
               write(6,*)'No info. on restart cell origin.'
               stop
            endif
         endif

c*** leave these 6 for the moment
         read (NCinpt,*) velinit
         read (NCinpt,*) rescalelog
         read (NCinpt,*) displace
         read (NCinpt,*) moveions
         read (NCinpt,*) dynam
         read (NCinpt,*) relaxconfig

c periodic output
         read(NCinpt,*)npereng
         read(NCinpt,*)npervel
         read(NCinpt,*)nperfri
	 read(NCinpt,*)npercell

c permanent periodic calls.
         read(NCinpt,*)nperrdf
         read(NCinpt,*)npersk

c DoF monitoring.
         read(NCinpt,*)nummon
         do 150 i=1,nummon
            read(NCinpt,*)nmon(i)
  150    continue

c energy calculation.
cc         read(NCinpt,*)kmax
         read(NCinpt,*)etainpt
         read(NCinpt,*)rcut
         read(NCinpt,*)conv

c thermostat parameters.
         read(NCinpt,*)nth
         if (nth) then
            read(NCinpt,*)relax
         endif

         read(NCinpt,*)nrscalelog
         if (nrscalelog) then
            read(NCinpt,*)nrscale
         else
            nrscale=10000000
         endif

         read(NCinpt,*)nanneallog
         if (nanneallog) then
            read(NCinpt,*)nperanneal
         else
            nperanneal=10000000
         endif

c barostat parameters
         read(NCinpt,*)nib
         if (nib) then
            read(NCinpt,*)taub
            read(NCinpt,*)pext
         endif
         read(NCinpt,*)nab
         if (nab) read(NCinpt,*)taub2
         read(NCinpt,*)ortho

c T/p ramping logical
         read(NCinpt,*)tramplog
         if (tramplog) then
            read(NCinpt,*)nsteptramp
            read(NCinpt,*)deltatramp
            read(NCinpt,*)tramprescalelog
         endif
         read(NCinpt,*)pramplog
         if (pramplog) then
            read(NCinpt,*)nsteppramp
            read(NCinpt,*)deltapramp
         endif
         read(NCinpt,*)fglog

         read(NCinpt,*)fixedlog
         if (fixedlog) then
         read(NCinpt,*)nionsfixed
            open(50,file='fixed.inpt',status='old')
               write(6,*)'fixed'
               do 654 i=1,nionsfixed
                  read(50,*)nfixed(i)
  654          continue
            close(50)
         endif

c********
c DIPPIM READ.
cc         read (NCinpt,*) assign






      close(NCinpt)
c
c Input Fumi-Tosi parameters (in atomic units).
c
      open(NFTinpt,file=FTfile,status='old')

         read (NFTinpt,'(a)')pottype

         if (pottype(1:2).eq.'FT') then

            do 20 i=1,nspec

               do 30 j=i,nspec

                  read(NFTinpt,*)ftalp(i,j)
                  read(NFTinpt,*)ftb(i,j)
                  read(NFTinpt,*)ftc(i,j)
                  read(NFTinpt,*)ftd(i,j)
                  read(NFTinpt,*)dddamp(i,j)
                  read(NFTinpt,*)dqdamp(i,j)

                  ftalp(j,i)=ftalp(i,j)
                  ftb(j,i)=ftb(i,j)
                  ftc(j,i)=ftc(i,j)
                  ftd(j,i)=ftd(i,j)
                  dddamp(j,i)=dddamp(i,j)
                  dqdamp(j,i)=dqdamp(i,j)

   30          continue

   20       continue

         endif

c RVP type.
         rvplog=.false.    
         if (pottype(1:3).eq.'RVP') then
            rvplog=.true.    

            read(NFTinpt,*)rvpmax
            do 201 i=1,nspec

               do 301 j=i,nspec

                  read(NFTinpt,*)rvph(i,j)
                  read(NFTinpt,*)rvpn(i,j)
                  read(NFTinpt,*)rvpr4(i,j)
                  read(NFTinpt,*)ftc(i,j)

                  ftc(j,i)=ftc(i,j)
                  rvph(j,i)=rvph(i,j)
                  rvpn(j,i)=rvpn(i,j)
                  rvpr4(j,i)=rvpr4(i,j)

  301          continue

  201       continue

         endif
c Gaussian core type.
         gcpotlog=.false.
         if (pottype(1:2).eq.'GC') then
            gcpotlog=.true.

            do 531 i=1,nspec

               do 431 j=i,nspec

                  read(NFTinpt,*)epsilon(i,j)
                  read(NFTinpt,*)sigma(i,j)

                  epsilon(j,i)=epsilon(i,j)
                  sigma(j,i)=sigma(i,j)

  431          continue

  531       continue

         endif

c LJ type.
         ljpotlog=.false.
         if (pottype(1:2).eq.'LJ') then
            ljpotlog=.true.

            read(NFTinpt,*)ljcutoff
            do 501 i=1,nspec

               do 401 j=i,nspec

                  read(NFTinpt,*)epsilon(i,j)
                  read(NFTinpt,*)sigma(i,j)

                  epsilon(j,i)=epsilon(i,j)
                  sigma(j,i)=sigma(i,j)

  401          continue

  501       continue

         endif

c DW type.
         dwpotlog=.false.
         if (pottype(1:2).eq.'DW') then
            dwpotlog=.true.

            read(NFTinpt,*)ljcutoff
            do 521 i=1,nspec

               do 421 j=i,nspec

                  read(NFTinpt,*)epsilon(i,j)
                  read(NFTinpt,*)sigma(i,j)
                  read(NFTinpt,*)Adw
                  read(NFTinpt,*)alphadw
                  read(NFTinpt,*)r0dw

                  epsilon(j,i)=epsilon(i,j)
                  sigma(j,i)=sigma(i,j)

  421          continue

  521       continue

         endif

c DW2 type.
         dwpot2log=.false.
         if (pottype(1:3).eq.'2DW') then
            dwpot2log=.true.

            read(NFTinpt,*)ljcutoff
            do 522 i=1,nspec

               do 422 j=i,nspec

                  read(NFTinpt,*)epsilon(i,j)
                  read(NFTinpt,*)sigma(i,j)
                  read(NFTinpt,*)Adw
                  read(NFTinpt,*)alphadw
                  read(NFTinpt,*)r0dw

                  epsilon(j,i)=epsilon(i,j)
                  sigma(j,i)=sigma(i,j)

  422          continue

  522       continue

         endif

c BKS type.
         bkspotlog=.false.
         if (pottype(1:3).eq.'BKS') then
            bkspotlog=.true.

            do 701 i=1,nspec

               do 601 j=i,nspec


                  read(NFTinpt,*)ftalp(i,j)
                  read(NFTinpt,*)ftb(i,j)
                  read(NFTinpt,*)ftc(i,j)
                  read(NFTinpt,*)ftd(i,j)
                  read(NFTinpt,*)dddamp(i,j)
                  read(NFTinpt,*)dqdamp(i,j)

                  ftalp(j,i)=ftalp(i,j)
                  ftb(j,i)=ftb(i,j)
                  ftc(j,i)=ftc(i,j)
                  ftd(j,i)=ftd(i,j)
                  dddamp(j,i)=dddamp(i,j)
                  dqdamp(j,i)=dqdamp(i,j)

                  read(NFTinpt,*)epsilon(i,j)
                  read(NFTinpt,*)sigma(i,j)

                  epsilon(j,i)=epsilon(i,j)
                  sigma(j,i)=sigma(i,j)

  601          continue

  701       continue

         endif


c Tangney-Scandolo type.
         tslog=.false.    
         if (pottype(1:2).eq.'TS') then
            tslog=.true.    

            do 403 i=1,nspec

               do 503 j=i,nspec

                  read(NFTinpt,*)tsr(i,j)
                  read(NFTinpt,*)tsd(i,j)
                  read(NFTinpt,*)tsgamma(i,j)

                  tsr(j,i)=tsr(i,j)
                  tsd(j,i)=tsd(i,j)
                  tsgamma(j,i)=tsgamma(i,j)

                  dddamp(i,j)=1000000.0d0
                  dddamp(j,i)=1000000.0d0

  503          continue

  403       continue

         endif

c TMP - CIM assume 2 species only and that only the anion-cation term
c is CIM-like at the moment. Also assume self energy in the form
c of FH thesis.
         if (pottype(1:4).eq.'CIM1') then
            cim1log=.true.
            cim2log=.false.

c    CIM parameters....
            read(NFTinpt,*)ftalp(1,2)
            read(NFTinpt,*)ftbeta(1,2)
            read(NFTinpt,*)ftgamma(1,2)
            read(NFTinpt,*)ftb(1,2)
            read(NFTinpt,*)ftb2(1,2)
            read(NFTinpt,*)ftb3(1,2)
            read(NFTinpt,*)ftc(1,2)
            read(NFTinpt,*)ftd(1,2)
            read(NFTinpt,*)dddamp(1,2)
            read(NFTinpt,*)dqdamp(1,2)
            read(NFTinpt,*)selfB
            read(NFTinpt,*)selfgam
            read(NFTinpt,*)selfgam2
            read(NFTinpt,*)selfgam3

            ftalp(2,1)=ftalp(1,2)
            ftbeta(2,1)=ftbeta(1,2)
            ftgamma(2,1)=ftgamma(1,2)
            ftb(2,1)=ftb(1,2)
            ftb2(2,1)=ftb2(1,2)
            ftb3(2,1)=ftb3(1,2)
            ftc(2,1)=ftc(1,2)
            ftd(2,1)=ftd(1,2)
            dddamp(2,1)=dddamp(1,2)
            dqdamp(2,1)=dqdamp(1,2)

            do 21 i=1,nspec

               read(NFTinpt,*)ftalp(i,i)
               read(NFTinpt,*)ftb(i,i)
               read(NFTinpt,*)ftc(i,i)
               read(NFTinpt,*)ftd(i,i)
               read(NFTinpt,*)dddamp(i,i)
               read(NFTinpt,*)dqdamp(i,i)

   21       continue

         endif

c TMP - CIM assume 2 species only and that only the anion-cation term
c is CIM-like at the moment. Also assume self energy in the form
c of MW thesis.
c At the moment this is the CIM for which an AIM run can be combined.
         if ((pottype(1:4).eq.'CIM2').or.
     x       (pottype(1:4).eq.'DAIM').or.
     x       (pottype(1:5).eq.'QUAIM')) then
            cim2log=.true.
            cim1log=.false.

c    CIM parameters....
            read(NFTinpt,*)ftalp(1,2)
            read(NFTinpt,*)ftbeta(1,2)
            read(NFTinpt,*)ftgamma(1,2)
            read(NFTinpt,*)ftb(1,2)
            read(NFTinpt,*)ftb2(1,2)
            read(NFTinpt,*)ftb3(1,2)
            read(NFTinpt,*)ftc(1,2)
            read(NFTinpt,*)ftd(1,2)
            read(NFTinpt,*)dddamp(1,2)
            read(NFTinpt,*)dqdamp(1,2)
            read(NFTinpt,*)selfB
            read(NFTinpt,*)selfgam
            read(NFTinpt,*)extraalpha
            read(NFTinpt,*)extrab

            ftalp(2,1)=ftalp(1,2)
            ftbeta(2,1)=ftbeta(1,2)
            ftgamma(2,1)=ftgamma(1,2)
            ftb(2,1)=ftb(1,2)
            ftb2(2,1)=ftb2(1,2)
            ftb3(2,1)=ftb3(1,2)
            ftc(2,1)=ftc(1,2)
            ftd(2,1)=ftd(1,2)
            dddamp(2,1)=dddamp(1,2)
            dqdamp(2,1)=dqdamp(1,2)

            do 22 i=1,nspec

               read(NFTinpt,*)ftalp(i,i)
               read(NFTinpt,*)ftb(i,i)
               read(NFTinpt,*)ftc(i,i)
               read(NFTinpt,*)ftd(i,i)
               read(NFTinpt,*)dddamp(i,i)
               read(NFTinpt,*)dqdamp(i,i)

   22       continue

         endif

c other options - extended FT, LJ......

c Extended Fumi-Tosi (B/R and/or second exponential)     
c (i.e. see PbF2 potential)
         xftlog=.false.    
         if (pottype(1:3).eq.'XFT') then

            xftlog=.true.
            do 120 i=1,nspec

               do 130 j=i,nspec

                  read(NFTinpt,*)ftalp(i,j)
                  read(NFTinpt,*)ftb(i,j)
                  read(NFTinpt,*)ftalpx(i,j)
                  read(NFTinpt,*)ftbx(i,j)
                  read(NFTinpt,*)nrpower(i,j)
                  read(NFTinpt,*)ftc(i,j)
                  read(NFTinpt,*)ftd(i,j)
                  read(NFTinpt,*)dddamp(i,j)
                  read(NFTinpt,*)dqdamp(i,j)

                  ftalp(j,i)=ftalp(i,j)
                  ftb(j,i)=ftb(i,j)
                  ftalpx(j,i)=ftalpx(i,j)
                  ftbx(j,i)=ftbx(i,j)
                  ftc(j,i)=ftc(i,j)
                  ftd(j,i)=ftd(i,j)
                  dddamp(j,i)=dddamp(i,j)
                  dqdamp(j,i)=dqdamp(i,j)
                  nrpower(j,i)=nrpower(i,j)

  130          continue

  120       continue

         endif

c Dipolar aim.
         if ((pottype(1:4).eq.'DAIM').or.
     x       (pottype(1:5).eq.'QUAIM')) then
            read (NFTinpt,*) selfeps
            read (NFTinpt,*) selfC
         endif

c Quadrupolar AIM.
         if (pottype(1:5).eq.'QUAIM') then
            read (NFTinpt,*) selfquaim
            read (NFTinpt,*) selfH
         endif

c========================================================================
c dipole parameters. Search over species. If they are polarizable
c then read in polarizability, damping parameters, thermostat parameters.
         do i=1,nspec
            do j=1,nspec
               dampa(i,j)=0.0d0
               dampfac(i,j)=0.0d0
               nkdamp(i,j)=0
               fgb(i,j)=0.0d0
               fgc(i,j)=0.0d0
               nkfg(i,j)=0.0d0
            enddo
         enddo

         if (dippimlog) then
            do 15 i=1,nspec
               if (polarizablelog(i)) then

                  read (NFTinpt,*) alppolar(i)
                  read (NFTinpt,*) dipmass(i)
                  read (NFTinpt,*) pimtherm(i)
                  if (pimtherm(i)) then
                     read (NFTinpt,*) diptemp(i)
                     read (NFTinpt,*) diprelax(i)
                  endif
                  do 16 j=1,nspec
c could relax this if anion-anion functions made a come-back.
                     if (j.ne.i) then
                        read (NFTinpt,*) dampa(j,i)
                        read (NFTinpt,*) nkdamp(j,i)
                        read (NFTinpt,*) dampfac(j,i)
                        write(6,*)j,i,dampa(j,i),dampfac(j,i)
                     endif
   16             continue

               endif

  15        continue

         endif

         if (quadpimlog) then
            do 151 i=1,nspec
               if (polarizablelog(i)) then

                  read (NFTinpt,*) alppolar(i)
                  read (NFTinpt,*) Bpolar(i)
                  read (NFTinpt,*) Cpolar(i)
                  read (NFTinpt,*) dipmass(i)
                  read (NFTinpt,*) quadmass(i)
                  read (NFTinpt,*) pimtherm(i)
                  if (pimtherm(i)) then
                     read (NFTinpt,*) diptemp(i)
                     read (NFTinpt,*) quadtemp(i)
                     read (NFTinpt,*) diprelax(i)
                     read (NFTinpt,*) quadrelax(i)
                  endif
                  do 161 j=1,nspec
c could relax this if anion-anion functions made a come-back.
c NEED TO PUT IN QUADRUPOLE DAMPING TERMS HERE...
                     if (j.ne.i) then
                        read (NFTinpt,*) dampa(j,i)
                        read (NFTinpt,*) nkdamp(j,i)
                        read (NFTinpt,*) dampfac(j,i)
                        write(6,*)j,i,dampa(j,i),dampfac(j,i)
                     endif
  161             continue

               endif

  151       continue

         endif

c TMP - CIM
         if ((cimlog).or.(daimlog).or.(quaimlog)) then
            read(NFTinpt,*)selfmass
            read(NFTinpt,*)cimtherm
            if (cimtherm) then
               read(NFTinpt,*)deltatemp
               read(NFTinpt,*)deltarelax
            endif
         endif

         if ((daimlog).or.(quaimlog)) then
            read (NFTinpt,*) epsmass
            read (NFTinpt,*) aimtherm
            if (aimtherm) then
               read (NFTinpt,*) epstemp
               read (NFTinpt,*) epsrelax
            endif
         endif

         if (quaimlog) then
            read (NFTinpt,*) quaimmass
            read (NFTinpt,*) quaimtherm
            if (quaimtherm) then
               read (NFTinpt,*) quaimtemp
               read (NFTinpt,*) quaimrelax
            endif
         endif

         if ((fglog).or.(quadpimlog)) then
            do 67 i=1,nspec
               if (polarizablelog(i)) then
                  do 66 j=1,nspec
                     if (j.ne.i) then
                        read(NFTinpt,*)fgb(j,i)
                        read(NFTinpt,*)nkfg(j,i)
                        read(NFTinpt,*)fgc(j,i)
                     endif
   66             continue
               endif
   67       continue
         endif

      close(NFTinpt)

c
c Input correlation function matrix.
c
      open(ncfinpt,file=cffile,status='old')

         read(ncfinpt,*)ncfmat2
         read(ncfinpt,*)nkmod2
         read(ncfinpt,*)kmax_ds
         read(ncfinpt,*)nbin
         do 43 i=1,3
            read(ncfinpt,*)(xkvec2(i,j),j=1,ncfmat2)
   43    continue

         do 44 i=1,3
            do 45 j=1,ncfmat2
               xkvec2(i,j)=dble(xkvec2(i,j))
   45       continue
   44    continue

c more this as well
         read(ncfinpt,*)docfcalc
         if (docfcalc) then
            read(ncfinpt,*)npercf
            read(ncfinpt,*)ncorr
         endif
         read(ncfinpt,*)msdcalllog
         if (msdcalllog) then
            read(ncfinpt,*)nmsdcalltime
         endif

      close(ncfinpt)

c VARIABLE CELL CODE ADDITION:
c Read in the unit cell parameters.
c
      open (ncfinpt,file=cellfile,status='old')
c
c The cell vectors are read into the columns of h:
c
c     ( a_x  b_x  c_x )
c     ( a_y  b_y  c_y )
c     ( a_x  b_z  c_z )

	 if ((boxlenfrominptlog).or.(.not.restart)) then

cc            read (ncfinpt,*) h(1,1),h(1,2),h(1,3)
cc            read (ncfinpt,*) h(2,1),h(2,2),h(2,3)
cc            read (ncfinpt,*) h(3,1),h(3,2),h(3,3)
            read (ncfinpt,*) h(1,1),h(2,1),h(3,1)
            read (ncfinpt,*) h(1,2),h(2,2),h(3,2)
            read (ncfinpt,*) h(1,3),h(2,3),h(3,3)
            read (ncfinpt,*) boxlenx
            read (ncfinpt,*) boxleny
            read (ncfinpt,*) boxlenz

	 endif

         if (.not.restart) then
            read (ncfinpt,*) nunitcellx
            read (ncfinpt,*) nunitcelly
            read (ncfinpt,*) nunitcellz
            do 55 i=1,nspec
               read (ncfinpt,*) nunitcellpos(i)
               read (ncfinpt,'(a)') cellcoordfile(i)
   55       continue
            read (ncfinpt,*) a0
            read (ncfinpt,*) b0
            read (ncfinpt,*) c0

         endif

      close(ncfinpt)




c TMP
c  WHY IS THIS HERE? MOVE TO SETUP
c	if (.not.restart) then


c	obtain vol by calling dcell

        call dcell(h,b)

      vol3=boxlenx*boxleny*boxlenz*b(10)

      write(6,*)'vol3 ',vol3
      do 100 j=1,3
         do 110 k=1,3
            vg3(j,k)=0.0d0
	    hlab3(j,k)=h(j,k)
            hlab2(j,k)=h(j,k) 
 110     continue
 100  continue

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

            do 765 ii=1,5
               pzeta2(ii)=0.0d0
               pzeta3(ii)=0.0d0
               bzeta2(ii)=0.0d0
               bzeta3(ii)=0.0d0
               vpzeta1(ii)=0.0d0
               vpzeta2(ii)=0.0d0
               vpzeta3(ii)=0.0d0
               vbzeta1(ii)=0.0d0
               vbzeta2(ii)=0.0d0
               vbzeta3(ii)=0.0d0
  765       continue

            eps2=0.0d0
            eps3=0.0d0
            veps1=0.0d0
            veps2=0.0d0
            veps3=0.0d0

c	endif

c         print*,nth,nib,nab,ortho
c         print*,neqc,taub,taub2,pext

      if (dynam) then
c=====================================================================
c Input correlation function matrix.
         open(ncfinpt,file='dynmat.inpt',status='old')

            read(ncfinpt,*)amovefac
c 5.2.02 dummied to allow ftol and tol to be inputted in main runtime file.
            read(ncfinpt,*)toldu
            read(ncfinpt,*)ftoldu
            read(ncfinpt,*)nmat1
            read(ncfinpt,*)nmat2
            do  77 i=1,3
               read(ncfinpt,*)(amove(i,j),j=1,3)
   77       continue

         close(ncfinpt)

      endif

      write (6,*)
      write (6,*) '**** Run-time parameters read in ****'
      write (6,*)

      return

      end
