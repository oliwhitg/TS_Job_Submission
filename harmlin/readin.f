      subroutine readin

c************************************************************
c   Reads in the data from the .inpt file.
c************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'
      include 'lj.inc'
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
      character*10 temp_pottype
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
      write(*,*) "in readin file"

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
      write(*,*) "read filenames"
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
         if (runtype2(1:3).eq.'cim') cimlog=.true.
         if (runtype2(1:4).eq.'daim') daimlog=.true.
         if (runtype2(1:5).eq.'quaim') quaimlog=.true.
         if (runtype2(1:3).eq.'mix') then
            mixlog=.true.
            read(NCinpt,*) num_potentials
         end if
           write(*,*) "run type" , runtype, runtype2
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
         write(*,*)etainpt
         read(NCinpt,*)rcut
         write(*,*)rcut
         read(NCinpt,*)conv
         write(*,*)conv

c thermostat parameters.
         read(NCinpt,*)nth
         if (nth) then
            read(NCinpt,*)relax
         endif

         read(NCinpt,*)nrscalelog
         if (nrscalelog) then
            read(NCinpt,*)nrscale
         else
            nrscale=1000000000
         endif

         read(NCinpt,*)nanneallog
         if (nanneallog) then
            read(NCinpt,*)nperanneal
         else
            nperanneal=1000000000
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

c LJ subsystem?
         read(NCinpt,*)ljsubsystem
         if (ljsubsystem) then
            read(NCinpt,*)fadeljpotlog
            if (fadeljpotlog) then
               read(NCinpt,*)nstepfadelj
               read(NCinpt,*)deltalambdalj
               xlambdalj=0.0d0
            else
               xlambdalj=1.0d0
            endif
         endif
         read(NCinpt,*)fixedlog
         if (fixedlog) then
         read(NCinpt,*)nionsfixed
            open(50,file='fixed.inpt',status='old')
               do 654 i=1,nionsfixed
                  read(50,*)nfixed(i)
            write(*,*)" mark7"
  654          continue
            close(50)
         endif
         read(NCinpt,*)fixedzlog
         if (fixedzlog) then
         read(NCinpt,*)nionsfixedz
            open(50,file='fixedz.inpt',status='old')
               do 655 i=1,nionsfixedz
                  read(50,*)nfixedz(i)
            write(*,*)" mark7"
  655          continue
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
         If (.not.mixlog) then 
           num_potentials =1  
         end if 

        ftlog=.false.    
         rvplog=.false.    
         ljpotlog=.false.
         sw2potlog=.false.
         sw3potlog=.false.
        tf2potlog=.false.
        xftlog=.false.    
        write(*,*) "start_reading_ pots"
        write(*,*) num_potentials, "num  pots"
      do 36 kkk=1,num_potentials

            write(*,*) kkk
            read (NFTinpt,'(a)')pottype
            write(*,*) pottype
         if (pottype(1:2).eq.'FT') then
            
        write(*,*) ftlog 
        write(*,*) "start_reading FT"
            ftlog=.true.    
          
        write(*,*) ftlog 
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

        write(*,*) "end FT"
         endif

c RVP type.
         if (pottype(1:3).eq.'RVP') then
            rvplog=.true.    
            epplog=.true.

            read(NFTinpt,*)rvpmax
            write(*,*)rvpmax
            do 201 i=1,nspec

               do 301 j=i,nspec

                  read(NFTinpt,*)rvph(i,j)
                  read(NFTinpt,*)rvpn(i,j)
                  read(NFTinpt,*)rvpr4(i,j)
                  read(NFTinpt,*)ftc(i,j)
                write(*,*)i,j,rvph(i,j)
                  write(*,*)rvpn(i,j)
                 write(*,*)rvpr4(i,j)
                  write(*,*)ftc(i,j)

                  ftc(j,i)=ftc(i,j)
                  rvph(j,i)=rvph(i,j)
                  rvpn(j,i)=rvpn(i,j)
                  rvpr4(j,i)=rvpr4(i,j)

  301          continue

  201       continue

         endif

c LJ type.
         if (pottype(1:2).eq.'LJ') then
            ljpotlog=.true.
         write(*,*) "start_readinglj"
            do 501 i=1,nspec

               do 401 j=i,nspec

                  read(NFTinpt,*)epsilon(i,j)
                  read(NFTinpt,*)sigma(i,j)

                  epsilon(j,i)=epsilon(i,j)
                  sigma(j,i)=sigma(i,j)

  401          continue

  501       continue
       rkb=   3.16669D-6
      do 291 i=1,nspec
         do 292 j=1,nspec
          write(*,*)i,j,sigma(i,j),epsilon(i,j)
          sigma(i,j)=sigma(i,j)/0.529177d0
          epsilon(i,j)=epsilon(i,j)*rkb

         sig12(i,j)=sigma(i,j)**12.0d0
         sig6(i,j)=sigma(i,j)**6.0d0

         fourepsil(i,j)=4.0d0*epsilon(i,j)
  292    continue
  291    continue

         write(*,*) "end_readinglj"
         endif

c SW type.
         if ((pottype(1:3).eq.'SW2').or.(pottype(1:3).eq.'SW3')) then
            sw2potlog=.true.

            do 701 i=1,nspec

               do 601 j=i,nspec

                  read(NFTinpt,*)swepsilon(i,j)
                  read(NFTinpt,*)swsigma(i,j)
                  read(NFTinpt,*)swa(i,j)
                  read(NFTinpt,*)swbigA(i,j)
                  read(NFTinpt,*)swB(i,j)
                  read(NFTinpt,*)swp(i,j)
                  read(NFTinpt,*)swq(i,j)

                  swepsilon(j,i)=swepsilon(i,j)
                  swsigma(j,i)=swsigma(i,j)
                  swa(j,i)=swa(i,j)
                  swbigA(j,i)=swbigA(i,j)
                  swB(j,i)=swB(i,j)
                  swp(j,i)=swp(i,j)
                  swq(j,i)=swq(i,j)

  601          continue


  701       continue

         endif

         if (pottype(1:3).eq.'SW3') then
            sw3potlog=.true.

            do 901 i=1,nspec

               do 801 j=i,nspec

                  read(NFTinpt,*)swgamma(i,j)
                  read(NFTinpt,*)swlambda(i,j)

                  swgamma(j,i)=swgamma(i,j)
                  swlambda(j,i)=swlambda(i,j)

       write(6,*)'gg',swgamma,swlambda
  801          continue

  901       continue

         endif

         write(6,*)'dfsdfsdf',pottype(1:3),kkk
         if (pottype(1:3).eq.'TF2') then
            tf2potlog=.true.
            write(*,*) "START READ TF2",tf2potlog
            do 702 i=1,nspec
              do 703 j=i, nspec

                  read(NFTinpt,*) bigA_type(i,j)
                  read(NFTinpt,*) bigB_type(i,j)
                  read(NFTinpt,*) tf_lamda1_type(i,j)
                  read(NFTinpt,*) tf_lamda2_type(i,j)
                  read(NFTinpt,*) tf_lamda3_type(i,j)
                  read(NFTinpt,*) tf_beta_type(i,j)
                  read(NFTinpt,*) tf_c_type(i,j)
                  read(NFTinpt,*) tf_d_type(i,j)
                  read(NFTinpt,*) tf_h_type(i,j)
                  read(NFTinpt,*) tf_n_type(i,j)
                  read(NFTinpt,*) bigR_type(i,j)
                  read(NFTinpt,*) bigD_type(i,j)
                  read(NFTinpt,*) ftc(i,j)
                  read(NFTinpt,*)dddamp(i,j)

                  ftc(j,i)=ftc(i,j)
                  dddamp(j,i)=dddamp(i,j)


  703       continue
  702       continue

            write(*,*) "end READ TF2"
         endif

         if (pottype(1:4).eq.'harm') then
            harmpotlog=.true.
            write(*,*) "START READ Harm"
            do 706 i=1,nspec
               do 707 j=i, nspec

                  read(NFTinpt,*)r0harm(i,j)
                  read(NFTinpt,*)xkharm(i,j)
                  r0harm(j,i)=r0harm(i,j)
                  xkharm(j,i)=xkharm(i,j)

  707          continue
  706       continue
            open(777,file='harmpairs.dat',status='old')
               read(777,*)npairs
               do 777 i=1,npairs
                  read(777,*)npairident1(i),npairident2(i)
  777          continue
            close(777)
         endif
         if (pottype(1:7).eq.'harmrep') then
            harmpotlog2=.true.
            write(*,*) "START READ Harm rep"
            rkb=3.16669D-6
            do 806 i=1,nspec
               do 807 j=i,nspec

                  read(NFTinpt,*)epsilon(i,j)
                  read(NFTinpt,*)sigma(i,j)

                  epsilon(j,i)=epsilon(i,j)
                  sigma(j,i)=sigma(i,j)

  807          continue
  806       continue
      rkb=   3.16669D-6
      do 491 i=1,nspec
         do 492 j=1,nspec
cc          sigma(i,j)=sigma(i,j)/0.529177d0
          epsilon(i,j)=epsilon(i,j)*rkb

         sig24(i,j)=sigma(i,j)**24.0d0
         sig12(i,j)=sigma(i,j)**12.0d0
            cut(i,j)=(2.0d0**(1.0d0/12.0d0))*sigma(i,j)
            cut(j,i)=(2.0d0**(1.0d0/12.0d0))*sigma(j,i)

            write(6,*)'cut',i,j,cut(i,j),sigma(i,j)
         fourepsil(i,j)=4.0d0*epsilon(i,j)

  492    continue
  491    continue
         endif

cc 2021. Linear potential for bilayer relaxation.
         if (pottype(1:7).eq.'harmlin') then
            harmpotlog3=.true.
            write(*,*) "START READ Harm linear"
            do 816 i=1,nspec
               do 817 j=i,nspec

                  read(NFTinpt,*)u0lin(i,j)
                  read(NFTinpt,*)cut(i,j)

                  cut(j,i)=cut(i,j)
                  u0lin(j,i)=u0lin(i,j)

  817          continue
  816       continue
            do 493 i=1,nspec
               do 494 j=1,nspec

                  xklin(i,j)=-u0lin(i,j)/cut(i,j)

            write(6,*)'xklin',i,j,xklin(i,j)

  494          continue
  493       continue
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

   36  continue
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

      if (docfcalc) then
c=====================================================================
c Input light scattering correlation function matrix.
         open(ncfinpt,file=lcffile,status='old')

            read(ncfinpt,*)hyperB
            read(ncfinpt,*)hypergamma
            read(ncfinpt,*)polarundist
	    do 46 i=1,nspec
               read(ncfinpt,*) xlipol(i)
   46	    continue
            do 50 i=1,nspec

               do 60 j=i,nspec

                  read(ncfinpt,*)asr(i,j)
                  read(ncfinpt,*)bsr(i,j)
                  read(ncfinpt,*)csr(i,j)
                  read(ncfinpt,*)dsr(i,j)

                  asr(j,i)=asr(i,j)
                  bsr(j,i)=bsr(i,j)
                  csr(j,i)=csr(i,j)
                  dsr(j,i)=dsr(i,j)

   60          continue

   50       continue

            do 65 i=1,nspec
               read(ncfinpt,*)sig(i)
   65       continue

         close(ncfinpt)

      endif
c
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
