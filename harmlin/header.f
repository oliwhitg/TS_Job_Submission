      subroutine header

c************************************************************
c Generates a header file for each run containing
c all the information required to restart that run.
c************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'

      open(26,file=head,status='new')

      write(26,1040)
      write(26,1000)nsofar+1,nanion,ncation

      if (nsofar.ne.0) then
         write(26,1010)
         do 10 i=1,nsofar
            write(26,1020)nstpbrk(i)
   10    continue
      endif

      write(26,1040)
      write(26,1030)ntotstp,nrun-ntotstp,nrun
c
c Dump options selected.
c
      write(26,1040)
      if (restart) then

         write(26,1050)
 
      else
 
         if (nacl) then
            write(26,1060)
         else
            write(26,1070)
         endif
 
         if (randdis) then
            write(26,1080)dismag
         else
            write(26,1090)
         endif

         if (.not.conjrodlog) then
            write(26,1100)
         else
            write(26,1110)nrodstep
         endif

         if (rescalelog) then
            write(26,1130)
         endif

      endif
c
c Dump chemical parameters.
c
      write(26,1040)
      write(26,1200)
      write(26,1210)dens
      write(26,1220)rmm
      write(26,1222)rmm/dens
      write(26,1224)boxlen
      write(26,1226)cellvol
      do 20 i=1,nspec
         write(26,1230)i,chg(i),amass(i)
   20 continue
c
c Dump physical parameters.
c
      write(26,1300)
      write(26,1310)trantemp
      write(26,1330)diptemp

c=====================================================================
c Dump energy parameters (inc. FT-parameters).
      write(26,1400)
      write(26,1410)kmax
      write(26,1430)ftalp(1,1),ftalp(1,2)
     x             ,ftalp(2,1),ftalp(2,2)
      write(26,1440)
      write(26,1430)ftb(1,1),ftb(1,2)
     x             ,ftb(2,1),ftb(2,2)
      write(26,1450)
      write(26,1430)ftc(1,1),ftc(1,2)
     x             ,ftc(2,1),ftc(2,2)
      write(26,1460)
      write(26,1430)ftd(1,1),ftd(1,2)
     x             ,ftd(2,1),ftd(2,2)

c=====================================================================
c Dump thermostat data.
      write(26,1500)
      write(26,1510)relax
      write(26,1530)diprelax
      write(26,1540)deltarelax

c=====================================================================
c Dump anion data.
      write(26,1600)
      write(26,1615)alppolar
      write(26,1620)xk1

c=====================================================================
c Dump I/O information.
      write(26,1800)
      if (restart) then
         write(26,1810)rstin
      endif
      write(26,1820)rstout
      write(26,1860)crdout
      write(26,1870)velout
      write(26,1880)accnout
      write(26,1890)qout
      write(26,1900)npervel
      write(26,1910)
      write(26,1920)engout
      write(26,1930)engout2
      write(26,1940)npereng
      write(26,1950)zetout
      write(26,1960)fluxout
      write(26,1970)tempout
      write(26,1900)nperfri
	write(26,1900)npercell
      write(26,1980)

 1000 format(/,' Run number',i3,' for the system of:',/
     x       ,T30,i8,' anions ',/,T30,i8,' cations ',/)
 1010 format(/,'Previous runs were of length:')
 1020 format(T30,i8)
 1030 format(/,'Total number of previous steps = ',T30,i8,/
     x       ,'Number of steps in this run    = ',T30,i8,/
     x       ,'Final step count should be     = ',T30,i8)
 1040 format('=================================================')

c=====================================================================
c Run-time option format statements.
 1045 format(/,/,'The following options were selected:'
     x      ,'====================================')
 1050 format(/,'** Run was started from a restart file. **')
 1060 format(/,'** The run was started from the NaCl startup. **')
 1070 format(/,'** The run was started from the CdCl2 startup. **')
 1080 format(/,'** With a random displacement of ',f16.12,' a.u. **')
 1090 format(/,'** With the particles on the lattice. **')
 1100 format(/,'** There was no pseudo-conjugate gradient annealing **')
 1110 format(/,'** Pseudo-conjugate gradient annealing with: **'
     x      ,/,T30,i5,' steps.')
 1120 format(/,'** with  ',T30,i5,' steps of a rigid pre-run **')
 1130 format(/,'** Velocity rescaling prior to the main run **')

c=====================================================================
c Chemical parameters.
 1200 format(/,'Chemical parameters',/,'===================')
 1210 format(/,'Density.....................',T30,f16.12,' gcm^-3')
 1220 format('Relative molecular mass.....',T30,f16.12,' gmol^-1')
 1222 format('Thus, the molar volume is...',T30,f16.12,' cm^3mol^-1')
 1224 format('Boxlength...................',T30,f16.12,' a.u.')
 1226 format('Cell volume.................',T30,f16.6,' a.u.')
 1230 format('Ion',I2,' has:',/,T20,'Charge....',T30,f16.12,' a.u.',/
     x      ,T20,'mass......',T30,f14.6,' a.u.')

c=====================================================================
c System parameters.
 1300 format(/'System parameters',/,'=================')
 1310 format(/,'Translational temperature....',T30,f9.4,' K')
 1330 format('Charge temperature...........',T30,f9.4,' K')
 1340 format(/,'Timestep.....................',T30,f9.4,' K')

c=====================================================================
c Energy parameters.
 1400 format(/,'Energy parameters',/,'=================')
 1410 format(/,'kmax.........................',T30,i3)
 1420 format(/,'Fumi-Tosi parameters: (all a.u.)'
     x           ,/,'====================',/
     x           ,'i).Alpha = ',f9.4)
 1430 format(/,T5,'============================',/,T5
     x           ,'|   |    1           2     |',/,T5
     x           ,'============================',/,T5
     x           ,'|   |                      |',/,T5
     x           ,'| 1 |',T11,f9.4,T21,f9.4,T32,'|',/,T5
     x           ,'|   |                      |',/,T5
     x           ,'|   |                      |',/,T5
     x           ,'| 2 |',T11,f9.4,T21,f9.4,T32,'|',/,T5
     x           ,'|   |                      |',/,T5
     x           ,'============================')
 1440 format(/,'ii) B (short-range) terms.')
 1450 format(/,'iii) C (Van der Waals) terms.')
 1460 format(/,'iv) D (Van der Waals) terms.')

c=====================================================================
c Thermostat formats.
 1500 format(/,'Thermostat parameters',/,'=====================')
 1510 format(/,'Translational relaxation time.....',T35,f16.6,' a.u.')
 1530 format('Charge relaxation time............',T35,f16.6,' a.u.')
 1540 format('Delta relaxation time............',T35,f16.6,' a.u.')

c=====================================================================
c Anion formats.
 1600 format(/,'Anion parameters',/,'================')
 1615 format('Polarisalility....................',T35,f9.4,' a.u.')
 1620 format('Elastic constant..................',T35,f9.4,' a.u.')
 1630 format('Charge mass.......................',T35,f9.4,' a.u.')

c=====================================================================
c I/O formats.
 1800 format(/,'Input/Output Procedures',/
     x        ,'=======================')
 1810 format('Restart information read from..............',T46,A20)
 1820 format('Restart information written to.............',T46,A20)
 1830 format('Conjugate gradient information read from...',T46,A20)
 1840 format('Conjugate gradient information written to..',T46,A20)
 1860 format('Coordinates written out to.................',T46,A20)
 1870 format('Velocities written out to..................',T46,A20)
 1880 format('Accelerations written to...................',T46,A20)
 1890 format('Partial charge values written to...........',T46,A20)
 1900 format(T5,'The following are written every',I5,' steps')
 1910 format(T40,'Coordinates',/,T40,'Velocities',/,T40
     x          ,'Accelerations.')
 1920 format('Potential,translational & kinetic energies.',T46,A20)
 1930 format('Ficticious k.e,restoring energy,total E....',T46,A20)
 1940 format(T5,'Energies written every',I5,' steps')
 1950 format('Friction coefficients written to...........',T46,A20)
 1960 format('Heat fluxes written to.....................',T46,A20)
 1970 format('Kinetic temperatures written to............',T46,A20)
 1980 format(T40,'Friction coefficients',/,T40,'Heat fluxes',/,T40
     x          ,'Kinetic temperatures.',/,T40,'Partial charges.')

      close(26)

      write(6,*)
      write(6,*)'**** Header file written out and closed ****'
      write(6,*)

      return

      end
