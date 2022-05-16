      subroutine setup

c************************************************************
c   Sets up:  i)  auxillary variables (put everything into
c                atomic units.
c            ii)  I/O channels
c************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'

c     local characters
      CHARACTER*9 indi
      CHARACTER*8 namex,namey,namez
      CHARACTER*10 namesrx,namesry,namesrz
      CHARACTER*10 namexx,nameyy,namezz
     x            ,namexy,namexz,nameyz
      CHARACTER*12 namesrxx,namesryy,namesrzz
     x            ,namesrxy,namesrxz,namesryz
      CHARACTER*60 filenamex,filenamey,filenamez
      CHARACTER*60 filenamesrx,filenamesry,filenamesrz
      CHARACTER*60 filenamexx,filenameyy,filenamezz
     x            ,filenamexy,filenamexz,filenameyz
      CHARACTER*60 filenamesrxx,filenamesryy,filenamesrzz
     x            ,filenamesrxy,filenamesrxz,filenamesryz

      indi = '123456789'

      if ((dippimlog).or.(quadpimlog)) then

         namex = 'dipx.out'
         namey = 'dipy.out'
         namez = 'dipz.out'
         namesrx = 'dipsrx.out'
         namesry = 'dipsry.out'
         namesrz = 'dipsrz.out'

         do 963 i=1,nummon

            nnn=3*(i-1)

            filenamex=namex//indi(i:i)
            filenamey=namey//indi(i:i)
            filenamez=namez//indi(i:i)

            filenamesrx=namesrx//indi(i:i)
            filenamesry=namesry//indi(i:i)
            filenamesrz=namesrz//indi(i:i)

            open(77+nnn,file=filenamex,status='new')
            open(78+nnn,file=filenamey,status='new')
            open(79+nnn,file=filenamez,status='new')

c****            open(901+nnn,file=filenamesrx,status='new')
c****            open(902+nnn,file=filenamesry,status='new')
c****            open(903+nnn,file=filenamesrz,status='new')

  963    continue

c****         open(502,file='stress_pol.out',status='new')
c****         open(503,file='stress_pol2.out',status='new')
c****         open(504,file='stress_polsr.out',status='new')

      endif

      write(6,*) 'opening output files'
c
c Open all files for data dumping.
c
      open (21,file=engout,status='new')
      open (22,file=velout,status='new',form='unformatted')
      open (23,file=crdout,status='new',form='unformatted')
      open (24,file=accnout,status='new')
      open (25,file=engout2,status='new')
      open (27,file=qout,status='new')
      open (28,file=zetout,status='new')
      open (29,file=fluxout,status='new')
      open (30,file=tempout,status='new')
      open (31,file=fgout,status='new')

      open (32,file=cellengout,status='new')
      open (33,file=engtotout,status='new')
      open (34,file=pzetaout,status='new')

      open (36,file=momout,status='new')

      open (44,file=bzetaout,status='new')

      open (38,file=polengout,status='new')
      open (39,file=tdipout,status='new')
      open (40,file=fkeout,status='new')
      open(48,file='kvectors.out',status='new')
      open (59,file=eleout,status='new')
      open (98,file='disp.out',status='new',form='unformatted')

      open (51,file=diagpresout,status='new')
      open (52,file=presout,status='new')

      open (53,file=celllenout,status='new')
      open (54,file=cellvolout,status='new')

      open (55,file=cellanglesout,status='new')

      open (56,file=diagstressout,status='new')
      open (57,file=xxyyzzstressout,status='new')
      open (58,file=xyxzyzstressout,status='new')
      open (60,file=poscartout,status='new')
      open (61,file=cellboxout,status='new')

c      open (62,file=polstressout,status='new')
c      open (63,file=coulsrstressout,status='new')

c      open(500,file='stress_sr.out',status='new')
c      open(501,file='stress_coul.out',status='new')

      onethird=1.0d0/3.0d0
c
c Initialize random no. generator.
c
      dummy=ran1(-1)
c
c Atomic unit of mass (electron mass) in grams.
c
      emass=9.109534d-28
c
c Bohr radii in Angstroms.
c
      bohr=0.529177d0
c
c Boltzmann in a.u's (K^-1)
c
      boltz=3.166664d-006
      avo=6.022045d023
c
c Factors of pi to be used in later loops.
c
      pi=4.0d0*atan(1.0d0)
      pisq=pi*pi
      sqrpi=dsqrt(pi)
      twopi=2.0d0*pi
      fourpi=2.0d0*twopi
      eightpi=2.0d0*fourpi
      fourpisq=4.0d0*pi*pi
c
c  Calculate the charge 'force constant' from the polarisability
c  Modified for the new formalism.
c  Modified to contain quadrupole terms.
c

      do 550 i=1,nspec
         if ((dippimlog).or.(quadpimlog))
     x              xk1(i)=1.0d0/(2.0d0*alppolar(i))
  550 continue


      nskstep=nrun
      num=0
      nanion=0
      ncation=0
      l=0

      trantkb=boltz*trantemp
      chgsum=0.0d0

c TMP
      nanion=nsp(1)
      do 936 i=2,nspec
         ncation=ncation+nsp(i)
 936  continue

      do 20 i=1,nspec

         num=num+nsp(i)

c
c Convert masses from .inpt (gmol^-1) into internal (atomic) units.
c
         amass(i)=amass(i)/(avo*emass)
         recamass(i)=1.0d0/amass(i)
         hmass(i)=amass(i)/2.0d0

         do 30 j=1,nsp(i)

            l=l+1
            ntype(l)=i
            q(l)=chg(i)
            chgsum=chgsum+(chg(i)*chg(i))

   30    continue
c
c Calculate the variance for the Gaussian velocity set up.
c vartrans in atomic units.
c
         vartrans(i)=dsqrt(trantkb/amass(i))
   20 continue

c
c Set up the known damping coefficients.
c dampa(i,1) is the parameter that is read in.
c
cc      do i=1,ndampparmax
cc         do j=2,ndamptermmax
cc            dampa(i,j)=dampa(i,j-1)*dampa(i,1)/dble(j)
cc         enddo
cc      enddo
c
c Set correlation function pointers.
c
      ncorrtime=1
      ncorrcall=0
c
c Setup the startup configuration if required.
c
      if (.not.restart) then
c
c Read in the position matrices.
c
         do 7 i=1,nspec
            open (7,file=cellcoordfile(i),status='old')
               do 8 j=1,nunitcellpos(i)
                  read (7,*) bx(i,j),by(i,j),bz(i,j)
    8          continue
            close(5)
    7    continue
c
c Generate the coordinates.
c
c the .mat files that contain the bx,y,z arrays are in the
c cell-frame and to be visualised, should be transformed into
c the lab-frame.
c
         m=1
         do 11 i=1,nspec
            do 21 j=1,nunitcellpos(i)
               do 31 nx=1,nunitcellx
                  do 41 ny=1,nunitcelly
                     do 51 nz=1,nunitcellz

                        x(m)=a0*(nx-1)+bx(i,j)*a0
                        y(m)=b0*(ny-1)+by(i,j)*b0
                        z(m)=c0*(nz-1)+bz(i,j)*c0

                        m=m+1

   51                continue
   41             continue
   31          continue
   21       continue
   11    continue
c
c Output cell-frame coordinates.
c
         open (96,file='cellframe.out',status='new')
            do 94 i=1,num
               write (96,*) x(i),y(i),z(i)
   94       continue
         close (96)
c
c Setup displaced sublattice for phonon investigations:
c
         if (displace) then

            write (6,*) 'Randomly displacing ions off lattice.'

            do 70 i=1,num

               x(i)=x(i)+0.1d0*(ran1(idum)-0.5d0)
               y(i)=y(i)+0.1d0*(ran1(idum)-0.5d0)
               z(i)=z(i)+0.1d0*(ran1(idum)-0.5d0)

   70       continue

         endif
c
c Enforce periodic boundary conditions.
c
         do i=1,num

            if (x(i).lt.0.0d0) x(i)=x(i)+boxlenx
            if (y(i).lt.0.0d0) y(i)=y(i)+boxleny
            if (z(i).lt.0.0d0) z(i)=z(i)+boxlenz

            if (x(i).gt.boxlenx) x(i)=x(i)-boxlenx
            if (y(i).gt.boxleny) y(i)=y(i)-boxleny
            if (z(i).gt.boxlenz) z(i)=z(i)-boxlenz

         enddo
c
c Output coordinates, in lab. frame:
c
         open (50,file='setuppos.out',status='new')

            do 75 i=1,num
               write (50,*) h(1,1)*x(i)+h(1,2)*y(i)+h(1,3)*z(i),
     x                      h(2,1)*x(i)+h(2,2)*y(i)+h(2,3)*z(i),
     x                      h(3,1)*x(i)+h(3,2)*y(i)+h(3,3)*z(i)
   75       continue

         close (50)

      endif
c
c Set up the known dispersion damping terms:
c
      do 333 i=1,nspec
         do 444 j=1,nspec

            dddamp2(i,j)=dddamp(i,j)*dddamp(i,j)/2.0d0
            dddamp3(i,j)=dddamp2(i,j)*dddamp(i,j)/3.0d0
            dddamp4(i,j)=dddamp3(i,j)*dddamp(i,j)/4.0d0
            dddamp5(i,j)=dddamp4(i,j)*dddamp(i,j)/5.0d0
            dddamp6(i,j)=dddamp5(i,j)*dddamp(i,j)/6.0d0

            dqdamp2(i,j)=dqdamp(i,j)*dqdamp(i,j)/2.0d0
            dqdamp3(i,j)=dqdamp2(i,j)*dqdamp(i,j)/3.0d0
            dqdamp4(i,j)=dqdamp3(i,j)*dqdamp(i,j)/4.0d0
            dqdamp5(i,j)=dqdamp4(i,j)*dqdamp(i,j)/5.0d0
            dqdamp6(i,j)=dqdamp5(i,j)*dqdamp(i,j)/6.0d0
            dqdamp7(i,j)=dqdamp6(i,j)*dqdamp(i,j)/7.0d0
            dqdamp8(i,j)=dqdamp7(i,j)*dqdamp(i,j)/8.0d0

  444    continue
  333 continue
c
c Setup timestep and a few dependant variables. If we are
c to use this for more than just changing the timestep for
c annealing, the thermostat stuff should be done in timeset.f
c
      call timeset(dmaintime)
c
c Setup boxlength related variables:
c
      call boxreset

      call dcell(fullh,bh)

      halfboxminsq=(0.5d0*dmin1(bh(7),bh(8),bh(9)))**2

      rsqmax=rcut**2.0d0

c rsqmax now directly from .inpt. halfboxminsq from
c the minimum lengths in the simulation cell.
      if (rsqmax.gt.halfboxminsq) then
         write(6,*)'Cutoff violation',halfboxminsq,rsqmax
         stop
      endif

c
c Calculate rksqmax, the 1.1 is just to avoid rounding errors.
c
      do i=1,3
         do j=1,3
            fullhit(i,j)=fullhi(j,i)
         enddo
      enddo
      call dcell(fullhit,bh)


c
c Factors of the eta (the Ewald 'smearing parameter') to be used later.
c
      ksqmax=kmax*kmax

cc      eta=etainpt*boxlenrec
      eta=etainpt/(2.0d0*rcut)
      etasq=eta*eta
      etapi=2.0d0*eta/sqrpi
      etaconst=-1.0d0/(4.0d0*eta*eta)
      fouretacubpi=etapi*etasq*2.0d0

c modified 24.10.2000 to do variable kmax.
      if (.not.boxlenfromrstlog) call kset

c chgsum has already been calculated in setup
c
      chgcorrec=chgsum*etapi*0.5d0

c***** end of moved from boxreset
c
c Thermostat parameters.
c
      relaxsq=relax*relax
      relaxsqrec=1.0d0/relaxsq
      tfluxtran=0.0d0
      gtran=3.0d0*dble(num-1)*trantkb
      gtranrec=1.0d0/gtran
      gtrantim=gtran*dtime/trantemp
      gtrankin=2.0d0/(3.0d0*dble(num-1)*boltz)

      if ((dippimlog).or.(quadpimlog)) then

         do 440 i=1,nspec

            diptkb(i)=boltz*diptemp(i)
            diprlxsq(i)=diprelax(i)*diprelax(i)
            diprlxsqrec(i)=1.0d0/diprlxsq(i)
            tfluxdip(i)=0.0d0
            dipmassrec(i)=1.0d0/dipmass(i)
            gdip(i)=3.0d0*dble(nsp(i))*diptkb(i)
            gdiprec(i)=1.0d0/gdip(i)
            gdiptim(i)=gdip(i)*dtime/diptemp(i)
            gdipkin(i)=1.0d0/(3.0d0*dble(nsp(i))*boltz)

  440    continue

      endif

c
c Initialize mean-squared displacement output file.
c
      write (98) nmsdcalltime,dtime,nrun
      do 78 i=1,num
         write(98) ntype(i)
   78 continue
c
c Zero all extra degrees of freedom - just to be sure!
c This is done before any of them are read in.
c
      do i=1,num

         xmu(i)=0.0d0
         ymu(i)=0.0d0
         zmu(i)=0.0d0

         qvelx(i)=0.0d0
         qvely(i)=0.0d0
         qvelz(i)=0.0d0

         xmun(i)=0.0d0
         ymun(i)=0.0d0
         zmun(i)=0.0d0

         qvelnx(i)=0.0d0
         qvelny(i)=0.0d0
         qvelnz(i)=0.0d0

      enddo

      if (assign) then

         write (6,*) 'setting up anion dipoles'

c Assignment of anion dipoles:
c
         do 80 i=1,num
            xmu(i)=ran1(idum)-0.5d0
            ymu(i)=ran1(idum)-0.5d0
            zmu(i)=ran1(idum)-0.5d0
   80    continue

         open (228,file='randip.out')

            do 85 i=1,num
               write (228,*) real(xmu(i)),real(ymu(i)),real(zmu(i))
   85       continue

         close (228)

      endif

      write (6,*)
      write (6,*) '**** Internal variables set up ****'
      write (6,*)

c set up variables for the D-S S(k) calculation.
c.....sk and norm need zeroing in setup, dk2_ds, kmax_ds, rksqmax_ds
c.....and nbin need initializing
c...........................................

      do 600 i=1,nbinmax

         norm_ds(i)=0

         do 610 j=1,nspairs

            sk_ds(j,i)=0.0d0

  610    continue

  600 continue

      rksqmax_ds=float(kmax_ds)*float(kmax_ds)*pisq/halfboxminsq

c taken from `trans_chains.f'
      free=3.0d0*float(num-1)

c Set anisotropic barostat parameter equal to isotropic b. parameter
      relaxb=taub
      relaxb2=taub2
      W=(free+3.0d0)*trantkb*(relaxb**2.0d0)
      Wgo=(free+3.0d0)*trantkb*(relaxb2**2.0d0)/3.0d0

C.......THIS SEEMS TO BE A VERY BAD WAY OF CHOOSING THE BAROSTAT INERTIA
      Wrec=1.0d0/W
      write(6,*)'W,Wrec',W,Wrec,relaxb
      Wgorec=1.0d0/Wgo

      if (.not.nib) then
         W=0.0d0
         Wrec=0.0d0
         veps1=0.0d0
         veps2=0.0d0
         eps2=0.0d0
         eps3=0.0d0
      endif

      if (.not.nab) then
         Wgo=0.0d0
         Wgorec=0.0d0
         do 211 j=1,3
            do 311 k=1,3
               vg2(j,k)=0.0d0
 311        continue
 211     continue
      endif

      dom=0.0d0

      if (nib) then
         dom=1.0d0
      endif

      if (nab) then
         dom=9.0d0
         if (ortho) then
            dom=3.0d0
            do 22 j=1,3
               do 32 k=1,3
                  if(j.ne.k) then
                     vg2(j,k)=0.0d0
                     h2(j,k)=0.0d0
                     hi2(j,k)=0.0d0
                  endif
 32            continue
 22         continue
         endif
      endif

c............particle thermostat masses
 
      CUEp=gtran/relaxsqrec
      CUEp2=CUEp/free
      CUEprec=1.0d0/CUEp
      CUEp2rec=1.0d0/CUEp2
      if (.not.nth) then
         CUEprec=0.0d0
         CUEp2rec=0.0d0
      endif
 
c............barostat thermostat masses
 
      CUEb=dom*trantkb/relaxsqrec
      CUEb2=CUEb/dom
      CUEbrec=1.0d0/CUEb
      CUEb2rec=1.0d0/CUEb2

      if (.not.nth) then
         CUEbrec=0.0d0
         CUEb2rec=0.0d0
      endif
      if (.not.nib) then
         CUEbrec=0.0d0
         CUEb2rec=0.0d0
         CUEb2=0.0d0
         CUEb=0.0d0
      endif

      fkeself=0.0d0
      fkeeps=0.0d0
      fkequaim=0.0d0

c
c Zero all extra degrees of freedom - just to be sure!
c This is done before any of them are read in.
c
      do i=1,num

         xmu(i)=0.0d0
         ymu(i)=0.0d0
         zmu(i)=0.0d0

         qvelx(i)=0.0d0
         qvely(i)=0.0d0
         qvelz(i)=0.0d0

         qvelnx(i)=0.0d0
         qvelny(i)=0.0d0
         qvelnz(i)=0.0d0

      enddo

      do 653 i=1,num
         fixedarray(i)=.false.
  653 continue
      if (fixedlog) then
         do 654 i=1,nionsfixed
            fixedarray(nfixed(i))=.true.
  654    continue
      endif

      do 432 i=1,num
        write(76,*)fixedarray(i)
  432 continue

      rkb=   3.16669D-6
      if ((ljpotlog).or.(dwpotlog).or.(dwpot2log).or.(gcpotlog)) then

         do 500 i=1,nspec

            do 510 j=1,nspec

               epsilon(i,j)=epsilon(i,j)*rkb
               fourepsil(i,j)=4.0d0*epsilon(i,j)

               sigma(i,j)=sigma(i,j)/0.529177d0

               sig12(i,j)=sigma(i,j)**12.0d0
               sig6(i,j)=sigma(i,j)**6.0d0

         write(6,*)i,j,sig6(i,j)
  510       continue

  500    continue

      endif

      if (bkspotlog) then

         do 800 i=1,nspec

            do 810 j=1,nspec

cc               epsilon(i,j)=epsilon(i,j)*rkb
               fourepsil(i,j)=4.0d0*epsilon(i,j)

cc               sigma(i,j)=sigma(i,j)/0.529177d0

               sig30(i,j)=sigma(i,j)**30.0d0
               sig6(i,j)=sigma(i,j)**6.0d0

         write(6,*)i,j,sig6(i,j)
  810       continue

  800    continue

      endif

      return
      end
