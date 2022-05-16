cc      implicit double precision(a-h,o-z)

      parameter (nummax=5000,nspmax=2,nksqtot=1000,nkmodmax2=9,
     x           ncfmatmax2=18,ncorrtimemax=1,nkmodmax=6,
     x           ncfmatmax=9,nruns=100,nktot=35,ndampparmax=3,
     x           ndamptermmax=6,nunimax=200,nchains=5)
c................for light scattering
      parameter (nsymmax=6)
      parameter (nindexmax=10)
      parameter (ncontribmax=4)
c   New D-S S(k) varaibles.
      parameter (nbinmax=1000)
      parameter (nktot_ds=20)
      parameter (nspairs =6)
c*********nspairs=(nspmax*(nspmax+1))/2


      logical cdiplog,cdprev
      logical nth
      logical ortho
      logical nib
      logical nab
      logical outp
      logical forfl
      logical backsteplog,boxlog,docfcalc,annealrad,dipcalclog,
     x        assign,randdis,displace,veldumplog,crddumplog,chgdumplog,
     x        fulldumplog,chgveldulog,endrun,ioniden,irreadlog,
     x        shutdown,pereng,perfric,pervel,rdfcall,restart,nacl,cdcl2,
     x        mcl3,pressurerst,boxlenkeep,rim,pim,cim,aim,qim,quaim,
     x        rescalelog,velinit,conjstart,
     x        debug,quaddamp,moveions,iontherm,pimtherm,qimtherm,
     x        cimtherm,aimtherm,dynam,ddpot,relaxconfig,quaimtherm,
     x        percell
c
      logical conjrodlog,conjgradlog,dippimlog,quadpimlog
     x             ,cimlog,rimlog,aimlog,epplog,cim1log,cim2log
     x             ,msdcalllog,rvplog
      logical boxlenfromdenslog
     x       ,boxlenfromrstlog
     x       ,boxlenfrominptlog
      logical nrscalelog,nanneallog
      logical polarizablelog
      logical tramplog,pramplog,tramprescalelog
      logical xftlog
      logical ftlog
      logical daimlog,quaimlog
      logical fglog

      character*20 engout,engout2,velout,crdout,accnout,rstin,crstin,
     x             rstout,crstout,head,qout,zetout,fluxout,tempout,
     x             eleout,untout,cfout,cfin,prerstin,prerstout,fileout,
     x             cellcoordfile
      character*20  diagstressout,xxyyzzstressout,xyxzyzstressout
     x             ,poscartout,cellboxout,polstressout
     x             ,coulsrstressout,momout
      character*20  cellengout,engtotout,pzetaout,bzetaout
     x             ,fkeout,polengout,tdipout,diagpresout
     x             ,presout,celllenout,cellvolout,cellanglesout
      character*200 fgout

c      double precision ar(nspmax,ncfmatmax2,nkmodmax2)
c     x                ,ai(nspmax,ncfmatmax2,nkmodmax2)

      common sk_ds(nspairs,0:nbinmax),rksqmax_ds,dk2_ds,
     x      norm_ds(0:nbinmax),kmax_ds,nbin

c===============================================
c Double precisions in common.
cc      common/didtm2/txxli(nummax),tyyli(nummax),tzzli(nummax)
cc     x             ,txyli(nummax),txzli(nummax),tyzli(nummax)
c.........contains DID polarizabilities

      common/elcflgu/elecxu(nummax)
     x             ,elecyu(nummax)
     x             ,eleczu(nummax)
c.........contains electric fields without the sr terms

      common/elegrd/exx(nummax),eyy(nummax),ezz(nummax)
     x             ,exy(nummax),exz(nummax),eyz(nummax)
      common/egrdsr/exxsr(nummax),eyysr(nummax),ezzsr(nummax)
     x             ,exysr(nummax),exzsr(nummax),eyzsr(nummax)
c.........contains electric field gradients

      common/srtrms/srxx(nummax),sryy(nummax),srzz(nummax)
     x             ,srxy(nummax),srxz(nummax),sryz(nummax)
c.........contains short-range polarizabilities
 
      common/light1/cohacc(nsymmax,ncontribmax)
      common/light2/lnorm(nsymmax,nindexmax,0:ncorrtimemax) 
      common/light3/cohav(nsymmax,ncontribmax)
     x             ,scf(nsymmax,nindexmax,0:ncorrtimemax)
      common/light4/storecohacc(nsymmax,ncontribmax,ncorrtimemax)

      common/light5/term(nsymmax,ncontribmax)

      common/isoarrays/cisokr(5,ncfmatmax2),cisoki(5,ncfmatmax2)

      common/isocf/scfiso(15,ncfmatmax2,0:ncorrtimemax)

c================================================
c Double precisions in common.
      common/hyperp/hyperB,hypergamma,polarundist 

      common/lgttms/root3,root12,root1twodiv3gam,twodiv3gamma
     x             ,root3div3gam,gammadiv3,bhyperpol3

      common/radius/sig(nspmax)
      common/srligt/asr(nspmax,nspmax),bsr(nspmax,nspmax)
     x             ,csr(nspmax,nspmax),dsr(nspmax,nspmax)

      common/polarr/xlipol(nspmax)

c................Integers in lightscattering

      common/avpolar/normav

      common/ffffff/storeisokr(5,ncfmatmax2,0:ncorrtimemax),
     c          storeisoki(5,ncfmatmax2,0:ncorrtimemax)

c.................end of light scattering

c Double precision arrays in common.
      common/FTpara/ftalp(nspmax,nspmax),ftb(nspmax,nspmax)
     x             ,ftc(nspmax,nspmax),ftd(nspmax,nspmax)
      common/FTpar2/ftbeta(nspmax,nspmax),ftb2(nspmax,nspmax)
      common/FTpar3/ftgamma(nspmax,nspmax),ftb3(nspmax,nspmax)
      common/acceln/ax(nummax),ay(nummax),az(nummax)
      common/accum/sreng(nspmax,nspmax)
     x            ,ddeng(nspmax,nspmax)
     x            ,dqeng(nspmax,nspmax)
      common/asquad/asquadxx(nummax),asquadyy(nummax),asquadzz(nummax)
     x             ,asquadxy(nummax),asquadxz(nummax),asquadyz(nummax)
      common/assdip/asdipx(nummax),asdipy(nummax),asdipz(nummax)
      common/cellm/h(3,3),hi(3,3),fullh(3,3),fullhi(3,3),fullhit(3,3),
     x             dcellinfo(10),bh(10)
      common/charge/q(nummax),chg(nspmax)
      common/chgvex/qvelx(nummax),qvelnx(nummax)
      common/chgvey/qvely(nummax),qvelny(nummax)
      common/chgvez/qvelz(nummax),qvelnz(nummax)
      common/coords/x(nummax),y(nummax),z(nummax)
      common/dampCD/dampa(ndampparmax,ndampparmax)
     x             ,dampfac(ndampparmax,ndampparmax)

      common/dipupd/xmun(nummax),ymun(nummax),zmun(nummax)
      common/dipval/xmu(nummax),ymu(nummax),zmu(nummax)
      common/dispda/dddamp(nspmax,nspmax)
     x             ,dddamp2(nspmax,nspmax)
     x             ,dddamp3(nspmax,nspmax)
     x             ,dddamp4(nspmax,nspmax)
     x             ,dddamp5(nspmax,nspmax)
     x             ,dddamp6(nspmax,nspmax)
      common/dispdb/dqdamp(nspmax,nspmax)
     x             ,dqdamp2(nspmax,nspmax)
     x             ,dqdamp3(nspmax,nspmax)
     x             ,dqdamp4(nspmax,nspmax)
     x             ,dqdamp5(nspmax,nspmax)
     x             ,dqdamp6(nspmax,nspmax)
     x             ,dqdamp7(nspmax,nspmax)
     x             ,dqdamp8(nspmax,nspmax)
      common/dynam/amove(3,3),amovefac
      common/elcflg/elecx(nummax),elecy(nummax),elecz(nummax)
 
      common/engftdots/engft1dotx(nummax),engft1doty(nummax),
     x                 engft1dotz(nummax),engft2dotx(nummax),
     x                 engft2doty(nummax),engft2dotz(nummax),
     x                 engft3dotx(nummax),engft3doty(nummax),
     x                 engft3dotz(nummax)
      common/engftdot2/engftdotxx(nummax,3),engftdotyy(nummax,3),
     x                 engftdotzz(nummax,3),engftdotxy(nummax,3),
     x                 engftdotxz(nummax,3),engftdotyz(nummax,3)
      common/engftdotsoo/engft1dotxoo(nummax),engft1dotyoo(nummax),
     x                 engft1dotzoo(nummax)
      common/engset/ak(nksqtot)
      common/epsilon/epsilonx(nummax),epsilony(nummax),epsilonz(nummax)
      common/epsilonn/epsilonxn(nummax),epsilonyn(nummax)
     x               ,epsilonzn(nummax)
      common/epsvel/epsvelx(nummax),epsvely(nummax),epsvelz(nummax)
      common/epsveln/epsvelnx(nummax),epsvelny(nummax),epsvelnz(nummax)
      common/forces/frrx(nummax),frry(nummax),frrz(nummax)
      common/ircorr/cfaa(ncfmatmax,nkmodmax,0:ncorrtimemax)
     x             ,cfcc(ncfmatmax,nkmodmax,0:ncorrtimemax)
     x             ,cfca(ncfmatmax,nkmodmax,0:ncorrtimemax)
      common/irstor/storeara(ncfmatmax,nkmodmax,ncorrtimemax)
     x             ,storeaia(ncfmatmax,nkmodmax,ncorrtimemax)
     x             ,storearc(ncfmatmax,nkmodmax,ncorrtimemax)
     x             ,storeaic(ncfmatmax,nkmodmax,ncorrtimemax)
      common/kvects/xkvec(3,ncfmatmax)
      common/kvect2/xkvec2(3,ncfmatmax2)
      common/masses/amass(nspmax),hmass(nspmax),recamass(nspmax)
      common/msdarray/xdisp(nummax),ydisp(nummax),zdisp(nummax)
      common/newpos/rxn(nummax),ryn(nummax),rzn(nummax)
      common/rdfarr/rdftot(0:1000),rdfpart(0:1000,nspmax,nspmax)
      common/selfde/delta(nummax),deltan(nummax)
     x             ,deltavel(nummax),deltaveln(nummax)
      common/selfe2/engft1(nummax),engft2(nummax),engft3(nummax)
      common/selfe3/engft1oo(nummax),engft2oo(nummax),engft3oo(nummax)
      common/sk1/arsk(nspmax,nksqtot),aisk(nspmax,nksqtot)
      common/sk2/sk(nspmax,nspmax,nksqtot)
c      common/sk2/skaa(nksqtot),skac(nksqtot),skcc(nksqtot)
c      common/skcalc/skanan(ncfmatmax2,nkmodmax2)
c     x             ,skcatcat(ncfmatmax2,nkmodmax2)
c     x             ,skancat(ncfmatmax2,nkmodmax2)

      common/skcalc/skc(nspmax,nspmax,ncfmatmax2,nkmodmax2)
      common/srdip/srdipx(nummax),srdipy(nummax),srdipz(nummax)
      common/srquad/srquadxx(nummax),srquadyy(nummax),srquadzz(nummax)
     x             ,srquadxy(nummax),srquadxz(nummax),srquadyz(nummax)
      common/quadru/quadxx(nummax),quadyy(nummax),quadzz(nummax)
     x             ,quadxy(nummax),quadxz(nummax),quadyz(nummax)
      common/quadrn/quadxxn(nummax),quadyyn(nummax),quadzzn(nummax)
     x             ,quadxyn(nummax),quadxzn(nummax),quadyzn(nummax)
      common/quadve/quadvelxx(nummax),quadvelyy(nummax)
     x             ,quadvelzz(nummax),quadvelxy(nummax)
     x             ,quadvelxz(nummax),quadvelyz(nummax)
      common/quadvn/quadvelxxn(nummax),quadvelyyn(nummax)
     x             ,quadvelzzn(nummax),quadvelxyn(nummax)
     x             ,quadvelxzn(nummax),quadvelyzn(nummax)
      common/quaim/quaimxx(nummax),quaimyy(nummax),quaimzz(nummax)
     x            ,quaimxy(nummax),quaimxz(nummax),quaimyz(nummax)
      common/quaimn/quaimxxn(nummax),quaimyyn(nummax),quaimzzn(nummax)
     x             ,quaimxyn(nummax),quaimxzn(nummax),quaimyzn(nummax)
      common/quaimvel/quaimvelxx(nummax),quaimvelyy(nummax)
     x               ,quaimvelzz(nummax),quaimvelxy(nummax)
     x               ,quaimvelxz(nummax),quaimvelyz(nummax)
      common/quaimveln/quaimvelxxn(nummax),quaimvelyyn(nummax)
     x                ,quaimvelzzn(nummax),quaimvelxyn(nummax)
     x                ,quaimvelxzn(nummax),quaimvelyzn(nummax)
      common/untcll/bx(nspmax,nunimax)
     x             ,by(nspmax,nunimax)
     x             ,bz(nspmax,nunimax)
      common/varian/vartrans(nspmax)
      common/veloc2/vxn(nummax),vyn(nummax),vzn(nummax)
      common/veloci/vx(nummax),vy(nummax),vz(nummax),vmag(nummax)
      common/workspace/work(nummax,4)

c       new dp arrays for variable cell

	common/newforces/frrx3(nummax),frry3(nummax),frrz3(nummax)
	common/newacclns/ax3(nummax),ay3(nummax),az3(nummax)
      	common/cellprev/h2(3,3),hi2(3,3),vg2(3,3),b2(10),b(10),bee(10)
	common/cell/h3(3,3),hi3(3,3),vg3(3,3),b3(10)
c	common/hcell/h(3,3),hi(3,3),deth
    	common/extr/hih(3,3)
	common/extras/com1(3,3),com2(3,3),com3(3,3),Fgo2(3,3),Fgo3(3,3)
	common/extras2/vgsq2(3,3),vgsq3(3,3),comexp(3,3),hhi(3,3)
      	common/extras3/hcom(3,3)
      	common/extras4/vg2sq(3,3),vg3sq(3,3)
      	common/extras5/hh(3,3),hu(3,3),hui(3,3)
      	common/ptensor/pint2(3,3),pint3(3,3)
      	common/iden/iden(3,3)
	common/hlab/hlab2(3,3),hlab3(3,3)
	common/clare/clare_old(8,8),clare_new(8,8)

        common/cellcor/xcell(8),ycell(8),zcell(8)
        common/cellcor2/xcellab(8),ycellab(8),zcellab(8)
c
c Integer arrays in common.
      common/damppar/nkdamp(ndampparmax,ndamptermmax)
      common/runstp/nstpbrk(nruns)
      common/sknorm/nskcount(nksqtot)
      common/typsp1/nsp(nspmax)
      common/typsp2/ntype(nummax)
      common/untcl2/nunitcellpos(nspmax)
c
c Character arrays in common.
      common/cellcd/cellcoordfile(nspmax)
c
c Double precisions in common.
      common/annealdt/dtanneal,dmaintime
      common/boltzm/boltz,trantkb,avo
      common/boxdim/boxlen,halfbox,halfboxrec,boxlenrec
      common/boxdi1/boxlenx,halfboxx,halfboxxrec
      common/boxdi2/boxleny,halfboxy,halfboxyrec
      common/boxdi3/boxlenz,halfboxz,halfboxzrec
      common/boxmov/boxmovemax,boxmovefac,boxbeta,boxtemp
      common/celvar/cellvol,a0,b0,c0

c	changes for catdip
      common/chgres/xk1(nspmax),xk2(nspmax),xk3(nspmax)
     x             ,xk4(nspmax),alppolar(nspmax),Bpolar(nspmax)
     x             ,Cpolar(nspmax)
      common/chgsmm/chgcorrec,chgsum
      common/conjto/conjtol
      common/delzet/deltazetanew,deltazeta,gdeltakin
     x             ,gdeltarec,deltarlxsqrec,deltarelax
     x             ,deltatemp,dtdelta2,tdelta
      common/degfre/gtran,gtranrec
      common/degfre2/gdip(nspmax),gdiprec(nspmax)
     x              ,gquad(nspmax),gquadrec(nspmax)
      common/degkin/gtrankin
      common/degkin2/gdipkin(nspmax),gquadkin(nspmax)
      common/degtim/gtrantim
      common/degtim/gdiptim(nspmax),gquadtim(nspmax)
      common/dipdat/dipmass(nspmax),dipmassrec(nspmax)
     x             ,quadmass(nspmax),quadmassrec(nspmax)
      common/dipthm/dipzeta(nspmax),dtdip2(nspmax)
     x             ,diprlxsq(nspmax),diprelax(nspmax)
     x             ,diptkb(nspmax),diptemp(nspmax)
     x            ,tfluxdip(nspmax),tdip(nspmax)
     x             ,diprlxsqrec(nspmax)
      common/fakeke/fketot
      common/fakeke2/fke(nspmax),fkequad(nspmax),reseng(nspmax),
     x              quadeng(nspmax),dipquadeng(nspmax)
     x              ,dipsqeng(nspmax)
      common/quadthm/quadzeta(nspmax),dtquad2(nspmax)
     x              ,quadrlxsq(nspmax),quadrelax(nspmax),
     x               quadtkb(nspmax),quadtemp(nspmax)
     x              ,tfluxquad(nspmax),tquad(nspmax)
     x              ,quadrlxsqrec(nspmax)
      common/zechne/dipzetanew(nspmax),quadzetanew(nspmax)
      common/zetrne/zetanew
      common/monitr/nmon(nspmax),nummon
      common/sizeof2/numb

      common/displc/dismag
      common/energy/engpetot,tranke,encoul,tranko
      common/energy2/qmueng,xmumueng,qquadeng,cpeacc
      common/epseng/fkeeps,epsselfengtot
      common/epspar/selfeps,selfC,epsmass,epsmassrec
      common/epszeta/gepskin,teps,epstemp,epszeta,epszetanew,dteps2,
     x               gepsrec,epsrelax,epsrlxsqrec
      common/etainn/etainpt

      common/fluxxx/fluxtran,tfluxtran
      common/forcfa/fac
      common/frictn/zeta,dtzeta2,relax,relaxsq,relaxsqrec
      common/msdisp/dmsdan,dmsdcat
      common/numbrs/onethird
      common/outpt1/tkin,tempfac,tkinres
      common/phondisp/xmove,ymove,zmove
      common/physic/dens,rmm
      common/picons/pi,sqrpi,twopi,fourpi,fourpisq,eightpi,tworootpi
      common/picon2/twopibox,fourpicell,cellvol3rec
      common/picon3/hzero1,hzero2,hzero3
      common/picon4/twopiboxx,twopiboxy,twopiboxz
      common/prespar1/erfcracc,qdipacc,dipdipacc,qquadacc,dipquadacc,
     x                quadquadacc,toteng,srdipeng1,srdipeng2,qqacc
      common/pressu/pressreq,accratreq
      common/pressure2/accrat

      common/quaimeng/fkequaim,quaimselfengtot
      common/quaimpar/selfquaim,selfH,quaimmass,quaimmassrec
      common/quaimzeta/gquaimkin,tquaim,quaimtemp,quaimzeta,
     x                 quaimzetanew,dtquaim2,gquaimrec,quaimrelax,
     x                 quaimrlxsqrec
      common/random/dummy
      common/selfpr/selfgam,selfgam2,selfgam3,selfB,selfmass,selfmassrec
      common/selffk/fkeself,selfengtot
      common/smeary/eta,etasq,etapi,etaconst,fouretacubpi,rsqmax,rksqmax
      common/stress/stsrxx,stsrxy,stsrxz,stsryy,stsryz,stsrzz,
     x              stcxx,stcxy,stcxz,stcyy,stcyz,stczz,
     x              stpxx,stpxy,stpxz,stpyy,stpyz,stpzz,
     x              stp2xx,stp2xy,stp2xz,stp2yy,stp2yz,stp2zz,
     x              stpsrxx,stpsrxy,stpsrxz,stpsryy,stpsryz,stpsrzz
      common/strqud/stpqquadxx,stpqquadyy,stpqquadzz
     x             ,stpqquadxy,stpqquadxz,stpqquadyz
     x             ,stpdipquadxx,stpdipquadyy,stpdipquadzz
     x             ,stpdipquadxy,stpdipquadxz,stpdipquadyz
     x             ,stpquadquadxx,stpquadquadyy,stpquadquadzz
     x             ,stpquadquadxy,stpquadquadxz,stpquadquadyz
      common/temkin/tke,trke
      common/temper/trantemp
      common/timest/dtime,dtimerec,dtime2
      common/velave/vansumx,vansumy,vansumz
     x             ,vcatsumx,vcatsumy,vcatsumz
      common/virial/xviraccsrdip,yviraccsrdip,zviraccsrdip,
     x              xviraccftsrc,yviraccftsrc,zviraccftsrc
      common/virial2/virtot,virtotx,virtoty,virtotz,srdipvir,srvirft
      common/virial3/virx,viry,virz
      common/xtraft/extraalpha,extrab

c	New dps in common

      common/pfrictn/eps1,eps2,eps3,veps1,veps2,veps3,vol2,vol3

      common/baro/presfinal,presfinalx,presfinaly,presfinalz
      common/pcont/taub,taub2,pext,free,relaxb,relaxb2
      common/deter/deth
      common/thermop/vpzeta1(nchains),vpzeta2(nchains)
     x                ,vpzeta3(nchains),pzeta2(nchains)
     x                ,pzeta3(nchains)
      common/thermob/vbzeta1(nchains),vbzeta2(nchains)
     x                ,vbzeta3(nchains),bzeta2(nchains)
     x                ,bzeta3(nchains)

c
c Integers in common.
      common/accept/naccept,nupaccept
      common/anneal/nperanneal
      common/celll3/nunitcellx,nunitcelly,nunitcellz
      common/corrf1/ncorrtime,ncorr,ncorrcall
      common/corrpo/npercf
      common/kvecs2/nkmod,ncfmat,nkmod2,ncfmat2

      common/msd/nmsdcalltime
      common/npres/nboxmove
      common/npres2/nperpres
      common/presse/nprescall
      common/prestp/nrodstep
      common/rdfcal/nperrdf,npperrdf,npersk
      common/rdfnrm/nrdfcall,nskcall
      common/rscale/nrscale
      common/runste/nsofar,ntotstp
      common/runtim/nrun,npereng,nstep,npervel,nperfri
      common/runtim2/npercell

      common/sizeof/num,nanion,ncation,nionunit

      common/sk3/nskstep
      common/specie/nspec
      common/vector/conv,kmax,ksqmax,kmaxx,kmaxy,kmaxz

c       new integers for variable cell
      common/upflag/flagup 
      common/celleq/neqc
c
c Logicals in common.
      common/backsp/backsteplog
      common/boxexp/boxlog
      common/cfcalc/docfcalc
      common/cimlol/annealrad
      common/diplog/dipcalclog
      common/dipset/assign
      common/displm/randdis,displace
      common/dumlog/veldumplog,crddumplog
     x             ,chgdumplog,fulldumplog
     x             ,chgveldulog
      common/endrun/endrun
      common/ionlog/ioniden(nspmax)
      common/irread/irreadlog
      common/outlog/shutdown,pereng,perfric,pervel,percell
      common/rdflog/rdfcall
      common/rstart/restart,nacl,cdcl2,mcl3
      common/presanneal/conjstart
      common/presrst/pressurerst,boxlenkeep
      common/runlogs/rim,pim,cim,aim,qim,quaim,debug,quaddamp,moveions,
     x               dynam,ddpot,relaxconfig
      common/thermlog/iontherm,pimtherm(nspmax),qimtherm
     x    ,cimtherm,aimtherm,
     x                quaimtherm
      common/strtup/rescalelog,velinit

c       new logicals for variable cell
      common/forceflag/forfl
      common/orthobaro/nth,ortho,nib,nab
      common/writepress/outp

c	logicals for catdip

      common/catdiplog/cdiplog,cdprev


c Characters in common.
      common/files/engout,engout2,velout,crdout,accnout,rstout,rstin,
     x             crstout,crstin,head,qout,zetout,fluxout,tempout,
     x             eleout,untout,cfout,cfin,prerstin,prerstout,fileout



      common/CG/tol,ftol,nmat1,nmat2

      common/cutoff/rcut


c   reviewed variables...
c Logicals
      common/runlog/conjrodlog,conjgradlog,dippimlog,quadpimlog
     x             ,cimlog,rimlog,aimlog,epplog,cim1log,cim2log
     x             ,msdcalllog,rvplog

      common/boxlog/boxlenfromdenslog
     x             ,boxlenfromrstlog
     x             ,boxlenfromcrdslog

      common/polarl/polarizablelog(nspmax)

      common/perlog/nrscalelog,nanneallog


c Variables for the conserved quantity for the variable cell.
      common/conser/tcell,tvol,tbzeta,tpzeta,PeeVee,pzeta,bzeta
     x             ,Pcomx,Pcomy,Pcomz

      common/displl/dxsav(10,10)
     x             ,dysav(10,10)
     x             ,dzsav(10,10),drsav(10,10)
      common/CUEppp/CUEp,CUEp2,CUEprec,CUEp2rec,dom
      common/CUEbbb/CUEb,CUEb2,CUEbrec,CUEb2rec
      common/WgoWWW/Wgo,W,Wrec,Wgorec

      common/srdipe/srdipeng(nspmax)
      common/srquda/srquadeng(nspmax)

c filenames added 18.4.2000
      common/filess/cellengout,engtotout,pzetaout,bzetaout
     x             ,fkeout,polengout,tdipout,diagpresout
     x             ,presout,celllenout,cellvolout,cellanglesout

      common/files2/diagstressout,xxyyzzstressout,xyxzyzstressout
     x             ,poscartout,cellboxout,polstressout
     x             ,coulsrstressout,momout

c T/p ramping variables
      common/ramplog/tramplog,pramplog,tramprescalelog

      common/rampva/nsteptramp,nsteppramp,deltatramp,deltapramp

c RVP variables
      common/rvprvp/rvph(nspmax,nspmax)
     x             ,rvpn(nspmax,nspmax)
     x             ,rvpr4(nspmax,nspmax)

c FG variables.
      common/fgfgfg/fgout
      common/fgfgf2/fgb(nspmax,nspmax)
     x             ,fgc(nspmax,nspmax)
     x             ,nkfg(nspmax,nspmax)
      common/fgfgf3/fglog


c AIM vairables

      common/aimaim/daimlog,quaimlog



      common/xftxft/xftlog
      common/xftxf3/ftlog
      common/xftxf2/ftalpx(nspmax,nspmax)                            
     x             ,ftbx(nspmax,nspmax)                              
     x             ,nrpower(nspmax,nspmax)          

c LJ variables.
      common/lj0001/epsilon(nspmax,nspmax),sigma(nspmax,nspmax)
      common/lj0002/fourepsil(nspmax,nspmax),sig12(nspmax,nspmax)
     x             ,sig6(nspmax,nspmax)
      common/lj0003/nochargelog,ljpotlog
      common/lj0004/sig24(nspmax,nspmax)
      common/lj0005/cut(nspmax,nspmax)

      logical nochargelog,ljpotlog

c SW variables
      common/sw0001/swepsilon(nspmax,nspmax),swsigma(nspmax,nspmax)
     x             ,swa(nspmax,nspmax),swp(nspmax,nspmax)
     x             ,swq(nspmax,nspmax),swbigA(nspmax,nspmax)
     x             ,swB(nspmax,nspmax),swgamma(nspmax,nspmax)
     x             ,swlambda(nspmax,nspmax)
      common/sw0002/sw2potlog,sw3potlog

      logical sw2potlog,sw3potlog

c Tersoff variables 
      common/cwcwcw/bigR, bigD, tf_lamda1, tf_lamda2, tf_lamda3
      common/tf2var/tf_beta, tf_n, tf_h, tf_c , tf_d
     x             , bigA, bigB 
      common/tf2tf2/tf2potlog
 
      logical tf2potlog
      common/cwcwc1/bigR_type(nspmax,nspmax), 
     x              bigD_type(nspmax,nspmax)
      common/cwcwc2/tf_lamda1_type(nspmax,nspmax),
     x                     tf_lamda2_type(nspmax,nspmax), 
     x                     tf_lamda3_type(nspmax,nspmax)
      common/cwcwc3/tf_beta_type(nspmax,nspmax),
     x               tf_n_type(nspmax,nspmax),
     x               tf_h_type(nspmax,nspmax),
     x               tf_c_type(nspmax,nspmax) ,
     x               tf_d_type(nspmax,nspmax)
     x             , bigA_type(nspmax,nspmax)
     x             , bigB_type(nspmax,nspmax)

c mixlog variables 
      common/mixmix/mixlog

      logical mixlog

c Morse variables
      common/morse1/amorse,remorse,demorse
      common/morse2/npairsmorse,idenmorse(nummax,2)
      common/morse3/morselog
      common/morse4/engmorse

      logical morselog
      logical fixedlog,fixedarray
      logical fixedzlog,fixedzarray
c fixed ions
      common/fixed/nionsfixed,nfixed(nummax)
     x            ,fixedlog,fixedarray(nummax)
     x , engTF2,eng_ion,engtot_lj
      common/fixed2/nionsfixedz,nfixedz(nummax)
     x            ,fixedzlog,fixedzarray(nummax)

c Harmonic potential variables
      common/harmp1/harmpotlog,harmpotlog2,harmpotlog3
      logical harmpotlog,harmpotlog2,harmpotlog3
      common/harmp2/r0harm(nspmax,nspmax),xkharm(nspmax,nspmax)
      common/harmp2a/xklin(nspmax,nspmax),u0lin(nspmax,nspmax)
      common/harmp3/npairs
      common/harmp4/npairident2(18000)
      common/harmp4/npairident1(18000)
