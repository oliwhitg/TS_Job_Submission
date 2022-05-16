      subroutine output
c
c************************************************************
c All (i.e non-temporary) diagnostic and generated output
c data should go to disc via this routine.
c************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'
c
c Periodic dump.
c
      polengtot=0.0d0

      if (pereng) then

         pereng=.false.

         if (dippimlog) then
            do 111 i=1,nspec
               polengtot=polengtot+reseng(i)
  111       continue
         endif

         if (quadpimlog) then
            do 112 i=1,nspec
               polengtot=polengtot+reseng(i)+dipsqeng(i)
     x                            +dipquadeng(i)+quadeng(i)
  112       continue
         endif

         fketot=fketot*0.5d0
c energies
         write(21,*)nstep,real(engpetot),real(tranke)
c         write(32,*)nstep,real(tcell),real(tvol)
c         write(34,*)nstep,real(tpzeta),real(pzeta)
c         write(44,*)nstep,real(tbzeta),real(bzeta)
         write(33,*)nstep,real(PeeVee)
     x        ,real(engpetot+tranke+fketot+polengtot
     x             +tcell+tvol+PeeVee+pzeta+tpzeta+tbzeta+bzeta)
c         write(40,*)nstep,(fke(i),i=1,nspec)
c         write(38,*)nstep,(reseng(i),i=1,nspec)
         write(86,*)nstep, real(eng_ion), real(engTF2), real(engtot_lj)  

         if (quadpimlog) then
            write(950,*)nstep,(fkequad(i),i=1,nspec)
            write(951,*)nstep,(dipquadeng(i),i=1,nspec)
            write(952,*)nstep,(dipsqeng(i),i=1,nspec)
            write(953,*)nstep,(quadeng(i),i=1,nspec)
         endif

         if (cimlog) then
        write(995,*)nstep,real(selfengtot),real(fkeself), real(engpetot)
     x                       ,real(tdelta)
            write(25,*)nstep,real(fketot),real(polengtot)
     x    ,real(engpetot+tranke+fketot+polengtot+selfengtot+fkeself)
         else if (daimlog) then
            write(995,*)nstep,real(selfengtot),real(fkeself)
     x                       ,real(tdelta)
            write(996,*)nstep,real(epsselfengtot),real(fkeeps)
     x                       ,real(teps)
            write(25,*)nstep,real(fketot),real(polengtot)
     x    ,real(engpetot+tranke+fketot+polengtot+selfengtot+fkeself
     x         +epsselfengtot+fkeeps)
         else if (quaimlog) then
            write(995,*)nstep,real(selfengtot),real(fkeself)
     x                       ,real(tdelta)
            write(996,*)nstep,real(epsselfengtot),real(fkeeps)
     x                       ,real(teps)
            write(997,*)nstep,real(quaimselfengtot),real(fkequaim)
     x                       ,real(tquaim)
            write(25,*)nstep,real(fketot),real(polengtot)
     x    ,real(engpetot+tranke+fketot+polengtot+selfengtot+fkeself
     x         +epsselfengtot+fkeeps+quaimselfengtot+fkequaim)
         else
            write(25,*)nstep,real(fketot),real(polengtot)
     x             ,real(engpetot+tranke+fketot+polengtot)
         endif

c cell momentum
c         write(36,*)nstep,real(Pcomx),real(Pcomy),real(Pcomz)

c stress tensor
c
c total term:
c
c         write (56,*) nstep,(stpxx+stcxx+stsrxx+stpyy+stcyy+
c     x                 stsryy+stpzz+stczz+stsrzz+stp2xx+
c     x                 stp2yy+stp2zz+stpsrxx+stpsryy+stpsrzz)*
c     x                 cellvol3rec
c
c components: (xx, yy, zz)
c
c         write (57,*) nstep,real((stpxx+stcxx+stsrxx+
c     x                 stp2xx+stpsrxx)*cellvol3rec),
c     x                 real((stpyy+stcyy+stsryy+stp2yy+stpsryy)*
c     x                 cellvol3rec),real((stpzz+stczz+stsrzz+
c     x                 stp2zz+stpsrzz)*cellvol3rec)
c
c components: (xy, xz, yz)
c
c         write (58,*) nstep,real((stpxy+stcxy+stsrxy+
c     x                 stp2xy+stpsrxy)*cellvol3rec),
c     x                 real((stpxz+stcxz+stsrxz+stp2xz+stpsrxz)*
c     x                 cellvol3rec),real((stpyz+stcyz+stsryz+
c     x                 stp2yz+stpsryz)*cellvol3rec)
c
c charge-dipole + s-r dipole + dipole-dipole term:
c
c         write (62,*) nstep,real((stpxx+stpyy+stpzz)*cellvol3rec),
c     x                 real((stpsrxx+stpsryy+stpsrzz)*cellvol3rec),
c     x                 real((stp2xx+stp2yy+stp2zz)*cellvol3rec)
c
c charge-charge + short-range term:
c
c         write (63,*) nstep,(stcxx+stcyy+stczz)*cellvol3rec,
c     x                 (stsrxx+stsryy+stsrzz)*cellvol3rec

c****         write(500,*)stsrxx,stsrxy,stsrxz
c****         write(500,*)stsryy,stsryz,stsrzz
c****         write(501,*)stcxx,stcxy,stcxz
c****         write(501,*)stcyy,stcyz,stczz
c****         write(502,*)stpxx,stpxy,stpxz
c****         write(502,*)stpyy,stpyz,stpzz
c****         write(503,*)stp2xx,stp2xy,stp2xz
c****         write(503,*)stp2yy,stp2yz,stp2zz
c****         write(504,*)stpsrxx,stpsrxy,stpsrxz
c****         write(504,*)stpsryy,stpsryz,stpsrzz


cc         write (173,*) nstep, (srdipeng(i),i=1,nspec)

cc         write(174,*)real(erfcracc),real(qdipacc),real(dipdipacc)
cc         write (121,*) real(qmueng),real(xmumueng),real(srdipengtot)


cc         write (280,*) sreng(1,1),sreng(2,1),sreng(2,2)
cc         write (281,*) ddeng(1,1),ddeng(2,1),ddeng(2,2)
cc         write (282,*) dqeng(1,1),dqeng(2,1),dqeng(2,2)
cc         write (283,*) ddeng(1,1)+ddeng(2,1)+ddeng(2,2)+
cc     x                 dqeng(1,1)+dqeng(2,1)+dqeng(2,2)

      endif

c=================================================================
c from trans_chains

      if (percell) then

         percell=.false.
         write(51,*) nstep,pint2(1,1),pint2(2,2),pint2(3,3)
         write(52,990) nstep,pinst3,((pint2(1,1)+pint2(2,2)+pint2(3,3))
     x                                    /3.0d0)

 990         format(I10, F4.1, E13.4)
         write(55,*) nstep,real(bee(4)),real(bee(5)),real(bee(6))
         write(53,*) nstep,real(boxlenx),real(boxleny)
     x                     ,real(boxlenz)
         write(54,930) nstep,vol3,(float(num)/vol3)
 930         format(I10, F12.2, E13.4)
      endif

c=================================================================
      if (pervel) then

         pervel=.false.

         write(61,*)hlab3(1,1),hlab3(1,2),hlab3(1,3)
         write(61,*)hlab3(2,1),hlab3(2,2),hlab3(2,3)
         write(61,*)hlab3(3,1),hlab3(3,2),hlab3(3,3)
         write(61,*)boxlenx,boxleny,boxlenz

         do 30 i=1,num
            write(22)vx(i),vy(i),vz(i)
            write(23)x(i),y(i),z(i)
            if (fglog) write(31,*)exx(i),eyy(i),ezz(i)
     x                           ,exy(i),exz(i),eyz(i)
     x                           ,exxsr(i),eyysr(i),ezzsr(i)
     x                           ,exysr(i),exzsr(i),eyzsr(i)

c       write out positions in lab. frame
            write(60,830)
     x                  hlab3(1,1)*x(i)+hlab3(1,2)*y(i)
     x                   +hlab3(1,3)*z(i),
     x                  hlab3(2,1)*x(i)+hlab3(2,2)*y(i)
     x                  +hlab3(2,3)*z(i),
     x                  hlab3(3,1)*x(i)+hlab3(3,2)*y(i)
     x                  +hlab3(3,3)*z(i)
 830         format (f10.5, f10.5, f10.5)
c            write (146,*) nstep,delta(i)
c            write (147,*) nstep,real(epsilonx(i)),real(epsilony(i)),
c     x                    real(epsilonz(i))
c            write (148,*) nstep,real(quaimxx(i)),real(quaimyy(i)),
c     x                    real(quaimzz(i))
c            write (149,*) nstep,real(quaimxy(i)),real(quaimxz(i)),
c     x                    real(quaimyz(i))
   30    continue

      endif

      if (perfric) then
 
         perfric=.false.

c         write(27,*)nstep,real(dtzeta2),real(dtdip2),real(dtquad2)
c         write(28,*)nstep,real(zeta),real(tfluxtran)
c         write(29,*)nstep,real(tfluxtran),real(tfluxdip),real(tfluxquad)
         write(30,*)nstep,real(tke)
c         write(39,*)nstep,(tdip(i),i=1,nspec)

 101     format(i7,f11.4,5f11.5)
c
c Write out the 3 dipole vector components and the five quadrupole
c tensor components. Note that these include the various symmetry
c relationships (i.e the delta_{alphagamma,betadelta} kind of
c things) as in Buckingham, Adv. Chem Phys. See theory book 3
c for the details.
c
         if (dippimlog) then
            nchan=77
            nnchan=901
            do 100 i=1,nummon
               nnn=3*(i-1)
               nt=ntype(nmon(i))
               write(nchan+nnn,*)nstep,real(xmu(nmon(i)))
     x          ,real(elecx(nmon(i))*alppolar(nt))
     x          ,real(elecx(nmon(i)))
               write(nchan+1+nnn,*)nstep,real(ymu(nmon(i)))
     x          ,real(elecy(nmon(i))*alppolar(nt))
     x          ,real(elecy(nmon(i)))
               write(nchan+2+nnn,*)nstep,real(zmu(nmon(i)))
     x          ,real(elecz(nmon(i))*alppolar(nt))
     x          ,real(elecz(nmon(i)))
c              write(nnchan+nnn,*)nstep,real(asdipx(nmon(i)))
c    x          ,real(srdipx(nmon(i)))
c              write(nnchan+1+nnn,*)nstep,real(asdipy(nmon(i)))
c    x          ,real(srdipy(nmon(i)))
c              write(nnchan+2+nnn,*)nstep,real(asdipz(nmon(i)))
c    x          ,real(srdipz(nmon(i)))
  100       continue

         endif

         if (quadpimlog) then
            nchan=80
            nnchan=801
            nnnchan=901
            nnnnchan=921
            do 105 i=1,nummon
               nnn=3*(i-1)
               nnnn=6*(i-1)
               nt=ntype(nmon(i))
               write(nchan+nnn,*)nstep,real(xmu(nmon(i)))
     x          ,real(elecx(nmon(i))*alppolar(nt))
     x          ,real(elecx(nmon(i))*alppolar(nt)
     x       +(1.0d0/3.0d0)*Bpolar(nt)*(elecx(nmon(i))*exx(nmon(i))
     x       +(1.5d0*elecy(nmon(i))*exy(nmon(i)))
     x       +(1.5d0*elecz(nmon(i))*exz(nmon(i)))))
               write(nchan+1+nnn,*)nstep,real(ymu(nmon(i)))
     x          ,real(elecy(nmon(i))*alppolar(nt))
     x        ,real(elecy(nmon(i))*alppolar(nt)
     x       +(1.0d0/3.0d0)*Bpolar(nt)*(elecy(nmon(i))*eyy(nmon(i))
     x       +(1.5d0*elecx(nmon(i))*exy(nmon(i)))
     x       +(1.5d0*elecz(nmon(i))*eyz(nmon(i)))))
               write(nchan+2+nnn,*)nstep,real(zmu(nmon(i)))
     x          ,real(elecz(nmon(i))*alppolar(nt))
     x        ,real(elecz(nmon(i))*alppolar(nt)
     x       +(1.0d0/3.0d0)*Bpolar(nt)*(elecz(nmon(i))*ezz(nmon(i))
     x       +(1.5d0*elecy(nmon(i))*eyz(nmon(i)))
     x       +(1.5d0*elecx(nmon(i))*exz(nmon(i)))))

               write(nnchan+nnnn,*)nstep,real(quadxx(nmon(i)))
     x        ,real((2.0d0/3.0d0)*Cpolar(nt)*exx(nmon(i))
     x        -(1.0d0/3.0d0)*Cpolar(nt)*eyy(nmon(i))
     x        -(1.0d0/3.0d0)*Cpolar(nt)*ezz(nmon(i)))
     x        ,real((2.0d0/3.0d0)*Cpolar(nt)*exx(nmon(i))
     x        -(1.0d0/3.0d0)*Cpolar(nt)*eyy(nmon(i))
     x        -(1.0d0/3.0d0)*Cpolar(nt)*ezz(nmon(i))
     x        +0.5d0*Bpolar(nt)*elecx(nmon(i))*elecx(nmon(i))
     x        -0.25d0*Bpolar(nt)*elecy(nmon(i))*elecy(nmon(i))
     x        -0.25d0*Bpolar(nt)*elecz(nmon(i))*elecz(nmon(i)))
               write(nnchan+1+nnnn,*)nstep,real(quadyy(nmon(i)))
     x        ,real((2.0d0/3.0d0)*Cpolar(nt)*eyy(nmon(i))
     x        -(1.0d0/3.0d0)*Cpolar(nt)*exx(nmon(i))
     x        -(1.0d0/3.0d0)*Cpolar(nt)*ezz(nmon(i)))
     x        ,real((2.0d0/3.0d0)*Cpolar(nt)*eyy(nmon(i))
     x        -(1.0d0/3.0d0)*Cpolar(nt)*exx(nmon(i))
     x        -(1.0d0/3.0d0)*Cpolar(nt)*ezz(nmon(i))
     x        +0.5d0*Bpolar(nt)*elecy(nmon(i))*elecy(nmon(i))
     x        -0.25d0*Bpolar(nt)*elecx(nmon(i))*elecx(nmon(i))
     x        -0.25d0*Bpolar(nt)*elecz(nmon(i))*elecz(nmon(i)))
               write(nnchan+2+nnnn,*)nstep,real(quadzz(nmon(i)))
     x        ,real((2.0d0/3.0d0)*Cpolar(nt)*ezz(nmon(i))
     x        -(1.0d0/3.0d0)*Cpolar(nt)*eyy(nmon(i))
     x        -(1.0d0/3.0d0)*Cpolar(nt)*exx(nmon(i)))
     x        ,real((2.0d0/3.0d0)*Cpolar(nt)*ezz(nmon(i))
     x        -(1.0d0/3.0d0)*Cpolar(nt)*eyy(nmon(i))
     x        -(1.0d0/3.0d0)*Cpolar(nt)*exx(nmon(i))
     x        +0.5d0*Bpolar(nt)*elecz(nmon(i))*elecz(nmon(i))
     x        -0.25d0*Bpolar(nt)*elecy(nmon(i))*elecy(nmon(i))
     x        -0.25d0*Bpolar(nt)*elecx(nmon(i))*elecx(nmon(i)))
               write(nnchan+3+nnnn,*)nstep,real(quadxy(nmon(i)))
     x        ,real(Cpolar(nt)*exy(nmon(i)))
     x        ,real((0.75d0*Bpolar(nt)*elecx(nmon(i))*elecy(nmon(i)))
     x        +(Cpolar(nt)*exy(nmon(i))))
               write(nnchan+4+nnnn,*)nstep,real(quadxz(nmon(i)))
     x        ,real(Cpolar(nt)*exz(nmon(i)))
     x        ,real((0.75d0*Bpolar(nt)*elecx(nmon(i))*elecz(nmon(i)))
     x        +(Cpolar(nt)*exz(nmon(i))))
               write(nnchan+5+nnnn,*)nstep,real(quadyz(nmon(i)))
     x        ,real(Cpolar(nt)*eyz(nmon(i)))
     x        ,real((0.75d0*Bpolar(nt)*elecy(nmon(i))*elecz(nmon(i)))
     x        +(Cpolar(nt)*eyz(nmon(i))))

               write(nnnchan+nnn,*)nstep,real(asdipx(nmon(i)))
     x          ,real(srdipx(nmon(i)))
               write(nnnchan+1+nnn,*)nstep,real(asdipy(nmon(i)))
     x          ,real(srdipy(nmon(i)))
               write(nnnchan+2+nnn,*)nstep,real(asdipz(nmon(i)))
     x          ,real(srdipz(nmon(i)))

               write(nnnnchan+nnnn,*)nstep,real(asquadxx(nmon(i)))
     x          ,real(srquadxx(nmon(i)))
               write(nnnnchan+1+nnnn,*)nstep,real(asquadyy(nmon(i)))
     x          ,real(srquadyy(nmon(i)))
               write(nnnnchan+2+nnnn,*)nstep,real(asquadzz(nmon(i)))
     x          ,real(srquadzz(nmon(i)))
               write(nnnnchan+3+nnnn,*)nstep,real(asquadxy(nmon(i)))
     x          ,real(srquadxy(nmon(i)))
               write(nnnnchan+4+nnnn,*)nstep,real(asquadxz(nmon(i)))
     x          ,real(srquadxz(nmon(i)))
               write(nnnnchan+5+nnnn,*)nstep,real(asquadyz(nmon(i)))
     x          ,real(srquadyz(nmon(i)))

  105       continue

         endif

         if ((cimlog).or.(daimlog).or.(quaimlog)) then

            nchan=701
            do 110 i=1,nummon
               nnn=(i-1)
               write(nchan+nnn,*)nstep,real(delta(nmon(i)))
  110       continue

         endif

         if ((daimlog).or.(quaimlog)) then

            nchan=711
            do 210 i=1,nummon
               nnn=(i-1)
               write (nchan+nnn,*)nstep,real(epsilonx(nmon(i)))
     x           ,real(epsilony(nmon(i))),real(epsilonz(nmon(i)))
  210       continue

         endif 

         if (quaimlog) then

            nchan=721
            do 310 i=1,nummon
               nnn=2*(i-1)
               write (nchan+nnn,*)nstep,real(quaimxx(nmon(i)))
     x                                ,real(quaimyy(nmon(i)))
     x                                ,real(quaimzz(nmon(i)))
               write (nchan+nnn+1,*)nstep,real(quaimxy(nmon(i)))
     x                                ,real(quaimxz(nmon(i)))
     x                                ,real(quaimyz(nmon(i)))

  310       continue

         endif

      endif

      return

      end
