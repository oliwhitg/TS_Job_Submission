      subroutine trans_vv

c************************************************************
c  Translates the ions via a velocity-Verlet algorithm.
c************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'

      double precision pint3c(3,3), vtx(nummax),vty(nummax),
     x vtz(nummax), vxni(nummax),vyni(nummax),vzni(nummax),
     x vxold(nummax),vyold(nummax),vzold(nummax),
     x tx(nummax),ty(nummax),tz(nummax)

      double precision gpzeta2(nchains),gpzeta3(nchains)
     x                ,gbzeta2(nchains),gbzeta3(nchains)

      real*8 riden(3,3)

c     update forces
      if (forfl) then
         do 10 i=1,num
            frrx(i)=frrx3(i)
            frry(i)=frry3(i)
            frrz(i)=frrz3(i)
 10      continue
      endif

      forfl=.true.

C.............thermostat for particles
      do ir=1,5
        vpzeta1(ir)=vpzeta2(ir)
        vpzeta2(ir)=vpzeta3(ir)
        pzeta2(ir)=pzeta3(ir)
        gpzeta2(ir)=0.0d0
        gpzeta3(ir)=0.0d0
        if (.not.nth) then
           vpzeta1(ir)=0.0d0
           vpzeta2(ir)=0.0d0
           vpzeta3(ir)=0.0d0
           pzeta2(ir)=0.0d0
         endif
      enddo

C.............thermostat for barostat
      do ir=1,5
        vbzeta1(ir)=vbzeta2(ir)
        vbzeta2(ir)=vbzeta3(ir)
        bzeta2(ir)=bzeta3(ir)
        gbzeta2(ir)=0.0d0
        gbzeta3(ir)=0.0d0
        if ((.not.nth).or.(.not.nib)) then
           vbzeta1(ir)=0.0d0
           vbzeta2(ir)=0.0d0
           vbzeta3(ir)=0.0d0
           bzeta2(ir)=0.0d0
        endif
      enddo
c............volume variables
      veps1=veps2
      veps2=veps3
      eps2=eps3
   
      vol2=vol3

      if (.not.nib) then
         veps1=0.0d0
         veps2=0.0d0
         veps3=0.0d0
      endif

c.............cell shape variables
      do 20 j=1,3
         do 30 k=1,3
            vg2(j,k)=vg3(j,k)
            riden(j,k)=0.0
            iden(j,k)=0
 30      continue
         riden(j,j)=1.0d0
         iden(j,j)=1
 20   continue

      if (.not.nab) then
         do 21 j=1,3
            do 31 k=1,3
               vg2(j,k)=0.0d0
               vg3(j,k)=0.0d0
   31       continue
   21    continue
      endif

      do 40 j=1,3
         do 51 k=1,3
            h2(j,k)=h3(j,k)
            hi2(j,k)=hi3(j,k)
            hlab2(j,k)=hlab3(j,k)
 51      continue
 40   continue
  
c 'xkesq2' is the k.e. at time t.
      xkesq2=0.0d0
      do 60 i=1,num

         ionpoint=ntype(i)

c vtx,y,z are the velocities in the Cartesian frame.
         vtx(i)=hlab2(1,1)*vx(i)+hlab2(1,2)*vy(i)+hlab2(1,3)*vz(i)
         vty(i)=hlab2(2,1)*vx(i)+hlab2(2,2)*vy(i)+hlab2(2,3)*vz(i)
         vtz(i)=hlab2(3,1)*vx(i)+hlab2(3,2)*vy(i)+hlab2(3,3)*vz(i)

         xkesq2=xkesq2+(amass(ionpoint)*(vtx(i)*vtx(i)+
     x               vty(i)*vty(i)+vtz(i)*vtz(i)))         

c positions in the Cartesian frame.
         tx(i)=hlab2(1,1)*x(i)+hlab2(1,2)*y(i)+hlab2(1,3)*z(i)
         ty(i)=hlab2(2,1)*x(i)+hlab2(2,2)*y(i)+hlab2(2,3)*z(i)
         tz(i)=hlab2(3,1)*x(i)+hlab2(3,2)*y(i)+hlab2(3,3)*z(i)

 60   continue

c     Obtain pressure at time t=0

      do 70 j=1,3
         do 80 k=1,3
            pint2(j,k)=0.0d0
 80      continue
 70   continue

      do 84 i=1,num
         ionpoint=ntype(i)
         rmassvol=amass(ionpoint)/cellvol
         pint2(1,1)=pint2(1,1)+ rmassvol*vtx(i)*vtx(i)
         pint2(1,2)=pint2(1,2)+ rmassvol*vtx(i)*vty(i)
         pint2(1,3)=pint2(1,3)+ rmassvol*vtx(i)*vtz(i)
         pint2(2,2)=pint2(2,2)+ rmassvol*vty(i)*vty(i)
         pint2(2,3)=pint2(2,3)+ rmassvol*vty(i)*vtz(i)
         pint2(3,3)=pint2(3,3)+ rmassvol*vtz(i)*vtz(i)
   84 continue

      pint2(1,1)=pint2(1,1)+(((3.0d0*(stpxx+
     x                 stcxx+stsrxx+
     x                 stp2xx+stpsrxx))*cellvol3rec))

      pint2(2,2)=pint2(2,2)+(((3.0d0*(stpyy+
     x                 stcyy+stsryy+
     x                 stp2yy+stpsryy))*cellvol3rec))

      pint2(3,3)=pint2(3,3)+(((3.0d0*(stpzz+
     x                 stczz+stsrzz+
     x                 stp2zz+stpsrzz))*cellvol3rec))

      pint2(1,2)=pint2(1,2)+(((3.0d0*(stpxy+
     x                 stcxy+stsrxy+
     x                 stp2xy+stpsrxy))*cellvol3rec))

      pint2(1,3)=pint2(1,3)+(((3.0d0*(stpxz+
     x                 stcxz+stsrxz+
     x                 stp2xz+stpsrxz))*cellvol3rec))

      pint2(2,3)=pint2(2,3)+(((3.0d0*(stpyz+
     x                 stcyz+stsryz+
     x                 stp2yz+stpsryz))*cellvol3rec))

      pint2(2,1)=pint2(1,2)
      pint2(3,1)=pint2(1,3)
      pint2(3,2)=pint2(2,3)

      pinst2=(pint2(1,1)+pint2(2,2)+pint2(3,3))/3.0d0

      if (ortho) then
         do 85 j=1,3
            do 86 k=1,3

	       if (j.ne.k) then
                  pint2(j,k)=0.0d0
   	       endif

 86         continue
 85      continue

      endif

      trvg2=0.0d0

c Calculate trace of vg**2
      do 90 j=1,3
         do 100 k=1,3
            trvg2=trvg2+(vg2(j,k)*vg2(k,j))
 100     continue
 90   continue

c.......update particle zetas
      gpzeta2(1)=(CUEprec)*((xkesq2) 
     x      -(free*trantkb))-vpzeta2(1)*vpzeta2(2)
      gpzeta2(2)=(CUEp2rec)*(CUEp*vpzeta2(1)**2 
     x      - trantkb)-vpzeta2(2)*vpzeta2(3)
      gpzeta2(3)=(CUEp2rec)*(CUEp2*vpzeta2(2)**2 
     x      - trantkb)-vpzeta2(3)*vpzeta2(4)
      gpzeta2(4)=(CUEp2rec)*(CUEp2*vpzeta2(3)**2 
     x      - trantkb)-vpzeta2(4)*vpzeta2(5)
      gpzeta2(5)=(CUEp2rec)*(CUEp2*vpzeta2(4)**2 
     x      - trantkb) 

      pzeta3(1)=pzeta2(1)+vpzeta2(1)*dtime+gpzeta2(1)*(dtime**2)/2.0d0
      pzeta3(2)=pzeta2(2)+vpzeta2(2)*dtime+gpzeta2(2)*(dtime**2)/2.0d0
      pzeta3(3)=pzeta2(3)+vpzeta2(3)*dtime+gpzeta2(3)*(dtime**2)/2.0d0
      pzeta3(4)=pzeta2(4)+vpzeta2(4)*dtime+gpzeta2(4)*(dtime**2)/2.0d0
      pzeta3(5)=pzeta2(5)+vpzeta2(5)*dtime+gpzeta2(5)*(dtime**2)/2.0d0

c.......update barostat zetas
      gbzeta2(1)=(CUEbrec)*( W*veps2**2
     x      +Wgo*trvg2  -dom*trantkb)-vbzeta2(1)*vbzeta2(2)
      gbzeta2(2)=(CUEb2rec)*(CUEb*vbzeta2(1)**2 
     x      - trantkb)-vbzeta2(2)*vbzeta2(3)
      gbzeta2(3)=(CUEb2rec)*(CUEb2*vbzeta2(2)**2 
     x      - trantkb)-vbzeta2(3)*vbzeta2(4)
      gbzeta2(4)=(CUEb2rec)*(CUEb2*vbzeta2(3)**2 
     x      - trantkb)-vbzeta2(4)*vbzeta2(5)
      gbzeta2(5)=(CUEb2rec)*(CUEb2*vbzeta2(4)**2 
     x      - trantkb) 

      bzeta3(1)=bzeta2(1)+vbzeta2(1)*dtime+gbzeta2(1)*(dtime**2)/2.0d0
      bzeta3(2)=bzeta2(2)+vbzeta2(2)*dtime+gbzeta2(2)*(dtime**2)/2.0d0
      bzeta3(3)=bzeta2(3)+vbzeta2(3)*dtime+gbzeta2(3)*(dtime**2)/2.0d0
      bzeta3(4)=bzeta2(4)+vbzeta2(4)*dtime+gbzeta2(4)*(dtime**2)/2.0d0
      bzeta3(5)=bzeta2(5)+vbzeta2(5)*dtime+gbzeta2(5)*(dtime**2)/2.0d0
 
c...........update the volume
      Feps2=((3.0d0)*vol2*(pinst2-pext))+
     x       ((3.0d0/free)*xkesq2)
   
      eps3=eps2+(dtime*veps2)+(((dtime**2)/2.0d0)*
     x       ((Feps2*Wrec)-(veps2*vbzeta2(1))))

      volfac=exp(eps3-eps2)
      vol3=vol2*volfac**3
 
c.......invert h0(t=0)

      call matinv(h2,hi2,detmh)

c........update the cell matrix

      do 110 j=1,3
         do 120 k=1,3

            Fgo2(j,k)=vol2*(pint2(j,k)-(pext*riden(j,k)))
     x        -((vol2/3.0d0)*(((pint2(1,1)-(riden(1,1)*pext))+
     x        (pint2(2,2)-(riden(2,2)*pext))+(pint2(3,3)-
     x        (riden(3,3)*pext)))*riden(j,k)))

 120     continue
 110  continue

      do 170 j=1,3
         do 180 k=1,3
            vg2sq(j,k)=0.0d0
            do 190 l=1,3
               vg2sq(j,k)=vg2sq(j,k)+(vg2(j,l)*vg2(l,k))
 190        continue

            hcom(j,k)=riden(j,k) + (vg2(j,k)*dtime) + (((dtime**2
     x        /2.0d0)*(Fgo2(j,k)*Wgorec)+(vg2sq(j,k))-(vg2(j,k)
     x        *vbzeta2(1))))

 180     continue
 170  continue

c............update h3, which is Martyna's h_0(Delta_t)
      do 200 j=1,3
         do 210 k=1,3
            h3(j,k)=0.0d0
            do 220 l=1,3

               h3(j,k)=h3(j,k)+(hcom(j,l)*h2(l,k))

c     Constrain cell to be orthorhombic

               if (ortho) then

                  if (j.ne.k) then
                     h3(j,k)=0.0d0
                  endif

                endif

 220        continue
 
 210     continue
 200  continue

c Impose constraint that |h(t+dt)| is equal to 1
      call shake

c.......invert h(t+dt)

      call matinv(h3,hi3,detmh)

      norvol3=detmh

c............hh corresponds to the matrix of cell-side vectors,
c............it is also known as fullh elsewhere in the program
      do 270 j=1,3
         do 280 k=1,3
            hh(j,k)=(vol3**(1.0d0/3.0d0))*h3(j,k)
 280     continue
 270  continue

      call dcell(hh,bee)

      bee(4)=acos(bee(4))
      bee(5)=acos(bee(5))
      bee(6)=acos(bee(6))

      boxlenx=bee(1)
      boxleny=bee(2)
      boxlenz=bee(3)

c..........hlab is the matrix which transforms from  cell to cartesian
c.......i.e. it is the matrix of unit vectors along the cell sides
c........it is the same as the matrix h which appears elsewhere in the program
      hlab3(1,1)=hh(1,1)/boxlenx
      hlab3(2,1)=hh(2,1)/boxlenx
      hlab3(3,1)=hh(3,1)/boxlenx
      hlab3(1,2)=hh(1,2)/boxleny
      hlab3(2,2)=hh(2,2)/boxleny
      hlab3(3,2)=hh(3,2)/boxleny
      hlab3(1,3)=hh(1,3)/boxlenz
      hlab3(2,3)=hh(2,3)/boxlenz
      hlab3(3,3)=hh(3,3)/boxlenz

      h(1,1)=hlab3(1,1)
      h(2,1)=hlab3(2,1)
      h(3,1)=hlab3(3,1)
      h(1,2)=hlab3(1,2)
      h(2,2)=hlab3(2,2)
      h(3,2)=hlab3(3,2)
      h(1,3)=hlab3(1,3)
      h(2,3)=hlab3(2,3)
      h(3,3)=hlab3(3,3)

c Update number density
c boxreset redetermines all the cell-dependent parameters....
c....AND ALSO INVERTS THE CELL MATRIX h -> hi as WELL AS fullh -> fullhi
c...i.e. the above is a highly circuitous route (but does some useful things!)
      call boxreset
      call kset

      do 320 l=1,3
         do 330 m=1,3

            comexp(l,m)=volfac*
     x        ((h3(l,1)*hi2(1,m))+(h3(l,2)*hi2(2,m))+
     x        (h3(l,3)*hi2(3,m)))
                  
            if (.not.nib) then
               comexp(l,m)=0.0d0
            endif
 330     continue

         if (.not.nib) then
            comexp(l,l)=1.0d0
         endif

 320  continue

      do 310 i=1,num

         ionpoint=ntype(i)

         ax(i)=frrx(i)*recamass(ionpoint)
         ay(i)=frry(i)*recamass(ionpoint)
         az(i)=frrz(i)*recamass(ionpoint)

c     calculate positions at t+dt

         rxn(i)=tx(i)+(dtime*vtx(i))+((0.5d0*dtime**2)*
     x   (ax(i)-(vpzeta2(1)*vtx(i))
     x   -(2.0d0*((vg2(1,1)*vtx(i))+(vg2(1,2)*vty(i))+
     x   (vg2(1,3)*vtz(i))))
     x   -((2.0d0+(3.0d0/free))*
     x   vtx(i)*veps2)))

         ryn(i)=ty(i)+(dtime*vty(i))+((0.5d0*dtime**2)*
     x   (ay(i)-(vpzeta2(1)*vty(i))
     x   -(2.0d0*((vg2(2,1)*vtx(i))+(vg2(2,2)*vty(i))+
     x   (vg2(2,3)*vtz(i))))
     x   -((2.0d0+(3.0d0/free))*
     x   vty(i)*veps2)))

         rzn(i)=tz(i)+(dtime*vtz(i))+((0.5d0*dtime**2)*
     x   (az(i)-(vpzeta2(1)*vtz(i))
     x   -(2.0d0*((vg2(3,1)*vtx(i))+(vg2(3,2)*vty(i))+
     x   (vg2(3,3)*vtz(i))))
     x   -((2.0d0+(3.0d0/free))*
     x   vtz(i)*veps2)))

         rxnn=rxn(i)
         rynn=ryn(i)
         rznn=rzn(i)

         rxn(i)=(comexp(1,1)*rxnn)+(comexp(1,2)*rynn)
     x  +(comexp(1,3)*rznn)
         ryn(i)=(comexp(2,1)*rxnn)+(comexp(2,2)*rynn)
     x  +(comexp(2,3)*rznn)
         rzn(i)=(comexp(3,1)*rxnn)+(comexp(3,2)*rynn)
     x  +(comexp(3,3)*rznn)

         rxnn=rxn(i)
         rynn=ryn(i)
         rznn=rzn(i)

c displacements for 'disp.out'. In Cartesiam frame.
         xdisp(i)=xdisp(i)+(rxn(i)-tx(i))
         ydisp(i)=ydisp(i)+(ryn(i)-ty(i))
         zdisp(i)=zdisp(i)+(rzn(i)-tz(i))

c.........rotate new ionic positions back to the cell frame ready for
c..... the update step
         rxn(i)=rxnn*hi(1,1)+rynn*hi(1,2)+rznn*hi(1,3)
         ryn(i)=rxnn*hi(2,1)+rynn*hi(2,2)+rznn*hi(2,3)
         rzn(i)=rxnn*hi(3,1)+rynn*hi(3,2)+rznn*hi(3,3)

 310  continue

c     forces at time t+dt are required 

      call dofmoveselect

      flagup=0

      call update_positions

      call ener

      do 350 i=1,num

         frrx3(i)=frrx(i)
         frry3(i)=frry(i)
         frrz3(i)=frrz(i)

 350  continue

      do 360 i=1,num

         ionpoint=ntype(i)

         ax3(i)=frrx3(i)*recamass(ionpoint)
         ay3(i)=frry3(i)*recamass(ionpoint)
         az3(i)=frrz3(i)*recamass(ionpoint)

 360  continue

c       Initial guess for vzeta3 and veps3 (from Verlet)

      vpzeta3(1)=vpzeta1(1)+(2.0d0*gpzeta2(1)*dtime)
      vpzeta3(2)=vpzeta1(2)+(2.0d0*gpzeta2(2)*dtime)
      vpzeta3(3)=vpzeta1(3)+(2.0d0*gpzeta2(3)*dtime)
      vpzeta3(4)=vpzeta1(4)+(2.0d0*gpzeta2(4)*dtime)
      vpzeta3(5)=vpzeta1(5)+(2.0d0*gpzeta2(5)*dtime)

      vbzeta3(1)=vbzeta1(1)+(2.0d0*gbzeta2(1)*dtime)
      vbzeta3(2)=vbzeta1(2)+(2.0d0*gbzeta2(2)*dtime)
      vbzeta3(3)=vbzeta1(3)+(2.0d0*gbzeta2(3)*dtime)
      vbzeta3(4)=vbzeta1(4)+(2.0d0*gbzeta2(4)*dtime)
      vbzeta3(5)=vbzeta1(5)+(2.0d0*gbzeta2(5)*dtime)

      veps3=veps1+(2.0d0*dtime)*((Feps2*Wrec)-(veps2*vbzeta2(1)))

c....update the configurational part of the stress tensor
 
      pint3c(1,1)=((3.0d0*(stpxx+
     x                 stcxx+stsrxx+
     x                 stp2xx+stpsrxx))*cellvol3rec) 

      pint3c(2,2)=((3.0d0*(stpyy+
     x                 stcyy+stsryy+
     x                 stp2yy+stpsryy))*cellvol3rec)

      pint3c(3,3)=((3.0d0*(stpzz+
     x                 stczz+stsrzz+
     x                 stp2zz+stpsrzz))*cellvol3rec)

      pint3c(1,2)=((3.0d0*(stpxy+
     x                 stcxy+stsrxy+
     x                 stp2xy+stpsrxy))*cellvol3rec)

      pint3c(1,3)=((3.0d0*(stpxz+
     x                 stcxz+stsrxz+
     x                 stp2xz+stpsrxz))*cellvol3rec)

      pint3c(2,3)=((3.0d0*(stpyz+
     x                 stcyz+stsryz+
     x                 stp2yz+stpsryz))*cellvol3rec)

      pint3c(2,1)=pint3c(1,2)
      pint3c(3,1)=pint3c(1,3)
      pint3c(3,2)=pint3c(2,3)
 
      xkesq=0.0d0
      do 530 i=1,num

         ionpoint=ntype(i)

c..........parts of velocity updates at dt=0 -> vxni, vyni etc

         comx=vtx(i)+((0.5d0*dtime)*(ax(i)-(vtx(i)*vpzeta2(1))
     x       -(2.0d0*((vg2(1,1)*vtx(i))+(vg2(1,2)*vty(i))
     x       +(vg2(1,3)*vtz(i))))
     x       -((2.0d0+(3.0d0/free))*veps2*vtx(i))))
         comy=vty(i)+((0.5d0*dtime)*(ay(i)-(vty(i)*vpzeta2(1))
     x       -(2.0d0*((vg2(2,1)*vtx(i))+(vg2(2,2)*vty(i))
     x       +(vg2(2,3)*vtz(i))))
     x       -((2.0d0+(3.0d0/free))*veps2*vty(i))))
         comz=vtz(i)+((0.5d0*dtime)*(az(i)-(vtz(i)*vpzeta2(1))
     x       -(2.0d0*((vg2(3,1)*vtx(i))+(vg2(3,2)*vty(i))
     x       +(vg2(3,3)*vtz(i))))
     x       -((2.0d0+(3.0d0/free))*veps2*vtz(i))))

c.........initial guess at the remaining parts of velocities --
c.........uses vg and v from t=0 rather than delta t -> vxn etc

         comx2=(0.5d0*dtime)*((ax3(i))-(vtx(i)*vpzeta3(1))
     x       -(2.0d0*((vg2(1,1)*vtx(i))+(vg2(1,2)*vty(i))
     x       +(vg2(1,3)*vtz(i))))
     x       -((2.0d0+(3.0d0/free))*vtx(i)*veps3))

         comy2=(0.5d0*dtime)*((ay3(i))-(vtx(i)*vpzeta3(1))
     x       -(2.0d0*((vg2(2,1)*vtx(i))+(vg2(2,2)*vty(i))
     x       +(vg2(2,3)*vtz(i))))
     x       -((2.0d0+(3.0d0/free))*vty(i)*veps3))
 
         comz2=(0.5d0*dtime)*((az3(i))-(vtz(i)*vpzeta3(1))
     x       -(2.0d0*((vg2(3,1)*vtx(i))+(vg2(3,2)*vty(i))
     x       +(vg2(3,3)*vtz(i))))
     x       -((2.0d0+(3.0d0/free))*vtz(i)*veps3))

         vxni(i)=((comexp(1,1)*comx)+(comexp(1,2)*comy)+
     x        (comexp(1,3)*comz))
              vxn(i)=vxni(i)+comx2

         vyni(i)=((comexp(2,1)*comx)+(comexp(2,2)*comy)+
     x        (comexp(2,3)*comz))
              vyn(i)=vyni(i)+comy2

         vzni(i)=((comexp(3,1)*comx)+(comexp(3,2)*comy)+
     x        (comexp(3,3)*comz))
              vzn(i)=vzni(i)+comz2

c
c Calculate kinetic energy:
c
 
         xkesq=xkesq+amass(ionpoint)*(vxn(i)*vxn(i)+
     x          vyn(i)*vyn(i)+vzn(i)*vzn(i))


 530  continue

      do 414 j=1,3
         do 424 k=1,3
            hhi(j,k)=0.0d0
            do 434 l=1,3
               hhi(j,k)=hhi(j,k)+(h2(j,l)*hi3(l,k))
               com1(j,k)=(Fgo2(j,k)*Wgorec)+vg2sq(j,k)
     x                       -(vbzeta2(1)*vg2(j,k))
               com2(j,k)=vg2(j,k)+(dtime*com1(j,k)/2.0d0)
 434        continue
 424     continue
 414  continue
c.......obtain all the contributions to vg3 up to t=0
 
      do 48 j=1,3
         do 49 k=1,3
            com3(j,k)=0.0d0
            do 50 l=1,3
               com3(j,k)=com3(j,k)+(com2(j,l)*hhi(l,k))
 50         continue
 49      continue
 48   continue

c.......make an initial guess at vg3
      trace=0.0d0
      do 42 j=1,3
         do 41 k=1,3
            vg3(j,k)=com3(j,k)+(dtime*com1(j,k)/2.0d0)

c     Constrain cell to be orthorhombic
            if (ortho) then           
               if (j.ne.k) then
                  vg3(j,k)=0.0d0
               endif
            endif

 41      continue
         trace=trace+vg3(j,j)     
 42   continue
      do 43 j=1,3
         vg3(j,j)=vg3(j,j)-(1.0d0/3.0d0)*trace
 43   continue

      if (.not.nib) then
         do 501 j=1,3
            do 502 k=1,3
               vg3(j,k)=0.0d0
 502        continue
 501     continue

         veps2=0.0d0
         veps3=0.0d0

      endif

c======================================================
c     begin iteration loop for velocities
      icnt=0
  380 continue
      icnt=icnt+1
      do 381 i=1,num
         vxold(i)=vxn(i)
         vyold(i)=vyn(i)
         vzold(i)=vzn(i)
  381 continue 

c     obtain pressure at time t+dt
 
      do 770 j=1,3
         do 880 k=1,3
            pint3(j,k)=0.0d0
 880     continue
 770  continue

      do 451 i=1,num
         ionpoint=ntype(i)
         rmassvol=amass(ionpoint)/cellvol
         pint3(1,1)=pint3(1,1)+ rmassvol*vxn(i)*vxn(i)
         pint3(1,2)=pint3(1,2)+ rmassvol*vxn(i)*vyn(i)
         pint3(1,3)=pint3(1,3)+ rmassvol*vxn(i)*vzn(i)
         pint3(2,2)=pint3(2,2)+ rmassvol*vyn(i)*vyn(i)
         pint3(2,3)=pint3(2,3)+ rmassvol*vyn(i)*vzn(i)
         pint3(3,3)=pint3(3,3)+ rmassvol*vzn(i)*vzn(i)
 451  continue
  
      do 771 j=1,3
         do 881 k=1,3
            pint3(j,k)=pint3(j,k)+pint3c(j,k)
 881     continue
 771  continue

      pint3(2,1)=pint3(1,2)
      pint3(3,1)=pint3(1,3)
      pint3(3,2)=pint3(2,3)
 
      if (ortho) then
         do 455 j=1,3
            do 456 k=1,3
               if (j.ne.k) then
                  pint3(j,k)=0.0d0
	       endif
 456        continue
 455     continue
      endif

      pinst3=0.0d0
      do k=1,3
         pinst3=pinst3+pint3(k,k)
      enddo
      pinst3=pinst3/3.0d0

c......get new force on cell shape

      do 410 j=1,3
         do 420 k=1,3
            vg3sq(j,k)=0.0d0
            do 430 l=1,3
               vg3sq(j,k)=vg3sq(j,k)+(vg3(j,l)*vg3(l,k))
 430        continue
 420     continue
 410  continue

      do 460 j=1,3
         do 470 k=1,3
                  
            Fgo3(j,k)=vol3*(pint3(j,k)-(pext*riden(j,k)))
     x        -((vol3/3.0d0)*(((pint3(1,1)-(riden(1,1)*pext))
     x                        +(pint3(2,2)-(riden(2,2)*pext))
     x                        +(pint3(3,3)-(riden(3,3)*pext)))
     x                        *riden(j,k)))

 470     continue
 460  continue

c.......update cell shape velocity 
      trace=0.0d0
      do 481 j=1,3
         do 491 k=1,3
 
            vg3(j,k)=com3(j,k)+(Fgo3(j,k)*Wgorec+
     x          vg3sq(j,k)-vg3(j,k)*vbzeta3(1))*(dtime/2.0d0)
                 

c     Constrain cell to be orthorhombic

            if (ortho) then
             
               if (j.ne.k) then
                  vg3(j,k)=0.0d0
               endif

            endif

 491     continue
         trace=trace+vg3(j,j)     
 481  continue
      do 482 j=1,3
         vg3(j,j)=vg3(j,j)-(1.0d0/3.0d0)*trace
 482  continue

      if (.not.nib) then
         do 510 j=1,3
            do 520 k=1,3
               vg3(j,k)=0.0d0
 520        continue
 510     continue

         veps2=0.0d0
         veps3=0.0d0

      endif

c...........now update the velocities
      xkesq=0.0d0

      Pcomx=0.0d0
      Pcomy=0.0d0
      Pcomz=0.0d0
      do 531 i=1,num

         ionpoint=ntype(i)

         comx2=(0.5d0*dtime)*((ax3(i))-(vxn(i)*vpzeta3(1))
     x       -(2.0d0*((vg3(1,1)*vxn(i))+(vg3(1,2)*vyn(i))
     x       +(vg3(1,3)*vzn(i))))
     x       -((2.0d0+(3.0d0/free))*vxn(i)*veps3))
         
         comy2=(0.5d0*dtime)*((ay3(i))-(vyn(i)*vpzeta3(1))
     x       -(2.0d0*((vg3(2,1)*vxn(i))+(vg3(2,2)*vyn(i))
     x       +(vg3(2,3)*vzn(i))))
     x       -((2.0d0+(3.0d0/free))*vyn(i)*veps3))

         comz2=(0.5d0*dtime)*((az3(i))-(vzn(i)*vpzeta3(1))
     x       -(2.0d0*((vg3(3,1)*vxn(i))+(vg3(3,2)*vyn(i))
     x       +(vg3(3,3)*vzn(i))))
     x       -((2.0d0+(3.0d0/free))*vzn(i)*veps3))

         vxn(i)=vxni(i)+comx2

         vyn(i)=vyni(i)+comy2

         vzn(i)=vzni(i)+comz2
 
         xkesq=xkesq+amass(ionpoint)*(vxn(i)*vxn(i)+
     x               vyn(i)*vyn(i)+vzn(i)*vzn(i))
       
         Pcomx=Pcomx+vxn(i)*amass(ionpoint)
         Pcomy=Pcomy+vyn(i)*amass(ionpoint)
         Pcomz=Pcomz+vzn(i)*amass(ionpoint)

 531  continue

c     Update the friction coefficient data

      trvg3=0.0d0
      do 600 j=1,3
         do 610 k=1,3
            vg3sq(j,k)=0.0d0
            do 620 l=1,3
               vg3sq(j,k)=vg3sq(j,k)+(vg3(j,l)*vg3(l,k))
 620        continue
 610     continue
         trvg3=trvg3+vg3sq(j,j)
 600  continue

c............update particle  zeta
 
      gpzeta3(1)=(CUEprec)*((xkesq) 
     x      -(free*trantkb))-vpzeta3(1)*vpzeta3(2)
      vpzeta3(1)=vpzeta2(1)+0.5d0*dtime*(gpzeta2(1)+gpzeta3(1))
      gpzeta3(2)=(CUEp2rec)*(CUEp*vpzeta3(1)**2 
     x      - trantkb)-vpzeta3(2)*vpzeta3(3)
      vpzeta3(2)=vpzeta2(2)+0.5d0*dtime*(gpzeta2(2)+gpzeta3(2))
      gpzeta3(3)=(CUEp2rec)*(CUEp2*vpzeta3(2)**2 
     x      - trantkb)-vpzeta3(3)*vpzeta3(4)
      vpzeta3(3)=vpzeta2(3)+0.5d0*dtime*(gpzeta2(3)+gpzeta3(3))
      gpzeta3(4)=(CUEp2rec)*(CUEp2*vpzeta3(3)**2 
     x      - trantkb)-vpzeta3(4)*vpzeta3(5)
      vpzeta3(4)=vpzeta2(4)+0.5d0*dtime*(gpzeta2(4)+gpzeta3(4))
      gpzeta3(5)=(CUEp2rec)*(CUEp2*vpzeta3(4)**2 
     x      - trantkb) 
      vpzeta3(5)=vpzeta2(5)+0.5d0*dtime*(gpzeta2(5)+gpzeta3(5))

c............  update barostat zeta 
 
      gbzeta3(1)=(CUEbrec)*(W*veps3**2
     x      +Wgo*trvg3-dom*trantkb)-vbzeta3(1)*vbzeta3(2)
      vbzeta3(1)=vbzeta2(1)+0.5d0*dtime*(gbzeta2(1)+gbzeta3(1))
      gbzeta3(2)=(CUEb2rec)*(CUEb*vbzeta3(1)**2 
     x      - trantkb)-vbzeta3(2)*vbzeta3(3)
      vbzeta3(2)=vbzeta2(2)+0.5d0*dtime*(gbzeta2(2)+gbzeta3(2))
      gbzeta3(3)=(CUEb2rec)*(CUEb2*vbzeta3(2)**2 
     x      - trantkb)-vbzeta3(3)*vbzeta3(4)
      vbzeta3(3)=vbzeta2(3)+0.5d0*dtime*(gbzeta2(3)+gbzeta3(3))
      gbzeta3(4)=(CUEb2rec)*(CUEb2*vbzeta3(3)**2 
     x      - trantkb)-vbzeta3(4)*vbzeta3(5)
      vbzeta3(4)=vbzeta2(4)+0.5d0*dtime*(gbzeta2(4)+gbzeta3(4))
      gbzeta3(5)=(CUEb2rec)*(CUEb2*vbzeta3(4)**2 
     x      - trantkb) 
      vbzeta3(5)=vbzeta2(5)+0.5d0*dtime*(gbzeta2(5)+gbzeta3(5))


c.............update epsilon
      Feps3=3.0d0*vol3*(pinst3-pext)+((3.0d0/free)*xkesq)

      veps3=veps2+(dtime/2.0d0)*(Feps2*Wrec-
     x        veps2*vbzeta2(1))+(dtime/2.0d0)*(Feps3*Wrec
     x        -veps3*vbzeta3(1))

 

c.......this is a pretty crude convergence indicator
c.......might be better to test the individual velocities
      converge=0.0d0 
      do i=1,num
         ion=ntype(i)
         converge=converge+amass(ion)*((vxn(i)-vxold(i))**2
     x           +(vyn(i)-vyold(i))**2+(vzn(i)-vzold(i))**2)
      enddo
      if (converge.gt.(free*trantkb*1.0d-17)) then
         write(6,*)icnt,converge,free*trantkb*1.0d-17
         if (icnt.gt.50) then
            write(6,*)'Convergence not in under 50 steps. Abort.'
            stop
         endif
         goto 380
      endif

      totmass =0.0d0

c  TMP - external
      do i=1,nspec
         totmass=totmass+amass(i)*nsp(i)
      enddo

      COMke=Pcomx**2+Pcomy**2+Pcomz**2
      COMke=COMke/(2.0d0*totmass)
      
      if(COMke.gt.free*trantkb*1.0d-30) then

         Vcomx=Pcomx/totmass
         Vcomy=Pcomy/totmass
         Vcomz=Pcomz/totmass

         Pcomx=0.0d0
         Pcomy=0.0d0
         Pcomz=0.0d0
         do 950 i=1,num
            ipoint=ntype(i)
            vxn(i)=vxn(i)-Vcomx
            vyn(i)=vyn(i)-Vcomy
            vzn(i)=vzn(i)-Vcomz
            Pcomx=Pcomx+vxn(i)*amass(ipoint)
            Pcomy=Pcomy+vyn(i)*amass(ipoint)
            Pcomz=Pcomz+vzn(i)*amass(ipoint)
  950    continue
 
         write(6,*) 'COM reset',COMke,free*trantkb*1.0d-30
     x     ,Pcomx,Pcomy,Pcomz

cc          goto 380
      endif         

 
c====================================================================
c Calculate the kinetic energy at time t. (i.e vmag() is the mean
c of the leapfroged v & vn velocities.)
c


c....  Note, these velocities are now in the cartesian frame
c.....and current with the positions and dipoles. i.e. do the energy
c......conservation at this point
      do 700 i=1,num

         vmag(i)=(vxn(i)*vxn(i))
     x          +(vyn(i)*vyn(i))
     x          +(vzn(i)*vzn(i))
         vtx(i)=vxn(i)
         vty(i)=vyn(i)
         vtz(i)=vzn(i)
         vxn(i)=hi(1,1)*vtx(i)+hi(1,2)*vty(i)+hi(1,3)*vtz(i)
         vyn(i)=hi(2,1)*vtx(i)+hi(2,2)*vty(i)+hi(2,3)*vtz(i)
         vzn(i)=hi(3,1)*vtx(i)+hi(3,2)*vty(i)+hi(3,3)*vtz(i)
 700  continue

C??????????????VXN NOW CONTAIN NEW VELOCITIES IN THE CELL FRAME --
 
      l=1
      tranko=tranke
      tranke=0.0d0

      do 710 i=1,nspec

         do 720 j=1,nsp(i)

            tranke=tranke+(hmass(i)*vmag(l))
            l=l+1

 720     continue

 710  continue

      tcell =  (Wgo/2.0d0)*trvg3

      tvol = (W/2.0d0)*veps3**2

      tpzeta =(CUEp/2.0d0)*vpzeta3(1)**2
      tbzeta =(CUEb/2.0d0)*vbzeta3(1)**2
      do i=2,5
        tpzeta=tpzeta+(CUEp2/2.0d0)*vpzeta3(i)**2
        tbzeta=tbzeta+(CUEb2/2.0d0)*vbzeta3(i)**2
      enddo

      PeeVee= pext*vol3

      pzeta=(free)*trantkb*pzeta3(1)
      bzeta=dom*trantkb*bzeta3(1)
      do i=2,5
        pzeta=pzeta+trantkb*pzeta3(i)
        bzeta=bzeta+trantkb*bzeta3(i)
      enddo

c...just add the PE to these variables -- should give the conserved quantity

c......Pcomx, Pcomy, Pcomz are available to monitor the total momentum
 
      tke=tranke*gtrankin

      if (relaxconfig) then
         if (tranke.le.tranko) then
            call velkill
            tranke=0.0d0
         endif
      endif

      return

      end
