      subroutine realE_quadpim

c************************************************************
c   Calculates the real energy and forces at each site.
c   The calculated energy is the real part of the Ewald
c  summation plus the short-range and Van der Waals
c  corrections from the Fumi-Tosi potential.
c
c  April 2000 - corrected Ewald.
c             - review loops to do multi-polarzable.
c
c  Nov. 2000  - quadrupoles re-introduced.
c************************************************************

      implicit double precision (a-h,o-z)

      include 'common.inc'
 
c================================================
c
c  Local double precision arrays.
      double precision elecxsr(nummax),elecysr(nummax)
     x                ,eleczsr(nummax),
     x         tempvirx(nspmax,nspmax),tempviry(nspmax,nspmax),
     x         tempvirz(nspmax,nspmax),frrxlab(nummax),
     x         frrylab(nummax),frrzlab(nummax)

      do 5 i=1,num

         elecxsr(i)=0.0d0
         elecysr(i)=0.0d0
         eleczsr(i)=0.0d0

         elecxu(i)=0.0d0
         elecyu(i)=0.0d0
         eleczu(i)=0.0d0

         srxx(i)=0.0d0
         sryy(i)=0.0d0
         srzz(i)=0.0d0
         srxy(i)=0.0d0
         srxz(i)=0.0d0
         sryz(i)=0.0d0

         exxsr(i)=0.0d0
         eyysr(i)=0.0d0
         ezzsr(i)=0.0d0
         exysr(i)=0.0d0
         exzsr(i)=0.0d0
         eyzsr(i)=0.0d0

    5 continue
c
c Zero energy accumulators.
c
      eng=0.0d0

      do 6 i=1,nspec

         srdipeng(i)=0.0d0
         srquadeng(i)=0.0d0
         do 7 j=1,nspec

            sreng(i,j)=0.0d0
            ddeng(i,j)=0.0d0
            dqeng(i,j)=0.0d0
            tempvirx(i,j)=0.0d0
            tempviry(i,j)=0.0d0
            tempvirz(i,j)=0.0d0

    7    continue

    6 continue
c
c Zero pressure calculation accumulators
c
      srvirft=0.0d0
      srdipvir=0.0d0
      virtot=0.0d0

      virtotx=0.0d0
      virtoty=0.0d0
      virtotz=0.0d0

      xviraccftsrc=0.0d0
      xviraccsrdip=0.0d0
      yviraccftsrc=0.0d0
      yviraccsrdip=0.0d0
      zviraccftsrc=0.0d0
      zviraccsrdip=0.0d0

      encoul=0.0d0

      cpeacc=0.0d0
      toteng=engpetot
      erfcracc=0.0d0
      qdipacc=0.0d0
      dipdipacc=0.0d0
      qquadacc=0.0d0
      dipquadacc=0.0d0
      quadquadacc=0.0d0
c
c zero the short range part, the other bits are initialized in recipE:
c
      stpsrxx=0.0d0
      stpsrxy=0.0d0
      stpsrxz=0.0d0
      stpsryy=0.0d0
      stpsryz=0.0d0
      stpsrzz=0.0d0

      stsrxx=0.0d0
      stsrxy=0.0d0
      stsrxz=0.0d0
      stsryy=0.0d0
      stsryz=0.0d0
      stsrzz=0.0d0

c zero the quadrupole-multipole contributions to the stress tensor.
c May have to move the `stpqquadxx'-like terms when the k-space part is done.
      stpqquadxx=0.0d0
      stpqquadxy=0.0d0
      stpqquadxz=0.0d0
      stpqquadyy=0.0d0
      stpqquadyz=0.0d0
      stpqquadzz=0.0d0

      stpdipquadxx=0.0d0
      stpdipquadxy=0.0d0
      stpdipquadxz=0.0d0
      stpdipquadyy=0.0d0
      stpdipquadyz=0.0d0
      stpdipquadzz=0.0d0

      stpquadquadxx=0.0d0
      stpquadquadxy=0.0d0
      stpquadquadxz=0.0d0
      stpquadquadyy=0.0d0
      stpquadquadyz=0.0d0
      stpquadquadzz=0.0d0
c
c Loop over anion-anion pairs calculating the permanent charge-permanent
c charge Fumi-Tosi potentials and forces.
c
      do 30 i=2,num

         ipoint=ntype(i)
         qirec=1.0d0/q(i)

         do 40 j=1,i-1

            jpoint=ntype(j)
            qjrec=1.0d0/q(j)

c
c The arrays x(), y() and z() now contain the system coordinates
c in the cell frame (hence the cf in the variable name).
c
            dxcf=x(i)-x(j)
            dycf=y(i)-y(j)
            dzcf=z(i)-z(j)
c
c Apply minimum image convention:
c
            dxcf=dxcf-boxlenx*int(dxcf*halfboxxrec)
            dycf=dycf-boxleny*int(dycf*halfboxyrec)
            dzcf=dzcf-boxlenz*int(dzcf*halfboxzrec)
c
c Transform separations into Cartesian coordinates.
c ie: the lab frame.
c
cc            dxsav(i,j)=hlab2(1,1)*dxcf
cc     x                +hlab2(1,2)*dycf+hlab2(1,3)*dzcf
cc            dysav(i,j)=hlab2(2,1)*dxcf
cc     x                +hlab2(2,2)*dycf+hlab2(2,3)*dzcf
cc            dzsav(i,j)=hlab2(3,1)*dxcf
cc     x                +hlab2(3,2)*dycf+hlab2(3,3)*dzcf

cc            dxsav(j,i)=-dxsav(i,j)
cc            dysav(j,i)=-dysav(i,j)
cc            dzsav(j,i)=-dzsav(i,j)

cc            drsq=dxsav(i,j)*dxsav(i,j)
cc cc    x          +dysav(i,j)*dysav(i,j)
cc     x          +dzsav(i,j)*dzsav(i,j)
c
c implement cut-off:
cc             drsav(i,j)=0.0d0
cc            drsav(j,i)=0.0d0
cc            if (drsq.ge.rsqmax) goto 40

            drsqrec=1.0d0/drsq
            drsq3rec=drsqrec*drsqrec*drsqrec
            dr=dsqrt(drsq)
cc            drsav(i,j)=dr
cc            drsav(j,i)=dr
            drrec=1.0d0/dr
            dr3rec=drrec*drsqrec
            dr4rec=drsqrec*drsqrec
            dr5rec=dr4rec*drrec
            dr7rec=dr5rec*drsqrec
            dr9rec=dr7rec*drsqrec
            dr11rec=dr9rec*drsqrec

            chgtem=q(i)*q(j)
            qq2api=etapi*chgtem
            expewld=exp(-etasq*drsq)
            erfc=erfunc(eta*dr)
            erfcr=chgtem*erfc*drrec
            egam=(erfcr+qq2api*expewld)*drsqrec
c
c Calculate terms for the dispersion (dd) damping function.
c
            dampsum=1.0d0+(dddamp(ipoint,jpoint)*dr)
     x               +(dddamp2(ipoint,jpoint)*drsq)
     x               +(dddamp3(ipoint,jpoint)*drsq*dr)
     x               +(dddamp4(ipoint,jpoint)*drsq*drsq)
     x               +(dddamp5(ipoint,jpoint)*drsq*drsq*dr)
     x               +(dddamp6(ipoint,jpoint)*drsq*drsq*drsq)
            dddampexp=dexp(-dddamp(ipoint,jpoint)*dr)
            f6=1.0d0-(dddampexp*dampsum)
            dddampforce=dddamp(ipoint,jpoint)
     x                 *dddamp6(ipoint,jpoint)
     x                 *drsq*drsq*drsq*dddampexp
c
c Calculate terms for the dispersion (dq) damping function.
c
            dampsumdq=1.0d0+(dqdamp(ipoint,jpoint)*dr)
     x               +(dqdamp2(ipoint,jpoint)*drsq)
     x               +(dqdamp3(ipoint,jpoint)*drsq*dr)
     x               +(dqdamp4(ipoint,jpoint)*drsq*drsq)
     x               +(dqdamp5(ipoint,jpoint)*drsq*drsq*dr)
     x               +(dqdamp6(ipoint,jpoint)*drsq*drsq*drsq)
     x               +(dqdamp7(ipoint,jpoint)*drsq*drsq*drsq*dr)
     x               +(dqdamp8(ipoint,jpoint)*drsq*drsq*drsq*drsq)
            dqdampexp=dexp(-dqdamp(ipoint,jpoint)*dr)
            f8=1.0d0-(dqdampexp*dampsumdq)
            dqdampforce=dqdamp(ipoint,jpoint)
     x                 *dqdamp8(ipoint,jpoint)
     x                 *drsq*drsq*drsq*drsq*dqdampexp

c
c Updated equation of motion to include new terms:
c

            sgam=(
     x          -(6.0d0*ftc(ipoint,jpoint)*f6
     x          -(ftc(ipoint,jpoint)*dddampforce*dr)
     x          +8.0d0*ftd(ipoint,jpoint)*f8*drsqrec
     x          -(ftd(ipoint,jpoint)*dqdampforce*drrec))
     x          *drsq3rec)*drsqrec

            fttx=sgam*dxsav(i,j)
            ftty=sgam*dysav(i,j)
            fttz=sgam*dzsav(i,j)

            cpe=-(ftc(ipoint,jpoint)*f6
     x             +ftd(ipoint,jpoint)*f8*drsqrec)*drsq3rec
            cpeacc=cpeacc+cpe
c
c energy accumulators for analysis purposes only:
c
            sreng(ipoint,jpoint)=sreng(ipoint,jpoint)+expft
            ddeng(ipoint,jpoint)=ddeng(ipoint,jpoint)
     x      -ftc(ipoint,jpoint)
     x      *f6*drsq3rec
            dqeng(ipoint,jpoint)=dqeng(ipoint,jpoint)-
     x      ftd(ipoint,jpoint)*f8*drsq3rec*drsqrec
c
c Don't include sgam here, as we have separated off the s-r contribtuion.
c
            gamtot=egam+sgam

            ettx=egam*dxsav(i,j)
            etty=egam*dysav(i,j)
            ettz=egam*dzsav(i,j)

            xviraccftsrc=xviraccftsrc+(sgam*dxsav(i,j)*dxsav(i,j))
            yviraccftsrc=yviraccftsrc+(sgam*dysav(i,j)*dysav(i,j))
            zviraccftsrc=zviraccftsrc+(sgam*dzsav(i,j)*dzsav(i,j))

            elecx(j)=elecx(j)-(ettx*qjrec)
            elecy(j)=elecy(j)-(etty*qjrec)
            elecz(j)=elecz(j)-(ettz*qjrec)

            elecx(i)=elecx(i)+(ettx*qirec)
            elecy(i)=elecy(i)+(etty*qirec)
            elecz(i)=elecz(i)+(ettz*qirec)
c
c Start of dipole handling section:
c
            xmuidotrij=xmu(i)*dxsav(i,j)
     x                +ymu(i)*dysav(i,j)
     x                +zmu(i)*dzsav(i,j)
            xmujdotrij=xmu(j)*dxsav(i,j)
     x                +ymu(j)*dysav(i,j)
     x                +zmu(j)*dzsav(i,j)
c
c Permanent charge - dipole energy:
c
            efact=etapi*expewld*dr+erfc
            chgdipengij=(q(j)*xmuidotrij)*dr3rec
            chgdipengji=-(q(i)*xmujdotrij)*dr3rec
            chgdipeng=(chgdipengij+chgdipengji)*efact

c
c Permanent charge - dipole forces:
c
            qdotmuxij=q(j)*xmu(i)
            qdotmuxji=-q(i)*xmu(j)
            qdotmuyij=q(j)*ymu(i)
            qdotmuyji=-q(i)*ymu(j)
            qdotmuzij=q(j)*zmu(i)
            qdotmuzji=-q(i)*zmu(j)

            efact2=2.0d0*etasq*drsqrec*etapi*expewld
c
c the sign difference in the efact2 term comes from the behaviour
c of this term as you invert from ion i to ion j, it does not
c change sign:
c
            frrxdipqij=(qdotmuxij*dr3rec
     x                 -3.0d0*dxsav(i,j)*drsqrec*chgdipengij)*efact
     x                 -dxsav(i,j)*xmuidotrij*q(j)*efact2
            frrxdipqji=(qdotmuxji*dr3rec
     x                 -3.0d0*dxsav(i,j)*drsqrec*chgdipengji)*efact
     x                 +dxsav(i,j)*xmujdotrij*q(i)*efact2
            frrydipqij=(qdotmuyij*dr3rec
     x                 -3.0d0*dysav(i,j)*drsqrec*chgdipengij)*efact
     x                 -dysav(i,j)*xmuidotrij*q(j)*efact2
            frrydipqji=(qdotmuyji*dr3rec
     x                 -3.0d0*dysav(i,j)*drsqrec*chgdipengji)*efact
     x                 +dysav(i,j)*xmujdotrij*q(i)*efact2
            frrzdipqij=(qdotmuzij*dr3rec
     x                 -3.0d0*dzsav(i,j)*drsqrec*chgdipengij)*efact
     x                 -dzsav(i,j)*xmuidotrij*q(j)*efact2
            frrzdipqji=(qdotmuzji*dr3rec
     x                 -3.0d0*dzsav(i,j)*drsqrec*chgdipengji)*efact
     x                 +dzsav(i,j)*xmujdotrij*q(i)*efact2

            elecx(j)=elecx(j)-frrxdipqij*qjrec
            elecy(j)=elecy(j)-frrydipqij*qjrec
            elecz(j)=elecz(j)-frrzdipqij*qjrec

            elecx(i)=elecx(i)+frrxdipqji*qirec
            elecy(i)=elecy(i)+frrydipqji*qirec
            elecz(i)=elecz(i)+frrzdipqji*qirec

            frrxdipq=frrxdipqij+frrxdipqji
            frrydipq=frrydipqij+frrydipqji
            frrzdipq=frrzdipqij+frrzdipqji

            egamchg=egam/chgtem

c
c Dipole-dipole interaction tensor:
c
            dxunitsq=dxsav(i,j)*dxsav(i,j)*drsqrec
            dyunitsq=dysav(i,j)*dysav(i,j)*drsqrec
            dzunitsq=dzsav(i,j)*dzsav(i,j)*drsqrec

            sxx=1.0d0-3.0d0*dxunitsq
            syy=1.0d0-3.0d0*dyunitsq
            szz=1.0d0-3.0d0*dzunitsq
            sxy=-3.0d0*dxsav(i,j)*dysav(i,j)*drsqrec
            sxz=-3.0d0*dxsav(i,j)*dzsav(i,j)*drsqrec
            syz=-3.0d0*dzsav(i,j)*dysav(i,j)*drsqrec

            txx=-sxx*dr3rec
            tyy=-syy*dr3rec
            tzz=-szz*dr3rec
            txy=-sxy*dr3rec
            txz=-sxz*dr3rec
            tyz=-syz*dr3rec

            xmuidotT=xmu(i)*sxx
     x              +ymu(i)*sxy 
     x              +zmu(i)*sxz 
            ymuidotT=xmu(i)*sxy 
     x              +ymu(i)*syy 
     x              +zmu(i)*syz 
            zmuidotT=xmu(i)*sxz 
     x              +ymu(i)*syz 
     x              +zmu(i)*szz 
c
c Dipole-dipole energy (ignoring the pbcs):
c
            dipdipeng=xmuidotT*xmu(j)
     x               +ymuidotT*ymu(j)
     x               +zmuidotT*zmu(j)

c**** Saves -T tensor for use in sr dip dip term.
            fdipdipdamp=dipdipeng

c
c second term corrects for pbcs:
c
            dipdipeng=dipdipeng*efact*dr3rec-
     x                 xmuidotrij*xmujdotrij*efact2
c

c========================================================
c Dipole-quadrupole interaction tensor (see Buckingham &
c theory book 3, p 135 -).
c There are 10 of these.
c
c this part is not 'ewalded':
c
            threedr7=3.0d0*dr7rec
            txxx=(5.0d0*dxsav(i,j)*dxsav(i,j)*dxsav(i,j)
     x          -3.0d0*dxsav(i,j)*drsq)*threedr7
            tyyy=(5.0d0*dysav(i,j)*dysav(i,j)*dysav(i,j)
     x          -3.0d0*dysav(i,j)*drsq)*threedr7
            tzzz=(5.0d0*dzsav(i,j)*dzsav(i,j)*dzsav(i,j)
     x          -3.0d0*dzsav(i,j)*drsq)*threedr7

            txxy=(5.0d0*dxsav(i,j)*dxsav(i,j)*dysav(i,j)
     x          -drsq*dysav(i,j))*threedr7
            txyy=(5.0d0*dxsav(i,j)*dysav(i,j)*dysav(i,j)
     x          -drsq*dxsav(i,j))*threedr7
            txxz=(5.0d0*dxsav(i,j)*dxsav(i,j)*dzsav(i,j)
     x          -drsq*dzsav(i,j))*threedr7
            txzz=(5.0d0*dxsav(i,j)*dzsav(i,j)*dzsav(i,j)
     x          -drsq*dxsav(i,j))*threedr7
            tyyz=(5.0d0*dysav(i,j)*dysav(i,j)*dzsav(i,j)
     x          -drsq*dzsav(i,j))*threedr7
            tyzz=(5.0d0*dysav(i,j)*dzsav(i,j)*dzsav(i,j)
     x          -drsq*dysav(i,j))*threedr7

            txyz=15.0d0*dr7rec*dxsav(i,j)*dysav(i,j)*dzsav(i,j)
c
c Dipole-dipole forces (see theory book 3, p92-). N.B. P.B.Cs are NOW
c INCLUDED -- NOTE CORRESPONDING CHANGE IN RECIPE
 
            fmumux=(xmu(j)*(xmu(i)*txxx
     x                     +ymu(i)*txxy
     x                     +zmu(i)*txxz))
     x            +(ymu(j)*(xmu(i)*txxy
     x                     +ymu(i)*txyy
     x                     +zmu(i)*txyz))
     x            +(zmu(j)*(xmu(i)*txxz
     x                     +ymu(i)*txyz
     x                     +zmu(i)*txzz))

            fmumuy=(xmu(j)*(xmu(i)*txxy
     x                     +ymu(i)*txyy
     x                     +zmu(i)*txyz))
     x            +(ymu(j)*(xmu(i)*txyy
     x                     +ymu(i)*tyyy
     x                     +zmu(i)*tyyz))
     x            +(zmu(j)*(xmu(i)*txyz
     x                     +ymu(i)*tyyz
     x                     +zmu(i)*tyzz))

            fmumuz=(xmu(j)*(xmu(i)*txxz
     x                     +ymu(i)*txyz
     x                     +zmu(i)*txzz))
     x            +(ymu(j)*(xmu(i)*txyz
     x                     +ymu(i)*tyyz
     x                     +zmu(i)*tyzz))
     x            +(zmu(j)*(xmu(i)*txzz
     x                     +ymu(i)*tyzz
     x                     +zmu(i)*tzzz))


            efact3= 2.0d0*efact2*(drsqrec+etasq)

            fmumux = -fmumux*efact
     x               +fdipdipdamp*efact2*dxsav(i,j)
     x               +xmu(i)*xmujdotrij*efact2
     x               +xmu(j)*xmuidotrij*efact2
     x               -xmuidotrij*xmujdotrij*dxsav(i,j)*efact3

            fmumuy = -fmumuy*efact
     x               +fdipdipdamp*efact2*dysav(i,j)
     x               +ymu(i)*xmujdotrij*efact2
     x               +ymu(j)*xmuidotrij*efact2
     x               -xmuidotrij*xmujdotrij*dysav(i,j)*efact3

            fmumuz = -fmumuz*efact
     x               +fdipdipdamp*efact2*dzsav(i,j)
     x               +zmu(i)*xmujdotrij*efact2
     x               +zmu(j)*xmuidotrij*efact2
     x               -xmuidotrij*xmujdotrij*dzsav(i,j)*efact3


c===============================================
c  accumulate all the real-space coulombic contributions to the
c  energy (not short-range dipoles)

            eng=eng+erfcr+cpe-chgdipeng+dipdipeng

c=======================================================
c
c update accumulators:
c
            erfcracc=erfcracc+erfcr
            qdipacc=qdipacc-chgdipeng
            dipdipacc=dipdipacc+dipdipeng

            stsrxx=stsrxx+fttx*dxsav(i,j)
            stsrxy=stsrxy+fttx*dysav(i,j)
            stsrxz=stsrxz+fttx*dzsav(i,j)
            stsryy=stsryy+ftty*dysav(i,j)
            stsryz=stsryz+ftty*dzsav(i,j)
            stsrzz=stsrzz+fttz*dzsav(i,j)

c
c charge-charge term:
c
            stcxx=stcxx+egam*dxsav(i,j)*dxsav(i,j)
            stcxy=stcxy+egam*dxsav(i,j)*dysav(i,j)
            stcxz=stcxz+egam*dxsav(i,j)*dzsav(i,j)
            stcyy=stcyy+egam*dysav(i,j)*dysav(i,j)
            stcyz=stcyz+egam*dysav(i,j)*dzsav(i,j)
            stczz=stczz+egam*dzsav(i,j)*dzsav(i,j)
c
c charge-dipole term:
c
            stpxx=stpxx+frrxdipq*dxsav(i,j)
            stpxy=stpxy+0.5d0*(frrxdipq*dysav(i,j)
     x                        +frrydipq*dxsav(i,j))
            stpxz=stpxz+0.5d0*(frrxdipq*dzsav(i,j)
     x                        +frrzdipq*dxsav(i,j))
            stpyy=stpyy+frrydipq*dysav(i,j)
            stpyz=stpyz+0.5d0*(frrydipq*dzsav(i,j)
     x                        +frrzdipq*dysav(i,j))
            stpzz=stpzz+frrzdipq*dzsav(i,j)
c
c dipole-dipole term:
c
            stp2xx=stp2xx+fmumux*dxsav(i,j)
            stp2xy=stp2xy+0.5d0*(fmumux*dysav(i,j)
     x                          +fmumuy*dxsav(i,j))
            stp2xz=stp2xz+0.5d0*(fmumux*dzsav(i,j)
     x                          +fmumuz*dxsav(i,j))
            stp2yy=stp2yy+fmumuy*dysav(i,j)
            stp2yz=stp2yz+0.5d0*(fmumuy*dzsav(i,j)
     x                          +fmumuz*dysav(i,j))
            stp2zz=stp2zz+fmumuz*dzsav(i,j)
c
c virial due to fumi-tosi s-r forces + dispersion
c
cc            tempvirx(1,1)=tempvirx(1,1)+fttx*dx
cc            tempviry(1,1)=tempviry(1,1)+ftty*dy
cc            tempvirz(1,1)=tempvirz(1,1)+fttz*dz
c            xviraccftsrc=xviraccftsrc+fttx*dx
c            yviraccftsrc=yviraccftsrc+ftty*dy
c            zviraccftsrc=zviraccftsrc+fttz*dz
c
c Sum together all forces.
c
c            fttx=fttx+egam*dx+frrxdipq+fmumux
c            ftty=ftty+egam*dy+frrydipq+fmumuy
c            fttz=fttz+egam*dz+frrzdipq+fmumuz

            fttx=gamtot*dxsav(i,j)+frrxdipq+fmumux
            ftty=gamtot*dysav(i,j)+frrydipq+fmumuy
            fttz=gamtot*dzsav(i,j)+frrzdipq+fmumuz

            frrx(j)=frrx(j)-fttx
            frry(j)=frry(j)-ftty
            frrz(j)=frrz(j)-fttz

            frrx(i)=frrx(i)+fttx
            frry(i)=frry(i)+ftty
            frrz(i)=frrz(i)+fttz

c***********************************************************
c***********************************************************
c***********************************************************
c***********************************************************
c End of dipole only interaction calcs.
c Start of quadrupole handling section.
c
            dx=dxsav(i,j)
            dy=dysav(i,j)
            dz=dzsav(i,j)

            dxx=dx*dx
            dyy=dy*dy
            dzz=dz*dz
            dxy=dx*dy
            dxz=dx*dz
            dyz=dy*dz
c
c quadrupole-charge interaction:
c   theta_alpha_beta(i)*T_alpha_beta*q(j)+
c   q(i)*T_alpha_beta*theta_alpha_beta(j)
c
            chgquadeng=((txx*quadxx(i)
     x                  +txy*quadxy(i)
     x                  +txy*quadxy(i)
     x                  +txz*quadxz(i)
     x                  +txz*quadxz(i)
     x                  +tyy*quadyy(i)
     x                  +tyz*quadyz(i)
     x                  +tzz*quadzz(i)
     x                  +tyz*quadyz(i))*q(j))
     x                  +((txx*quadxx(j)
     x                  +txy*quadxy(j)
     x                  +txy*quadxy(j)
     x                  +txz*quadxz(j)
     x                  +txz*quadxz(j)
     x                  +tyy*quadyy(j)
     x                  +tyz*quadyz(j)
     x                  +tzz*quadzz(j)
     x                  +tyz*quadyz(j))*q(i))

            chgquad2=((dxx*quadxx(i)
     x                  +dxy*quadxy(i)
     x                  +dxy*quadxy(i)
     x                  +dxz*quadxz(i)
     x                  +dxz*quadxz(i)
     x                  +dyy*quadyy(i)
     x                  +dyz*quadyz(i)
     x                  +dzz*quadzz(i)
     x                  +dyz*quadyz(i))*q(j))
     x                  +((dxx*quadxx(j)
     x                  +dxy*quadxy(j)
     x                  +dxy*quadxy(j)
     x                  +dxz*quadxz(j)
     x                  +dxz*quadxz(j)
     x                  +dyy*quadyy(j)
     x                  +dyz*quadyz(j)
     x                  +dzz*quadzz(j)
     x                  +dyz*quadyz(j))*q(i))

            chgquadeng=(chgquadeng*efact)+chgquad2*efact2
c
c factor of one-third from Buckingham.
c
            chgquadeng=chgquadeng*onethird
c
c Quadrupole-charge forces:
c  -theta_beta_gamma(i)*T_alpha_beta_gamma*q(j)
c
            fqquadxij=(q(j)*(quadxx(i)*txxx
     x                    +quadxy(i)*txxy
     x                    +quadxz(i)*txxz
     x                    +quadxy(i)*txxy
     x                    +quadyy(i)*txyy
     x                    +quadyz(i)*txyz
     x                    +quadxz(i)*txxz
     x                    +quadyz(i)*txyz
     x                    +quadzz(i)*txzz))
c
c  -q(i)*T_alpha_beta_gamma*theta_beta_gamma(j)
c
            fqquadxji=(q(i)*(quadxx(j)*txxx
     x                    +quadxy(j)*txxy
     x                    +quadxz(j)*txxz
     x                    +quadxy(j)*txxy
     x                    +quadyy(j)*txyy
     x                    +quadyz(j)*txyz
     x                    +quadxz(j)*txxz
     x                    +quadyz(j)*txyz
     x                    +quadzz(j)*txzz))

            fqquadyij=(q(j)*(quadxx(i)*txxy
     x                    +quadxy(i)*txyy
     x                    +quadxz(i)*txyz
     x                    +quadxy(i)*txyy
     x                    +quadyy(i)*tyyy
     x                    +quadyz(i)*tyyz
     x                    +quadxz(i)*txyz
     x                    +quadyz(i)*tyyz
     x                    +quadzz(i)*tyzz))
            fqquadyji=(q(i)*(quadxx(j)*txxy
     x                    +quadxy(j)*txyy
     x                    +quadxz(j)*txyz
     x                    +quadxy(j)*txyy
     x                    +quadyy(j)*tyyy
     x                    +quadyz(j)*tyyz
     x                    +quadxz(j)*txyz
     x                    +quadyz(j)*tyyz
     x                    +quadzz(j)*tyzz))

            fqquadzij=(q(j)*(quadxx(i)*txxz
     x                    +quadxy(i)*txyz
     x                    +quadxz(i)*txzz
     x                    +quadxy(i)*txyz
     x                    +quadyy(i)*tyyz
     x                    +quadyz(i)*tyzz
     x                    +quadxz(i)*txzz
     x                    +quadyz(i)*tyzz
     x                    +quadzz(i)*tzzz))
            fqquadzji=(q(i)*(quadxx(j)*txxz
     x                    +quadxy(j)*txyz
     x                    +quadxz(j)*txzz
     x                    +quadxy(j)*txyz
     x                    +quadyy(j)*tyyz
     x                    +quadyz(j)*tyzz
     x                    +quadxz(j)*txzz
     x                    +quadyz(j)*tyzz
     x                    +quadzz(j)*tzzz))
c
c Factor of onethird from Buckingham.
c
            fqquadxij=fqquadxij*onethird
            fqquadyij=fqquadyij*onethird
            fqquadzij=fqquadzij*onethird
            fqquadxji=fqquadxji*onethird
            fqquadyji=fqquadyji*onethird
            fqquadzji=fqquadzji*onethird

            elecx(j)=elecx(j)-fqquadxij*qjrec
            elecy(j)=elecy(j)-fqquadyij*qjrec
            elecz(j)=elecz(j)-fqquadzij*qjrec

            elecx(i)=elecx(i)+fqquadxji*qirec
            elecy(i)=elecy(i)+fqquadyji*qirec
            elecz(i)=elecz(i)+fqquadzji*qirec

            fqquadx=fqquadxij+fqquadxji
            fqquady=fqquadyij+fqquadyji
            fqquadz=fqquadzij+fqquadzji

            stpqquadxx=stpqquadxx+fqquadx*dxsav(i,j)
            stpqquadxy=stpqquadxy+0.5d0*(fqquadx*dysav(i,j)
     x                        +fqquady*dxsav(i,j))
            stpqquadxz=stpqquadxz+0.5d0*(fqquadx*dzsav(i,j)
     x                        +fqquadz*dxsav(i,j))
            stpqquadyy=stpqquadyy+fqquady*dysav(i,j)
            stpqquadyz=stpqquadyz+0.5d0*(fqquadz*dysav(i,j)
     x                        +fqquady*dzsav(i,j))
            stpqquadzz=stpqquadzz+fqquadz*dzsav(i,j)

c
c Dipole - quadrupole energy.
c
            quaddipengij=xmu(i)*(txxx*quadxx(j)
     x                +(2.0d0*txxy*quadxy(j))
     x                +(2.0d0*txxz*quadxz(j))
     x                +(txyy*quadyy(j))
     x                +(txzz*quadzz(j))
     x                +(2.0d0*txyz*quadyz(j)))
     x                +ymu(i)*(txxy*quadxx(j)
     x                +(2.0d0*txyy*quadxy(j))
     x                +(2.0d0*txyz*quadxz(j))
     x                +(tyyy*quadyy(j))
     x                +(tyzz*quadzz(j))
     x                +(2.0d0*tyyz*quadyz(j)))
     x                +zmu(i)*(txxz*quadxx(j)
     x                +(2.0d0*txyz*quadxy(j))
     x                +(2.0d0*txzz*quadxz(j))
     x                +(tyyz*quadyy(j))
     x                +(tzzz*quadzz(j))
     x                +(2.0d0*tyzz*quadyz(j)))

            quaddipengji=xmu(j)*(txxx*quadxx(i)
     x                +(2.0d0*txxy*quadxy(i))
     x                +(2.0d0*txxz*quadxz(i))
     x                +(txyy*quadyy(i))
     x                +(txzz*quadzz(i))
     x                +(2.0d0*txyz*quadyz(i)))
     x                +ymu(j)*(txxy*quadxx(i)
     x                +(2.0d0*txyy*quadxy(i))
     x                +(2.0d0*txyz*quadxz(i))
     x                +(tyyy*quadyy(i))
     x                +(tyzz*quadzz(i))
     x                +(2.0d0*tyyz*quadyz(i)))
     x                +zmu(j)*(txxz*quadxx(i)
     x                +(2.0d0*txyz*quadxy(i))
     x                +(2.0d0*txzz*quadxz(i))
     x                +(tyyz*quadyy(i))
     x                +(tzzz*quadzz(i))
     x                +(2.0d0*tyzz*quadyz(i)))

            quaddipeng=(quaddipengij-quaddipengji)*onethird
c
c Quadrupole-quadrupole interaction tensor (see MW's theory book 3
c p 137 - checked with MVdH so it must be right !).
c There are 15 of these.
c
            txxxx=(105.0d0*dx*dx*dx*dx*dr9rec)
     x           +(9.0d0*dr5rec)
     x           -(90.0d0*dx*dx*dr7rec)
            tyyyy=(105.0d0*dy*dy*dy*dy*dr9rec)
     x           +(9.0d0*dr5rec)
     x           -(90.0d0*dy*dy*dr7rec)
            tzzzz=(105.0d0*dz*dz*dz*dz*dr9rec)
     x           +(9.0d0*dr5rec)
     x           -(90.0d0*dz*dz*dr7rec)

            txxxy=(105.0d0*dx*dx*dx*dy*dr9rec)
     x           -(45.0d0*dx*dy*dr7rec)
            txxxz=(105.0d0*dx*dx*dx*dz*dr9rec)
     x           -(45.0d0*dx*dz*dr7rec)
            txyyy=(105.0d0*dx*dy*dy*dy*dr9rec)
     x           -(45.0d0*dx*dy*dr7rec)
            txzzz=(105.0d0*dx*dz*dz*dz*dr9rec)
     x           -(45.0d0*dx*dz*dr7rec)
            tyyyz=(105.0d0*dz*dy*dy*dy*dr9rec)
     x           -(45.0d0*dz*dy*dr7rec)
            tyzzz=(105.0d0*dz*dz*dz*dy*dr9rec)
     x           -(45.0d0*dz*dy*dr7rec)

            txxyy=(105.0d0*dx*dx*dy*dy*dr9rec)
     x           +(3.0d0*dr5rec)
     x           -(15.0d0*((dx*dx)+(dy*dy))*dr7rec)
            txxzz=(105.0d0*dx*dx*dz*dz*dr9rec)
     x           +(3.0d0*dr5rec)
     x           -(15.0d0*((dx*dx)+(dz*dz))*dr7rec)
            tyyzz=(105.0d0*dz*dz*dy*dy*dr9rec)
     x           +(3.0d0*dr5rec)
     x           -(15.0d0*((dz*dz)+(dy*dy))*dr7rec)

            txxyz=(105.0d0*dx*dx*dy*dz*dr9rec)
     x           -(15.0d0*dy*dz*dr7rec)
            txyyz=(105.0d0*dx*dy*dy*dz*dr9rec)
     x           -(15.0d0*dx*dz*dr7rec)
            txyzz=(105.0d0*dx*dz*dy*dz*dr9rec)
     x           -(15.0d0*dy*dx*dr7rec)
c
c Dipole - quadrupole forces.
c
            fdipquadijx=xmu(i)*(txxxx*quadxx(j)
     x                  +(2.0d0*txxxy*quadxy(j))
     x                  +(2.0d0*txxxz*quadxz(j))
     x                  +(txxyy*quadyy(j))
     x                  +(txxzz*quadzz(j))
     x                  +(2.0d0*txxyz*quadyz(j)))
     x          +ymu(i)*((txxxy*quadxx(j))
     x                  +(2.0d0*txxyy*quadxy(j))
     x                  +(2.0d0*txxyz*quadxz(j))
     x                         +txyyy*quadyy(j)
     x                  +(txyzz*quadzz(j))
     x                  +(2.0d0*txyyz*quadyz(j)))
     x          +zmu(i)*((txxxz*quadxx(j))
     x                  +(2.0d0*txxyz*quadxy(j))
     x                  +(2.0d0*txxzz*quadxz(j))
     x                  +(txyyz*quadyy(j))
     x                         +txzzz*quadzz(j)
     x                  +(2.0d0*txyzz*quadyz(j)))

            fdipquadjix=xmu(j)*(txxxx*quadxx(i)
     x                  +(2.0d0*txxxy*quadxy(i))
     x                  +(2.0d0*txxxz*quadxz(i))
     x                  +(txxyy*quadyy(i))
     x                  +(txxzz*quadzz(i))
     x                  +(2.0d0*txxyz*quadyz(i)))
     x          +ymu(j)*((txxxy*quadxx(i))
     x                  +(2.0d0*txxyy*quadxy(i))
     x                  +(2.0d0*txxyz*quadxz(i))
     x                         +txyyy*quadyy(i)
     x                  +(txyzz*quadzz(i))
     x                  +(2.0d0*txyyz*quadyz(i)))
     x          +zmu(j)*((txxxz*quadxx(i))
     x                  +(2.0d0*txxyz*quadxy(i))
     x                  +(2.0d0*txxzz*quadxz(i))
     x                  +(txyyz*quadyy(i))
     x                         +txzzz*quadzz(i)
     x                  +(2.0d0*txyzz*quadyz(i)))

            fdipquadijy=xmu(i)*(txxxy*quadxx(j)
     x                  +(2.0d0*txxyy*quadxy(j))
     x                  +(2.0d0*txxyz*quadxz(j))
     x                  +(txyyy*quadyy(j))
     x                  +(txyzz*quadzz(j))
     x                  +(2.0d0*txyyz*quadyz(j)))
     x          +ymu(i)*((txxyy*quadxx(j))
     x                  +(txyyy*quadxy(j))
     x                  +(2.0d0*txyyz*quadxz(j))
     x                         +tyyyy*quadyy(j)
     x                  +(tyyzz*quadzz(j))
     x                  +(2.0d0*tyyyz*quadyz(j)))
     x          +zmu(i)*((txxyz*quadxx(j))
     x                  +(2.0d0*txyyz*quadxy(j))
     x                  +(2.0d0*txyzz*quadxz(j))
     x                  +(tyyyz*quadyy(j))
     x                         +tyzzz*quadzz(j)
     x                  +(2.0d0*tyyzz*quadyz(j)))

            fdipquadjiy=xmu(j)*(txxxy*quadxx(i)
     x                  +(2.0d0*txxyy*quadxy(i))
     x                  +(2.0d0*txxyz*quadxz(i))
     x                  +(txyyy*quadyy(i))
     x                  +(txyzz*quadzz(i))
     x                  +(2.0d0*txyyz*quadyz(i)))
     x          +ymu(j)*((txxyy*quadxx(i))
     x                  +(2.0d0*txyyy*quadxy(i))
     x                  +(2.0d0*txyyz*quadxz(i))
     x                         +tyyyy*quadyy(i)
     x                  +(tyyzz*quadzz(i))
     x                  +(2.0d0*tyyyz*quadyz(i)))
     x          +zmu(j)*((txxyz*quadxx(i))
     x                  +(2.0d0*txyyz*quadxy(i))
     x                  +(2.0d0*txyzz*quadxz(i))
     x                  +(tyyyz*quadyy(i))
     x                         +tyzzz*quadzz(i)
     x                  +(2.0d0*tyyzz*quadyz(i)))

            fdipquadijz=xmu(i)*(txxxz*quadxx(j)
     x                  +(2.0d0*txxyz*quadxy(j))
     x                  +(2.0d0*txxzz*quadxz(j))
     x                  +(txyyz*quadyy(j))
     x                  +(txzzz*quadzz(j))
     x                   +(2.0d0*txyzz*quadyz(j)))
     x          +ymu(i)*((txxyz*quadxx(j))
     x                  +(2.0d0*txyyz*quadxy(j))
     x                  +(2.0d0*txyzz*quadxz(j))
     x                         +tyyyz*quadyy(j)
     x                  +(tyzzz*quadzz(j))
     x                  +(2.0d0*tyyzz*quadyz(j)))
     x          +zmu(i)*((txxzz*quadxx(j))
     x                  +(2.0d0*txyzz*quadxy(j))
     x                  +(2.0d0*txzzz*quadxz(j))
     x                  +(tyyzz*quadyy(j))
     x                         +tzzzz*quadzz(j)
     x                  +(2.0d0*tyzzz*quadyz(j)))

            fdipquadjiz=xmu(j)*(txxxz*quadxx(i)
     x                  +(2.0d0*txxyz*quadxy(i))
     x                  +(2.0d0*txxzz*quadxz(i))
     x                  +(txyyz*quadyy(i))
     x                  +(txzzz*quadzz(i))
     x                         +(2.0d0*txyzz*quadyz(i)))
     x          +ymu(j)*((txxyz*quadxx(i))
     x                  +(2.0d0*txyyz*quadxy(i))
     x                         +(2.0d0*txyzz*quadxz(i))
     x                         +tyyyz*quadyy(i)
     x                  +(tyzzz*quadzz(i))
     x                  +(2.0d0*tyyzz*quadyz(i)))
     x          +zmu(j)*((txxzz*quadxx(i))
     x                         +(2.0d0*txyzz*quadxy(i))
     x                  +(2.0d0*txzzz*quadxz(i))
     x                  +(tyyzz*quadyy(i))
     x                         +tzzzz*quadzz(i)
     x                  +(2.0d0*tyzzz*quadyz(i)))

            fdipquadx=onethird*(fdipquadijx-fdipquadjix)
            fdipquady=onethird*(fdipquadijy-fdipquadjiy)
            fdipquadz=onethird*(fdipquadijz-fdipquadjiz)

            stpdipquadxx=stpdipquadxx+fdipquadx*dxsav(i,j)
            stpdipquadxy=stpdipquadxy+0.5d0*(fdipquadx*dysav(i,j)
     x                        +fdipquady*dxsav(i,j))
            stpdipquadxz=stpdipquadxz+0.5d0*(fdipquadx*dzsav(i,j)
     x                        +fdipquadz*dxsav(i,j))
            stpdipquadyy=stpdipquadyy+fdipquady*dysav(i,j)
            stpdipquadyz=stpdipquadyz+0.5d0*(fdipquadz*dysav(i,j)
     x                        +fdipquady*dzsav(i,j))
            stpdipquadzz=stpdipquadzz+fdipquadz*dzsav(i,j)

c
c Quadrupole-quadrupole energy.
c
            quadquadengxx=quadxx(i)*(txxxx*quadxx(j)
     x                    +(2.0d0*txxxy*quadxy(j))
     x                    +(2.0d0*txxxz*quadxz(j))
     x                    +(2.0d0*txxyz*quadyz(j))
     x                    +txxyy*quadyy(j)
     x                    +txxzz*quadzz(j))
            quadquadengyy=quadyy(i)*(txxyy*quadxx(j)
     x                    +(2.0d0*txyyy*quadxy(j))
     x                    +(2.0d0*txyyz*quadxz(j))
     x                    +(2.0d0*tyyyz*quadyz(j))
     x                    +tyyyy*quadyy(j)
     x                    +tyyzz*quadzz(j))
            quadquadengzz=quadzz(i)*(txxzz*quadxx(j)
     x                    +(2.0d0*txyzz*quadxy(j))
     x                    +(2.0d0*txzzz*quadxz(j))
     x                    +(2.0d0*tyzzz*quadyz(j))
     x                    +tyyzz*quadyy(j)
     x                    +tzzzz*quadzz(j))
            quadquadengxy=quadxy(i)*((2.0d0*txxxy*quadxx(j))
     x                    +(4.0d0*txxyy*quadxy(j))
     x                    +(4.0d0*txxyz*quadxz(j))
     x                    +(2.0d0*txyyy*quadyy(j))
     x                    +(4.0d0*txyyz*quadyz(j))
     x                    +(2.0d0*txyzz*quadzz(j)))
            quadquadengxz=quadxz(i)*((2.0d0*txxxz*quadxx(j))
     x                    +(4.0d0*txxyz*quadxy(j))
     x                    +(4.0d0*txxzz*quadxz(j))
     x                    +(2.0d0*txyyz*quadyy(j))
     x                    +(4.0d0*txyzz*quadyz(j))
     x                    +(2.0d0*txzzz*quadzz(j)))
            quadquadengyz=quadyz(i)*((2.0d0*txxyz*quadxx(j))
     x                    +(4.0d0*txyyz*quadxy(j))
     x                    +(4.0d0*txyzz*quadxz(j))
     x                    +(2.0d0*tyyyz*quadyy(j))
     x                    +(4.0d0*tyyzz*quadyz(j))
     x                    +(2.0d0*tyzzz*quadzz(j)))

            quadquadeng=quadquadengxx+quadquadengyy
     x                 +quadquadengzz+quadquadengxy
     x                 +quadquadengxz+quadquadengyz

            quadquadeng=quadquadeng*onethird*onethird
c
c accumlate total real-space energies:
c
            eng=eng+chgquadeng+quadquadeng-quaddipeng 
c
c update accumlators:
c
            qquadacc=qquadacc+chgquadeng
            dipquadacc=dipquadacc-quaddipeng
            quadquadacc=quadquadacc+quadquadeng
c
c Quadrupole - ? interaction tensor.
c (21 in all).
c (You are really not going to believe this !)
c Txxxxx,Tyyyyy,Tzzzzz.
c
            txxxxx=-945.0d0*dx*dx*dx*dx*dx*dr11rec
     x            +1050d0*dx*dx*dx*dr9rec
     x            -225.0d0*dx*dr7rec
            tyyyyy=-945.0d0*dy*dy*dy*dy*dy*dr11rec
     x            +1050d0*dy*dy*dy*dr9rec
     x            -225.0d0*dy*dr7rec
            tzzzzz=-945.0d0*dz*dz*dz*dz*dz*dr11rec
     x            +1050d0*dz*dz*dz*dr9rec
     x            -225.0d0*dz*dr7rec
c
c Txxxxy,Txxxxz,Txyyyy,Txzzzz,Tyyyyz,Tyzzzz.
c
            txxxxy=-945.0d0*dx*dx*dx*dx*dy*dr11rec
     x            -45.0d0*dy*dr7rec
     x            +630.0d0*dx*dx*dy*dr9rec
            txxxxz=-945.0d0*dx*dx*dx*dx*dz*dr11rec
     x            -45.0d0*dz*dr7rec
     x            +630.0d0*dx*dx*dz*dr9rec
            tyyyyz=-945.0d0*dy*dy*dy*dy*dz*dr11rec
     x            -45.0d0*dz*dr7rec
     x            +630.0d0*dy*dy*dz*dr9rec
            txyyyy=-945.0d0*dx*dy*dy*dy*dy*dr11rec
     x            -45.0d0*dx*dr7rec
     x            +630.0d0*dx*dy*dy*dr9rec
            txzzzz=-945.0d0*dx*dz*dz*dz*dz*dr11rec
     x            -45.0d0*dx*dr7rec
     x            +630.0d0*dx*dz*dz*dr9rec
            tyzzzz=-945.0d0*dy*dz*dz*dz*dz*dr11rec
     x            -45.0d0*dy*dr7rec
     x            +630.0d0*dy*dz*dz*dr9rec
c
c Txxxyy,Txxxzz,Txxyyy,Txxzzz,Tyyyzz,Tyyzzz.
c
            txxxyy=-945.0d0*dx*dx*dx*dy*dy*dr11rec
     x            +105.0d0*dx*dx*dx*dr9rec
     x            +315.0d0*dx*dy*dy*dr9rec
     x            -45.0d0*dx*dr7rec
            txxxzz=-945.0d0*dx*dx*dx*dz*dz*dr11rec
     x            +105.0d0*dx*dx*dx*dr9rec
     x            +315.0d0*dx*dz*dz*dr9rec
     x            -45.0d0*dx*dr7rec
            tyyyzz=-945.0d0*dy*dy*dy*dz*dz*dr11rec
     x            +105.0d0*dy*dy*dy*dr9rec
     x            +315.0d0*dy*dz*dz*dr9rec
     x            -45.0d0*dy*dr7rec
            txxyyy=-945.0d0*dx*dx*dy*dy*dy*dr11rec
     x            +105.0d0*dy*dy*dy*dr9rec
     x            +315.0d0*dx*dx*dy*dr9rec
     x            -45.0d0*dy*dr7rec
            txxzzz=-945.0d0*dx*dx*dz*dz*dz*dr11rec
     x            +105.0d0*dz*dz*dz*dr9rec
     x            +315.0d0*dx*dx*dz*dr9rec
     x            -45.0d0*dz*dr7rec
            tyyzzz=-945.0d0*dy*dy*dz*dz*dz*dr11rec
     x            +105.0d0*dz*dz*dz*dr9rec
     x            +315.0d0*dy*dy*dz*dr9rec
     x            -45.0d0*dz*dr7rec
c
c Txxxyz,Txyyyz,Txyzzz.
c
            txxxyz=-945.0d0*dx*dx*dx*dy*dz*dr11rec
     x            +315.0d0*dx*dy*dz*dr9rec
            txyyyz=-945.0d0*dy*dy*dy*dx*dz*dr11rec
     x            +315.0d0*dx*dy*dz*dr9rec
            txyzzz=-945.0d0*dz*dz*dz*dx*dy*dr11rec
     x            +315.0d0*dx*dy*dz*dr9rec
c
c Txxyyz,Txxyzz,Txyyzz.
c
            txxyyz=-945.0d0*dx*dx*dy*dy*dz*dr11rec
     x            -15.0d0*dz*dr7rec
     x            +105.0d0*dz*dr9rec*((dx*dx)+(dy*dy))
            txxyzz=-945.0d0*dx*dx*dy*dz*dz*dr11rec
     x            -15.0d0*dy*dr7rec
     x            +105.0d0*dy*dr9rec*((dx*dx)+(dz*dz))
            txyyzz=-945.0d0*dx*dy*dy*dz*dz*dr11rec
     x            -15.0d0*dx*dr7rec
     x            +105.0d0*dx*dr9rec*((dy*dy)+(dz*dz))
c
c Quadrupole-quadrupole forces.
c xx
c
            fquadquadxxbyrx=quadxx(i)*(txxxxx*quadxx(j)
     x                     +(2.0d0*txxxxy*quadxy(j))
     x                     +(2.0d0*txxxxz*quadxz(j))
     x                     +(2.0d0*txxxyz*quadyz(j))
     x                     +txxxyy*quadyy(j)
     x                     +txxxzz*quadzz(j))
            fquadquadxxbyry=quadxx(i)*(txxxxy*quadxx(j)
     x                     +(2.0d0*txxxyy*quadxy(j))
     x                     +(2.0d0*txxxyz*quadxz(j))
     x                     +(2.0d0*txxyyz*quadyz(j))
     x                     +txxyyy*quadyy(j)
     x                     +txxyzz*quadzz(j))
            fquadquadxxbyrz=quadxx(i)*(txxxxz*quadxx(j)
     x                     +(2.0d0*txxxyz*quadxy(j))
     x                     +(2.0d0*txxxzz*quadxz(j))
     x                     +(2.0d0*txxyzz*quadyz(j))
     x                     +txxyyz*quadyy(j)
     x                     +txxzzz*quadzz(j))
c
c yy.
c
            fquadquadyybyrx=quadyy(i)*(txxxyy*quadxx(j)
     x                     +(2.0d0*txxyyy*quadxy(j))
     x                     +(2.0d0*txxyyz*quadxz(j))
     x                     +(2.0d0*txyyyz*quadyz(j))
     x                     +txyyyy*quadyy(j)
     x                     +txyyzz*quadzz(j))
            fquadquadyybyry=quadyy(i)*(txxyyy*quadxx(j)
     x                     +(2.0d0*txyyyy*quadxy(j))
     x                     +(2.0d0*txyyyz*quadxz(j))
     x                     +(2.0d0*tyyyyz*quadyz(j))
     x                     +tyyyyy*quadyy(j)
     x                     +tyyyzz*quadzz(j))
            fquadquadyybyrz=quadyy(i)*(txxyyz*quadxx(j)
     x                     +(2.0d0*txyyyz*quadxy(j))
     x                     +(2.0d0*txyyzz*quadxz(j))
     x                     +(2.0d0*tyyyzz*quadyz(j))
     x                     +tyyyyz*quadyy(j)
     x                     +tyyzzz*quadzz(j))
c
c zz.
c
            fquadquadzzbyrx=quadzz(i)*(txxxzz*quadxx(j)
     x                     +(2.0d0*txxyzz*quadxy(j))
     x                     +(2.0d0*txxzzz*quadxz(j))
     x                     +(2.0d0*txyzzz*quadyz(j))
     x                     +txyyzz*quadyy(j)
     x                     +txzzzz*quadzz(j))
            fquadquadzzbyry=quadzz(i)*(txxyzz*quadxx(j)
     x                     +(2.0d0*txyyzz*quadxy(j))
     x                     +(2.0d0*txyzzz*quadxz(j))
     x                     +(2.0d0*tyyzzz*quadyz(j))
     x                     +tyyyzz*quadyy(j)
     x                     +tyzzzz*quadzz(j))
            fquadquadzzbyrz=quadzz(i)*(txxzzz*quadxx(j)
     x                     +(2.0d0*txyzzz*quadxy(j))
     x                     +(2.0d0*txzzzz*quadxz(j))
     x                     +(2.0d0*tyzzzz*quadyz(j))
     x                     +tyyzzz*quadyy(j)
     x                     +tzzzzz*quadzz(j))
c
c xy.
c
            fquadquadxybyrx=quadxy(i)*((2.0d0*txxxxy*quadxx(j))
     x                    +(4.0d0*txxxyy*quadxy(j))
     x                    +(4.0d0*txxxyz*quadxz(j))
     x                    +(2.0d0*txxyyy*quadyy(j))
     x                    +(4.0d0*txxyyz*quadyz(j))
     x                    +(2.0d0*txxyzz*quadzz(j)))
            fquadquadxybyry=quadxy(i)*((2.0d0*txxxyy*quadxx(j))
     x                    +(4.0d0*txxyyy*quadxy(j))
     x                    +(4.0d0*txxyyz*quadxz(j))
     x                    +(2.0d0*txyyyy*quadyy(j))
     x                    +(4.0d0*txyyyz*quadyz(j))
     x                    +(2.0d0*txyyzz*quadzz(j)))
            fquadquadxybyrz=quadxy(i)*((2.0d0*txxxyz*quadxx(j))
     x                    +(4.0d0*txxyyz*quadxy(j))
     x                    +(4.0d0*txxyzz*quadxz(j))
     x                    +(2.0d0*txyyyz*quadyy(j))
     x                    +(4.0d0*txyyzz*quadyz(j))
     x                    +(2.0d0*txyzzz*quadzz(j)))
c
c xz.
c
            fquadquadxzbyrx=quadxz(i)*((2.0d0*txxxxz*quadxx(j))
     x                    +(4.0d0*txxxyz*quadxy(j))
     x                    +(4.0d0*txxxzz*quadxz(j))
     x                    +(2.0d0*txxyyz*quadyy(j))
     x                    +(4.0d0*txxyzz*quadyz(j))
     x                    +(2.0d0*txxzzz*quadzz(j)))
            fquadquadxzbyry=quadxz(i)*((2.0d0*txxxyz*quadxx(j))
     x                    +(4.0d0*txxyyz*quadxy(j))
     x                    +(4.0d0*txxyzz*quadxz(j))
     x                    +(2.0d0*txyyyz*quadyy(j))
     x                    +(4.0d0*txyyzz*quadyz(j))
     x                    +(2.0d0*txyzzz*quadzz(j)))
            fquadquadxzbyrz=quadxz(i)*((2.0d0*txxxzz*quadxx(j))
     x                    +(4.0d0*txxyzz*quadxy(j))
     x                    +(4.0d0*txxzzz*quadxz(j))
     x                    +(2.0d0*txyyzz*quadyy(j))
     x                    +(4.0d0*txyzzz*quadyz(j))
     x                    +(2.0d0*txzzzz*quadzz(j)))
c
c yz.
c
            fquadquadyzbyrx=quadyz(i)*((2.0d0*txxxyz*quadxx(j))
     x                    +(4.0d0*txxyyz*quadxy(j))
     x                    +(4.0d0*txxyzz*quadxz(j))
     x                    +(2.0d0*txyyyz*quadyy(j))
     x                    +(4.0d0*txyyzz*quadyz(j))
     x                    +(2.0d0*txyzzz*quadzz(j)))
            fquadquadyzbyry=quadyz(i)*((2.0d0*txxyyz*quadxx(j))
     x                    +(4.0d0*txyyyz*quadxy(j))
     x                    +(4.0d0*txyyzz*quadxz(j))
     x                    +(2.0d0*tyyyyz*quadyy(j))
     x                    +(4.0d0*tyyyzz*quadyz(j))
     x                    +(2.0d0*tyyzzz*quadzz(j)))
            fquadquadyzbyrz=quadyz(i)*((2.0d0*txxyzz*quadxx(j))
     x                    +(4.0d0*txyyzz*quadxy(j))
     x                    +(4.0d0*txyzzz*quadxz(j))
     x                    +(2.0d0*tyyyzz*quadyy(j))
     x                    +(4.0d0*tyyzzz*quadyz(j))
     x                    +(2.0d0*tyzzzz*quadzz(j)))

            fquadquadx=fquadquadxxbyrx+fquadquadyybyrx
     x                +fquadquadzzbyrx+fquadquadxybyrx
     x                +fquadquadxzbyrx+fquadquadyzbyrx
            fquadquady=fquadquadxxbyry+fquadquadyybyry
     x                +fquadquadzzbyry+fquadquadxybyry
     x                +fquadquadxzbyry+fquadquadyzbyry
            fquadquadz=fquadquadxxbyrz+fquadquadyybyrz
     x                +fquadquadzzbyrz+fquadquadxybyrz
     x                +fquadquadxzbyrz+fquadquadyzbyrz

            fquadquadx=fquadquadx*onethird*onethird
            fquadquady=fquadquady*onethird*onethird
            fquadquadz=fquadquadz*onethird*onethird

            stpquadquadxx=stpquadquadxx+fquadquadx*dxsav(i,j)
            stpquadquadxy=stpquadquadxy+0.5d0*(fquadquadx*dysav(i,j)
     x                        +fquadquady*dxsav(i,j))
            stpquadquadxz=stpquadquadxz+0.5d0*(fquadquadx*dzsav(i,j)
     x                        +fquadquadz*dxsav(i,j))
            stpquadquadyy=stpquadquadyy+fquadquady*dysav(i,j)
            stpquadquadyz=stpquadquadyz+0.5d0*(fquadquadz*dysav(i,j)
     x                        +fquadquady*dzsav(i,j))
            stpquadquadzz=stpquadquadzz+fquadquadz*dzsav(i,j)

c
c Sum together all forces.
c
            fttx=fqquadx-fdipquadx-fquadquadx
            ftty=fqquady-fdipquady-fquadquady
            fttz=fqquadz-fdipquadz-fquadquadz

            frrx(j)=frrx(j)-fttx
            frry(j)=frry(j)-ftty
            frrz(j)=frrz(j)-fttz

            frrx(i)=frrx(i)+fttx
            frry(i)=frry(i)+ftty
            frrz(i)=frrz(i)+fttz



c*************************************************
c*************************************************
c*************************************************
c*************************************************
c*************************************************
c*************************************************

c damping (enhancement) moved here 17.4.2000
c TMP i-j term - i is defined as at the centre of the cluster.
cc            dampa1=dampa(jpoint,ipoint)
cc            dampa2=dampa1*dampa1/2.0d0
cc            dampa3=dampa1*dampa1*dampa1/6.0d0
cc            dampa4=dampa1*dampa1*dampa1*dampa1/24.0d0

cc            dampsumfi=1+(dampa1*dr)+(dampa2*drsq)+(dampa3*drsq*dr)
cc     x             +(dampa4*drsq*drsq)

            dampsumfi=1.0d0
            xf=1.0d0
            factorial=1.0d0
            do 18 kk=1,nkdamp(jpoint,ipoint)
               xf=(xf*dampa(jpoint,ipoint)*dr)
               factorial=factorial*float(kk)
               dampsumfi=dampsumfi+(xf/factorial)
   18       continue

            dampaexpi=dexp(-dampa(jpoint,ipoint)*dr)
            dampfunci=-dampsumfi*dampaexpi*dampfac(jpoint,ipoint)
cc            dampfuncdiffi=dampa(jpoint,ipoint)*dampa4*drsq*drsq*dampaexpi
            dampfuncdiffi=dampa(jpoint,ipoint)*dampaexpi
     x             *dampfac(jpoint,ipoint)
     x             *((dampa(jpoint,ipoint)*dr)**nkdamp(jpoint,ipoint))
     x             /(factorial)

            r3dampi=dr3rec*dampfunci*q(j)

            elecsrxi=dxsav(i,j)*r3dampi
            elecsryi=dysav(i,j)*r3dampi
            elecsrzi=dzsav(i,j)*r3dampi

            elecxsr(i)=elecxsr(i)+elecsrxi
            elecysr(i)=elecysr(i)+elecsryi
            eleczsr(i)=eleczsr(i)+elecsrzi

c TMP j-i term - j now at the centre of the cluster.
cc            dampa1=dampa(ipoint,jpoint)
cc            dampa2=dampa1*dampa1/2.0d0
cc            dampa3=dampa1*dampa1*dampa1/6.0d0
cc            dampa4=dampa1*dampa1*dampa1*dampa1/24.0d0

cc            dampsumfj=1+(dampa1*dr)+(dampa2*drsq)+(dampa3*drsq*dr)
cc     x             +(dampa4*drsq*drsq)

            dampsumfj=1.0d0
            xf=1.0d0
            factorial=1.0d0
            do 16 kk=1,nkdamp(ipoint,jpoint)
               xf=(xf*dampa(ipoint,jpoint)*dr)
               factorial=factorial*float(kk)
               dampsumfj=dampsumfj+(xf/factorial)
   16       continue

            dampaexpj=dexp(-dampa(ipoint,jpoint)*dr)
            dampfuncj=-dampsumfj*dampaexpj*dampfac(ipoint,jpoint)
cc            dampfuncdiffj=dampa1*dampa4*drsq*drsq*dampaexpj
            dampfuncdiffj=dampa(ipoint,jpoint)*dampaexpj
     x             *dampfac(ipoint,jpoint)
     x             *((dampa(ipoint,jpoint)*dr)**nkdamp(ipoint,jpoint))
     x             /(factorial)

            r3dampj=dr3rec*dampfuncj*q(i)

            elecsrxj=-dxsav(i,j)*r3dampj
            elecsryj=-dysav(i,j)*r3dampj
            elecsrzj=-dzsav(i,j)*r3dampj

            elecxsr(j)=elecxsr(j)+elecsrxj
            elecysr(j)=elecysr(j)+elecsryj
            eleczsr(j)=eleczsr(j)+elecsrzj

c
c Calculate contribution to ion forces:
c
            xmuidotrij=xmu(i)*dxsav(i,j)
     x                +ymu(i)*dysav(i,j)
     x                +zmu(i)*dzsav(i,j)
            xmujdotrij=xmu(j)*dxsav(i,j)
     x                +ymu(j)*dysav(i,j)
     x                +zmu(j)*dzsav(i,j)

            xmuidotT=(3.0d0*dxsav(i,j)*xmuidotrij
     x              -drsq*xmu(i))*dr5rec
            ymuidotT=(3.0d0*dysav(i,j)*xmuidotrij
     x              -drsq*ymu(i))*dr5rec
            zmuidotT=(3.0d0*dzsav(i,j)*xmuidotrij
     x              -drsq*zmu(i))*dr5rec

            xmujdotT=(3.0d0*dxsav(i,j)*xmujdotrij
     x              -drsq*xmu(j))*dr5rec
            ymujdotT=(3.0d0*dysav(i,j)*xmujdotrij
     x              -drsq*ymu(j))*dr5rec
            zmujdotT=(3.0d0*dzsav(i,j)*xmujdotrij
     x              -drsq*zmu(j))*dr5rec

c
c Calculate the short range energy contribution for the dipoles:
c
            srdipeng(ipoint)=srdipeng(ipoint)
     x                         -(elecsrxi*xmu(i)
     x                          +elecsryi*ymu(i)
     x                          +elecsrzi*zmu(i))
            srdipeng(jpoint)=srdipeng(jpoint)
     x                         -(elecsrxj*xmu(j)
     x                          +elecsryj*ymu(j)
     x                          +elecsrzj*zmu(j))

            dipi=xmu(i)*dxsav(i,j)
     x         +ymu(i)*dysav(i,j)+zmu(i)*dzsav(i,j)
            dipj=xmu(j)*dxsav(i,j)
     x         +ymu(j)*dysav(i,j)+zmu(j)*dzsav(i,j)

            diprfunci=dipi*drrec*q(j)
     x              *((-3.0d0*dampfunci*dr3rec*drrec)
     x              +(dampfuncdiffi*dr3rec))
            diprfuncj=dipj*drrec*q(i)
     x              *((-3.0d0*dampfuncj*dr3rec*drrec)
     x              +(dampfuncdiffj*dr3rec))

            frrxsri=(xmu(i)*r3dampi)+(dxsav(i,j)*diprfunci)
            frrysri=(ymu(i)*r3dampi)+(dysav(i,j)*diprfunci)
            frrzsri=(zmu(i)*r3dampi)+(dzsav(i,j)*diprfunci)
            frrxsrj=(xmu(j)*r3dampj)+(dxsav(i,j)*diprfuncj)
            frrysrj=(ymu(j)*r3dampj)+(dysav(i,j)*diprfuncj)
            frrzsrj=(zmu(j)*r3dampj)+(dzsav(i,j)*diprfuncj)

c
c Add these other s-r terms to the virial:
c
c need to be checked think 1st set right.
c Also need to add dip-dip contribution somewhere.

cc            xviraccsrdip=xviraccsrdip+((frrxsr)*dxsav(i,j))
cc            yviraccsrdip=yviraccsrdip+((frrysr)*dysav(i,j))
cc            zviraccsrdip=zviraccsrdip+((frrzsr)*dzsav(i,j))

c
c Stress tensor accumulation:
c charge-dipole damping term:
c TMP - COME BACK TO THIS - MISSING TERMS....
            stpsrxx=stpsrxx+frrxsri*dxsav(i,j)
     x                     -frrxsrj*dxsav(i,j)
            stpsrxy=stpsrxy+0.5d0*(frrxsri*dysav(i,j)
     x                            +frrysri*dxsav(i,j))
     x                     -0.5d0*(frrxsrj*dysav(i,j)
     x                            +frrysrj*dxsav(i,j))
            stpsrxz=stpsrxz+0.5d0*(frrxsri*dzsav(i,j)
     x                            +frrzsri*dxsav(i,j))
     x                     -0.5d0*(frrxsrj*dzsav(i,j)
     x                            +frrzsrj*dxsav(i,j))
            stpsryy=stpsryy+frrysri*dysav(i,j)
     x                     -frrysrj*dysav(i,j)
            stpsryz=stpsryz+0.5d0*(frrysri*dzsav(i,j)
     x                            +frrzsri*dysav(i,j))
     x                     -0.5d0*(frrysrj*dzsav(i,j)
     x                            +frrzsrj*dysav(i,j))
            stpsrzz=stpsrzz+frrzsri*dzsav(i,j)
     x                     -frrzsrj*dzsav(i,j)

c
c Accumulate forces:
c
            frrx(i)=frrx(i)+frrxsri-frrxsrj
            frry(i)=frry(i)+frrysri-frrysrj
            frrz(i)=frrz(i)+frrzsri-frrzsrj

            frrx(j)=frrx(j)-frrxsri+frrxsrj
            frry(j)=frry(j)-frrysri+frrysrj
            frrz(j)=frrz(j)-frrzsri+frrzsrj

c Field gradients. (from AJR via MW)
c 1) Permanent charge contribution -T_2q.
c

            dxx=dxsav(i,j)*dxsav(i,j)
            dyy=dysav(i,j)*dysav(i,j)
            dzz=dzsav(i,j)*dzsav(i,j)
            dxy=dxsav(i,j)*dysav(i,j)
            dxz=dxsav(i,j)*dzsav(i,j)
            dyz=dysav(i,j)*dzsav(i,j)

            txxt=txx*efact+efact2*dxx
            tyyt=tyy*efact+efact2*dyy
            tzzt=tzz*efact+efact2*dzz
            txyt=txy*efact+efact2*dxy
            txzt=txz*efact+efact2*dxz
            tyzt=tyz*efact+efact2*dyz

            exx(j)=exx(j)-txxt*q(i)
            eyy(j)=eyy(j)-tyyt*q(i)
            ezz(j)=ezz(j)-tzzt*q(i)
            exy(j)=exy(j)-txyt*q(i)
            exz(j)=exz(j)-txzt*q(i)
            eyz(j)=eyz(j)-tyzt*q(i)

            exx(i)=exx(i)-txxt*q(j)
            eyy(i)=eyy(i)-tyyt*q(j)
            ezz(i)=ezz(i)-tzzt*q(j)
            exy(i)=exy(i)-txyt*q(j)
            exz(i)=exz(i)-txzt*q(j)
            eyz(i)=eyz(i)-tyzt*q(j)

c
c 2) Dipole contribution +T_3 mu.
c
            exxdipi=txxx*xmu(i)
     x             +txxy*ymu(i)
     x             +txxz*zmu(i)
            eyydipi=tyyy*ymu(i)
     x             +txyy*xmu(i)
     x             +tyyz*zmu(i)
            ezzdipi=tzzz*zmu(i)
     x             +tyzz*ymu(i)
     x             +txzz*xmu(i)
            exydipi=txxy*xmu(i)
     x             +txyy*ymu(i)
     x             +txyz*zmu(i)
            exzdipi=txxz*xmu(i)
     x             +txyz*ymu(i)
     x             +txzz*zmu(i)
            eyzdipi=txyz*xmu(i)
     x             +tyyz*ymu(i)
     x             +tyzz*zmu(i)

            exx(j)=exx(j)+exxdipi
            eyy(j)=eyy(j)+eyydipi
            ezz(j)=ezz(j)+ezzdipi
            exy(j)=exy(j)+exydipi
            exz(j)=exz(j)+exzdipi
            eyz(j)=eyz(j)+eyzdipi

            exxdipj=-txxx*xmu(j)
     x             -txxy*ymu(j)
     x             -txxz*zmu(j)
            eyydipj=-tyyy*ymu(j)
     x             -txyy*xmu(j)
     x             -tyyz*zmu(j)
            ezzdipj=-tzzz*zmu(j)
     x             -tyzz*ymu(j)
     x             -txzz*xmu(j)
            exydipj=-txxy*xmu(j)
     x             -txyy*ymu(j)
     x             -txyz*zmu(j)
            exzdipj=-txxz*xmu(j)
     x             -txyz*ymu(j)
     x             -txzz*zmu(j)
            eyzdipj=-txyz*xmu(j)
     x             -tyyz*ymu(j)
     x             -tyzz*zmu(j)

            exx(i)=exx(i)+exxdipj
            eyy(i)=eyy(i)+eyydipj
            ezz(i)=ezz(i)+ezzdipj
            exy(i)=exy(i)+exydipj
            exz(i)=exz(i)+exzdipj
            eyz(i)=eyz(i)+eyzdipj

c**********************************************
c
c 3) Quadrupole contribution - T_4 quad.
c
            exxquadi=-txxxx*quadxx(i)
     x               -2.0d0*txxxy*quadxy(i)
     x               -2.0d0*txxxz*quadxz(i)
     x               -txxyy*quadyy(i)
     x               -txxzz*quadzz(i)
     x               -2.0d0*txxyz*quadyz(i)
            eyyquadi=-txxyy*quadxx(i)
     x               -2.0d0*txyyy*quadxy(i)
     x               -2.0d0*txyyz*quadxz(i)
     x               -tyyyy*quadyy(i)
     x               -tyyzz*quadzz(i)
     x               -2.0d0*tyyyz*quadyz(i)
            ezzquadi=-txxzz*quadxx(i)
     x               -2.0d0*txyzz*quadxy(i)
     x               -2.0d0*txzzz*quadxz(i)
     x               -tyyzz*quadyy(i)
     x               -tzzzz*quadzz(i)
     x               -2.0d0*tyzzz*quadyz(i)
            exyquadi=-txxxy*quadxx(i)
     x               -2.0d0*txxyy*quadxy(i)
     x               -2.0d0*txxyz*quadxz(i)
     x               -txyyy*quadyy(i)
     x               -txyzz*quadzz(i)
     x               -2.0d0*txyyz*quadyz(i)
            exzquadi=-txxxz*quadxx(i)
     x               -2.0d0*txxyz*quadxy(i)
     x               -2.0d0*txxzz*quadxz(i)
     x               -txyyz*quadyy(i)
     x               -txzzz*quadzz(i)
     x               -2.0d0*txyzz*quadyz(i)
            eyzquadi=-txxyz*quadxx(i)
     x               -2.0d0*txyyz*quadxy(i)
     x               -2.0d0*txyzz*quadxz(i)
     x               -tyyyz*quadyy(i)
     x               -tyzzz*quadzz(i)
     x               -2.0d0*tyyzz*quadyz(i)

            exx(j)=exx(j)+exxquadi*onethird
            eyy(j)=eyy(j)+eyyquadi*onethird
            ezz(j)=ezz(j)+ezzquadi*onethird
            exy(j)=exy(j)+exyquadi*onethird
            exz(j)=exz(j)+exzquadi*onethird
            eyz(j)=eyz(j)+eyzquadi*onethird

            exxquadj=-txxxx*quadxx(j)
     x               -2.0d0*txxxy*quadxy(j)
     x               -2.0d0*txxxz*quadxz(j)
     x               -txxyy*quadyy(j)
     x               -txxzz*quadzz(j)
     x               -2.0d0*txxyz*quadyz(j)
            eyyquadj=-txxyy*quadxx(j)
     x               -2.0d0*txyyy*quadxy(j)
     x               -2.0d0*txyyz*quadxz(j)
     x               -tyyyy*quadyy(j)
     x               -tyyzz*quadzz(j)
     x               -2.0d0*tyyyz*quadyz(j)
            ezzquadj=-txxzz*quadxx(j)
     x               -2.0d0*txyzz*quadxy(j)
     x               -2.0d0*txzzz*quadxz(j)
     x               -tyyzz*quadyy(j)
     x               -tzzzz*quadzz(j)
     x               -2.0d0*tyzzz*quadyz(j)
            exyquadj=-txxxy*quadxx(j)
     x               -2.0d0*txxyy*quadxy(j)
     x               -2.0d0*txxyz*quadxz(j)
     x               -txyyy*quadyy(j)
     x               -txyzz*quadzz(j)
     x               -2.0d0*txyyz*quadyz(j)
            exzquadj=-txxxz*quadxx(j)
     x               -2.0d0*txxyz*quadxy(j)
     x               -2.0d0*txxzz*quadxz(j)
     x               -txyyz*quadyy(j)
     x               -txzzz*quadzz(j)
     x               -2.0d0*txyzz*quadyz(j)
            eyzquadj=-txxyz*quadxx(j)
     x               -2.0d0*txyyz*quadxy(j)
     x               -2.0d0*txyzz*quadxz(j)
     x               -tyyyz*quadyy(j)
     x               -tyzzz*quadzz(j)
     x               -2.0d0*tyyzz*quadyz(j)

            exx(i)=exx(i)+exxquadj*onethird
            eyy(i)=eyy(i)+eyyquadj*onethird
            ezz(i)=ezz(i)+ezzquadj*onethird
            exy(i)=exy(i)+exyquadj*onethird
            exz(i)=exz(i)+exzquadj*onethird
            eyz(i)=eyz(i)+eyzquadj*onethird

c**********************************************
c Field grad. SR contribution.
cc            dampa1=fgb(jpoint,ipoint)
cc            dampa2=dampa1*dampa1/2.0d0
cc            dampa3=dampa1*dampa1*dampa1/6.0d0
cc            dampa4=dampa1*dampa1*dampa1*dampa1/24.0d0

cc            dampsumfi=1+(dampa1*dr)+(dampa2*drsq)+(dampa3*drsq*dr)
cc     x             +(dampa4*drsq*drsq)

cc            dampaexpi=dexp(-dampa1*dr)
cc            dampfunci=-dampsumfi*dampaexpi*fgc(jpoint,ipoint)
cc            dampfuncdiffi=dampa1*dampa4*drsq*drsq*dampaexpi

            dampsumfi=1.0d0
            xf=1.0d0
            factorial=1.0d0
            do 15 kk=1,nkfg(jpoint,ipoint)
               xf=(xf*fgb(jpoint,ipoint)*dr)
               factorial=factorial*float(kk)
               dampsumfi=dampsumfi+(xf/factorial)
   15       continue

            dampaexpi=dexp(-fgb(jpoint,ipoint)*dr)
            dampfunci=-dampsumfi*dampaexpi*fgc(jpoint,ipoint)
cc            dampfuncdiffj=dampa1*dampa4*drsq*drsq*dampaexpj
            dampfuncdiffi=fgb(jpoint,ipoint)*dampaexpi
     x             *((fgb(jpoint,ipoint)*dr)**nkfg(jpoint,ipoint))
     x             /(factorial)


c short-range damping contribution to electric fields.
c
            exxsrg=txx*dampfunci*q(j)
            eyysrg=tyy*dampfunci*q(j)
            ezzsrg=tzz*dampfunci*q(j)
            exysrg=txy*dampfunci*q(j)
            exzsrg=txz*dampfunci*q(j)
            eyzsrg=tyz*dampfunci*q(j)

            exxsr(i)=exxsr(i)-exxsrg
            eyysr(i)=eyysr(i)-eyysrg
            ezzsr(i)=ezzsr(i)-ezzsrg
            exysr(i)=exysr(i)-exysrg
            exzsr(i)=exzsr(i)-exzsrg
            eyzsr(i)=eyzsr(i)-eyzsrg

            srquadeng(ipoint)=srquadeng(ipoint)
     x                       +onethird*(exxsrg*quadxx(i)
     x                                  +eyysrg*quadyy(i)
     x                                  +ezzsrg*quadzz(i)
     x                           +2.0d0*(exysrg*quadxy(i)
     x                                  +exzsrg*quadxz(i)
     x                                  +eyzsrg*quadyz(i)))

c TMP j-i term - j now at the centre of the cluster.
cc            dampa1=fgb(ipoint,jpoint)
cc            dampa2=dampa1*dampa1/2.0d0
cc            dampa3=dampa1*dampa1*dampa1/6.0d0
cc            dampa4=dampa1*dampa1*dampa1*dampa1/24.0d0

cc            dampsumfj=1+(dampa1*dr)+(dampa2*drsq)+(dampa3*drsq*dr)
cc    x             +(dampa4*drsq*drsq)

cc            dampaexpj=dexp(-dampa1*dr)
cc            dampfuncj=-dampsumfj*dampaexpj*fgc(ipoint,jpoint)
cc            dampfuncdiffj=dampa1*dampa4*drsq*drsq*dampaexpj

            dampsumfj=1.0d0
            xf=1.0d0
            factorial=1.0d0
            do 17 kk=1,nkfg(ipoint,jpoint)
               xf=(xf*fgb(ipoint,jpoint)*dr)
               factorial=factorial*float(kk)
               dampsumfj=dampsumfj+(xf/factorial)
   17       continue

            dampaexpj=dexp(-fgb(ipoint,jpoint)*dr)
            dampfuncj=-dampsumfj*dampaexpj*fgc(ipoint,jpoint)
            dampfuncdiffj=fgb(ipoint,jpoint)*dampaexpj
     x             *((fgb(ipoint,jpoint)*dr)**nkfg(ipoint,jpoint))
     x             /(factorial)

            exxsrg=txx*dampfuncj*q(i)
            eyysrg=tyy*dampfuncj*q(i)
            ezzsrg=tzz*dampfuncj*q(i)
            exysrg=txy*dampfuncj*q(i)
            exzsrg=txz*dampfuncj*q(i)
            eyzsrg=tyz*dampfuncj*q(i)

            exxsr(j)=exxsr(j)-exxsrg
            eyysr(j)=eyysr(j)-eyysrg
            ezzsr(j)=ezzsr(j)-ezzsrg
            exysr(j)=exysr(j)-exysrg
            exzsr(j)=exzsr(j)-exzsrg
            eyzsr(j)=eyzsr(j)-eyzsrg

            srquadeng(jpoint)=srquadeng(jpoint)
     x                       +onethird*(exxsrg*quadxx(j)
     x                                  +eyysrg*quadyy(j)
     x                                  +ezzsrg*quadzz(j)
     x                           +2.0d0*(exysrg*quadxy(j)
     x                                  +exzsrg*quadxz(j)
     x                                  +eyzsrg*quadyz(j)))

            T2dotquadj=quadxx(j)*txx+quadyy(j)*tyy
     x               +quadzz(j)*tzz+2.0d0*(quadxy(j)*txy
     x               +quadxz(j)*txz+quadyz(j)*tyz)
            T2dotquadi=quadxx(i)*txx+quadyy(i)*tyy
     x                +quadzz(i)*tzz+2.0d0*(quadxy(i)*txy
     x                +quadxz(i)*txz+quadyz(i)*tyz)

            frrxsrqji=-((-txxx*quadxx(j)-txyy*quadyy(j)-txzz*quadzz(j)
     x               -2.0d0*txxy*quadxy(j)
     x               -2.0d0*txxz*quadxz(j)
     x               -2.0d0*txyz*quadyz(j))*dampfunc
     x               +(T2dotquadj*dampfuncdiffj*dx*drrec))
     x               *onethird*q(i)

            frrxsrqij=-((-txxx*quadxx(i)-txyy*quadyy(i)-txzz*quadzz(i)
     x               -2.0d0*txxy*quadxy(i)
     x               -2.0d0*txxz*quadxz(i)
     x               -2.0d0*txyz*quadyz(i))*dampfunc
     x               +(T2dotquadi*dampfuncdiffi*dx*drrec))
     x               *onethird*q(j)

            frrysrqji=-((-txxy*quadxx(j)-tyyy*quadyy(j)-tyzz*quadzz(j)
     x               -2.0d0*txyy*quadxy(j)
     x               -2.0d0*txyz*quadxz(j)
     x               -2.0d0*tyyz*quadyz(j))*dampfunc
     x               +(T2dotquadj*dampfuncdiffj*dy*drrec))
     x               *onethird*q(i)

            frrysrqij=-((-txxy*quadxx(i)-tyyy*quadyy(i)-tyzz*quadzz(i)
     x               -2.0d0*txyy*quadxy(i)
     x               -2.0d0*txyz*quadxz(i)
     x               -2.0d0*tyyz*quadyz(i))*dampfunc
     x               +(T2dotquadi*dampfuncdiffi*dy*drrec))
     x               *onethird*q(j)

            frrzsrqji=-((-txxz*quadxx(j)-tyyz*quadyy(j)-tzzz*quadzz(j)
     x               -2.0d0*txyz*quadxy(j)
     x               -2.0d0*txzz*quadxz(j)
     x               -2.0d0*tyzz*quadyz(j))*dampfunc
     x               +(T2dotquadj*dampfuncdiffj*dz*drrec))
     x               *onethird*q(i)

            frrzsrqij=-((-txxz*quadxx(i)-tyyz*quadyy(i)-tzzz*quadzz(i)
     x               -2.0d0*txyz*quadxy(i)
     x               -2.0d0*txzz*quadxz(i)
     x               -2.0d0*tyzz*quadyz(i))*dampfunc
     x               +(T2dotquadi*dampfuncdiffi*dz*drrec))
     x               *onethird*q(j)

cc         write(779,*)i,j,frrysrqji
cc     x     ,frrysrqijj,fqquady
            frrx(i)=frrx(i)+frrxsrqij+frrxsrqji
            frry(i)=frry(i)+frrysrqij+frrysrqji
            frrz(i)=frrz(i)+frrzsrqij+frrzsrqji

            frrx(j)=frrx(j)-frrxsrqij-frrxsrqji
            frry(j)=frry(j)-frrysrqij-frrysrqji
            frrz(j)=frrz(j)-frrzsrqij-frrzsrqji

   40    continue

   30 continue

      do 140 i=1,num

         ipoint=ntype(i)

         asdipx(i)=alppolar(ipoint)*elecx(i)
         asdipy(i)=alppolar(ipoint)*elecy(i)
         asdipz(i)=alppolar(ipoint)*elecz(i)

         srdipx(i)=alppolar(ipoint)*elecxsr(i)
         srdipy(i)=alppolar(ipoint)*elecysr(i)
         srdipz(i)=alppolar(ipoint)*eleczsr(i)

         elecxu(i)=elecx(i)
         elecyu(i)=elecy(i)
         eleczu(i)=elecz(i)

         elecx(i)=elecx(i)+elecxsr(i)
         elecy(i)=elecy(i)+elecysr(i)
         elecz(i)=elecz(i)+eleczsr(i)

         asquadxx(i)=Cpolar(ipoint)*exx(i)
         asquadyy(i)=Cpolar(ipoint)*eyy(i)
         asquadzz(i)=Cpolar(ipoint)*ezz(i)
         asquadxy(i)=Cpolar(ipoint)*exy(i)
         asquadxz(i)=Cpolar(ipoint)*exz(i)
         asquadyz(i)=Cpolar(ipoint)*eyz(i)

         srquadxx(i)=Cpolar(ipoint)*exxsr(i)
         srquadyy(i)=Cpolar(ipoint)*eyysr(i)
         srquadzz(i)=Cpolar(ipoint)*ezzsr(i)
         srquadxy(i)=Cpolar(ipoint)*exysr(i)
         srquadxz(i)=Cpolar(ipoint)*exzsr(i)
         srquadyz(i)=Cpolar(ipoint)*eyzsr(i)

         exx(i)=exx(i)+exxsr(i)
         eyy(i)=eyy(i)+eyysr(i)
         ezz(i)=ezz(i)+ezzsr(i)
         exy(i)=exy(i)+exysr(i)
         exz(i)=exz(i)+exzsr(i)
         eyz(i)=eyz(i)+eyzsr(i)

  140 continue

      write(970,*)real(qquadacc),real(qquadeng)
     x     ,real(dipquadacc),real(quadquadacc)

c
c energy accumulator for pressure calculation:
c multiply qdipacc by 2 as it is an r^-2 term
c similarly for dipdipacc (r^-3)
c toteng contains reciprocal space q-q energy.
c
      qqacc=erfcracc+toteng
      qdipacc=(qdipacc+qmueng)*2.0d0
      dipdipacc=(dipdipacc+xmumueng)*3.0d0
 
      qquadacc=(qquadacc+qquadeng)*3.0d0
      dipquadacc=dipquadacc*4.0d0
      quadquadacc=quadquadacc*5.0d0
    
      srvirft=xviraccftsrc+yviraccftsrc+zviraccftsrc
      srdipvir=xviraccsrdip+yviraccsrdip+zviraccsrdip
      virtot=srvirft+srdipvir

      virtotx=xviraccftsrc+xviraccsrdip
      virtoty=yviraccftsrc+yviraccsrdip
      virtotz=zviraccftsrc+zviraccsrdip

c..........estimate of the coulomb contribution to the pressure
      encoul=engpetot+erfcracc+qdipacc+dipdipacc 
c
c Final energy summation:
c
      srdipengtot=0.0d0
      do 963 i=1,nspec
         srdipengtot=srdipengtot+srdipeng(i)+srquadeng(i)
  963 continue

      engpetot=engpetot+qmueng+xmumueng+eng+srdipengtot+qquadeng

      write(980,*)exx(1),eyy(1),ezz(1)

      return
      end
