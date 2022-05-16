      subroutine realE_dippim

c************************************************************
c   Calculates the real energy and forces at each site.
c   The calculated energy is the real part of the Ewald
c  summation plus the short-range and Van der Waals
c  corrections from the Fumi-Tosi potential.
c
c  April 2000 - corrected Ewald.
c             - review loops to do multi-polarzable.
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

    5 continue
c
c Zero energy accumulators.
c
      eng=0.0d0

      do 6 i=1,nspec

         srdipeng(i)=0.0d0
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
cc     x          +dysav(i,j)*dysav(i,j)
cc     x          +dzsav(i,j)*dzsav(i,j)
c
c implement cut-off:
cc             drsav(i,j)=0.0d0
cc            drsav(j,i)=0.0d0
c            if (drsq.ge.rsqmax) goto 40

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

c damping (enhancement) moved here 17.4.2000
c TMP i-j term - i is defined as at the centre of the cluster.
            dampa1=dampa(jpoint,ipoint)
            dampa2=dampa1*dampa1/2.0d0
            dampa3=dampa1*dampa1*dampa1/6.0d0
            dampa4=dampa1*dampa1*dampa1*dampa1/24.0d0

            dampsumfi=1+(dampa1*dr)+(dampa2*drsq)+(dampa3*drsq*dr)
     x             +(dampa4*drsq*drsq)

            dampaexpi=dexp(-dampa1*dr)
            dampfunci=-dampsumfi*dampaexpi*dampfac(jpoint,ipoint)
            dampfuncdiffi=dampa1*dampa4*drsq*drsq*dampaexpi

            r3dampi=dr3rec*dampfunci*q(j)

            elecsrxi=dxsav(i,j)*r3dampi
            elecsryi=dysav(i,j)*r3dampi
            elecsrzi=dzsav(i,j)*r3dampi

            elecxsr(i)=elecxsr(i)+elecsrxi
            elecysr(i)=elecysr(i)+elecsryi
            eleczsr(i)=eleczsr(i)+elecsrzi

c TMP j-i term - j now at the centre of the cluster.
            dampa1=dampa(ipoint,jpoint)
            dampa2=dampa1*dampa1/2.0d0
            dampa3=dampa1*dampa1*dampa1/6.0d0
            dampa4=dampa1*dampa1*dampa1*dampa1/24.0d0

            dampsumfj=1+(dampa1*dr)+(dampa2*drsq)+(dampa3*drsq*dr)
     x             +(dampa4*drsq*drsq)

            dampaexpj=dexp(-dampa1*dr)
            dampfuncj=-dampsumfj*dampaexpj*dampfac(ipoint,jpoint)
            dampfuncdiffj=dampa1*dampa4*drsq*drsq*dampaexpj

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

  140 continue

c
c energy accumulator for pressure calculation:
c multiply qdipacc by 2 as it is an r^-2 term
c similarly for dipdipacc (r^-3)
c toteng contains reciprocal space q-q energy.
c
      qqacc=erfcracc+toteng
      qdipacc=(qdipacc+qmueng)*2.0d0
      dipdipacc=(dipdipacc+xmumueng)*3.0d0
 
    
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
         srdipengtot=srdipengtot+srdipeng(i)
  963 continue

      engpetot=engpetot+qmueng+xmumueng+eng+srdipengtot

      return
      end
