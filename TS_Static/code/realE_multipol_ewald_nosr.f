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

cc         write(6,*)i
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
            dx=hlab2(1,1)*dxcf
     x                +hlab2(1,2)*dycf+hlab2(1,3)*dzcf
            dy=hlab2(2,1)*dxcf
     x                +hlab2(2,2)*dycf+hlab2(2,3)*dzcf
            dz=hlab2(3,1)*dxcf
     x                +hlab2(3,2)*dycf+hlab2(3,3)*dzcf

            drsq=dx*dx
     x          +dy*dy
     x          +dz*dz
c
c implement cut-off:
cc            if (drsq.ge.rsqmax) goto 40

            drsqrec=1.0d0/drsq
            drsq3rec=drsqrec*drsqrec*drsqrec
            dr=dsqrt(drsq)
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

            fttx=sgam*dx
            ftty=sgam*dy
            fttz=sgam*dz

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

            ettx=egam*dx
            etty=egam*dy
            ettz=egam*dz

            elecx(j)=elecx(j)-(ettx*qjrec)
            elecy(j)=elecy(j)-(etty*qjrec)
            elecz(j)=elecz(j)-(ettz*qjrec)

            elecx(i)=elecx(i)+(ettx*qirec)
            elecy(i)=elecy(i)+(etty*qirec)
            elecz(i)=elecz(i)+(ettz*qirec)
c
c Start of dipole handling section:
c
            xmuidotrij=xmu(i)*dx
     x                +ymu(i)*dy
     x                +zmu(i)*dz
            xmujdotrij=xmu(j)*dx
     x                +ymu(j)*dy
     x                +zmu(j)*dz
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
     x                 -3.0d0*dx*drsqrec*chgdipengij)*efact
     x                 -dx*xmuidotrij*q(j)*efact2
            frrxdipqji=(qdotmuxji*dr3rec
     x                 -3.0d0*dx*drsqrec*chgdipengji)*efact
     x                 +dx*xmujdotrij*q(i)*efact2
            frrydipqij=(qdotmuyij*dr3rec
     x                 -3.0d0*dy*drsqrec*chgdipengij)*efact
     x                 -dy*xmuidotrij*q(j)*efact2
            frrydipqji=(qdotmuyji*dr3rec
     x                 -3.0d0*dy*drsqrec*chgdipengji)*efact
     x                 +dy*xmujdotrij*q(i)*efact2
            frrzdipqij=(qdotmuzij*dr3rec
     x                 -3.0d0*dz*drsqrec*chgdipengij)*efact
     x                 -dz*xmuidotrij*q(j)*efact2
            frrzdipqji=(qdotmuzji*dr3rec
     x                 -3.0d0*dz*drsqrec*chgdipengji)*efact
     x                 +dz*xmujdotrij*q(i)*efact2

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
            dxunitsq=dx*dx*drsqrec
            dyunitsq=dy*dy*drsqrec
            dzunitsq=dz*dz*drsqrec

            sxx=1.0d0-3.0d0*dxunitsq
            syy=1.0d0-3.0d0*dyunitsq
            szz=1.0d0-3.0d0*dzunitsq
            sxy=-3.0d0*dx*dy*drsqrec
            sxz=-3.0d0*dx*dz*drsqrec
            syz=-3.0d0*dz*dy*drsqrec

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
            txxx=(5.0d0*dx*dx*dx
     x          -3.0d0*dx*drsq)*threedr7
            tyyy=(5.0d0*dy*dy*dy
     x          -3.0d0*dy*drsq)*threedr7
            tzzz=(5.0d0*dz*dz*dz
     x          -3.0d0*dz*drsq)*threedr7

            txxy=(5.0d0*dx*dx*dy
     x          -drsq*dy)*threedr7
            txyy=(5.0d0*dx*dy*dy
     x          -drsq*dx)*threedr7
            txxz=(5.0d0*dx*dx*dz
     x          -drsq*dz)*threedr7
            txzz=(5.0d0*dx*dz*dz
     x          -drsq*dx)*threedr7
            tyyz=(5.0d0*dy*dy*dz
     x          -drsq*dz)*threedr7
            tyzz=(5.0d0*dy*dz*dz
     x          -drsq*dy)*threedr7

            txyz=15.0d0*dr7rec*dx*dy*dz
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
     x               +fdipdipdamp*efact2*dx
     x               +xmu(i)*xmujdotrij*efact2
     x               +xmu(j)*xmuidotrij*efact2
     x               -xmuidotrij*xmujdotrij*dx*efact3

            fmumuy = -fmumuy*efact
     x               +fdipdipdamp*efact2*dy
     x               +ymu(i)*xmujdotrij*efact2
     x               +ymu(j)*xmuidotrij*efact2
     x               -xmuidotrij*xmujdotrij*dy*efact3

            fmumuz = -fmumuz*efact
     x               +fdipdipdamp*efact2*dz
     x               +zmu(i)*xmujdotrij*efact2
     x               +zmu(j)*xmuidotrij*efact2
     x               -xmuidotrij*xmujdotrij*dz*efact3


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

            stsrxx=stsrxx+fttx*dx
            stsrxy=stsrxy+fttx*dy
            stsrxz=stsrxz+fttx*dz
            stsryy=stsryy+ftty*dy
            stsryz=stsryz+ftty*dz
            stsrzz=stsrzz+fttz*dz

c
c charge-charge term:
c
            stcxx=stcxx+egam*dx*dx
            stcxy=stcxy+egam*dx*dy
            stcxz=stcxz+egam*dx*dz
            stcyy=stcyy+egam*dy*dy
            stcyz=stcyz+egam*dy*dz
            stczz=stczz+egam*dz*dz
c
c charge-dipole term:
c
            stpxx=stpxx+frrxdipq*dx
            stpxy=stpxy+0.5d0*(frrxdipq*dy
     x                        +frrydipq*dx)
            stpxz=stpxz+0.5d0*(frrxdipq*dz
     x                        +frrzdipq*dx)
            stpyy=stpyy+frrydipq*dy
            stpyz=stpyz+0.5d0*(frrydipq*dz
     x                        +frrzdipq*dy)
            stpzz=stpzz+frrzdipq*dz
c
c dipole-dipole term:
c
            stp2xx=stp2xx+fmumux*dx
            stp2xy=stp2xy+0.5d0*(fmumux*dy
     x                          +fmumuy*dx)
            stp2xz=stp2xz+0.5d0*(fmumux*dz
     x                          +fmumuz*dx)
            stp2yy=stp2yy+fmumuy*dy
            stp2yz=stp2yz+0.5d0*(fmumuy*dz
     x                          +fmumuz*dy)
            stp2zz=stp2zz+fmumuz*dz
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

            fttx=gamtot*dx+frrxdipq+fmumux
            ftty=gamtot*dy+frrydipq+fmumuy
            fttz=gamtot*dz+frrzdipq+fmumuz

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

            elecsrxi=dx*r3dampi
            elecsryi=dy*r3dampi
            elecsrzi=dz*r3dampi

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

            elecsrxj=-dx*r3dampj
            elecsryj=-dy*r3dampj
            elecsrzj=-dz*r3dampj

            elecxsr(j)=elecxsr(j)+elecsrxj
            elecysr(j)=elecysr(j)+elecsryj
            eleczsr(j)=eleczsr(j)+elecsrzj

c
c Calculate contribution to ion forces:
c
            xmuidotrij=xmu(i)*dx
     x                +ymu(i)*dy
     x                +zmu(i)*dz
            xmujdotrij=xmu(j)*dx
     x                +ymu(j)*dy
     x                +zmu(j)*dz

            xmuidotT=(3.0d0*dx*xmuidotrij
     x              -drsq*xmu(i))*dr5rec
            ymuidotT=(3.0d0*dy*xmuidotrij
     x              -drsq*ymu(i))*dr5rec
            zmuidotT=(3.0d0*dz*xmuidotrij
     x              -drsq*zmu(i))*dr5rec

            xmujdotT=(3.0d0*dx*xmujdotrij
     x              -drsq*xmu(j))*dr5rec
            ymujdotT=(3.0d0*dy*xmujdotrij
     x              -drsq*ymu(j))*dr5rec
            zmujdotT=(3.0d0*dz*xmujdotrij
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

            dipi=xmu(i)*dx
     x         +ymu(i)*dy+zmu(i)*dz
            dipj=xmu(j)*dx
     x         +ymu(j)*dy+zmu(j)*dz

            diprfunci=dipi*drrec*q(j)
     x              *((-3.0d0*dampfunci*dr3rec*drrec)
     x              +(dampfuncdiffi*dr3rec))
            diprfuncj=dipj*drrec*q(i)
     x              *((-3.0d0*dampfuncj*dr3rec*drrec)
     x              +(dampfuncdiffj*dr3rec))

            frrxsri=(xmu(i)*r3dampi)+(dx*diprfunci)
            frrysri=(ymu(i)*r3dampi)+(dy*diprfunci)
            frrzsri=(zmu(i)*r3dampi)+(dz*diprfunci)
            frrxsrj=(xmu(j)*r3dampj)+(dx*diprfuncj)
            frrysrj=(ymu(j)*r3dampj)+(dy*diprfuncj)
            frrzsrj=(zmu(j)*r3dampj)+(dz*diprfuncj)

c
c Add these other s-r terms to the virial:
c
c need to be checked think 1st set right.
c Also need to add dip-dip contribution somewhere.

cc            xviraccsrdip=xviraccsrdip+((frrxsr)*dx)
cc            yviraccsrdip=yviraccsrdip+((frrysr)*dy)
cc            zviraccsrdip=zviraccsrdip+((frrzsr)*dz)

c
c Stress tensor accumulation:
c charge-dipole damping term:
c TMP - COME BACK TO THIS - MISSING TERMS....
            stpsrxx=stpsrxx+frrxsri*dx
     x                     -frrxsrj*dx
            stpsrxy=stpsrxy+0.5d0*(frrxsri*dy
     x                            +frrysri*dx)
     x                     -0.5d0*(frrxsrj*dy
     x                            +frrysrj*dx)
            stpsrxz=stpsrxz+0.5d0*(frrxsri*dz
     x                            +frrzsri*dx)
     x                     -0.5d0*(frrxsrj*dz
     x                            +frrzsrj*dx)
            stpsryy=stpsryy+frrysri*dy
     x                     -frrysrj*dy
            stpsryz=stpsryz+0.5d0*(frrysri*dz
     x                            +frrzsri*dy)
     x                     -0.5d0*(frrysrj*dz
     x                            +frrzsrj*dy)
            stpsrzz=stpsrzz+frrzsri*dz
     x                     -frrzsrj*dz

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
