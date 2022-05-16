      subroutine rgdrealE

c************************************************************
c   Calculates the real energy and forces at each site.
c   The calculated energy is the real part of the Ewald
c  summation plus the short-range and Van der Waals
c  corrections from the Fumi-Tosi potential.
c
c  April 2000 - corrected Ewald.
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

         srxx(i)=0.0d0
         sryy(i)=0.0d0
         srzz(i)=0.0d0
         srxy(i)=0.0d0
         srxz(i)=0.0d0
         sryz(i)=0.0d0

         elecxu(i)=0.0d0
         elecyu(i)=0.0d0
         eleczu(i)=0.0d0

    5 continue
c
c Zero energy accumulators.
c
      eng=0.0d0

      do 6 i=1,nspec

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
c
c zero the short range part, the other bits are initialized in recipE:
c
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
            dr=0.0d0
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

            xviraccftsrc=xviraccftsrc+(sgam*dx*dx)
            yviraccftsrc=yviraccftsrc+(sgam*dy*dy)
            zviraccftsrc=zviraccftsrc+(sgam*dz*dz)

            elecx(j)=elecx(j)-(ettx*qjrec)
            elecy(j)=elecy(j)-(etty*qjrec)
            elecz(j)=elecz(j)-(ettz*qjrec)

            elecx(i)=elecx(i)+(ettx*qirec)
            elecy(i)=elecy(i)+(etty*qirec)
            elecz(i)=elecz(i)+(ettz*qirec)

c===============================================
c  accumulate all the real-space coulombic contributions to the
c  energy (not short-range dipoles)

            eng=eng+erfcr+cpe

c=======================================================
c
c update accumulators:
c
            erfcracc=erfcracc+erfcr

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

            fttx=gamtot*dx
            ftty=gamtot*dy
            fttz=gamtot*dz

            frrx(j)=frrx(j)-fttx
            frry(j)=frry(j)-ftty
            frrz(j)=frrz(j)-fttz

            frrx(i)=frrx(i)+fttx
            frry(i)=frry(i)+ftty
            frrz(i)=frrz(i)+fttz

   40    continue

   30 continue

      do 140 i=1,num

         elecxu(i)=elecx(i)
         elecyu(i)=elecy(i)
         eleczu(i)=elecz(i)

  140 continue

c
c energy accumulator for pressure calculation:
c multiply qdipacc by 2 as it is an r^-2 term
c similarly for dipdipacc (r^-3)
c toteng contains reciprocal space q-q energy.
c
      qqacc=erfcracc+toteng
    
      srvirft=xviraccftsrc+yviraccftsrc+zviraccftsrc
      srdipvir=xviraccsrdip+yviraccsrdip+zviraccsrdip
      virtot=srvirft+srdipvir

      virtotx=xviraccftsrc+xviraccsrdip
      virtoty=yviraccftsrc+yviraccsrdip
      virtotz=zviraccftsrc+zviraccsrdip

c..........estimate of the coulomb contribution to the pressure
      encoul=engpetot+erfcracc
c
c Final energy summation:
c
      engpetot=engpetot+eng

      return

      end
