      subroutine sr_energy

c************************************************************
c   Calculates the SR energy and forces at each site.
c
c   Possible models are the i) EPP -FT.
c                          ii) CIM for anion-cation + EPP for others
c                         iii) EPP -XFT.
c                          iv) DAIM
c                           v) QUAIM
c             
c************************************************************

      implicit double precision (a-h,o-z)

      include 'common.inc'
c
c Zero energy accumulators.
c
      engsr=0.0d0
      engsrat=0.0d0
      engsrrep=0.0d0

c Select SR energy to be calculated.
      if (epplog) then

         if ((.not.(xftlog)).and.(.not.(rvplog))
     x       .and.(.not.(ljpotlog))
     x       .and.(.not.(bkspotlog))) then
c Loop over all pairs calculating the
c standard Fumi-Tosi potentials and forces.
            do 30 i=2,num

               ipoint=ntype(i)

               do 40 j=1,i-1

                  jpoint=ntype(j)

            dxcf=x(i)-x(j)
            dycf=y(i)-y(j)
            dzcf=z(i)-z(j)

            dxcf=dxcf-boxlenx*int(dxcf*halfboxxrec)
            dycf=dycf-boxleny*int(dycf*halfboxyrec)
            dzcf=dzcf-boxlenz*int(dzcf*halfboxzrec)

            dx=hlab2(1,1)*dxcf
     x             +hlab2(1,2)*dycf+hlab2(1,3)*dzcf
            dy=hlab2(2,1)*dxcf
     x             +hlab2(2,2)*dycf+hlab2(2,3)*dzcf
            dz=hlab2(3,1)*dxcf
     x             +hlab2(3,2)*dycf+hlab2(3,3)*dzcf
                  drsq=dx*dx
     x                +dy*dy
     x                +dz*dz
c
c implement cut-off:
cc            if (drsq.ge.rsqmax) goto 40

                  drsqrec=1.0d0/drsq
                  dr=dsqrt(drsq)
                  drrec=1.0d0/dr

                  expft=ftb(ipoint,jpoint)*exp(-ftalp(ipoint,jpoint)
     x      *dr)

c Updated equation of motion to include new terms:
c
                  sgam=(expft*ftalp(ipoint,jpoint)*drrec)

                  fttx=sgam*dx
                  ftty=sgam*dy
                  fttz=sgam*dz

                  engsr=engsr+expft

c
c Stress tensor accumulation:
c short-range forces:
c
                  stsrxx=stsrxx+fttx*dx
                  stsrxy=stsrxy+fttx*dy
                  stsrxz=stsrxz+fttx*dz
                  stsryy=stsryy+ftty*dy
                  stsryz=stsryz+ftty*dz
                  stsrzz=stsrzz+fttz*dz

                  frrx(j)=frrx(j)-fttx
                  frry(j)=frry(j)-ftty
                  frrz(j)=frrz(j)-fttz

                  frrx(i)=frrx(i)+fttx
                  frry(i)=frry(i)+ftty
                  frrz(i)=frrz(i)+fttz

   40          continue

   30       continue

         endif

         if (rvplog) then

c Loop over all pairs calculating the
c standard RVP potentials and forces.
            do 301 i=2,num

               ipoint=ntype(i)

               do 401 j=1,i-1

                  jpoint=ntype(j)

            dxcf=x(i)-x(j)
            dycf=y(i)-y(j)
            dzcf=z(i)-z(j)

            dxcf=dxcf-boxlenx*int(dxcf*halfboxxrec)
            dycf=dycf-boxleny*int(dycf*halfboxyrec)
            dzcf=dzcf-boxlenz*int(dzcf*halfboxzrec)

            dx=hlab2(1,1)*dxcf
     x             +hlab2(1,2)*dycf+hlab2(1,3)*dzcf
            dy=hlab2(2,1)*dxcf
     x             +hlab2(2,2)*dycf+hlab2(2,3)*dzcf
            dz=hlab2(3,1)*dxcf
     x             +hlab2(3,2)*dycf+hlab2(3,3)*dzcf

            drsq=dx*dx+dy*dy+dz*dz
c
c implement cut-off:
cc            if (drsq.ge.rsqmax) goto 40

                  drsqrec=1.0d0/drsq
                  dr=dsqrt(drsq)
                  drrec=1.0d0/dr

c RVP energy. Note that the C_6 term (I-I only for AgI)
c is handled in the long-range energy routines.
c As a result, we consider only the 1/(R^n) SR term
c and 1/r^4 term here.
                  expft=(rvph(ipoint,jpoint)/(dr**rvpn(ipoint,jpoint)))
     x                 -(rvpr4(ipoint,jpoint)/(dr**4.0d0))

                  sgam=((rvpn(ipoint,jpoint)*rvph(ipoint,jpoint)
     x             /(dr**(rvpn(ipoint,jpoint)+1.0d0)))
     x             -(4.0d0*rvpr4(ipoint,jpoint)/(dr**5.0d0)))*drrec

                  fttx=sgam*dx
                  ftty=sgam*dy
                  fttz=sgam*dz

                  engsr=engsr+expft

c
c Stress tensor accumulation:
c short-range forces:
c
                  stsrxx=stsrxx+fttx*dx
                  stsrxy=stsrxy+fttx*dy
                  stsrxz=stsrxz+fttx*dz
                  stsryy=stsryy+ftty*dy
                  stsryz=stsryz+ftty*dz
                  stsrzz=stsrzz+fttz*dz

                  frrx(j)=frrx(j)-fttx
                  frry(j)=frry(j)-ftty
                  frrz(j)=frrz(j)-fttz

                  frrx(i)=frrx(i)+fttx
                  frry(i)=frry(i)+ftty
                  frrz(i)=frrz(i)+fttz

  401          continue

  301       continue

         endif

c Tangney-Scandalo
         if (tslog) then

c Loop over all pairs calculating the
c standard RVP potentials and forces.
            do 501 i=2,num

               ipoint=ntype(i)

               do 601 j=1,i-1

                  jpoint=ntype(j)

            dxcf=x(i)-x(j)
            dycf=y(i)-y(j)
            dzcf=z(i)-z(j)

            dxcf=dxcf-boxlenx*int(dxcf*halfboxxrec)
            dycf=dycf-boxleny*int(dycf*halfboxyrec)
            dzcf=dzcf-boxlenz*int(dzcf*halfboxzrec)

            dx=hlab2(1,1)*dxcf
     x             +hlab2(1,2)*dycf+hlab2(1,3)*dzcf
            dy=hlab2(2,1)*dxcf
     x             +hlab2(2,2)*dycf+hlab2(2,3)*dzcf
            dz=hlab2(3,1)*dxcf
     x             +hlab2(3,2)*dycf+hlab2(3,3)*dzcf
                  drsq=dx*dx
     x                +dy*dy
     x                +dz*dz

                  drsqrec=1.0d0/drsq
                  dr=dsqrt(drsq)
                  drrec=1.0d0/dr

                  term1=1.0d0-(dr/tsr(ipoint,jpoint))
                  term2=dexp(tsgamma(ipoint,jpoint)*term1)
                  term3=2.0d0*dexp(0.5d0*tsgamma(ipoint,jpoint)*term1)

                  expft=tsd(ipoint,jpoint)*(term2-term3)

                  sgam1=-1.0d0*term2*tsd(ipoint,jpoint)
     x                *tsgamma(ipoint,jpoint)/tsr(ipoint,jpoint)
                  sgam2=0.5d0*term3*tsd(ipoint,jpoint)
     x                *tsgamma(ipoint,jpoint)/tsr(ipoint,jpoint)
                  sgam=-(sgam1+sgam2)*drrec

                  fttx=sgam*dx
                  ftty=sgam*dy
                  fttz=sgam*dz

                  engsr=engsr+expft
c
c Stress tensor accumulation:
c short-range forces:
c
                  stsrxx=stsrxx+fttx*dx
                  stsrxy=stsrxy+fttx*dy
                  stsrxz=stsrxz+fttx*dz
                  stsryy=stsryy+ftty*dy
                  stsryz=stsryz+ftty*dz
                  stsrzz=stsrzz+fttz*dz

                  frrx(j)=frrx(j)-fttx
                  frry(j)=frry(j)-ftty
                  frrz(j)=frrz(j)-fttz

                  frrx(i)=frrx(i)+fttx
                  frry(i)=frry(i)+ftty
                  frrz(i)=frrz(i)+fttz

  601          continue

  501       continue

         endif

c Bexp(-ar)/r + B'exp(-a'r^n) potential (i.e PbF2 work).

         if (xftlog) then

            do 300 i=2,num

               ipoint=ntype(i)

               do 400 j=1,i-1

                  jpoint=ntype(j)

            dxcf=x(i)-x(j)
            dycf=y(i)-y(j)
            dzcf=z(i)-z(j)

            dxcf=dxcf-boxlenx*int(dxcf*halfboxxrec)
            dycf=dycf-boxleny*int(dycf*halfboxyrec)
            dzcf=dzcf-boxlenz*int(dzcf*halfboxzrec)

            dx=hlab2(1,1)*dxcf
     x             +hlab2(1,2)*dycf+hlab2(1,3)*dzcf
            dy=hlab2(2,1)*dxcf
     x             +hlab2(2,2)*dycf+hlab2(2,3)*dzcf
            dz=hlab2(3,1)*dxcf
     x             +hlab2(3,2)*dycf+hlab2(3,3)*dzcf

            drsq=dx*dx+dy*dy+dz*dz
c
c implement cut-off:
cc            if (drsq.ge.rsqmax) goto 40

                  drsqrec=1.0d0/drsq
                  dr=dsqrt(drsq)
                  drrec=1.0d0/dr

                  drnrec=drrec**nrpower(ipoint,jpoint)

                  expft=ftb(ipoint,jpoint)*exp(-ftalp(ipoint,jpoint)
     x      *dr)*drnrec

                  if (ftbx(ipoint,jpoint).gt.0.0d0) then 
                     expftx=ftbx(ipoint,jpoint)
     x                     *exp(-ftalpx(ipoint,jpoint)*drsq)
                  else
                     expftx=0.0d0
                  endif

                  sgam=expft*ftalp(ipoint,jpoint)*drrec
     x          +(expft*nrpower(ipoint,jpoint)*drsqrec)
     x          +expftx*2.0d0*ftalpx(ipoint,jpoint)

c Updated equation of motion to include new terms:
c

                  fttx=sgam*dx
                  ftty=sgam*dy
                  fttz=sgam*dz

                  engsr=engsr+expftx+expft
c
c Stress tensor accumulation:
c short-range forces:
c
                  stsrxx=stsrxx+fttx*dx
                  stsrxy=stsrxy+fttx*dy
                  stsrxz=stsrxz+fttx*dz
                  stsryy=stsryy+ftty*dy
                  stsryz=stsryz+ftty*dz
                  stsrzz=stsrzz+fttz*dz

                  frrx(j)=frrx(j)-fttx
                  frry(j)=frry(j)-ftty
                  frrz(j)=frrz(j)-fttz

                  frrx(i)=frrx(i)+fttx
                  frry(i)=frry(i)+ftty
                  frrz(i)=frrz(i)+fttz

  400          continue

  300       continue

         endif

      endif

c anion-anion EPP
c TMP - CIM, AIM.
      if ((cimlog).or.(daimlog).or.(quaimlog)) then

         do 70 i=2,nanion

            ipoint=ntype(i)

            do 80 j=1,i-1

               jpoint=ntype(j)

            dxcf=x(i)-x(j)
            dycf=y(i)-y(j)
            dzcf=z(i)-z(j)

            dxcf=dxcf-boxlenx*int(dxcf*halfboxxrec)
            dycf=dycf-boxleny*int(dycf*halfboxyrec)
            dzcf=dzcf-boxlenz*int(dzcf*halfboxzrec)

            dx=hlab2(1,1)*dxcf
     x             +hlab2(1,2)*dycf+hlab2(1,3)*dzcf
            dy=hlab2(2,1)*dxcf
     x             +hlab2(2,2)*dycf+hlab2(2,3)*dzcf
            dz=hlab2(3,1)*dxcf
     x             +hlab2(3,2)*dycf+hlab2(3,3)*dzcf

            drsq=dx*dx+dy*dy+dz*dz
c
c implement cut-off:
cc            if (drsq.ge.rsqmax) goto 40

               drsqrec=1.0d0/drsq
               dr=dsqrt(drsq)
               drrec=1.0d0/dr

               expft=ftb(ipoint,jpoint)*exp(-ftalp(ipoint,jpoint)
     x      *dr)

c Updated equation of motion to include new terms:
c
               sgam=(expft*ftalp(ipoint,jpoint)*drrec)

               fttx=sgam*dx
               ftty=sgam*dy
               fttz=sgam*dz

               engsr=engsr+expft

c
c Stress tensor accumulation:
c short-range forces:
c
               stsrxx=stsrxx+fttx*dx
               stsrxy=stsrxy+fttx*dy
               stsrxz=stsrxz+fttx*dz
               stsryy=stsryy+ftty*dy
               stsryz=stsryz+ftty*dz
               stsrzz=stsrzz+fttz*dz

               frrx(j)=frrx(j)-fttx
               frry(j)=frry(j)-ftty
               frrz(j)=frrz(j)-fttz

               frrx(i)=frrx(i)+fttx
               frry(i)=frry(i)+ftty
               frrz(i)=frrz(i)+fttz

   80       continue

   70    continue

c cation-cation EPP
c TMP - CIM
         do 90 i=nanion+2,num

            ipoint=ntype(i)

            do 100 j=1+nanion,i-1

               jpoint=ntype(j)

            dxcf=x(i)-x(j)
            dycf=y(i)-y(j)
            dzcf=z(i)-z(j)

            dxcf=dxcf-boxlenx*int(dxcf*halfboxxrec)
            dycf=dycf-boxleny*int(dycf*halfboxyrec)
            dzcf=dzcf-boxlenz*int(dzcf*halfboxzrec)

            dx=hlab2(1,1)*dxcf
     x             +hlab2(1,2)*dycf+hlab2(1,3)*dzcf
            dy=hlab2(2,1)*dxcf
     x             +hlab2(2,2)*dycf+hlab2(2,3)*dzcf
            dz=hlab2(3,1)*dxcf
     x             +hlab2(3,2)*dycf+hlab2(3,3)*dzcf

            drsq=dx*dx+dy*dy+dz*dz
c
c implement cut-off:
cc            if (drsq.ge.rsqmax) goto 40

               drsqrec=1.0d0/drsq
               dr=dsqrt(drsq)
               drrec=1.0d0/dr

               expft=ftb(ipoint,jpoint)*exp(-ftalp(ipoint,jpoint)
     x      *dr)

c Updated equation of motion to include new terms:
c
               sgam=(expft*ftalp(ipoint,jpoint)*drrec)

               fttx=sgam*dx
               ftty=sgam*dy
               fttz=sgam*dz

               engsr=engsr+expft
c
c Stress tensor accumulation:
c short-range forces:
c
               stsrxx=stsrxx+fttx*dx
               stsrxy=stsrxy+fttx*dy
               stsrxz=stsrxz+fttx*dz
               stsryy=stsryy+ftty*dy
               stsryz=stsryz+ftty*dz
               stsrzz=stsrzz+fttz*dz

               frrx(j)=frrx(j)-fttx
               frry(j)=frry(j)-ftty
               frrz(j)=frrz(j)-fttz

               frrx(i)=frrx(i)+fttx
               frry(i)=frry(i)+ftty
               frrz(i)=frrz(i)+fttz

  100       continue

   90    continue

      endif

      if (ljpotlog) then

         stsrxx=0.0d0
         stsrxy=0.0d0
         stsrxz=0.0d0
         stsryy=0.0d0
         stsryz=0.0d0
         stsrzz=0.0d0

         engpetot=0.0d0

         do 290 i=1,num

            frrx(i)=0.0d0
            frry(i)=0.0d0
            frrz(i)=0.0d0

  290    continue

         do 305 i=2,num

            ipoint=ntype(i)

            do 310 j=1,i-1

               jpoint=ntype(j)

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
     x             +dy*dy
     x             +dz*dz

               dr=dsqrt(drsq)
               if (dr.gt.ljcutoff) goto 310

               drrec=1.0d0/dr
               drsqrec=drrec*drrec
               dr6rec=drrec**6.0d0
               dr12rec=drrec**12.0d0

               ulj=fourepsil(ipoint,jpoint)
     x                     *((sig12(ipoint,jpoint)*dr12rec)
     x                     -(sig6(ipoint,jpoint)*dr6rec))

               gam=12.0d0*fourepsil(ipoint,jpoint)*sig12(ipoint,jpoint)
     x         *dr12rec*drsqrec
     x         -6.0d0*fourepsil(ipoint,jpoint)*sig6(ipoint,jpoint)
     x         *dr6rec*drsqrec

               engsr=engsr+ulj
               engsrrep=engsrrep+fourepsil(ipoint,jpoint)
     x                     *sig12(ipoint,jpoint)*dr12rec
               engsrat=engsrat-fourepsil(ipoint,jpoint)
     x                     *sig6(ipoint,jpoint)*dr6rec

               frrx(j)=frrx(j)-dx*gam
               frry(j)=frry(j)-dy*gam
               frrz(j)=frrz(j)-dz*gam
               frrx(i)=frrx(i)+dx*gam
               frry(i)=frry(i)+dy*gam
               frrz(i)=frrz(i)+dz*gam

               stsrxx=stsrxx+dx*gam*dx
               stsrxy=stsrxy+dx*gam*dy
               stsrxz=stsrxz+dx*gam*dz
               stsryy=stsryy+dy*gam*dy
               stsryz=stsryz+dy*gam*dz
               stsrzz=stsrzz+dz*gam*dz

  310       continue

  305    continue

      endif

cc double well potential.
      if (dwpotlog) then

         stsrxx=0.0d0
         stsrxy=0.0d0
         stsrxz=0.0d0
         stsryy=0.0d0
         stsryz=0.0d0
         stsrzz=0.0d0

         engpetot=0.0d0

         do 291 i=1,num

            frrx(i)=0.0d0
            frry(i)=0.0d0
            frrz(i)=0.0d0

  291    continue

         do 355 i=2,num

            ipoint=ntype(i)

            do 356 j=1,i-1

               jpoint=ntype(j)

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
     x             +dy*dy
     x             +dz*dz

               dr=dsqrt(drsq)
               if (dr.gt.ljcutoff) goto 356

               drrec=1.0d0/dr
               drsqrec=drrec*drrec
               dr6rec=drrec**6.0d0
               dr12rec=drrec**12.0d0

c LJ potential.
               ulj=fourepsil(ipoint,jpoint)
     x                     *((sig12(ipoint,jpoint)*dr12rec)
     x                     -(sig6(ipoint,jpoint)*dr6rec))

c ggg is the gaussian function.
               ggg=Adw*dexp(-alphadw*((dr-r0dw)**2.0d0))
               term4=2.0d0*alphadw*(dr-r0dw)
               u2=ggg*term4

               utot=ulj*(u2+1.0d0)

cc               gam=12.0d0*fourepsil(ipoint,jpoint)*sig12(ipoint,jpoint)
cc     x         *dr12rec*drsqrec
cc     x         -6.0d0*fourepsil(ipoint,jpoint)*sig6(ipoint,jpoint)
cc     x         *dr6rec*drsqrec

               uljderiv=12.0d0*fourepsil(ipoint,jpoint)
     x                 *sig12(ipoint,jpoint)
     x         *dr12rec*drrec
     x         -6.0d0*fourepsil(ipoint,jpoint)*sig6(ipoint,jpoint)
     x         *dr6rec*drrec

               uderiv1=uljderiv*(u2+1.0d0)

               term5=2.0d0*alphadw*ggg
               term6=4.0d0*alphadw*alphadw*((dr-r0dw)**2.0d0)*ggg

               du2dr=term5-term6

               dudr=uderiv1-(ulj*du2dr)

               gam=dudr*drrec

               engsr=engsr+utot

               frrx(j)=frrx(j)-dx*gam
               frry(j)=frry(j)-dy*gam
               frrz(j)=frrz(j)-dz*gam
               frrx(i)=frrx(i)+dx*gam
               frry(i)=frry(i)+dy*gam
               frrz(i)=frrz(i)+dz*gam

               stsrxx=stsrxx+dx*gam*dx
               stsrxy=stsrxy+dx*gam*dy
               stsrxz=stsrxz+dx*gam*dz
               stsryy=stsryy+dy*gam*dy
               stsryz=stsryz+dy*gam*dz
               stsrzz=stsrzz+dz*gam*dz

  356       continue

  355    continue

      endif

cc double well potential (type II).
      if (dwpot2log) then

       write(6,*)'dwdwdwdwdw2'
         stsrxx=0.0d0
         stsrxy=0.0d0
         stsrxz=0.0d0
         stsryy=0.0d0
         stsryz=0.0d0
         stsrzz=0.0d0

         engpetot=0.0d0

         do 292 i=1,num

            frrx(i)=0.0d0
            frry(i)=0.0d0
            frrz(i)=0.0d0

  292    continue

         do 357 i=2,num

            ipoint=ntype(i)

            do 358 j=1,i-1

               jpoint=ntype(j)

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
     x             +dy*dy
     x             +dz*dz

               dr=dsqrt(drsq)

               if (dr.gt.ljcutoff) goto 358

               drrec=1.0d0/dr
               drsqrec=drrec*drrec
               dr6rec=drrec**6.0d0
               dr12rec=drrec**12.0d0

c LJ potential.
               ulj=fourepsil(ipoint,jpoint)
     x                     *((sig12(ipoint,jpoint)*dr12rec)
     x                     -(sig6(ipoint,jpoint)*dr6rec))

c ggg is the gaussian function.
               u2=Adw*dexp(-alphadw*((dr-r0dw)**2.0d0))
cc This was the original formalism, later changed.
cc               u2=1.0d0-ggg

               utot=ulj-u2

cc               gam=12.0d0*fourepsil(ipoint,jpoint)*sig12(ipoint,jpoint)
cc     x         *dr12rec*drsqrec
cc     x         -6.0d0*fourepsil(ipoint,jpoint)*sig6(ipoint,jpoint)
cc     x         *dr6rec*drsqrec

               uljderiv=12.0d0*fourepsil(ipoint,jpoint)
     x                 *sig12(ipoint,jpoint)
     x         *dr12rec*drrec
     x         -6.0d0*fourepsil(ipoint,jpoint)*sig6(ipoint,jpoint)
     x         *dr6rec*drrec

               term4=2.0d0*alphadw*(dr-r0dw)

               du2dr=term4*u2

               dudr=uljderiv-du2dr

               gam=dudr*drrec

               engsr=engsr+utot

               frrx(j)=frrx(j)-dx*gam
               frry(j)=frry(j)-dy*gam
               frrz(j)=frrz(j)-dz*gam
               frrx(i)=frrx(i)+dx*gam
               frry(i)=frry(i)+dy*gam
               frrz(i)=frrz(i)+dz*gam

               stsrxx=stsrxx+dx*gam*dx
               stsrxy=stsrxy+dx*gam*dy
               stsrxz=stsrxz+dx*gam*dz
               stsryy=stsryy+dy*gam*dy
               stsryz=stsryz+dy*gam*dz
               stsrzz=stsrzz+dz*gam*dz

  358       continue

  357    continue

      endif

cc gaussian core potential.
      if (gcpotlog) then

         stsrxx=0.0d0
         stsrxy=0.0d0
         stsrxz=0.0d0
         stsryy=0.0d0
         stsryz=0.0d0
         stsrzz=0.0d0

         engpetot=0.0d0

         do 294 i=1,num

            frrx(i)=0.0d0
            frry(i)=0.0d0
            frrz(i)=0.0d0

  294    continue

         do 359 i=2,num

            ipoint=ntype(i)

            do 361 j=1,i-1

               jpoint=ntype(j)

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
     x             +dy*dy
     x             +dz*dz

               dr=dsqrt(drsq)
cc               if (dr.gt.ljcutoff) goto 356

               drrec=1.0d0/dr
               drsqrec=drrec*drrec

c GC potential.
               ugc=epsilon(ipoint,jpoint)
     x      *dexp(-drsq/(sigma(ipoint,jpoint)*sigma(ipoint,jpoint)))

         dudr=(-2.0d0*dr/(sigma(ipoint,jpoint)*sigma(ipoint,jpoint)))
     x       *ugc

               gam=-dudr*drrec

               engsr=engsr+ugc

               frrx(j)=frrx(j)-dx*gam
               frry(j)=frry(j)-dy*gam
               frrz(j)=frrz(j)-dz*gam
               frrx(i)=frrx(i)+dx*gam
               frry(i)=frry(i)+dy*gam
               frrz(i)=frrz(i)+dz*gam

               stsrxx=stsrxx+dx*gam*dx
               stsrxy=stsrxy+dx*gam*dy
               stsrxz=stsrxz+dx*gam*dz
               stsryy=stsryy+dy*gam*dy
               stsryz=stsryz+dy*gam*dz
               stsrzz=stsrzz+dz*gam*dz

  361       continue

  359    continue

      endif

      if (bkspotlog) then

cc         stsrxx=0.0d0
cc         stsrxy=0.0d0
cc         stsrxz=0.0d0
cc         stsryy=0.0d0
cc         stsryz=0.0d0
cc         stsrzz=0.0d0

cc         engpetot=0.0d0

cc         do 290 i=1,num

cc            frrx(i)=0.0d0
cc            frry(i)=0.0d0
cc            frrz(i)=0.0d0

cc  290    continue

         engsrft=0.0d0
         do 505 i=2,num

            ipoint=ntype(i)

            do 511 j=1,i-1

               jpoint=ntype(j)

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
     x             +dy*dy
     x             +dz*dz

               dr=dsqrt(drsq)
cc               if (dr.gt.ljcutoff) goto 310

               drrec=1.0d0/dr
               drsqrec=drrec*drrec
               dr6rec=drrec**6.0d0
               dr30rec=drrec**30.0d0

               ulj=fourepsil(ipoint,jpoint)
     x                     *((sig30(ipoint,jpoint)*dr30rec)
     x                     -(sig6(ipoint,jpoint)*dr6rec))

               gam=30.0d0*fourepsil(ipoint,jpoint)*sig30(ipoint,jpoint)
     x         *dr30rec*drsqrec
     x         -6.0d0*fourepsil(ipoint,jpoint)*sig6(ipoint,jpoint)
     x         *dr6rec*drsqrec

               engsr=engsr+ulj
               engsrrep=engsrrep+fourepsil(ipoint,jpoint)
     x                     *sig30(ipoint,jpoint)*dr30rec
               engsrat=engsrat-fourepsil(ipoint,jpoint)
     x                     *sig6(ipoint,jpoint)*dr6rec

               frrx(j)=frrx(j)-dx*gam
               frry(j)=frry(j)-dy*gam
               frrz(j)=frrz(j)-dz*gam
               frrx(i)=frrx(i)+dx*gam
               frry(i)=frry(i)+dy*gam
               frrz(i)=frrz(i)+dz*gam

               stsrxx=stsrxx+dx*gam*dx
               stsrxy=stsrxy+dx*gam*dy
               stsrxz=stsrxz+dx*gam*dz
               stsryy=stsryy+dy*gam*dy
               stsryz=stsryz+dy*gam*dz
               stsrzz=stsrzz+dz*gam*dz

                  drsqrec=1.0d0/drsq
                  drrec=1.0d0/dr

                  expft=ftb(ipoint,jpoint)*exp(-ftalp(ipoint,jpoint)
     x      *dr)

c Updated equation of motion to include new terms:
c
                  sgam=(expft*ftalp(ipoint,jpoint)*drrec)

                  fttx=sgam*dx
                  ftty=sgam*dy
                  fttz=sgam*dz

                  engsr=engsr+expft
                  engsrft=engsrft+expft

c
c Stress tensor accumulation:
c short-range forces:
c
                  stsrxx=stsrxx+fttx*dx
                  stsrxy=stsrxy+fttx*dy
                  stsrxz=stsrxz+fttx*dz
                  stsryy=stsryy+ftty*dy
                  stsryz=stsryz+ftty*dz
                  stsrzz=stsrzz+fttz*dz

                  frrx(j)=frrx(j)-fttx
                  frry(j)=frry(j)-ftty
                  frrz(j)=frrz(j)-fttz

                  frrx(i)=frrx(i)+fttx
                  frry(i)=frry(i)+ftty
                  frrz(i)=frrz(i)+fttz

  511       continue

  505    continue

      endif

      write(87,*)engpetot,engsr
      write(86,*)engsrat,engsrrep,engsrft
      engpetot=engpetot+engsr
cc      write(88,*)engmorse

      return

      end
