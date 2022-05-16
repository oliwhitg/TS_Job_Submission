      subroutine sr_energy

c************************************************************
c   Calculates the SR energy and forces at each site.
c
c   Possible models are the i) EPP -FT.
c                          ii) CIM for anion-cation + EPP for others
c                         iii) EPP -XFT.
c                          iv) 
c                           v) QUAIM
c             
c************************************************************

      implicit double precision (a-h,o-z)

      include 'common.inc'
c
c Zero energy accumulators.
c
      engsr=0.0d0
c         write(*,*) "begin sr_energy" 
c         write(*,*) "engpetot", engpetot
         write(*,*) "engsr", engsr
      do 10 i=1,num
         engft1(i)=0.0d0
         engft2(i)=0.0d0
         engft3(i)=0.0d0
         engft1dotx(i)=0.0d0
         engft1doty(i)=0.0d0
         engft1dotz(i)=0.0d0
         engft2dotx(i)=0.0d0
         engft2doty(i)=0.0d0
         engft2dotz(i)=0.0d0
         engft3dotx(i)=0.0d0
         engft3doty(i)=0.0d0
         engft3dotz(i)=0.0d0

         do j=1,3
            engftdotxx(i,j)=0.0d0
            engftdotyy(i,j)=0.0d0
            engftdotzz(i,j)=0.0d0
            engftdotxy(i,j)=0.0d0
            engftdotxz(i,j)=0.0d0
            engftdotyz(i,j)=0.0d0
         enddo

   10 continue
c         write(*,*) ftlog, epplog
c Select SR energy to be calculated.
      if ((epplog).or.(ftlog)) then
           write(*,*) "inner"
         if ((.not.(xftlog)).and.(.not.(rvplog))) then
           write(*,*) "outer"
c Loop over all pairs calculating the
c standard Fumi-Tosi potentials and forces.
            do 30 i=2,num

               ipoint=ntype(i)

               do 40 j=1,i-1

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
               drsq=dxcf*dxcf
     x             +dycf*dycf
     x             +dzcf*dzcf
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

                  fttx=sgam*dxsav(i,j)
                  ftty=sgam*dysav(i,j)
                  fttz=sgam*dzsav(i,j)
cc
                  engsr=engsr+expft

c
c Stress tensor accumulation:
c short-range forces:
c
cc                  stsrxx=stsrxx+fttx*dxsav(i,j)
cc                  stsrxy=stsrxy+fttx*dysav(i,j)
cc                  stsrxz=stsrxz+fttx*dzsav(i,j)
cc                  stsryy=stsryy+ftty*dysav(i,j)
cc                  stsryz=stsryz+ftty*dzsav(i,j)
cc                  stsrzz=stsrzz+fttz*dzsav(i,j)

                  frrx(j)=frrx(j)-fttx
                  frry(j)=frry(j)-ftty
                  frrz(j)=frrz(j)-fttz

                  frrx(i)=frrx(i)+fttx
                  frry(i)=frry(i)+ftty
                  frrz(i)=frrz(i)+fttz

   40          continue

   30       continue
         write(*,*) "end 1 sr_energy"
         write(*,*) "engpetot", engpetot
         write(*,*) "engsr", engsr

         eng_ion=engpetot+engsr
         endif

      write(6,*)'rvp1'
         if (rvplog) then
      write(6,*)'rvp2'

c Loop over all pairs calculating the
c standard RVP potentials and forces.
            do 301 i=2,num

               ipoint=ntype(i)

               do 401 j=1,i-1

                  jpoint=ntype(j)

cc                  drsq=dxsav(i,j)*dxsav(i,j)
cc     x                +dysav(i,j)*dysav(i,j)
cc     x                +dzsav(i,j)*dzsav(i,j)
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

c                  write(*,*) "rvph",rvph(ipoint,jpoint)
c              write(*,*) "dr",dr , "rvpn", rvpn(ipoint,jpoint)
c                  write(*,*)"rvpr4" ,rvpr4(ipoint,jpoint)

                  sgam=((rvpn(ipoint,jpoint)*rvph(ipoint,jpoint)
     x             /(dr**(rvpn(ipoint,jpoint)+1.0d0)))
     x             -(4.0d0*rvpr4(ipoint,jpoint)/(dr**5.0d0)))*drrec

cc                  fttx=sgam*dxsav(i,j)
cc                  ftty=sgam*dysav(i,j)
cc                  fttz=sgam*dzsav(i,j)
c                  write(*,*) engsr, "1 sr"
                  engsr=engsr+expft
c                  write(*,*) engsr, "2 sr"

         eng_ion=engpetot+engsr
c
c Stress tensor accumulation:
c short-range forces:
c
cc                  stsrxx=stsrxx+fttx*dxsav(i,j)
cc                  stsrxy=stsrxy+fttx*dysav(i,j)
cc                  stsrxz=stsrxz+fttx*dzsav(i,j)
cc                  stsryy=stsryy+ftty*dysav(i,j)
cc                  stsryz=stsryz+ftty*dzsav(i,j)
cc                  stsrzz=stsrzz+fttz*dzsav(i,j)

                  frrx(j)=frrx(j)-fttx
                  frry(j)=frry(j)-ftty
                  frrz(j)=frrz(j)-fttz

                  frrx(i)=frrx(i)+fttx
                  frry(i)=frry(i)+ftty
                  frrz(i)=frrz(i)+fttz

  401          continue

  301       continue
             write(*,*)"ion", eng_ion
         eng_ion=engpetot+engsr
         endif

c Bexp(-ar)/r + B'exp(-a'r^n) potential (i.e PbF2 work).

         if (xftlog) then

            do 300 i=2,num

               ipoint=ntype(i)

               do 400 j=1,i-1

                  jpoint=ntype(j)

cc                  drsq=dxsav(i,j)*dxsav(i,j)
cc     x                +dysav(i,j)*dysav(i,j)
cc     x                +dzsav(i,j)*dzsav(i,j)
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

cc                  fttx=sgam*dxsav(i,j)
cc                  ftty=sgam*dysav(i,j)
cc                  fttz=sgam*dzsav(i,j)

                  engsr=engsr+expftx+expft
c
c Stress tensor accumulation:
c short-range forces:
c
cc                  stsrxx=stsrxx+fttx*dxsav(i,j)
cc                  stsrxy=stsrxy+fttx*dysav(i,j)
cc                  stsrxz=stsrxz+fttx*dzsav(i,j)
cc                  stsryy=stsryy+ftty*dysav(i,j)
cc                  stsryz=stsryz+ftty*dzsav(i,j)
cc                  stsrzz=stsrzz+fttz*dzsav(i,j)

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

c CIM energy + force calc.
      if (cimlog) then

c TMP - CIM1
         if (cim1log) then 

            do 50 i=1,nsp(1)

               ipoint=ntype(i)

               do 60 j=nsp(1)+1,num

                  jpoint=ntype(j)

cc                  drsq=dxsav(i,j)*dxsav(i,j)
cc     x             +dysav(i,j)*dysav(i,j)
cc     x             +dzsav(i,j)*dzsav(i,j)
c
c implement cut-off:
cc            if (drsq.ge.rsqmax) goto 40

                  drsqrec=1.0d0/drsq
                  dr=dsqrt(drsq)
                  drrec=1.0d0/dr

                  expft=(ftb(ipoint,jpoint)/dr)
     x              *exp(-ftalp(ipoint,jpoint)*
     x              (dr-delta(i)))
                  expft2=(ftb2(ipoint,jpoint)/dr)
     x              *exp(-ftbeta(ipoint,jpoint)*
     x              (dr-delta(i)))
                  expft3=(ftb3(ipoint,jpoint)/dr)
     x              *exp(-ftgamma(ipoint,jpoint)*
     x              (dr-delta(i)))

                  engft1(i)=engft1(i)+expft
                  engft2(i)=engft2(i)+expft2
                  engft3(i)=engft3(i)+expft3

                  sgam=(ftalp(ipoint,jpoint)*expft)+(expft*drrec)
     x          +(ftbeta(ipoint,jpoint)*expft2)+(expft2*drrec)
     x          +(ftgamma(ipoint,jpoint)*expft3)+(expft3*drrec)
cc                  fttx=sgam*(dxsav(i,j)*drrec)
cc                  ftty=sgam*(dysav(i,j)*drrec)
cc                  fttz=sgam*(dzsav(i,j)*drrec)

                  engsr=engsr+expft+expft2+expft3
c
c Stress tensor accumulation:
c short-range forces:
c
cc                  stsrxx=stsrxx+fttx*dxsav(i,j)
cc                  stsrxy=stsrxy+fttx*dysav(i,j)
cc                  stsrxz=stsrxz+fttx*dzsav(i,j)
cc                  stsryy=stsryy+ftty*dysav(i,j)
cc                  stsryz=stsryz+ftty*dzsav(i,j)
cc                  stsrzz=stsrzz+fttz*dzsav(i,j)

                  frrx(j)=frrx(j)-fttx
                  frry(j)=frry(j)-ftty
                  frrz(j)=frrz(j)-fttz

                  frrx(i)=frrx(i)+fttx
                  frry(i)=frry(i)+ftty
                  frrz(i)=frrz(i)+fttz

   60          continue

   50       continue

         endif

c TMP - CIM2
         if (cim2log) then 

            do 55 i=1,nsp(1)

               ipoint=ntype(i)

               do 65 j=nsp(1)+1,num

                  jpoint=ntype(j)

cc                  drsq=dxsav(i,j)*dxsav(i,j)
cc     x             +dysav(i,j)*dysav(i,j)
cc     x             +dzsav(i,j)*dzsav(i,j)
c
c implement cut-off:
cc            if (drsq.ge.rsqmax) goto 40

                  drsqrec=1.0d0/drsq
                  dr=dsqrt(drsq)
                  drrec=1.0d0/dr

                  expft=(ftb(ipoint,jpoint))
     x              *exp(-ftalp(ipoint,jpoint)*
     x              (dr-delta(i)))
                  expft2=(ftb2(ipoint,jpoint))
     x              *exp(-ftbeta(ipoint,jpoint)*
     x              (dr-delta(i)))
                  expft3=(ftb3(ipoint,jpoint))
     x              *exp(-ftgamma(ipoint,jpoint)*
     x              (dr-delta(i)))
                  expft4=extrab*exp(-extraalpha*dr)

                  engft1(i)=engft1(i)+expft
                  engft2(i)=engft2(i)+expft2
                  engft3(i)=engft3(i)+expft3

                  sgam=(ftalp(ipoint,jpoint)*expft)
     x          +(ftbeta(ipoint,jpoint)*expft2)
     x          +(ftgamma(ipoint,jpoint)*expft3)
     x          +expft4*extraalpha

cc                  fttx=sgam*(dxsav(i,j)*drrec)
cc                  ftty=sgam*(dysav(i,j)*drrec)
cc                  fttz=sgam*(dzsav(i,j)*drrec)

                  engsr=engsr+expft+expft2+expft3+expft4
c
c Stress tensor accumulation:
c short-range forces:
c
cc                  stsrxx=stsrxx+fttx*dxsav(i,j)
cc                  stsrxy=stsrxy+fttx*dysav(i,j)
cc                  stsrxz=stsrxz+fttx*dzsav(i,j)
cc                  stsryy=stsryy+ftty*dysav(i,j)
cc                  stsryz=stsryz+ftty*dzsav(i,j)
cc                  stsrzz=stsrzz+fttz*dzsav(i,j)

                  frrx(j)=frrx(j)-fttx
                  frry(j)=frry(j)-ftty
                  frrz(j)=frrz(j)-fttz

                  frrx(i)=frrx(i)+fttx
                  frry(i)=frry(i)+ftty
                  frrz(i)=frrz(i)+fttz

   65          continue

   55       continue

         endif

      endif

c DAIM energy + force calc.
      if (daimlog) then

         do 500 j=1,nsp(1)

            ipoint=ntype(j)

            do 600 i=nsp(1)+1,num

               jpoint=ntype(i)

cc               drsq=dxsav(i,j)*dxsav(i,j)
cc     x             +dysav(i,j)*dysav(i,j)
cc     x             +dzsav(i,j)*dzsav(i,j)
c
c implement cut-off:
cc            if (drsq.ge.rsqmax) goto 40

               drsqrec=1.0d0/drsq
               dr=dsqrt(drsq)
               drrec=1.0d0/dr

cc               dx=dxsav(i,j)
cc               dy=dysav(i,j)
cc               dz=dzsav(i,j)

               repsilon=epsilonx(j)*dx+epsilony(j)*dy+epsilonz(j)*dz
c
c divide by dr to make r in r.epsilon a unit vector:
c
               repsilon=repsilon*drrec

               expft=ftb(ipoint,jpoint)
     x    *exp(-ftalp(ipoint,jpoint)*(dr-delta(j)-repsilon))
               expft2=ftb2(ipoint,jpoint)
     x    *exp(-ftbeta(ipoint,jpoint)*(dr-delta(j)-repsilon))
               expft3=ftb3(ipoint,jpoint)
     x    *exp(-ftgamma(ipoint,jpoint)*(dr-delta(j)-repsilon))

               expft4=extrab*exp(-extraalpha*dr)

               engft1(j)=engft1(j)+expft
               engft2(j)=engft2(j)+expft2
               engft3(j)=engft3(j)+expft3

c calculate eps force components for use in epsilonmv:
c
               dxunit=dx*drrec
               dyunit=dy*drrec
               dzunit=dz*drrec

               engft1dotx(j)=engft1dotx(j)+expft*dxunit
               engft1doty(j)=engft1doty(j)+expft*dyunit
               engft1dotz(j)=engft1dotz(j)+expft*dzunit
               engft2dotx(j)=engft2dotx(j)+expft2*dxunit
               engft2doty(j)=engft2doty(j)+expft2*dyunit
               engft2dotz(j)=engft2dotz(j)+expft2*dzunit
               engft3dotx(j)=engft3dotx(j)+expft3*dxunit
               engft3doty(j)=engft3doty(j)+expft3*dyunit
               engft3dotz(j)=engft3dotz(j)+expft3*dzunit

               sgam=ftalp(1,2)*expft+ftbeta(1,2)*expft2
     x               +ftgamma(1,2)*expft3
               sgam2x=drrec*epsilonx(j)-dx*drsqrec*drrec*repsilon
               sgam2y=drrec*epsilony(j)-dy*drsqrec*drrec*repsilon
               sgam2z=drrec*epsilonz(j)-dz*drsqrec*drrec*repsilon

               sgam3=extraalpha*drrec*expft4

               fttx=sgam*(dx*drrec-sgam2x)+(sgam3)*dx
               ftty=sgam*(dy*drrec-sgam2y)+(sgam3)*dy
               fttz=sgam*(dz*drrec-sgam2z)+(sgam3)*dz

               engsr=engsr+expft+expft2+expft3+expft4

c Stress tensor accumulation:
c short-range forces:
c
cc               stsrxx=stsrxx+fttx*dxsav(i,j)
cc               stsrxy=stsrxy+fttx*dysav(i,j)
cc               stsrxz=stsrxz+fttx*dzsav(i,j)
cc               stsryy=stsryy+ftty*dysav(i,j)
cc               stsryz=stsryz+ftty*dzsav(i,j)
cc               stsrzz=stsrzz+fttz*dzsav(i,j)

               frrx(j)=frrx(j)-fttx
               frry(j)=frry(j)-ftty
               frrz(j)=frrz(j)-fttz

               frrx(i)=frrx(i)+fttx
               frry(i)=frry(i)+ftty
               frrz(i)=frrz(i)+fttz

  600       continue

  500    continue

      endif

c QAIM energy + force calc.
      if (quaimlog) then

         do 510 j=1,nsp(1)

            ipoint=ntype(j)

            do 610 i=nsp(1)+1,num

               jpoint=ntype(i)

cc               drsq=dxsav(i,j)*dxsav(i,j)
cc    x             +dysav(i,j)*dysav(i,j)
cc     x             +dzsav(i,j)*dzsav(i,j)
c
c implement cut-off:
cc            if (drsq.ge.rsqmax) goto 40

               drsqrec=1.0d0/drsq
               dr=dsqrt(drsq)
               drrec=1.0d0/dr

cc               dx=dxsav(i,j)
cc               dy=dysav(i,j)
cc               dz=dzsav(i,j)

               repsilon=epsilonx(j)*dx+epsilony(j)*dy+epsilonz(j)*dz
c
c divide by dr to make r in r.epsilon a unit vector:
c
               repsilon=repsilon*drrec

c quadrupolar distortion code:
c
               dxx=dx*dx
               dyy=dy*dy
               dzz=dz*dz
               dxy=dx*dy
               dxz=dx*dz
               dyz=dy*dz
c
c txy, etc.
c
               txx=3.0d0*dxx*drsqrec-1.0d0
               tyy=3.0d0*dyy*drsqrec-1.0d0
               tzz=3.0d0*dzz*drsqrec-1.0d0
               txy=3.0d0*dxy*drsqrec
               txz=3.0d0*dxz*drsqrec
               tyz=3.0d0*dyz*drsqrec

               rrtheta=quaimxx(j)*txx+quaimyy(j)*tyy+quaimzz(j)*tzz+
     x              2.0d0*(quaimxy(j)*txy+quaimxz(j)*txz+quaimyz(j)*tyz)

               expft=ftb(ipoint,jpoint)
     x    *exp(-ftalp(ipoint,jpoint)*(dr-delta(j)-repsilon-rrtheta))
               expft2=ftb2(ipoint,jpoint)
     x    *exp(-ftbeta(ipoint,jpoint)*(dr-delta(j)-repsilon-rrtheta))
               expft3=ftb3(ipoint,jpoint)
     x    *exp(-ftgamma(ipoint,jpoint)*(dr-delta(j)-repsilon-rrtheta))

               expft4=extrab*exp(-extraalpha*dr)

               engft1(j)=engft1(j)+expft
               engft2(j)=engft2(j)+expft2
               engft3(j)=engft3(j)+expft3

c calculate eps force components for use in epsilonmv:
c
               dxunit=dx*drrec
               dyunit=dy*drrec
               dzunit=dz*drrec

               engft1dotx(j)=engft1dotx(j)+expft*dxunit
               engft1doty(j)=engft1doty(j)+expft*dyunit
               engft1dotz(j)=engft1dotz(j)+expft*dzunit
               engft2dotx(j)=engft2dotx(j)+expft2*dxunit
               engft2doty(j)=engft2doty(j)+expft2*dyunit
               engft2dotz(j)=engft2dotz(j)+expft2*dzunit
               engft3dotx(j)=engft3dotx(j)+expft3*dxunit
               engft3doty(j)=engft3doty(j)+expft3*dyunit
               engft3dotz(j)=engft3dotz(j)+expft3*dzunit

c calculate quaim derivatives for use in quaimmv:
c
               engftdotxx(j,1)=engftdotxx(j,1)+expft*txx
               engftdotyy(j,1)=engftdotyy(j,1)+expft*tyy
               engftdotzz(j,1)=engftdotzz(j,1)+expft*tzz
               engftdotxy(j,1)=engftdotxy(j,1)+2.0d0*expft*txy
               engftdotxz(j,1)=engftdotxz(j,1)+2.0d0*expft*txz
               engftdotyz(j,1)=engftdotyz(j,1)+2.0d0*expft*tyz

               engftdotxx(j,2)=engftdotxx(j,2)+expft2*txx
               engftdotyy(j,2)=engftdotyy(j,2)+expft2*tyy
               engftdotzz(j,2)=engftdotzz(j,2)+expft2*tzz
               engftdotxy(j,2)=engftdotxy(j,2)+2.0d0*expft2*txy
               engftdotxz(j,2)=engftdotxz(j,2)+2.0d0*expft2*txz
               engftdotyz(j,2)=engftdotyz(j,2)+2.0d0*expft2*tyz

               engftdotxx(j,3)=engftdotxx(j,3)+expft3*txx
               engftdotyy(j,3)=engftdotyy(j,3)+expft3*tyy
               engftdotzz(j,3)=engftdotzz(j,3)+expft3*tzz
               engftdotxy(j,3)=engftdotxy(j,3)+2.0d0*expft3*txy
               engftdotxz(j,3)=engftdotxz(j,3)+2.0d0*expft3*txz
               engftdotyz(j,3)=engftdotyz(j,3)+2.0d0*expft3*tyz

               sgam=ftalp(1,2)*expft+ftbeta(1,2)*expft2
     x               +ftgamma(1,2)*expft3
               sgam2x=drrec*epsilonx(j)-dx*drsqrec*drrec*repsilon
               sgam2y=drrec*epsilony(j)-dy*drsqrec*drrec*repsilon
               sgam2z=drrec*epsilonz(j)-dz*drsqrec*drrec*repsilon

c Quaim contributions:
c
                sgam2xx=quaimxx(j)
     x             *(6.0d0*dx*drsqrec-6.0d0*dxx*dx*dr4rec)-
     x              quaimyy(j)*6.0d0*dx*dyy*dr4rec-
     x              quaimzz(j)*6.0d0*dx*dzz*dr4rec+
     x              quaimxy(j)*(6.0d0*dy*drsqrec-12.0d0*dy*dxx*dr4rec)+
     x              quaimxz(j)*(6.0d0*dz*drsqrec-12.0d0*dz*dxx*dr4rec)-
     x              quaimyz(j)*12.0d0*dx*dy*dz*dr4rec

                sgam2yy=quaimyy(j)*(6.0d0*dy*drsqrec
     x             -6.0d0*dyy*dy*dr4rec)-
     x              quaimxx(j)*6.0d0*dy*dxx*dr4rec-
     x              quaimzz(j)*6.0d0*dy*dzz*dr4rec+
     x              quaimxy(j)*(6.0d0*dx*drsqrec-12.0d0*dx*dyy*dr4rec)+
     x              quaimyz(j)*(6.0d0*dz*drsqrec-12.0d0*dz*dyy*dr4rec)-
     x              quaimxz(j)*12.0d0*dx*dy*dz*dr4rec

                sgam2zz=quaimzz(j)*(6.0d0*dz*drsqrec
     x             -6.0d0*dzz*dz*dr4rec)-
     x              quaimxx(j)*6.0d0*dz*dxx*dr4rec-
     x              quaimyy(j)*6.0d0*dz*dyy*dr4rec+
     x              quaimxz(j)*(6.0d0*dx*drsqrec-12.0d0*dx*dzz*dr4rec)+
     x              quaimyz(j)*(6.0d0*dy*drsqrec-12.0d0*dy*dzz*dr4rec)-
     x              quaimxy(j)*12.0d0*dx*dy*dz*dr4rec

               sgam3=extraalpha*drrec*expft4

               fttx=sgam*(dx*drrec-sgam2x-sgam2xx)+(sgam3)*dx
               ftty=sgam*(dy*drrec-sgam2y-sgam2yy)+(sgam3)*dy
               fttz=sgam*(dz*drrec-sgam2z-sgam2zz)+(sgam3)*dz

               engsr=engsr+expft+expft2+expft3+expft4

               frrx(j)=frrx(j)-fttx
               frry(j)=frry(j)-ftty
               frrz(j)=frrz(j)-fttz

               frrx(i)=frrx(i)+fttx
               frry(i)=frry(i)+ftty
               frrz(i)=frrz(i)+fttz

c Stress tensor accumulation:
c short-range forces:
c
cc               stsrxx=stsrxx+fttx*dxsav(i,j)
cc               stsrxy=stsrxy+fttx*dysav(i,j)
cc               stsrxz=stsrxz+fttx*dzsav(i,j)
cc               stsryy=stsryy+ftty*dysav(i,j)
cc               stsryz=stsryz+ftty*dzsav(i,j)
cc               stsrzz=stsrzz+fttz*dzsav(i,j)

  610       continue

  510    continue

      endif

c anion-anion EPP
c TMP - CIM, AIM.
      if ((cimlog).or.(daimlog).or.(quaimlog)) then

         do 70 i=2,nanion

            ipoint=ntype(i)

            do 80 j=1,i-1

               jpoint=ntype(j)

cc               drsq=dxsav(i,j)*dxsav(i,j)
cc     x             +dysav(i,j)*dysav(i,j)
cc     x             +dzsav(i,j)*dzsav(i,j)
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

cc               fttx=sgam*dxsav(i,j)
cc               ftty=sgam*dysav(i,j)
cc               fttz=sgam*dzsav(i,j)

               engsr=engsr+expft

c
c Stress tensor accumulation:
c short-range forces:
c
cc               stsrxx=stsrxx+fttx*dxsav(i,j)
cc               stsrxy=stsrxy+fttx*dysav(i,j)
cc               stsrxz=stsrxz+fttx*dzsav(i,j)
cc               stsryy=stsryy+ftty*dysav(i,j)
cc               stsryz=stsryz+ftty*dzsav(i,j)
cc               stsrzz=stsrzz+fttz*dzsav(i,j)

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

cc               drsq=dxsav(i,j)*dxsav(i,j)
cc     x             +dysav(i,j)*dysav(i,j)
cc     x             +dzsav(i,j)*dzsav(i,j)
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

cc               fttx=sgam*dxsav(i,j)
cc               ftty=sgam*dysav(i,j)
cc               fttz=sgam*dzsav(i,j)

               engsr=engsr+expft
c
c Stress tensor accumulation:
c short-range forces:
c
cc               stsrxx=stsrxx+fttx*dxsav(i,j)
cc               stsrxy=stsrxy+fttx*dysav(i,j)
cc               stsrxz=stsrxz+fttx*dzsav(i,j)
cc               stsryy=stsryy+ftty*dysav(i,j)
cc               stsryz=stsryz+ftty*dzsav(i,j)
cc               stsrzz=stsrzz+fttz*dzsav(i,j)

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
       engtot_lj=0
c         write(*,*) "begin lj"
c         write(*,*) "engpetot", engpetot
c	 write(*,*) "lj eng",engtot_lj
c         write(*,*) "engsr", engsr
       if (.not.mixlog) then
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

       end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
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
cc               dxsav(i,j)=hlab2(1,1)*dxcf
cc     x                +hlab2(1,2)*dycf+hlab2(1,3)*dzcf
cc               dysav(i,j)=hlab2(2,1)*dxcf
cc     x                +hlab2(2,2)*dycf+hlab2(2,3)*dzcf
cc               dzsav(i,j)=hlab2(3,1)*dxcf
cc     x                +hlab2(3,2)*dycf+hlab2(3,3)*dzcf

cc               dxsav(j,i)=-dxsav(i,j)
cc               dysav(j,i)=-dysav(i,j)
cc               dzsav(j,i)=-dzsav(i,j)

cc               drsq=dxsav(i,j)*dxsav(i,j)
cc     x             +dysav(i,j)*dysav(i,j)
cc     x             +dzsav(i,j)*dzsav(i,j)

               dr=dsqrt(drsq)
cc               drsav(i,j)=dr
cc               drsav(j,i)=dr

               drrec=1.0d0/dr
               drsqrec=drrec*drrec
               dr6rec=drrec**6.0d0
               dr12rec=drrec**12.0d0

               ulj=fourepsil(ipoint,jpoint)
     x                     *((sig12(ipoint,jpoint)*dr12rec)
     x                     -(sig6(ipoint,jpoint)*dr6rec))
                if((j==775).and.(i==1890)) then
c               write(*,*) ulj 
c               write(*,*) fourepsil(ipoint,jpoint)
c     x              ,sig12(ipoint,jpoint)
c     x              ,sig6(ipoint,jpoint)
c     x               , ipoint, jpoint
c                if(ulj.le.(0-10)) then
c               write(*,*) i,j,ipoint, jpoint, ulj, "ulj"
c               write(*,*) dr
               end if
               gam=12.0d0*fourepsil(ipoint,jpoint)*sig12(ipoint,jpoint)
     x         *dr12rec*drsqrec
     x         -6.0d0*fourepsil(ipoint,jpoint)*sig6(ipoint,jpoint)
     x         *dr6rec*drsqrec

               engsr=engsr+ulj
               engtot_lj=engtot_lj+ulj
cc               frrx(j)=frrx(j)-dxsav(i,j)*gam
cc               frry(j)=frry(j)-dysav(i,j)*gam
cc               frrz(j)=frrz(j)-dzsav(i,j)*gam
cc               frrx(i)=frrx(i)+dxsav(i,j)*gam
cc               frry(i)=frry(i)+dysav(i,j)*gam
cc               frrz(i)=frrz(i)+dzsav(i,j)*gam

cc               stsrxx=stsrxx+dxsav(i,j)*gam*dxsav(i,j)
cc               stsrxy=stsrxy+dxsav(i,j)*gam*dysav(i,j)
cc               stsrxz=stsrxz+dxsav(i,j)*gam*dzsav(i,j)
cc               stsryy=stsryy+dysav(i,j)*gam*dysav(i,j)
cc               stsryz=stsryz+dysav(i,j)*gam*dzsav(i,j)
cc               stsrzz=stsrzz+dzsav(i,j)*gam*dzsav(i,j)

  310       continue

  305    continue
c        write(*,*) "end lj"
c        write(*,*) "engpetot", engpetot
c	 write(*,*) "lj eng",engtot_lj
c         write(*,*) "engsr", engsr

       endif

      if ((sw2potlog).or.(sw3potlog)) then

       if (.not.mixlog) then
         stsrxx=0.0d0
         stsrxy=0.0d0
         stsrxz=0.0d0
         stsryy=0.0d0
         stsryz=0.0d0
         stsrzz=0.0d0

         engpetot=0.0d0

         do 390 i=1,num

            frrx(i)=0.0d0
            frry(i)=0.0d0
            frrz(i)=0.0d0

  390    continue

        end if

         do 405 i=2,num

            ipoint=ntype(i)

            do 410 j=1,i-1

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
cc               dxsav(i,j)=hlab2(1,1)*dxcf
cc     x                +hlab2(1,2)*dycf+hlab2(1,3)*dzcf
cc               dysav(i,j)=hlab2(2,1)*dxcf
cc     x                +hlab2(2,2)*dycf+hlab2(2,3)*dzcf
cc               dzsav(i,j)=hlab2(3,1)*dxcf
cc     x                +hlab2(3,2)*dycf+hlab2(3,3)*dzcf

c               dxsav(j,i)=-dxsav(i,j)
cc               dysav(j,i)=-dysav(i,j)
cc               dzsav(j,i)=-dzsav(i,j)

cc               drsq=dxsav(i,j)*dxsav(i,j)
cc     x             +dysav(i,j)*dysav(i,j)
cc     x             +dzsav(i,j)*dzsav(i,j)

               dr=dsqrt(drsq)
cc               drsav(i,j)=dr
cc               drsav(j,i)=dr

c reduced separation.
               drsig=dr/swsigma(ipoint,jpoint)

               if (drsig.gt.swa(ipoint,jpoint)) goto 410

               swexp=1.0d0/(drsig-swa(ipoint,jpoint))
               t1=exp(swexp)

               t2=swbigA(ipoint,jpoint)
     x           *((swB(ipoint,jpoint)*drsig**(-swp(ipoint,jpoint)))
     x                               -(drsig**(-swq(ipoint,jpoint))))

               usw=swepsilon(ipoint,jpoint)*t1*t2

c now the forces.
               term1=swbigA(ipoint,jpoint)*((-swp(ipoint,jpoint)
     x  *swB(ipoint,jpoint)*(drsig**(-(swp(ipoint,jpoint)+1.0d0))))
     x +(swq(ipoint,jpoint)*(drsig**(-(swq(ipoint,jpoint)+1.0d0)))))
     x              *t1
               term2=t1*t2*swexp*swexp

               force=-term1+term2

               gam=swepsilon(ipoint,jpoint)*force
     x            /swsigma(ipoint,jpoint)/dr

               engsr=engsr+usw

cc               frrx(j)=frrx(j)-dxsav(i,j)*gam
cc               frry(j)=frry(j)-dysav(i,j)*gam
c               frrz(j)=frrz(j)-dzsav(i,j)*gam
c               frrx(i)=frrx(i)+dxsav(i,j)*gam
c               frry(i)=frry(i)+dysav(i,j)*gam
c               frrz(i)=frrz(i)+dzsav(i,j)*gam

c               stsrxx=stsrxx+dxsav(i,j)*gam*dxsav(i,j)
c               stsrxy=stsrxy+dxsav(i,j)*gam*dysav(i,j)
c               stsrxz=stsrxz+dxsav(i,j)*gam*dzsav(i,j)
c               stsryy=stsryy+dysav(i,j)*gam*dysav(i,j)
c               stsryz=stsryz+dysav(i,j)*gam*dzsav(i,j)
c               stsrzz=stsrzz+dzsav(i,j)*gam*dzsav(i,j)

  410       continue

  405    continue

      endif

c Stillinger-Weber 3-body term.
      if ((sw3potlog)) then

         engsw3=0.0d0
         n3=0
         do 700 i=1,num

            ipoint=ntype(i)

            do 710 j=2,num

               jpoint=ntype(j)

               if (i.eq.j) goto 710

c               drsigij=drsav(i,j)/swsigma(ipoint,jpoint)
               if (drsigij.gt.swa(ipoint,jpoint)) goto 710

               do 720 k=1,j

                  kpoint=ntype(k)

                  if ((i.eq.k).or.(j.eq.k)) goto 720

c                  drsigik=drsav(i,k)/swsigma(ipoint,jpoint)
                  if (drsigik.gt.swa(ipoint,kpoint)) goto 720

                  n3=n3+1
c calculate energy.
c                  dxij=dxsav(i,j)/swsigma(ipoint,jpoint)
c                  dyij=dysav(i,j)/swsigma(ipoint,jpoint)
c                  dzij=dzsav(i,j)/swsigma(ipoint,jpoint)

                  drsq=(dxij*dxij)+(dyij*dyij)+(dzij*dzij)
                  drij=dsqrt(drsq)

c                  dxik=dxsav(i,k)/swsigma(ipoint,jpoint)
c                  dyik=dysav(i,k)/swsigma(ipoint,jpoint)
c                  dzik=dzsav(i,k)/swsigma(ipoint,jpoint)

                  drsq=(dxik*dxik)+(dyik*dyik)+(dzik*dzik)
                  drik=dsqrt(drsq)

                  cstheta=(dxij*dxik+dyij*dyik+dzij*dzik)/(drik*drij)

                  angle=acos(cstheta)

c h1 is for theta_jik
                  h1a=swlambda(ipoint,jpoint)
     x     *dexp((swgamma(ipoint,jpoint)/(drij-swa(ipoint,jpoint)))
     x         +(swgamma(ipoint,kpoint)/(drik-swa(ipoint,kpoint))))

                  term1b=(cstheta+(1.0d0/3.0d0))**2.0d0

                  h1=h1a*term1b

                  engsr=engsr+(swepsilon(ipoint,jpoint)*h1)
                  engsw3=engsw3+(swepsilon(ipoint,jpoint)*h1)


c now attempt to work out the forces.
c dh1/drij
cc                  term1a=(-swlambda(ipoint,jpoint)
cc     x     *swgamma(ipoint,jpoint)*dexp(swgamma(ipoint,jpoint)
cc     x     /(drij-swa(ipoint,jpoint))))
cc     x     /((drij-swa(ipoint,jpoint))**2.0d0)
                  term1a=-swgamma(ipoint,jpoint)*h1a
     x     /((drij-swa(ipoint,jpoint))**2.0d0)

c dh1/drik
cc                  term1a2=(-swlambda(ipoint,kpoint)
cc     x     *swgamma(ipoint,kpoint)*dexp(swgamma(ipoint,kpoint)
cc     x     /(drik-swa(ipoint,kpoint))))
cc     x     /((drik-swa(ipoint,kpoint))**2.0d0)
                  term1a2=-swgamma(ipoint,kpoint)*h1a
     x     /((drik-swa(ipoint,kpoint))**2.0d0)

                  term1arijx=term1a*dxij/drij
                  term1arijy=term1a*dyij/drij
                  term1arijz=term1a*dzij/drij

                  term1arikx=term1a2*dxik/drik
                  term1ariky=term1a2*dyik/drik
                  term1arikz=term1a2*dzik/drik

                  term1rijx=term1arijx*term1b
                  term1rijy=term1arijy*term1b
                  term1rijz=term1arijz*term1b

                  term1rikx=term1arikx*term1b
                  term1riky=term1ariky*term1b
                  term1rikz=term1arikz*term1b

                  u=dxij*dxik+dyij*dyik+dzij*dzik
                  v=1.0d0/drij
                  w=1.0d0/drik

                  dudrijx=dxik
                  dudrijy=dyik
                  dudrijz=dzik
                  dudrikx=dxij
                  dudriky=dyij
                  dudrikz=dzij

                  dvdrijx=-dxij/(drij**3.0d0)
                  dvdrijy=-dyij/(drij**3.0d0)
                  dvdrijz=-dzij/(drij**3.0d0)
                  dwdrikx=-dxik/(drik**3.0d0)
                  dwdriky=-dyik/(drik**3.0d0)
                  dwdrikz=-dzik/(drik**3.0d0)

                  dcosthetajikdrijx=(u*w*dvdrijx)+(v*w*dudrijx)
                  dcosthetajikdrijy=(u*w*dvdrijy)+(v*w*dudrijy)
                  dcosthetajikdrijz=(u*w*dvdrijz)+(v*w*dudrijz)
                  dcosthetajikdrikx=(u*v*dwdrikx)+(v*w*dudrikx)
                  dcosthetajikdriky=(u*v*dwdriky)+(v*w*dudriky)
                  dcosthetajikdrikz=(u*v*dwdrikz)+(v*w*dudrikz)

c derivatives of the (cos(theta)+0.333)^2 term:
                  dfdrijx=2.0d0*(cstheta+0.3333333333333)
     x               *dcosthetajikdrijx
                  dfdrijy=2.0d0*(cstheta+0.3333333333333)
     x               *dcosthetajikdrijy
                  dfdrijz=2.0d0*(cstheta+0.3333333333333)
     x               *dcosthetajikdrijz
                  dfdrikx=2.0d0*(cstheta+0.3333333333333)
     x               *dcosthetajikdrikx
                  dfdriky=2.0d0*(cstheta+0.3333333333333)
     x               *dcosthetajikdriky
                  dfdrikz=2.0d0*(cstheta+0.3333333333333)
     x               *dcosthetajikdrikz

c put together into the complete second term:
                  term2a=swlambda(ipoint,jpoint)
     x       *dexp((swgamma(ipoint,jpoint)
     x       /(drij-swa(ipoint,jpoint)))
     x            +(swgamma(ipoint,kpoint)
     x       /(drik-swa(ipoint,kpoint))))

                  term2rijx=term2a*dfdrijx
                  term2rijy=term2a*dfdrijy
                  term2rijz=term2a*dfdrijz
                  term2rikx=term2a*dfdrikx
                  term2riky=term2a*dfdriky
                  term2rikz=term2a*dfdrikz

                  fulldxij=term1rijx+term2rijx
                  fulldyij=term1rijy+term2rijy
                  fulldzij=term1rijz+term2rijz
                  fulldxik=term1rikx+term2rikx
                  fulldyik=term1riky+term2riky
                  fulldzik=term1rikz+term2rikz

                  frrxij=(fulldxij*swepsilon(ipoint,jpoint)
     x            /swsigma(ipoint,jpoint))
                  frryij=(fulldyij*swepsilon(ipoint,jpoint)
     x            /swsigma(ipoint,jpoint))
                  frrzij=(fulldzij*swepsilon(ipoint,jpoint)
     x            /swsigma(ipoint,jpoint))

                  frrx(j)=frrx(j)+frrxij
                  frry(j)=frry(j)+frryij
                  frrz(j)=frrz(j)+frrzij
                  frrx(i)=frrx(i)-frrxij
                  frry(i)=frry(i)-frryij
                  frrz(i)=frrz(i)-frrzij

                  frrxik=(fulldxik*swepsilon(ipoint,kpoint)
     x            /swsigma(ipoint,jpoint))
                  frryik=(fulldyik*swepsilon(ipoint,kpoint)
     x            /swsigma(ipoint,jpoint))
                  frrzik=(fulldzik*swepsilon(ipoint,kpoint)
     x            /swsigma(ipoint,jpoint))

                  frrx(k)=frrx(k)+frrxik
                  frry(k)=frry(k)+frryik
                  frrz(k)=frrz(k)+frrzik
                  frrx(i)=frrx(i)-frrxik
                  frry(i)=frry(i)-frryik
                  frrz(i)=frrz(i)-frrzik

c                  stsrxx=stsrxx-dxsav(i,j)*frrxij
c                  stsrxy=stsrxy-dxsav(i,j)*frryij
c                  stsrxz=stsrxz-dxsav(i,j)*frrzij
c                  stsryy=stsryy-dysav(i,j)*frryij
c                  stsryz=stsryz-dysav(i,j)*frrzij
c                  stsrzz=stsrzz-dzsav(i,j)*frrzij
c                  stsrxx=stsrxx-dxsav(i,k)*frrxik
c                  stsrxy=stsrxy-dxsav(i,k)*frryik
c                  stsrxz=stsrxz-dxsav(i,k)*frrzik
c                  stsryy=stsryy-dysav(i,k)*frryik
c                  stsryz=stsryz-dysav(i,k)*frrzik
c                  stsrzz=stsrzz-dzsav(i,k)*frrzik
c
  720          continue

  710       continue

  700    continue

      endif

      if ((morselog)) then

         engmorse=0.0d0
         do 800 i=2,nanion

            do 810 j=1,i-1

c               drsq=dxsav(i,j)*dxsav(i,j)
c     x             +dysav(i,j)*dysav(i,j)
c     x             +dzsav(i,j)*dzsav(i,j)

               dr=dsqrt(drsq)

               v=demorse*((1.0d0-exp(-amorse*(dr-remorse)))**2.0d0)

               dvdr=2.0d0*demorse*(1.0d0-exp(-amorse*(dr-remorse)))
     x          *amorse*exp(-amorse*(dr-remorse))

c               frrxi=dvdr*dxsav(i,j)/dr
c               frryi=dvdr*dysav(i,j)/dr
c               frrzi=dvdr*dzsav(i,j)/dr

            frrx(i)=frrx(i)-frrxi
            frry(i)=frry(i)-frryi
            frrz(i)=frrz(i)-frrzi
            frrx(j)=frrx(j)+frrxi
            frry(j)=frry(j)+frryi
            frrz(j)=frrz(j)+frrzi

            engmorse=engmorse+v-demorse
            engsr=engsr+v-demorse

  810       continue
  800    continue

      endif

cc      write(88,*)engpetot,engsr,engsw3,n3
cc      write(88,*)engmorse

!!!!!!!!!!!!!!!!!!!!!!!!!!! TERSOFF II !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


          if ((tf2potlog))then 
cc            call separations
         pi=4.0d0*atan(1.0d0)

           if (.not.mixlog) then

         stsrxx=0.0d0
         stsrxy=0.0d0
         stsrxz=0.0d0
         stsryy=0.0d0
         stsryz=0.0d0
         stsrzz=0.0d0
 
         engpetot=0.0d0

          do 890 i=1,num
            frrx(i)=0.0d0
            frry(i)=0.0d0
            frrz(i)=0.0d0
  890     continue

           end if
c         write(*,*) "begin tf"
c         write(*,*) "engpetot", engpetot
c         write(*,*) "engsr", engsr

             tf_xi=0
             engTF2=0.0d0
             engTF2_2=0.0d0
             engTF2_3=0.0d0
             engdisp=0.0d0

cc         write(6,*)'nnnnnnnnnnnnnnnnnnn',num
          do 900 i=1,num           !!!!!!!loop over all atoms

               ipoint=ntype(i)
cc         write(6,*)'nnnnnnnnnnnnnnnnnnn',ntype(i)

          test_interaction=bigA_type(ipoint,ipoint)
     x                    +bigB_type(ipoint,ipoint)


          if(test_interaction==0)goto 900   !!!!!!!!not a carbon if A and B are zero
!!!!!!!!let the parameters be for carbon-carbon!!!!!!!!!!!!!!!!!!!!!!!!!!

         bigA      =bigA_type(ipoint,ipoint)
         bigB      =bigB_type(ipoint,ipoint)
         tf_lamda1 =tf_lamda1_type(ipoint,ipoint)
         tf_lamda2 =tf_lamda2_type(ipoint,ipoint)
         tf_lamda3 =tf_lamda3_type(ipoint,ipoint)
         tf_beta   =tf_beta_type(ipoint,ipoint)
         tf_c      =tf_c_type(ipoint,ipoint)
         tf_d      =tf_d_type(ipoint,ipoint)
         tf_h      =tf_h_type(ipoint,ipoint)
         tf_n      =tf_n_type(ipoint,ipoint)
         bigR      =bigR_type(ipoint,ipoint)
         bigD      =bigD_type(ipoint,ipoint)
c       write(*,*) bigA, bigB ,tf_lamda1, tf_lamda2, tf_lamda3
c       write(*,*) tf_beta, tf_c, tf_d, tf_h, tf_n, bigR, bigD

c        write(*,*)test_interaction


             do 910 j=1,num                 !!!!loop over all pairs 
cc         write(6,*)'nnnnnnnnnnnnnnnnnnn',ntype(j)
cc          if(ntype(i).ne.ntype(j)) goto 910 !!!!!!there is no interaction, i is carbon but j isn't

c           write(*,*)ntype(i),ntype(j)
               jpoint=ntype(j)

                if (i.eq.j) goto 910         !!!!!!!if i and j are the same, not a pair

              tf_xi_sum=0 ! wipe running total
              dxi_drijx=0.0d0 
              dxi_drijy=0.0d0 
              dxi_drijz=0.0d0 
              dxi_drikx=0.0d0 
              dxi_driky=0.0d0 
              dxi_drikz=0.0d0 

c              distance between i and j 

c                 drij=(dxsav(i,j)**2.0d0)
c     x                +(dysav(i,j)**2.0d0)
c     x                +(dzsav(i,j)**2.0d0)
                 drij=drij**(0.5d0)       !dist between i and j

ccccc DISPERSION - FEB. 2011... REALE... NOT CALLED FOR TERSOFF...
            dr=drij
            drsq=dr*dr
            drsqrec=1.0d0/(drij*drij)
            drsq3rec=drsqrec*drsqrec*drsqrec
            drrec=1.0d0/drij
            dr3rec=drrec*drsqrec

cc            sgam=-(6.0d0*ftc(ipoint,jpoint)
cc     x          +8.0d0*ftd(ipoint,jpoint)*drsqrec)
cc     x          *drsq3rec*drsqrec

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

            sgam=(
     x          -(6.0d0*ftc(ipoint,jpoint)*f6
     x          -(ftc(ipoint,jpoint)*dddampforce*dr)
     x          +8.0d0*ftd(ipoint,jpoint)*f8*drsqrec
     x          -(ftd(ipoint,jpoint)*dqdampforce*drrec))
     x          *drsq3rec)*drsqrec

c            fttx=sgam*dxsav(i,j)
c            ftty=sgam*dysav(i,j)
c            fttz=sgam*dzsav(i,j)

            cpe=-(ftc(ipoint,jpoint)*f6)*drsq3rec
            engsr=engsr+cpe

            frrx(j)=frrx(j)-fttx
            frry(j)=frry(j)-ftty
            frrz(j)=frrz(j)-fttz

            frrx(i)=frrx(i)+fttx
            frry(i)=frry(i)+ftty
            frrz(i)=frrz(i)+fttz

cc            write(200,*)i,j,fttx,ftty,fttz
cccccc end 2011



c              fc(rij) 
                 If (drij.le.(bigR-bigD))then
                 cutoffij=1.0d0
                 else if (drij.le.(bigR+bigD)) then
                 cutoffij=0.5d0*pi*(drij-bigR)/bigD
                 cutoffij=0.5d0-(0.5d0*dsin(cutoffij))
                  else 
                 cutoffij=0.0d0
                 goto 910
                 Endif


c               fR(rij)=repulse and Fa(rij)=attract
 
                repulse=bigA*dexp(-1.0d0*tf_lamda1*drij)
                attract=-1.0d0*bigB*dexp(-1.0d0*tf_lamda2*drij)

                                            
c                generate all triplets by do look over all atoms k
c               calculate xi(i,j) which is a sum over all k triplets


                do 920 k=1,num

               kpoint=ntype(k)

cc          if(ntype(k).ne.ntype(i)) goto 920 !there is no interaction, i is carbon but k isn't

!!!!!!!!!!!!!!!!!!ENERGY TERMS - depends on k!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if ((i.eq.k).or.(j.eq.k)) goto 920
c                  drik=(dxsav(i,k)**2.0d0)
c     x                +(dysav(i,k)**2.0d0)
c     x                +(dzsav(i,k)**2.0d0)
                 drik=drik**(0.5d0)      !distance between i and k 

c              fc(rik)


                If (drik.le.(bigR-bigD))then
                 cutoffik=1.0d0
                 else if (drik.le.(bigR+bigD)) then
                 cutoffik=0.5d0*pi*(drik-bigR)/bigD
                 cutoffik=0.5d0-(0.5d0*dsin(cutoffik))
                 else 
                      cutoffik=0.0d0
                      goto 920
                 Endif

c               g(theta)


c                cstheta= (dxsav(i,j)*dxsav(i,k))
c     x                  +(dysav(i,j)*dysav(i,k))
c     x                  +(dzsav(i,j)*dzsav(i,k))
                cstheta=cstheta/(drij*drik)

                goftheta=1+((tf_c**2.0d0)/(tf_d**2.0d0))
     x                  -( (tf_c**2.0d0)/((tf_d**2.0d0)
     x                  +((tf_h-cstheta)**2.0d0)) )

c               xi(i,j)
                
                tf_xi=(cutoffik*goftheta*
     x           dexp((tf_lamda3**3.0d0)*((drij-drik)**3.0d0)))
                tf_xi_sum=tf_xi_sum+tf_xi
c                 write(*,*) tf_xi_sum, tf_xi , "XI"
  920          continue     ! ends 1st summations over triplets 

!!!!!!!!!!!!!!!!!CALCULATE REMAINING ENERGY TERMS!!!!!!!!!!!!!!!!!!!

c          calculate bij

              bij=(1.0d0+( (tf_beta**tf_n)
     x           *(tf_xi_sum**tf_n) ))
     x           **(-1.0d0/(2.0d0*tf_n))
c          calculate the energy
              engsr=engsr+(0.5d0*( cutoffij*(repulse+(bij*attract)) ))  
c               write(*,*)cutoffij, repulse, bij, attract 
              engTF2=engTF2+(0.5d0*( cutoffij*(repulse+(bij*attract)) ))
              engTF2_2=engTF2_2+(0.5d0*(cutoffij*(repulse)))
              engTF2_3=engTF2_3+(0.5d0*(cutoffij*(bij*attract)))
cc               write(*,*)i,j, engTF2
!!!!!!!!!!!!!!!!!!!!CALCULATE U_b (depends on k) FORCE TERMS w.r.t R(i,j)  !!!!!!

                do 903 k=1,num  !!!!!!!!!!!!!2nd loop over triplets

            kpoint=ntype(k)

cc          if(ntype(i).ne.ntype(k)) goto 903 !there is no interaction, i is carbon but k isn't


                if ((i.eq.k).or.(j.eq.k)) goto 903       !same atom, so not a triplet

c                  drik=(dxsav(i,k)**2.0d0)
c     x                +(dysav(i,k)**2.0d0)
c     x                +(dzsav(i,k)**2.0d0)
                 drik=drik**(0.5d0)      !distance between i and k

             If (drik.le.(bigR-bigD))then
                 cutoffik=1.0d0
                 else if (drik.le.(bigR+bigD)) then
                 cutoffik=0.5d0*pi*(drik-bigR)/bigD
                 cutoffik=0.5d0-(0.5d0*dsin(cutoffik))
                 else
                      cutoffik=0.0d0
                      goto 903
                 Endif

c          rij_dot_rik=  (dxsav(i,j)*dxsav(i,k))
c     x                  +(dysav(i,j)*dysav(i,k))
c     x                  +(dzsav(i,j)*dzsav(i,k))
c                cstheta= (dxsav(i,j)*dxsav(i,k))
c     x                  +(dysav(i,j)*dysav(i,k))
c     x                  +(dzsav(i,j)*dzsav(i,k))
                cstheta=cstheta/(drij*drik)
          goftheta=1+((tf_c**2.0d0)/(tf_d**2.0d0))
     x                  -( (tf_c**2.0d0)/((tf_d**2.0d0)
     x                  +((tf_h-cstheta)**2.0d0)) )


c           dcstheta_drijx=(dxsav(i,k)/(drij*drik))
c     x             -((dxsav(i,j)/drik)*(rij_dot_rik/(drij**3.d0)))


           dcstheta_drijy=(dysav(i,k)/(drij*drik))
     x          -((dysav(i,j)/drik)*(rij_dot_rik/(drij**3.d0)))



           dcstheta_drijz=(dzsav(i,k)/(drij*drik))
     x         -((dzsav(i,j)/drik)*(rij_dot_rik/(drij**3.d0)))



              !  dg(theta)/dr(i,j)

                dgoftheta_drijx=2.0d0*(tf_c**2.0d0)
     x          *(tf_h-cstheta)
     x          *(((tf_d**2.0d0)+((tf_h-cstheta)**2.0d0))**-2.0d0)
     x          *(-1.0d0*dcstheta_drijx)


                dgoftheta_drijy=2.0d0*(tf_c**2.0d0)
     x          *(tf_h-cstheta)
     x          *(((tf_d**2.0d0)+((tf_h-cstheta)**2.0d0))**-2.0d0)
     x          *(-1.0d0*dcstheta_drijy)

                dgoftheta_drijz=2.0d0*(tf_c**2.0d0)
     x          *(tf_h-cstheta)
     x          *(((tf_d**2.0d0)+((tf_h-cstheta)**2.0d0))**-2.0d0)
     x          *(-1.0d0*dcstheta_drijz)

       !          d xi/dr(i,j)

c                dxi_drijx=
c     x              (cutoffik*dgoftheta_drijx
c     x              *dexp((tf_lamda3**3.0d0)*((drij-drik)**3.0d0)))
c     x              +(cutoffik*goftheta*
c     x              dexp((tf_lamda3**3.0d0)*((drij-drik)**3.0d0)) 
c     x              *(tf_lamda3**3.0d0)*((drij-drik)**2.0d0)
c     x              *3.0d0*(dxsav(i,j)/drij))

c              dxi_drijy=
c     x             (cutoffik*dgoftheta_drijy
c     x             *dexp((tf_lamda3**3.0d0)*((drij-drik)**3.0d0)))
c     x             +(cutoffik*goftheta*
c     x             dexp((tf_lamda3**3.0d0)*((drij-drik)**3.0d0)) 
c     x             *(tf_lamda3**3.0d0)*((drij-drik)**2.0d0)
c     x             *3.0d0*(dysav(i,j)/drij))

c               dxi_drijz=
c     x            (cutoffik*dgoftheta_drijz
c     x            *dexp((tf_lamda3**3.0d0)*((drij-drik)**3.0d0)))
c     x            +(cutoffik*goftheta*
c     x            dexp((tf_lamda3**3.0d0)*((drij-drik)**3.0d0)) 
c     x            *(tf_lamda3**3.0d0)*((drij-drik)**2.0d0)
c     x             *3.0d0*(dzsav(i,j)/drij))


c  !        d(b(i,j))/dr(ij)
          if(tf_xi_sum==0)then 
          dbij_drijx=0.0d0
          dbij_drijy=0.0d0
          dbij_drijz=0.0d0
          else 
              dbij_drijx= (-1.0d0/(2.0d0*tf_n))*
     x        ( (1.0d0+((tf_beta**tf_n)
     x           *(tf_xi_sum**tf_n)))
     x        **((-1.0d0-(2.0d0*tf_n))/(2.0d0*tf_n)) )
     x        * (tf_beta**tf_n)
     x        * (tf_n*(tf_xi_sum**(tf_n-1.0d0)))
     x        * dxi_drijx


        dbij_drijy= (-1.0d0/(2.0d0*tf_n))*      !-1/2n
     x        ( (1.0d0+((tf_beta**tf_n)      
     x           *(tf_xi_sum**tf_n)))       
     x        **((-1.0d0-(2.0d0*tf_n))/(2.0d0*tf_n)) )
     x        * (tf_beta**tf_n)
     x        * (tf_n*(tf_xi_sum**(tf_n-1.0d0)))
     x        * dxi_drijy

        dbij_drijz= (-1.0d0/(2.0d0*tf_n))*
     x        ( (1.0d0+((tf_beta**tf_n)
     x           *(tf_xi_sum**tf_n)))
     x        **((-1.0d0-(2.0d0*tf_n))/(2.0d0*tf_n)) )
     x        * (tf_beta**tf_n)
     x        * (tf_n*(tf_xi_sum**(tf_n-1.0d0)))
     x        * dxi_drijz

           endif


               dU_drijx_b=cutoffij*dbij_drijx*attract
               dU_drijy_b=cutoffij*dbij_drijy*attract
               dU_drijz_b=cutoffij*dbij_drijz*attract


c update forces
               frrx(j)=frrx(j)+(0.5d0*dU_drijx_b)
               frry(j)=frry(j)+(0.5d0*dU_drijy_b)
               frrz(j)=frrz(j)+(0.5d0*dU_drijz_b)
               frrx(i)=frrx(i)-(0.5d0*dU_drijx_b)
               frry(i)=frry(i)-(0.5d0*dU_drijy_b)
               frrz(i)=frrz(i)-(0.5d0*dU_drijz_b)
c               stsrxx=stsrxx-dxsav(i,j)*(0.5d0*dU_drijx_b)
c               stsrxy=stsrxy-dxsav(i,j)*(0.5d0*dU_drijy_b)
c               stsrxz=stsrxz-dxsav(i,j)*(0.5d0*dU_drijz_b)
c               stsryy=stsryy-dysav(i,j)*(0.5d0*dU_drijy_b)
c               stsryz=stsryz-dysav(i,j)*(0.5d0*dU_drijz_b)
c               stsrzz=stsrzz-dzsav(i,j)*(0.5d0*dU_drijz_b)
c





!!!!!!!!!!!!!!!!!!!!CALCULATE REMAINING FORCE TERMS w.r.t R(i,k)  !!!!!!

c          d(b(i,j))/dr(i,k)

c!!!!!!!!!!!!!!!!!!! dU/dr(i,k)!!!!!!!!!!!!!!!!!!!!!!!!!!


c           d cos(theta)/dr(i,k)

c           dcstheta_drikx=(dxsav(i,j)/(drij*drik))
c     x           -((dxsav(i,k)/drij)*(rij_dot_rik/(drik**3.0d0)))

           dcstheta_driky=(dysav(i,j)/(drij*drik))
     x           -((dysav(i,k)/drij)*(rij_dot_rik/(drik**3.0d0)))

           dcstheta_drikz=(dzsav(i,j)/(drij*drik))
     x           -((dzsav(i,k)/drij)*(rij_dot_rik/(drik**3.0d0)))

c       !           dg(theta)/dr(i,k)

               dgoftheta_drikx=2.0d0*(tf_c**2.0d0)
     x          *(tf_h-cstheta)
     x          *(((tf_d**2.0d0)+((tf_h-cstheta)**2.0d0))**-2.0d0)
     x          *(-1.0d0*dcstheta_drikx)

               dgoftheta_driky=2.0d0*(tf_c**2.0d0)
     x          *(tf_h-cstheta)
     x          *(((tf_d**2.0d0)+((tf_h-cstheta)**2.0d0))**-2.0d0)
     x          *(-1.0d0*dcstheta_driky)

               dgoftheta_drikz=2.0d0*(tf_c**2.0d0)
     x          *(tf_h-cstheta)
     x          *(((tf_d**2.0d0)+((tf_h-cstheta)**2.0d0))**-2.0d0)
     x          *(-1.0d0*dcstheta_drikz)


c                 d fc(rik)/drik

               If (drik.le.(bigR-bigD))then
                 dcutoffik_drikx=0.0d0
                 else if (drik.le.(bigR+bigD)) then
                 dcutoffik_drikx_t1=0.5d0*pi*(drik-bigR)/bigD
c                 dcutoffik_drikx=-1.0d0*(dcos(dcutoffik_drikx_t1))
c     x                          *pi*dxsav(i,k)/(4.0d0*drik*bigD)
                 else
                   dcutoffik_drikx=0.0d0
                 Endif
c               write(*,*) "DDDik", nstep, dcutoffik_drikx, i, k,"x"


               If (drik.le.(bigR-bigD))then
                 dcutoffik_driky=0.0d0
                 else if (drik.le.(bigR+bigD)) then
                 dcutoffik_driky=0.5d0*pi*(drik-bigR)/bigD
                 dcutoffik_driky=-1.0d0*(dcos(dcutoffik_driky))
     x                          *pi*dysav(i,k)/(4.0d0*drik*bigD)
                 else
                    dcutoffik_driky=0.0d0
                 Endif

c               write(*,*) "DDDik", nstep, dcutoffik_driky, i, k,"y"

               If (drik.le.(bigR-bigD))then
                 dcutoffik_drikz=0.0d0
                 else if (drik.le.(bigR+bigD)) then
                 dcutoffik_drikz=0.5d0*pi*(drik-bigR)/bigD
                 dcutoffik_drikz=-1.d0*dcos((dcutoffik_drikz))
     x                          *pi*dzsav(i,k)/(4.0d0*drik*bigD)
                 else
                 dcutoffik_drikz=0.0d0
                 Endif

c              write(*,*) "DDDik", nstep, dcutoffik_drikz, i, k,"z"

c                 d xi/dr(i,k)
c              dxi_drikx=( dcutoffik_drikx*goftheta*
c     x        dexp((tf_lamda3**3.0d0)*((drij-drik)**3.0d0)) )
c     x        +(cutoffik*dgoftheta_drikx*
c     x        dexp((tf_lamda3**3.0d0)*((drij-drik)**3.0d0)) )
c     x         -(cutoffik*goftheta*
c     x         (dexp((tf_lamda3**3.0d0)*((drij-drik)**3.0d0)))
c     x         *(tf_lamda3**3.00)*(dxsav(i,k)/drik)
c     x        *3.0d0*((drij-drik)**2.0d0) )

c                dxi_driky=( dcutoffik_driky*goftheta*
c     x           dexp((tf_lamda3**3.0d0)*((drij-drik)**3.0d0)) )
c     x            +(cutoffik*dgoftheta_driky*
c     x            dexp((tf_lamda3**3.0d0)*((drij-drik)**3.0d0)) )
c     x            -(cutoffik*goftheta*
c     x            (dexp((tf_lamda3**3.0d0)*((drij-drik)**3.0d0)))
c     x            *(tf_lamda3**3.0d0)*(dysav(i,k)/drik)
c     x            *3.0d0*((drij-drik)**2.0d0) )

c                dxi_drikz=( dcutoffik_drikz*goftheta*
c     x           dexp((tf_lamda3**3.0d0)*((drij-drik)**3.0d0)) )
c     x           +(cutoffik*dgoftheta_drikz*
c     x           dexp((tf_lamda3**3.0d0)*((drij-drik)**3.0d0)) )
c     x           -(cutoffik*goftheta*
c     x           (dexp((tf_lamda3**3.0d0)*((drij-drik)**3.0d0)))
c     x           *(tf_lamda3**3.0d0)*(dzsav(i,k)/drik)
c     x           *3.0d0*((drij-drik)**2.0d0) )


c          d(b(i,j))/dr(i,k)
         if(tf_xi_sum==0)then
          dbij_drikx=0.0d0
          dbij_driky=0.0d0
          dbij_drikz=0.0d0
          else

            dbij_drikx=(-1.0d0/(2.0d0*tf_n))*
     x        ( (1.0d0+((tf_beta**tf_n)
     x           *(tf_xi_sum**tf_n)))
     x        **((-1.0d0-(2.0d0*tf_n))/(2.0d0*tf_n)) )
     x        * (tf_beta**tf_n)
     x        * (tf_n*(tf_xi_sum**(tf_n-1.0d0)))
     x        * dxi_drikx

            dbij_driky=(-1.0d0/(2.0d0*tf_n))*
     x        ( (1.0d0+((tf_beta**tf_n)
     x           *(tf_xi_sum**tf_n)))
     x        **((-1.0d0-(2.0d0*tf_n))/(2.0d0*tf_n)) )
     x        * (tf_beta**tf_n)
     x        * (tf_n*(tf_xi_sum**(tf_n-1.0d0)))
     x        * dxi_driky

            dbij_drikz=(-1.0d0/(2.0d0*tf_n))*
     x        ( (1.0d0+((tf_beta**tf_n)
     x           *(tf_xi_sum**tf_n)))
     x        **((-1.0d0-(2.0d0*tf_n))/(2.0d0*tf_n)) )
     x        * (tf_beta**tf_n)
     x        * (tf_n*(tf_xi_sum**(tf_n-1.0d0)))
     x        * dxi_drikz
          end if




c differential
               dU_drikx=cutoffij*attract*dbij_drikx
               dU_driky=cutoffij*attract*dbij_driky
               dU_drikz=cutoffij*attract*dbij_drikz
c update forces
               frrx(k)=frrx(k)+(0.5d0*dU_drikx)
               frry(k)=frry(k)+(0.5d0*dU_driky)
              frrz(k)=frrz(k)+(0.5d0*dU_drikz)
               frrx(i)=frrx(i)-(0.5d0*dU_drikx)
               frry(i)=frry(i)-(0.5d0*dU_driky)
               frrz(i)=frrz(i)-(0.5d0*dU_drikz)

c update stress
c               stsrxx=stsrxx-dxsav(i,k)*(0.5d0*dU_drikx)
c              stsrxy=stsrxy-dxsav(i,k)*(0.5d0*dU_driky)
c              stsrxz=stsrxz-dxsav(i,k)*(0.5d0*dU_drikz)
c               stsryy=stsryy-dysav(i,k)*(0.5d0*dU_driky)
c               stsryz=stsryz-dysav(i,k)*(0.5d0*dU_drikz)
c               stsrzz=stsrzz-dzsav(i,k)*(0.5d0*dU_drikz)
 
  903         continue !!!!!!!!!!finish 2nd loop over triplets

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!dU_a/drij!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c     !       d(fc(i,j))/d(rij)

                If (drij.le.(bigR-bigD))then
                 dcutoffij_drijx=0.0d0
                 else if (drij.le.(bigR+bigD)) then
                 dcutoffij_drijx=0.5d0*pi*(drij-bigR)/bigD
c                 dcutoffij_drijx=-(dcos(dcutoffij_drijx))
c     x                          *pi*dxsav(i,j)/(4.0d0*drij*bigD)
                 else
                 dcutoffij_drijx=0.0d0
                 Endif


               If (drij.le.(bigR-bigD))then
                 dcutoffij_drijy=0.0d0
                 else if (drij.le.(bigR+bigD)) then
                 dcutoffij_drijy=0.5d0*pi*((drij)-bigR)/bigD
                 dcutoffij_drijy=-(dcos(dcutoffij_drijy))
     x                         *pi*dysav(i,j)/(4.0d0*drij*bigD)
                  else
                 dcutoffij_drijy=0.0d0
                 Endif

               If (drij.le.(bigR-bigD))then
                 dcutoffij_drijz=0.0d0
                 else if (drij.le.(bigR+bigD)) then
                 dcutoffij_drijz=0.5d0*pi*(drij-bigR)/bigD
                 dcutoffij_drijz=-(dcos(dcutoffij_drijz))
     x                          *pi*dzsav(i,j)/(4.0d0*drij*bigD)
                 else
                 dcutoffij_drijz=0.0d0
                 Endif

c     !      d(fR(i,j))/dr(ij)
c               drepulse_drijx=-1.0d0*bigA*tf_lamda1*(dxsav(i,j)/drij)*
c     x                  dexp(-1.0d0*tf_lamda1*drij)
c               drepulse_drijy=-1.0d0*bigA*tf_lamda1*(dysav(i,j)/drij)*
c     x                  dexp(-1.0d0*tf_lamda1*drij)
c               drepulse_drijz=-1.0d0*bigA*tf_lamda1*(dzsav(i,j)/drij)*
c     x                   dexp(-1.0d0*tf_lamda1*drij)

c  !       d(fA(i,j))/dr(ij)
c                dattract_drijx=bigB*tf_lamda2*(dxsav(i,j)/drij)
c     x                         *dexp(-1.0d0*tf_lamda2*drij)
c                dattract_drijy=bigB*tf_lamda2*(dysav(i,j)/drij)
c     x                         *dexp(-1.0d0*tf_lamda2*drij)
c                dattract_drijz=bigB*tf_lamda2*(dzsav(i,j)/drij)
c     x                         *dexp(-1.0d0*tf_lamda2*drij)


c differential
               dU_drijx_a=( dcutoffij_drijx*(repulse+(bij*attract)) )
     x         +cutoffij*(drepulse_drijx+(bij*dattract_drijx))
               dU_drijy_a=( dcutoffij_drijy*(repulse+(bij*attract)) )
     x         +cutoffij*(drepulse_drijy+(bij*dattract_drijy))
               dU_drijz_a=( dcutoffij_drijz*(repulse+(bij*attract)) )
     x         +cutoffij*(drepulse_drijz+(bij*dattract_drijz))


c update forces
               frrx(j)=frrx(j)+(0.5d0*dU_drijx_a)
               frry(j)=frry(j)+(0.5d0*dU_drijy_a)
               frrz(j)=frrz(j)+(0.5d0*dU_drijz_a)
               frrx(i)=frrx(i)-(0.5d0*dU_drijx_a)
               frry(i)=frry(i)-(0.5d0*dU_drijy_a)
               frrz(i)=frrz(i)-(0.5d0*dU_drijz_a)
c update stress
c               stsrxx=stsrxx-dxsav(i,j)*(0.5d0*dU_drijx_a)
c               stsrxy=stsrxy-dxsav(i,j)*(0.5d0*dU_drijy_a)
c               stsrxz=stsrxz-dxsav(i,j)*(0.5d0*dU_drijz_a)
c               stsryy=stsryy-dysav(i,j)*(0.5d0*dU_drijy_a)
c               stsryz=stsryz-dysav(i,j)*(0.5d0*dU_drijz_a)
c               stsrzz=stsrzz-dzsav(i,j)*(0.5d0*dU_drijz_a)






  910         continue
            
  900         continue
         write(*,*) "end tf"

         write(*,*) "engpetot", engpetot
         write(89,*)engTF2,engTF2_2,engTF2_3
  
               endif

cc Harmonic potential
      if ((harmpotlog))then
cc         call separations
         pi=4.0d0*atan(1.0d0)

         if (.not.mixlog) then

            stsrxx=0.0d0
            stsrxy=0.0d0
            stsrxz=0.0d0
            stsryy=0.0d0
            stsryz=0.0d0
            stsrzz=0.0d0

            engpetot=0.0d0

            do 995 i=1,num
               frrx(i)=0.0d0
               frry(i)=0.0d0
               frrz(i)=0.0d0
  995       continue

         end if

cc loop over pre-defined pairs only.
         do 996 k=1,npairs

            i=npairident1(k)
            j=npairident2(k)

            ipoint=ntype(i)
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

            drijsq=(dxcf**2.0d0)
     x          +(dycf**2.0d0)
     x          +(dzcf**2.0d0)
            drij=drijsq**(0.5d0)

            u=xkharm(ipoint,jpoint)
     x       *((drij-r0harm(ipoint,jpoint))**2.0d0)

            engsr=engsr+u

            dudr=2.0d0*xkharm(ipoint,jpoint)
     x       *(drij-r0harm(ipoint,jpoint))
            frrxi=dudr*dxcf/drij
            frryi=dudr*dycf/drij
            frrzi=dudr*dzcf/drij

            frrx(i)=frrx(i)-frrxi
            frry(i)=frry(i)-frryi
            frrz(i)=frrz(i)-frrzi
            frrx(j)=frrx(j)+frrxi
            frry(j)=frry(j)+frryi
            frrz(j)=frrz(j)+frrzi

            stsrxx=stsrxx+frrxi*dxcf
            stsrxy=stsrxy+frrxi*dycf
            stsrxz=stsrxz+frrxi*dzcf
            stsryy=stsryy+frryi*dycf
            stsryz=stsryz+frryi*dzcf
            stsrzz=stsrzz+frrzi*dzcf

  996 continue

      endif

cc shifted 1/r^24 potential with harmonic term?
      if (harmpotlog2) then
       engrep=0
       do 605 i=2,num

            ipoint=ntype(i)

            do 611 j=1,i-1

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

               if (dr.gt.cut(ipoint,jpoint)) goto 611

               drrec=1.0d0/dr
               drsqrec=drrec*drrec
               dr24rec=drrec**24.0d0
               dr12rec=drrec**12.0d0

               ulj=fourepsil(ipoint,jpoint)
     x                     *((sig24(ipoint,jpoint)*dr24rec)
     x                     -(sig12(ipoint,jpoint)*dr12rec))
               ulj=ulj+epsilon(ipoint,jpoint)
               gam=24.0d0*fourepsil(ipoint,jpoint)*sig24(ipoint,jpoint)
     x         *dr24rec*drsqrec
     x         -12.0d0*fourepsil(ipoint,jpoint)*sig12(ipoint,jpoint)
     x         *dr12rec*drsqrec

cc            write(666,*)i,j,dr,ulj
cc            write(666,*)sig24(ipoint,jpoint),dr24rec
cc     x                 ,sig12(ipoint,jpoint),dr12rec
cc     x                 ,fourepsil(ipoint,jpoint),dr12rec
               engsr=engsr+ulj
               engrep=engrep+ulj
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

  611       continue

  605    continue

      endif

cc linear potential with harmonic term.
      if (harmpotlog3) then
       engrep=0
       do 606 i=2,num

            ipoint=ntype(i)

            do 612 j=1,i-1

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

               if (dr.gt.cut(ipoint,jpoint)) goto 612

               ulj=(xklin(ipoint,jpoint)*dr)+u0lin(ipoint,jpoint)

               gam=-xklin(ipoint,jpoint)/dr

               engsr=engsr+ulj
               engrep=engrep+ulj
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

  612       continue

  606    continue

      endif

         write(*,*) "engsr,engpetot",engsr,engpetot
         write(*,*) "engrep",engrep
      engpetot=engpetot+engsr
         write(*,*) "engpetotfinal", engpetot
                  return

         end

