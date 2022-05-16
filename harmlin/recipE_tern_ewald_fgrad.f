      subroutine recipE_dippim

c************************************************************
c Calculates the reciprocal space summation part of the
c Ewald sum. (Dipolar version -INCLUDING THE DIPOLE-DIPOLE FORCES)
c************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'
c
c Local double precision arrays:
      double precision elcall(nummax,0:nktot),emcall(nummax,0:nktot),
     x                 encall(nummax,0:nktot),elsall(nummax,0:nktot),
     x                 emsall(nummax,0:nktot),ensall(nummax,0:nktot)
      double precision clmall(nummax),slmall(nummax),
     x                 ckcnoqall(nummax),cksnoqall(nummax)
c
c Zero all force and electric field arrays:
c
      dipsum=0.0d0

      do 10 i=1,num

         frrx(i)=0.0d0
         frry(i)=0.0d0
         frrz(i)=0.0d0

         elecx(i)=0.0d0
         elecy(i)=0.0d0
         elecz(i)=0.0d0

         exx(i)=0.0d0
         eyy(i)=0.0d0
         ezz(i)=0.0d0
         exy(i)=0.0d0
         exz(i)=0.0d0
         eyz(i)=0.0d0

         dipsum=dipsum+(xmu(i)*xmu(i))
     x                +(ymu(i)*ymu(i))
     x                +(zmu(i)*zmu(i))

   10 continue

c
c zero structure factor arrays:
c
      do 20 i=1,ksqmax

         do 30 j=1,num

            jpoint=ntype(j)
            arsk(jpoint,i)=0.0d0
            aisk(jpoint,i)=0.0d0

   30    continue

   20 continue
c
c zero energy accumulators:
c
      engpetot=0.0d0
      qmueng=0.0d0
      xmumueng=0.0d0
c
c zero the coulomb and polarization parts of stress tensor:
c the other bits are initialized in realE.
c
      stp2xx=0.0d0
      stp2xy=0.0d0
      stp2xz=0.0d0
      stp2yy=0.0d0
      stp2yz=0.0d0
      stp2zz=0.0d0

      stpxx=0.0d0
      stpxy=0.0d0
      stpxz=0.0d0
      stpyy=0.0d0
      stpyz=0.0d0
      stpzz=0.0d0

      stcxx=0.0d0
      stcxy=0.0d0
      stcxz=0.0d0
      stcyy=0.0d0
      stcyz=0.0d0
      stczz=0.0d0

cc      goto 222
c
c Setup cos and sin arrays for addition rule.
c
      do 50 i=1,num
       
          elcall(i,0)=1.0d0
          emcall(i,0)=1.0d0
          encall(i,0)=1.0d0
          elsall(i,0)=0.0d0
          emsall(i,0)=0.0d0
          ensall(i,0)=0.0d0

          elcall(i,1)=cos(twopiboxx*x(i))
          emcall(i,1)=cos(twopiboxy*y(i))
          encall(i,1)=cos(twopiboxz*z(i))
          elsall(i,1)=sin(twopiboxx*x(i))
          emsall(i,1)=sin(twopiboxy*y(i))
          ensall(i,1)=sin(twopiboxz*z(i))

   50 continue
c
c Calculate all cosines and sines
c
      do 60 l=2,kmax

         do 90 i=1,num

            elcall(i,l)=elcall(i,l-1)*elcall(i,1)
     x                -elsall(i,l-1)*elsall(i,1)
            emcall(i,l)=emcall(i,l-1)*emcall(i,1)
     x                -emsall(i,l-1)*emsall(i,1)
            encall(i,l)=encall(i,l-1)*encall(i,1)
     x                -ensall(i,l-1)*ensall(i,1)
            elsall(i,l)=elsall(i,l-1)*elcall(i,1)
     x                +elcall(i,l-1)*elsall(i,1)
            emsall(i,l)=emsall(i,l-1)*emcall(i,1)
     x                +emcall(i,l-1)*emsall(i,1)
            ensall(i,l)=ensall(i,l-1)*encall(i,1)
     x                +encall(i,l-1)*ensall(i,1)
   90    continue
   60 continue

      mmin=0
      nmin=1
c
c Primary vector loop.
c
      do 110 ll=0,kmax

         l=iabs(ll)
         rl=dble(ll)*twopi
c
c Generate the k-vectors corresponding to the noncubic box.
c These are the cartesian components of k, which we require since
c we are taking cartesian derivatives. fullhi derives from the
c inverse of the matrix containing the full cell vectors. We are
c actually using the transpose of fullhi.
c
c (fullhi)^T =  ( fullhi(1,1)  fullhi(2,1)  fullhi(3,1) )
c               ( fullhi(1,2)  fullhi(2,2)  fullhi(3,2) )
c               ( fullhi(1,3)  fullhi(2,3)  fullhi(3,3) )
c
c The full rkx3 (which is k_x) is 2pi*n_1*fullhi(1,1) + 2pi*n_2*fullhi(1,2) +
c                                 2pi*n_3*fullhi(1,3).
c
         rkx1=rl*fullhi(1,1)
         rky1=rl*fullhi(1,2)
         rkz1=rl*fullhi(1,3)
c
c Secondary vector loop.
c
         do 120 mm=mmin,kmax

            m=iabs(mm)
            rm=dble(mm)*twopi
c
c Generate the k-vectors corresponding to the noncubic box.
c
            rkx2=rkx1+rm*fullhi(2,1)
            rky2=rky1+rm*fullhi(2,2)
            rkz2=rkz1+rm*fullhi(2,3)

            if (mm.ge.0) then
c
c Utilise addition rule to generate new values.
c
               do 160 i=1,num

                  clmall(i)=elcall(i,l)*emcall(i,m)
     x                    -elsall(i,l)*emsall(i,m)                     
                  slmall(i)=elsall(i,l)*emcall(i,m)
     x                    +emsall(i,m)*elcall(i,l)

  160          continue
            else

               do 190 i=1,num

                  clmall(i)=elcall(i,l)*emcall(i,m)
     x                    +elsall(i,l)*emsall(i,m)                     
                  slmall(i)=elsall(i,l)*emcall(i,m)
     x                    -emsall(i,m)*elcall(i,l)

  190          continue

            endif
c
c Tertiary vector loop.
c
            do 200 nn=nmin,kmax

               n=iabs(nn)
               rn=dble(nn)*twopi
c
c Generate the k-vectors corresponding to the noncubic box.
c
               rkx3=rkx2+rn*fullhi(3,1)
               rky3=rky2+rn*fullhi(3,2)
               rkz3=rkz2+rn*fullhi(3,3)
c
c xkk is k^2, generated from the cartesian components of k.
c
               xkk=rkx3*rkx3+rky3*rky3+rkz3*rkz3

               kk=ll*ll+mm*mm+nn*nn
               if (kk.eq.0) goto 200
               if (xkk.le.rksqmax) then
c
c Generate total cosines and sines.
c
                  if (nn.ge.0) then

                     do 230 i=1,num

                        ckcnoqall(i)=clmall(i)*encall(i,n)
     x                             -slmall(i)*ensall(i,n)
                        cksnoqall(i)=slmall(i)*encall(i,n)
     x                             +clmall(i)*ensall(i,n)

  230                continue

                  else

                     do 260 i=1,num

                        ckcnoqall(i)=clmall(i)*encall(i,n)
     x                             +slmall(i)*ensall(i,n)
                        cksnoqall(i)=slmall(i)*encall(i,n)
     x                             -clmall(i)*ensall(i,n)       

  260                continue

                  endif

                  ckcs=0.0d0
                  ckss=0.0d0
                  ckcsaf=0.0d0
                  ckssaf=0.0d0
                  cx=0.0d0
                  cy=0.0d0
                  cz=0.0d0
                  sx=0.0d0
                  sy=0.0d0
                  sz=0.0d0

c	changes for ternary systems

                  do 290 i=1,num

                   ipoint=ntype(i)

                     ckcsaf=ckcsaf+ckcnoqall(i)
                     ckssaf=ckssaf+cksnoqall(i)

                     ckcs=ckcs+ckcnoqall(i)*q(i)
                     ckss=ckss+cksnoqall(i)*q(i)

                     cx=cx+ckcnoqall(i)*xmu(i)
                     cy=cy+ckcnoqall(i)*ymu(i)
                     cz=cz+ckcnoqall(i)*zmu(i)


                     sx=sx+cksnoqall(i)*xmu(i)
                     sy=sy+cksnoqall(i)*ymu(i)
                     sz=sz+cksnoqall(i)*zmu(i)

                     arsk(ipoint,kk)=arsk(ipoint,kk)+ckcnoqall(i)
                     aisk(ipoint,kk)=aisk(ipoint,kk)+cksnoqall(i)

  290             continue


c
c structure factor calculation: (Needs modifying.)
c

                 do 291 i=1,nspec
                  do 292 j=i,nspec
        sk(i,j,kk)=sk(i,j,kk)+(arsk(i,kk)*arsk(j,kk))
     x                              +(aisk(i,kk)*aisk(j,kk))
292            continue
291            continue
               do 293 i=1,nspec
                  arsk(i,kk)=0.0d0
                  aisk(i,kk)=0.0d0
293            continue


c                  skaa(kk)=skaa(kk)+(arsk(1,kk)*arsk(1,kk))
c     x                              +(aisk(1,kk)*aisk(1,kk))
c                  skac(kk)=skac(kk)+(arsk(1,kk)*arsk(2,kk))
c     x                              +(aisk(1,kk)*aisk(2,kk))
c                  skcc(kk)=skcc(kk)+(arsk(2,kk)*arsk(2,kk))
c     x                              +(aisk(2,kk)*aisk(2,kk))
c                  arsk(1,kk)=0.0d0
c                  arsk(2,kk)=0.0d0
c                  aisk(1,kk)=0.0d0
c                  aisk(2,kk)=0.0d0
c
c Setup akk. Note this is different to the cubic code which
c defined an ak() array externally (calculated in setup).
c The variable akk is equivalent to this but now calculated
c in the main recip. code in order to clarify the role
c of the L^2 term that appears in it.
c
                  akk=exp(etaconst*xkk)/xkk
c
c charge-charge energy:
c
                  efac=akk*(ckss*ckss+ckcs*ckcs)
                  engpetot=engpetot+efac
c
c N.B. a self-interaction term needs to be subtracted at the end:
c this is chgcorrec, calculated in setup.f
c Use the cartesian components of k:
c
                  arl=akk*rkx3
                  arm=akk*rky3
                  arn=akk*rkz3
                  arll=arl*rkx3
                  armm=arm*rky3
                  arnn=arn*rkz3
                  arlm=arl*rky3
                  arln=arl*rkz3
                  armn=arm*rkz3
c
c k-space stress tensor calculation:
c
cc (rk2 is now xkk)  rk2=rl*rl+rm*rm+rn*rn
cc                  stfac=(1.0d0+xkk/(4.0d0*etasq))/xkk

                  stfac=(4.0d0*etasq+xkk)/(4.0d0*etasq*xkk)

                  clx=cx*rkx3
                  cly=cy*rky3
                  clz=cz*rkz3
                  slx=sx*rkx3
                  sly=sy*rky3
                  slz=sz*rkz3

                  rdotmuc=clx+cly+clz
                  rdotmus=slx+sly+slz
                  ddfac=rdotmuc*rdotmuc+rdotmus*rdotmus
c
c charge-charge contributions:
c
                  stcxx=stcxx+efac*(1.0d0-2.0d0*rkx3*rkx3*stfac)
                  stcyy=stcyy+efac*(1.0d0-2.0d0*rky3*rky3*stfac)
                  stczz=stczz+efac*(1.0d0-2.0d0*rkz3*rkz3*stfac)
                  stcxy=stcxy-efac*2.0d0*rkx3*rky3*stfac
                  stcxz=stcxz-efac*2.0d0*rkx3*rkz3*stfac
                  stcyz=stcyz-efac*2.0d0*rky3*rkz3*stfac
c
c dipole-dipole contributions:
c
                  stp2xx=stp2xx+akk*(ddfac*(1.0d0-2.0d0*rkx3*rkx3*stfac)
     x                  +2.0d0*(clx*rdotmuc+slx*rdotmus))
                  stp2yy=stp2yy+akk*(ddfac*(1.0d0-2.0d0*rky3*rky3*stfac)
     x                  +2.0d0*(cly*rdotmuc+sly*rdotmus))
                  stp2zz=stp2zz+akk*(ddfac*(1.0d0-2.0d0*rkz3*rkz3*stfac)
     x                  +2.0d0*(clz*rdotmuc+slz*rdotmus))

                  stp2xy=stp2xy+akk*(-ddfac*2.0d0*rkx3*rky3*stfac
     x                  +(rkx3*cy+rky3*cx)*rdotmuc
     x                  +(rkx3*sy+rky3*sx)*rdotmus)
                  stp2xz=stp2xz+akk*(-ddfac*2.0d0*rkx3*rkz3*stfac
     x                  +(rkx3*cz+rkz3*cx)*rdotmuc
     x                  +(rkx3*sz+rkz3*sx)*rdotmus)
                  stp2yz=stp2yz+akk*(-ddfac*2.0d0*rky3*rkz3*stfac
     x                  +(rky3*cz+rkz3*cy)*rdotmuc
     x                  +(rky3*sz+rkz3*sy)*rdotmus)
c
c charge-dipole contributions:
c
                  dqfac=-2.0d0*akk*(ckcs*rdotmus-ckss*rdotmuc)

                  stpxx=stpxx+2.0d0*arl*(cx*ckss-sx*ckcs)
     x                 +dqfac*(1.0d0-2.0d0*rkx3*rkx3*stfac)
                  stpyy=stpyy+2.0d0*arm*(cy*ckss-sy*ckcs)
     x                 +dqfac*(1.0d0-2.0d0*rky3*rky3*stfac)
                  stpzz=stpzz+2.0d0*arn*(cz*ckss-sz*ckcs)
     x                 +dqfac*(1.0d0-2.0d0*rkz3*rkz3*stfac)

                  stpxy=stpxy+arl*(cy*ckss-sy*ckcs)
     x                 +arm*(cx*ckss-sx*ckcs)
     x                 -dqfac*2.0d0*rkx3*rky3*stfac
                  stpxz=stpxz+arl*(cz*ckss-sz*ckcs)
     x                 +arn*(cx*ckss-sx*ckcs)
     x                 -dqfac*2.0d0*rkx3*rkz3*stfac
                  stpyz=stpyz+arm*(cz*ckss-sz*ckcs)
     x                 +arn*(cy*ckss-sy*ckcs)
     x                 -dqfac*2.0d0*rky3*rkz3*stfac

                  do 320 i=1,num
c
c force and electric field due to charges:
c
                     elefld=cksnoqall(i)*ckcs
     x                     -ckcnoqall(i)*ckss
                     qforce=elefld*q(i)

                     frrx(i)=frrx(i)+arl*qforce
                     frry(i)=frry(i)+arm*qforce
                     frrz(i)=frrz(i)+arn*qforce

                     elecx(i)=elecx(i)+arl*elefld
                     elecy(i)=elecy(i)+arm*elefld
                     elecz(i)=elecz(i)+arn*elefld
c
c -mu.field gives charge-dipole energy:
c extra factor of 2 at end over and above q-q and mu-mu energies
c since there are two of these terms in the energy equation:
c
                     qmueng=qmueng-arl*elefld*xmu(i)
     x                 -ymu(i)*arm*elefld-zmu(i)*arn*elefld
c
c this is the qTmu term:
c force and electric field due to dipoles:
c
                     elefld=-rkx3*((ckcnoqall(i)*cx
     x                      +cksnoqall(i)*sx))
     x                      -rky3*((ckcnoqall(i)*cy
     x                      +cksnoqall(i)*sy))
     x                      -rkz3*((ckcnoqall(i)*cz
     x                      +cksnoqall(i)*sz))

                     elecx(i)=elecx(i)+arl*elefld
                     elecy(i)=elecy(i)+arm*elefld
                     elecz(i)=elecz(i)+arn*elefld

                     qforce=elefld*q(i)

                     frrx(i)=frrx(i)+arl*qforce
                     frry(i)=frry(i)+arm*qforce
                     frrz(i)=frrz(i)+arn*qforce
c
c -mu.field gives dipole-dipole energy:
c multiplied by same factor as q-q energy at end.
c
                     xmumueng=xmumueng-xmu(i)*arl*elefld
     x                   -ymu(i)*arm*elefld-zmu(i)*arn*elefld
c
c N.B. there is a self-interaction term here as well, REMOVED
C AT THE BOTTOM
c
c muTq term:
c this only interacts with the ion dipoles to create a force
c on the ions, the electric field gradient is irrelevant to the dipole
c equations of motion.
c
                     rkdotmui= rkx3*xmu(i)+rky3*ymu(i)+rkz3*zmu(i)
                     elefld= rkdotmui
     x                        *(ckcnoqall(i)*ckcs
     x                        +cksnoqall(i)*ckss)

                     frrx(i)=frrx(i)+arl*elefld
                     frry(i)=frry(i)+arm*elefld
                     frrz(i)=frrz(i)+arn*elefld
c
c energy part is as above, so multiply that term by two.

C..........dipole-dipole forces
                     prod= rkdotmui*
     x                       (cksnoqall(i)*rdotmuc
     x                       -ckcnoqall(i)*rdotmus)

                     frrx(i)=frrx(i)+arl*prod
                     frry(i)=frry(i)+arm*prod
                     frrz(i)=frrz(i)+arn*prod

c added FG.
                     exx(i)=exx(i)+arll*egrad
                     eyy(i)=eyy(i)+armm*egrad
                     ezz(i)=ezz(i)+arnn*egrad
                     exy(i)=exy(i)+arlm*egrad
                     exz(i)=exz(i)+arln*egrad
                     eyz(i)=eyz(i)+armn*egrad

  320             continue

               endif

  200       continue
            nmin=-kmax

  120    continue
         mmin=-kmax

  110 continue
c
c chgcorrec term (charge-charge self-interaction)is calculated in setup.
c
      engpetot=fourpicell*engpetot-chgcorrec
      qmueng=qmueng*2.0d0*fourpicell
      xmumueng=xmumueng*fourpicell

c...dipole-dipole self-interaction
      xmumuengfac=2.0d0*(eta**3.0d0)/(3.0d0*sqrpi)*dipsum
      xmumueng=xmumueng-xmumuengfac

      stcxx=stcxx*fourpicell
      stcxy=stcxy*fourpicell
      stcxz=stcxz*fourpicell
      stcyy=stcyy*fourpicell
      stcyz=stcyz*fourpicell
      stczz=stczz*fourpicell

      stpxx=stpxx*fourpicell
      stpxy=stpxy*fourpicell
      stpxz=stpxz*fourpicell
      stpyy=stpyy*fourpicell
      stpyz=stpyz*fourpicell
      stpzz=stpzz*fourpicell

      stp2xx=stp2xx*fourpicell
      stp2xy=stp2xy*fourpicell
      stp2xz=stp2xz*fourpicell
      stp2yy=stp2yy*fourpicell
      stp2yz=stp2yz*fourpicell
      stp2zz=stp2zz*fourpicell

c         write (703,*) (stcxx+stcyy+stczz)*cellvol3rec,
c     x                 (stpxx+stpyy+stpzz)*cellvol3rec,
c     x                 (stp2xx+stp2yy+stp2zz)*cellvol3rec

c         write (704,*) engpetot*cellvol3rec,2.0d0*qmueng*
c     x                 cellvol3rec,3.0d0*xmumueng*cellvol3rec
c
c multiply forces and fields by 4pi/vol (double A&T factor due to half
c counting in k-loop)
c
c........fac is 2 times fourpicell (two forces from differentiating double sum)
      do 350 i=1,num

         frrx(i)=frrx(i)*fac
         frry(i)=frry(i)*fac
         frrz(i)=frrz(i)*fac

         elecx(i)=elecx(i)*fac
         elecy(i)=elecy(i)*fac
         elecz(i)=elecz(i)*fac

c added FG.
         exx(i)=exx(i)*fac
         eyy(i)=eyy(i)*fac
         ezz(i)=ezz(i)*fac
         exy(i)=exy(i)*fac
         exz(i)=exz(i)*fac
         eyz(i)=eyz(i)*fac

c............correct the electric field for the mu-mu self-interaction
         fac2=4.0d0*(eta**3.0d0)/(3.0d0*sqrpi)

         elecx(i)=elecx(i)+(fac2*xmu(i))
         elecy(i)=elecy(i)+(fac2*ymu(i))
         elecz(i)=elecz(i)+(fac2*zmu(i))

  350 continue

cc  222 continue

      return
      end
