      subroutine rgdrecipE

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
      do 10 i=1,num

         frrx(i)=0.0d0
         frry(i)=0.0d0
         frrz(i)=0.0d0

         elecx(i)=0.0d0
         elecy(i)=0.0d0
         elecz(i)=0.0d0

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
c
c zero the coulomb and polarization parts of stress tensor:
c the other bits are initialized in realE.
c
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

      do 60 l=2,kmaxx

         do 90 i=1,num

            elcall(i,l)=elcall(i,l-1)*elcall(i,1)
     x                -elsall(i,l-1)*elsall(i,1)
            elsall(i,l)=elsall(i,l-1)*elcall(i,1)
     x                +elcall(i,l-1)*elsall(i,1)
   90    continue
   60 continue

      do 61 l=2,kmaxy

         do 91 i=1,num

            emcall(i,l)=emcall(i,l-1)*emcall(i,1)
     x                -emsall(i,l-1)*emsall(i,1)
            emsall(i,l)=emsall(i,l-1)*emcall(i,1)
     x                +emcall(i,l-1)*emsall(i,1)
   91    continue
   61 continue

      do 62 l=2,kmaxz

         do 92 i=1,num

            encall(i,l)=encall(i,l-1)*encall(i,1)
     x                -ensall(i,l-1)*ensall(i,1)
            ensall(i,l)=ensall(i,l-1)*encall(i,1)
     x                +encall(i,l-1)*ensall(i,1)
   92    continue
   62 continue

      mmin=0
      nmin=1
c
c Primary vector loop.
c
      do 110 ll=0,kmaxx

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
         do 120 mm=mmin,kmaxy

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
            do 200 nn=nmin,kmaxz

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

c	changes for ternary systems

                  do 290 i=1,num

                   ipoint=ntype(i)

                     ckcsaf=ckcsaf+ckcnoqall(i)
                     ckssaf=ckssaf+cksnoqall(i)

                     ckcs=ckcs+ckcnoqall(i)*q(i)
                     ckss=ckss+cksnoqall(i)*q(i)

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
c
c k-space stress tensor calculation:
c
                  stfac=(4.0d0*etasq+xkk)/(4.0d0*etasq*xkk)
c
c charge-charge contributions:
c
                  stcxx=stcxx+efac*(1.0d0-2.0d0*rkx3*rkx3*stfac)
                  stcyy=stcyy+efac*(1.0d0-2.0d0*rky3*rky3*stfac)
                  stczz=stczz+efac*(1.0d0-2.0d0*rkz3*rkz3*stfac)
                  stcxy=stcxy-efac*2.0d0*rkx3*rky3*stfac
                  stcxz=stcxz-efac*2.0d0*rkx3*rkz3*stfac
                  stcyz=stcyz-efac*2.0d0*rky3*rkz3*stfac

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

  320             continue

               endif

  200       continue
            nmin=-kmaxz

  120    continue
         mmin=-kmaxy

  110 continue
c
c chgcorrec term (charge-charge self-interaction)is calculated in setup.
c
      engpetot=fourpicell*engpetot-chgcorrec

      stcxx=stcxx*fourpicell
      stcxy=stcxy*fourpicell
      stcxz=stcxz*fourpicell
      stcyy=stcyy*fourpicell
      stcyz=stcyz*fourpicell
      stczz=stczz*fourpicell

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

  350 continue

cc  222 continue

      return
      end
