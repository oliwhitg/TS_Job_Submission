      subroutine debye_scherer

c************************************************************
c Calculates the  structure factor
c  .. input is the ion positions in the cell frame and
c... the cell matrix h. Output should be the structure factor
c .. i.e. this should handle the definition of the appropriate
c...cartesian length of the wavevector. We will also need some 
c...form factors.
c************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'
c
c Local double precision arrays:
      double precision elcall(nummax,0:nktot_ds)
     x                ,emcall(nummax,0:nktot_ds)
     x                ,encall(nummax,0:nktot_ds)
     x                ,elsall(nummax,0:nktot_ds)
     x                ,emsall(nummax,0:nktot_ds)
     x                ,ensall(nummax,0:nktot_ds)

c..............nktot_ds should be bigger than kmax_ds
      double precision clmall(nummax),slmall(nummax),
     x                 ckcnoqall(nummax),cksnoqall(nummax),
     x                 arsk_ds(nspmax),aisk_ds(nspmax)

c.....sk and norm need zeroing in setup, dk2_ds, kmax_ds, rksqmax_ds
c.....and nbin need initializing
c...........................................

 
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
      do 60 l=2,kmax_ds

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
      do 110 ll=0,kmax_ds

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
         do 120 mm=mmin,kmax_ds

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
            do 200 nn=nmin,kmax_ds

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
               if (xkk.le.rksqmax_ds) then
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

  

c	changes for ternary systems

                  do ipoint=1,nspec
                     arsk_ds(ipoint)=0.0d0
                     aisk_ds(ipoint)=0.0d0
                  enddo

                  do 290 i=1,num

                   ipoint=ntype(i)
 
                     arsk_ds(ipoint)=arsk_ds(ipoint)+ckcnoqall(i)
                     aisk_ds(ipoint)=aisk_ds(ipoint)+cksnoqall(i)

  290             continue

c
c structure factor calculation: (Needs modifying.)
c

                 kbin = int((xkk/rksqmax_ds)*float(nbin))
                 norm_ds(kbin)=norm_ds(kbin)+1
               
                 lk=1
               do 291 i=1,nspec
                 do 292 j=i,nspec
                    sk_ds(lk,kbin)= sk_ds(lk,kbin) 
     x                             +(arsk_ds(i)*arsk_ds(j))
     x                              +(aisk_ds(i)*aisk_ds(j))
                    lk=lk+1
292              continue
291            continue

              endif

200           continue
              nmin=-kmax_ds
120        continue
           mmin=-kmax_ds
110     continue
  
      return

      end
