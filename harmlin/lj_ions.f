      subroutine lj_ions

c************************************************************
c
c  Calculate the energy and forces due to the lennard-jones
c fluid present interacting with the ions. 
c  The potential used is a straight 12-6.
c
c************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'
      include 'lj.inc'

      engpetotlj=0.0d0

      engljionsr=0.0d0
      engljion=0.0d0

c=====================================================================
c  Loop over all the Lennard-Jones particles in the system.
      do 1000 i=1,numlj

         do 1010 j=1,num

            ipoint=ntype(j)

            dxcf=xlj(i)-x(j)
            dycf=ylj(i)-y(j)
            dzcf=zlj(i)-z(j)

            dxcf=dxcf-boxlenx*int(dxcf*halfboxxrec)
            dycf=dycf-boxleny*int(dycf*halfboxyrec)
            dzcf=dzcf-boxlenz*int(dzcf*halfboxzrec)

            dx=hlab2(1,1)*dxcf
     x                +hlab2(1,2)*dycf+hlab2(1,3)*dzcf
            dy=hlab2(2,1)*dxcf
     x                +hlab2(2,2)*dycf+hlab2(2,3)*dzcf
            dz=hlab2(3,1)*dxcf
     x                +hlab2(3,2)*dycf+hlab2(3,3)*dzcf

            drsq=dx*dx+dy*dy+dz*dz

            dr=dsqrt(drsq)

            drrec=1.0d0/dr
            drsqrec=drrec*drrec
            dr6rec=drrec**6.0d0
            dr12rec=drrec**12.0d0

            ulj=xlambdalj
     x         *fourepsilcyl(ipoint)*((sig12(ipoint)*dr12rec)
     x                     -(sig6(ipoint)*dr6rec))

            gam=12.0d0*fourepsilcyl(ipoint)*sig12(ipoint)
     x         *dr12rec*drsqrec
     x         -6.0d0*fourepsilcyl(ipoint)*sig6(ipoint)
     x         *dr6rec*drsqrec

            gam=gam*xlambdalj

            engljion=engljion+ulj

            frrx(j)=frrx(j)-dx*gam
            frry(j)=frry(j)-dy*gam
            frrz(j)=frrz(j)-dz*gam

            stsrxx=stsrxx+dx*gam*dx
            stsrxy=stsrxy+dx*gam*dy
            stsrxz=stsrxz+dx*gam*dz
            stsryy=stsryy+dy*gam*dy
            stsryz=stsryz+dy*gam*dz
            stsrzz=stsrzz+dz*gam*dz

 1010    continue

 1000 continue

      engpetotlj=engpetotlj+engljion

      write(88,*)engpetotlj

      return

      end
