      subroutine rearrange

c************************************************************
c   Calculates the SR energy and forces at each site.
c
c   Possible models are the i) EPP
c                          ii) CIM for anion-cation + EPP for others
c             
c************************************************************

      implicit double precision (a-h,o-z)

      include 'common.inc'
c
c Zero energy accumulators.
c
      engsr=0.0d0

      do 10 i=1,num
         engft1(i)=0.0d0
         engft2(i)=0.0d0
         engft3(i)=0.0d0
   10 continue

c TMP - CIM
      if (cim1log) then

         do 50 i=1,nsp(1)

            ipoint=ntype(i)

            do 60 j=nsp(1)+1,num

               jpoint=ntype(j)
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
               dxsav(i,j)=hlab2(1,1)*dxcf
     x                +hlab2(1,2)*dycf+hlab2(1,3)*dzcf
               dysav(i,j)=hlab2(2,1)*dxcf
     x                +hlab2(2,2)*dycf+hlab2(2,3)*dzcf
               dzsav(i,j)=hlab2(3,1)*dxcf
     x                +hlab2(3,2)*dycf+hlab2(3,3)*dzcf

               dxsav(j,i)=-dxsav(i,j)
               dysav(j,i)=-dysav(i,j)
               dzsav(j,i)=-dzsav(i,j)
               drsq=dxsav(i,j)*dxsav(i,j)
     x          +dysav(i,j)*dysav(i,j)
     x          +dzsav(i,j)*dzsav(i,j)
c
c implement cut-off:
cc            if (drsq.ge.rsqmax) goto 40

               drsqrec=1.0d0/drsq
               dr=dsqrt(drsq)
               drrec=1.0d0/dr

               expft=(ftb(ipoint,jpoint)/dr)
     x           *exp(-ftalp(ipoint,jpoint)*
     x           (dr-delta(i)))
               expft2=(ftb2(ipoint,jpoint)/dr)
     x           *exp(-ftbeta(ipoint,jpoint)*
     x           (dr-delta(i)))
               expft3=(ftb3(ipoint,jpoint)/dr)
     x           *exp(-ftgamma(ipoint,jpoint)*
     x           (dr-delta(i)))

               engft1(i)=engft1(i)+expft
               engft2(i)=engft2(i)+expft2
               engft3(i)=engft3(i)+expft3

   60       continue

   50    continue

      endif

      if (cim2log) then

         do 55 i=1,nsp(1)

            ipoint=ntype(i)

            do 65 j=nsp(1)+1,num

               jpoint=ntype(j)
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
               dxsav(i,j)=dxcf
               dysav(i,j)=dycf
               dzsav(i,j)=dzcf

               dxsav(j,i)=-dxsav(i,j)
               dysav(j,i)=-dysav(i,j)
               dzsav(j,i)=-dzsav(i,j)
               drsq=dxsav(i,j)*dxsav(i,j)
     x          +dysav(i,j)*dysav(i,j)
     x          +dzsav(i,j)*dzsav(i,j)
c
c implement cut-off:
cc            if (drsq.ge.rsqmax) goto 40

               drsqrec=1.0d0/drsq
               dr=dsqrt(drsq)
               drrec=1.0d0/dr

               expft=(ftb(ipoint,jpoint))
     x           *exp(-ftalp(ipoint,jpoint)*
     x           (dr-delta(i)))
               expft2=(ftb2(ipoint,jpoint))
     x           *exp(-ftbeta(ipoint,jpoint)*
     x           (dr-delta(i)))
               expft3=(ftb3(ipoint,jpoint))
     x           *exp(-ftgamma(ipoint,jpoint)*
     x           (dr-delta(i)))

               engft1(i)=engft1(i)+expft
               engft2(i)=engft2(i)+expft2
               engft3(i)=engft3(i)+expft3

   65       continue

   55    continue

      endif

      return

      end
