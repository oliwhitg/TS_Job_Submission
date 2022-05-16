      subroutine rescalen

c************************************************************
c  Performs velocity rescaling to the desired temperature.
c************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'

      n1=1
      n2=1
      xke2=0.0d0

      do 10 i=1,nspec

         vxtot=0.0d0
         vytot=0.0d0
         vztot=0.0d0

         if (nsp(i).ne.1) then

         do 20 j=1,nsp(i)

            vxtot=vxtot+vxn(n1)
            vytot=vytot+vyn(n1)
            vztot=vztot+vzn(n1)

            n1=n1+1

   20    continue

         sprec=1.0d0/dble(nsp(i))
         vxmean=vxtot*sprec
         vymean=vytot*sprec
         vzmean=vztot*sprec

         do 30 j=1,nsp(i)

            vxn(n2)=vxn(n2)-vxmean
            vyn(n2)=vyn(n2)-vymean
            vzn(n2)=vzn(n2)-vzmean

c
c Noncubic box code:
c
            vttx=hlab2(1,1)*vxn(n2)+hlab2(1,2)*vyn(n2)
     x          +hlab2(1,3)*vzn(n2)
            vtty=hlab2(2,1)*vxn(n2)+hlab2(2,2)*vyn(n2)
     x          +hlab2(2,3)*vzn(n2)
            vttz=hlab2(3,1)*vxn(n2)+hlab2(3,2)*vyn(n2)
     x          +hlab2(3,3)*vzn(n2)

            xke2=xke2+(vttx*vttx+vtty*vtty+vttz*vttz)*hmass(i)

            n2=n2+1

   30    continue

         else
c
c Cope with the case that we only have one
c ion of a given species. (The mean velocity is
c then the actual velocity!)
c
            vttx=hlab2(1,1)*vxn(n1)+hlab2(1,2)*vyn(n1)
     x          +hlab2(1,3)*vzn(n1)
            vtty=hlab2(2,1)*vxn(n1)+hlab2(2,2)*vyn(n1)
     x          +hlab2(2,3)*vzn(n1)
            vttz=hlab2(3,1)*vxn(n1)+hlab2(3,2)*vyn(n1)
     x          +hlab2(3,3)*vzn(n1)

            xke2=xke2+(vttx*vttx+vtty*vtty+vttz*vttz)*hmass(i)



c
c need to update both counters.
c
           n1=n1+1
           n2=n2+1

         endif

   10 continue
c
c Calculate the kinetic temperature and the rescaling factor.
c
      tkin=xke2*gtrankin
      tempfac=dsqrt(trantemp/tkin)

      n=1
      tini=0.0d0

      do 70 i=1,nspec
         do 80 j=1,nsp(i)

            vxn(n)=vxn(n)*tempfac
            vyn(n)=vyn(n)*tempfac
            vzn(n)=vzn(n)*tempfac

            tini=tini+(vxn(n)*vxn(n)+vyn(n)*vyn(n)+
     x                 vzn(n)*vzn(n))*hmass(i)
            n=n+1

   80    continue
   70 continue

      tkinres=tini*gtrankin

      write(6,*)
      write(6,*)'**** Velocity rescaling complete ****'
      write(6,*)

      return

      end
