      subroutine rescale

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

            vxtot=vxtot+vx(n1)
            vytot=vytot+vy(n1)
            vztot=vztot+vz(n1)

            n1=n1+1

   20    continue

         sprec=1.0d0/dble(nsp(i))
         vxmean=vxtot*sprec
         vymean=vytot*sprec
         vzmean=vztot*sprec

         do 30 j=1,nsp(i)

            vx(n2)=vx(n2)-vxmean
            vy(n2)=vy(n2)-vymean
            vz(n2)=vz(n2)-vzmean

c
c Noncubic box code:
c
            vttx=hlab2(1,1)*vx(n2)+hlab2(1,2)*vy(n2)+hlab2(1,3)*vz(n2)
            vtty=hlab2(2,1)*vx(n2)+hlab2(2,2)*vy(n2)+hlab2(2,3)*vz(n2)
            vttz=hlab2(3,1)*vx(n2)+hlab2(3,2)*vy(n2)+hlab2(3,3)*vz(n2)

            xke2=xke2+(vttx*vttx+vtty*vtty+vttz*vttz)*hmass(i)


            n2=n2+1

   30    continue

         else
c
c Cope with the case that we only have one
c ion of a given species. (The mean velocity is
c then the actual velocity!)
c
            vttx=hlab2(1,1)*vx(n1)+hlab2(1,2)*vy(n1)+hlab2(1,3)*vz(n1)
            vtty=hlab2(2,1)*vx(n1)+hlab2(2,2)*vy(n1)+hlab2(2,3)*vz(n1)
            vttz=hlab2(3,1)*vx(n1)+hlab2(3,2)*vy(n1)+hlab2(3,3)*vz(n1)

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

            vx(n)=vx(n)*tempfac
            vy(n)=vy(n)*tempfac
            vz(n)=vz(n)*tempfac

            tini=tini+(vx(n)*vx(n)+vy(n)*vy(n)+
     x                 vz(n)*vz(n))*hmass(i)
            n=n+1

   80    continue
   70 continue

      tkinres=tini*gtrankin

      write(6,*)
      write(6,*)'**** Velocity rescaling complete ****'
      write(6,*)

      return

      end
