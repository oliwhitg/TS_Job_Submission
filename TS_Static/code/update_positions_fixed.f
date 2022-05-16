      subroutine update_positions

c************************************************************
c   Updates:  i) ion positions.
c            ii) ion velocities.
c           iii) dipole components.
c            iv) dipole component velocities.
c
c   Note that velocities are at the half timesteps.
c************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'

      if (moveions) then
c
c Update ion positions and velocities.
c
      do 10 i=1,num

         if (fixedarray(i)) then
            goto 10
         endif

cc         if (fixedzarray(i)) then
cc            write(6,*)'fixed ',i,x(i)
cc            z(i)=rzn(i)
cc            if (z(i).lt.0.0d0) z(i)=z(i)+boxlenz
cc            if (z(i).gt.boxlenz) z(i)=z(i)-boxlenz
cc            goto 10
cc         endif

cc         xdisp(i)=xdisp(i)+(rxn(i)-x(i))
cc         ydisp(i)=ydisp(i)+(ryn(i)-y(i))
cc         zdisp(i)=zdisp(i)+(rzn(i)-z(i)) 

         x(i)=rxn(i)
         y(i)=ryn(i)
         z(i)=rzn(i)
c
c Enforce periodic boundary conditions.
c
         if (x(i).lt.0.0d0) x(i)=x(i)+boxlenx
         if (y(i).lt.0.0d0) y(i)=y(i)+boxleny
         if (z(i).lt.0.0d0) z(i)=z(i)+boxlenz

         if (x(i).gt.boxlenx) x(i)=x(i)-boxlenx
         if (y(i).gt.boxleny) y(i)=y(i)-boxleny
         if (z(i).gt.boxlenz) z(i)=z(i)-boxlenz

   10 continue

      endif
c
c Update dipole components and velocities.
c
      if ((dippimlog).or.(quadpimlog)) then

         do 20 i=1,num

            xmu(i)=xmun(i)
            ymu(i)=ymun(i)
            zmu(i)=zmun(i)

            qvelx(i)=qvelnx(i)
            qvely(i)=qvelny(i)
            qvelz(i)=qvelnz(i)

  20     continue

      endif

      return

      end
