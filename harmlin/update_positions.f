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

         if (fixedzarray(i)) then
            write(6,*)'fixed ',i,x(i)
            z(i)=rzn(i)
            if (z(i).lt.0.0d0) z(i)=z(i)+boxlenz
            if (z(i).gt.boxlenz) z(i)=z(i)-boxlenz
            goto 10
         endif

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

      if (quadpimlog) then

         do 25 i=1,num

            quadxx(i)=quadxxn(i)
            quadyy(i)=quadyyn(i)
            quadzz(i)=quadzzn(i)
            quadxy(i)=quadxyn(i)
            quadxz(i)=quadxzn(i)
            quadyz(i)=quadyzn(i)

            quadvelxx(i)=quadvelxxn(i)
            quadvelyy(i)=quadvelyyn(i)
            quadvelzz(i)=quadvelzzn(i)
            quadvelxy(i)=quadvelxyn(i)
            quadvelxz(i)=quadvelxzn(i)
            quadvelyz(i)=quadvelyzn(i)

  25     continue

      endif

      if ((cimlog).or.(daimlog).or.(quaimlog)) then

         do 30 i=1,num

            delta(i)=deltan(i)
            deltavel(i)=deltaveln(i)

  30    continue

      endif

      if ((daimlog).or.(quaimlog)) then

         do 300 i=1,num

            epsilonx(i)=epsilonxn(i)
            epsilony(i)=epsilonyn(i)
            epsilonz(i)=epsilonzn(i)
            epsvelx(i)=epsvelnx(i)
            epsvely(i)=epsvelny(i)
            epsvelz(i)=epsvelnz(i)

 300     continue

      endif

      if (quaimlog) then

         do 330 i=1,num

            quaimxx(i)=quaimxxn(i)
            quaimyy(i)=quaimyyn(i)
            quaimzz(i)=quaimzzn(i)
            quaimxy(i)=quaimxyn(i)
            quaimxz(i)=quaimxzn(i)
            quaimyz(i)=quaimyzn(i)

            quaimvelxx(i)=quaimvelxxn(i)
            quaimvelyy(i)=quaimvelyyn(i)
            quaimvelzz(i)=quaimvelzzn(i)
            quaimvelxy(i)=quaimvelxyn(i)
            quaimvelxz(i)=quaimvelxzn(i)
            quaimvelyz(i)=quaimvelyzn(i)

 330     continue

      endif

      return

      end
