      subroutine conjrod

c************************************************************
c Performs an energy minimilisation on the startup configuration
c via a method analogous to the conjugant gradient method.
c Splits the Pim and Cim parts of the minimisation up, as they are
c independent of each other.
c************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc' 

      logical updatepim,updatecim
      integer prestep

      updatepim=.true.
      updatecim=.true.

      prestep=1
      fkepimprev=0.0d0
      fkecimprev=0.0d0
c      open(99,file='anneal1.out',status='new')
c      open(100,file='anneal2.out',status='new')
c      open(97,file='anneal3.out',status='new')
  10  if (annealrad) then
cc         call rearrange
         call separations
         call sr_energy
      else
         call ener
      endif

c
c Stop thermostating for the annealing procedure.
c
      do 112 i=1,nspec
         dipzeta(i)=0.0d0
         dipzetanew(i)=0.0d0
         dtdip2(i)=0.0d0
 112  continue

      deltazeta=0.0d0
      dtdelta2=0.0d0
      epszeta=0.0d0
      dteps2=0.0d0
      quaimzeta=0.0d0
      dtquaim2=0.0d0

      call dofmoveselect

      write(6,*)prestep,fketot
      if ((fketot).ge.fkepimprev) then

         fkepimprev=fketot
         updatepim=.true.

      else
c
c Quench pim-based fake velocities:
c
         do 30 i=1,num

            qvelx(i)=0.0d0
            qvely(i)=0.0d0
            qvelz(i)=0.0d0
            quadvelxx(i)=0.0d0
            quadvelyy(i)=0.0d0
            quadvelzz(i)=0.0d0
            quadvelxy(i)=0.0d0
            quadvelxz(i)=0.0d0
            quadvelyz(i)=0.0d0
            fkepimprev=0.0d0

   30    continue

         updatepim=.false.

      endif

c CIM anneal if required.
      if ((cimlog).or.(daimlog).or.(quaimlog)) then

         if ((fkeself+fkeeps+fkequaim).ge.fkecimprev) then

            fkecimprev=fkeself+fkeeps+fkequaim
            updatecim=.true.

         else
c
c Quench cim-based fake velocities:
c
            do i=1,nsp(1)
               deltavel(i)=0.0d0
               epsvelx(i)=0.0d0
               epsvely(i)=0.0d0
               epsvelz(i)=0.0d0
               quaimvelxx(i)=0.0d0
               quaimvelyy(i)=0.0d0
               quaimvelzz(i)=0.0d0
               quaimvelxy(i)=0.0d0
               quaimvelxz(i)=0.0d0
               quaimvelyz(i)=0.0d0
            enddo

            fkecimprev=0.0d0
            updatecim=.false.

         endif

      endif

      if (updatepim) then

         do 50 i=1,num

            qvelx(i)=qvelnx(i)
            qvely(i)=qvelny(i)
            qvelz(i)=qvelnz(i)
            xmu(i)=xmun(i)
            ymu(i)=ymun(i)
            zmu(i)=zmun(i)

            quadvelxx(i)=quadvelxxn(i)
            quadvelyy(i)=quadvelyyn(i)
            quadvelzz(i)=quadvelzzn(i)
            quadvelxy(i)=quadvelxyn(i)
            quadvelxz(i)=quadvelxzn(i)
            quadvelyz(i)=quadvelyzn(i)
            quadxx(i)=quadxxn(i)
            quadyy(i)=quadyyn(i)
            quadzz(i)=quadzzn(i)
            quadxy(i)=quadxyn(i)
            quadxz(i)=quadxzn(i)
            quadyz(i)=quadyzn(i)

  50     continue

      endif

      if ((cimlog).or.(daimlog).or.(quaimlog)) then

         if (updatecim) then

            do 40 i=1,nsp(1)

               delta(i)=deltan(i)
               deltavel(i)=deltaveln(i)

               epsilonx(i)=epsilonxn(i)
               epsilony(i)=epsilonyn(i)
               epsilonz(i)=epsilonzn(i)
               epsvelx(i)=epsvelnx(i)
               epsvely(i)=epsvelny(i)
               epsvelz(i)=epsvelnz(i)

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

  40        continue

         endif

      endif

      prestep=prestep+1
      write(6,*)prestep,(fke(i),i=1,nspec)
      write(6,*)prestep,real(fke(1)),real(fkeself),
     x             real(fkeeps),real(fkequaim)

      write(6,*) prestep,(selfengtot),(delta(1))
     x   ,(engft1(1)),(engft2(1)),(engft3(1))

      if (prestep.le.nrodstep) goto 10

c      do 111 i=1,num
c         write(15,*)i,real(delta(i))
c         write(17,*)i,real(xmu(i)),real(ymu(i)),real(zmu(i))
c         write(19,*)i,real(quadxx(i)),real(quadyy(i)),real(quadzz(i))
c         write(19,*)i,real(quadxy(i)),real(quadxz(i)),real(quadyz(i))
c         write(18,*)i,real(elecx(i)),real(elecy(i))
c     x                 ,real(elecz(i))
c         write(16,*)i,real(epsilonx(i)),real(epsilony(i))
c     x                ,real(epsilonz(i))
c         write(13,*)real(quaimxx(i)),real(quaimyy(i)),real(quaimzz(i))
c         write(14,*)real(quaimxy(i)),real(quaimxz(i)),real(quaimyz(i))
c  111 continue

c      close(99)
c      close(100)
c      close(97)

      write (6,*)
      write (6,*) '**** Anion annealing completed ****'
      write (6,*)

      return

      end
