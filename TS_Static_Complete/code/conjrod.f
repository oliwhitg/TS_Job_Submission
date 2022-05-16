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
cc         call separations
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
            fkepimprev=0.0d0

   30    continue

         updatepim=.false.

      endif

      if (updatepim) then

         do 50 i=1,num

            qvelx(i)=qvelnx(i)
            qvely(i)=qvelny(i)
            qvelz(i)=qvelnz(i)
            xmu(i)=xmun(i)
            ymu(i)=ymun(i)
            zmu(i)=zmun(i)

  50     continue

      endif

      prestep=prestep+1
      write(6,*)prestep,(fke(i),i=1,nspec)
      write(6,*)prestep,real(fke(1)),real(fkeself)

      if (prestep.le.nrodstep) goto 10

      write (6,*)
      write (6,*) '**** Anion annealing completed ****'
      write (6,*)

      return

      end
