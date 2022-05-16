      subroutine readin_lj

c************************************************************
c   Reads in the data from the .inpt file.
c************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'
      include 'lj.inc'

      character*80 rstinlj

c=====================================================================
c Input I/O filenames.
      open(20,file='lj.inpt',status='old')

         read(20,*)numlj
         do 10 i=1,nspec
            read(20,*)epsilon(i)
            read(20,*)sigma(i)
   10    continue
         read(20,'(a)')rstinlj

      close(20)

c  convert to au.
      rkb=   3.16669D-6
      do 20 i=1,nspec
         sigma(i)=sigma(i)/0.529177d0

         epsilon(i)=epsilon(i)*rkb

         sig12(i)=sigma(i)**12.0d0
         sig6(i)=sigma(i)**6.0d0

         fourepsilcyl(i)=4.0d0*epsilon(i)

   20 continue

      write(6,*)
      write(6,*)'**** LJ Run-time parameters read in ****'
      write(6,*)

      return

      end
