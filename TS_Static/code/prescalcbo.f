      subroutine prescalc

c***************************************************************
c calculates the pressure at each step in the PIM
c***************************************************************

      implicit double precision (a-h,o-z)

      include 'common.inc'

      xnumdens=dble(num)/cellvol
      xkbtterm=xnumdens*boltz*tke

c      write(75,*)'virials',srvirft,srdipvir,virtot
c      write(75,*)
c      write(75,*)'erf',erfcracc
c      write(75,*)'qdip',qdipacc
c      write(75,*)'dipdip',dipdipacc
c      write(75,*)'qquad',qquadacc
c      write(75,*)'dipquad',dipquadacc
c      write(75,*)'quadquad',quadquadacc
c      write(75,*)'toteng',toteng
c      write(75,*)

        if (outp) then

      srvirft=srvirft/(cellvol*3.0d0)
      srdipvir=srdipvir/(cellvol*3.0d0)

         endif

      rtotal=(erfcracc+qdipacc+dipdipacc+qquadacc+dipquadacc+
     x        quadquadacc+toteng)*cellvol3rec

c      write(76,*)
c      write(76,*)'ideal term =',xkbtterm
c      write(76,*)
c      write(76,*)'s-r f-t contrib. =',srvirft
c      write(76,*)'s-r dip contrib. =',srdipvir
c      write(76,*)'total s-r contrib. =',srvirft+srdipvir
c      write(76,*)
c      write(76,*)'erf contrib. =',qqacc*cellvol3rec
c      write(76,*)'qdip contrib. =',qdipacc*cellvol3rec
c      write(76,*)'dipdip contrib. =',dipdipacc*cellvol3rec
c      write(76,*)'qquad l-r contrib. =',qquadacc*cellvol3rec
c      write(76,*)'dipquad l-r contrib. =',dipquadacc*cellvol3rec
c      write(76,*)'quadquad l-r contrib. =',quadquadacc*cellvol3rec
c      write(76,*)'engpetot contrib. =',toteng*cellvol3rec
c      write(76,*)'total l-r contrib. =',rtotal
c      write(76,*)

        if (outp) then 

        if (mod(nstep,npercell).eq.0) then

      write(587,*) nstep,xkbtterm,srvirft

      write(589,*) nstep,(srdipvir+(qdipacc/(cellvol*3.0d0)))
      write(590,*) nstep,dipdipacc/(cellvol*3.0d0),
     x ((erfcracc+toteng)/(cellvol*3.0d0))

	endif

        endif



      encoul=erfcracc+qdipacc+dipdipacc+qquadacc+dipquadacc+
     x       quadquadacc+toteng

      srcoulterm1=virtot*cellvol3rec

      srcoulterm1x=virtotx/(cellvol*3.0d0)
      srcoulterm1y=virtoty/(cellvol*3.0d0)
      srcoulterm1z=virtotz/(cellvol*3.0d0)


      srcoulterm2=encoul*cellvol3rec

c      write(75,*)'step number',nstep
c      write(75,*)'coulombic energy',encoul
c      write(75,*)'kbt',xkbtterm

c      write(75,*)'virial pressure',srcoulterm1
c      write(75,*)'coulombic pressure',srcoulterm2

      presfinal=xkbtterm+srcoulterm1+srcoulterm2

      presfinalx=(xkbtterm)+srcoulterm2+
     x (3.0d0*srcoulterm1x)

      presfinaly=xkbtterm+srcoulterm2+
     x  (3.0d0*srcoulterm1y)

      presfinalz=xkbtterm+srcoulterm2+
     x  (3.0d0*srcoulterm1z)




	if (outp) then

        if (mod(nstep,npercell).eq.0) then

c      write(75,*)'final pressure',presfinal

      write(74,*) nstep,presfinal



      write(6,*) 'pressure=',presfinal

	endif

	endif

      return
      end
