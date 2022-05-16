      subroutine boxreset

c***********************************************************
c
c Sets/resets all box length dependent variables.
c
c***********************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'
c
c Local double precision arrays.
      double precision lengths(3)

      lengths(1)=boxlenx
      lengths(2)=boxleny
      lengths(3)=boxlenz

      call dcell(h,dcellinfo)

      cellvol=dcellinfo(10)*boxlenx*boxleny*boxlenz
      halfboxx=boxlenx/2.0d0
      halfboxxrec=1.0d0/halfboxx
      halfboxy=boxleny/2.0d0
      halfboxyrec=1.0d0/halfboxy
      halfboxz=boxlenz/2.0d0
      halfboxzrec=1.0d0/halfboxz
      twopiboxx=twopi/boxlenx
      twopiboxy=twopi/boxleny
      twopiboxz=twopi/boxlenz
      fourpicell=fourpi/cellvol
      cellvol3rec=onethird/cellvol
c
c Find the minimum boxlen:
c
      boxlen=boxlenx
      if (boxleny.lt.boxlen) boxlen=boxleny
      if (boxlenz.lt.boxlen) boxlen=boxlenz

      halfbox=boxlen*0.5d0
      boxlenrec=1.0d0/boxlen
      twopibox=twopi/boxlen
c
c Now invert the cell matrix for use with the force calculations
c and the Ewald summation.
c
      call invert(h,hi,det)

      do 10 i=1,3
         do 10 j=1,3
            fullh(i,j)=h(i,j)*lengths(j)
 10   continue

      call invert(fullh,fullhi,fulldet)

c
c Calculate factor for reciprocal space forces.
c
      fac=eightpi/cellvol

      call dcell(fullh,bh)

      halfboxminsq=(0.5d0*dmin1(bh(7),bh(8),bh(9)))**2

      rsqmax=rcut**2.0d0

c rsqmax now directly from .inpt. halfboxminsq from
c the minimum lengths in the simulation cell.
      if (rsqmax.gt.halfboxminsq) write(6,*)halfboxminsq
     x          ,rsqmax

      return

      end
