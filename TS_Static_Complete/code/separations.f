      subroutine separations

c************************************************************
c Performs an energy minimilisation on the startup configuration
c via a method analogous to the conjugant gradient method.
c Splits the Pim and Cim parts of the minimisation up, as they are
c independent of each other.
c************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc' 

      do 10 i=2,num

         do 20 j=1,i-1
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
     x             +hlab2(1,2)*dycf+hlab2(1,3)*dzcf
            dysav(i,j)=hlab2(2,1)*dxcf
     x             +hlab2(2,2)*dycf+hlab2(2,3)*dzcf
            dzsav(i,j)=hlab2(3,1)*dxcf
     x             +hlab2(3,2)*dycf+hlab2(3,3)*dzcf

            dxsav(j,i)=-dxsav(i,j)
            dysav(j,i)=-dysav(i,j)
            dzsav(j,i)=-dzsav(i,j)

   20    continue
 
   10 continue

      return

      end

