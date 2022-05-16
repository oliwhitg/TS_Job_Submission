      subroutine dcell(h,b)

c************************************************************
c
c  Calculate the properies of a cell specified by the 
c matrix h(3,3). 
c  The results are returned in the array b, with :
c     b(1 to 3) - lengths of cell vectors
c     b(4 to 6) - cosines of cell angles
c     b(7 to 9) - perpendicular cell widths
c     b(10)     - cell volume
c
c  Modified from CCP5 code.
c
c************************************************************
c
      implicit double precision (a-h,o-z)
      double precision h(3,3),b(10)

c=====================================================================
c Calculate lengths of cell vectors.
      b(1)=sqrt(h(1,1)*h(1,1)+h(2,1)*h(2,1)+h(3,1)*h(3,1))
      b(2)=sqrt(h(1,2)*h(1,2)+h(2,2)*h(2,2)+h(3,2)*h(3,2))
      b(3)=sqrt(h(1,3)*h(1,3)+h(2,3)*h(2,3)+h(3,3)*h(3,3))
 
c=====================================================================
c Calculate the cosines of cell angles.
      b(4)=(h(1,1)*h(1,2)+h(2,1)*h(2,2)+h(3,1)*h(3,2))/(b(1)*b(2))
      b(5)=(h(1,1)*h(1,3)+h(2,1)*h(2,3)+h(3,1)*h(3,3))/(b(1)*b(3))
      b(6)=(h(1,2)*h(1,3)+h(2,2)*h(2,3)+h(3,2)*h(3,3))/(b(2)*b(3))
 
c=====================================================================
c Calculate the vector products of cell vectors....
      axb1=h(2,1)*h(3,2)-h(3,1)*h(2,2)
      axb2=h(3,1)*h(1,2)-h(1,1)*h(3,2)
      axb3=h(1,1)*h(2,2)-h(2,1)*h(1,2)
      bxc1=h(2,2)*h(3,3)-h(3,2)*h(2,3)
      bxc2=h(3,2)*h(1,3)-h(1,2)*h(3,3)
      bxc3=h(1,2)*h(2,3)-h(2,2)*h(1,3)
      cxa1=h(2,3)*h(3,1)-h(2,1)*h(3,3)
      cxa2=h(1,1)*h(3,3)-h(3,1)*h(1,3)
      cxa3=h(2,1)*h(1,3)-h(1,1)*h(2,3)
 
c=====================================================================
c ...and hence the volume of the cell....
      b(10)=abs(h(1,1)*bxc1+h(2,1)*bxc2+h(3,1)*bxc3)
 
c=====================================================================
c ...and finally the perpendicular widths
      b(7)=b(10)/sqrt(bxc1*bxc1+bxc2*bxc2+bxc3*bxc3)
      b(8)=b(10)/sqrt(cxa1*cxa1+cxa2*cxa2+cxa3*cxa3)
      b(9)=b(10)/sqrt(axb1*axb1+axb2*axb2+axb3*axb3)

      return
      end
