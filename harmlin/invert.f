      subroutine invert(h,hi,d)

c************************************************************
c
c  Amazingly unsubtle matrix inversion algorithm
c removing the need for nag f01abf in the inversion
c of the short-range force matrix.
c
c  h(3,3) is the cell vector matrix.
c  hi(3,3) is the inverse of the above.
c  d is the determinant of h().
c
c  Modified from CCP5 code.
c
c************************************************************

      implicit double precision (a-h,o-z)
      double precision h(3,3),hi(3,3)

c	write(999,*) 'h',h
c
c Calculate adjoint matrix.
c
      hi(1,1)=h(2,2)*h(3,3)-h(3,2)*h(2,3)
      hi(2,1)=h(3,1)*h(2,3)-h(2,1)*h(3,3)
      hi(3,1)=h(2,1)*h(3,2)-h(3,1)*h(2,2)
      hi(1,2)=h(3,2)*h(1,3)-h(1,2)*h(3,3)
      hi(2,2)=h(1,1)*h(3,3)-h(3,1)*h(1,3)
      hi(3,2)=h(3,1)*h(1,2)-h(1,1)*h(3,2)
      hi(1,3)=h(1,2)*h(2,3)-h(2,2)*h(1,3)
      hi(2,3)=h(2,1)*h(1,3)-h(1,1)*h(2,3)
      hi(3,3)=h(1,1)*h(2,2)-h(2,1)*h(1,2)
c
c Calculate determinant.
c
      d=h(1,1)*hi(1,1)+h(1,2)*hi(2,1)+h(1,3)*hi(3,1)
      r=0.d0
      if(abs(d).gt.0.d0)r=1.d0/d
c
c Complete inverse matrix.
c
      hi(1,1)=r*hi(1,1)
      hi(2,1)=r*hi(2,1)
      hi(3,1)=r*hi(3,1)
      hi(1,2)=r*hi(1,2)
      hi(2,2)=r*hi(2,2)
      hi(3,2)=r*hi(3,2)
      hi(1,3)=r*hi(1,3)
      hi(2,3)=r*hi(2,3)
      hi(3,3)=r*hi(3,3)

      return
      end
