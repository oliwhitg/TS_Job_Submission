      subroutine matinv(h,hi,deterh)
      implicit double precision (a-h,o-z)
      double precision h(3,3),hi(3,3)

c      include 'common.inc'

c....calculate adjoint matrix

      cofh11=h(2,2)*h(3,3)-h(2,3)*h(3,2)
      cofh12=h(2,1)*h(3,3)-h(2,3)*h(3,1)
      cofh13=h(2,1)*h(3,2)-h(2,2)*h(3,1)
      cofh21=h(1,2)*h(3,3)-h(1,3)*h(3,2)
      cofh22=h(1,1)*h(3,3)-h(1,3)*h(3,1)
      cofh23=h(1,1)*h(3,2)-h(1,2)*h(3,1)
      cofh31=h(1,2)*h(2,3)-h(1,3)*h(2,2)
      cofh32=h(1,1)*h(2,3)-h(1,3)*h(2,1)
      cofh33=h(1,1)*h(2,2)-h(1,2)*h(2,1)

      deterh=(h(1,1)*cofh11)-(h(1,2)*cofh12)+(h(1,3)*cofh13)
      r=0.0d0

      if (abs(deterh).gt.0.0d0) r=1.0d0/deterh

      hi(1,1)=r*cofh11
      hi(1,2)=-r*cofh21
      hi(1,3)=r*cofh31
      hi(2,1)=-r*cofh12
      hi(2,2)=r*cofh22
      hi(2,3)=-r*cofh32
      hi(3,1)=r*cofh13
      hi(3,2)=-r*cofh23
      hi(3,3)=r*cofh33

      return
      end

