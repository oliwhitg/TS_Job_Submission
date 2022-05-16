      subroutine quaimmv
c
c updates the variables describing the quadrupolar distortions
c of the anions according to the derived equations of motion.
c
      implicit double precision (a-h,o-z)

      include 'common.inc'

      fkequaim=0.0d0
      quaimselfengtot=0.0d0

      do 10 i=1,nsp(1)

         thetasq=quaimxx(i)*quaimxx(i)+
     x           quaimyy(i)*quaimyy(i)+
     x           quaimzz(i)*quaimzz(i)+
     x           2.0d0*(quaimxy(i)*quaimxy(i)+
     x                  quaimxz(i)*quaimxz(i)+
     x                  quaimyz(i)*quaimyz(i))

cc         quaimselfeng=selfH*selfquaim*selfquaim*thetasq
         selfexp=selfH*exp(selfquaim*selfquaim*thetasq)
         quaimselfeng=selfexp-selfH
         quaimselfengtot=quaimselfengtot+quaimselfeng
         selfengderiv=2.0d0*selfH*selfquaim*selfquaim*selfexp

         quaimaccxx=-(ftalp(1,2)*engftdotxx(i,1)
     x               +ftbeta(1,2)*engftdotxx(i,2)
     x               +ftgamma(1,2)*engftdotxx(i,3)
     x               +quaimxx(i)*selfengderiv)*quaimmassrec
         quaimaccyy=-(ftalp(1,2)*engftdotyy(i,1)
     x               +ftbeta(1,2)*engftdotyy(i,2)
     x               +ftgamma(1,2)*engftdotyy(i,3)
     x               +quaimyy(i)*selfengderiv)*quaimmassrec
         quaimacczz=-(ftalp(1,2)*engftdotzz(i,1)
     x               +ftbeta(1,2)*engftdotzz(i,2)
     x               +ftgamma(1,2)*engftdotzz(i,3)
     x               +quaimzz(i)*selfengderiv)*quaimmassrec

         quaimaccxy=-(ftalp(1,2)*engftdotxy(i,1)
     x               +ftbeta(1,2)*engftdotxy(i,2)
     x               +ftgamma(1,2)*engftdotxy(i,3)
     x               +2.0d0*quaimxy(i)*selfengderiv)*quaimmassrec
         quaimaccxz=-(ftalp(1,2)*engftdotxz(i,1)
     x               +ftbeta(1,2)*engftdotxz(i,2)
     x               +ftgamma(1,2)*engftdotxz(i,3)
     x               +2.0d0*quaimxz(i)*selfengderiv)*quaimmassrec
         quaimaccyz=-(ftalp(1,2)*engftdotyz(i,1)
     x               +ftbeta(1,2)*engftdotyz(i,2)
     x               +ftgamma(1,2)*engftdotyz(i,3)
     x               +2.0d0*quaimyz(i)*selfengderiv)*quaimmassrec

         quaimvelxxn(i)=quaimvelxx(i)+dtime*quaimaccxx
         quaimvelyyn(i)=quaimvelyy(i)+dtime*quaimaccyy
         quaimvelzzn(i)=quaimvelzz(i)+dtime*quaimacczz

         quaimvelxyn(i)=quaimvelxy(i)+dtime*quaimaccxy
         quaimvelxzn(i)=quaimvelxz(i)+dtime*quaimaccxz
         quaimvelyzn(i)=quaimvelyz(i)+dtime*quaimaccyz

         quaimxxn(i)=quaimxx(i)+dtime*quaimvelxxn(i)
         quaimyyn(i)=quaimyy(i)+dtime*quaimvelyyn(i)
         quaimzzn(i)=quaimzz(i)+dtime*quaimvelzzn(i)

         quaimxyn(i)=quaimxy(i)+dtime*quaimvelxyn(i)
         quaimxzn(i)=quaimxz(i)+dtime*quaimvelxzn(i)
         quaimyzn(i)=quaimyz(i)+dtime*quaimvelyzn(i)
c
c Apply the constraint that the quadrupole tensor trace is zero.
c
         tracetheta3=onethird*(quaimxxn(i)+quaimyyn(i)+quaimzzn(i))
         quaimxxn(i)=quaimxxn(i)-tracetheta3
         quaimyyn(i)=quaimyyn(i)-tracetheta3
         quaimzzn(i)=quaimzzn(i)-tracetheta3
c
c Modify the velocities as well.
c
         quaimvelxxn(i)=(quaimxxn(i)-quaimxx(i))*dtimerec
         quaimvelyyn(i)=(quaimyyn(i)-quaimyy(i))*dtimerec
         quaimvelzzn(i)=(quaimzzn(i)-quaimzz(i))*dtimerec

         xxvelsq=0.25d0*(quaimvelxxn(i)+quaimvelxx(i))*
     x                  (quaimvelxxn(i)+quaimvelxx(i))
         yyvelsq=0.25d0*(quaimvelyyn(i)+quaimvelyy(i))*
     x                  (quaimvelyyn(i)+quaimvelyy(i))
         zzvelsq=0.25d0*(quaimvelzzn(i)+quaimvelzz(i))*
     x                  (quaimvelzzn(i)+quaimvelzz(i))
         xyvelsq=0.25d0*(quaimvelxyn(i)+quaimvelxy(i))*
     x                  (quaimvelxyn(i)+quaimvelxy(i))
         xzvelsq=0.25d0*(quaimvelxzn(i)+quaimvelxz(i))*
     x                  (quaimvelxzn(i)+quaimvelxz(i))
         yzvelsq=0.25d0*(quaimvelyzn(i)+quaimvelyz(i))*
     x                  (quaimvelyzn(i)+quaimvelyz(i))

         fkequaim=fkequaim+(xxvelsq+yyvelsq+zzvelsq
     x      +xyvelsq+xzvelsq+yzvelsq)

  10  continue

      fkequaim=fkequaim*quaimmass*0.5d0
      tquaim=fkequaim*2.0d0*gquaimkin

      return

      end
