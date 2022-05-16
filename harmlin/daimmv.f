      subroutine daimmv
c
c updates the variables describing the dipolar distortions
c of the anions according to the derived equations of motion.
c
      implicit double precision (a-h,o-z)

      include 'common.inc'

      fkeeps=0.0d0
      fkehalf=0.0d0
      epsselfengtot=0.0d0
c
c Calculate the denominator and numerator for the velocity leapfrog
c for the anion radius thermostating:
c
      epszetadenom=1.0d0/(1.0d0+dteps2)
      epszetanum=1.0d0-dteps2

      do 10 i=1,nsp(1)

         epsilonsq=epsilonx(i)*epsilonx(i)+
     x      epsilony(i)*epsilony(i)+epsilonz(i)*epsilonz(i)
         selfexp=selfC*dexp(selfeps*selfeps*epsilonsq)
cc         epsselfeng=selfC*selfeps*selfeps*epsilonsq
         epsselfeng=selfexp-selfC
         epsselfengtot=epsselfengtot+epsselfeng
cc         selfengderiv=2.0d0*selfC*selfeps*selfeps
         selfengderiv=2.0d0*selfC*selfeps*selfeps*selfexp

         epsaccx=-(ftalp(1,2)*engft1dotx(i)
     x             +ftbeta(1,2)*engft2dotx(i)
     x             +ftgamma(1,2)*engft3dotx(i)
     x             +epsilonx(i)*selfengderiv)*epsmassrec
         epsaccy=-(ftalp(1,2)*engft1doty(i)
     x             +ftbeta(1,2)*engft2doty(i)
     x             +ftgamma(1,2)*engft3doty(i)
     x             +epsilony(i)*selfengderiv)*epsmassrec
         epsaccz=-(ftalp(1,2)*engft1dotz(i)
     x             +ftbeta(1,2)*engft2dotz(i)
     x             +ftgamma(1,2)*engft3dotz(i)
     x             +epsilonz(i)*selfengderiv)*epsmassrec

         if (aimtherm) then
            epsvelnx(i)=(epsvelx(i)*epszetanum+dtime*epsaccx)*
     x                   epszetadenom
            epsvelny(i)=(epsvely(i)*epszetanum+dtime*epsaccy)*
     x                   epszetadenom
            epsvelnz(i)=(epsvelz(i)*epszetanum+dtime*epsaccz)*
     x                   epszetadenom
         else
            epsvelnx(i)=epsvelx(i)+dtime*epsaccx
            epsvelny(i)=epsvely(i)+dtime*epsaccy
            epsvelnz(i)=epsvelz(i)+dtime*epsaccz
         endif

         epsilonxn(i)=epsilonx(i)+dtime*epsvelnx(i)
         epsilonyn(i)=epsilony(i)+dtime*epsvelny(i)
         epsilonzn(i)=epsilonz(i)+dtime*epsvelnz(i)

         epsvelavx=(epsvelnx(i)+epsvelx(i))*0.5d0
         epsvelavy=(epsvelny(i)+epsvely(i))*0.5d0
         epsvelavz=(epsvelnz(i)+epsvelz(i))*0.5d0

         fkeeps=fkeeps+(epsvelavx*epsvelavx+epsvelavy*epsvelavy+
     x                   epsvelavz*epsvelavz)
         fkehalf=fkehalf+(epsvelnx(i)*epsvelnx(i)+
     x                    epsvelny(i)*epsvelny(i)+
     x                    epsvelnz(i)*epsvelnz(i))

  10  continue

      fkeeps=fkeeps*epsmass*0.5d0
      fkehalf=fkehalf*epsmass*0.5d0
      teps=fkeeps*2.0d0*gepskin
c
c Epsilon thermostating code:
c
      epszetanew=epszeta+dtime*(fkehalf*2.0d0*gepsrec-1.0d0)*epsrlxsqrec
      dteps2=epszetanew*dtime2

      if (debug) then
         write (110,*) nstep,real(fkeeps),real(epsselfengtot)
      endif

      return

      end
