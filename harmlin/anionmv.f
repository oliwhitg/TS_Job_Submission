      subroutine anionmv

c*************************************************************
c Updates the anion radii according to the derived 
c derived equations of motion.
c*************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'
c
c Zero the energy accumulators:
c
      fkeself=0.0d0
      fkehalf=0.0d0
      selfengtot=0.0d0
c
c Calculate the denominator and numerator for the velocity leapfrog
c for the anion radius thermostating:
c
      deltazetadenom=1.0d0/(1.0d0+dtdelta2)
      deltazetanum=1.0d0-dtdelta2

c TMP - CIM
      do 10 i=1,nsp(1)

         selfeng=(selfB*(delta(i)**2.0d0))
     x          +(selfgam*(delta(i)**3.0d0))
     x          +(selfgam2*(delta(i)**4.0d0))
     x          +(selfgam3*(delta(i)**5.0d0))
         deltasq=delta(i)*delta(i)
         selfengtot=selfengtot+selfeng
c
c This is wrong if we have more than a single exponential for
c the O-O interaction.
c
         deltaaccn=-(ftalp(1,2)*engft1(i)
     x             +ftbeta(1,2)*engft2(i)
     x             +ftgamma(1,2)*engft3(i)
     x		+((2.0d0*selfB*(delta(i)))
     x          +(3.0d0*selfgam*(delta(i)**2.0d0))
     x          +(4.0d0*selfgam2*(delta(i)**3.0d0))
     x          +(5.0d0*selfgam3*(delta(i)**4.0d0))))
     x		*selfmassrec

         if (cimtherm) then
            deltaveln(i)=((deltavel(i)*deltazetanum)
     x                  +(deltaaccn*dtime))*deltazetadenom
         else
            deltaveln(i)=deltavel(i)+deltaaccn*dtime
         endif

         deltan(i)=delta(i)+dtime*deltaveln(i)

         deltavelav=(deltavel(i)+deltaveln(i))*0.5d0

         fkeself=fkeself+(deltavelav*deltavelav)
         fkehalf=fkehalf+deltaveln(i)*deltaveln(i)

  10  continue

      fkeself=fkeself*selfmass*0.5d0
      fkehalf=fkehalf*selfmass*0.5d0
c
c Anion radius thermostating code:
c
      tdelta=fkeself*2.0d0*gdeltakin

      deltazetanew=deltazeta+dtime
     x            *(((fkehalf*2.0d0*gdeltarec)-1.0d0)
     x            *deltarlxsqrec)
      dtdelta2=deltazetanew*dtime2

      return

      end
