      subroutine dipolemove

c***********************************************************
c Updates the dipolar degress of freedom according
c to the derived equations of motion.
c***********************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'

      double precision fkehalf(nspmax)
      double precision zetadenom(nspmax)
      double precision zetanum(nspmax)

c local array for polarization energy at t+dt.
      double precision resengdt(nspmax)
      double precision fkedt(nspmax)
c
c Zero the ficticious k.e. and restoring energy accumulators.
c

      do 100 i=1,nspec
         fke(i)=0.0d0
         fkedt(i)=0.0d0
         fkehalf(i)=0.0d0
         resengdt(i)=0.0d0
         reseng(i)=0.0d0
c
c Calculate the denominator and numerator for the velocity leapfrog.
c
         zetadenom(i)=1.0d0/(1.0d0+dtdip2(i))
         zetanum(i)=1.0d0-dtdip2(i)
  100 continue

      l=1
      do 30 j=1,nspec

         if (polarizablelog(j)) then

            do 10 i=l,nsp(j)+l-1

               ipoint=ntype(i)

               dipsq=(xmu(i)*xmu(i))
     x              +(ymu(i)*ymu(i))
     x              +(zmu(i)*zmu(i))

               reseng(ipoint)=reseng(ipoint)+dipsq

               dipaccnx=((-2.0d0*xk1(ipoint)*xmu(i))
     x               +elecx(i))*dipmassrec(ipoint)
               dipaccny=((-2.0d0*xk1(ipoint)*ymu(i))
     x               +elecy(i))*dipmassrec(ipoint)
               dipaccnz=((-2.0d0*xk1(ipoint)*zmu(i))
     x               +elecz(i))*dipmassrec(ipoint)

               if (pimtherm(ipoint)) then

                  qvelnx(i)=((qvelx(i)*zetanum(ipoint))
     x      		+(dipaccnx*dtime))*zetadenom(ipoint)
                  qvelny(i)=((qvely(i)*zetanum(ipoint))
     x		      +(dipaccny*dtime))*zetadenom(ipoint)
                  qvelnz(i)=((qvelz(i)*zetanum(ipoint))
     x       		+(dipaccnz*dtime))*zetadenom(ipoint)

               else

                  qvelnx(i)=qvelx(i)+dipaccnx*dtime
                  qvelny(i)=qvely(i)+dipaccny*dtime
                  qvelnz(i)=qvelz(i)+dipaccnz*dtime

c guess the dipole velocities at t+(3/2)dt by assuming
c that the acceleration at time t+dt = that at time t.
cc            qvelnnx=qvelnx(i)+dipaccnx*dtime
cc            qvelnny=qvelny(i)+dipaccnx*dtime
cc            qvelnnz=qvelnz(i)+dipaccnx*dtime
                  qvelnnx=qvelnx(i)+dipaccnx*dtime*0.5d0
                  qvelnny=qvelny(i)+dipaccnx*dtime*0.5d0
                  qvelnnz=qvelnz(i)+dipaccnx*dtime*0.5d0

               endif
c
c Get fke at t+half dt for the friction coeff:
c
               fkehalf(ipoint)=fkehalf(ipoint)
     x                  +(qvelnx(i)*qvelnx(i))
     x                  +(qvelny(i)*qvelny(i))
     x                  +(qvelnz(i)*qvelnz(i))

               xmun(i)=xmu(i)+dtime*qvelnx(i)
               ymun(i)=ymu(i)+dtime*qvelny(i)
               zmun(i)=zmu(i)+dtime*qvelnz(i)

c polarization energy at t+dt
               dipsq=(xmun(i)*xmun(i))
     x              +(ymun(i)*ymun(i))
     x              +(zmun(i)*zmun(i))

               resengdt(ipoint)=resengdt(ipoint)+dipsq

               xvelsq=0.25d0*(qvelnx(i)+qvelx(i))
     x                      *(qvelnx(i)+qvelx(i))
               yvelsq=0.25d0*(qvelny(i)+qvely(i))
     x                      *(qvelny(i)+qvely(i))
               zvelsq=0.25d0*(qvelnz(i)+qvelz(i))
     x                      *(qvelnz(i)+qvelz(i))
               fke(ipoint)=fke(ipoint)+xvelsq+yvelsq+zvelsq

c approximate fke at t+dt 
               xvelsq=0.25d0*(qvelnx(i)+qvelnnx)
     x                      *(qvelnx(i)+qvelnnx)
               yvelsq=0.25d0*(qvelny(i)+qvelnny)
     x                      *(qvelny(i)+qvelnny)
               zvelsq=0.25d0*(qvelnz(i)+qvelnnz)
     x                      *(qvelnz(i)+qvelnnz)
               xvelsq=(qvelnnx)*(qvelnnx)
               yvelsq=(qvelnny)*(qvelnny)
               zvelsq=(qvelnnz)*(qvelnnz)
               fkedt(ipoint)=fkedt(ipoint)+xvelsq+yvelsq+zvelsq

   10       continue

         endif

c update counter for ion positions for inner loop.
         l=l+nsp(j)

   30 continue

      fketot=0.0d0

      do 20 i=1,nspec

         if (polarizablelog(i)) then

            fke(i)=fke(i)*dipmass(i)
            fkehalf(i)=fkehalf(i)*dipmass(i)
            reseng(i)=reseng(i)*xk1(i)

            fketot=fketot+fke(i)
c
c Thermostating code:
c
            tdip(i)=fke(i)*gdipkin(i)

            tfluxdip(i)=tfluxdip(i)+(tdip(i)*gdiptim(i)*dipzeta(i))

            dipzetanew(i)=dipzeta(i)+dtime
     x           *(((fkehalf(i)*gdiprec(i))-1.0d0)*diprlxsqrec(i))
            dtdip2(i)=dipzetanew(i)*dtime2

            fke(i)=fke(i)*0.5d0
            fkehalf(i)=fkehalf(i)*0.5d0

         endif

         if (.not.pimtherm(i)) then
            dipzeta(i)=0.0d0
            dipzetanew(i)=0.0d0
            dtdip2(i)=0.0d0
         endif

   20 continue

      return

      end
