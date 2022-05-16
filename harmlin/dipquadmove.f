      subroutine dipquadmove

c***********************************************************
c Updates the dipolar degress of freedom according
c to the derived equations of motion.
c***********************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'

      double precision fkehalf(nspmax)
      double precision fkequadhalf(nspmax)
      double precision zetadenom(nspmax)
      double precision zetanum(nspmax)
      double precision quaddenom(nspmax)
      double precision quadnum(nspmax)

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
         fkequad(i)=0.0d0
         fkequadhalf(i)=0.0d0
         dipquadeng(i)=0.0d0
         dipsqeng(i)=0.0d0
         quadeng(i)=0.0d0

         resengdt(i)=0.0d0
         reseng(i)=0.0d0
c
c Calculate the denominator and numerator for the velocity leapfrog.
c
         zetadenom(i)=1.0d0/(1.0d0+dtdip2(i))
         zetanum(i)=1.0d0-dtdip2(i)
         quaddenom(i)=1.0d0/(1.0d0+dtquad2(i))
         quadnum(i)=1.0d0-dtquad2(i)

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

c
c Dipole equations of motion including cross (B) terms.
c
               dipaccnx=((-2.0d0*xk1(ipoint)*xmu(i))
     x           -2.0d0*xk2(ipoint)*((xmu(i)*quadxx(i))
     x           +ymu(i)*quadxy(i)
     x           +zmu(i)*quadxz(i))
     x           -4.0d0*xk4(ipoint)*xmu(i)*dipsq
     x         +elecx(i))*dipmassrec(ipoint)
               dipaccny=((-2.0d0*xk1(ipoint)*ymu(i))
     x           -2.0d0*xk2(ipoint)*(xmu(i)*quadxy(i)
     x           +(ymu(i)*quadyy(i))
     x           +zmu(i)*quadyz(i))
     x           -4.0d0*xk4(ipoint)*ymu(i)*dipsq
     x         +elecy(i))*dipmassrec(ipoint)
               dipaccnz=((-2.0d0*xk1(ipoint)*zmu(i))
     x           -2.0d0*xk2(ipoint)*(xmu(i)*quadxz(i)
     x           +ymu(i)*quadyz(i)
     x           +(zmu(i)*quadzz(i)))
     x           -4.0d0*xk4(ipoint)*zmu(i)*dipsq
     x         +elecz(i))*dipmassrec(ipoint)

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

               dipsqeng(ipoint)=dipsqeng(ipoint)+dipsq*dipsq

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

c
c Quadrupole equations of motion including cross (B) terms.
c
               quadxxaccn=(xk2(ipoint)*(xmu(i)*xmu(i)
     x                   -dipsq*onethird)
     x             +2.0d0*xk3(ipoint)*quadxx(i)
     x             -onethird*exx(i))*quadmassrec(ipoint)
               quadyyaccn=(xk2(ipoint)*(ymu(i)*ymu(i)
     x                   -dipsq*onethird)
     x             +2.0d0*xk3(ipoint)*quadyy(i)
     x             -onethird*eyy(i))*quadmassrec(ipoint)
               quadzzaccn=(xk2(ipoint)*(zmu(i)*zmu(i)
     x                   -dipsq*onethird)
     x             +2.0d0*xk3(ipoint)*quadzz(i)
     x             -onethird*ezz(i))*quadmassrec(ipoint)
               quadxyaccn=(xk2(ipoint)*xmu(i)*ymu(i)
     x             +2.0d0*xk3(ipoint)*quadxy(i)
     x             -onethird*exy(i))*quadmassrec(ipoint)
               quadxzaccn=(xk2(ipoint)*xmu(i)*zmu(i)
     x             +2.0d0*xk3(ipoint)*quadxz(i)
     x             -onethird*exz(i))*quadmassrec(ipoint)
               quadyzaccn=(xk2(ipoint)*zmu(i)*ymu(i)
     x             +2.0d0*xk3(ipoint)*quadyz(i)
     x             -onethird*eyz(i))*quadmassrec(ipoint)
c
c why is it quadvelnext=quadvel-quadacc*dt?
c
               if (pimtherm(ipoint)) then
c
c thermostated equations of motion:
c
                  quadvelxxn(i)=((quadvelxx(i)*quadnum(ipoint))-
     x               (quadxxaccn*dtime))*quaddenom(ipoint)
                  quadvelyyn(i)=((quadvelyy(i)*quadnum(ipoint))-
     x               (quadyyaccn*dtime))*quaddenom(ipoint)
                  quadvelzzn(i)=((quadvelzz(i)*quadnum(ipoint))-
     x               (quadzzaccn*dtime))*quaddenom(ipoint)
                  quadvelxyn(i)=((quadvelxy(i)*quadnum(ipoint))-
     x               (quadxyaccn*dtime))*quaddenom(ipoint)
                  quadvelxzn(i)=((quadvelxz(i)*quadnum(ipoint))-
     x               (quadxzaccn*dtime))*quaddenom(ipoint)
                  quadvelyzn(i)=((quadvelyz(i)*quadnum(ipoint))-
     x               (quadyzaccn*dtime))*quaddenom(ipoint)

               else

                  quadvelxxn(i)=quadvelxx(i)-quadxxaccn*dtime
                  quadvelyyn(i)=quadvelyy(i)-quadyyaccn*dtime
                  quadvelzzn(i)=quadvelzz(i)-quadzzaccn*dtime
                  quadvelxyn(i)=quadvelxy(i)-quadxyaccn*dtime
                  quadvelxzn(i)=quadvelxz(i)-quadxzaccn*dtime
                  quadvelyzn(i)=quadvelyz(i)-quadyzaccn*dtime

               endif

               quadxxn(i)=quadxx(i)+dtime*quadvelxxn(i)
               quadyyn(i)=quadyy(i)+dtime*quadvelyyn(i)
               quadzzn(i)=quadzz(i)+dtime*quadvelzzn(i)
               quadxyn(i)=quadxy(i)+dtime*quadvelxyn(i)
               quadxzn(i)=quadxz(i)+dtime*quadvelxzn(i)
               quadyzn(i)=quadyz(i)+dtime*quadvelyzn(i)


               fkequadhalf(ipoint)=fkequadhalf(ipoint)
     x               +(quadvelxxn(i)*quadvelxxn(i)
     x               +quadvelyyn(i)*quadvelyyn(i)
     x               +quadvelzzn(i)*quadvelzzn(i)
     x               +quadvelxyn(i)*quadvelxyn(i)
     x               +quadvelxzn(i)*quadvelxzn(i)
     x               +quadvelyzn(i)*quadvelyzn(i))
c
c dipole-quadrupole energy (k2 term):
c
               dipquadeng(ipoint)=dipquadeng(ipoint)
     x           +(xmu(i)*quadxx(i)*xmu(i)
     x           +ymu(i)*quadyy(i)*ymu(i)
     x           +zmu(i)*quadzz(i)*zmu(i)
     x           +2.0d0*(xmu(i)*quadxy(i)*ymu(i)
     x                  +xmu(i)*quadxz(i)*zmu(i)
     x                  +ymu(i)*quadyz(i)*zmu(i)))
c
c quadrupole energy.
c
               quadsq=quadxx(i)*quadxx(i)
     x         +quadyy(i)*quadyy(i)
     x         +quadzz(i)*quadzz(i)
     x         +2.0d0*(quadxy(i)*quadxy(i)
     x                +quadxz(i)*quadxz(i)
     x                +quadyz(i)*quadyz(i))

               quadeng(ipoint)=quadeng(ipoint)+quadsq
c
c Apply the contraint that the quadrupole tensor trace is zero.
c Note that this 'brute force' method is entirely equivalent to
c the lagrange multiplier method as the constraint is linear.
c
         tracetheta3=onethird*(quadxxn(i)+quadyyn(i)+quadzzn(i))
         quadxxn(i)=quadxxn(i)-tracetheta3
         quadyyn(i)=quadyyn(i)-tracetheta3
         quadzzn(i)=quadzzn(i)-tracetheta3
c
c Modify the velocities as well.
c
         quadvelxxn(i)=(quadxxn(i)-quadxx(i))*dtimerec
         quadvelyyn(i)=(quadyyn(i)-quadyy(i))*dtimerec
         quadvelzzn(i)=(quadzzn(i)-quadzz(i))*dtimerec

               xxvelsq=0.25d0*(quadvelxxn(i)+quadvelxx(i))*
     x                  (quadvelxxn(i)+quadvelxx(i))
               yyvelsq=0.25d0*(quadvelyyn(i)+quadvelyy(i))*
     x                  (quadvelyyn(i)+quadvelyy(i))
               zzvelsq=0.25d0*(quadvelzzn(i)+quadvelzz(i))*
     x                  (quadvelzzn(i)+quadvelzz(i))
               xyvelsq=0.25d0*(quadvelxyn(i)+quadvelxy(i))*
     x                  (quadvelxyn(i)+quadvelxy(i))
               xzvelsq=0.25d0*(quadvelxzn(i)+quadvelxz(i))*
     x                  (quadvelxzn(i)+quadvelxz(i))
               yzvelsq=0.25d0*(quadvelyzn(i)+quadvelyz(i))*
     x                  (quadvelyzn(i)+quadvelyz(i))

               fkequad(ipoint)=fkequad(ipoint)
     x      +(xxvelsq+yyvelsq+zzvelsq
     x      +xyvelsq+xzvelsq+yzvelsq)

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

            fkequad(i)=fkequad(i)*quadmass(i)
            fkequadhalf(i)=fkequadhalf(i)*quadmass(i)

            dipsqeng(i)=dipsqeng(i)*xk4(i)
            dipquadeng(i)=dipquadeng(i)*xk2(i)
            quadeng(i)=quadeng(i)*xk3(i)

c
c Thermostating code:
c
            tdip(i)=fke(i)*gdipkin(i)
            tquad(i)=fkequad(i)*gquadkin(i)

            tfluxdip(i)=tfluxdip(i)+(tdip(i)*gdiptim(i)*dipzeta(i))
cc            tfluxquad=tfluxquad+(tquad*gquadtim*quadzeta)

            dipzetanew(i)=dipzeta(i)+dtime
     x           *(((fkehalf(i)*gdiprec(i))-1.0d0)*diprlxsqrec(i))
            dtdip2(i)=dipzetanew(i)*dtime2

            quadzetanew(i)=quadzeta(i)+dtime
     x        *(((fkequadhalf(i)*gquadrec(i))-1.0d0)*diprlxsqrec(i))
            dtquad2(i)=quadzetanew(i)*dtime2

            fke(i)=fke(i)*0.5d0
            fkehalf(i)=fkehalf(i)*0.5d0

            fkequad(i)=fkequad(i)*0.5d0
            fkequadhalf(i)=fkequadhalf(i)*0.5d0

            fketot=fketot+fke(i)+fkequad(i)

         endif

         if (.not.pimtherm(i)) then
            dipzeta(i)=0.0d0
            dipzetanew(i)=0.0d0
            dtdip2(i)=0.0d0
         endif

   20 continue

      return

      end
