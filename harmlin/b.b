cc Harmonic potential
      if ((harmpotlog))then
         call separations
         pi=4.0d0*atan(1.0d0)

         if (.not.mixlog) then

            stsrxx=0.0d0
            stsrxy=0.0d0
            stsrxz=0.0d0
            stsryy=0.0d0
            stsryz=0.0d0
            stsrzz=0.0d0

            engpetot=0.0d0

            do 995 i=1,num
               frrx(i)=0.0d0
               frry(i)=0.0d0
               frrz(i)=0.0d0
  995       continue

         end if

cc loop over pre-defined pairs only.
         do 996 k=1,npairs

            i=npairident(k,1)
            j=npairident(k,2)

            ipoint=ntype(i)
            jpoint=ntype(j)

            drijsq=(dxsav(i,j)**2.0d0)
     x          +(dysav(i,j)**2.0d0)
     x          +(dzsav(i,j)**2.0d0)
            drij=drijsq**(0.5d0)

            u=xkharm(ipoint,jpoint)
     x       *((drij-r0harm(ipoint,jpoint))**2.0d0)

            engsr=engsr+u

            dudr=2.0d0*xkharm(ipoint,jpoint)
     x       *(drij-r0harm(ipoint,jpoint))
            frrxi=dudr*dxsav(i,j)/drij
            frryi=dudr*dysav(i,j)/drij
            frrzi=dudr*dzsav(i,j)/drij

            frrx(i)=frrx(i)-frrxi
            frry(i)=frry(i)-frryi
            frrz(i)=frrz(i)-frrzi
            frrx(j)=frrx(j)+frrxi
            frry(j)=frry(j)+frryi
            frrz(j)=frrz(j)+frrzi

  996 continue

      endif
