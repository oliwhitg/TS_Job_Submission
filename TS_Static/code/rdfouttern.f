      subroutine rdfout
 
c************************************************************
c Calculates partial and total rdf's and outputs them.
c************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'

      dimension xion(6)

c     local characters
      CHARACTER*9 indi
      CHARACTER*7 rdfname
      CHARACTER*60 filename

      indi = '123456789'
      rdfname = 'rdf.out'

      const=dble(num)*fourpi*cellvol3rec
      const=const*dble(nrdfcall*num)
c      drrdf=halfboxx*0.01d0

      rlmax=(cellvol**(1.0d0/3.0d0))/2.0d0
      drrdf=0.04
      nbins =int(rlmax/drrdf)

      if (nbins.gt.1000) nbins=1000

      xtot = 0.0d0
      do 20 i=1, nspec
      xtot = float(nsp(i))/float(nionunit) + xtot
      xion(i) = float(nsp(i))/float(nionunit)
20    continue 

      open(64,file='rdftot.out',status='new')
    
       
         do 10 i=1,nbins
 
            rlower=dble(i)*drrdf
            rupper=rlower+drrdf
            xnorm=const*((rupper**3.0d0)-(rlower**3.0d0))
            xnormrec=1.0d0/xnorm
 
            rdftot(i)=rdftot(i)*xnormrec

            do 11 ii=1,nspec      
               do 12 jj=ii,nspec      
                  rdfpart(i,ii,jj)=
     x                  0.5d0*xtot*xtot/((xion(ii))*(xion(jj)))
     x                   *rdfpart(i,ii,jj)/xnorm
                  if (ii.eq.jj) then
                     rdfpart(i,ii,jj)=rdfpart(i,ii,jj)*2.0d0
                  endif

12             continue
11          continue     

         write(64,*) rlower,rdftot(i)

   10    continue

      close(64)

      do 100 i=1,nspec

         do 110 j=i,nspec

            filename=rdfname//indi(i:i)//indi(j:j)

            open(65,file=filename,status='new')

               do 120 k=1,nbins

                  rlower=dble(k)*drrdf
                  rlower=rlower+(0.5d0*drrdf)

                  write(65,*) rlower,rdfpart(k,i,j)

  120          continue

            close(65)

  110    continue

  100 continue

      write (6,*)
      write (6,*) '**** Radial distribution function written out ****'
      write (6,*)

      return

      end

