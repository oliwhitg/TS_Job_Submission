      subroutine debye_scherer_out

c**********************************************************
c Outputs the structure factor calculated in recipE.
c**********************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'

      double precision xijrec(nspmax,nspmax)

c     local characters
      CHARACTER*9 indi
      CHARACTER*6 skname
      CHARACTER*60 filename

      indi = '123456789'
      skname = 'sk.out'

c	change for ternary systems

      do 9 i=1,nspec
         do 8 j=i,nspec
            xijrec(i,j)=1.0d0/dsqrt(dble((nsp(i))*(nsp(j))))
8        continue
9     continue

      runrec=1.0d0/dble(nskcall)
      fnbin=float(nbin)

      lk=1
      do 11 ii=1,nspec
         do 12 jj=ii,nspec

            filename=skname//indi(ii:ii)//indi(jj:jj)

            open(65,file=filename,status='new')

               do 100 ibin=1,nbin
                  rk=sqrt(rksqmax_ds*(float(ibin)+0.5d0)/fnbin)
                  binwidth=sqrt(float(ibin+1)*rksqmax_ds) -
     x            sqrt(float(ibin)*rksqmax_ds)
                  binwidth=rk**2*binwidth

                  count=float(norm_ds(ibin))
c.......NB the 5 here is arbitrary  -- to let in only those
c       bins which are appreciably sampled
                  icount =int(5.0d0*count *runrec)

                  if(icount.ne.0) then

                     skk=sk_ds(lk,ibin)*xijrec(ii,jj)
     x                                 /count
cc     x                                 *runrec/count

	             write(65,*)rk,skk 
                  endif
  100          continue

            close(65)

            lk=lk+1
   12    continue
   11 continue

            open(65,file='sk_norm.out',status='new')

               do 110 ibin=1,nbin
                  rk=sqrt(rksqmax_ds*(float(ibin)+0.5d0)/fnbin)
                  binwidth=sqrt(float(ibin+1)*rksqmax_ds) -
     x            sqrt(float(ibin)*rksqmax_ds)
                  binwidth=rk**2*binwidth

                  count=float(norm_ds(ibin))
c.......NB the 5 here is arbitrary  -- to let in only those
c       bins which are appreciably sampled
                  icount =int(5.0d0*count *runrec)

                  if(icount.ne.0) then
                     write(65,*)rk,binwidth,count
                  endif
  110          continue

            close(65)

      return

      end
