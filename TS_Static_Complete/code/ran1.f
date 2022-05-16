      function ran1(idum)

c************************************************************
c   From `Numerical Recipes'.
c   Produces a uniform random number in the 0-1 range.
c   Called with idum = -1 to seed or re-seed the sequence.
c************************************************************

      implicit double precision (a-h,o-z)

      parameter (m1=259200,ia1=7141,ic1=54773,rm1=1.0d0/m1)
      parameter (m2=134456,ia2=8121,ic2=28411,rm2=1.0d0/m2)
      parameter (m3=243000,ia3=4561,ic3=51349)

      save r,ix1,ix2,ix3

      double precision r(100)

c if statement executed on initializing the sequence.
   
      if (idum.lt.0.) then

         ix1=mod(ic1-idum,m1)
         ix1=mod(ia1*ix1+ic1,m1)
         ix2=mod(ix1,m2)
         ix1=mod(ia1*ix1+ic1,m1)
         ix3=mod(ix1,m3)

         do 10 j=1,97
            ix1=mod(ia1*ix1+ic1,m1)
            ix2=mod(ia2*ix2+ic2,m2)
            r(j)=(dble(ix1)+dble(ix2)*rm2)*rm1
   10    continue

      endif
 
      ix1=mod(ia1*ix1+ic1,m1)
      ix2=mod(ia2*ix2+ic2,m2)
      ix3=mod(ia3*ix3+ic3,m3)

      j=1+(97*ix3)/m3
      ran1=r(j)
      r(j)=(dble(ix1)+dble(ix2)*rm2)*rm1

      return

      end
