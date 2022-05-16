      subroutine rdf

c************************************************************
c  Radial distribution function code.
c************************************************************

      implicit double precision(a-h,o-z)
 
      include 'common.inc'

      nrdfcall=nrdfcall+1

c      dr=halfbox/100.0d0
      dr=0.04
c      dr=(cellvol**(1.0d0/3.0d0))/400.0d0

      do 10 i=1,num-1
 
         ipoint=ntype(i)
 
 
         do 20 j=i+1,num
 
            jpoint=ntype(j)
            nparttype=ipoint+jpoint-1
c
c The arrays x(), y() and z() now contain the system coordinates
c in the cell frame (hence the cf in the variable name).
c
cc            dist=drsav(i,j)
 
            nbinrdf=int((dist)/dr)

            if(nbinrdf.gt.0.and.nbinrdf.le.1000) then
                rdftot(nbinrdf)=rdftot(nbinrdf)+2.0d0

c	changes for ternary system

                rdfpart(nbinrdf,ipoint,jpoint)
     x   =rdfpart(nbinrdf,ipoint,jpoint)+2.0d0
            endif


c                rdfpart(nbinrdf,nparttype)=rdfpart(nbinrdf,nparttype)+2.0d0
 
   20    continue

   10 continue

      return
 
      end
