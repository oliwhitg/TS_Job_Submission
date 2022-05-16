      subroutine rdf_lj

c************************************************************
c
c  Radial distribution function code.
c
c************************************************************

      implicit double precision(a-h,o-z)
 
      include 'common.inc'
      include 'lj.inc'

c  LJ-ion rdf (into rdfpart(4+5))
      do 20 i=1,num
 
         ipoint=ntype(i)
 
         xi=x(i)
         yi=y(i)
         zi=z(i)
 
         do 10 j=1,numlj
 
            jpoint=3

            dx=xlj(j)-xi
            dy=ylj(j)-yi
            dz=zlj(j)-zi

            dist=dx*dx+dy*dy+dz*dz
            dist=dsqrt(dist)
 
            nbinrdf=int((dist)/drrdf)

            if(nbinrdf.ge.0.and.nbinrdf.le.400) then
                rdftot(nbinrdf)=rdftot(nbinrdf)+2.0d0
                rdfpart(nbinrdf,ipoint,jpoint)
     x               =rdfpart(nbinrdf,ipoint,jpoint)+2.0d0
            endif
 
   10    continue

   20 continue

      do 40 i=1,numlj
 
         xi=xlj(i)
         yi=ylj(i)
         zi=zlj(i)
 
         do 30 j=1,numlj
 
            if (i.eq.j) goto 30
            dx=xlj(j)-xi
            dy=ylj(j)-yi
            dz=zlj(j)-zi

            dist=dx*dx+dy*dy+dz*dz
            dist=dsqrt(dist)
 
            nbinrdf=int((dist)/drrdf)

            if(nbinrdf.ge.0.and.nbinrdf.le.400) then
                rdftot(nbinrdf)=rdftot(nbinrdf)+2.0d0
                rdfpart(nbinrdf,3,3)
     x               =rdfpart(nbinrdf,3,3)+2.0d0
            endif
 
   30    continue

   40 continue

      return
 
      end
