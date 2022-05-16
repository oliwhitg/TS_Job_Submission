      program msd_recover

      implicit double precision(a-h,o-z)

      parameter(nummax=500)

      double precision xold(nummax)
     x                ,yold(nummax)
     x                ,zold(nummax)

      character*80 filein,fileout

      open(10,file='msd_recover3.inpt',status='old')

         read(10,*)num
         read(10,*)nanion
         read(10,*)ncation
         read(10,*)nconfigs
         read(10,*)boxlen
         read(10,'(a)')filein
         read(10,'(a)')fileout

      close(10)

      halfbox=0.5d0*boxlen

      none=1
      ntwo=2

      open(10,file=filein,status='old')
      open(99,file=fileout,form='unformatted')

c Read in intitial positions of each ion

      do 10 i=1,num

          read(10)xold(i),yold(i),zold(i)

10    continue

c Calculate displacement during each step

      do 20 i=1,nconfigs-1
          do 30 j=1,num

              read(10)x,y,z
 
              dx=x-xold(j)
              dy=y-yold(j)
              dz=z-zold(j)

              dx=dx-boxlen*int(dx/halfbox)
              dy=dy-boxlen*int(dy/halfbox)
              dz=dz-boxlen*int(dz/halfbox)

              xold(j)=x
              yold(j)=y
              zold(j)=z

              write(99)dx,dy,dz

30        continue

20    continue

      do 40 i=1,num
          write(99)dx,dy,dz
40    continue

      close(99)
      close(10)

      end
