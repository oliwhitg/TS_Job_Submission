      program aaa

      open(10,file='crys.crds',status='old')

         do 10 i=1,4800
            read(10,*)x,y,z
            write(99,*)x/.529177,y/.529177,z/.529177
   10    continue

      close(10)

      end
