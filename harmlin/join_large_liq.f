         Program Join 


         Implicit None
         
         Integer, Parameter::K6=selected_real_kind(12, 10)
         Real (k6)   delta, a_length, b_length, c_length, c
         Real (k6), dimension(1:1000):: x, y, z
         INTEGER a,b, no_atoms
         Character d

         Open (4, file= "join.inpt", status = "old")
         read(4,*) no_atoms
         read(4,*) delta
         Close (4)

         Open (1, file= "liq_a", status = "old")
         b=0

         Do a=1,5
         Read(1,*)d 
         End do 

         Do a = 1, no_atoms
         b=b+1
         Read (1,*) x(b),y(b),z(b)
         End Do
         
         Do a = 1, no_atoms*3 !1x velocities, 1x dipole positions, 1 x dipole vel
         Read (1,*) c 
         End Do

         !Do a=1,32 !Read past 2xsim length 30 x barostat data 
         !Read(1,*) c
         !End Do
         
         Do a=1,6  !Read past cell matrix
         Read(1,*) c
         End do         
         
         Read(1,*) a_length
         write(*,*) a_length
         Read(1,*) b_length
         write(*,*) b_length
         Read(1,*) c_length
         write(*,*) c_length
         
         Close (1)
          
        
         Open (2, file= "liquid_nano.crds", status = "new")
         
         Write(2,*) "T" 
         Do a=1,4
          Write(2,*) "F"
         End Do

         Do a=1,no_atoms/2
         Write(2,*) x(a),y(a),z(a)
         End do

         Do a=1,no_atoms/2
         Write(2,*) x(a),y(a)+b_length+delta,z(a)
         End do

         Do a=1,no_atoms/2
         Write(2,*) x(a)+a_length+delta,y(a),z(a)
         End do

         Do a=1,no_atoms/2
         Write(2,*) x(a)+a_length+delta,y(a)+b_length+delta,z(a)
         End do

         Do a=1+(no_atoms/2), no_atoms
         Write(2,*) x(a),y(a),z(a)
         End do

         Do a=1+(no_atoms/2), no_atoms
         Write(2,*) x(a),y(a)+b_length+delta,z(a)
         End do

         Do a=1+(no_atoms/2), no_atoms
         Write(2,*) x(a)+a_length+delta,y(a),z(a)
         End do

         Do a=1+(no_atoms/2), no_atoms
         Write(2,*) x(a)+a_length+delta,y(a)+b_length+delta,z(a)
         End do




         write(2,*) "1. 0. 0."
         write(2,*) "0. 1. 0."
         write(2,*) "0. 0. 1."

           write (2,*) (a_length*2)+delta
           write (2,*) (b_length*2)+delta
           write (2,*) c_length

         Close (2) 

        End Program Join
