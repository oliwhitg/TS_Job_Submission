      function erfunc(x)
c
c Calculate the error function:
c
      implicit double precision (a-h,o-z)

      erfunc=1.0d0/((1.0d0+x*(0.0705230745d0+x*(0.0422820123d0
     x +x*(0.0092705272d0+x*(0.0001520143d0+x*(0.0002765672d0+x*
     x 0.0000430638d0))))))**16.0d0)

      return
      end
