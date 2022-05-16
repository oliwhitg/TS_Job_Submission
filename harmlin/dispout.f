      subroutine dispout

c*************************************************************
c Writes out mean-squared displacement data for later use.
c*************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'
         running_tot=0 
      do 10 i=1,num

         write(98)xdisp(i),ydisp(i),zdisp(i)
c         write(*,*)xdisp(i)*frrx(i),ydisp(i)*frry(i)
c         write(*,*)"fdr",(xdisp(i)*frrx(i))+(ydisp(i)*frry(i))
         running_tot=running_tot+(xdisp(i)*frrx(i))+(ydisp(i)*frry(i))
         xdisp(i)=0.0d0
         ydisp(i)=0.0d0
         zdisp(i)=0.0d0

   10 continue
c         write(*,*) "tot_fdr", running_tot

      return

      end
