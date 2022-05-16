      subroutine kset

c************************************************************
c Reads in ion positions,velocities and dipolar values and
c 'velocities' and thermostating information
c from a file (filename in the filenames.inpt file) for restart.
c************************************************************

      implicit double precision(a-h,o-z)

      include 'common.inc'

      do i=1,3
         do j=1,3
            fullhit(i,j)=fullhi(j,i)
         enddo
      enddo

      call dcell(fullhit,bh)

      rkmin=(twopi*dmin1(dabs(bh(7)),dabs(bh(8))
     x                  ,dabs(bh(9))))
      rksqmax=-log(conv)*4.0d0*etasq/(rkmin*rkmin)
      kmax=int(dsqrt(rksqmax))

      ksqmax=kmax*kmax

      kmaxx=int((rkmin/(twopi*bh(7)))*float(kmax))+1
      kmaxy=int((rkmin/(twopi*bh(8)))*float(kmax))+1
      kmaxz=int((rkmin/(twopi*bh(9)))*float(kmax))+1

     


      write(48,*)'(in kset )bh7,8,9',bh(7),bh(8),bh(9)
      write(48,*)'kmax ',kmax
      write(48,*)'kmaxx,y,z ',kmaxx,kmaxy,kmaxz
      write(48,*)twopi*bh(7),rkmin,twopi*bh(7)/rkmin

      
      if ((kmaxx.gt.nktot).or.(kmaxy.gt.nktot)
     x          .or.(kmaxz.gt.nktot)) then
         write(6,*)'kmaxx(y,z) > nktot'
      endif

      return

      end
