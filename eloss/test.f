      implicit none
      real*8 length, dens, aeff, zeff, epart, mpart, Eloss
      real*8 deloss, eloss_mean, eloss_mp, lambda, k, shift, xi
      integer typeflag, nevents
      integer i

      common /results/eloss_mp, eloss_mean, deloss, 
     >     lambda, k, shift, xi

      namelist /output/eloss_mp, eloss_mean, deloss , typeflag, 
     >     Eloss, lambda, k, shift, xi

      print *, 'enter length, dens, aeff, zeff, epart, mpart, typflag :'
      read (*,*) length,
     >      dens,
     >      aeff,
     >      zeff,
     >      epart,
     >      mpart,
     >      typeflag
      print *, 'number of events : '
      read (*,*) nevents

c init random numbers      
      call recuin( 287364, 782394)
      open (16,file = 'test.data')

      write (16,*) '# lenght = ', length
      write (16,*) '# density = ', dens
      write (16,*) '# aeff = ', aeff
      write (16,*) '# zeff = ', zeff
      write (16,*) '# epart = ', epart
      write (16,*) '# mpart = ', mpart
      write (16,*) '# typeflag = ', typeflag
      write (16,*) '# '
      write (16,*)'#! demp[f,0]/ demean[f,1]/ de[f,2]/ type[i,3]/ '//
     >     'eloss[f,4]/ lam[f,5]/ k[f,6]/ shift[f,7]/ xi[f,8]/ ' 

      do i = 1, nevents
         call enerloss_new(length,dens,zeff,aeff,epart,mpart,
     >        typeflag, Eloss)
         write (16,*)eloss_mp, eloss_mean, deloss , typeflag, 
     >     Eloss, lambda, k, shift, xi
      enddo
      stop
      end
