c for testing
      real*8 function get_eloss(length,dens,zeff,aeff,epart,mpart,
     >     typeflag)
      real*8 length, dens, aeff, zeff, epart, mpart, Eloss
      integer typeflag

      call enerloss_new(length,dens,zeff,aeff,epart,mpart,
     >     typeflag, Eloss)
      get_eloss = Eloss 
      return 
      end

c the real routine
      subroutine enerloss_new(length,dens,zeff,aeff,epart,mpart,
     >     typeflag, Eloss)

      implicit none

      real*8 length, dens, aeff, zeff
      integer typeflag
      real*8 pi, euler, mec2, coef
      real*8 p1, p2, p3, p4, p5, p6
      parameter (pi=3.14159265358979324D0)
      parameter (euler=0.577215D0)
      parameter (mec2=0.51099906D0)
      parameter (p1=.60715D0,p2=.67794D0,p3=.052382D0,p4=.94753D0,
     +     p5=.74442D0,p6=1.1934D0)

      parameter (coef = 0.1535374581D0) ! 2*pi*N_a*(r_e)^2*m_e*c^2
      integer zpart
      real*8 ppart, epart, mpart, thick
      real*8 beta, gamma, beta2, eta
      real*8 ratio, zx2, xi, f1, f2
      real*8 emax, k, shift, lambda, sigma
      real*8 Eloss, eloss_mp, eloss_mean, deloss

      real*8 straggling, eloss_hadron, eloss_electron
      namelist /eloss_nml/length,dens,zeff,aeff,epart,mpart,
     >     typeflag,Eloss

      common /eloss_variables/thick, gamma, ppart, beta, beta2,
     > eta, ratio, zx2, f1, f2, emax

      common /eloss_results/eloss_mp, eloss_mean, deloss, lambda,
     >     k, shift, xi

      thick = length*dens
      if (thick .eq.0d0) then
         eloss = 0.d0
         return
      endif
      gamma = epart/mpart
      ppart = sqrt(epart**2 - mpart**2)
      beta = ppart/epart

c 1=normal eloss (picked from distribution)
c 2=min eloss
c 3=max eloss
c 4=most probable eloss

c check if it an electron or a heavier particle
      if (dabs(mpart - mec2).lt.0.01) then
         beta2 = beta*beta
         xi = coef*thick*zeff/(aeff*beta2)
         k = xi/epart
         eloss_mean = eloss_electron(beta,zeff, aeff,dens,thick)
      else
c assume hadron with charge 1       
         zpart = 1

c parameter for MP energyloss
         beta2 = beta*beta
         eta = beta*gamma
         ratio = mec2/mpart
         zx2 = dfloat(zpart*zpart)
         xi = coef*thick*zx2*zeff/(aeff*beta2)
         f1 = 2.D0*mec2*eta*eta
         f2 = 1.D0 + 2.D0*ratio*gamma + ratio*ratio
         emax = f1/f2
         k = xi/emax
c calculate mean energy loss
         eloss_mean = eloss_hadron(zpart,beta,zeff, aeff,dens,thick)
      endif
      eloss_mp = eloss_mean
c calculate most probable energy loss and k to get the proper deviate
      if (k.le.10.D0) then
         shift = xi*(beta2+dlog(k)-euler+1.D0)
         eloss_mp = eloss_mean + shift
      else
         shift = 0.d0
         lambda = 0.d0
         sigma = eloss_mean*emax*(1.D0-beta2/2.D0)
         sigma = sqrt(sigma)
      endif
c random deviate from mean value either landau, vavilov or gauss depending on k        
      if (typeflag .eq. 4) then
c     most probable eloss
         deloss = shift
      elseif (typeflag .eq. 3) then
c maximum eloss see straggling.f
         if (k .gt. 10) then
            deloss = 5.d0*sigma ! gaussian distribution
         else
            lambda = 15.
            deloss = xi*lambda  + shift ! now mean value = 0
         endif
      elseif (typeflag .eq. 2) then
c     minimum eloss
         if (k .gt. 10) then
            deloss = -5.d0*sigma ! gaussian distribution
         else
            lambda = -3.5
            deloss = xi*lambda  + shift ! now mean value = 0
         endif
      elseif (typeflag .eq. 1) then
c     normal eloss
         deloss =  straggling(k, xi, beta, eloss_mean, emax, lambda)
      endif
      eloss = eloss_mean + deloss
 
c for debugging     
c      write (6, nml = eloss_nml)
      return 
      end

