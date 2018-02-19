      real*8 function straggling(k, xi, beta, eloss_mean, emax, dlambda)

       implicit none

       integer n,option,ew 

       real*8 argu,emax
       real*8 e0,fin
       real*8 beta, be2
       real*8 euler,p1,p2,p3,p4,p5,p6,xmean,landmax,k
       real*8 elossv,sigma,fi,pi,xi
       real*8 shift,eloss_mp_tmp
       real*8 eloss_mean,eloss_tot
       real*8 dlambda

       real rndm(2),lambda

       parameter (euler=0.577215D0)
       parameter (p1=.60715D0,p2=.67794D0,p3=.052382D0,p4=.94753D0,
     +            p5=.74442D0,p6=1.1934D0)
       parameter (pi=3.14159265358979324D0)

       be2 = beta*beta
c------------------------------------------------------------------------------
c  "known" is in MeV/g cm^2 ; =4 pi Nav re^2 me c^2/ A (A=1g/mol);
c  from Physical Review D, Particle and Fields, July 1996
c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
c
c     Compute STRAGGLING and MOST PROBABLE ENERGY LOSS
c
c     Both of these will depend on the type of distribution.
c
c     Note that as an expedient one might use the formula found in
c       - W.R. Leo, Techniques for Nuclear and
c         Particle Physics Experiments, 2nd rev. Ed., Springer-Verlag.
c         (see Eq. 2.99):
c
c         log_epsilon = log( (1.-be2)*imep*imep/(2.*mec2*be2) ) + be2
c         epsilon  = exp(log_epsilon)
c         eloss_mp_tmp = xi * ( log(xi/epsilon)+0.198-delta(n) )
c
c     However, to be consistent with the straggling, I do it case by
c     case instead.
c
c     Note that the straggling is about the mean energy loss,
c     not the most probable energy loss.  Therefore, for instance in the
c     Landau case, the values lambda returned by glandg have a peak
c     very close to zero.  The quantity elossv is this distribution
c     plus  xi*(+be2+log(k)-euler+1.D0), the latter being a negative
c     quantity.  Thus this "addition" moves the landau distribution leftward
c     so that it's mean value is now zero instead of its most probable
c     value being zero.  To get the total energy loss, one just adds
c     this shifted distribution to the mean energy loss computed above.
c     The most probable energy loss is just the difference between the
c     mean value and the magnitude of the shift.  Simple!
c
c------------------------------------------------------------------------------

       if (k.le.0.01D0) then
c------------------------------------------------------------------------------
c if k<=0.01 Landau distribution
c------------------------------------------------------------------------------
          xmean = -be2-log(k)+euler-1.D0
          landmax = p1+p6*xmean+(p2+p3*xmean)*exp(p4*xmean+p5)  
 25       call glandg(lambda)    
          dlambda = dble(lambda)
          if (dlambda.gt.landmax) go to 25
          shift  = xi*(be2+log(k)-euler+1.D0)   ! negative quantity
          elossv = xi*dlambda  + shift     ! now mean value = 0
          eloss_mp_tmp = eloss_mean + shift     ! most probable value

       else if (k.gt.0.01D0.and.k.le.10.D0) then
c------------------------------------------------------------------------------
c if 0.01<k<=10.0 Vavilov distribution
c------------------------------------------------------------------------------
          call RANECU(rndm,1)
          call gvaviv(lambda,real(k),real(be2),rndm(1))
          dlambda = dble(lambda)
          shift  = xi*(be2+log(k)-euler+1.D0)   ! negative quantity
          elossv = xi*dlambda  + shift     ! now mean value = 0
          eloss_mp_tmp = eloss_mean + shift     ! most probable value

       else if (k.gt.10.D0) then
c------------------------------------------------------------------------------
c if k>10.0 Gaussian distribution
c------------------------------------------------------------------------------
          sigma = eloss_mean*emax*(1.D0-be2/2.D0)
          sigma = sqrt(sigma)
 30       call RANECU(rndm,2)
          if (dble(rndm(1)).le.0.D0) go to 30
          fi = -2.D0*log(dble(rndm(1)))
          elossv = sigma*sqrt(fi)*cos(2.D0*pi*dble(rndm(2)))
          eloss_mp_tmp = eloss_mean            ! gaussian is symmetric
       end if 

c------------------------------------------------------------------------------
c  Total energy loss = (mean) + (straggling about the mean)
c------------------------------------------------------------------------------

       eloss_tot = eloss_mean + elossv
       if (eloss_tot.le.0.001D0) then   ! energy loss can't be negative
            eloss_tot = 0.D0
       end if
       straggling = elossv
       return 
       end
