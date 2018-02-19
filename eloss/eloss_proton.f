       subroutine eloss_h(path,betah,eloss_tot,zproj,option,hw)
 
C------------------------------------------------------------------------------
c subroutine: eloss_h
c     author: L.Todor
c       date: July 1999
c    purpose: Calculate the energy loss for hadrons
c             in a layer of thickness thick; first we
c             calculate an average energy loss; 
c             according to thick layer definition we 
c             then try to simulate a quantity Gaussian, Landau or
c             Vavilov distributed with that average 
c             
c  reference: William R. Leo "Techniques for Nuclear Particle"
c             C. Grupen "Particle Detectors", Cambridge University Press,1996
c       used: apriori we suppose that particle is proton; this limitation
c             can be easily removed later
c             option allows to choose if average energy loss is used
c             or will be decided in the routine what distribution will
c             be used on top of the average energy loss and consequently
c             set the option
c             option=0 -> only average eloss calculation
c             option=1(else) -> it will be decided in the routine what 
c             distribution to be used; as return code than
c             path=layer thickness in m
c             rho=target density in g/cm^3
c             imep=mean excitation potential of the target
c
c   Modifications:
c
c     P.E. Ulmer - Aug. 8, 2005
c     Add calculation of most probable energy loss, fixed bug and many
c     stylistic changes.
c
C------------------------------------------------------------------------------

       implicit none

       integer zproj,option,hw
 
       double precision thick,thick_equiv,betah,zp2,path,fin
       double precision cfixed,dedx,beta2,lncalc,gamma,eta
       double precision wmax,mec2,k,pi,euler,p1,p2,p3,p4,p5,p6
       double precision elossv,sigma,fi,xmean
       double precision emax,xi,f1,f2,ratio,mproj,landmax
       double precision shift,eloss_mp_tmp
       double precision eloss_window,eloss_tgt,eloss_mean,eloss_tot

       real rndm(2),lambda

       parameter (pi=3.14159265358979324D0)
       parameter (euler=0.577215D0)
       parameter (mec2=0.51099906D0)
       parameter (p1=.60715D0,p2=.67794D0,p3=.052382D0,p4=.94753D0,
     +            p5=.74442D0,p6=1.1934D0)

c------------------------------------------------------------------------------
c  "cfixed" is in MeV/g cm^2 ; =4 pi Nav re^2 me c^2/ A (A=1g/mol);
c  from Physical Review D, Particle and Fields, July 1996
c------------------------------------------------------------------------------

       parameter (cfixed=0.307075D0)

       INCLUDE 'masses.cmn'
       INCLUDE 'eloss.cmn'

       beta2 = betah*betah
       gamma = 1.d0/sqrt(1.d0-beta2)
       eta   = betah*gamma
       mproj = EJECT_MASS
       zp2   = dble(zproj*zproj)

c------------------------------------------------------------------------------
c Maximum kinetic energy which can be imparted to a free electron in
c a single collision, assuming Mproton>>me
c------------------------------------------------------------------------------

       wmax   = 2.D0*mec2*eta*eta    

       lncalc = log(wmax/imep)
       fin    = 2.D0*lncalc-2.D0*beta2-delta(3)-2.D0*cshell(3)/dble(ztg)
       dedx   = ( (cfixed/2.D0)*(zp2*dble(ztg)/atg)/beta2 ) * fin

c------------------------------------------------------------------------------
c  Assume path in meters and rho in g/cm^3 -> thickness in g/cm^2
c------------------------------------------------------------------------------

       thick = path*100.D0*rho

c------------------------------------------------------------------------------
c     eloss_tgt    = mean energy loss (not most probable)
c     eloss_window = mean energy loss (not most probable)
c------------------------------------------------------------------------------

       eloss_tgt = thick * dedx
       if (hw.eq.4.or.hw.eq.5) then
          eloss_window = elwalk(hw)
       else
          eloss_window = 0.d0
       end if
       eloss_mean = eloss_tgt + eloss_window   ! total mean energy loss

c------------------------------------------------------------------------------
c     Now assign the values for performing the mean energy loss correction
c     (implemented in mceep.f).
c------------------------------------------------------------------------------

       celoss(3)  =  eloss_mean                  ! 3=ejectile

       if (option.eq.0) then       ! only mean value requested
          eloss_tot = eloss_mean
          return
       endif

c------------------------------------------------------------------------------
c     Compute the "equivalent thickness" used to define xi.
c
c         The equiv. thickness is the actual thickness of target
c         material traversed, scaled up to account for the windows:
c
c             thick_equiv = thickness * ( eloss_mean / eloss_tgt)
c
c     The equivalent thickness (i.e. xi) will be used to calculate 
c     both the most probable energy loss and the straggling.
c     So these will both effectively include the windows.
c------------------------------------------------------------------------------

       if (eloss_tgt .ne. 0.) then
          thick_equiv = thick * ( eloss_mean / eloss_tgt )
          xi = (dedx/fin) * thick_equiv
       endif

       if (eloss_tgt .eq. 0) then           ! this should be VERY rare
          if (eloss_window .ne. 0.) then
             xi = eloss_window/fin
          else
             eloss_tot   = 0.d0
             eloss_mp(3) = 0.d0             ! 3=ejectile
             return
          endif
       endif

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
c         log_epsilon = log( (1.-beta2)*imep*imep/(2.*mec2*beta2) ) + beta2
c         epsilon  = exp(log_epsilon)
c         eloss_mp_tmp = xi * ( log(xi/epsilon)+0.198-delta(3) )
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

       ratio = mec2/mproj
       f1    = 2.D0*mec2*eta*eta
       f2    = 1.D0+2.D0*ratio*gamma+ratio*ratio
       emax  = f1/f2  
       k     = xi/emax

       if (k.le.0.01D0) then
c------------------------------------------------------------------------------
c if k<=0.01 Landau distribution
c------------------------------------------------------------------------------
          xmean = -beta2-log(k)+euler-1.d0
          landmax = p1+p6*xmean+(p2+p3*xmean)*exp(p4*xmean+p5)  
 25       call glandg(lambda)    
          if (dble(lambda).gt.landmax) go to 25
          shift  = xi*(beta2+log(k)-euler+1.D0)   ! negative quantity
          elossv = xi*dble(lambda)  + shift       ! now mean value = 0
          eloss_mp_tmp = eloss_mean + shift       ! most probable value

       else if (k.gt.0.01D0.and.k.le.10.D0) then
c------------------------------------------------------------------------------
c if 0.01<k<=10.0 Vavilov distribution
c------------------------------------------------------------------------------
          call RANECU(rndm,1)
          call gvaviv(lambda,real(k),real(beta2),rndm(1))
          shift  = xi*(beta2+log(k)-euler+1.D0)   ! negative quantity
          elossv = xi*dble(lambda)  + shift       ! now mean value = 0
          eloss_mp_tmp = eloss_mean + shift       ! most probable value

       else if (k.gt.10.D0) then
c------------------------------------------------------------------------------
c if k>10.0 Gaussian distribution
c------------------------------------------------------------------------------
          sigma = eloss_mean*wmax*(1.D0-beta2/2.D0)
          sigma = sqrt(sigma)
 30       call RANECU(rndm,2)
          if (dble(rndm(1)).le.0.D0) go to 30
          fi = -2.D0*log(dble(rndm(1)))
cpeu      elossv = sigma*sqrt(f1)*cos(2.D0*pi*dble(rndm(2)))  ! bug
          elossv = sigma*sqrt(fi)*cos(2.D0*pi*dble(rndm(2)))
          eloss_mp_tmp = eloss_mean            ! gaussian is symmetric
       end if 

c------------------------------------------------------------------------------
c  Total energy loss = (mean) + (straggling about the mean)
c------------------------------------------------------------------------------

       eloss_tot = eloss_mean + elossv
       if (eloss_tot.le.0.001D0) then  ! energy loss can't be negative
            eloss_tot = 0.D0
       endif

c------------------------------------------------------------------------------
c     Now assign the values for performing the most probable energy loss 
c     correction (implemented in mceep.f).
c------------------------------------------------------------------------------

       eloss_mp(3) = eloss_mp_tmp    ! 3=ejectile

       return 
       end

