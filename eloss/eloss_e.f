       subroutine eloss_e(path,e0,n,eloss_tot,option,ew)
 
C------------------------------------------------------------------------------
c subroutine: eloss_e
c     author: L.Todor
c       date: July 1999
c    purpose: Calculate the electron energy loss by collision
c             while going through a path of thickness t
c  reference: Physical Review D, Particle and Fields, 1 Jul 1996
c             William R. Leo "Techniques for Nuclear Particle" (basic formula)
c             C. Grupen "Particle Detectors", Cambridge University Press,1996
c   comment:  In the energy range of the electron for the Hall A experiments
c             (momentum >350MeV/c) the external Bremsstrahlung dominates
c             over the electron energy loss; energy loss manifests a kind 
c             of saturation
c             same dedx regardless the momentum value.
c             option=0 -> only average eloss calculation
c             option=1(else) -> it will be decided in the routine what      
c             distribution to be used;             
c
c   Modifications:
c
c     P.E. Ulmer - Aug. 8, 2005
c     Add calculation of most probable energy loss, fixed bug and many
c     stylistic changes.
c
C------------------------------------------------------------------------------

       implicit none

       integer n,option,ew 

       double precision argu,emax
       double precision path,e0,known,tau,ftau,mec2,fin,dedx
       double precision betae,tau2,be2,imepmec2,thick,thick_equiv
       double precision euler,p1,p2,p3,p4,p5,p6,xmean,landmax,k
       double precision elossv,sigma,fi,pi,xi
       double precision shift,eloss_mp_tmp
       double precision eloss_window,eloss_tgt,eloss_mean,eloss_tot

       real rndm(2),lambda

       parameter (euler=0.577215D0)
       parameter (mec2=0.51099906D0)
       parameter (p1=.60715D0,p2=.67794D0,p3=.052382D0,p4=.94753D0,
     +            p5=.74442D0,p6=1.1934D0)
       parameter (pi=3.14159265358979324D0)

c------------------------------------------------------------------------------
c  "known" is in MeV/g cm^2 ; =4 pi Nav re^2 me c^2/ A (A=1g/mol);
c  from Physical Review D, Particle and Fields, July 1996
c------------------------------------------------------------------------------

       parameter (known=0.307075D0)

       INCLUDE 'eloss.cmn'
      
       tau      = e0/mec2-1.d0    ! kinetic energy in units of mec^2
       tau2     = tau*tau
       betae    = (sqrt(e0*e0-mec2*mec2))/e0        
       be2      = betae*betae 
       imepmec2 = imep*imep/(mec2*mec2) 

       ftau = 1.d0-be2+(tau2/8.d0-(2.d0*tau+1.d0)*log(2.D0))/
     &               (tau2+2.d0*tau+1.d0)
       argu = (tau2/2.D0)*(tau+2.d0)/imepmec2
       fin  = ftau+log(argu)-delta(n)-2.D0*cshell(n)/dble(ztg)
       dedx = ( (known/2.D0)*(dble(ztg)/atg)/be2 ) * fin

c------------------------------------------------------------------------------
c     Assume path in meters and rho in g/cm^3 -> thickness in g/cm^2
c------------------------------------------------------------------------------

       thick = path*100.D0*rho

c------------------------------------------------------------------------------
c     eloss_tgt    = mean energy loss (not most probable)
c     eloss_window = mean energy loss (not most probable)
c------------------------------------------------------------------------------

       eloss_tgt = thick * dedx
       if (ew.ne.0.and.ew.le.3) then
          eloss_window = elwalk(ew)
       else
          eloss_window = 0.d0
       endif
       eloss_mean = eloss_tgt + eloss_window   ! total mean energy loss

c------------------------------------------------------------------------------
c     Now assign the values for performing the mean energy loss correction
c     (implemented in mceep.f).
c------------------------------------------------------------------------------

       if(ew.eq.1) then
          celoss(1)   = eloss_mean               ! beam
       else
          celoss(2)   = eloss_mean               ! scattered electron
       end if

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
             eloss_tot = 0.d0
             if(ew.eq.1) then
                eloss_mp(1) = 0.d0          ! beam
             else
                eloss_mp(2) = 0.d0          ! scattered electron
             end if
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

       emax  = e0    ! maximum energy transfer to atomic electron
       k     = xi/emax

       if (k.le.0.01D0) then
c------------------------------------------------------------------------------
c if k<=0.01 Landau distribution
c------------------------------------------------------------------------------
          xmean = -be2-log(k)+euler-1.D0
          landmax = p1+p6*xmean+(p2+p3*xmean)*exp(p4*xmean+p5)  
 25       call glandg(lambda)    
          if (dble(lambda).gt.landmax) go to 25
          shift  = xi*(be2+log(k)-euler+1.D0)   ! negative quantity
          elossv = xi*dble(lambda)  + shift     ! now mean value = 0
          eloss_mp_tmp = eloss_mean + shift     ! most probable value

       else if (k.gt.0.01D0.and.k.le.10.D0) then
c------------------------------------------------------------------------------
c if 0.01<k<=10.0 Vavilov distribution
c------------------------------------------------------------------------------
          call RANECU(rndm,1)
          call gvaviv(lambda,real(k),real(be2),rndm(1))
          shift  = xi*(be2+log(k)-euler+1.D0)   ! negative quantity
          elossv = xi*dble(lambda)  + shift     ! now mean value = 0
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
cpeu      elossv = sigma*sqrt(f1)*cos(2.D0*pi*dble(rndm(2)))  ! bug
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

c------------------------------------------------------------------------------
c     Now assign the values for performing the most probable energy loss 
c     correction (implemented in mceep.f).
c------------------------------------------------------------------------------

       if(ew.eq.1) then
          eloss_mp(1) = eloss_mp_tmp               ! beam
       else
          eloss_mp(2) = eloss_mp_tmp               ! scattered electron
       end if

       return 
       end
