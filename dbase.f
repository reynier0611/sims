	subroutine dbase_read(H)

! dbase_read reads the input file and sets the run-time flags for the run.
! Note that there are four INDEPENDENT ways to run SIMC.
!
! 1. doing_eep: (e,e'p) subcases:doing_hyd_elast, doing_deuterium, doing_heavy, doing_bound
!
! 2. doing_kaon:(e,e'K) subcases:doing_hydkaon, doing_deutkaon,doing_hekaon.
!	which_kaon= 0/ 1/ 2 for Lambda/Sigam0/Sigma- quasifree.
!	which_kaon=10/11/12 for Lambda/Sigam0/Sigma- coherent production.
!	For bound hypernuclear states, set doing_hydkaon true for ALL
!	targets (i.e. treat as a heavy proton) but with targ.Mtar_struck
!	and targ.Mrec_struck set appropriatly.
!
! 3. doing_pion:(e,e'p) subcases:doing_hydpi, doing_deutpi,doing_hepi.
!	which_pion =  0( 1) gives pi+ (pi-) quasifree production.
!	which_pion = 10(11) gives pi+ (pi-) coherent production.
!	This sets doing_hydpi true for ALL targets (i.e. treat as
!	a heavy proton) but with targ.Mtar_struck and targ.Mrec_struck
!	set appropriatly.
!
! 4. doing_phsp:Generate acceptance with radiation and cross section disabled,
!	use doing_kaon or doing_pion to set hadron mass, then
!	set doing_eep,doing_kaon and doing_pion to FALSE.
!	May or may not be properly  implemented.
!
! 5. doing_delta:(e,e'p)pi subcases:doing_hyddelta,doing_deutdelta,doing_hedelta.
!	Identical to doing_pion, except for changes in mass parameters and
!	cross section model.  Initial implementation will be for proton target
!	only (generation well be general, but cross section model will be
!	a 'shortcut' version to start with.
! 6. doing_semi: H(e,e'pi)X (doing_semipi) and H(e,e'k)X (doing_semika)
! 7. doing_rho: H(e,e'rho)p

	implicit none
	include 'radc.inc'
	include 'histograms.inc'
	include 'simulate.inc'
	include 'g_dump_all_events.inc'

	real*8 dum1,dum2,dum3,dum4,dum5,dum6,dum7
	real*8 mrec_guess

	integer*4 iread, iq2, iang
	integer*4 i, j, k, ii
	integer*4 ierr, thload, thbook, read_parameters
	logical success
	character filename*80,tmpfile*80
	character dbase_file*60 !needs to be shorter than filename


	type (histograms):: H

	! for debugging
	character*80 dump1, dump2, dump3
	data dump1/'dump1.dat'/, dump2/'dump2.dat'/,
     >       dump3/'simc_parameter_dump.dat'/

! ... Register CTP vars

	call get_defaults
!	call dump_parameters(dump1)

! ... read the dbase file.

	dbase_file=' '
	extra_dbase_file=' '
	write(6,*) ' '
	write(6,*) '************************************************'
	write(6,*) '*                                              *'
	write(6,*) '*                --- SIMC ---                  *'
	write(6,*) '*                                              *'
	write(6,*) '*          Welcome to Jefferson Lab            *'
	write(6,*) '*         Last edited on: 9/27/2017           *'
	write(6,*) '*                                              *'
	write(6,*) '************************************************'
	write(6,*) 'Enter the input filename'
	write(6,*) '(assumed to be in infiles directory):'
	read(5,'(a)') dbase_file

	write(6,*) '*************************************************'

	j=index(dbase_file,'/')
	if (j.eq.0) then					!no path given
	  write(filename,'(a)') 'infiles/'//dbase_file
	else
	  write(filename,'(a)') dbase_file
	endif
	i=index(filename,'.')
	j=index(filename,'/')
	if(i.eq.0) then		!add .inp if not included in filename
	  i=index(filename,' ')
	  if(i+2.le.len(filename)) write(filename(i:),'(''.inp'')')
	endif
	write(6,'(a10,a69)')'filename=',filename
	if (i.gt.1) base=filename(j+1:i-1)

! ... load and book input file

!	ierr = thload(filename)
	ierr = read_parameters(filename)
	if (ierr.ne.0) stop 'Loading problem!  Not going to try again...wouldnt be prudent.'
!	ierr = thbook()
	if (ierr.ne.0) stop ' Booking problem!  Not going to try again...wouldnt be prudent'
!	call dump_parameters(dump1)


! ... read the secondary dbase file.

	if (extra_dbase_file.ne.' ') then	!new filename from dbase file
	  write(filename,'(a)') 'infiles/'//extra_dbase_file
	  i=index(filename,'.')
	  j=index(filename,'/')
	  if(i.eq.0) then		!add .inp if not included in filename
	    i=index(filename,' ')
	    i f(i+2.le.len(filename)) write(filename(i:),'(''.inp'')')
	  endif
	  write(6,'(a10,a69)')'filename=',filename
!	  ierr = thload(filename)
	  ierr = read_parameters(filename)
	  if (ierr.ne.0) stop ' Loading problem!  Not going to try again...wouldnt be prudent.'
!	  ierr = thbook()
	  if (ierr.ne.0) stop ' Booking problem!  Not going to try again...wouldnt be prudent'
	endif	!extra dbase input file

! all has been read in, dump the parameters to check
!	call dump_parameters(dump2)

C DJG: Ugly hack! This must come before the test on doing_pion
	if(doing_pion .and. doing_semi) then
	   doing_semipi=.true.
	   doing_semika=.false.
	   doing_pion=.false.
	endif
	if(doing_kaon .and. doing_semi) then
	   doing_semika=.true.
	   doing_semipi=.false.
	   doing_kaon=.false.
	endif
C DJG:

! ... dbase field experiment.
	if (doing_pion) then
	  Mh=Mpi
	  if (nint(targ%A).eq.1 .and. which_pion.eq.1)
     >		stop 'Pi- production from Hydrogen not allowed!'
	  if (nint(targ%A).le.2 .and. which_pion.ge.10)
     >		stop 'Coherent production from Hydrogen/Deuterium not allowed!'
	  if (nint(targ%A).eq.3 .and. which_pion.eq.11)
     >		stop 'Coherent Pi- production from 3He not allowed!'
	  if (nint(targ%A).eq.4 .and. which_pion.ge.10)
     >		stop 'Coherent production from 4He not allowed!'
	  doing_hydpi = (nint(targ%A).eq.1)
	  doing_deutpi = (nint(targ%A).eq.2)
	  doing_hepi = (nint(targ%A).ge.3)
	  doing_eep = .false.

* quasifree production is default (which_pion=0,1 for pi+/pi0)
* Treat (A+pi) final state as production from heavy proton (which_pion=10,11)
	  if (which_pion.ge.10) then
	    doing_hydpi = .true.
	    doing_deutpi = .false.
	    doing_hepi = .false.
	  endif

	else if (doing_kaon) then
	  Mh=Mk
	  if (nint(targ%A).eq.1 .and. which_kaon.eq.2)
     >		stop 'Sigma- production from Hydrogen not allowed!'
	  if (nint(targ%A).le.2 .and. which_kaon.ge.10)
     >		stop 'Coherent production from Hydrogen/Deuterium not allowed!'
	  doing_hydkaon = (nint(targ%A).eq.1)
	  doing_deutkaon = (nint(targ%A).eq.2)
	  doing_hekaon = (nint(targ%A).ge.3)
	  doing_eep = .false.

* quasifree production is default (which_kaon=0,1,2 for lambda/sigma/sigma-).
* Treat (A+pi) final state as production from heavy proton (+10 to which_kaon)
	  if (which_kaon.ge.10) then
	    doing_hydkaon = .true.
	    doing_deutkaon = .false.
	    doing_hekaon = .false.
	  endif

	else if (doing_delta) then
	  Mh=Mp
	  if (nint(targ%A).ge.2)
     >      write(6,*) 'WARNING: Delta cross section model only set up for proton target!'
	  doing_hyddelta = (nint(targ%A).eq.1)
	  doing_deutdelta = (nint(targ%A).eq.2)
	  doing_hedelta = (nint(targ%A).ge.3)
	  doing_eep = .false.

	else if (doing_semi) then
	   if(doing_semipi) then
	      Mh=Mpi
	   elseif(doing_semika) then
	      Mh=Mk
	   endif
	   doing_hydsemi = (nint(targ%A).eq.1)
	   doing_deutsemi = (nint(targ%A).eq.2)
	   if(doing_hydsemi.and.do_fermi) then
	      write(6,*) 'WARNING: Cannot do Fermi motion for Hydrogen!'
	      write(6,*) 'Do you mean to be running deuterium?'
	      do_fermi = .false.
	   endif
	   doing_eep=.false.

        else if (doing_rho) then
	   Mh = Mrho
	   doing_hydrho = (nint(targ%A).eq.1)
	   doing_deutrho = (nint(targ%A).eq.2)
	   doing_herho = (nint(targ%A).eq.3)
	   doing_eep=.false.

	else		!doing_eep if nothing else set.
	  Mh=Mp
	  doing_eep = .true.
	  doing_hyd_elast = (nint(targ%A).eq.1)
	  doing_deuterium = (nint(targ%A).eq.2)
	  doing_heavy = (nint(targ%A).ge.3)
	endif

	Mh2 = Mh*Mh
	Mh2_final = Mh2

	if (doing_phsp) then
	  rad_flag=0
	  doing_eep=.false.	!need to set one of these to determine Mh,
	  doing_pion=.false.	!but doing_phsp is independent of others.
	  doing_kaon=.false.
	  doing_delta=.false.
	  doing_rho=.false.
	endif

! ... dbase field kinematics_main

	dEbeam=Ebeam*dEbeam/100.
	spec%e%theta=abs(spec%e%theta)/degrad
	spec%e%cos_th=cos(spec%e%theta)
	spec%e%sin_th=sin(spec%e%theta)
	spec%p%theta=abs(spec%p%theta)/degrad
	spec%p%cos_th=cos(spec%p%theta)
	spec%p%sin_th=sin(spec%p%theta)

! ... Get phi for spectrometers.  HMS/HRSr are on right, phi=3*pi/2
	if (electron_arm.eq.1 .or. electron_arm.eq.3 .or. electron_arm.eq.7) then
	  spec%e%phi = 3*pi/2.
	else if (electron_arm.eq.2 .or. electron_arm.eq.4 .or.
     >		 electron_arm.eq.5 .or. electron_arm.eq.6 .or. electron_arm.eq.8) then
	  spec%e%phi = pi/2.
	else
	  stop 'I dont know what phi should be for the electron arm'
	endif

	if (hadron_arm.eq.1 .or. hadron_arm.eq.3) then
	  spec%p%phi = 3*pi/2.
	else if (hadron_arm.eq.2 .or. hadron_arm.eq.4 .or.
     >		 hadron_arm.eq.5 .or. hadron_arm.eq.6) then
	  spec%p%phi = pi/2.
	else
	  stop 'I dont know what phi should be for the hadron arm'
	endif

! ... dbase field target

	targ%N = targ%A - targ%Z
	targ%M = targ%mass_amu*amu
	targ%Mrec = targ%mrec_amu*amu

! ... improve precision for hydrogen and deuterium targets.  For deuterium,
! ... need to wait to fill targ.Mrec until we know which nucleon was struck.

	if (nint(targ%A).eq.1) then	 !proton target
	  if (abs(targ%M-Mp).gt.0.1) write(6,*) 'WARNING: forcing target mass to Mp!!!'
	  if (abs(targ%Mrec).gt.0.1) write(6,*) 'WARNING: forcing recoil mass to zero!!!'
	  targ%M = Mp
	  targ%Mrec = 0.
	else if (nint(targ%A).eq.2) then !deuterium target
	  if (abs(targ%M-Md).gt.0.1) write(6,*) 'WARNING: forcing target mass to Md!!!'
	  targ%M = Md
	else 				!heavy target
	  Mrec_guess = targ%M - Mp	!approx A-1 mass, don't know mtar_struck yet
	  if (abs(targ%Mrec-Mrec_guess).gt.100.) then !assume OK if within 100MeV.
	    write(6,*) 'targ.Mrec=',targ%Mrec,' I was expecting closer to',Mrec_guess
	    write(6,*) 'I AM REPLACING YOUR TARG.MREC WITH MY EXPECTATION!!!!!'
	    targ%Mrec = Mrec_guess
	  endif
	endif


! ... set target/recoil masses.
	if (doing_eep) then		!Strike proton, no 'recoil' particle
	  targ%Mtar_struck = Mp
	  targ%Mrec_struck = 0.0
	  sign_hadron=1.0
	else if (doing_delta) then		!Strike (and detect) proton, pion 'recoil'
	  targ%Mtar_struck = Mp
	  targ%Mrec_struck = Mpi
	  sign_hadron=1.0
	else if(doing_semi) then         ! For now, just assuming proton mass
	   targ%Mtar_struck = Mp
	   targ%Mrec_struck = Mp    !must have at LEAST a recoiling proton
	                            !Probably more...
	   if(doing_deutsemi) then
	      targ%Mtar_struck = (Mp+Mn)/2.0
	      targ%Mrec_struck = (Mp+Mn)/2.0
	   endif
	   if(doing_hplus) then
	      sign_hadron=1.0
	   else
	      sign_hadron=-1.0
	   endif

	else if(doing_rho) then
	   targ%Mtar_struck = Mp
	   targ%Mrec_struck = Mp
	   if(doing_hplus) then
	      sign_hadron=1.0
	   else
	      sign_hadron=-1.0
	   endif


! ... for normal production, Strike p (n), recoil partile is n(p).
! ... for bound final state, use targ.Mrec if it appears to be OK (same A
! ... A as target, with one n->p or p->n.

	else if (doing_pion) then
	  if (which_pion .eq. 0) then
	    targ%Mtar_struck = Mp      ! H(e,e'pi+)n or D(e,e'pi+)nn
	    targ%Mrec_struck = Mn
	    sign_hadron=1.0
	  else if (which_pion .eq. 1) then
	    targ%Mtar_struck = Mn      ! D(e,e'pi-) pp
	    targ%Mrec_struck = Mp
	    sign_hadron=-1.0
	  else if (which_pion .eq. 10) then
	    targ%Mtar_struck = targ%M  ! A(e,e'pi+)A'
	    targ%Mrec_struck = targ%Mrec
	    Mrec_guess = targ%M - Mp + Mn
	    sign_hadron=1.0
	  else if (which_pion .eq. 11) then
	    targ%Mtar_struck = targ%M  ! A(e,e'pi-)A'
	    targ%Mrec_struck = targ%Mrec
	    Mrec_guess = targ%M - Mn + Mp
	    sign_hadron=-1.0
	  else
	    stop 'Bad value for which_pion'
	  endif

	  if (which_pion .ge. 10) then	!Coherent pion production.
	    if (abs(targ%Mrec_struck-Mrec_guess).gt.100.) then
	      targ%Mrec_struck = Mrec_guess
	      write(6,*) 'For coherent pion production, targ.mrec_amu should'
	      write(6,*) 'be the mass of the final nucleus.'
	      write(6,*) 'I AM REPLACING THE MASS FROM THE INPUT FILE WITH MY GUESS:'
	      write(6,*) 'M_final = M_A + M_N(struck) - M_N(recoil) =',Mrec_guess
	    endif
	    targ%Mrec = 0.	!no 'A-1' recoil system for coherent production.
	  endif


! ... for normal production, Strike p(n), recoil is Lambda/Sigma/Sigma-.
! ... for bound final state, use targ.Mrec if it appears to be OK (same A
! ... A as target, with one n->p or p->n.

	else if (doing_kaon) then	!Strike p(n), recoil is L0/S0(S-)
	  if (which_kaon .eq. 0) then		! p(e,eK+)Lambda0
	    targ%Mtar_struck = Mp
	    targ%Mrec_struck = Mlambda
	    sign_hadron=1.0
	  else if (which_kaon .eq. 1) then	! p(e,eK+)Sigma0
	    targ%Mtar_struck = Mp
	    targ%Mrec_struck = Msigma0
	    sign_hadron=1.0
	  else if (which_kaon .eq. 2) then	! n(e,eK+)Sigma-
	    targ%Mtar_struck = Mn
	    targ%Mrec_struck = Msigma_minus
	    sign_hadron=1.0
	  else if (which_kaon .eq. 10) then	!A(e,eK+)A_Lambda0 (bound hypernucleon)
	    targ%Mtar_struck = targ%M
	    targ%Mrec_struck = targ%Mrec
	    Mrec_guess = targ%M - Mp + Mlambda
	    sign_hadron=1.0
	  else if (which_kaon .eq. 11) then	!A(e,eK+)A_Sigma0 (bound hypernucleon)
	    targ%Mtar_struck = targ%M
	    targ%Mrec_struck = targ%Mrec
	    Mrec_guess = targ%M - Mp + Msigma0
	    sign_hadron=1.0
	  else if (which_kaon .eq. 12) then	!A(e,eK+)A_Sigma- (bound hypernucleon)
	    targ%Mtar_struck = targ%M
	    targ%Mrec_struck = targ%Mrec
	    Mrec_guess = targ%M - Mn + Msigma_minus
	    sign_hadron=1.0
	  else
	    stop 'Bad value for which_kaon'
	  endif

	  if (which_kaon .ge. 10) then	!Coherent hypernucleus final state.
	    if (abs(targ%Mrec_struck-Mrec_guess).gt.100.) then
	      targ%Mrec_struck = Mrec_guess
	      write(6,*) 'For coherent kaon production, targ.mrec_amu should'
	      write(6,*) 'be the mass of the final hypernucleus.'
	      write(6,*) 'I AM REPLACING THE MASS FROM THE INPUT FILE WITH MY GUESS:'
	      write(6,*) 'M_hypernucleus = M_A - M_N + M_Y =',Mrec_guess
	    endif
	    targ%Mrec = 0.	!no 'A-1' recoil system for coherent production.
	  endif
	endif

! ... Now that we know targ.Mtar_struck, we can improve the precision of
! ... targ.Mrec for deuterium target.

	if (nint(targ%A).eq.2) then			!force Mrec.
	  Mrec_guess = Mp + Mn - targ%Mtar_struck
	  if (abs(targ%Mrec-Mrec_guess).gt.0.1) write(6,*)
     >		'WARNING: forcing recoil mass to M_Nucleon'
	  targ%Mrec = Mrec_guess
	endif


	targ%thick=targ%thick/1000.
	targ%length=targ%thick/targ%rho
	targ%angle=targ%angle/degrad
	targ_Bangle = targ_Bangle/degrad
	targ_Bphi = targ_Bphi/degrad
	if(targ%Z.lt.2.4) then		!liquid target
	  if (abs(targ%angle) .gt. 0.0001) then
	    targ%angle = 0.0
	    write(6,*) 'Forcing target angle to zero for cryotarget.'
	  endif
	  if (targ%can.ne.1 .and. targ%can.ne.2) stop 'bad targ.can value'
	endif
	if(sin(targ%angle) .gt. 0.85) then
	  write(6,*) 'BAD targ.angle (0 is perp. to beam, +ve is rotated towards SOS)'
	  write(6,*) 'If the target really is >60 degrees from normal to the beam, then'
	  write(6,*) 'comment out this error checking code.'
	  stop
	endif

! ... dbase field spect_offset

	spec%e%offset%xptar=spec%e%offset%xptar/1000.
	spec%e%offset%yptar=spec%e%offset%yptar/1000.
	spec%p%offset%xptar=spec%p%offset%xptar/1000.
	spec%p%offset%yptar=spec%p%offset%yptar/1000.

! ... dbase fields e_arm_accep, p_arm_accep.  Convert angles to radians.

	if (SPedge%e%delta%min.le.-100.0) SPedge%e%delta%min=-99.99
	if (SPedge%p%delta%min.le.-100.0) SPedge%p%delta%min=-99.99
	SPedge%e%yptar%min = SPedge%e%yptar%min/1000.
	SPedge%e%yptar%max = SPedge%e%yptar%max/1000.
	SPedge%e%xptar%min = SPedge%e%xptar%min/1000.
	SPedge%e%xptar%max = SPedge%e%xptar%max/1000.
	SPedge%p%yptar%min = SPedge%p%yptar%min/1000.
	SPedge%p%yptar%max = SPedge%p%yptar%max/1000.
	SPedge%p%xptar%min = SPedge%p%xptar%min/1000.
	SPedge%p%xptar%max = SPedge%p%xptar%max/1000.

! ... dbase field simulate

! ... Set flags for e,e',hadron radiation tails.

	doing_tail(1) = one_tail.eq.0 .or. one_tail.eq.1 .or.
     >		one_tail.eq.-2 .or. one_tail.eq.-3
	doing_tail(2) = one_tail.eq.0 .or. one_tail.eq.2 .or.
     >		one_tail.eq.-3 .or. one_tail.eq.-1
	doing_tail(3) = one_tail.eq.0 .or. one_tail.eq.3 .or.
     >		one_tail.eq.-1 .or. one_tail.eq.-2

	if (.not.doing_tail(1) .or. .not.doing_tail(2)) then
	  write(6,*) '**WARNING: Shutting off tails 1 or 2 not well defined'
	endif

	if (abs(one_tail).gt.3 .and. using_rad)
     >     stop 'Moron! one_tail>3 turns radiation off, but using_rad wants it on.'

	if (.not.using_rad) then
	  do i = 1, 3
	    doing_tail(i) = .false.
	  enddo
	endif

	if (Egamma_gen_max.gt.0.01) then !use hardwired limits if gen_max > 0
	  hardwired_rad = .true.
	else
	  hardwired_rad = .false.
	endif

! ... change angles to radians
	gen%e%yptar%min=gen%e%yptar%min/1000.
	gen%e%yptar%max=gen%e%yptar%max/1000.
	gen%p%yptar%min=gen%p%yptar%min/1000.
	gen%p%yptar%max=gen%p%yptar%max/1000.
	gen%e%xptar%min=gen%e%xptar%min/1000.
	gen%e%xptar%max=gen%e%xptar%max/1000.
	gen%p%xptar%min=gen%p%xptar%min/1000.
	gen%p%xptar%max=gen%p%xptar%max/1000.

! ... set some flags
	using_E_arm_montecarlo = spect_mode.ne.1 .and. spect_mode.ne.-1
	using_P_arm_montecarlo = spect_mode.ne.1 .and. spect_mode.ne.-2

	if (doing_pion .or. doing_kaon .or. doing_delta .or.
     >		(cuts%Em%min.eq.cuts%Em%max) ) then
	  cuts%Em%min = -1.d6
	  cuts%Em%max =  1.d6
	endif

	if (abs(deForest_flag).gt.1) stop 'Idiot! check setting of deForest_flag'

	if (correct_Eloss .and. .not. using_Eloss) then
	  write(6,*) 'using_Eloss is false, SO WE ARE FORCING CORRECT_ELOSS TO BE FALSE!!!'
	  correct_Eloss = .false.
	endif

	if (nint(targ%Z).eq.1) using_Coulomb=.false.  !no coulomb corr. for Z=1

! ... target

	call target_init(using_Eloss, using_Coulomb, spec%e%theta, spec%p%theta,
     >		spec%p%P, Mh2, Ebeam, spec%e%P)

! ... Read in the theory file (for A(e,e'p))
	if (doing_deuterium .or. doing_heavy) then
	  call theory_init(success)
	  if (.not.success) stop 'THEORY_INIT failed!'
	endif

! ... Read in (and normalize) momentum distribution.  Normalizing to one
! ... may not be correct if there is strength beyond p=p_max (MeV/c)

	if(doing_deutpi.or.doing_hepi.or.doing_deutkaon.or.doing_hekaon.or.doing_deutsemi) then
	  if(doing_deutpi .or. doing_deutkaon .or. doing_deutsemi) open(1,file='deut.dat',status='old',form='formatted')
	  if(doing_hepi .or. doing_hekaon) then
	    if (nint(targ%A).eq.3) then
	      open(1,file='he3.dat',status='old',form='formatted')
	    else if (nint(targ%A).eq.4) then
	      open(1,file='he4.dat',status='old',form='formatted')
	    else if (nint(targ%A).eq.12) then
	      open(1,file='c12.dat',status='old',form='formatted')
	    else
	      write(6,*) 'No Momentum Distribution for A = ',targ%A
	      write(6,*) 'Defaulting to carbon momentum distribution'
	      open(1,file='c12.dat',status='old',form='formatted')
	    endif
	  endif
	  do ii=1,2000
	    read(1,*,end=999) pval(ii),mprob(ii)
	    nump=ii
	  enddo
999	  continue
	  do ii=1,nump
	    mprob(ii)=mprob(ii)/mprob(nump)
	  enddo
	  close(1)
	endif

! ... Read in (and normalize) spectral function for A>2 electroproduction.
! ... Normalizing to one may not be correct if there is strength beyond
! ... Em_max,Pm_max.

	if(doing_hepi.or.doing_hekaon .or. (doing_heavy.and.use_sf)) then
	  if (nint(targ%A).eq.3) then
	    if (sf_version.eq.0) then
	      write(6,*) 'Using the Benhar version of 3He S.F.'
	    else if (sf_version.eq.1) then
	      write(6,*) 'Using the Kaptari version of the 3He S.F.'
	    else
	      write(6,*) 'No proper S.F. specified. Defaulting to Benhar S.F.'
	    endif
! ============================================================
! RCT 10/26/2016 - using different spectral function versions
!           Benhar spectral function
	    if (sf_version.eq.0) then
	      tmpfile='benharsf_3mod.dat'
!           Kaptari spectral function
	    else if (sf_version.eq.1) then
	      tmpfile='kaptarisf_3par.dat'
	    else
	      write(6,*) 'Wrong value for sf_version. Defaulting to Benhar'
	      tmpfile='benharsf_3mod.dat'
	    endif
! ============================================================
	  else if (nint(targ%A).eq.4) then
	    tmpfile='benharsf_4.dat'
	  else if (nint(targ%A).eq.12) then
	    tmpfile='benharsf_12.dat'
	  else if (nint(targ%A).eq.56) then
	    tmpfile='benharsf_56.dat'
	  else if (nint(targ%A).eq.197) then
	    tmpfile='benharsf_197.dat'
	  else
	    write(6,*) 'No spectral function for A = ',targ%A
	    write(6,*) 'Defaulting to Carbon S.F.'
	    tmpfile='benharsf_12.dat'
	  endif

! *****************************************************************************
! RCT 9/13/2016 Added this little piece of code that takes the neutron
! 3He distribution and uses it for the proton distribution in tritium.
	  if ((nint(targ%A).eq.3).and.(nint(targ%Z).eq.2)) then
	    call sf_lookup_init(tmpfile,.true.)                 !the true flag calls the proton S.F.
	  else if ((nint(targ%A).eq.3).and.(nint(targ%Z).eq.1)) then
	    call sf_lookup_init(tmpfile,.false.)                !the false flag calls the neutron S.F.
	  endif
! Choos proton or neutron spectral function based on targ.Mtar_struck
!	  if (abs(targ%Mtar_struck-Mp).le.1.d-6) then
!	    call sf_lookup_init(tmpfile,.true.)			!proton S.F.
!	  else if (abs(targ%Mtar_struck-Mn).le.1.d-6) then
!	    call sf_lookup_init(tmpfile,.false.)		!neutron S.F.
!	  else
!	    write(6,*) 'targ%Mtar_struck = ',targ%Mtar_struck
!	    write(6,*) 'targ%Mtar_struck not equal to Mp or Mn'
!	    write(6,*) 'DEFAULTING TO PROTON SPECTRAL FUNCTION!!!'
!	    call sf_lookup_init(tmpfile,.true.)			!proton S.F.
!	   endif
	endif
! *****************************************************************************

! ... Read in input file for kaon cross section model (saghai model).
	if (doing_kaon) then		!read in model xsec info.
	  open(1,file='weightc.dat',status='old',form='formatted')
	  do i=1,50
	    do j=1,20
	      read(1,*) weightc(j,i)
	    enddo
	  enddo
	  close(1)

	  open(1,file='weightd.dat',status='old',form='formatted')
	  do i=1,30
	    do j=1,40
	      do k=1,8
	        read(1,*) weightd(k,j,i)
	      enddo
	    enddo
	  enddo
	  close(1)

	  open(3,file='saghai_proton.dat',status='old',form='formatted')
	  do iread=1,10
	    do iq2=1,11
	      read(3,*) dum1,dum2
	      do iang=1,19
		read(3,4005) dum3,dum4,dum5,dum6,dum7
		read(3,4005) zrff1(iread,iq2,iang),ziff1(iread,iq2,iang),
     1			zrff2(iread,iq2,iang),ziff2(iread,iq2,iang),
     2			zrff3(iread,iq2,iang),ziff3(iread,iq2,iang)
		read(3,4005) zrff4(iread,iq2,iang),ziff4(iread,iq2,iang),
     1			zrff5(iread,iq2,iang),ziff5(iread,iq2,iang),
     2			zrff6(iread,iq2,iang),ziff6(iread,iq2,iang)
	      enddo
	    enddo
	  enddo
4005	  format(6e12.4)
	  close(3)

	  open(3,file='saghai_sigma0.dat',status='old',form='formatted')
	  do iread=1,20
	    do iq2=1,10
	      do iang=1,19
		read(3,*) dum1,dum2
		read(3,4005) dum3,dum4,dum5,dum6,dum7
		read(3,4005) zsrff1(iread,iq2,iang),zsiff1(iread,iq2,iang),
     1			zsrff2(iread,iq2,iang),zsiff2(iread,iq2,iang),
     2			zsrff3(iread,iq2,iang),zsiff3(iread,iq2,iang)
		read(3,4005) zsrff4(iread,iq2,iang),zsiff4(iread,iq2,iang),
     1			zsrff5(iread,iq2,iang),zsiff5(iread,iq2,iang),
     2			zsrff6(iread,iq2,iang),zsiff6(iread,iq2,iang)
	      enddo
	    enddo
	  enddo
	  close(3)

	endif		!kaon input files.



! ... initialize limits on generation, cuts, edges

	call limits_init(H)

! ... some announcements

	if (doing_eep) then
	  if (doing_hyd_elast) then
	    write(6,*) ' ****--------  H(e,e''p)  --------****'
	  else if (doing_deuterium) then
	    write(6,*) ' ****--------  D(e,e''p)  --------****'
	  else if (doing_heavy) then
	    write(6,*) ' ****--------  A(e,e''p)  --------****'
	  else
	    stop 'I don''t have ANY idea what (e,e''p) we''re doing!!!'
	  endif

	else if (doing_semi) then
           if (doing_semipi) then
	      if(doing_hydsemi) then
		 if(doing_hplus) then
		    write(6,*) ' ****--------  H(e,e''pi+)X  --------****'
		 else
		    write(6,*) ' ****--------  H(e,e''pi-)X  --------****'
		 endif
	      elseif (doing_deutsemi) then
		 if(doing_hplus) then
		    write(6,*) ' ****--------  D(e,e''pi+)X  --------****'
		 else
		    write(6,*) ' ****--------  D(e,e''pi-)X  --------****'
		 endif
	      endif

	   else if (doing_semika) then
	      if (doing_hydsemi) then
		 if(doing_hplus) then
		    write(6,*) ' ****--------  H(e,e''K+)X  --------****'
		 else
		    write(6,*) ' ****--------  H(e,e''K-)X  --------****'
		 endif
	      elseif(doing_deutsemi) then
		 if(doing_hplus) then
		    write(6,*) ' ****--------  D(e,e''K+)X  --------****'
		 else
		    write(6,*) ' ****--------  D(e,e''K-)X  --------****'
		 endif
	      endif
	   else
	      stop 'I don''t have ANY idea what (e,e''m)X we''re doing!!!'
	   endif

	else if (doing_delta) then
	  if (doing_hyddelta) then
	    write(6,*) ' ****--------  H(e,e''p)pi  --------****'
	  else if (doing_deutpi) then
	    write(6,*) ' ****--------  D(e,e''p)pi  --------****'
	  else if (doing_hepi) then
	    write(6,*) ' ****--------  A(e,e''p)pi  --------****'
	  else
	    stop 'I don''t have ANY idea what (e,e''p)pi we''re doing!!!'
	  endif
	else if (doing_pion) then
	  if (doing_hydpi) then
	    if (targ%A .eq. 1) then
	      write(6,*) ' ****--------  H(e,e''pi)  --------****'
	    else if (targ%A .ge. 3) then
	      write(6,*) ' ****--------  A(e,e''pi)  --------****'
	    endif
	  else if (doing_deutpi) then
	    write(6,*) ' ****--------  D(e,e''pi)  --------****'
	  else if (doing_hepi) then
	    write(6,*) ' ****--------  A(e,e''pi)  --------****'
	  else
	    stop 'I don''t have ANY idea what (e,e''pi) we''re doing!!!'
	  endif
	  if (which_pion.eq.0 .or. which_pion.eq.10) then
	    write(6,*) ' ****-------  pi+ production  -------****'
	  else if (which_pion.eq.1 .or. which_pion.eq.11) then
	    write(6,*) ' ****-------  pi- production  -------****'
	  endif
	  if (which_pion.eq.0 .or. which_pion.eq.1) then
	    write(6,*) ' ****---- Quasifree Production ----****'
	  else if (which_pion.eq.10 .or. which_pion.eq.11) then
	    write(6,*) ' ****----  Coherent Production ----****'
	  endif
	else if (doing_kaon) then
	  if (doing_hydkaon) then
	    if (targ%A .eq. 1) then
	      write(6,*) ' ****--------  H(e,e''K)  --------****'
	    else if (targ%A .ge. 3) then
	      write(6,*) ' ****--------  A(e,e''K)  --------****'
	    endif
	  else if (doing_deutkaon) then
	    write(6,*) ' ****--------  D(e,e''K)  --------****'
	  else if (doing_hekaon) then
	    write(6,*) ' ****--------  A(e,e''K)  --------****'
	  else
	    stop 'I don''t have ANY idea what (e,e''K) we''re doing!!!'
	  endif
	  if (which_kaon.eq.0) then
	    write(6,*) ' ****----- LAMBDA Quasifree -----****'
	  else if (which_kaon.eq.1) then
	    write(6,*) ' ****----- SIGMA0 Quasifree -----****'
	  else if (which_kaon.eq.2) then
	    write(6,*) ' ****----- SIGMA- Quasifree -----****'
	  else if (which_kaon.eq.10) then
	    write(6,*) ' ****----- LAMBDA Coherent  -----****'
	  else if (which_kaon.eq.11) then
	    write(6,*) ' ****----- SIGMA0 Coherent  -----****'
	  else if (which_kaon.eq.12) then
	    write(6,*) ' ****----- SIGMA- Coherent  -----****'
	  else
	    stop 'I don''t have ANY idea what (e,e''K) we''re doing!!!'
	  endif
	else if (doing_rho) then
	   if (doing_hydrho) then
	      if (targ%A .eq. 1) then
		 write(6,*) ' ****--------  H(e,e''rho)  --------****'
	      else if (targ%A .ge. 3) then
		 write(6,*) ' ****--------  A(e,e''rho)  --------****'
		 write(6,*) ' **** ---- Not yet implemented -----****'
		 stop
	      endif
	   else if (doing_deutrho) then
	      write(6,*) ' ****--------  D(e,e''rho)  --------****'
	      write(6,*) ' **** ---- Not yet implemented -----****'
	      stop
	   else if (doing_herho) then
	      write(6,*) ' ****--------  A(e,e''rho)  --------****'
	      write(6,*) ' **** ---- Not yet implemented -----****'
	      stop
	   else
	      stop 'I don''t have ANY idea what (e,e''rho) we''re doing!!!'
	   endif

	else if (doing_phsp) then
	  write(6,*) ' ****--------  PHASE SPACE - NO physics, NO radiation  --------****'
	else
	  stop 'I don''t have ANY idea what we''re doing!!!'
	endif

	if (electron_arm.eq.1) then
	  write(6,*) '   HMS  is detecting electrons'
	else if (electron_arm.eq.2) then
	  write(6,*) '   SOS  is detecting electrons'
	else if (electron_arm.eq.3) then
	  write(6,*) '   HRSR is detecting electrons'
	else if (electron_arm.eq.4) then
	  write(6,*) '   HRSL is detecting electrons'
	else if (electron_arm.eq.5) then
	  write(6,*) '   SHMS is detecting electrons (SSA TUNE)'
	else if (electron_arm.eq.6) then
	  write(6,*) '   SHMS is detecting electrons (LSA TUNE)'
	else if (electron_arm.eq.7) then
	  write(6,*) '   BIGCAL is detecting electrons (HMS side of beam)'
	else if (electron_arm.eq.8) then
	  write(6,*) '   BIGCAL is detecting electrons (SOS side of beam)'
	endif

	if (hadron_arm.eq.1) then
	  write(6,*) '   HMS  is detecting hadrons'
	else if (hadron_arm.eq.2) then
	  write(6,*) '   SOS  is detecting hadrons'
	else if (hadron_arm.eq.3) then
	  write(6,*) '   HRSR is detecting hadrons'
	else if (hadron_arm.eq.4) then
	  write(6,*) '   HRSL is detecting hadrons'
	else if (hadron_arm.eq.5) then
	  write(6,*) '   SHMS is detecting hadrons (SSA TUNE)'
	else if (hadron_arm.eq.6) then
	  write(6,*) '   SHMS is detecting hadrons (LSA TUNE)'
	else if (hadron_arm.eq.7 .or. hadron_arm.eq.8) then
	  write(6,*) ' Cannot use Bigcal for hadrons'
	  stop
	endif

	if (hadron_arm.eq.electron_arm) then
	  write(6,*) '**WARNING: SIMC works best with two DIFFERENT spectrometers!!!'
	else if ((hadron_arm+electron_arm).eq.3) then
	  write(6,*) 'Rolf welcomes you to Hall C at Jefferson Lab - Get to work!'
	else if (electron_arm.eq.5 .or. hadron_arm.eq.5 .or.
     >		 electron_arm.eq.6 .or. hadron_arm.eq.6) then
	  write(6,*) 'Grumpy Dutch Giant welcomes you to Hall C++ at JLab12 - Get to work!'
	else if ((hadron_arm+electron_arm).eq.7) then
	  write(6,*) '------------------------------------------------------'
	  write(6,*) 'Kees welcomes you to Hall A at JLab - Enjoy your stay!'
	  write(6,*) '------------------------------------------------------'
	else if ( electron_arm.eq.7 .and. hadron_arm .eq. 5) then
	  write(6,*) ' Bigcal and SHMS'
	else if ( electron_arm.eq.8 .and. hadron_arm .eq. 1) then
	  write(6,*) ' Bigcal and HMS'
	else
	  write(6,*) '**WARNING: SIMC works best when both spectrometers are in the same hall!!'
	endif

	if (doing_phsp) write(6,*) '**WARNING: phsp option is not really implemented - buyer beware!!'
	if (doing_kaon .and. .not. doing_decay) write(6,*) 'NOTE: not doing decay, so a decay weight will be applied to WEIGHT'
	if (doing_semika .and. .not. doing_decay) write(6,*) 'NOTE: not doing decay, so a decay weight will be applied to WEIGHT'
	if(doing_deutsemi.and.do_fermi) write(6,*) 'NOTE: Fermi motion enabled for semi-incusive production from deuterium'
	if (.not.using_rad) write(6,*) 'NOTE: Will NOT be applying radiative corrections'
	if (.not.using_E_arm_montecarlo) write(6,*) 'NOTE: Will NOT be running events through the E arm Monte Carlo'
	if (.not.using_P_arm_montecarlo) write(6,*) 'NOTE: Will NOT be running events through the P arm Monte Carlo'
	if (.not.using_Eloss) write(6,*) 'NOTE: Will NOT be calculating energy loss in the target'
	if (.not.standard_Eloss) write(6,*) 'NOTE: Will be calculating full energy loss in the target for recon'
	if (.not.using_Coulomb) write(6,*) 'NOTE: Will NOT be calculating Coulomb correction (default for Hydrogen target)'
	if (using_Coulomb) write(6,*) 'NOTE: WILL be calculating Coulomb corrections.
     >     (Implemented for beam and scattered electron only)'
	if (.not.correct_Eloss) write(6,*) 'NOTE: Will NOT correct reconstructed data for energy loss'
	if (.not.correct_raster) write(6,*) 'NOTE: Will NOT use raster terms in reconstruction'
	write(6,*) ''
	write(6,*) ' *************************************************'
! everything has been set now dump the parameters
	call dump_parameters(dump3)

	return
	end

!---------------------------------------------------------------

	subroutine get_defaults
	implicit none
	include 'dbase_namelists.inc'
! get the initialization values
! the rewind statments make sure that the sequence of namelists do not matter
!
! these namlists MUST occur in the nml_default.data file for proper
! initialization of the program

	if (debug(2)) write(6,*)'default namelist: entering'
! open defaults
	open(io_nl, file = 'nml_default.data', status = 'old')
	read( io_nl, NML = RESTSW, iostat = iret)
	rewind(io_nl)
	if (iret .ne. 0) print *, 'get_defaults: read error :  RESTSW', iret
	read( io_nl, NML = EXPERIMENT, iostat = iret)
	rewind(io_nl)
	if (iret .ne. 0) print *, 'get_defaults: read error :  EXPERIMENT', iret
	read( io_nl, NML = DEBUG_PARM, iostat = iret)
	rewind(io_nl)
	if (iret .ne. 0) print *, 'get_defaults: read error :  DEBUG_PARM', iret
	read( io_nl, NML = TARGET, iostat = iret)
	rewind(io_nl)
	if (iret .ne. 0) print *, 'get_defaults: read error :  TARGET', iret
	read( io_nl, NML = E_ARM_MAIN, iostat = iret)
	rewind(io_nl)
	if (iret .ne. 0) print *, 'get_defaults: read error : E_ARM_MAIN ', iret
	read( io_nl, NML = P_ARM_MAIN, iostat = iret)
	rewind(io_nl)
	if (iret .ne. 0) print *, 'get_defaults: read error :  P_ARM_MAIN', iret
	read( io_nl, NML = E_ARM_OFFSET, iostat = iret)
	rewind(io_nl)
	if (iret .ne. 0) print *, 'get_defaults: read error :  E_ARM_OFFSET', iret
	read( io_nl, NML = P_ARM_OFFSET, iostat = iret)
	rewind(io_nl)
	if (iret .ne. 0) print *, 'get_defaults: read error :  P_ARM_OFFSET', iret
	read( io_nl, NML = MISC2INT, iostat = iret)
	rewind(io_nl)
	if (iret .ne. 0) print *, 'get_defaults: read error :  MISC2INT', iret
	read( io_nl, NML = SIMULATE, iostat = iret)
	rewind(io_nl)
	if (iret .ne. 0) print *, 'get_defaults: read error :  SIMULATE', iret
	read( io_nl, NML = E_ARM_ACCEPT, iostat = iret)
	rewind(io_nl)
	if (iret .ne. 0) print *, 'get_defaults: read error :  E_ARM_ACCEPT', iret
	read( io_nl, NML = P_ARM_ACCEPT, iostat = iret)
	rewind(io_nl)
	if (iret .ne. 0) print *, 'get_defaults: read error :  P_ARM_ACCEPT', iret
	read( io_nl, NML = THEORY_CTRL, iostat = iret)
	rewind(io_nl)
	if (iret .ne. 0) print *, 'get_defaults: read error :  THEORY', iret

	close(io_nl)

	if (debug(2)) write(6,*)'user namelist: entering'
	return
	end
!---------------------------------------------------------------
	integer function read_parameters(filename)

! find which namelists are in the parameter (dbase) file and
! in what order. Then read the corresponding namelists
!
! for this to work the names stores in dbase_namelists.inc MUST agree
! with the namelists that have been defined

	implicit none
	character*80 line, filename, current_name
	character *80 dump
	logical always
	integer nl_counter
	include 'dbase_namelists.inc'
!       loop indices
	integer i, j
! it is assumed that each namelist occurs only once in the input field
	integer nl_array(nnames)


	data dump/'dump_rp.dat'/
	read_parameters = -1

! find the namelist sequence in the input file
	open(io_nl, file = filename)
	nl_counter = 0
! endless loop
	input:do
	   read(io_nl, '(a)', iostat = iret) line
	   if (iret .ne. 0) then
	      exit
	   endif
	   i = index(line,'&')! namelists start with an &
	   if (i .ne. 0) then ! found a namelist
	      current_name = trim(line(i+1:))
	      do j = 1, nnames
		 if (nl_names(j) == current_name) then
		    nl_counter = nl_counter + 1
		    ! store the name index
		    if (nl_counter .gt. nnames) then
		       write (6, *) 'too many namelists in data file :' // filename
		       return
		    endif
		    nl_array(nl_counter) = j
		 endif
	      enddo
	   endif
	enddo input
 999	rewind(io_nl)
	if (nl_counter .eq. 0) then
	   write (6,*) 'read_parameters: nothing to read from : '//filename
	   return
	endif
! now read the namelists according to what is in the file
	do i = 1, nl_counter
	   select case( nl_array(i) )
	   case (1)
	      read( io_nl, NML = RESTSW, iostat = iret)
	      if (iret .ne. 0) print *, 'read_parameters: read error : RESTSW', iret
	   case (2)
	      read( io_nl, NML = EXPERIMENT, iostat = iret)
	      if (iret .ne. 0) print *, 'read_parameters: read error :  EXPERIMENT', iret
	   case (3)
	      read( io_nl, NML = DEBUG_PARM, iostat = iret)
	      if (iret .ne. 0) print *, 'read_parameters: read error :  DEBUG_PARM', iret
	   case (4)
	      read( io_nl, NML = TARGET, iostat = iret)
	      if (iret .ne. 0) print *, 'read_parameters: read error : TARGET ', iret
	   case (5)
	      read( io_nl, NML = E_ARM_MAIN, iostat = iret)
	      if (iret .ne. 0) print *, 'read_parameters: read error :  E_ARM_MAIN', iret
	   case (6)
	      read( io_nl, NML = P_ARM_MAIN, iostat = iret)
	      if (iret .ne. 0) print *, 'read_parameters: read error :  P_ARM_MAIN ', iret
	   case (7)
	      read( io_nl, NML = E_ARM_OFFSET, iostat = iret)
	      if (iret .ne. 0) print *, 'read_parameters: read error :   E_ARM_OFFSET', iret
	   case (8)
	      read( io_nl, NML = P_ARM_OFFSET, iostat = iret)
	      if (iret .ne. 0) print *, 'read_parameters: read error :  P_ARM_OFFSET', iret
	   case (9)
	      read( io_nl, NML = MISC2INT, iostat = iret)
	      if (iret .ne. 0) print *, 'read_parameters: read error :  MISC2INT', iret
	   case (10)
	      read( io_nl, NML = SIMULATE, iostat = iret)
	      if (iret .ne. 0) print *, 'read_parameters: read error :  SIMULATE', iret
	   case (11)
	      read( io_nl, NML = E_ARM_ACCEPT, iostat = iret)
	      if (iret .ne. 0) print *, 'read_parameters: read error : E_ARM_ACCEPT', iret
	   case (12)
	      read( io_nl, NML = P_ARM_ACCEPT, iostat = iret)
	      if (iret .ne. 0) print *, 'read_parameters: read error :  P_ARM_ACCEPT', iret
	   case (13)
	      read( io_nl, NML = KINEMATICS_MAIN, iostat = iret)
	      if (iret .ne. 0) print *, 'read_parameters: read error :  KINEMATICS_MAIN', iret
	   case (14)
	      read( io_nl, NML = BEAM_AND_TARGET_INFO, iostat = iret)
	      if (iret .ne. 0) print *, 'read_parameters: read error :  BEAM_AND_TARGET_INFO', iret
	   case (15)
	      read( io_nl, NML = THEORY_CTRL, iostat = iret)
	      if (iret .ne. 0) print *, 'read_parameters: read error :  THEORY', iret
	   end select
	enddo
	close( io_nl)
	if (iret .ne. 0) then
	   read_parameters = -1
	else
!	   call dump_parameters(dump)
	   read_parameters = 0
	endif
	return
	end


!----------------------------------------------------------------------
	subroutine dump_parameters( dump_file)
	implicit none
	include 'dbase_namelists.inc'
	character*80 dump_file

! open dump file
	open(io_nl, file = dump_file)
	write( io_nl, NML = RESTSW, iostat = iret)
	write( io_nl, NML = EXPERIMENT, iostat = iret)
	write( io_nl, NML = DEBUG_PARM, iostat = iret)
	write( io_nl, NML = TARGET, iostat = iret)
	write( io_nl, NML = E_ARM_MAIN, iostat = iret)
	write( io_nl, NML = P_ARM_MAIN, iostat = iret)
	write( io_nl, NML = E_ARM_OFFSET, iostat = iret)
	write( io_nl, NML = P_ARM_OFFSET, iostat = iret)
	write( io_nl, NML = MISC2INT, iostat = iret)
	write( io_nl, NML = SIMULATE, iostat = iret)
	write( io_nl, NML = E_ARM_ACCEPT, iostat = iret)
	write( io_nl, NML = P_ARM_ACCEPT, iostat = iret)
	write( io_nl, NML = KINEMATICS_MAIN, iostat = iret)
	write( io_nl, NML = BEAM_AND_TARGET_INFO, iostat = iret)
	write( io_nl, NML = THEORY_CTRL, iostat = iret)

	close(io_nl)
	return
	end
