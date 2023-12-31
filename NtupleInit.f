	subroutine NtupleInit(filename)
	implicit none
	save

	include  'hbook.inc'
	include  'simulate.inc'

	character*80 filename,directory
	character*30 name,title
	integer*4 m,io,recl,bank,id,status, i
	parameter(recl = 1024)
	parameter(bank = 8000)
	parameter(io = 29)
	parameter(name = 'SNT')
	parameter(title = 'SIMTUPLE')
	
	data NtupleTag /80*' '/

	NtupleID = defaultID
	id = NtupleID
	NtupleIO = io
	NtupleName = name

!	call HCDIR(directory,'R') !CERNLIB read current directory
!	call HROPEN(io,name,filename,'N',recl,status)  !CERNLIB
	open (io, file=filename, status='UNKNOWN', iostat=status)
						!directory set to "//TUPLE"
	if (status.ne.0)then
	  write(6,*) 'NtupleInit: OPEN error: istat=',status
	  stop
	endif

	m = 0
	if(electron_arm.eq.1 .or. electron_arm.eq.3.or. electron_arm.eq.7)then !electron = right side.
	   m = m+1
	   NtupleTag(m) = 'e_delta' !  1
	   m = m+1
	   NtupleTag(m) = 'e_yptar' !  2
	   m = m+1
	   NtupleTag(m) = 'e_xptar' !  3
	   m = m+1
	   NtupleTag(m) = 'e_ytar' !  4
	   m = m+1
	   NtupleTag(m) = 'e_xfp' !  5
	   m = m+1
	   NtupleTag(m) = 'e_xpfp' !  6
	   m = m+1
	   NtupleTag(m) = 'e_yfp' !  7
	   m = m+1
	   NtupleTag(m) = 'e_ypfp' !  8
	   m = m+1
	   NtupleTag(m) = 'e_deltai' !  9
	   m = m+1
	   NtupleTag(m) = 'e_yptari' ! 10
	   m = m+1
	   NtupleTag(m) = 'e_xptari' ! 11
	   m = m+1
	   NtupleTag(m) = 'e_ytari' ! 12
	   m = m+1
	   NtupleTag(m) = 'h_delta' ! 13
	   m = m+1
	   NtupleTag(m) = 'h_yptar' ! 14
	   m = m+1
	   NtupleTag(m) = 'h_xptar' ! 15
	   m = m+1
	   NtupleTag(m) = 'h_ytar' ! 16
	   m = m+1
	   NtupleTag(m) = 'h_xfp' ! 17
	   m = m+1
	   NtupleTag(m) = 'h_xpfp' ! 18
	   m = m+1
	   NtupleTag(m) = 'h_yfp' ! 19
	   m = m+1
	   NtupleTag(m) = 'h_ypfp' ! 20
	   m = m+1
	   NtupleTag(m) = 'h_deltai' ! 21
	   m = m+1
	   NtupleTag(m) = 'h_yptari' ! 22
	   m = m+1
	   NtupleTag(m) = 'h_xptari' ! 23
	   m = m+1
	   NtupleTag(m) = 'h_ytari' ! 24
	else if (electron_arm.eq.2 .or. electron_arm.eq.4 .or.
     >		 electron_arm.eq.5 .or. electron_arm.eq.6.or. electron_arm.eq.8) then  !e- = left.
	   m = m+1
	   NtupleTag(m) = 'h_delta' !  1
	   m = m+1
	   NtupleTag(m) = 'h_yptar' !  2
	   m = m+1
	   NtupleTag(m) = 'h_xptar' !  3
	   m = m+1
	   NtupleTag(m) = 'h_ytar' !  4
	   m = m+1
	   NtupleTag(m) = 'h_xfp' !  5
	   m = m+1
	   NtupleTag(m) = 'h_xpfp' !  6
	   m = m+1
	   NtupleTag(m) = 'h_yfp' !  7
	   m = m+1
	   NtupleTag(m) = 'h_ypfp' !  8
	   m = m+1
	   NtupleTag(m) = 'h_deltai' !  9
	   m = m+1
	   NtupleTag(m) = 'h_yptari' ! 10
	   m = m+1
	   NtupleTag(m) = 'h_xptari' ! 11
	   m = m+1
	   NtupleTag(m) = 'h_ytari' ! 12

	   m = m+1
	   NtupleTag(m) = 'e_delta' ! 13
	   m = m+1
	   NtupleTag(m) = 'e_yptar' ! 14
	   m = m+1
	   NtupleTag(m) = 'e_xptar' ! 15
	   m = m+1
	   NtupleTag(m) = 'e_ytar' ! 16
	   m = m+1
	   NtupleTag(m) = 'e_xfp' ! 17
	   m = m+1
	   NtupleTag(m) = 'e_xpfp' ! 18
	   m = m+1
	   NtupleTag(m) = 'e_yfp' ! 19
	   m = m+1
	   NtupleTag(m) = 'e_ypfp' ! 20
	   m = m+1
	   NtupleTag(m) = 'e_deltai' ! 21
	   m = m+1
	   NtupleTag(m) = 'e_yptari' ! 22
	   m = m+1
	   NtupleTag(m) = 'e_xptari' ! 23
	   m = m+1
	   NtupleTag(m) = 'e_ytari' ! 24
	endif
	m = m+1
	NtupleTag(m) = 'q'	! 25
	m = m+1
	NtupleTag(m) = 'nu'		! 26
	m = m+1
	NtupleTag(m) = 'Q2'		! 27
	m = m+1
	NtupleTag(m) = 'W'		! 28
	m = m+1
	NtupleTag(m) = 'epsilon'	! 29
	m = m+1
	NtupleTag(m) = 'Em'		! 30
	m = m+1
	NtupleTag(m) = 'Pm'		! 31
	m = m+1
	NtupleTag(m) = 'theta_pq'	! 32
	m = m+1
	NtupleTag(m) = 'phi_pq'		! 33

	if (doing_pion .or. doing_kaon .or. doing_delta) then
	  m = m+1
	  NtupleTag(m) = 'missmass'	! 34
	  m = m+1
	  NtupleTag(m) = 'mmnuc'	! 35
	  m = m+1
	  NtupleTag(m) = 'phad'		! 36
	  m = m+1
	  NtupleTag(m) = 't'		! 37
	  m = m+1
	  NtupleTag(m) = 'pmpar'	! 38
	  m = m+1
	  NtupleTag(m) = 'pmper'	! 39
	  m = m+1
	  NtupleTag(m) = 'pmoop'	! 40
	  m = m+1
	  NtupleTag(m) = 'fry'		! 41		!+y is up.
	  m = m+1
	  NtupleTag(m) = 'radphot'	! 42
	  m = m+1
	  NtupleTag(m) = 'pfermi'	! 43
	  m = m+1
	  NtupleTag(m) = 'siglab'	! 44
	  m = m+1
	  NtupleTag(m) = 'sigcm'	! 45
	  m = m+1
	  NtupleTag(m) = 'Weight'	! 46
	  m = m+1
	  NtupleTag(m) = 'decdist'	! 47
	  m = m+1
	  NtupleTag(m) = 'Mhadron'	! 48
	  m = m+1
	  NtupleTag(m) = 'pdotqhat'	! 49
	  m = m+1
	  NtupleTag(m) = 'Q2i'		! 50
	  m = m+1
	  NtupleTag(m) = 'Wi'		! 51
	  m = m+1
	  NtupleTag(m) = 'ti'		! 52
	  m = m+1
	  NtupleTag(m) = 'phipqi'	! 53
	  if(using_tgt_field) then
	     m = m+1
	     NtupleTag(m) = 'th_tarq' ! 54
	     m = m+1 
	     NtupleTag(m) = 'phitarq' ! 55 
	     m = m+1
	     NtupleTag(m) = 'beta' ! 56
	     m = m+1
	     NtupleTag(m) = 'phis' ! 57
	     m = m+1
	     NtupleTag(m) = 'phic' ! 58
	     m = m+1
	     NtupleTag(m) = 'betai' ! 59
	     m = m+1
	     NtupleTag(m) = 'phisi' ! 60
	     m = m+1
	     NtupleTag(m) = 'phici' ! 61
	  endif
	  if (doing_kaon) then
	    m = m+1
	    NtupleTag(m) = 'saghai'	! 54
	    m = m+1
	    NtupleTag(m) = 'factor'	! 55
	  endif
	else if (doing_semi.or.doing_rho) then
	  m = m+1
	  NtupleTag(m) = 'missmass'	! 34 <- Wprime for semi-inclusive folks
	  m = m+1
	  NtupleTag(m) = 'ppi'		! 35
	  m = m+1
	  NtupleTag(m) = 't'		! 36
	  m = m+1
	  NtupleTag(m) = 'fry'		! 37		!+y is up.
	  m = m+1
	  NtupleTag(m) = 'radphot'	! 38
	  m = m+1
	  NtupleTag(m) = 'siglab'	! 39
	  m = m+1
	  NtupleTag(m) = 'sigcent'	! 40
	  m = m+1
	  NtupleTag(m) = 'Weight'	! 41
	  m = m+1
	  NtupleTag(m) = 'decdist'	! 42
	  m = m+1
	  NtupleTag(m) = 'Mhadron'	! 43
	  m = m+1
	  NtupleTag(m) = 'z' 	        ! 44
	  m = m+1
	  NtupleTag(m) = 'zi' 	        ! 45
	  m = m+1
	  NtupleTag(m) = 'pt2' 	        ! 46
	  m = m+1
	  NtupleTag(m) = 'pt2i' 	! 47
	  m = m+1
	  NtupleTag(m) = 'xbj' 	        ! 48
	  m = m+1
	  NtupleTag(m) = 'xbji' 	! 49
	  m = m+1
	  NtupleTag(m) = 'thqi' 	! 50	  
	  m = m+1
	  NtupleTag(m) = 'sighad' 	! 51	  
	  m = m+1
	  NtupleTag(m) = 'jacobian' 	! 52	  
	  m = m+1
	  NtupleTag(m) = 'centjac' 	! 53
	  m = m+1
	  NtupleTag(m) = 'pfermi'       ! 54
	  m = m+1
	  NtupleTag(m) = 'xfermi'       ! 55
	  m = m+1
	  NtupleTag(m) = 'phipqi'       ! 56
	  if(using_tgt_field) then
	     m = m+1
	     NtupleTag(m) = 'th_tarq' ! 57
	     m = m+1 
	     NtupleTag(m) = 'phitarq' ! 58 
	     m = m+1
	     NtupleTag(m) = 'beta' ! 59
	     m = m+1
	     NtupleTag(m) = 'phis' ! 60
	     m = m+1
	     NtupleTag(m) = 'phic' ! 61
	     m = m+1
	     NtupleTag(m) = 'betai' ! 62
	     m = m+1
	     NtupleTag(m) = 'phisi' ! 63
	     m = m+1
	     NtupleTag(m) = 'phici' ! 64
	  endif
	  if(doing_rho) then
	     m = m+1
	     NtupleTag(m) = 'Mrho' ! 57 or 65
	     m = m+1
	     NtupleTag(m) = 'Thrho' ! 58 or 66
	  endif
	    
	else if (doing_hyd_elast .or. doing_deuterium .or. doing_heavy) then
	  m = m+1
	  NtupleTag(m) = 'corrsing'	! 34
	  m = m+1
	  NtupleTag(m) = 'Pmx'		! 35		!for Heepcheck
	  m = m+1
	  NtupleTag(m) = 'Pmy'		! 36		!for Heepcheck
	  m = m+1
	  NtupleTag(m) = 'Pmz'		! 37		!for Heepcheck
	  m = m+1
	  NtupleTag(m) = 'PmPar'	! 38
	  m = m+1
	  NtupleTag(m) = 'PmPer'	! 39
	  m = m+1
	  NtupleTag(m) = 'PmOop'	! 40
	  m = m+1
	  NtupleTag(m) = 'fry'		! 41		!+y is up.
	  m = m+1
	  NtupleTag(m) = 'radphot'	! 42
	  m = m+1
	  NtupleTag(m) = 'sigcc'	! 43
	  m = m+1
	  NtupleTag(m) = 'Weight'	! 44
	  m = m+1
	  NtupleTag(m) = 'Jacobian'	! 45          ! WB added 
	  m = m+1
	  NtupleTag(m) = 'theta_e'	! 46          ! WB added 
	  m = m+1
	  NtupleTag(m) = 'theta_p'	! 47          ! WB added 
	  m = m+1
	  NtupleTag(m) = 'tar_x'	! 48          ! WB added 
	  m = m+1
	  NtupleTag(m) = 'tar_y'	! 49          ! WB added 
	  m = m+1
	  NtupleTag(m) = 'tar_z'	! 50          ! WB added 
	  m = m+1
	  NtupleTag(m) = 'Genweight'	! 51 general weight
	  m = m+1
	  NtupleTag(m) = 'SF_weight'	! 52  spectral function   
	  m = m+1
	  NtupleTag(m) = 'Jacobian_corr' ! 53   correction to jac.        
	  m = m+1
	  NtupleTag(m) = 'sig'	        ! 54    full cross section          
	  m = m+1
	  NtupleTag(m) = 'sig_recon'	! 55   full cross section from reconstruction (Laget only) 
	  m = m+1
	  NtupleTag(m) = 'sigcc_recon'	! 56   sigma cc (el. cross section) from recon.
	  m = m+1
	  NtupleTag(m) = 'coul_corr'	! 57   coul. corr
 	  m = m+1
	  NtupleTag(m) = 'h_zv'	! 58 recon. vertex p
	  m = m+1
	  NtupleTag(m) = 'h_yv'	! 59 recon. vertex p
	  m = m+1
	  NtupleTag(m) = 'e_zv'	! 60 recon. vertex e
	  m = m+1
	  NtupleTag(m) = 'e_yv'	! 61 recon. vertex e
	  m = m+1
	  NtupleTag(m) = 'h_pf'	! 62 recon. proton final momentum
	  m = m+1
	  NtupleTag(m) = 'e_pf'	! 63 recon. electron final momentum
	  m = m+1
	  NtupleTag(m) = 'h_pfi'! 64 recon. electron final momentum
	  m = m+1
	  NtupleTag(m) = 'e_pfi'! 65 recon. electron final momentum
	  m = m+1
	  NtupleTag(m) = 'Ein'	! 66 recon. beam energy
	  m = m+1
	  NtupleTag(m) = 'theta_rq'	! 67 recoil angle
	  m = m+1
	  NtupleTag(m) = 'SF_weight_recon' ! 68 spectral function with reconstructed quantities
	  m = m+1 
	  NtupleTag(m) = 'theta_ei'        ! 69 RCT 5/26/2017 outgoing electron in-plane angle generated
	  m = m+1
	  NtupleTag(m) = 'theta_pi'        ! 70 RCT 5/26/2017 outgoing proton n in-plane angle generated
	  m = m+1
	  NtupleTag(m) = 'e_spec_p'        ! 71 RCT 5/26/2017 electron spectrometer central momentum
	  m = m+1
	  NtupleTag(m) = 'h_spec_p'        ! 72 RCT 5/26/2017 proton   spectrometer central momentum
	  m = m+1
	  NtupleTag(m) = 'e_spec_th'       ! 73 RCT 5/26/2017 electron spectrometer central angle
	  m = m+1
	  NtupleTag(m) = 'h_spec_th'       ! 74 RCT 5/26/2017 proton   spectrometer central angle  
	  m = m+1
	  NtupleTag(m) = 'xB'              ! 75 RCT 5/26/2017 Bjorken x
	endif

!	else		!used to be the if (doing_phsp) option.
!	 m=m+1
!	 NtupleTag(m)='gd'
!	 m=m+1
!	 NtupleTag(m)='gt'
!	 m=m+1
!	 NtupleTag(m)='gp'
!	 m=m+1
!	 NtupleTag(m)='gy'
!	 m=m+1
!	 NtupleTag(m)='rd'
!	 m=m+1
!	 NtupleTag(m)='rt'
!	 m=m+1
!	 NtupleTag(m)='rp'
!	 m=m+1
!	 NtupleTag(m)='ry'
!	 m=m+1
!	 NtupleTag(m)='w'
!	endif

	NtupleSize = m

!	call HBOOKN(id,title,NtupleSize,name,bank,NtupleTag) !create Ntuple

!	call HCDIR(NtupleDirectory,'R') !record Ntuple directory

!	call HCDIR(directory,' ')       !reset CERNLIB directory
! create the output header
	write (NtupleIO, '(a)') title
	write (NtupleIO, '(a)') name
	write (NtupleIO, *) m
	write (NtupleIO,  '(a)') (NtupleTag(i), i = 1, m)
	return
	end
