	include 'simulate.inc'
	include 'radc.inc'
	include 'g_dump_all_events.inc'

	integer*4 ierr, io_nl, iret
	integer*4 nnames
	integer*4 regparmint, regparmdouble
	integer*4 regparmintarray, regparmdoublearray
	integer*4 regparmstring

! default unit number for io
	parameter(io_nl = 56)
! 15 namelists are defined currently
	parameter(nnames = 15)
	character*80 nl_names(nnames)

!	integer*4 regparmreal
*	RESTSW

	namelist /RESTSW/
     >	   mc_smear,
     >	   electron_arm,
     >	   hadron_arm,
     >	   use_first_cer,
     >	   extra_dbase_file

*	EXPERIMENT
	namelist /EXPERIMENT/
     >	   ngen,
     >	   Ebeam,
     >	   dEbeam,
     >	   EXPER,
!      EXPER%charge,
     >     doing_bound,
     >     doing_phsp,
     >	   doing_kaon,
     >	   which_kaon,
     >	   doing_pion,
     >	   which_pion,
     >	   doing_delta,
     >	   doing_semi,
     >	   doing_hplus,
     >	   doing_rho,
     >	   doing_decay,
     >	   extra_dbase_file,
     >	   ctau,
     >     use_sf,
     >     sf_version,
     >     transparency

*	DEBUG
	namelist /DEBUG_PARM/
     >	    debug

*	TARGET
	namelist /TARGET/
     >	   targ,
 !    targ%A,
 !    targ%Z,
 !    targ%mass_amu,
 !    targ%mrec_amu,
 !    targ%rho,
 !    targ%thick,
 !    targ%xoffset,
 !    targ%yoffset,
 !    targ%fr_pattern,
 !    targ%fr1,
 !    targ%fr2,
 !    targ%zoffset,
 !    targ%angle,
 !    targ%abundancy,
 !    targ%can,
     >     targ_Bangle,
     >     targ_Bphi,
     >     targ_pol

*	E_ARM_MAIN
	namelist /E_ARM_MAIN/
     >     spec
!       spec%e%P,
!       spec%e%theta

*	P_ARM_MAIN
	namelist /P_ARM_MAIN/
     >	   spec
!       spec%p%theta

*	E_ARM_OFFSET
	namelist /E_ARM_OFFSET/
     >	   gen,
     >	   spec

!    gen%ywid,
!    gen%xwid,
!    spec%e%offset%x,
!    spec%e%offset%y,
!    spec%e%offset%z,
!    spec%e%offset%xptar,
!    spec%e%offset%yptar

*	P_ARM_OFFSET
	namelist /P_ARM_OFFSET/
     >     spec
! 	   spec%p%offset%x,
! 	   spec%p%offset%y,
! 	   spec%p%offset%z,
! 	   spec%p%offset%xptar,
! 	   spec%p%offset%yptar

*	MISC2INT
	namelist /MISC2INT/
     >	   use_expon,
     >	   one_tail,
     >     dump_all_in_ntuple,
     >	   intcor_mode

*	SIMULATE
	namelist /SIMULATE/
     >	   ngen,
     >	   mc_smear,
     >	   mc_smear_addl,
     >	   hard_cuts,
     >	   using_rad,
     >	   spect_mode,
     >	   doing_phsp,
     >     cuts,
!   cuts%Em%min,
!   cuts%Em%max,
     >	   using_Eloss,
     >	   correct_Eloss,
     >	   standard_Eloss,
     >	   correct_raster,
     >	   deForest_flag,
     >	   rad_flag,
     >	   extrad_flag,
     >	   lambda,
     >	   using_cit_generation,
     >	   Nntu,
     >	   using_Coulomb,
     >	   dE_edge_test,
     >	   use_offshell_rad,
     >	   Egamma_gen_max,
     >	   do_fermi,
     >	   pt_b_param,
     >	   sigc_flag,
     >	   sigc_nbin,
     >	   sigc_kin_min,
     >	   sigc_kin_max,
     >	   sigc_kin_ind,
     >	   using_tgt_field,
     >	   tgt_field_file,
     >     col_flag,
     >	   drift_to_cal

*	E_ARM_ACCEPT
	namelist /E_ARM_ACCEPT/
     >	   SPedge
! 	   SPedge%e%delta%min,
! 	   SPedge%e%delta%max,
! 	   SPedge%e%yptar%min,
! 	   SPedge%e%yptar%max,
! 	   SPedge%e%xptar%min,
! 	   SPedge%e%xptar%max

*	P_ARM_ACCEPT
	namelist /P_ARM_ACCEPT/
     >	   SPedge
!	   SPedge%p%delta%min,
!	   SPedge%p%delta%max,
!	   SPedge%p%yptar%min,
!	   SPedge%p%yptar%max,
!	   SPedge%p%xptar%min,
!	   SPedge%p%xptar%max
! first read the defaults


*	KINEMATICS_MAIN
	namelist /KINEMATICS_MAIN/
     >	   Ebeam,
     >	   dEbeam,
     >	   EXPER,
     >	   electron_arm,
     >	   hadron_arm,
     >	   use_first_cer,
     >     spec

*	BEAM_AND_TARGET_INFO
	namelist /BEAM_AND_TARGET_INFO/
     >      targ,
     >      gen

* parameters to control theory models
        namelist /THEORY_CTRL/theory_par

! the read the user parameters

	data nl_names(1) /"RESTSW"/
	data nl_names(2) /"EXPERIMENT" /
	data nl_names(3) /"DEBUG_PARM"/
	data nl_names(4) /"TARGET"/
	data nl_names(5) /"E_ARM_MAIN"/
	data nl_names(6) /"P_ARM_MAIN"/
	data nl_names(7) /"E_ARM_OFFSET"/
	data nl_names(8) /"P_ARM_OFFSET"/
	data nl_names(9) /"MISC2INT"/
	data nl_names(10) /"SIMULATE"/
	data nl_names(11) /"E_ARM_ACCEPT"/
	data nl_names(12) /"P_ARM_ACCEPT"/
	data nl_names(13) /"KINEMATICS_MAIN"/
	data nl_names(14) /"BEAM_AND_TARGET_INFO"/
	data nl_names(15) /"THEORY_CTRL"/
