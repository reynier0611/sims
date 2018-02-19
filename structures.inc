! STRUCTURES.INC 
!   (note that all initializations of fields have been moved
!    into the file "structures_init.inc")

! Define some BASIC record structures and associated parameters

! ... generic cut -->  initialized with MAX large and MIN small
	type cutstype
                sequence
		real*8		min, max
	end type

! ... generic range (rather than a cut) --> initialized with HI small and LO large
	type rangetype
                sequence
		real*8		lo, hi
	end type

! ... generic Cartesian vector
	type Cartesian
                sequence
		real*8		x,y,z
	end type

! ... minimal description of scattered particles, or cuts/ranges on these qties

! declare structures arm and arm2 first
!
! xv, yv, zv are reconstructed vertices
        type arm
                sequence
		real*8		delta, yptar, xptar, z, xv, yv, zv
	end type 		

        type arm2
                sequence
		real*8		delta, yptar, xptar, z, xv, yv, zv
	end type 		


	type double_arm
                sequence
		type(arm)::e
		type(arm2)::p
	end type

	type arm_cuts
             sequence
	     type(cutstype):: delta, yptar, xptar, z
	end type

	type arm_cuts2
             sequence
	     type(cutstype):: delta, yptar, xptar, z
	end type

	type double_arm_cuts
                sequence
		type(arm_cuts)::e
		type (arm_cuts2)::p
	end type

	type arm_range
             sequence
	     type(rangetype):: delta,xptar,yptar,z
        end type

	type arm_range2
             sequence
	     type(rangetype):: delta,xptar,yptar,z
        end type

	type double_arm_range
                sequence
		type(arm_range)::e
		type(arm_range2)::p
	end type

! ... generic focal plane vectors (transport convention, in both spectrometers)

       type arm_FP
           sequence
           real*8		x, dx, y, dy, path
       end type

       type arm_FP2
           sequence
           real*8		x, dx, y, dy, path
       end type

	type double_arm_FP
                sequence
		type(arm_FP)::e
		type(arm_FP2)::p
	end type

! ... full description of a given particle
	type arm_full
                sequence
		real*8	delta,xptar,yptar,z
		real*8	theta,phi,E,P
	end type

	type cuts_true
                sequence
		type(cutstype)::	Em, Pm
	end type

! EVENT structures

! ... description of event -- both actual and reconstructed (measured) versions calculated

	type event
                sequence
		real*8	Ein
		real*8	Em, Pm
		real*8	Emiss, Pmiss
		real*8	Pmx, Pmy, Pmz
		real*8	PmPar, PmPer, PmOop
		real*8	nu, q, Q2, Trec, W, Mrec
		real*8	epsilon, theta_pq, theta_rq, theta_tarq,phi_pq,phi_targ
		real*8  beta, phi_s, phi_c
		real*8  zhad,pt2,xbj
		type (arm_full):: e, p
		type (Cartesian):: ue, up, uq
	end type

! ... target-specific event quantities
	type event_target
                sequence
		real*8	x, y, z, rasterx, rastery
		real*8	teff(3), Eloss(3), Coulomb
	end type

! ... quantities that are determined only once per event
! WB 7/9/2016 added SF_weight_recon
	type event_main
                sequence
		real*8 weight, SF_weight, SF_weight_recon, gen_weight, jacobian, jacobian_corr
		real*8 Ein_shift, Ee_shift
		real*8 sigcc, sigcc_recon, sigcent, sig, sig_recon,coul_corr
                real*8 epsilon,theta_pq,theta_tarq,phi_pq,phi_targ,beta
		real*8 w,t,tmin,q2
                real*8 pcm,thetacm,phicm,wcm
		real*8 davejac,johnjac
		type(event_target)::target
		type(double_arm):: SP, RECON
		type(double_arm_FP)::FP
		real*8 Trec
	end type

! ... a gross structure that serves merely to catch all interesting qties for a 'central' event
        type event_central_rad
                sequence
		real*8 hardcorfac, etatzai, frac(3), lambda(3), bt(2)
		real*8 c_int(0:3), c_ext(0:3), c(0:3), g_int, g_ext, g(0:3)
        end type

	type event_central
                sequence
		real*8		sigcc, nu, q, Q2
		real*8		Em, Pm, W, MM
		type(event_central_rad)::rad
		type (arm_full):: e, p
	end type

! OTHER stuff

! ... spectrometer settings and specifications
      type spec_offset
              sequence
	      real*8	x,y,z,xptar,yptar
      end type

      type spec_offset2
              sequence
	      real*8	x,y,z,xptar,yptar
      end type

       type spectrometer
           sequence
           real*8	P,theta,cos_th,sin_th,phi
	   real*8       sig_musc(2)
	   type(spec_offset)::offset
       end type

       type spectrometer2
           sequence
           real*8	P,theta,cos_th,sin_th,phi
	   real*8       sig_musc(2)
	   type(spec_offset2)::offset
       end type


	type both_spec
                sequence
		type(spectrometer)::e
		type(spectrometer2)::p
	end type

! ... acceptance edges for TRUE and VERTEX quantities, both BEFORE reconstruction
        type edge_arm
           sequence
	   type(cutstype)::		E, yptar, xptar
        end type

        type edge_arm2
           sequence
	   type(cutstype)::		E, yptar, xptar
        end type

	type edge_true
                sequence
		type(edge_arm)::e
		type(edge_arm2)::p
		type (cutstype)::		Em, Pm, Mrec, Trec, Trec_struck
	end type

! ... pieces of the EXP dbase field that we'll need
	type EXP_field
                sequence
		real*8	charge
	end type

! ... generic description of a histogram axis
	type axis
                sequence
		real*8		min,max,bin
		integer		n
	end type

! ... ranges for the quantities that are generated randomly for each event / edges on quantities at the GENERATE stage
       type arm_limits
          sequence
          type (cutstype)::	delta, yptar, xptar, E
       end type

       type arm_limits2
          sequence
          type (cutstype)::	delta, yptar, xptar, E
       end type

	type gen_limits
                sequence
		type(arm_limits)::e
		type(arm_limits2)::p
		type (cutstype)::	sumEgen, Trec
		real*8		xwid, ywid
	end type

! ... ranges of event qties which actually contributed
        type contrib_gen
                  sequence
		  type (arm_range):: e, p
		  type (rangetype):: Trec, sumEgen
        end type

	type contrib_arm
                sequence
          	type (rangetype):: E, yptar, xptar
        end type


	type contrib_arm2
                sequence
          	type (rangetype):: E, yptar, xptar
        end type

        type contrib_true
                  sequence
		  type(contrib_arm)::e
		  type(contrib_arm2)::p
		  type (rangetype):: Em, Pm, Trec
        end type

	type contrib_vertex
            sequence
            type (rangetype):: Trec, Em, Pm
        end type

	type contrib_rad
             sequence
             type (rangetype):: Egamma(3), Egamma_total
        end type

	type contribtype
                sequence
		type(contrib_gen)::gen
		type(contrib_true)::tru
		type(double_arm_range)::SP
		type(contrib_vertex)::vertex
		type(contrib_rad)::rad
	end type

! ... values, and ranges of values which actually contributed, for useful slops (some are local to limits_init)
	type slop_item
                sequence
		real*8		lo, hi, used
	end type
	
	type slop_total
           sequence
           type (slop_item)::	Em, Pm
        end type

	type slop_MC_arm
             sequence	
             type (slop_item):: delta, yptar, xptar
        end type

	type slop_MC_arm2
             sequence
             type (slop_item):: delta, yptar, xptar
        end type

	type slop_MC
                 sequence
                 type(slop_MC_arm)::e
		 type(slop_MC_arm2)::p
        end type

	type sloptype
                sequence
		type(slop_total)::total
		type(slop_MC)::MC
	end type


! ... sum and sum**2 of reconstruction errors (needed to get resolutions)
       type sums_electron
           sequence
           real*8 delta,xptar,yptar,ytar
       end type

       type sums_proton
           sequence
           real*8 delta,xptar,yptar,ytar
       end type

	type sums_twoarm
                sequence
		type(sums_electron)::e
		type(sums_proton)::p
	end type

! ... control of theory models
        type theory_param
                sequence
                character(80) model
                character(80) data_file	        
                real, dimension(10) :: parameter
                logical, dimension(10) :: flags
        end type theory_param

