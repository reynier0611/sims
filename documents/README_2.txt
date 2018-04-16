******************************************
              ---  SIMC  ---
	    Changes to the code
       Notes by Reynier Cruz Torres
******************************************
This simc version I got from Werner Boeglin (boeglinw@fiu.edu)
and then implemented the following changes:

***************************
            1
ADDING ANOTHER PARAMETER
TO THE INPUT FILES
***************************
1) In the file: "dbase_namelists.inc"	, added "> doing_bound,"
2) In the file: "simulate.inc"       	, added "logical doing_bound"
					, added "> doing_bound"
3) In the file: "nml_default.data"	, added "doing_bound = F"
This last change maked the default of said parameter be false



***************************
            2
RESTORING ALL THE PIECES
RELATED TO THE SPECTRAL
FUNCTION
***************************
Werner had gotten rid of the parts related to the Benhar Spectral
function. I put everything back. I followed what was being done in
the default version of SIMC. I also made sure that the code
understands the new input parameter: doing_bound and put back the
benhar and transparency input parameters.



***************************
            3
MODIFIED THE FILE
SF_LOOKUP.F
***************************
This file contains the function that reads in and spits out values
from the spectral function external files. It was reading the
spectral function as a whole, when, in reality, the spectral
function contains information for the 2-body breakup (when the
nucleus is split into two components) and the 3-body breakup (when
the nucleus is split into three components) kinematics. These
kinematics are very different, and treating them as equal is a big
error. I modified the function to account for this difference.
doing_bound is the flag in the input file that allows for us to
select the kinematics:

** doing_bound = T ---> 2-body breakup
** doing_bound = F ---> 3-body breakup



***************************
            4
CHANGED LOGICAL STATEMENTS
DOING_DEUTERIUM AND
DOING_HEAVY
***************************
Went through the code and changed (in most cases, but not in all
cases) the flags:
doing_deuterium	-> changed to:	doing_deuterium.or.(doing_heavy.and.doing_bound)
doing_heavy	-> changed to:	doing_heavy.and.(.not.doing_bound)
Now spectral function is separated into 2- and 3-body breakups. 

* For a full list of all the parameters that I changed and those that were
  not changed, see ./0_REYNIER/list.txt



***************************
            5
MINOR CHANGES WITH RESPECT
TO THE PREVIOUS VERSION (4)
***************************
Added a few lines of comments, and organized the code a little bit. No major
changes. The only real change was the step where we "de-integrate" the
spectral function (couldn't find a better way to say that). That is, in the
SF file, the probability values "S" are given already integrated as follows:

S = 4*Pi [s] Pm^2 dPm dEm

where s is the actual ("un-integrated") probability. In other words, we are
given S and we actually need s. Thus, we have to do:

s = S / (4*Pi Pm^2 dPm dEm)

This step was done in a step that was not very convenient. We moved it and
now it is done as soon as the normalization is calculated (Note: this MUST
be done after the normalization is computed, otherwise the results are
simply wrong!) 





***************************
            6
ADDED A NEW BRANCH TO THE
OUTPUT TREE
***************************
I was looking at the in-plane angle of the outgoing proton, and I was doing
the conversion from h_yptar (ssyptar) to "theta" by hand. I decided to add
a new branch that outputs theta in degrees directly. To do this, I had to
implement two small changes:

1) In: NtupleInit.f I added two lines:
	From:
	m = m+1
        NtupleTag(m) = 'SF_weight_recon'      ! 66 spectral function with reconstructed quantities

	To:
	m = m+1
        NtupleTag(m) = 'SF_weight_recon'      ! 66 spectral function with reconstructed quantities
        m = m+1                     ! 67 RCT 8/9/2016 outgoing proton in-plane angle
        NtupleTag(m) = 'h_Thf'      ! 67 RCT 8/9/2016 outgoing proton in-plane angle

2) In: results_write.f I added one line
	From:
        ntu(67) = (spec%p%theta+recon%p%yptar)*180./3.1415926536     !  RCT 8/9/2016 outgoing reconstructed proton in-plane angle 

	To:
        ntu(66) = main%SF_weight                     !  Spectral Function fro recon. quantities
        ntu(67) = (spec%p%theta+recon%p%yptar)*180./3.1415926536     !  RCT 8/9/2016 outgoing reconstructed proton in-plane angle 

I did this following Werner's email:
"Here is (in short) how to add a variable to the output tree. This is for e,e’p only:
In NtupleInit.f add another NtupleTag (the last once currently is SF_weight_recon)
and make sure you increment the counter variable m.
Then in results_write.f you store the corresponding value in the ntu array. The way
it is programmed is a bit dangerous as here you need to give the actual array index
.e.g. 67, in other words you need to increment it yourself adn hardwire it. This could
be changed but I was too lazy to change the original code.
This is basically all you have to do. The branches in the tree are created automatically,
so nothing else needs to be changed."


***************************
            7
REMOVED DUPLICATED CODE
***************************
The missing momentum distribution was giving a yield twice as big as we expected.
I emailed Werner about this and he replied back explaining that he had accidentally
duplicated a piece of code in event.f that effectively multiplied the cross section
times 2. I removed that piece of code from this version.


***************************
            8
ATTEMPTED TO ADD TRITIUM
SPECTRAL FUNCTION
***************************
The tritium proton spectral function is the same as the 3He neutron spectral
function. I added a small piece of code in dbase.f. Here it is:

! RCT 9/13/2016 Added this little piece of code that takes the neutron
! 3He distribution and uses it for the proton distribution in tritium.
          if ((nint(targ%A).eq.3).and.(nint(targ%Z).eq.2)) then
            call sf_lookup_init(tmpfile,.true.)                 !proton S.F.
          else if ((nint(targ%A).eq.3).and.(nint(targ%Z).eq.1)) then
            call sf_lookup_init(tmpfile,.false.)                !neutron S.F.
          endif
! Choos proton or neutron spectral function based on targ.Mtar_struck
!         if (abs(targ%Mtar_struck-Mp).le.1.d-6) then
!           call sf_lookup_init(tmpfile,.true.)                 !proton S.F.
!         else if (abs(targ%Mtar_struck-Mn).le.1.d-6) then
!           call sf_lookup_init(tmpfile,.false.)                !neutron S.F.
!         else
!           write(6,*) 'targ%Mtar_struck = ',targ%Mtar_struck
!           write(6,*) 'targ%Mtar_struck not equal to Mp or Mn'
!           write(6,*) 'DEFAULTING TO PROTON SPECTRAL FUNCTION!!!'
!           call sf_lookup_init(tmpfile,.true.)                 !proton S.F.
!          endif
        endif


***************************
            9
INCORPORATING KAPTARI SF
*************************** 
1) renamed the input parameter use_benhar_sf to use_sf in every single file.
                  
2) added a new input parameter: sf_version
*** In the file: "dbase_namelists.inc"   , added "> sf_version,"
*** In the file: "simulate.inc"          , added "integer*4 sf_version"
                                        , added "> sf_version"
*** In the file: "nml_default.data"      , added "sf_version = 0"
This last change maked the default of said parameter be 0

3) added Kaptari as a version in the code: in dbase.f:

! ============================================================
! RCT 10/26/2016 - using different spectral function versions
!           Benhar spectral function
            if(use_sf.eq.0) then
              tmpfile='benharsf_3mod.dat'
!           Kaptari spectral function
            else if(sf_version.eq.1) then
              tmpfile='kaptarisf_3par.dat'
            else
              write(6,*) 'Wrong value for sf_version. Defaulting to Benhar'
              tmpfile='benharsf_3mod.dat'
            endif
! ============================================================









