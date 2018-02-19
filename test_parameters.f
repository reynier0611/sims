! test the input parameter file
!---------------------------------------------------------------
	implicit none
	character*80 line, filename, current_name
	logical always
	integer nl_counter
	include 'dbase_namelists.inc'
!       loop indices
	integer i, j
! it is assumed that each namelist occurs only once in the input file
!
! determine which namelist is in the file and determine also the sequence in which
! they have to be read
 
	integer nl_array(nnames)

! get default values first
	call get_defaults

! find the namelist sequence in the input file
	open(io_nl, file = 'test_file.data')
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
	      print *, 'found namelist : ', current_name
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
 999	close(io_nl)
	print *, 'found ', nl_counter, ' namelist entries' 
	if (nl_counter .eq. 0) then
	   write (6,*) 'read_parameters: nothing to read from : '//filename
	   return 
	endif
! now read the namelists according to the sequence in the file
!	re-open
	open(io_nl, file = 'test_file.data')
	do i = 1, nl_counter
	   print *, 'i = ', i, 'nl_arra(i) = ', nl_array(i)
	   select case( nl_array(i) ) ! select which namelist to read
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
	   end select
	enddo
	print *, 'dumping parameters'
	close(io_nl)
	call dump_parameters
	stop 'all done'
	end

	
!----------------------------------------------------------------------
	subroutine dump_parameters
	implicit none
	include 'dbase_namelists.inc'

! open dump file
	open(io_nl, file = 'dump_par.data')
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

	close(io_nl)
	return
	end
!---------------------------------------------------------------

	subroutine get_defaults
	implicit none
	include 'dbase_namelists.inc'

! the rewind statments make sure that the sequence of namelists
! do not matter, but one should only appear once 

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

	close(io_nl)

	if (debug(2)) write(6,*)'user namelist: entering'
	return
	end
