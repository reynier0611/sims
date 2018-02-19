	subroutine dump_parameters
	implicit none
	include 'dbase_namelists.inc'

! open dump file
	open(io_nl, file = 'nml_dump.data')
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

	close(io_nl)
	return
	end
