        implicit none
        include 'dbase_namelists.inc'
        open(15, file = 'nml_default.data')
	read( 15, NML = RESTSW)         !RESTSW       
        rewind(15)
	read( 15, NML = E_ARM_MAIN)     !E_ARM_MAIN   
        rewind(15)
	read( 15, NML = EXPERIMENT)     !EXPERIMENT   
        rewind(15)
	read( 15, NML = DEBUG_PARM)     !DEBUG_PARM   
        rewind(15)
	read( 15, NML = TARGET)         !TARGET       
        rewind(15)
	read( 15, NML = P_ARM_MAIN)     !P_ARM_MAIN   
        rewind(15)
	read( 15, NML = MISC2INT)       !MISC2INT     
        rewind(15)
	read( 15, NML = E_ARM_OFFSET)   !E_ARM_OFFSET 
        rewind(15)
	read( 15, NML = P_ARM_OFFSET)   !P_ARM_OFFSET 
        rewind(15)
	read( 15, NML = SIMULATE)       !SIMULATE     
        rewind(15)
	read( 15, NML = E_ARM_ACCEPT)   !E_ARM_ACCEPT 
        rewind(15)
	read( 15, NML = P_ARM_ACCEPT)   !P_ARM_ACCEPT 

       write (*, NML=RESTSW)
	write( *, NML = EXPERIMENT)
	write( *, NML = DEBUG_PARM)
	write( *, NML = TARGET)
	write( *, NML = E_ARM_MAIN)
	write( *, NML = P_ARM_MAIN)
	write( *, NML = E_ARM_OFFSET)
	write( *, NML = P_ARM_OFFSET)
	write( *, NML = MISC2INT)
	write( *, NML = SIMULATE)
	write( *, NML = E_ARM_ACCEPT)
	write( *, NML = P_ARM_ACCEPT)
        close(15)
        stop
	end
