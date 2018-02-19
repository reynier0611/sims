	subroutine sf_lookup_init(infile,protonflag)

! Read in spectral function.  First line has to have # of Pm,Em bins.
! Subsequent lines (number of lines = numPm*numEm) are of the form:
!
! <Pmvalue> <Emvalue> <SF_PROB> <dPm> <dEm>
!
! where Pmvalue/Emvalue are the center of the Pm/Em bins, dPm/dEm are widths,
! and <SF_PROB> is the integrated spectral function within that bin,
! (SF_PROB = SF * dE * P**2 * dP * 4 * pi), and is the probability within
! the Em/Pm bin, IF the spectral function is normalized to 1.
!
	implicit none

	include 'sf_lookup.inc'
	include 'target.inc'
	include 'constants.inc'

	real*8 tmpPm,tmpEm,tmpsfp,tmpsfn,tmpdEm,tmpdPm
	integer*4 iPm,iEm
	character infile*80
	logical protonflag, doing_bound	

	open(unit=55,file=infile,status='old')

	read (55,*) numPm , numEm

! *************************************************************
! Loop over the SF file to get all the values
	do iPm = 1 , numPm
	  do iEm = 1 , numEm
	    read (55,*) tmpPm,tmpEm,tmpsfp,tmpsfn,tmpdPm,tmpdEm

! Choose proton or neutron spectral function.
	    if (protonflag) then
	      sfval(iEm,iPm)=tmpsfp 
	    else
	      sfval(iEm,iPm)=tmpsfn 
	    endif

! Check that Pm value is consistent with others in same bin
	    if (iEm.eq.1) then
	      Pmval(iPm)=tmpPm
	      dPm(iPm)=tmpdPm
	    else
	      if (abs(tmpPm-Pmval(iPm)).gt.0.01) then
	        write(6,*) 'Bad Pm value, line number is approx',(numEm*(iPm-1)+iEm)
	        write(6,*) 'iPm,oldPm,newPm=',iPm,Pmval(iPm),tmpPm
	      endif
	    endif

! Check that Em value is consistent with others in same bin
	    if (iPm.eq.1) then
	      Emval(iEm)=tmpEm
	      dEm(iEm)=tmpdEm
	    else
	      if (abs(tmpEm-Emval(iEm)).gt.0.01) then
	        write(6,*) 'Bad Em value, line number is approx',(numEm*(iPm-1)+iEm)
	        write(6,*) 'iEm,oldEm,newEm=',iEm,Emval(iEm),tmpEm
	      endif
	    endif

	  enddo
	enddo

	sftotnorm = 0.0
	do iPm = 1 , numPm
	  sfnorm(iPm) = 0.0
	  do iEm = 1 , numEm
	    sfnorm(iPm) = sfnorm(iPm) + sfval(iEm,iPm)
	    sftotnorm = sftotnorm + sfval(iEm,iPm)
	  enddo

	enddo
	write(6,*) '------------------------------------------------------'
	write(6,*) 'sf_lookup: Uncorrected S.F. has normalization = ',sftotnorm

	do iPm = 1 , numPm
	  do iEm = 1 , numEm
	    sfval(iEm,iPm) = sfval(iEm,iPm) / sftotnorm
	    sfval(iEm,iPm) = sfval(iEm,iPm) /4./3.1415926535
	    sfval(iEm,iPm) = sfval(iEm,iPm) /(Pmval(iPm)*Pmval(iPm))
	    sfval(iEm,iPm) = sfval(iEm,iPm) /dPm(iPm)/dEm(iEm)
	  enddo
	enddo
	write(6,*) 'sf_lookup: Uncorrected S.F. has been normalized to 1'
	write(6,*) '------------------------------------------------------'
	return
	end

!----------------------------------------------------------------------
! RCT 8/2/2016 This is the subroutine used to get the weighing
! for (for example) 3-Helium

        subroutine sf_lookup_diff(Em,Pm,SFd, doing_bound)
! Get value of fully differential spectral function at Em, Pm.
        implicit none
        real*8 Em, Pm
        real*8 SFd
	logical doing_bound

        real*8 SF
        call sf_lookup(Em, Pm, SF,doing_bound)
        !SFd = SF/4/3.1415926535/(Pm*Pm)/5.0/20.0
	! RCT 8/4/2016 in the line above we were dividing by 4*Pi,
	!Pm**2, dPm, and dEm. This process is done now right when we
	!import the spectral function. Furthermore, dEm and dPm are
	!not assumed to be constant (5.0 and 20.0 above)
        SFd = SF
	return
        end

!#######################################################################################
	subroutine sf_lookup(Em,Pm,SF,doing_bound)

! Get value of spectral function <SF_PROB> at Em, Pm.
!
! Where <SF_PROB> is the integrated spectral function within that bin,
! (SF_PROB = SF * dE * P**2 * dP * 4 * pi), and is the probability within
! the Em/Pm bin, once the spectral function has been normalized to 1.
!

	implicit none
	include 'sf_lookup.inc'

	integer*4 ind,iPm,iEm
	real*8 Em,Pm
	real*8 Em1,Em2,sf1,sf2,logsf,w1,w2
	real*8 SF
	logical doing_bound

! --------------------------------------------------------------
!find nearest Pm value
! ----------------------------------------------
! If Pm is greater than the greatest value of Pm in the SF file
! *** w1 is the weight from the right. Similarly,
! *** w2 is the weight from the left.

	if (Pm.ge.Pmval(numPm)) then
	  iPm=numPm-1
	  w1=0
	  w2=1
! --------------------------------------------------------------
! If Pm is smaller than the smallest value of Pm in the SF file
	else if (Pm.le.Pmval(1)) then
	  iPm=1
	  w1=1
	  w2=0
! --------------------------------------------------------------
! If Pm lies somewhere in the middle of Pm values in the SF file
	else
	  ind=1
	  do while(Pm.gt.Pmval(ind))
	    ind = ind + 1
	  enddo
	  iPm=ind-1
	      w2=(Pm-Pmval(iPm))/(Pmval(iPm+1)-Pmval(iPm))
	      w1=(Pmval(iPm+1)-Pm)/(Pmval(iPm+1)-Pmval(iPm))
	      if (abs(w1*Pmval(iPm)+w2*Pmval(iPm+1)-Pm).gt.0.0001) then
	        write(6,*) 'w1,w2,Pm,Pmval(iPm)=',w1,w2,Pm,Pmval(iPm)
	        stop
	      endif

	endif
! --------------------------------------------------------------
	if (abs(w1+w2-1).gt.0.0001) then
	  write(6,*) 'iPm,Pm,w1+w2=',iPm,Pm,w1+w2
	  stop
	endif

! --------------------------------------------------------------
! RCT 8/2/2016 Only take values that are above the 3-body break up threshold
! The first E miss bin represents the 2-body break up threshold. Therefore, if not
! doing bound should ignore anything at or below the first bin in E miss.

! -----------------------------------------------------
! RCT 8/4/2016 This corresponds to the 3-body break up
	if (.not.doing_bound) then

	  if (Em.gt.Emval(numEm)) then	!linear extrapolation of LOG(f(Em))
	    Em1=Emval(numEm-1)
	    Em2=Emval(numEm)
	    sf1=w1*sfval(numEm-1,iPm)+w2*sfval(numEm-1,iPm+1)
	    sf2=w1*sfval(numEm,iPm)+w2*sfval(numEm,iPm+1)
	  else if ((Em.lt.Emval(2))) then   !The first bin represents 2-body breakup!
	    SF=0
	    return
	  else
	    do iEm=2,numEm-1
	      if (Em.ge.Emval(iEm) .and. Em.lt.Emval(iEm+1)) then
	        Em1=Emval(iEm)
	        Em2=Emval(iEm+1)  
	        sf1 = w1*sfval(  iEm,iPm) + w2*sfval(  iEm,iPm+1)
	        sf2 = w1*sfval(iEm+1,iPm) + w2*sfval(iEm+1,iPm+1)
	      endif
	    enddo
	  endif

	  logsf=(sf1 + (Em-Em1)*(sf2-sf1)/(Em2-Em1))  
	  
	  SF=logsf
	  if (SF.lt.1.d-20) SF=0 !If the obtained value is extremely small, just return zero!

! -----------------------------------------------------
! RCT 8/4/2016 This corresponds to the 2-body break up
	else if (doing_bound) then
	! Return the value corresponding to the first bin!
	  if(w1.ge.w2) then
	    SF = sfval(1,iPm) 
	  else if (w1.lt.w2) then
	    SF = sfval(1,iPm+1)
	  endif
	  if (SF.lt.1.d-20) SF=0 !If the obtained value is extremely small, just return zero!
	endif

	return
	end


!#######################################################################################

	subroutine generate_em(Pm,Em, doing_bound)

! Generate a missing energy for an event with a given Pm.
!
! Have to get spectral function for several Em bins, at the given value
! of Pm.  Then, select Em according to the energy distribution at fixed Pm.
!

	implicit none

	include 'sf_lookup.inc'

	integer*4 ind,iEm
	real*8 Pm,Em			!input Pm value, generated Em value.
	real*8 x(nEmmax),y(nEmmax)
	real*8 ynorm
	real*8 ranprob
	logical doing_bound
	real*8 grnd

! Determine Em distribution at fixed Pm (integrated SF up to Em).
	ynorm = 0.
	x(1) = Emval(1)
	call sf_lookup(x(1),Pm,y(1),doing_bound)

	do iEm = 2 , numEm
	  x(iEm) = Emval(iEm)
	  call sf_lookup(x(iEm),Pm,y(iEm),doing_bound)
	  y(iEm) = y(iEm) + y(iEm-1)
	enddo

! Normalize the Em distribution.
	do iEm = 1 , numEm
	  y(iEm) = y(iEm) / y(numEm)
	enddo

! Generate Em.
	ranprob = grnd()
	ind = 1
	do while (ranprob.gt.y(ind))
	  ind = ind + 1
	enddo
	Em = Emval(ind) + dEm(ind)*(grnd()-0.5)

! CHECK NORMALIZATION, HAVE CODE SAVE EM/PM DISTRIBUTIONS TO CHECK.
! MAKE SURE EM IS NOT BELOW EM_MIN.

!	write(99,'(2f8.2)') Pm,Em

	return
	end
