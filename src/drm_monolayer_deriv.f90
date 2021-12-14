
! Main program
!-----------------------------------------------------------------------------------------
PROGRAM monolayer_main
use drm_monolayer_mod
use global
implicit none
integer :: ncpu, res
real(8) :: summarydata(100)
character*(128) :: infile, outfile, runfile
character*(64) :: travelfile = 'travel_time_dist.out'
integer :: status, nlen, cnt, i, inbuflen, outbuflen
integer :: jstep, hour, ntot, ncog, inflow, irun, i_hypoxia_cutoff,i_growth_cutoff, nsumm_interval
character*(128) :: b, c, progname
real :: vasc
real(8) :: t1, t2
!logical :: simulate_colony
integer :: idist, ndist = 40
real(8) :: PE, colony_days, dist(40), ddist = 50
integer, parameter :: npvary = 12
real(8) :: SF, SF0, dSF, delta_p, dSFave(npvary)
character*(12) :: paramStr(npvary)
integer :: phdist(0:8)
integer :: iprun
real(8) :: param

paramStr = ['baseRate','mitRate','Kaber','Klethal','K_ATM(2)','K_ATM(3)','K_ATM(4)','K_ATR(2)','K_ATR(3)','K_ATR(4)','Pcomplex','PHRsimple']
call disableTCP

call get_command (b, nlen, status)
if (status .ne. 0) then
    write (*,*) 'get_command failed with status = ', status
    stop
end if
!write (nfrun,*) 'command line = ', b(1:nlen)
call get_command_argument (0, c, nlen, status)
if (status .ne. 0) then
    write (*,*) 'Getting command name failed with status = ', status
    stop
end if
!write (*,*) 'command name = ', c(1:len)
progname = c(1:nlen)
cnt = command_argument_count ()
!write (*,*) 'number of command arguments = ', cnt
if (cnt < 3) then
    write(*,*) 'Use: ',trim(progname),' num_cpu input_file colony_days'
    write(*,*) '  simulate colony if colony_days > 0'
    stop
endif

do i = 1, cnt
    call get_command_argument (i, c, nlen, status)
    if (status .ne. 0) then
        write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', i
        stop
    end if
    if (i == 1) then
!        read(c(1:len),'(i)') ncpu
        read(c(1:nlen),*) ncpu															! --> ncpu
        write(*,*) 'Requested threads: ',ncpu
    elseif (i == 2) then
        infile = c(1:nlen)																! --> infile
        write(*,*) 'Input file: ',infile
    elseif (i == 3) then
        read(c(1:nlen),*) colony_days															! --> ncpu
        simulate_colony = (colony_days > 0)
    endif
    use_PEST = .false.
end do
if (cnt == 4) then
    call get_command_argument (4, c, nlen, status)
    PEST_outputfile = c(1:nlen)
    use_PEST = .true.
elseif (cnt /= 3) then
	write(*,*) 'Error: wrong number of arguments'
	stop
endif

if (use_PEST) then
    outfile = PEST_outputfile
else
    outfile = 'drm_monolayer_deriv.out'
endif

!call get_dimensions(NX,NY,NZ,nsteps,DELTA_T, MAX_CHEMO, cused);
i_hypoxia_cutoff = 3
i_growth_cutoff = 1
do iprun = 1,npvary
do irun = 1,2
	write(*,*) 'irun: ',irun
	inbuflen = len(infile)
	outbuflen = len(outfile)
	write(*,*) 'call execute'
!	write(nfrun,*) 'infile: ',infile
!	write(nfrun,*) 'outfile: ',outfile
	call execute(ncpu,infile,inbuflen,outfile,outbuflen,res)
	if (res /= 0) stop
	!call cpu_time(t1)
	t1 = wtime()
	write(*,*) 'did execute: nsteps, DELTA_T: ',nsteps, DELTA_T
	nsumm_interval = (60*60)/DELTA_T   ! number of time steps per hour
	write(*,*) 'nsumm_interval: ',nsumm_interval
!	call get_summary(summarydata,i_hypoxia_cutoff,i_growth_cutoff)
    if (irun == 2) then
    if (iprun == 1) then
        param = baseRate
!        paramStr(iprun) = 'baseRate'
        delta_p = 0.01*param
        baseRate = param + delta_p
    elseif (iprun == 2) then
        param = mitRate
!        paramStr(iprun) = 'mitRat'e
        delta_p = 0.01*param
        mitRate = param + delta_p
    elseif (iprun == 3) then
        param = Kaber
!        paramStr(iprun) = 'Kaber'
        delta_p = 0.01*param
        Kaber = param + delta_p
    elseif (iprun == 4) then
        param = Klethal
!        paramStr(iprun) = 'Klethal'
        delta_p = 0.01*param
        Klethal = param + delta_p
    elseif (iprun == 5) then
        param = K_ATM(2)
        paramStr(iprun) = 'K_ATM(2)'
        delta_p = 0.01*param
        K_ATM(2) = param + delta_p
    elseif (iprun == 6) then
        param = K_ATM(3)
!        paramStr(iprun) = 'K_ATM(3)'
        delta_p = 0.01*param
        K_ATM(3) = param + delta_p
    elseif (iprun == 7) then
        param = K_ATM(4)
!        paramStr(irun-1) = 'K_ATM(4)'
       delta_p = 0.01*param
        K_ATM(4) = param + delta_p
    elseif (iprun == 8) then
        param = K_ATR(2)
!        paramStr(irun-1) = 'K_ATR(2)'
        delta_p = 0.01*param
        K_ATR(2) = param + delta_p
    elseif (iprun == 9) then
        param = K_ATR(3)
!        paramStr(irun-1) = 'K_ATR(3)'
        delta_p = 0.01*param
        K_ATR(3) = param + delta_p
    elseif (iprun == 10) then
        param = K_ATR(4)
!        paramStr(irun-1) = 'K_ATR(4)'
        delta_p = 0.01*param
        K_ATR(4) = param + delta_p
    elseif (iprun == 11) then
        param = Pcomplex
!        paramStr(irun-1) = 'Pcomplex'
        delta_p = 0.01*param
        Pcomplex = param + delta_p
    elseif (iprun == 12) then
        param = PHRsimple
!        paramStr(irun-1) = 'PHRsimple'
        delta_p = 0.01*param
        PHRsimple = param + delta_p
    endif
    endif
	do jstep = 1,Nsteps+1
!		write(*,*) 'jstep: ',jstep
		call simulate_step(res)
		if (mod(jstep,nsumm_interval) == 0) then
!			call get_summary(summarydata,i_hypoxia_cutoff,i_growth_cutoff)
		endif
		if (res == 1) then
		    write(*,*) 'Success'
		    exit
		elseif (res == 2) then
			write(*,*) 'All cells have died'
			stop
		elseif (res /= 0) then
			write(*,*) 'Error exit: ',res
			stop
		endif
	enddo
	if (simulate_colony) then
	    call make_colony_distribution(colony_days,dist,ddist,ndist,PE)
	endif
	call getResults(SF,phdist)
	if (irun == 1) then
	    SF0 = SF
	else
	    dSF = SF - SF0
	    dSFave(iprun) = dSF/SF0
	endif
	call terminate_run(res)
	!call cpu_time(t2)
!	t2 = wtime()
!	write(*,*) 'time: ',t2-t1
    write(*,'(a,2e12.3)') 'SFave: ',SF,log10(SF)
    write(*,'(a,9i6)') 'phase_dist: ',phdist(0:3)
!    if (irun == 2) then
!        write(*,*)
!        write(*,'(a,e12.3)') 'derivative: ',dSF/delta_p
!        write(*,'(a,e12.3)') 'fractional change: ',dSF/SF0
!    endif
    write(*,*)
enddo
enddo
open(10,file='dSFave.out',status='replace')
do iprun = 1,npvary
    write(10,'(i4,2x,a12,e12.3)') iprun,paramStr(iprun),dSFave(iprun)
enddo
close(10)
end

