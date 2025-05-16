!-----------------------------------------------------------------------------------------
! Main program
!-----------------------------------------------------------------------------------------
PROGRAM monolayer_main
use drm_monolayer_mod
use global
implicit none
integer :: ncpu, res
character*(128) :: infile, outfile, runfile
integer :: status, nlen, cnt, i, inbuflen, outbuflen
integer :: jstep, irun
character*(128) :: b, c, progname
real(8) :: t1, t2
real(8) :: colony_days

real(8) :: progress(30)
integer :: Nph
logical :: use_single

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
!        outfile = c(1:nlen)																! --> outfile
!        write(*,*) 'Output file: ',outfile
        read(c(1:nlen),*) colony_days															! --> ncpu
!        simulate_colony = (colony_days > 0)
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
    outfile = 'drm_monolayer_main.out'
endif

! Synchronisation of cell IR
use_synchronise = .false.
use_single = .true. ! to simulate a cell (or cells) at specified synch_phase and synch_progress
synch_phase = S_phase   !G1 is 1 - 6, S is 7 - 15, G2 is 16 - 19
synch_fraction = 0.5    ! only used if use_single true (also synch_phase)
nph = 1
if (use_synchronise) then
    if (use_single) then
        nph = 1
        progress(1) = synch_fraction
    else
        nph = 10
        do i = 1,nph
            progress(i) = (i-1)*1.0/nph
        enddo
    endif
endif

!do synch_phase = 1,3
do irun = 1,nph
    if (use_synchronise) then
        synch_fraction = progress(irun)    !(irun-1)*0.2
        write(*,*)
    	write(*,'(a,2i4,f6.3)') 'main: irun, synch_phase, synch_fraction: ',irun,synch_phase,synch_fraction
    endif
	inbuflen = len(infile)
	outbuflen = len(outfile)
	write(*,*) 'call execute'
    res = irun
	call execute(ncpu,infile,inbuflen,outfile,outbuflen,res)
	if (res /= 0) stop
	t1 = wtime()
	write(*,*) 'did execute: nsteps, DELTA_T: ',nsteps, DELTA_T
	do jstep = 1,Nsteps+1
		call simulate_step(res)
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
    write(*,*) 'res: ',res
    if (res == 0) write(*,*) 'Exceeded nsteps, not all cells reached mitosis, increase ndays'
	call terminate_run(res)
	t2 = wtime()
	write(*,*) 'time: ',t2-t1
enddo
!enddo
end

