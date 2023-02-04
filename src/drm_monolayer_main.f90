subroutine syncher(Np,phase,progress)
real(8) :: progress(30), cyc(4), p, h
integer :: ih, iph, phase(30), Np

cyc(1:4) = [5.36, 9.51, 3.41, 0.43]

ih = 0
h = 0
do
    if (h < cyc(1)) then
        iph = 1
        p = h/cyc(1)
    elseif (h < (cyc(1) + cyc(2))) then
        iph = 2
        p = (h - cyc(1))/cyc(2)
    elseif (h < (cyc(1) + cyc(2) + cyc(3))) then
        iph = 3
        p = (h - cyc(1) - cyc(2))/cyc(3)
    else
        exit
    endif
    phase(ih+1) = iph
    progress(ih+1) = p
    ih = ih+1
    h = h+1
enddo
Np = ih
end subroutine
        
!-----------------------------------------------------------------------------------------
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

real(8) :: progress(30)
integer :: phase(30), Nph

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
    outfile = 'drm_monolayer_main.out'
endif

! Synchronisation of cell IR
use_synchronise = .false.
if (use_synchronise) then
    call syncher(Nph, phase, progress)
    do i = 1,Nph
        write(*,'(2i6,f6.3)') i-1, phase(i), progress(i)
    enddo
endif
synch_phase = S_phase   !G1 is 1 - 6, S is 7 - 15, G2 is 16 - 19
synch_fraction = 0.0
G2_katm3_factor = 1.0
G2_katm4_factor = 1.0
G2_katr3_factor = 1.0
G2_katr4_factor = 1.0

use_fixed_CP = .false.
compute_cycle = .true.
!call get_dimensions(NX,NY,NZ,nsteps,DELTA_T, MAX_CHEMO, cused);
i_hypoxia_cutoff = 3
i_growth_cutoff = 1
do irun = 1,1   ! 1,6   ! 16,19     ! 7,15
!    synch_fraction = progress(irun)    !(irun-1)*0.2
    if (use_synchronise) then
    	write(*,'(a,2i4,f6.3)') 'irun, synch_phase, synch_fraction: ',irun,synch_phase,synch_fraction
    endif
	inbuflen = len(infile)
	outbuflen = len(outfile)
	write(*,*) 'call execute'
!	write(nfrun,*) 'infile: ',infile
!	write(nfrun,*) 'outfile: ',outfile
    res = irun
	call execute(ncpu,infile,inbuflen,outfile,outbuflen,res)
	if (res /= 0) stop
	!call cpu_time(t1)
	t1 = wtime()
	write(*,*) 'did execute: nsteps, DELTA_T: ',nsteps, DELTA_T
	nsumm_interval = (60*60)/DELTA_T   ! number of time steps per hour
	write(*,*) 'nsumm_interval: ',nsumm_interval
	call get_summary(summarydata,i_hypoxia_cutoff,i_growth_cutoff)
	do jstep = 1,Nsteps+1
!		write(*,*) 'jstep: ',jstep
		call simulate_step(res)
		if (mod(jstep,nsumm_interval) == 0) then
			call get_summary(summarydata,i_hypoxia_cutoff,i_growth_cutoff)
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
	call terminate_run(res)
	!call cpu_time(t2)
	t2 = wtime()
	write(*,*) 'time: ',t2-t1
enddo
end

