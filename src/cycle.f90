module cycle_mod
use real_kind_mod
use global
use mc

implicit none

contains

!--------------------------------------------------------------------------
! Average phase durations have been computed to sum to the average
! cycle time, assuming no checkpoint delays.
! The average phase times are stored in ccp:
! ccp%T_G1, ccp%T_S, ccp%T_G2, ccp%T_M
! Each individual cell is randomly assigned a cycle time, initially and after
! cell division.  The cycle time is a random variate drawn from a lognormal
! distribution. This is done in subroutine set_divide_volume().
! At the same time, a multiplying factor fg is computed.  This is the basic
! factor multiplying the phase times for G1, S and G2 (T_M is fixed.), 
! assuming unconstrained growth, i.e. growth rate = max_growthrate.
! Tgrowth = average growth time = ccp%T_G1 + ccp%T_S + ccp%T_G2
! Tdiv = generated cycle time for the cell (lognormal RV)
! then cp%fg = (Tdiv-ccp%T_M)/Tgrowth
! The rate of volume growth computed for the cell not taking into account
! the cycle time, cp%dVdt, must be divided by this factor cp%fg.
! Immediately before entering a phase, the time to be spent in the phase is
! determined (for example) by: 
! cp%G1_time = ccp%T_G1*max_growthrate/(cp%dVdt/cp%fg))
! and the growth of the cell in a time step dt is given by:
! dV = dt*(cp%dVdt/cp%fg)
! Note that cp%dVdt is computed from the levels of oxygen and glucose:
! cp%dVdt = metab*max_growthrate
! where metab is the product of oxygen and glucose metabolism factors.
! This formulation ensures that the cell grows at a rate that results in
! (approximately) a constant volume = Vdiv at mitosis:
! Vdiv = Vdiv/2 + Tgrowth*max_growthrate,
! Vdiv/2 = Tgrowth*max_growthrate
!--------------------------------------------------------------------------
subroutine log_timestep(cp, ccp, dt, dies)
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: dt
logical :: dies
integer :: phase, ityp, kpar=0
real(REAL_KIND) :: R, tfactor
logical :: switch, S_switch, M_switch

S_switch = .false.
M_switch = .false.
if (cp%dVdt == 0) then
	write(nflog,*) 'dVdt=0, kcell: ',kcell_now,cp%phase
	stop
endif
ityp = cp%celltype
10 continue
phase = cp%phase
if (phase == G1_phase) then
    switch = (tnow > cp%G1_time)
    if (switch) then
        if (.not.is_radiation) then
            cp%phase = S_phase
            goto 10
        endif
        cp%phase = G1_checkpoint
        cp%G1_flag = .true.     ! no G1 checkpoint delay
        cp%G1S_time = tnow      ! no G1 checkpoint
        N_checkpoint = N_checkpoint + 1
        goto 10
    endif
elseif (phase == G1_checkpoint) then  ! this checkpoint combines the release from G1 delay and the G1S repair check
    cp%G1S_flag = (tnow >= cp%G1S_time)
    if (cp%G1_flag .and. cp%G1S_flag) then  ! switch to S-phase
        S_switch = .true.
        cp%phase = S_phase
        tfactor = (max_growthrate(ityp)/cp%dVdt)*cp%fg
		cp%S_time = tnow + tfactor*ccp%T_S     ! restoring previous mode - no arrest
        N_checkpoint = N_checkpoint - 1
        goto 10
    endif
elseif (phase == S_phase) then
    switch = (tnow >= cp%S_time)
    if (switch) then
        if (.not.is_radiation) then
            cp%phase = G2_phase
            goto 10
        endif
        cp%phase = S_checkpoint
        cp%S_flag = .true. ! no fixed checkpoint delay
        cp%SG2_time = tnow + S_checkpoint_time(cp)
        N_checkpoint = N_checkpoint + 1
        goto 10
    endif
elseif (phase == S_checkpoint) then
    cp%SG2_flag = (tnow > cp%SG2_time)
    if (cp%S_flag .and. cp%SG2_flag) then
        cp%phase = G2_phase
        tfactor = (max_growthrate(ityp)/cp%dVdt)*cp%fg
		cp%G2_time = tnow + tfactor*ccp%T_G2
        N_checkpoint = N_checkpoint - 1
        goto 10
    endif
elseif (phase == G2_phase) then
	switch = (tnow > cp%G2_time .and. cp%V > cp%divide_volume) ! to prevent volumes decreasing 
    if (switch) then
        cp%V = cp%divide_volume     ! correct for slight volume discrepancy here, to maintain correct cell volume
        if (.not.is_radiation) then
            cp%phase = M_phase
            goto 10
        endif
        cp%phase = G2_checkpoint
        cp%G2_flag = .false.
        cp%G2M_time = tnow + G2_checkpoint_time(cp)
        N_checkpoint = N_checkpoint + 1
        goto 10
    endif
elseif (phase == G2_checkpoint) then ! this checkpoint combines the release from G2 delay and the G2M repair check
    cp%G2_flag = .true.
    cp%G2M_flag = (tnow > cp%G2M_time)
    if (cp%G2_flag .and. cp%G2M_flag) then  ! switch to M-phase
        M_switch = .true.
        cp%phase = M_phase
        cp%M_time = tnow + ccp%T_M   
        N_checkpoint = N_checkpoint - 1
        goto 10
    endif
elseif (phase == M_phase) then
    ! do nothing - in new_growcells the phase is immediately changed to cp%dividing, and the mitosis timer starts
endif
dies = .false.
call updateRepair(cp, dt)
end subroutine

end module

