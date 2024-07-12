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
! At the same time, multiplying factors fg(:) are computed.  These are the basic
! factors multiplying the phase times for G1, S, G2, and M. 
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
subroutine log_timestep(cp, ccp, dt)
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: dt
real(REAL_KIND) :: tIR, Cdrug, totDSB(2)
integer :: Nwrite, jpp
logical :: switch

!if (cp%dVdt == 0) then
!	write(nflog,*) 'dVdt=0, kcell: ',kcell_now,cp%phase
!	stop
!endif
!if (kcell_now <= 10) write(nflog,'(a,2i4,f6.3)') 'kcell_now, phase, progress: ',kcell_now,cp%phase,cp%progress
tIR = (tnow - t_irradiation)/3600
10 continue
if (cp%phase == G1_phase) then
!    if (single_cell) then
!        write(*,'(a,5e12.3)') 'G1: progress,fp,dt,T_G1, inc: ',cp%progress,cp%fp,dt,ccp%T_G1,cp%fp*dt/ccp%T_G1
!        write(nflog,'(a,5e12.3)') 'G1: progress,fp,dt,T_G1, inc: ',cp%progress,cp%fp,dt,ccp%T_G1,cp%fp*dt/ccp%T_G1
!    endif
    cp%progress = cp%progress + cp%fp*dt/ccp%T_G1
    if (cp%progress >= 1) then
        if (use_G1_stop) then
            ! At start of CP, need to compute CP delay
            call get_CP_delay(cp)
!            if (kcell_now <= 10) write(*,'(a,i6,f8.1)') 'G1 CP_delay: ', kcell_now,cp%CP_delay/3600
            cp%phase = G1_checkpoint
            cp%progress = 0
!            write(*,'(a,i6)') 'G1 -> checkpoint: ',kcell_now
            goto 20
        endif
        cp%phase = S_phase
        cp%progress = 0
        cp%t_S_phase = tnow
        if (single_cell) then
            write(*,'(a,2f6.2)') 'S entry: tnow, N_DSB: ',tnow/3600,sum(cp%DSB(1:3,:))
            write(nflog,'(a,2f6.2)') 'S entry: tnow, N_DSB: ',tnow/3600,sum(cp%DSB(1:3,:))
        endif
!        goto 10    ! could continue in next phase with the remainder of the time step
    endif
elseif (cp%phase == S_phase) then
    cp%progress = cp%progress + cp%fp*dt/ccp%T_S
    if (cp%progress >= 1) then
        if (use_S_stop) then
            ! At start of CP, need to compute CP delay
            call get_CP_delay(cp)
            cp%phase = S_checkpoint
            cp%progress = 0
            goto 20
        endif
        cp%phase = G2_phase
        cp%progress = 0
        cp%t_start_G2 = istep*DELTA_T
        nSdelay = nSdelay + 1   ! only S doesn't use stops
        if (single_cell) then
            do jpp = 1,2
                totDSB(jpp) = sum(cp%DSB(:,jpp))
            enddo
            write(*,'(a,2f6.2)') 'G2 entry: tnow, N_DSB: ',tnow/3600,sum(cp%DSB(1:3,:))
            write(nflog,'(a,4f8.3,4x,3f8.1)') 'G2 entry: tIR,DSB, Nmis: ',tIR,totDSB,sum(totDSB),2*cp%Nmis(1),cp%Nmis(2),2*cp%Nmis(1)+cp%Nmis(2)
            write(nflog,'(a,2f8.3)') 'tnow, t_irradiation: ',tnow/3600, t_irradiation/3600
! Try this
!            cp%ATM_act = 0
!            cp%ATR_act = 0
        endif
    endif
elseif (cp%phase == G2_phase) then
    if (use_Jaiswal) then
        if (is_radiation) then      ! post-IR
            Nwrite = 0.1*3600/DELTA_T
 !           if (single_cell) write(nflog,'(f6.3,2f9.6)') istep*DELTA_T/3600.0,cp%CC_act,cp%dCC_act_dt
 !           if (single_cell .and. istep <= Nwrite) write(nflog,'(f6.3,2f9.6)') istep*DELTA_T/3600.0,cp%CC_act,cp%dCC_act_dt
 !           if (single_cell .and. mod(istep,Nwrite) == 0) write(nflog,'(a,7f8.4)') 'G2_phase: ',istep*DELTA_T/3600.0,cp%CC_act,cp%ATR_act,cp%ATM_act,cp%DSB(2:3),cp%dCC_act_dt
            tIR = (t_simulation - t_irradiation)/3600   ! time since IR, in hours
            switch = (cp%CC_act >= CC_threshold)
        else
            switch = (cp%CC_act >= CC_threshold)
        endif

        if (switch) then
!            write(nflog,*) 'cell enters M: kcell, CC_act: ',kcell_now, cp%CC_act
            cp%phase = M_phase
            cp%progress = 0
            cp%V = cp%divide_volume     ! set volume here, to maintain correct cell volume at cell division
            cp%t_mitosis = tIR
            if (single_cell .or. test_run) then
                write(*,*) 'Reached start of mitosis at (since IR): ',tIR
                write(nflog,*) 'kcell,reached start of mitosis at tIR,CC_act: ',kcell_now,tIR,cp%CC_act
            endif
            !write(*,*) 'stopping on mitosis'
            !write(nflog,*) 'stopping on mitosis: ',kcell_now
            !stop
!            if (istep == 0) write(nflog,'(a,i4,3f8.3)') 'Exit G2: CC_act, threshold, t: ',kcell_now,cp%CC_act,CC_threshold,t_simulation/3600
!            if (cp%generation == 1) then
!                npet = npet+1
!                phase_exit_time_sum = t_simulation + phase_exit_time_sum + dt
!            endif
        endif
    else
        cp%progress = cp%progress + cp%fp*dt/ccp%T_G2
    !    if (kcell_now < 10) write(*,*) 'phase, progress: ',kcell_now, cp%phase,cp%progress
        if (cp%progress >= 1) then
    !        if (kcell_now <= 10) write(nflog,'(a,i6,2e12.3)') 'To M_phase, V, divide_volume: ',kcell_now, cp%V, cp%divide_volume
            if (use_G2_stop) then
                ! At start of CP, need to compute CP delay
                call get_CP_delay(cp)
                cp%phase = G2_checkpoint
                cp%progress = 0
                goto 20
            endif
            cp%phase = M_phase
            cp%progress = 0
            cp%V = cp%divide_volume     ! correct for slight volume discrepancy here, to maintain correct cell volume
        endif
    endif
elseif (cp%phase == M_phase) then
    ! We never get here - in grower() %phase is changed to dividing
    ! do nothing - in new_growcells the phase is immediately changed to cp%dividing, and the mitosis timer starts
elseif (cp%phase == G1_checkpoint) then
    if (cp%CP_delay > 0) then
        cp%progress = cp%progress + dt/cp%CP_delay
    else
        cp%progress = 1
    endif
    if (cp%progress >= 1) then
        cp%phase = S_phase
        cp%progress = 0
!        write(*,'(a,i6)') 'G1 -> S: ',kcell_now
    endif 
elseif (cp%phase == S_checkpoint) then
    cp%progress = cp%progress + dt/cp%CP_delay
    if (cp%progress >= 1) then
        cp%phase = G2_phase
        cp%progress = 0
    endif 
elseif (cp%phase == G2_checkpoint) then
    cp%progress = cp%progress + dt/cp%CP_delay
    if (kcell_now == 1) write(*,'(a,3f10.2)') 'G2_checkpoint, dt,delay,progress: ',dt,cp%CP_delay,cp%progress
    if (cp%progress >= 1) then
        cp%phase = M_phase
        cp%progress = 0
        if (kcell_now == 1) write(*,'(a,i6)') 'G2 -> M: ',kcell_now
    endif 
endif
20  continue
! Note: if cp%phase = dividing, no repair
if (cp%phase < M_phase) then
    call updateRepair(cp, dt)
endif
end subroutine

#if 0
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
subroutine zlog_timestep(cp, ccp, dt, dies)
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
if (kcell_now == 1) write(nflog,'(a,2i4,2f6.3)') 'kcell, phase, progress, fp: ',kcell_now,cp%phase,cp%progress,cp%fp
ityp = cp%celltype
10 continue
phase = cp%phase
if (phase == G1_phase) then
    switch = (tnow > cp%G1_time)    ! G1_time = time to leave G1
    if (switch) then
        if (.not.is_radiation) then     ! before radiation, there is no G1 checkpoint
            cp%phase = S_phase
            tfactor = (max_growthrate(ityp)/cp%dVdt)*cp%fg(G1_phase)
!            if (tfactor > 1.1) then
!                write(*,'(a,i6,3e12.3)') 'tfactor: ',kcell_now,tfactor,(max_growthrate(ityp)/cp%dVdt),cp%fg
!                stop
!            endif
		    cp%S_time = tnow + tfactor*ccp%T_S     ! restoring previous mode - no arrest
!		    if (kcell_now < 10) write(*,'(3e12.3)') kcell_now,tnow,tfactor*ccp%T_S
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
        tfactor = (max_growthrate(ityp)/cp%dVdt)*cp%fg(G1_phase)
		cp%S_time = tnow + tfactor*ccp%T_S     ! restoring previous mode - no arrest
!		if (kcell_now == 1) write(*,'(3e12.3)') kcell_now,tnow,tfactor*ccp%T_S
        N_checkpoint = N_checkpoint - 1
        goto 10
    endif
elseif (phase == S_phase) then
    switch = (tnow >= cp%S_time)    ! S_time = time to leave S
    if (switch) then
        if (.not.is_radiation) then
            cp%phase = G2_phase
            tfactor = (max_growthrate(ityp)/cp%dVdt)*cp%fg(S_phase)
!            if (kcell_now == 1) then
!                write(nfphase,'(a,i4,2f6.3)') 'S_phase entry: V: ',kcell_now, cp%V/Vdivide0, (1 + (ccp%T_G1 + ccp%T_S)/cp%divide_time)/2
!                write(nfphase,'(a,i4,f6.3,2e12.3,2f6.3)') 'tnow, tfactor: ',kcell_now,tnow/3600,max_growthrate(ityp),cp%dVdt,cp%fg,tfactor
!            endif
		    cp%G2_time = tnow + tfactor*ccp%T_G2
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
        tfactor = (max_growthrate(ityp)/cp%dVdt)*cp%fg(S_phase)
		cp%G2_time = tnow + tfactor*ccp%T_G2
        N_checkpoint = N_checkpoint - 1
        goto 10
    endif
elseif (phase == G2_phase) then
	switch = (tnow > cp%G2_time .and. cp%V > cp%divide_volume) ! G2_time = time to leave G2 (V check to prevent volumes decreasing) 
    if (switch) then
        cp%V = cp%divide_volume     ! correct for slight volume discrepancy here, to maintain correct cell volume
        if (.not.is_radiation) then
            cp%phase = M_phase
!            cp%M_time = tnow + ccp%T_M   
            cp%M_time = tnow + cp%mitosis_duration   
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
!        cp%M_time = tnow + ccp%T_M   
        cp%M_time = tnow + cp%mitosis_duration   
        N_checkpoint = N_checkpoint - 1
        goto 10
    endif
elseif (phase == M_phase) then
    ! do nothing - in new_growcells the phase is immediately changed to cp%dividing, and the mitosis timer starts
endif
dies = .false.
call updateRepair(cp, dt)
end subroutine
#endif

end module

