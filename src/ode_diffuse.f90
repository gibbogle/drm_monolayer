!----------------------------------------------------------------------------------
! Note: The value of spcrad was first determined by writing out the value computed in rkc.
! Later it was just determined by trial, then made into a run parameter.
!----------------------------------------------------------------------------------
double precision function spcrad(neqn,t,y)
!DEC$ ATTRIBUTES DLLEXPORT :: spcrad
use global
integer :: neqn
double precision :: t, y(neqn)
spcrad = spcrad_value
end function

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
module ode_diffuse

use chemokine
use metabolism
use cycle_mod
use rkc_90

implicit none

integer :: ivdbug

!real(REAL_KIND) :: work_rkc(8+5*2*MAX_CHEMO)
real(REAL_KIND) :: work_rkc(8+5*NUTS*(N1D+1))
logical :: chemo_active(2*MAX_CHEMO)    ! flags necessity to solve for the constituent
real(REAL_KIND) :: CO2_rkc				! O2 concentration for f_rkc_drug
integer :: idrug_rkc					! drug number for f_rkc_drug

contains

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine CheckDrugConcs
integer :: ndrugs_present, drug_present(3*MAX_DRUGTYPES), drug_number(3*MAX_DRUGTYPES)
integer :: idrug, iparent, im, kcell, ichemo, i
type(cell_type), pointer :: cp

ndrugs_present = 0
drug_present = 0
do idrug = 1,ndrugs_used
	iparent = DRUG_A + 2*(idrug-1)
	if (chemo(iparent)%present) then		! simulation with this drug has started
	    do im = 0,1
	        ichemo = iparent + im
	        ndrugs_present = ndrugs_present + 1
	        drug_present(ndrugs_present) = ichemo
	        drug_number(ndrugs_present) = idrug
	    enddo
	endif
enddo

do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
    cp => cell_list(kcell)
	do i = 1,ndrugs_present
	    ichemo = drug_present(i)
	    idrug = drug_number(i)
	    if (cp%Cin(ichemo) > Cthreshold) drug_gt_cthreshold(idrug) = .true.
!	    if (cp%Cex(ichemo) > Cthreshold) drug_gt_cthreshold(idrug) = .true.
	enddo
enddo
do i = 1,ndrugs_present
    ichemo = drug_present(i)
    idrug = drug_number(i)
    if (Caverage(MAX_CHEMO + ichemo) > Cthreshold) drug_gt_cthreshold(idrug) = .true.
enddo
end subroutine



!----------------------------------------------------------------------------------
! Cells in checkpoint are excluded from Nmetabolisingcells(?) Check with Bill!
! Note: it is assumed that no cells are dying or dead.
! Note: O2 consumption resulting from drug metabolism is ignored!
!----------------------------------------------------------------------------------
subroutine f_rkc_drug(neqn,t,y,dydt,icase)
integer :: neqn, icase
real(REAL_KIND) :: t, y(neqn), dydt(neqn)
integer :: k, kk, i, ichemo, idrug, iparent, im, ict, Nmetabolisingcells
real(REAL_KIND) :: dCsum, dCdiff, dCreact, vol_cm3, Cex
real(REAL_KIND) :: decay_rate, C, membrane_kin, membrane_kout, membrane_flux, area_factor, n_O2(0:2)
logical :: metabolised(MAX_CELLTYPES,0:2)
real(REAL_KIND) :: metab, cell_flux, dMdt, KmetC, vcell_actual, dC, CO2, A, d, dX, dV, Kd, KdAVX, z1, z2
type(drug_type), pointer :: dp
real(REAL_KIND) :: average_volume = 1.2
logical :: use_average_volume = .true.
logical :: is_metab1

!write(nflog,'(42f8.4)') y(1:neqn)
!do k = 1,neqn
!    dydt(k) = -y(k)*chemo(DRUG_A)%decay_rate
!enddo
!return

ict = icase
CO2 = CO2_rkc
idrug = idrug_rkc
A = well_area
d = total_volume/A
dX = d/N1D
dV = A*dX
if (use_average_volume) then
    vol_cm3 = Vcell_cm3*average_volume	  ! not accounting for cell volume change
    area_factor = (average_volume)**(2./3.)
endif
!Nmetabolisingcells = Ncells - (Ndying(1) + Ndying(2))
Nmetabolisingcells = Ncells - N_checkpoint
iparent = DRUG_A + 2*(idrug-1)
dp => drug(idrug)
metabolised(:,:) = (dp%Kmet0(:,:) > 0)
n_O2(:) = dp%n_O2(ict,:)

k = 0
do im = 0,1
	! First process IC reactions
	k = k+1
	C = y(k)
	Cex = y(k+1)
	ichemo = iparent + im
	Kd = chemo(ichemo)%medium_diff_coef
    decay_rate = chemo(ichemo)%decay_rate
    membrane_kin = chemo(ichemo)%membrane_diff_in
    membrane_kout = chemo(ichemo)%membrane_diff_out
	membrane_flux = area_factor*(membrane_kin*Cex - membrane_kout*C)
	dCreact = 0
    if (im == 0) then
        if (metabolised(ict,0) .and. C > 0) then
		    KmetC = dp%Kmet0(ict,0)*C
		    if (dp%Vmax(ict,0) > 0) then
			    KmetC = KmetC + dp%Vmax(ict,0)*C/(dp%Km(ict,0) + C)
		    endif
		    dCreact = -(1 - dp%C2(ict,0) + dp%C2(ict,0)*dp%KO2(ict,0)**n_O2(0)/(dp%KO2(ict,0)**n_O2(0) + CO2**n_O2(0)))*KmetC
	    endif
!		write(nflog,'(a,3e12.3)') 'dCreact, flux, decay: ',dCreact,membrane_flux/vol_cm3,-C*decay_rate
	    dCreact = dCreact + membrane_flux/vol_cm3
    elseif (im == 1) then	! kk=1 is the PARENT drug
		kk = 1
	    if (metabolised(ict,0) .and. y(kk) > 0) then
		    dCreact = (1 - dp%C2(ict,0) + dp%C2(ict,0)*dp%KO2(ict,0)**n_O2(0)/(dp%KO2(ict,0)**n_O2(0) + CO2**n_O2(0)))*dp%Kmet0(ict,0)*y(kk)
	    endif
	    if (metabolised(ict,1) .and. C > 0) then
		    dCreact = dCreact - (1 - dp%C2(ict,1) + dp%C2(ict,1)*dp%KO2(ict,1)**n_O2(1)/(dp%KO2(ict,1)**n_O2(1) + CO2**n_O2(1)))*dp%Kmet0(ict,1)*C
	    endif
	    dCreact = dCreact + membrane_flux/vol_cm3
!    elseif (im == 2) then	! kk=N1D+2 is the METAB1
!		kk = N1D+2
!	    if (metabolised(ict,1) .and. y(kk) > 0) then
!		    dCreact = (1 - dp%C2(ict,1) + dp%C2(ict,1)*dp%KO2(ict,1)**n_O2(1)/(dp%KO2(ict,1)**n_O2(1) + CO2**n_O2(1)))*dp%Kmet0(ict,1)*y(kk)
!	    endif
!	    if (metabolised(ict,2) .and. C > 0) then
!		    dCreact = dCreact - (1 - dp%C2(ict,2) + dp%C2(ict,2)*dp%KO2(ict,2)**n_O2(2)/(dp%KO2(ict,2)**n_O2(2) + CO2**n_O2(2)))*dp%Kmet0(ict,2)*C
!	    endif
!	    dCreact = dCreact + membrane_flux/vol_cm3
    endif
	dydt(k) = dCreact - C*decay_rate
!	if (im == 0) write(nflog,'(a,f10.6,4e12.3)') 'IC: C,flux,dCreact,C*decay,dydt: ',C,membrane_flux,dCreact,C*decay_rate,dydt(k)
!	write(nflog,'(a,i4,e12.3)') 'dydt: ',im,dydt(k)
	if (isnan(dydt(k))) then
		write(nflog,*) 'f_rkc_drug: dydt isnan: ',im,dydt(k)
		write(*,*) 'f_rkc_drug: dydt isnan: ',im,dydt(k)
		stop
	endif
	
	! Next process grid cell next to the cell layer - note that membrane _flux has already been computed
	k = k+1
	C = y(k)
	z1 = (-Nmetabolisingcells*membrane_flux - Kd*A*(C - y(k+1))/dX)/dV
	dydt(k) = z1 - C*decay_rate
!	write(*,'(a,f10.6,3e12.3)') 'EC: C,z1,C*decay_rate,dydt: ',C,z1,C*decay_rate,dydt(k)
!    if (im == 0) write(nflog,'(i6,6e12.3)') Nmetabolisingcells, membrane_flux, Kd*A, C, y(k+1), (C - y(k+1)),dydt(k)
	
	! Next compute diffusion and decay on the FD grid
	KdAVX = Kd*A/(dV*dX)
	do i = 2,N1D
		k = k+1
		C = y(k)
		if (i < N1D) then
			dydt(k) = KdAVX*(y(k+1) - 2*C + y(k-1)) - C*decay_rate
		else
			dydt(k) = KdAVX*(-C + y(k-1)) - C*decay_rate
		endif
	enddo
enddo
end subroutine


!----------------------------------------------------------------------------------
! This version assumes a single metabolism solution for all cells (all phases)
! Cells in checkpoint are excluded from Nmetabolisingcells.
! Note: it is assumed that no cells are dying or dead.
!----------------------------------------------------------------------------------
subroutine f_rkc_OG(neqn,t,y,dydt,icase)
integer :: neqn, icase
real(REAL_KIND) :: t, y(neqn), dydt(neqn)
integer :: k, kk, i, ichemo, ict, Nmetabolisingcells
real(REAL_KIND) :: dCsum, dCdiff, dCreact, vol_cm3, Cex, Cin(NUTS+1), C_GlnEx
real(REAL_KIND) :: C, membrane_kin, membrane_kout, membrane_flux, area_factor, Cbnd
real(REAL_KIND) :: A, d, dX, dV, Kd, KdAVX  !, K1, K2
type(metabolism_type), pointer :: mp
type(cell_type), pointer :: cp
real(REAL_KIND) :: average_volume = 1.2
logical :: use_average_volume = .true.
integer :: res

!cp => cell_list(1)
cp => master_cell
mp => cp%metab
knt = knt+1
ict = icase
A = well_area
d = total_volume/A
dX = d/N1D
dV = A*dX
if (use_average_volume) then
    vol_cm3 = Vcell_cm3*average_volume	  ! not accounting for cell volume change
    area_factor = (average_volume)**(2./3.)
endif
!Nmetabolisingcells = Ncells - (Ndying(1) + Ndying(2))
Nmetabolisingcells = Ncells - N_checkpoint
do ichemo = 1,NUTS
    Cin(ichemo) = y((ichemo-1)*(N1D+1) + 1)
enddo
C_GlnEx = 0     ! not used
if (noSS) then
    Cin(NUTS+1) = y(NUTS*(N1D+1) + 1)
!    K1 = K_PL
!    K2 = K_LP
endif
!write(nflog,*)
!write(nflog,'(a,i4,f10.3,4e15.6)') 'f_rkc_OGL: knt, t, Cin: ',knt,t,Cin(1:4) 
if (knt > 10000) then
    write(nflog,*) 'ERROR: knt > 10000'
    write(*,*) 'ERROR: knt > 10000'
    stop
endif
call get_metab_rates(cp,Cin,C_GlnEx,res)     ! needs to be from y()
if (res /= 0) then
    write(nflog,*) 'Error: get_metab_rates: res: ',res
    stop
endif
k = 0
do ichemo = 1,NUTS     
	! First process IC reactions
	k = k+1
	C = y(k)
	Cex = y(k+1)
	Kd = chemo(ichemo)%medium_diff_coef
    membrane_kin = chemo(ichemo)%membrane_diff_in
    membrane_kout = chemo(ichemo)%membrane_diff_out
	membrane_flux = area_factor*(membrane_kin*Cex - membrane_kout*C)
    if (ichemo == OXYGEN) then
		dCreact = (-mp%O_rate + membrane_flux)/vol_cm3		! O_rate is rate of consumption
!		if (istep > 285 .and. istep < 290) then
!        write(nflog,'(a,2i5,2f8.5,3e14.6)') 'knt,ichemo,Cex,Cin,membrane_flux,r_O,dydt: ',knt,ichemo,Cex,C,membrane_flux,mp%O_rate,dCreact
!        endif
    elseif (ichemo == GLUCOSE) then	! 
		dCreact = (-mp%G_rate + membrane_flux)/vol_cm3		! G_rate is rate of consumption
    endif
!    write(nflog,'(a,2i5,2f8.5,2e14.6)') 'knt,ichemo,Cex,Cin,membrane_flux,dydt: ',knt,ichemo,Cex,C,membrane_flux,dCreact
	dydt(k) = dCreact
	if (isnan(dydt(k))) then
		write(nflog,'(a,i4,4e12.3)') 'f_rkc_OGL: ichemo, Cex, C, membrane_flux,dydt isnan: ',ichemo,Cex,C,membrane_flux,dydt(k)
		write(*,'(a,i4,4e12.3)') 'f_rkc_OGL: ichemo, Cex,C membrane_flux,dydt isnan: ',ichemo,Cex,C,membrane_flux,dydt(k)
		if (ichemo == OXYGEN) then
		    write(nflog,*) 'mp%O_rate: ',mp%O_rate
		    write(*,*) 'mp%O_rate: ',mp%O_rate
		endif
		stop
	endif
	
	! Next process grid cell next to the cell layer - note that membrane _flux has already been computed
	k = k+1
	C = y(k)
	dydt(k) = (-Nmetabolisingcells*membrane_flux - Kd*A*(C - y(k+1))/dX)/dV
!	write(*,'(a,i4,e12.3)') 'k,dydt: ',k,dydt(k)
	
	! Next compute diffusion and decay on the FD grid
	! Need special treatment for oxygen at air boundary
	KdAVX = Kd*A/(dV*dX)
	do i = 2,N1D
		k = k+1
		C = y(k)
		if (i < N1D) then
			dydt(k) = KdAVX*(y(k+1) - 2*C + y(k-1))
		else
			if (ichemo == OXYGEN) then
				Cbnd = chemo(OXYGEN)%bdry_conc
				dydt(k) = KdAVX*(Cbnd - 2*C + y(k-1))
			else
				dydt(k) = KdAVX*(-C + y(k-1))
			endif
		endif
	enddo
enddo
end subroutine

!----------------------------------------------------------------------------------
! Caverage(1:MAX_CHEMO) is the cell concentrations (all cells the same) = IC
! Caverage(MAX_CHEMO+1:2*MAX_CHEMO) is the concentrations close to the bottom = EC
! CmediumAve() is the average concentrations in the medium
!----------------------------------------------------------------------------------
subroutine Solver(it,tstart,dt,nc,ok)
integer :: it, nc
real(REAL_KIND) :: tstart, dt
integer :: kcell
logical :: ok
logical :: use_drugsolver = .true.

ok = .true.
call OGSolver(tstart,dt,ok)
if (.not.use_drugsolver) return
!if (DRUG_A_inhibiter) then
!    if (use_inhibiter) then     ! fixed drug concentration
!        do kcell = 1,ncells
!            cell_list(kcell)%Cin(DRUG_A) = event(1)%conc
!        enddo
!    endif
!    return
!endif
if (chemo(DRUG_A)%present) then
	call DrugSolver(DRUG_A,tstart,dt,1,ok)
endif
!if (chemo(DRUG_B)%present) then
!	call DrugSolver(DRUG_B,tstart,dt,2,ok)
!endif
end subroutine

!----------------------------------------------------------------------------------
! For phase-dependent drug, e.g. EDU, calls DrugPhaseSolver
!----------------------------------------------------------------------------------
subroutine DrugSolver(iparent,tstart,dt,idrug,ok)
integer :: iparent, idrug
real(REAL_KIND) :: tstart, dt
logical :: ok
integer :: ichemo, k, ict, neqn, i, kcell, im
real(REAL_KIND) :: t, tend
real(REAL_KIND) :: C(3*N1D+3), Csum
real(REAL_KIND) :: timer1, timer2
! Variables for RKC
integer :: info(4), idid
real(REAL_KIND) :: rtol, atol(1)
type(rkc_comm) :: comm_rkc(1)

if (drug(idrug)%phase_dependent) then
	call DrugPhaseSolver(iparent,tstart,dt,idrug,ok)
	return
endif

!write(nflog,*) 'DrugSolver: ',istep
ict = 1 ! for now just a single cell type
idrug_rkc = idrug
CO2_rkc = Caverage(OXYGEN)

k = 0
do im = 0,1
	ichemo = iparent + im
	if (.not.chemo(ichemo)%present) cycle
	k = k+1
	C(k) = Caverage(ichemo)		! IC 
	do i = 1,N1D
		k = k+1
!		C(k) = Cdrug(im,i)		! EC
		C(k) = chemo(ichemo)%Cmedium(i)
		if (im == 1 .and. C(k) > 1.0) then      ! bad value
		    write(*,'(a,2i4,e12.3)') 'bad chemo(ichemo)%Cmedium(i): ichemo,i,C: ',ichemo,i,C
		    write(nflog,'(a,2i4,e12.3)') 'bad chemo(ichemo)%Cmedium(i): ichemo,i,C: ',ichemo,i,C
		    stop
		endif
	enddo
enddo
!write(*,'(a,5f12.8)') 'pre DrugSolver: EC: ',(chemo(iparent+i)%Cmedium(1),i=0,1)

neqn = k
!write(nflog,*) 'neqn: ',neqn
info(1) = 1
info(2) = 1		! = 1 => use spcrad() to estimate spectral radius, != 1 => let rkc do it
info(3) = 1
info(4) = 0
rtol = 1d-5
atol = rtol

idid = 0
t = tstart
tend = t + dt
call rkc(comm_rkc(1),neqn,f_rkc_drug,C,t,tend,rtol,atol,info,work_rkc,idid,ict)
if (idid /= 1) then
	write(logmsg,*) 'Solver: Failed at t = ',t,' with idid = ',idid
	call logger(logmsg)
	ok = .false.
	return
endif
!write(*,'(a,3e12.3)') 'IC: ',C(1),C(N1D+2)

! This determines average cell concentrations, assumed the same for all cells
! Now put the concentrations into the cells 

k = 0
do im = 0,1
	ichemo = iparent + im
	if (.not.chemo(ichemo)%present) cycle
	k = k+1
	Caverage(ichemo) = C(k)		! IC 
    Caverage(MAX_CHEMO + ichemo) = C(k+1)	! not really average, this is medium at the cell layer, i.e. EC
	Csum = 0
	do i = 1,N1D
		k = k+1
		chemo(ichemo)%Cmedium(i) = C(k)
		Csum = Csum + C(k)
	enddo
	Cmediumave(ichemo) = Csum/N1D
    do kcell = 1,nlist
        if (cell_list(kcell)%state == DEAD) cycle
        cell_list(kcell)%Cin(ichemo) = Caverage(ichemo)
    enddo
enddo
!write(*,'(a,5f12.8)') 'EC: ',(Caverage(MAX_CHEMO+iparent+i),i=0,1)

!do im = 0,1
!    ichemo = iparent + im
!	if (.not.chemo(ichemo)%present) cycle
!    k = im*(N1D+1) + 1
!    Caverage(ichemo) = C(k)
!	Csum = 0
!    do i = 1,N1D
!		Csum = Csum + C(k+i)
!	enddo
!	Cmediumave(ichemo) = Csum/N1D
!    do kcell = 1,nlist
!        if (cell_list(kcell)%state == DEAD) cycle
!        cell_list(kcell)%Cin(ichemo) = Caverage(ichemo)
!    enddo
!    Caverage(MAX_CHEMO + ichemo) = C(k+1)	! not really average, this is medium at the cell layer, i.e. EC
!!	write(nflog,'(a,i3,5e12.3)') 'Cdrug: im: ',im,Cdrug(im,1:5)
!enddo
!write(*,'(a,3e12.3)') 'Cell drug conc: ',(Caverage(DRUG_A+k),k=0,2)

end subroutine


!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine OGSolver(tstart,dt,ok)
real(REAL_KIND) :: tstart, dt
logical :: ok
integer :: ichemo, k, ict, neqn, i, kcell, it, res, k0
real(REAL_KIND) :: t, tend
real(REAL_KIND) :: C(NUTS*(N1D+1)+1), Cin(NUTS+1), Csum, dCdt(NUTS*(N1D+1)+1), C_P, dtt
real(REAL_KIND) :: timer1, timer2
! Variables for RKC
integer :: info(4), idid
real(REAL_KIND) :: rtol, atol(1)
type(rkc_comm) :: comm_rkc(1)
type(metabolism_type), pointer :: mp
type(cell_type), pointer :: cp, cpm
real(REAL_KIND) :: Cic,Cex,area_factor,membrane_kin,membrane_kout,membrane_flux
integer :: nt = 1000
logical :: use_explicit = .false.		! The explicit approach is hopelessly unstable, even with nt = 1000
! Checking
real(REAL_KIND) :: dC_Pdt, vol_cm3, K1, K2
real(REAL_KIND) :: average_volume = 1.2

!write(nflog,*)
!write(nflog,*) 'OGSolver: ',istep
ict = selected_celltype ! for now just a single cell type 
cpm => master_cell
mp => cpm%metab
mp%recalcable = -1     ! This ensures full solution procedure at the start of each time step

k = 0
do ichemo = 1,NUTS
    Caverage(ichemo) = max(0.0,Caverage(ichemo))    ! try adding this to prevent -ve C_Gln
	k = k+1
	C(k) = Caverage(ichemo)		! IC 
!	write(nflog,'(a,2i4,e12.3)') 'ichemo,k,C(k): ',ichemo,k,C(k)
	k0 = k
	do i = 1,N1D
		k = k+1
		C(k) = C_OGL(ichemo,i)	! EC
	enddo
!	write(nflog,'(10e12.3)') C(k0+1:k0+N1D)
enddo

! Locations of Cin and Cex: 
! Cin = C(ichemo + (ichemo-1)*N1D) = Caverage(ichemo)
! Cex = C_OGL(ichemo,1)
!
!write(nflog,'(10e12.3)') C(1:N1D+1)
!write(nflog,'(a,7f8.3)') 'OGLsolver: Cin, tstart, dt (h): ',Caverage(1:NUTS),tstart/3600,dt/3600
!write(nflog,'(a,f6.3)') 'OGLsolver: glutamine: IC: ',C(3*(N1D+1) + 1)
!write(nflog,'(10f6.3)') C(3*(N1D+1)+2: 3*(N1D+1)+N1D+1)
!C_P = 0
!if (noSS) then
!    k = k+1
!    C_P = mp%C_P
!    C(k) = C_P
!endif
neqn = k

info(1) = 1
info(2) = 1		! = 1 => use spcrad() to estimate spectral radius, != 1 => let rkc do it
info(3) = 1
info(4) = 0
rtol = 5d-4		! was 5d-4
atol = rtol
idid = 0
t = tstart
tend = t + dt
knt = 0
call rkc(comm_rkc(1),neqn,f_rkc_OG,C,t,tend,rtol,atol,info,work_rkc,idid,ict)
if (idid /= 1) then
	write(logmsg,*) 'Solver: Failed at t = ',t,' with idid = ',idid
	call logger(logmsg)
	ok = .false.
	return
endif

! This determines average cell concentrations, assumed the same for all cells
! Now put the concentrations into the cells 
do ichemo = 1,NUTS
	if (.not.chemo(ichemo)%present) cycle
    k = (ichemo-1)*(N1D+1) + 1
    C(k) = max(0.0,C(k))
    Caverage(ichemo) = C(k)
    
    ! Try smoothing
    Caverage(ichemo) = (Caverage(ichemo) + master_cell%Cin(ichemo))/2
    
    Csum = 0
    do i = 1,N1D
        C(k+i) = max(0.0,C(k+i))
		C_OGL(ichemo,i) = C(k+i)
		Csum = Csum + C(k+i)
		chemo(ichemo)%Cmedium(i) = C(k+i)
	enddo
	Cmediumave(ichemo) = Csum/N1D
    do kcell = 1,nlist
        if (cell_list(kcell)%state == DEAD) cycle
        cell_list(kcell)%Cin(ichemo) = Caverage(ichemo)
    enddo
    master_cell%Cin(ichemo) = Caverage(ichemo)
    Caverage(MAX_CHEMO + ichemo) = C(k+1)	! not really average, this is medium at the cell layer, i.e. EC
                                            ! = chemo(ichemo)%Cmedium(1)
!	write(nflog,'(a,i3,5e12.3)') 'Cdrug: im: ',im,Cdrug(im,1:5)
enddo
do kcell = 1,nlist
    cp => cell_list(kcell)
    if (cp%state == DEAD) cycle
enddo
do ichemo = 1,NUTS
    Cin(ichemo) = C((ichemo-1)*(N1D+1) + 1)
enddo

! Update C_A in phase_metabolic(1)
!call update_C_A(dt,mp)

do kcell = 1,nlist
	cp => cell_list(kcell)
!    if (cp%state == DEAD .or. cp%state == DYING) cycle
    if (cp%state == DEAD) cycle
    ! First back up cell metabolism parameters that we need to preserve
    cp%metab = master_cell%metab
enddo
!write(nflog,'(a,5f10.6)') 'master_cell%Cin: ',master_cell%Cin(1:NUTS)
!write(*,'(a,2e12.3)') 'did OGLSolver: Grate, Orate: ',cell_list(1)%metab%G_rate, cell_list(1)%metab%O_rate
return

end subroutine


!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
function smoothstep(x) result(f)
real(REAL_KIND) :: x, f
f = x*x*(3 - 2*x)
end function

!----------------------------------------------------------------------------------
! Note: This computes a rate of change of concentration! mM/s
! Currently only for O2!!! 
! There are two options: use_Cex_Cin = true/false
!
! use_Cex_Cin = true
! ------------------
! The idea is that the speed of the intracellular reactions, compared with other
! processes, is so fast that effectively the intracellular concentration is always
! in equilibrium with the extracellular value.  This means that the rate of consumption
! in the cell matches the rate of transport across the cell membrane: both these rates 
! depend on Cin, therefore we can solve for Cin given Cex then deduce uptake rate
!
! use_Cex_Cin = false
! -------------------
! In this case we just use Cin = Cex to calculate the consumption rate - no
! dependence on chemo(OXYGEN)%membrane_diff
!----------------------------------------------------------------------------------
real(REAL_KIND) function UptakeRate(ichemo,Cex)
integer :: ichemo
real(REAL_KIND) :: Cex
real(REAL_KIND) :: vol, K1, Cin, flux
integer :: n, i

if (ichemo == OXYGEN) then
!    vol = Vsite_cm3
!    vol = Vsite_cm3 - Vextra_cm3	! this was used in the RKC solution
    vol = Vextra_cm3	! the current extracellular volume should be used I think !!!!!!!!!!!!!!!
	if (use_Cex_Cin) then
		Cin = getCin(ichemo,Cex)
!		flux = chemo(ichemo)%membrane_diff*(Cex - Cin)
		flux = (chemo(ichemo)%membrane_diff_in*Cex - chemo(ichemo)%membrane_diff_out*Cin)
	else	! 
		flux = O2_metab(Cex)*chemo(ichemo)%max_cell_rate
	endif
	if (dbug) write(nfout,'(a,2e12.4)') 'Cex, flux: ',Cex,flux
	UptakeRate = flux/vol	! concentration rate (mM/s)
else
	write(logmsg,*) 'ERROR: UptakeRate: currently only for OXYGEN'
	call logger(logmsg)
	stop
endif
end function

!----------------------------------------------------------------------------------
! Computes intracellular O2 concentration as a function of the extracellular level C,
! assuming equilibrium, i.e. rate of consumption = rate of membrane transport.
! Note that the cell's O2 uptake rate is taken to be independent of any other factors,
! e.g. independent of cell size.
! NOTE: Currently only for OXYGEN and GLUCOSE - OK because membrane_diff_in = membrane_diff_out
! Note: needs to be amended to account for HIF-1
!----------------------------------------------------------------------------------
!real(REAL_KIND) function getCinO2(C)
real(REAL_KIND) function getCin(ichemo,C)
integer :: ichemo
real(REAL_KIND) :: C
real(REAL_KIND) :: K1, K2, K2K1, C0, a, b, cc, D, r(3), Cin
integer :: i, n

if (ichemo >= DRUG_A) then
	write(logmsg,*) 'ERROR: getCin: currently only for OXYGEN, GLUCOSE, LACTATE'
	call logger(logmsg)
	stop
endif
!ichemo = OXYGEN
!K1 = chemo(OXYGEN)%membrane_diff*(Vsite_cm3 - Vextra_cm3)
K1 = chemo(ichemo)%membrane_diff_in
K2 = chemo(ichemo)%max_cell_rate
K2K1 = K2/K1
C0 = chemo(ichemo)%MM_C0
if (chemo(ichemo)%Hill_N == 2) then
	a = K2K1 - C
	b = C0*C0
	cc = -b*C
	call cubic_roots(a,b,cc,r,n)
	if (n == 1) then
		Cin = r(1)
	else
		n = 0
		do i = 1,3
			if (r(i) > 0) then
				n = n+1
				Cin = r(i)
			endif
		enddo
		if (n > 1) then
			write(nflog,*) 'getCin: two roots > 0: ',r
			stop
		endif
	endif
elseif (chemo(ichemo)%Hill_N == 1) then
	b = K2K1 + C0 - C
	cc = -C0*C
	D = sqrt(b*b - 4*cc)
	Cin = (D - b)/2
endif
getCin = Cin
end function





!-----------------------------------------------------------------------------------------
! Drug reactions and fluxes are solved for all cells separately.
! Currently only set up for labelling drugs like EDU, for which only the parent can
! exist in free form in the cell, and metabolite "concentration" in actuality
! represents metabolite that has been incorporated into DNA.
! For now only a single metabolite is simulated.
! Note: ichemo = iparent = parent drug
!-----------------------------------------------------------------------------------------
subroutine DrugPhaseSolver(ichemo,tstart,dt,idrug,ok)
integer :: ichemo, idrug
real(REAL_KIND) :: tstart, dt
logical :: ok
logical :: tagged, active
type(cell_type), pointer :: cp
type(drug_type), pointer :: dp
integer :: ict, n_O2, kcell, it, k, i, n_S_phase, n
real(REAL_KIND) :: dtt, decay_rate, membrane_kin, membrane_kout, membrane_flux, Cex, Cex0, vol_cm3, area_factor, R
real(REAL_KIND) :: CO2, C, Clabel, KmetC, dCreact, totalflux, F(N1D+1), A, d, dX, dV, Kd, t, PI_factor
real(REAL_KIND) :: average_volume = 1.2
real(REAL_KIND), dimension(:), pointer :: Cmedium
logical :: use_average_volume = .false.
integer :: nt = 20
integer :: ndt = 20
integer :: kpar = 0
real(REAL_KIND) :: cov = 0.002

Cex = Caverage(MAX_CHEMO+ichemo)
Cex0 = Cex
!if (Cex == 0 .and. chemo(ichemo)%present) then	! stop processing when the parent drug is removed 
!	Caverage(ichemo) = 0
!	chemo(ichemo)%present = .false.
!	return
!endif

dtt = (dt/nt)
dp => drug(idrug)
n_O2 = dp%n_O2(ict,0)
decay_rate = chemo(ichemo)%decay_rate
membrane_kin = chemo(ichemo)%membrane_diff_in
membrane_kout = chemo(ichemo)%membrane_diff_out
Cmedium => chemo(ichemo)%Cmedium
if (use_average_volume) then
    vol_cm3 = Vcell_cm3*average_volume	  ! not accounting for cell volume change 
    area_factor = (average_volume)**(2./3.)
endif

do it = 1,nt
	t = tstart + (it-1)*dtt
    ! Solve for each cell separately
    n_S_phase = 0
    totalflux = 0
    do kcell = 1,nlist
   	    cp => cell_list(kcell)
	    if (cp%state == DEAD) cycle
	    if (.not.use_average_volume) then
		    vol_cm3 = cp%V
		    area_factor = (vol_cm3/Vcell_cm3)**(2./3.)
	    endif
	    ict = cp%celltype
	    CO2 = cp%Cin(OXYGEN)
	    active = drug(idrug)%active_phase(cp%phase)
	    if (active) then	! .and. .not.tagged) then
		    n_S_phase = n_S_phase + 1
	    endif
    !	do it = 1,nt
!	    cellfluxsum = 0     ! Is this OK?????????
	    C = cp%Cin(ichemo)
	    Clabel = cp%Cin(ichemo+1)
	    membrane_flux = area_factor*(membrane_kin*Cex - membrane_kout*C)
	    if (active) then	! .and. .not.tagged) then
		    KmetC = dp%Kmet0(ict,0)*C
		    if (dp%Vmax(ict,0) > 0) then
			    KmetC = KmetC + dp%Vmax(ict,0)*C/(dp%Km(ict,0) + C)
		    endif
		    dCreact = -(1 - dp%C2(ict,0) + dp%C2(ict,0)*dp%KO2(ict,0)**n_O2/(dp%KO2(ict,0)**n_O2 + CO2**n_O2))*KmetC
		    if (trim(dp%name) == 'EDU') then
			    dCreact = dCreact*cp%dVdt/max_growthrate(ict)
		    endif
		    if (trim(dp%name) == 'PI') then
			    if (cp%phase < S_phase) then
				    PI_factor = 1
			    elseif (cp%phase > S_phase) then
				    PI_factor = 2
			    else
!					PI_factor = 1 + (t - cp%S_start_time)/(cp%S_time - cp%S_start_time)
                    PI_factor = 1 + cp%S_time/cp%S_duration
			    endif
!				dCreact = dCreact*cp%dVdt/max_growthrate(ict)
			    dCreact = 0.1*PI_factor*dCreact
		    endif
		    cp%dCdt(ichemo) = dCreact + membrane_flux/vol_cm3 - C*decay_rate
		    cp%dCdt(ichemo+1) = -dCreact
		    Clabel = Clabel + dtt*cp%dCdt(ichemo+1)
	    else
		    cp%dCdt(ichemo) = membrane_flux/vol_cm3 - C*decay_rate	
	    endif
!		R = par_uni(kpar)
	    C = C + dtt*cp%dCdt(ichemo)	!*(1 + (R-0.5)*cov)
!	    cellfluxsum = cellfluxsum + membrane_flux
    !	enddo
	    cp%Cin(ichemo) = C
	    cp%Cin(ichemo+1) = Clabel*(1 + (par_uni(kpar)-0.5)*cov)
!        cp%dMdt(ichemo) = -cellfluxsum/nt	! average flux of parent drug NOT USED ANYWHERE
!        totalflux = totalflux + cp%dMdt(ichemo)
        totalflux = totalflux + membrane_flux
    enddo
    	
    ! Next solve for ID concentrations of parent in medium, %Cmedium(:)
    Kd = chemo(ichemo)%medium_diff_coef
    A = well_area
    d = total_volume/A
    dX = d/N1D
    dV = A*dX
    !Cex = Caverage(MAX_CHEMO + ichemo)
    do k = 1,ndt
	    F(1) = totalflux
	    do i = 2,N1D
		    F(i) = Kd*A*(Cmedium(i-1) - Cmedium(i))/dX
	    enddo
	    F(N1D+1) = 0
	    do i = 1,N1D
		    Cmedium(i) = Cmedium(i)*(1 - decay_rate) + (F(i) - F(i+1))*(dt/ndt)/dV
	    enddo
	    Cex = Cmedium(1)
    enddo
    !write(*,*) 'istep,Cex: ',istep,Cex
enddo
C = 0
Clabel = 0
n = 0
do kcell = 1,nlist
   	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	n = n+1
    C = C + cp%Cin(ichemo)
    Clabel = Clabel + cp%Cin(ichemo+1)
enddo
Caverage(ichemo) = C/n
Caverage(ichemo+1) = Clabel/n
Caverage(MAX_CHEMO + ichemo) = Cex
if (Cex0 == 0) then	! stop processing when the parent drug is removed 
!	Caverage(ichemo) = 0
	chemo(ichemo)%present = .false.
endif
end subroutine

end module

