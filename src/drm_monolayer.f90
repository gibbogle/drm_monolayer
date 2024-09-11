!-----------------------------------------------------------------------------------------
! Units:
!     time				s = seconds
!     distance			cm
!     volume			cm^3
!     mass				micromole = 10^-6 mol = mumol
!     flux				mumol/s
!     concentration		mumol/cm^3 = mM
!-----------------------------------------------------------------------------------------
module drm_monolayer_mod
use global
use chemokine
use ode_diffuse
use cellstate
!use winsock  
use colony
use transfer
use metabolism
!use Tcp_mod

IMPLICIT NONE

real(REAL_KIND) :: phdist0(NP)

contains 

!-----------------------------------------------------------------------------------------
! This subroutine is called to initialize a simulation run. 
! ncpu = the number of processors to use
! infile = file with the input data
! outfile = file to hold the output 
!-----------------------------------------------------------------------------------------
subroutine Setup(ncpu,infile,outfile,ok)
integer :: ncpu
character*(*) :: infile, outfile
logical :: ok
character*(64) :: msg
integer :: ichemo, error, kcell, idrug, ityp
real(REAL_KIND) :: tgrowth(MAX_CELLTYPES)
type(cycle_parameters_type),pointer :: ccp
type(metabolism_type), pointer :: mp
logical :: isopen

ok = .true.
initialized = .false.
par_zig_init = .false.
colony_simulation = .false.

inputfile = infile
outputfile = outfile
call logger("ReadCellParams new")
call ReadCellParams(ok)
if (.not.ok) return
call logger("did ReadCellParams")
write(*,*) 'seed: ',seed 

!open(nfphase, file='phase.log',status='replace')
start_wtime = wtime()

if (ncpu == 0) then
	ncpu = ncpu_input
endif
Mnodes = ncpu
write(logmsg,*) 'ncpu: ',ncpu 
call logger(logmsg)

#if defined(OPENMP) || defined(_OPENMP)
    call logger("OPENMP defined")
    call omp_initialisation(ok)
    if (.not.ok) return
#else
    call logger("OPENMP NOT defined")
    if (Mnodes > 1) then
        write(logmsg,'(a)') 'No OpenMP, using one thread only'
        call logger(logmsg)
        Mnodes = 1
    endif
#endif

call ArrayInitialisation(ok)
if (.not.ok) return
call logger('did ArrayInitialisation')


use_cell_cycle = .true.
!if (use_metabolism) then
	chemo(OXYGEN)%controls_growth = .false.     
	chemo(OXYGEN)%controls_death = .false.
	chemo(GLUCOSE)%controls_growth = .false.
	chemo(GLUCOSE)%controls_death = .false.
!else
!	chemo(OXYGEN)%controls_growth = .true.     
!	chemo(OXYGEN)%controls_death = .true.
!	chemo(GLUCOSE)%controls_growth = .true.
!	chemo(GLUCOSE)%controls_death = .true.
!endif

call SetupChemo

! New cell cycle formulation - need a value for max (unconstrained) growth rate
!use_volume_method = .not.use_cell_cycle
!if (use_cell_cycle .and. .not.use_volume_based_transition) then
!    use_constant_growthrate = .true.
!endif
! Growth occurs during G1, S and G2, not in checkpoints
do ityp = 1,2
	ccp => cc_parameters(ityp)
	tgrowth(ityp) = ccp%T_G1 + ccp%T_S + ccp%T_G2
	max_growthrate(ityp) = Vdivide0/(2*tgrowth(ityp))
	write(nflog,*) 'ityp, Vdivide0, max_growthrate(ityp): ',ityp,Vdivide0, max_growthrate(ityp)
enddo

mp => phase_metabolic(1)
call SetupMetabolism(mp,ok)
if (.not.ok) stop

! Cell synchronisation
!use_synchronise = .false.  ! set in main
!synch_phase = G1_phase
!synch_fraction = 0.0    !0.7    ! 2.3/3.3

!write(nflog,*) 'before PlaceCells: npar_uni, npar_rnor = ',npar_uni,npar_rnor

call PlaceCells(ok)

!write(nflog,*) 'after PlaceCells: npar_uni, npar_rnor = ',npar_uni,npar_rnor

!call checkInitialPhaseDistribution
!write(*,*) 'stopping after PlaceCells'
!stop

if (.not.ok) return

call CreateMastercell

istep = 0
do ichemo = 1,NUTS
	if (chemo(ichemo)%used) then
		call InitConcs(ichemo)
		call SetupMedium(ichemo)
	endif
enddo
call UpdateChemomap
call AdjustMM
!call SetInitialGrowthRate
NATP_tag = 0
NGLN_tag = 0
Nradiation_tag = 0
Ndrug_tag = 0
Ndrug_tag = 0
Nradiation_dead = 0
Ndrug_dead = 0
NATP_dead = 0
NGLN_dead = 0
ndivided = 0
Ndying = 0
Ndead = 0

N_checkpoint = 0
ntphase = 0

ndoublings = 0
doubling_time_sum = 0
ncells_mphase = 0

!radiation_dosed = .false.
t_simulation = 0
total_dMdt = 0
chemo(:)%total_flux_prev = 0
t_lastmediumchange = 0
medium_change_step = .false.
limit_stop = .false.
kcell_dbug = 0
allocate(nphase(0:ndays*24,8))
!write(logmsg,'(a,i6)') 'Startup procedures have been executed: initial T cell count: ',Ncells0
!call logger(logmsg)
!call averages
!call GenerateSFlookup(1)
use_permute = .true.
rad_count = 0
! To compute average phase time
phase_exit_time_sum = 0
npet = 0
!call counter
!count_Nlethal = 0
!count_totDSB = 0
end subroutine

!----------------------------------------------------------------------------------------- 
!----------------------------------------------------------------------------------------- 
subroutine show_volume_data
integer :: kcell
real(REAL_KIND) :: Vsum, Vdivsum

write(nflog,*) 'Volume data:'
write(nflog,'(a,L)') 'use_divide_time_distribution: ',use_divide_time_distribution
write(nflog,'(a,L)') 'use_V_dependence: ',use_V_dependence
Vsum = 0
Vdivsum = 0
do kcell = 1,nlist
	write(nflog,'(i6,2f6.2)') kcell,cell_list(kcell)%V, cell_list(kcell)%divide_volume
	Vsum = Vsum + cell_list(kcell)%V
	Vdivsum = Vdivsum + cell_list(kcell)%divide_volume
enddo
write(nflog,*)
write(nflog,'(a,f6.2)') 'Average initial volume: ',Vsum/nlist
write(nflog,'(a,f6.2)') 'Average divide volume: ',Vdivsum/nlist
end subroutine

!----------------------------------------------------------------------------------------- 
! Initialise medium concentration
!-----------------------------------------------------------------------------------------
subroutine SetupMedium(ichemo)
integer :: ichemo, im

if (chemo(ichemo)%present) then
    Caverage(MAX_CHEMO+ichemo) = chemo(ichemo)%bdry_conc
    Cmediumave(ichemo) = chemo(ichemo)%bdry_conc
else
    Caverage(MAX_CHEMO+ichemo) = 0
    Cmediumave(ichemo) = 0
endif
if (ichemo <= NUTS) then
    C_OGL(ichemo,:) = chemo(ichemo)%bdry_conc
endif
if (ichemo <= NUTS) then
	chemo(ichemo)%Cmedium = chemo(ichemo)%bdry_conc
endif
if (ichemo > NUTS) then
	chemo(ichemo)%Cmedium = 0
endif
end subroutine

!----------------------------------------------------------------------------------------- 
!-----------------------------------------------------------------------------------------
subroutine omp_initialisation(ok)
logical :: ok
integer :: npr, nth

ok = .true.
#if defined(OPENMP) || defined(_OPENMP)
write(logmsg,'(a,i2)') 'Requested Mnodes: ',Mnodes
call logger(logmsg)
npr = omp_get_num_procs()
write(logmsg,'(a,i2)') 'Machine processors: ',npr
call logger(logmsg)

nth = omp_get_max_threads()
write(logmsg,'(a,i2)') 'Max threads available: ',nth
call logger(logmsg)
if (nth < Mnodes) then
    Mnodes = nth
    write(logmsg,'(a,i2)') 'Setting Mnodes = max thread count: ',nth
	call logger(logmsg)
endif

call omp_set_num_threads(Mnodes)
!$omp parallel
nth = omp_get_num_threads()
write(logmsg,*) 'Threads, max: ',nth,omp_get_max_threads()
call logger(logmsg)
!$omp end parallel
#endif

call logger('did omp_initialisation')
!call test_omp1

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine test_omp
integer, parameter :: n = 10
integer :: i

integer :: sum1, sum2
integer, allocatable :: y1(:)
integer :: y2(n)

allocate(y1(n))
y1 = 1
y2 = 1

sum1 = 0
sum2 = 0
!$omp parallel do
do i = 1,n
	sum1 = sum1 + y1(i)
	sum2 = sum2 + y2(i)
enddo
!$omp end parallel do
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine ArrayInitialisation(ok)
logical :: ok
integer :: x,y,z,k, ichemo
integer :: MAXX, z1, z2, nc0, inflow
integer :: cog_size
real(REAL_KIND) :: d, rr(3)

ok = .false.
call RngInitialisation

! These are deallocated here instead of in subroutine wrapup so that when a simulation run ends
! it will still be possible to view the cell distributions and chemokine concentration fields.
!if (allocated(occupancy)) deallocate(occupancy)
if (allocated(cell_list)) deallocate(cell_list)
!if (allocated(allstate)) deallocate(allstate)
!if (allocated(ODEdiff%ivar)) deallocate(ODEdiff%ivar)
if (allocated(gaplist)) deallocate(gaplist)
!if (allocated(Cslice)) deallocate(Cslice)
if (allocated(nphase)) deallocate(nphase)
if (allocated(Psurvive)) deallocate(Psurvive)
call logger('did deallocation')

!nsteps_per_min = 1.0/DELTA_T
ngaps = 0
nlist = 0

write(logmsg,*) 'Initial count, max_nlist: ',initial_count, max_nlist
call logger(logmsg)

! How big is cell_list?
!write(*,*) 'size of cell_list: ',sizeof(cell_list(1)),sizeof(cell_list(1))*max_nlist
allocate(cell_list(max_nlist))
allocate(gaplist(max_ngaps))

ok = .true.

end subroutine

!----------------------------------------------------------------------------------------- 
!-----------------------------------------------------------------------------------------
subroutine test_rnor
integer :: kpar=0, n = 100000, i
real(8) :: R, sum, sum2, ave, var

sum = 0
sum2 = 0
do i = 1,n
    R = par_rnor(kpar)
    sum = sum + R
    sum2 = sum2 + R*R
enddo
ave = sum/n
var = (sum2 - n*ave*ave)/(n-1)
write(*,'(a,3f10.4)') 'ave, var, std: ',ave,var,sqrt(var)
end subroutine

!----------------------------------------------------------------------------------------- 
!-----------------------------------------------------------------------------------------
subroutine RngInitialisation
integer, allocatable :: zig_seed(:)
integer :: i
integer :: npar, grainsize = 32

npar = Mnodes
write(logmsg,*) 'npar = ',npar,seed
call logger(logmsg)
allocate(zig_seed(0:npar-1))
do i = 0,npar-1
    zig_seed(i) = seed(1)*seed(2)*(i+1)
enddo
call par_zigset(npar,zig_seed,grainsize)
par_zig_init = .true.

!call test_rnor
!stop
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine ReadCellParams(ok)
logical :: ok
integer :: i, idrug, imetab, nmetab, im, itestcase, Nmm3, ichemo, itreatment, iuse_extra, iuse_relax, iuse_par_relax, iuse_FD
integer :: iuse_oxygen, iuse_glucose, iuse_lactate, iuse_glutamine, iuse_othernutrient, iuse_drug, iuse_metab, iV_depend
integer :: iV_random, iuse_gd_all, iuse_divide_dist, iuse_lognormal, ityp
!integer ::  idrug_decay, imetab_decay
integer :: ictype, idisplay, isconstant, ioxygengrowth, iglucosegrowth, ilactategrowth, ioxygendeath, iglucosedeath
integer :: iuse_drop, iconstant, isaveprofiledata, isaveslicedata, iusecellcycle, iusemetabolism, ifullymixed, isynchronise
logical :: use_metabolites
integer :: isaveFACSdata
real(REAL_KIND) :: bdry_conc, percent, d_n_limit
real(REAL_KIND) :: sigma(2), DXmm, anoxia_tag_hours, anoxia_death_hours, aglucosia_tag_hours, aglucosia_death_hours
character*(12) :: drug_name
character*(1) :: numstr
type(cycle_parameters_type),pointer :: ccp
logical :: write_hourly_results

ok = .true.
chemo(:)%used = .false.
Vsite_cm3 = 0.2000E-08  ! from spheroid - to retain same scaling 

open(nfcell,file=inputfile,status='old')
read(nfcell,'(a)') header
if (header(1:3) == 'GUI') then
	gui_run_version = header
	header = 'DD/MM/YYYY header_string'
else
	read(nfcell,*) gui_run_version				! program run version number
endif
read(nfcell,*) dll_run_version				! DLL run version number 
write(*,*) dll_run_version
!read(nfcell,*) NX							! size of grid
read(nfcell,*) initial_count				! initial number of tumour cells
if (greens) then
    initial_count = NgreenCells
endif
write(*,*) 'initial_count: ',initial_count
read(nfcell,*) iuse_lognormal
use_exponential_cycletime = (iuse_lognormal /= 1)
read(nfcell,*) divide_time_median(1)
read(nfcell,*) divide_time_shape(1)
read(nfcell,*) divide_time_median(2)
read(nfcell,*) divide_time_shape(2)
read(nfcell,*) iV_depend
read(nfcell,*) iV_random
read(nfcell,*) ndays							! number of days to simulate
!read(nfcell,*) d_n_limit					! possible limit on diameter or number of cells
!diam_count_limit = d_n_limit
read(nfcell,*) DELTA_T						! time step size (sec)
!read(nfcell,*) NXB							! size of coarse grid = NXB = NYB
!read(nfcell,*) NZB							! size of coarse grid = NZB
!read(nfcell,*) DXF							! fine grid spacing, off-lattice model (um)
!read(nfcell,*) a_separation
!read(nfcell,*) a_force
!read(nfcell,*) c_force
!read(nfcell,*) x0_force
!read(nfcell,*) x1_force
!read(nfcell,*) kdrag
!read(nfcell,*) frandom
read(nfcell,*) NT_CONC						! number of subdivisions of DELTA_T for diffusion computation
!read(nfcell,*) Nmm3							! number of cells/mm^3
!DXmm = 1.0/(Nmm3**(1./3))
!DELTA_X = DXmm/10							! mm -> cm
!Vsite_cm3 = DELTA_X*DELTA_X*DELTA_X			! total site volume (cm^3)
!read(nfcell,*) fluid_fraction				! fraction of the (non-necrotic) tumour that is fluid
read(nfcell,*) Vcell_pL                     ! nominal cell volume in pL
read(nfcell,*) well_area                    ! well bottom area (cm^2)
read(nfcell,*) medium_volume0				! initial total volume (cm^3)
read(nfcell,*) ifullymixed					! medium is fully mixed
fully_mixed = (ifullymixed == 1)
read(nfcell,*) Vdivide0						! nominal cell volume multiple for division
read(nfcell,*) dVdivide						! variation about nominal divide volume
read(nfcell,*) MM_THRESHOLD					! O2 concentration threshold Michaelis-Menten "soft-landing" (uM)
read(nfcell,*) anoxia_threshold			    ! O2 threshold for anoxia (uM)
read(nfcell,*) anoxia_tag_hours				! hypoxic time leading to tagging to die by anoxia (h)
read(nfcell,*) anoxia_death_hours			! time after tagging to death by anoxia (h)
read(nfcell,*) aglucosia_threshold			! O2 threshold for aglucosia (uM)
read(nfcell,*) aglucosia_tag_hours			! hypoxic time leading to tagging to die by aglucosia (h)
read(nfcell,*) aglucosia_death_hours		! time after tagging to death by aglucosia (h)
read(nfcell,*) itestcase                    ! test case to simulate
read(nfcell,*) seed(1)						! seed vector(1) for the RNGs
read(nfcell,*) seed(2)						! seed vector(2) for the RNGs
read(nfcell,*) ncpu_input					! for GUI just a placeholder for ncpu, used only when execute parameter ncpu = 0
read(nfcell,*) Ncelltypes					! maximum number of cell types in the spheroid
do ictype = 1,Ncelltypes
	read(nfcell,*) percent
	celltype_fraction(ictype) = percent/100
enddo
read(nfcell,*) NT_GUI_OUT					! interval between GUI outputs (timesteps)
read(nfcell,*) show_progeny                 ! if != 0, the number of the cell to show descendents of
read(nfcell,*) iuse_oxygen		! chemo(OXYGEN)%used
!read(nfcell,*) ioxygengrowth
!chemo(OXYGEN)%controls_growth = (ioxygengrowth == 1)
!read(nfcell,*) ioxygendeath
!chemo(OXYGEN)%controls_death = (ioxygendeath == 1)
read(nfcell,*) chemo(OXYGEN)%diff_coef
read(nfcell,*) chemo(OXYGEN)%medium_diff_coef
read(nfcell,*) chemo(OXYGEN)%membrane_diff_in
read(nfcell,*) chemo(OXYGEN)%membrane_diff_out
read(nfcell,*) chemo(OXYGEN)%bdry_conc
read(nfcell,*) iconstant
chemo(OXYGEN)%constant = (iconstant == 1)
read(nfcell,*) chemo(OXYGEN)%max_cell_rate
read(nfcell,*) chemo(OXYGEN)%MM_C0
read(nfcell,*) chemo(OXYGEN)%Hill_N
read(nfcell,*) iuse_glucose		!chemo(GLUCOSE)%used
!read(nfcell,*) iglucosegrowth
!chemo(GLUCOSE)%controls_growth = (iglucosegrowth == 1)
!read(nfcell,*) iglucosedeath
!chemo(GLUCOSE)%controls_death = (iglucosedeath == 1)
read(nfcell,*) C_G_base
read(nfcell,*) chemo(GLUCOSE)%diff_coef
read(nfcell,*) chemo(GLUCOSE)%medium_diff_coef
read(nfcell,*) chemo(GLUCOSE)%membrane_diff_in
read(nfcell,*) chemo(GLUCOSE)%membrane_diff_out
read(nfcell,*) chemo(GLUCOSE)%bdry_conc
if (chemo(GLUCOSE)%bdry_conc < 0) chemo(GLUCOSE)%bdry_conc = -chemo(GLUCOSE)%bdry_conc*C_G_base
chemo(GLUCOSE)%dose_conc = chemo(GLUCOSE)%bdry_conc
read(nfcell,*) iconstant
chemo(GLUCOSE)%constant = (iconstant == 1)
read(nfcell,*) chemo(GLUCOSE)%max_cell_rate
read(nfcell,*) chemo(GLUCOSE)%MM_C0
read(nfcell,*) chemo(GLUCOSE)%Hill_N



read(nfcell,*) LQ(1)%alpha_H
read(nfcell,*) LQ(1)%beta_H
read(nfcell,*) LQ(1)%OER_am
read(nfcell,*) LQ(1)%OER_bm
read(nfcell,*) LQ(1)%K_ms
read(nfcell,*) LQ(2)%alpha_H
read(nfcell,*) LQ(2)%beta_H
read(nfcell,*) LQ(2)%OER_am
read(nfcell,*) LQ(2)%OER_bm
read(nfcell,*) LQ(2)%K_ms
read(nfcell,*) iusecellcycle
use_cell_cycle = .true.
read(nfcell,*) isynchronise
synchronise = (isynchronise == 1)
call ReadCellCycleParameters(nfcell)
call ReadMcParameters(nfcell)
use_metabolism = .false.
read(nfcell,*) O2cutoff(1)
read(nfcell,*) O2cutoff(2)
read(nfcell,*) O2cutoff(3)
read(nfcell,*) hypoxia_threshold
write(nflog,*) 'hypoxia_threshold: ',hypoxia_threshold
read(nfcell,*) growthcutoff(1)
read(nfcell,*) growthcutoff(2)
read(nfcell,*) growthcutoff(3)
read(nfcell,*) Cthreshold
read(nfcell,*) Clabel_threshold
read(nfcell,*) spcrad_value
read(nfcell,*) Ndrugs_used
write(nflog,*) 'Ndrugs_used: ',Ndrugs_used
use_inhibiter = .false.
if (Ndrugs_used > 0) then
    call ReadDrugData(nfcell)
    chemo(DRUG_A)%used = .true.
    chemo(DRUG_A+1)%used = .true.
endif

is_radiation = .false.
if (use_events) then
	call ReadProtocol(nfcell)
	use_treatment = .false.
endif
close(nfcell)
call logger('Finished reading input data')

write(nflog,*) 'Hard-coded parameters:'
write(nflog,'(a,i3)') 'ATR_in_S: ',ATR_in_S
write(nflog,'(a,f6.1)') 'dose_threshold: ',dose_threshold

!call checkPD
!stop

! Try setting this for each cell unless use_cell_kcc2a_dependence
Kcc2a = get_Kcc2a(kmccp,CC_tot,CC_threshold_factor,cc_parameters(1)%T_G2/3600)
write(nflog,*) 'did get_Kcc2a: ',Kcc2a
single_cell = (initial_count==1)
write(nflog,*) 'single_cell: ',single_cell

! Rescale
chemo(OXYGEN)%membrane_diff_in = chemo(OXYGEN)%membrane_diff_in*Vsite_cm3/60		! /min -> /sec
chemo(OXYGEN)%membrane_diff_out = chemo(OXYGEN)%membrane_diff_out*Vsite_cm3/60		! /min -> /sec
chemo(OXYGEN)%max_cell_rate = chemo(OXYGEN)%max_cell_rate*1.0e6						! mol/cell/s -> mumol/cell/s
chemo(OXYGEN)%MM_C0 = chemo(OXYGEN)%MM_C0/1000										! uM -> mM
chemo(GLUCOSE)%membrane_diff_in = chemo(GLUCOSE)%membrane_diff_in*Vsite_cm3/60		! /min -> /sec
chemo(GLUCOSE)%membrane_diff_out = chemo(GLUCOSE)%membrane_diff_out*Vsite_cm3/60	! /min -> /sec
chemo(GLUCOSE)%max_cell_rate = chemo(GLUCOSE)%max_cell_rate*1.0e6					! mol/cell/s -> mumol/cell/s
chemo(GLUCOSE)%MM_C0 = chemo(GLUCOSE)%MM_C0/1000                                    ! uM -> mM

if (celltype_fraction(1) == 1.0) then
	write(nflog,*) 'Type 1 cells'
	selected_celltype = 1
elseif (celltype_fraction(2) == 1.0) then
	write(nflog,*) 'Type 2 cells'
	selected_celltype = 2
else
	write(logmsg,*) 'Error: cells must all be of the same type'
	call logger(logmsg)
	ok = .false.
	return
endif

if (chemo(OXYGEN)%Hill_N /= 1 .and. chemo(OXYGEN)%Hill_N /= 2) then
	call logger('Error: OXYGEN_HILL_N must be 1 or 2')
	ok = .false.
	return
endif
!if (chemo(GLUCOSE)%Hill_N /= 1 .and. chemo(GLUCOSE)%Hill_N /= 2) then
!	call logger('Error: GLUCOSE_HILL_N must be 1 or 2')
!	ok = .false.
!	return
!endif
!DXB = 4*DXF
if (S_phase_RR) then
!    S_phase_RR_progress = MM_THRESHOLD
else
    MM_THRESHOLD = MM_THRESHOLD/1000					! uM -> mM
endif
anoxia_threshold = anoxia_threshold/1000			! uM -> mM
aglucosia_threshold = aglucosia_threshold/1000		! uM -> mM
O2cutoff = O2cutoff/1000							! uM -> mM
hypoxia_threshold = hypoxia_threshold/1000			! uM -> mM
!relax = (iuse_relax == 1)
!use_parallel = (iuse_par_relax == 1)
!use_FD = (iuse_FD == 1)
chemo(OXYGEN)%used = (iuse_oxygen == 1)
chemo(GLUCOSE)%used = (iuse_glucose == 1)
if (.not.chemo(OXYGEN)%used) then
    chemo(OXYGEN)%controls_growth = .false.
    chemo(OXYGEN)%controls_death = .false.
endif
if (.not.chemo(GLUCOSE)%used) then
    chemo(GLUCOSE)%controls_growth = .false.
    chemo(GLUCOSE)%controls_death = .false.
endif

#if 0
if (use_exponential_cycletime) then
    divide_dist(1:2)%class = EXPONENTIAL_DIST
    do ityp = 1,2
        ccp => cc_parameters(ityp)
        divide_time_mean(ityp) = ccp%T_G1 + ccp%G1_mean_delay + ccp%T_S + ccp%S_mean_delay &
            + ccp%T_G2 + ccp%G2_mean_delay + ccp%T_M
    enddo
else
	write(nflog,*) 'divide_time_median: ',divide_time_median
	write(nflog,*) 'divide_time_shape: ',divide_time_shape
    divide_dist(1:2)%class = LOGNORMAL_DIST
    divide_time_median(1:2) = 60*60*divide_time_median(1:2)		! hours -> seconds
    sigma(1:2) = log(divide_time_shape(1:2))
    divide_dist(1:2)%p1 = log(divide_time_median(1:2))	
    divide_dist(1:2)%p2 = sigma(1:2)
    divide_time_mean(1:2) = exp(divide_dist(1:2)%p1 + 0.5*divide_dist(1:2)%p2**2)	! mean = median.exp(sigma^2/2)
    write(logmsg,'(a,24e12.4)') 'shape, sigma: ',divide_time_shape(1:2),sigma(1:2)
    call logger(logmsg)
    write(logmsg,'(a,4e12.4)') 'Median, mean divide time: ',divide_time_median(1:2)/3600,divide_time_mean(1:2)/3600
    call logger(logmsg)
!    call AdjustCycleTimes
endif
do ityp = 1,2
    ccp => cc_parameters(ityp)
    if (ccp%G1_mean_delay > 0) then
	    ccp%Pk_G1 = 1./ccp%G1_mean_delay    ! /sec
	else
	    ccp%Pk_G1 = 0
	endif
    if (ccp%S_mean_delay > 0) then
	    ccp%Pk_S = 1./ccp%S_mean_delay    ! /sec
	else
	    ccp%Pk_S = 0
	endif
    if (ccp%G2_mean_delay > 0) then
	    ccp%Pk_G2 = 1./ccp%G2_mean_delay    ! /sec
	else
	    ccp%Pk_G2 = 0
	endif
enddo
#endif

use_divide_time_distribution = (iuse_lognormal == 1)
use_V_dependence = (iV_depend == 1)
randomise_initial_volume = (iV_random == 1)
use_constant_divide_volume = (dVdivide == 0)
t_anoxia_limit = 60*60*anoxia_tag_hours				! hours -> seconds
anoxia_death_delay = 60*60*anoxia_death_hours		! hours -> seconds
t_aglucosia_limit = 60*60*aglucosia_tag_hours		! hours -> seconds
aglucosia_death_delay = 60*60*aglucosia_death_hours	! hours -> seconds
Vcell_cm3 = 1.0e-9*Vcell_pL							! nominal cell volume in cm3
Vdivide0 = Vdivide0*Vcell_cm3
dVdivide = dVdivide*Vcell_cm3

write(nflog,*) 'Vdivide0: ',Vdivide0, ' medium_volume0: ',medium_volume0

!write(logmsg,'(a,3e12.4)') 'DELTA_X, cell_radius: ',DELTA_X,cell_radius
!call logger(logmsg)
!write(logmsg,'(a,4e12.4)') 'Volumes: site, extra, cell (average, base): ',Vsite_cm3, Vextra_cm3, Vsite_cm3-Vextra_cm3, Vcell_cm3
!call logger(logmsg)

! Setup test_case
test_case = .false.
!if (itestcase /= 0) then
!    test_case(itestcase) = .true.
!endif
if (itestcase == 1) then
    use_constant_growthrate = .true.
    use_metabolism = .false.
    use_cell_cycle = .true.
endif

!if (mod(NX,2) /= 0) NX = NX+1					! ensure that NX is even
!NYB = NXB

if (use_PEST) then
    if (simulate_colony) then
        open(nfout,file=outputfile,status='replace')
    endif
else
    open(nfout,file=outputfile,status='replace')
!    write(nfout,'(a,a)') 'GUI version: ',gui_run_version
!    write(nfout,'(a,a)') 'DLL version: ',dll_run_version
!    write(nfout,*)
endif
write(nflog,*)
write(nflog,'(a,a)') 'GUI version: ',gui_run_version
write(nflog,'(a,a)') 'DLL version: ',dll_run_version
write(nflog,*)

if (.not.(use_PEST .and. simulate_colony)) then
    if (.not.use_PEST) then
	    open(nfres,file='drm_monolayer_ts.out',status='replace')
	    write(nflog,*) 'Opened drm_monolayer_ts.out'
    else
	    open(nfres,file=outputfile,status='replace')
	    write(nflog,*) 'Opened ',trim(outputfile)
    endif
    
if (.not.use_PEST .and. .false.) then
write(nfres,'(a)') 'date info GUI_version DLL_version &
istep hour Ncells(1) Ncells(2) Nviable Nnonviable &
NdrugA_dead(1) NdrugA_dead(2) &
Nradiation_dead(1) Nradiation_dead(2) Ntotal_dead &
Ntagged_drugA(1) Ntagged_drugA(2) Ntagged_radiation(1) Ntagged_radiation(2) &
f_viable f_hypoxic(1) f_hypoxic(2) f_hypoxic(3) f_clonohypoxic(1) f_clonohypoxic(2) f_clonohypoxic(3) f_growth(1) f_growth(2) f_growth(3) &
f_nogrow f_clonogenic plating_efficiency(1) plating_efficiency(2) &
EC_oxygen EC_glucose  EC_drugA EC_drugA_metab1 &
IC_oxygen IC_glucose  IC_drugA IC_drugA_metab1 &
medium_oxygen medium_glucose  medium_drugA medium_drugA_metab1 &
G1_phase G1_checkpoint S_phase G2_phase G2_checkpoint M_phase S_phase_nonarrest &
doubling_time oxygen_rate glycolysis_rate Ndivided'

write(logmsg,*) 'Wrote nfres header: '
call logger(logmsg)
endif
endif

Nsteps = ndays*24*60*60/DELTA_T		! DELTA_T in seconds
NT_DISPLAY = 2						! This is the updating interval (calls to get_summary) in the GUI version.  Not used by command-line version.
DT_DISPLAY = NT_DISPLAY*DELTA_T
write(logmsg,'(a,2i6,f6.0)') 'nsteps, NT_CONC, DELTA_T: ',nsteps,NT_CONC,DELTA_T
call logger(logmsg)
write(logmsg,'(a,i6,f6.0)') 'NT_DISPLAY, DT_DISPLAY: ',NT_DISPLAY, DT_DISPLAY
call logger(logmsg)

ok = .true.

end subroutine

!-----------------------------------------------------------------------------------------
! The cell cycle parameters include the parameters for radiation damage and repair, 
! and for the associated checkpoint duration limits Tcp(:).
! Time unit = hour
! Note: Now with the log-normal phase times there are no checkpoint delays except those
! created by cell damage (e.g. by radiation.)  Therefore now the mean cycle time is just the
! sum of G1, S, G2 and M times.
!-----------------------------------------------------------------------------------------
subroutine ReadCellCycleParameters(nf)
integer :: nf
type(cycle_parameters_type),pointer :: ccp
integer :: ityp, i
real(REAL_KIND) :: sigma, total

write(nflog,*) 'ReadCellCycleParameters:'
do ityp = 1,2
    write(nflog,*) 'ityp: ',ityp
    ccp => cc_parameters(ityp)

!read(nf,*) ccp%T_G1
!read(nf,*) ccp%T_S
!read(nf,*) ccp%T_G2
!read(nf,*) ccp%T_M
read(nf,*) ccp%f_G1
read(nf,*) ccp%f_S
read(nf,*) ccp%f_G2
read(nf,*) ccp%f_M
read(nf,*) ccp%Apoptosis_median
read(nf,*) ccp%Apoptosis_shape

divide_dist(ityp)%class = LOGNORMAL_DIST
divide_time_median(ityp) = 60*60*divide_time_median(ityp)		! hours -> seconds
sigma = log(divide_time_shape(ityp))
divide_dist(ityp)%p1 = log(divide_time_median(ityp))	
divide_dist(ityp)%p2 = sigma
divide_time_mean(ityp) = exp(divide_dist(ityp)%p1 + 0.5*divide_dist(ityp)%p2**2)	! mean = median.exp(sigma^2/2)
write(logmsg,'(a,i4,2e12.4)') 'ityp, shape, sigma: ',ityp,divide_time_shape(ityp),sigma
call logger(logmsg)
write(logmsg,'(a,2e12.4)') 'Median, mean divide time (hours): ',divide_time_median(ityp)/3600,divide_time_mean(ityp)/3600
call logger(logmsg)

call SteelMethod(ityp)

total = ccp%T_G1 + ccp%T_S + ccp%T_G2 + ccp%T_M
write(nflog,'(a,8f8.3)') 'T_G1,T_S,T_G2,T_M, total: ',ccp%T_G1,ccp%T_S,ccp%T_G2,ccp%T_M, total
ccp%T_G1 = 3600*ccp%T_G1                    ! hours -> seconds
ccp%T_S = 3600*ccp%T_S
ccp%T_G2 = 3600*ccp%T_G2
ccp%T_M = 3600*ccp%T_M
enddo

end subroutine

!-----------------------------------------------------------------------------------------
! Compute phase durations from phase fractions.  All corresponding to the average cycle time.
!-----------------------------------------------------------------------------------------
subroutine SteelMethod(ityp)
integer :: ityp
real(REAL_KIND) :: Tc, b
type(cycle_parameters_type),pointer :: ccp

ccp => cc_parameters(ityp)
Tc = divide_time_mean(ityp)/3600    ! seconds -> hours
b = log(2.0)/Tc
ccp%T_G1 = -(log(1-ccp%f_G1/2))/b
ccp%T_S = -(log(1-(ccp%f_G1+ccp%f_S)/2))/b - ccp%T_G1
ccp%T_M = log(1 + ccp%f_M)/b     ! Smith & Dendy
ccp%T_G2 = Tc - ccp%T_G1 - ccp%T_S - ccp%T_M
!ccp%T_G2 = -(log(1-(ccp%f_G1+ccp%f_S+ccp%f_G2)/2))/b - ccp%T_G1 - ccp%T_S
!ccp%T_M = Tc - ccp%T_G1 - ccp%T_S - ccp%T_G2

end subroutine

!-----------------------------------------------------------------------------------------
! NOT USED
!-----------------------------------------------------------------------------------------
subroutine AdjustCycleTimes
integer :: ityp
real(REAL_KIND) :: tmean, tsum, tfactor
type(cycle_parameters_type),pointer :: ccp

do ityp = 1,2
	ccp => cc_parameters(ityp)
	tmean = divide_time_mean(ityp)
	if (use_exponential_cycletime) then
	    tsum = ccp%T_G1 + ccp%T_S + ccp%T_G2 + ccp%T_M + ccp%G1_mean_delay + ccp%S_mean_delay + ccp%G2_mean_delay
	    tfactor = tmean/tsum
	    ccp%T_G1 = tfactor*ccp%T_G1
	    ccp%T_S = tfactor*ccp%T_S
	    ccp%T_G2 = tfactor*ccp%T_G2
	    ccp%T_M = tfactor*ccp%T_M
	    ccp%G1_mean_delay = tfactor*ccp%G1_mean_delay
	    ccp%S_mean_delay = tfactor*ccp%S_mean_delay
	    ccp%G2_mean_delay= tfactor*ccp%G2_mean_delay
    else    ! no checkpoint delays in log_timestep (unless there is radiation damage)
! Modified to keep T_M constant, scaling only T_G1, T_S, T_G2
!	    tsum = ccp%T_G1 + ccp%T_S + ccp%T_G2 + ccp%T_M
!	    tfactor = tmean/tsum
	    tsum = ccp%T_G1 + ccp%T_S + ccp%T_G2 ! + ccp%T_M
	    tfactor = (tmean - ccp%T_M)/tsum
	    ccp%T_G1 = tfactor*ccp%T_G1
	    ccp%T_S = tfactor*ccp%T_S
	    ccp%T_G2 = tfactor*ccp%T_G2
!	    ccp%T_M = tfactor*ccp%T_M
	    ccp%G1_mean_delay = 0
	    ccp%S_mean_delay = 0
	    ccp%G2_mean_delay= 0
	endif
!	ccp%Pk_G1 = 1./ccp%G1_mean_delay    ! /sec
!	ccp%Pk_S = 1./ccp%S_mean_delay    ! /sec
!	ccp%Pk_G2 = 1./ccp%G2_mean_delay    ! /sec
!	write(nflog,'(a,4e12.3)') 'Pk_G1, Pk_G2: ',ccp%Pk_G1,ccp%Pk_G2
enddo
end subroutine


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine ReadDrugData(nf)
integer :: nf
integer :: idrug, im, ictyp, ival
character*(16) :: drugname

write(logmsg,*) 'ReadDrugData'
call logger(logmsg)
if (allocated(drug)) then
	deallocate(drug)
endif
allocate(drug(Ndrugs_used))
do idrug = 1,Ndrugs_used
	read(nf,'(a)') drug(idrug)%classname        ! not used
	if (drug(idrug)%classname == 'TPZ') then
		drug(idrug)%drugclass = TPZ_CLASS
	elseif (drug(idrug)%classname == 'DNB') then
		drug(idrug)%drugclass = DNB_CLASS
	endif
	read(nf,*) drug(idrug)%nmetabolites		! currently all DNA-PK drugs have 1 metabolite
	drug(idrug)%use_metabolites = .true.	! currently simulate metabolites
	drug(idrug)%phase_dependent = .false.
	drug(idrug)%active_phase = .false.
    do im = 0,drug(idrug)%nmetabolites		! 0 = parent, 1 = metab_1, 2 = metab_2
		read(nf,'(a)') drugname
		if (im == 0) then
			drug(idrug)%name = drugname
		endif
		read(nf,*) drug(idrug)%diff_coef(im)
		read(nf,*) drug(idrug)%medium_diff_coef(im)
		read(nf,*) drug(idrug)%membrane_diff_in(im)
		read(nf,*) drug(idrug)%membrane_diff_out(im)
		read(nf,*) drug(idrug)%halflife(im)
		drug(idrug)%membrane_diff_in(im) = drug(idrug)%membrane_diff_in(im)*Vsite_cm3/60	! /min -> /sec
		drug(idrug)%membrane_diff_out(im) = drug(idrug)%membrane_diff_out(im)*Vsite_cm3/60	! /min -> /sec
		do ictyp = 1,ncelltypes
            read(nf,*) drug(idrug)%Kmet0(ictyp,im)
            read(nf,*) drug(idrug)%C2(ictyp,im)
            read(nf,*) drug(idrug)%KO2(ictyp,im)
            read(nf,*) drug(idrug)%Vmax(ictyp,im)
            read(nf,*) drug(idrug)%Km(ictyp,im)
            read(nf,*) drug(idrug)%n_O2(ictyp,im)
            drug(idrug)%Vmax(ictyp,im) = drug(idrug)%Vmax(ictyp,im)/60						! /min -> /sec
            drug(idrug)%Kmet0(ictyp,im) = drug(idrug)%Kmet0(ictyp,im)/60					! /min -> /sec
            drug(idrug)%KO2(ictyp,im) = 1.0e-3*drug(idrug)%KO2(ictyp,im)					! um -> mM
		enddo
	    if (drug(idrug)%name == 'EDU') then
			drug(idrug)%nmetabolites = 1
			drug(idrug)%phase_dependent = .true.
			drug(idrug)%active_phase(S_phase) = .true.
		endif
	    if (drug(idrug)%name == 'PI') then
			drug(idrug)%nmetabolites = 1
			drug(idrug)%phase_dependent = .true.
			drug(idrug)%active_phase(1:6) = .true.
		endif
    enddo
    write(nflog,*) 'drug: ',idrug,drug(idrug)%classname,'  ',drug(idrug)%name
    if (drug(1)%halflife(0) == 0) then
        Khalflife = 0
    else
        Khalflife = 0.693/drug(1)%halflife(0)
    endif
    use_drug_halflife = (Khalflife > 0)
    ! Repair inhibiting drug
    ! Now assume that any drug used is a DNA-PK inhibiter
    use_inhibiter = .true.
!    This is not needed now
!    if (drug(1)%SER_KO2(1,0) < 0) then ! drug 1 is a DNA repair inhibiter
!        DRUG_A_inhibiter = .true.
!        a_inhibit = drug(1)%SER_max(1,0)
!        b_inhibit = drug(1)%SER_Km(1,0)
!        use_inhibiter = drug(1)%sensitises(1,0)
!    endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Skip lines until the 'PROTOCOL' line
!-----------------------------------------------------------------------------------------
subroutine ReadProtocol(nf)
integer :: nf
integer :: itime, ntimes, kevent, ichemo, idrug, im, kevent_flush, drop_event
character*(64) :: line
character*(16) :: drugname
character*(1)  :: numstr
character*(1) :: fullstr
real(REAL_KIND) :: t, dt, vol, conc, O2conc, O2flush, dose, O2medium, glumedium
type(event_type) :: E, Eflush
logical :: flushing

write(logmsg,*) 'ReadProtocol'
call logger(logmsg)
total_volume = medium_volume0
chemo(DRUG_A:)%used = .false.
flushing = .false.
do
	read(nf,'(a)') line
	if (trim(line) == 'PROTOCOL') exit
enddo
read(nf,*) ntimes
if (ntimes == 0) then
	call logger('no events')
	Nevents = 0
	return
endif
Nevents = ntimes
if (allocated(event)) deallocate(event)
allocate(event(2*ntimes))
kevent = 0
do itime = 1,ntimes
	read(nf,'(a)') line
	write(nflog,'(a)') line
	if (trim(line) == 'DRUG') then
		kevent = kevent + 1
		event(kevent)%etype = DRUG_EVENT
		read(nf,'(a)') line
		write(nflog,'(a)') line
		drugname = trim(line)
		write(nflog,*) 'ndrugs_used: ',ndrugs_used
		do idrug = 1,ndrugs_used
			if (drugname == drug(idrug)%name) then
				ichemo = DRUG_A + 2*(idrug-1)
				exit
			endif
		enddo
		write(nflog,*) 'ichemo: ',ichemo
		! Need to copy drug(idrug) parameters to chemo(ichemo) 
		call CopyDrugParameters(idrug,ichemo)
		read(nf,*) t
		read(nf,*) dt
		read(nf,*) vol
		read(nf,*) O2conc
		read(nf,*) O2flush
		read(nf,*) conc
      !  if (CA_at_flushing) then
    		!CA_time = (t + dt)*60*60    ! assuming that CA occurs at drug flushing time
      !  else
      !      CA_time = 24*60*60          ! for Cho2 CDTD expts
      !  endif
		event(kevent)%time = t
		event(kevent)%ichemo = ichemo
		event(kevent)%idrug = idrug
		event(kevent)%volume = vol
		event(kevent)%conc = conc
        drug_conc0 = conc   ! used to compute effect of decay on Caverage, in grower().
		event(kevent)%O2conc = O2conc
		event(kevent)%dose = 0
		event(kevent)%full = .false.	
		chemo(ichemo)%used = .true.
		write(nflog,'(a,i3,2f8.3)') 'define DRUG_EVENT: volume, O2conc: ',kevent,event(kevent)%volume,event(kevent)%O2conc
		if (drug(idrug)%use_metabolites) then
			do im = 1,drug(idrug)%nmetabolites
				chemo(ichemo+im)%used = .true.
			enddo
		endif

!        if (.not.(idrug == 1 .and. use_inhibiter)) then     ! the flushing MEDIUM_EVENT is not added if the drug is an inhibiter
!		    kevent = kevent + 1
            flushing = .true.
            if (dt < 0) then    ! this signals a CDTD expt for which CA_time = flushing time, i.e. Cho1 only.  Otherwise CA_time takes the input parameter value.
                dt = -dt
                CA_time_h = dt
            endif
		    Eflush%etype = MEDIUM_EVENT
! assuming DRUG_EVENT is at time 0, and RADIATION_EVENT is at time > 0
! because RADIATION_EVENT t > 0 implies a delay of one time step (DELTA_T) before IR
! Must add an increment slightly less than DELTA_T to ensure washout is at the right time step.
            Eflush%time = 0.999*DELTA_T/3600 + dt        
		    Eflush%ichemo = 0
!		    event(kevent)%volume = medium_volume0
		    Eflush%volume = total_volume
		    Eflush%conc = 0
		    Eflush%O2medium = O2flush		
		    Eflush%glumedium = chemo(GLUCOSE)%dose_conc
		    Eflush%lacmedium = chemo(LACTATE)%dose_conc	
		    Eflush%full = .false.	
		    Eflush%dose = 0
		    write(nflog,'(a,i3,2f8.3)') 'define MEDIUM_EVENT: volume: ',kevent,Eflush%volume,Eflush%O2medium
!		endif
	elseif (trim(line) == 'MEDIUM') then
		kevent = kevent + 1
		event(kevent)%etype = MEDIUM_EVENT
		read(nf,*) t
		read(nf,*) vol
		read(nf,*) fullstr
		read(nf,*) O2medium
		read(nf,*) glumedium
		event(kevent)%time = t
		event(kevent)%volume = vol	
		event(kevent)%ichemo = 0
		event(kevent)%O2medium = O2medium
		event(kevent)%glumedium = glumedium
		event(kevent)%lacmedium = chemo(LACTATE)%dose_conc
		event(kevent)%full = (trim(fullstr) == 'Y' .or. trim(fullstr) == 'y')
		event(kevent)%dose = 0
		write(nflog,'(a,i3,2f8.3)') 'define MEDIUM_EVENT: volume: ',kevent,event(kevent)%volume,event(kevent)%O2medium
	elseif (trim(line) == 'RADIATION') then
!        is_radiation = .true.      ! moved to ProcessEvent
		kevent = kevent + 1
		event(kevent)%etype = RADIATION_EVENT
		read(nf,*) t
		if (use_synchronise) then
		    t = 0.001   ! might need to be > 0
        endif
        CA_time_h = CA_time_h + t
		read(nf,*) dose
		event(kevent)%time = t
		event(kevent)%dose = dose	
		event(kevent)%ichemo = 0
		event(kevent)%volume = 0
		event(kevent)%conc = 0
		if (use_PEST .and. nphase_hours > 0) then
		    phase_hour(1:nphase_hours) = t + phase_hour(1:nphase_hours)
		endif
		write(nflog,'(a,i3,a,f6.1,a,f6.2)') 'Radiation event: ', kevent,'  dose: ',dose,' hour: ',t
        IR_time_h = t
	endif
enddo
Nevents = kevent
! Check that the DRUG event has cdrug > 0
drop_event = 0
do kevent = 1,Nevents
    if (event(kevent)%etype == DRUG_EVENT) then
        if (event(kevent)%conc == 0) then
            drop_event = kevent
            exit
        endif
    endif
enddo
if (drop_event > 0) then    ! need to remove DRUG event from the list
    do kevent = drop_event+1, Nevents
        event(kevent-1) = event(kevent)
    enddo
    Nevents = Nevents - 1
!    CA_time = 18*60*60
    flushing = .false.
	do idrug = 1,ndrugs_used
		if (drugname == drug(idrug)%name) then
			ichemo = DRUG_A + 2*(idrug-1)
			exit
		endif
	enddo
	chemo(ichemo)%used = .false.
	if (drug(idrug)%use_metabolites) then
		do im = 1,drug(idrug)%nmetabolites
			chemo(ichemo+im)%used = .false.
		enddo
	endif
endif
! Set events not done
! convert time from hours to seconds
! convert volume from uL to cm^3  NO LONGER - now using cm^3 everywhere
write(logmsg,*) 'nevents: ',nevents
call logger(logmsg)
if (flushing) then  ! add the flushing event - for now assume only one
    kevent_flush = 0
    do kevent = 1,Nevents
        if (Eflush%time < event(kevent)%time) then
            kevent_flush = kevent
            exit
        endif
    enddo
    if (kevent_flush == 0) then     ! add event after all other events
        event(Nevents+1) = Eflush
    else                            ! insert event at kevent_flush
        do kevent = Nevents,kevent_flush,-1
            event(kevent+1) = event(kevent)
        enddo
        event(kevent_flush) = Eflush
    endif
    Nevents = Nevents+1
endif
            
do kevent = 1,Nevents
	event(kevent)%done = .false.
	event(kevent)%time = event(kevent)%time*60*60
!	event(kevent)%volume = event(kevent)%volume*1.0e-3
	E = event(kevent)
	write(nflog,'(a,i3,f8.3,2i3,3f8.4)') 'event: ',kevent,E%time/3600,E%etype,E%ichemo,E%volume,E%conc,E%dose
enddo
! Check that events are sequential
do kevent = 1,Nevents-1
	if (event(kevent)%time >= event(kevent+1)%time) then
		write(logmsg,'(a,2(i4,f6.2))') 'Error: non-sequential event: ',kevent,event(kevent)%time,kevent+1,event(kevent+1)%time
		call logger(logmsg)
		stop
	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine CopyDrugParameters(idrug,ichemo)
integer :: idrug,ichemo
integer :: im, im1, im2
character*(1) :: numstr

im1 = 0
chemo(ichemo)%name = drug(idrug)%name
if (drug(idrug)%use_metabolites) then
	do im = 1,drug(idrug)%nmetabolites
		chemo(ichemo+im)%used = .true.
		chemo(ichemo+im)%name = trim(chemo(ichemo)%name) // '_metab'
		write(numstr,'(i1)') im
		chemo(ichemo+im)%name = trim(chemo(ichemo+im)%name) // numstr
	enddo
	im2 = 1
else
	im2 = 0
endif
do im = im1, im2
	chemo(ichemo+im)%diff_coef = drug(idrug)%diff_coef(im)
	chemo(ichemo+im)%medium_diff_coef = drug(idrug)%medium_diff_coef(im)
	chemo(ichemo+im)%membrane_diff_in = drug(idrug)%membrane_diff_in(im)
	chemo(ichemo+im)%membrane_diff_out = drug(idrug)%membrane_diff_out(im)
	chemo(ichemo+im)%halflife = drug(idrug)%halflife(im)
!	chemo(ichemo+im)%medium_dlayer = d_layer
	chemo(ichemo+im)%decay = (chemo(ichemo+im)%halflife > 0)
	if (chemo(ichemo+im)%decay) then
		chemo(ichemo+im)%decay_rate = DecayRate(chemo(ichemo+im)%halflife)
	else
		chemo(ichemo+im)%decay_rate = 0
	endif
	write(nflog,'(a,i4,2e12.3)') 'CopyDrugParameters: im,halflife,decay_rate: ',im,chemo(ichemo+im)%halflife,chemo(ichemo+im)%decay_rate
enddo
end subroutine



!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine PlaceCells(ok)
logical :: ok
integer :: kcell, k, ichemo, ityp, site(3), phase, kpar = 0
real(REAL_KIND) :: rsite(3)
real(REAL_KIND) :: fract(0:8), total
real(REAL_KIND) :: kt2cc_min, kt2cc_max, kcc2a_sum
!integer :: phase_count(0:8), hour
integer :: counts(NP)
type(cell_type), pointer :: cp
type(cycle_parameters_type),pointer :: ccp
	
!call checkInitialPhaseDistribution
!stop
ccp => cc_parameters(1)
counts = 0
lastID = 0
kcell = 0
Ncells_type = 0
t_irradiation = -1
rsite = [0.,0.,0.]
kt2cc_min = 9999
kt2cc_max = 0
kcc2a_sum = 0
!do k = 1,20
!    write(nflog,*) 'k,par_uni: ',k,par_uni(kpar)
!enddo
do kcell = 1,initial_count
!    write(*,*) 'PlaceCells: kcell: ',kcell,initial_count
	call AddCell(kcell,rsite)
!	write(*,*) 'Added kcell: ',kcell
!    write(nflog,*) 'kcell,par_uni,npar-rnor: ',kcell,par_uni(kpar),npar_rnor
    cp => cell_list(kcell)
	kt2cc_min = min(kt2cc_min,cp%kt2cc)
	kt2cc_max = max(kt2cc_max,cp%kt2cc)
	phase = cp%phase
	phase = min(phase,M_phase)
	counts(phase) = counts(phase) + 1
	if (phase == G2_phase) then
	    cp%t_start_G2 = 0
	endif
	kcc2a_sum = kcc2a_sum + cp%kcc2a
!    write(*,*) 'Kcc2a: ',cp%Kcc2a
enddo
!write(nflog,'(a,2f8.3)') 'kcc2a_ave, calc: ',kcc2a_ave,kcc2a_sum/initial_count
!write(*,'(a,2f10.4)') 'kt2cc min, max: ',kt2cc_min,kt2cc_max
!write(nflog,'(a,2f10.4)') 'kt2cc min, max: ',kt2cc_min,kt2cc_max
!write(*,*)
write(*,'(a,5i6)') 'Initial phase counts: ',counts
write(*,'(a,5f8.3)') 'Initial phase %ages: ',100.0*real(counts)/sum(counts)
write(*,*)
write(nflog,'(a,5i6)') 'Initial phase counts: ',counts
write(nflog,'(a,5f7.2)') 'Initial phase %ages: ',100.0*real(counts)/sum(counts)
phdist0 = 100.0*real(counts)/sum(counts)
write(nflog,*) 'Ncells_Mphase: ',Ncells_Mphase
nlist = kcell-1
Ncells = nlist
Ncells0 = Ncells
Nviable = Ncells_type
!Nreuse = 0	

ok = .true.
!write(logmsg,*) 'idbug: ',idbug
!call logger(logmsg)
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine AddCell(kcell, rsite)
integer :: kcell
real(REAL_KIND) :: rsite(3)
integer :: ityp, k, kpar = 0
real(REAL_KIND) :: v(3), c(3), R1, R2, V0, Tdiv, Vdiv, p(3), R, gfactor, kfactor
real(8) :: kcc2a_std = 0.7
type(cell_type), pointer :: cp
type(cycle_parameters_type),pointer :: ccp
	
cp => cell_list(kcell)
cp%ID = kcell
cp%state = ALIVE
cp%generation = 1
cp%birthtime = 0
cp%rad_state = 0
cp%celltype = 1     !random_choice(celltype_fraction,Ncelltypes,kpar)
ityp = cp%celltype
ccp => cc_parameters(ityp)
Ncells_type(ityp) = Ncells_type(ityp) + 1
cp%Iphase = .true.
cp%totMis = 0
!cp%nspheres = 1
! Set cell's mitosis duration as a Gaussian rv
!R = par_rnor(kpar)	! N(0,1)
!cp%mitosis_duration = (1 + mitosis_std*R)*ccp%T_M
cp%mitosis_duration = get_mitosis_duration()
!write(nflog,'(a,i6,2f6.3)') 'mitosis_duration: ',kcell,ccp%T_M/3600,cp%mitosis_duration/3600
V0 = Vdivide0/2
kcell_now = kcell
!!!call set_divide_volume(cp, V0)  ! sets %divide_volume and %divide_time
!cp%dVdt = max_growthrate(ityp)
cp%metab = phase_metabolic(1)
!cp%metab%I_rate = r_Iu	! this is just to ensure that initial growth rate is not 0
cp%metab%C_GlnEx_prev = 0

! Jaiswal
!R = par_rnor(kpar)
!if (abs(R) > 2) R = R/2
!kfactor = max(0.0,1 + R*jaiswal_std)
R = par_uni(kpar)
kfactor = 1 + (R - 0.5)*jaiswal_std
if (single_cell .or. test_run .OR. use_no_random) kfactor = 1
cp%kt2cc = kt2cc*kfactor
R = par_uni(kpar)
kfactor = 1 + (R - 0.5)*jaiswal_std
if (single_cell .or. test_run .OR. use_no_random) kfactor = 1
cp%ke2cc = ke2cc*kfactor

!if (kcell <= 100) write(nflog,'(a,5f10.4)') 'Jaiswal R: ',R,kfactor,cp%kt2cc
cp%CC_act = 0   ! CC_act0
!R = par_uni(kpar)
!cp%CC_act = CC_act0 + R*(CC_threshold - CC_act0)
cp%ATR_act = 0
cp%ATM_act = 0
cp%G2_time = 0
!if (use_volume_method) then
!    !cp%divide_volume = Vdivide0
!    if (initial_count == 1) then
!	    cp%V = 0.9*cp%divide_volume
!    else
!    !	cp%V = (0.5 + 0.49*par_uni(kpar))*cp%divide_volume
!        if (randomise_initial_volume) then
!	        cp%V = cp%divide_volume*0.5*(1 + par_uni(kpar))
!        else
!	        cp%V = cp%divide_volume/1.6
!        endif
!    endif
!    cp%t_divide_last = 0    ! not correct
!else	! use cell cycle
!    cp%NL1 = 0
!    cp%NL2 = 0
!    cp%N_PL = 0
!    cp%N_Ch1 = 0
!    cp%N_Ch2 = 0
!    cp%irrepairable = .false.
    ! Need to assign phase, volume to complete phase, current volume
!    if (use_cell_kcc2a_dependence) then
!        cp%Kcc2a = get_Kcc2a(kmccp,CC_tot,CC_threshold_factor,cp%fg(3)*ccp%T_G2/3600)
!        R = par_rnor(kpar)
!        if (abs(R) > 2) R = R/2
!        kfactor = max(0.5,1 + R*kcc2a_std)
!        kfactor = 0.9*kfactor
!        cp%Kcc2a = kfactor*Kcc2a_ave
!        if (kcell <= 100) write(nflog,'(a,5f10.4)') 'Kcc2a R: ',R,kfactor,cp%kcc2a
!    endif

!!!    call SetInitialCellCycleStatus(kcell,cp)
!   added
    cp%phase = G1_phase
    cp%progress = 0
!!!--------------------------------------------------------    
    
    if (cp%phase == M_phase) then
        write(*,'(a,i8,f6.3)') 'M_phase, mitosis_duration: ',kcell,cp%mitosis_duration/3600
        stop
    endif
!endif
!cp%dVdt = max_growthrate(ityp)
!if (use_metabolism) then	! Fraction of I needed to divide = fraction of volume needed to divide
!	cp%metab%I2Divide = get_I2Divide(cp)
!	cp%metab%Itotal = cp%metab%I2Divide*(cp%V - V0)/(cp%divide_volume - V0)
!endif
!cp%radius(1) = (3*cp%V/(4*PI))**(1./3.)
!cp%centre(:,1) = rsite 
!cp%site = rsite/DELTA_X + 1
!cp%d = 0
!cp%birthtime = 0
!cp%growthrate = test_growthrate
!cp2%divide_volume = get_divide_volume()
!cp%d_divide = (3*cp%divide_volume/PI)**(1./3.)
!cp%mitosis = 0

!cp%ATP_tag = .false.
!cp%GLN_tag = .false.
cp%drug_tag = .false.
cp%radiation_tag = .false.
!cp%anoxia_tag = .false.
!cp%aglucosia_tag = .false.
!cp%growth_delay = .false.
cp%G2_M = .false.
cp%p_rad_death = 0

!cp%t_anoxia = 0
!cp%t_aglucosia = 0

!call get_random_vector3(v)	! set initial axis direction
!cp%d = 0.1*small_d
!c = cp%centre(:,1)
!cp%centre(:,1) = c + (cp%d/2)*v
!cp%centre(:,2) = c - (cp%d/2)*v
!cp%nbrs = 0
!cp%Cex = Caverage(1,1,1,:)	! initially the concentrations are the same everywhere
cp%Cin(OXYGEN) = chemo(OXYGEN)%bdry_conc
cp%Cin(GLUCOSE) = chemo(GLUCOSE)%bdry_conc
#if 0
cp%Cin(LACTATE) = chemo(LACTATE)%bdry_conc
cp%Cin(GLUTAMINE) = chemo(GLUTAMINE)%bdry_conc
cp%Cin(OTHERNUTRIENT) = chemo(OTHERNUTRIENT)%bdry_conc
cp%CFSE = generate_CFSE(1.d0)
#endif

!cp%growth_rate_factor = get_growth_rate_factor()
!cp%ATP_rate_factor = get_ATP_rate_factor()
!cp%ndt = ndt
cp%Psurvive = -1    ! flags Psurvive not yet computed
!if (cp%phase == G1_phase) write(*,'(a,i6,f8.3)') 'G1_phase cell: ',kcell,cp%progress

end subroutine

!--------------------------------------------------------------------------------------
! Steel: 
! Probability density function of progress through cell cycle: f(t) = 2b exp(-bt) 
! => cumulative distribution function F(t) = 2(1 - exp(-bt)) = fraction of cells less than t since creation
! To generate a variate from this CDF, first generate R from U(0,1)
! R = 2(1 - exp(-bt)), exp(-bt) = 1 - R/2, -bt = ln(1 - R/2)
! t = -(1/b)ln(1 - R/2)
! To do:
! set flags
!
! Note: assumes use of log-normal cycle time!
!
! Calculates synch_time for each cell.  This is the initial time point in the cell cycle 
! that all cells share, i.e. at time = 0 all cells are in the same phase and at the same 
! fractional point in the phase.
! For the radiation dose to occur while the cells are synchronised requires the dose time
! to be very close to the start, i.e. close to 0.
! NEED TO set initial cell volume!!!
!--------------------------------------------------------------------------------------
subroutine SetInitialCellCycleStatus(kcell,cp)
integer :: kcell
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
integer :: ityp, kpar = 0
real(REAL_KIND) :: Tc, Tmean, scale, b, t, R, tswitch(3), rVmax, V0, Vprev, fg(4), metab, f_CP, fp(4)
real(REAL_KIND) :: T_G1, T_S, T_G2, T_M, tleft, Vleft, Vphase, dth

if (use_exponential_cycletime) then
    write(*,*) 'SetInitialCellCycleStatus: Must use log-normal cycle time!'
    stop
endif
ityp = cp%celltype
ccp => cc_parameters(ityp)
!if (S_phase_RR) then
!if (use_synchronise .and. (initial_count == 1)) then
Tmean = divide_time_mean(ityp)
if (single_cell .or. test_run) then
    Tc = divide_time_mean(1)
else
    Tc = cp%divide_time         ! log-normal, implies %fg
endif
scale = Tc/Tmean
fg = cp%fg
metab = 1.0
f_CP = 1.0
fp(:) = metab*f_CP/fg(:)
T_G1 = ccp%T_G1/fp(1)
T_S = ccp%T_S/fp(2)
T_G2 = ccp%T_G2/fp(3)
if (test_run) then
    T_M = ccp%T_M/fp(4)
else
    T_M = cp%mitosis_duration
endif
if (use_cell_kcc2a_dependence) then
    cp%Kcc2a = get_Kcc2a(kmccp,CC_tot,CC_threshold_factor,T_G2/3600)
    cp%Kcc2a = min(cp%kcc2a, 0.9*CC_threshold)
!    write(nflog,*) 'T_G2, cp%Kcc2a: ',T_G2/3600,cp%Kcc2a
endif

if (use_synchronise) then
!    if (kcell == 1) then
!        write(nfphase,'(a,i4,f6.3)') 'synch_phase, synch_fraction: ',synch_phase, synch_fraction
!        write(nflog,'(a,i4,f6.3)') 'synch_phase, synch_fraction: ',synch_phase, synch_fraction
!        write(*,'(a,i4,f6.3)') 'synch_phase, synch_fraction: ',synch_phase, synch_fraction
!    endif
    if (synch_phase == G1_phase) then
        t = synch_fraction*T_G1
    elseif (synch_phase == S_phase) then
        t = T_G1 + synch_fraction*T_S
    elseif (synch_phase == G2_phase) then
        t = T_G1 + T_S + synch_fraction*T_G2
    else
        write(*,*) 'Error: SetInitialCellCycleStatus: bad synch_phase: ',synch_phase
        write(nflog,*) 'Error: SetInitialCellCycleStatus: bad synch_phase: ',synch_phase
        stop
    endif
    cp%progress = synch_fraction
    write(nflog,'(a,i4,f6.2,5f8.0)') 'synchronise: kcell,synch_fraction,T_G1,T_S,T_G2,T_M,t: ',kcell,synch_fraction,T_G1,T_S,T_G2,T_M,t
    if (write_nfres) write(nfres,'(a,i4,2f10.2)') 'use_synchronise: synch_phase, synch_fraction, t: ',synch_phase, synch_fraction, t
else
    b = log(2.0)/Tc
    R = par_uni(kpar)
    t = -(1/b)*log(1 - R/2)     ! cycle progression, log-normal r.v. (t/Tc = fractional progression)
endif
!if (kcell <= 10) write(*,*) 'kcell, t/Tc: ',kcell,t/Tc
tswitch(1) = T_G1 
tswitch(2) = tswitch(1) + T_S
tswitch(3) = tswitch(2) + T_G2

V0 = Vdivide0/2
rVmax = max_growthrate(ityp)

if (t < tswitch(1)) then
    cp%phase = G1_phase
    cp%fp = fp(1)
    cp%dVdt = cp%fp*rVmax
    cp%V = V0 + t*cp%dVdt
    cp%progress = t/T_G1
elseif (t <= tswitch(2)) then
    Vprev = V0 + fp(1)*rVmax*T_G1
    cp%phase = S_phase
    cp%fp = fp(2)
    cp%dVdt = cp%fp*rVmax
    cp%V = Vprev + (t - tswitch(1))*cp%dVdt 
    cp%progress = (t - tswitch(1))/T_S
elseif (t <= tswitch(3)) then
    Vprev = V0 + fp(1)*rVmax*T_G1 + fp(2)*rVmax*T_S
    cp%phase = G2_phase
    cp%fp = fp(3)
    cp%dVdt = cp%fp*rVmax
    cp%V = Vprev + (t - tswitch(2))*cp%dVdt
    cp%progress = (t - tswitch(2))/T_G2
    tleft = tswitch(3) - t
    Vleft = cp%dVdt*tleft
    Vphase = rVmax*T_G2*fp(3)
    !if (cp%progress < 0.01) write(*,'(a,i6,f6.3)') 'SetInitialCellCycleStatus: kcell, progress: ',kcell,cp%progress
!    if (kcell <= 10) then
!        write(*,*)
!        write(*,*) 'SetInitialCellCycleStatus'
!        write(*,'(a,2e12.3)') 'Total cycle volume: ',rVmax*(T_G1*fp(1) + T_S*fp(2) + T_G2*fp(3)),Vdivide0/2 
!        write(*,'(a,3f8.2)') 'phase times: ',T_G1/3600,T_S/3600,T_G2/3600
!        write(*,'(a,e12.3,3f8.3)') 'rVmax,fg: ',rVmax,fg(G1_phase),fg(S_phase),fg(G2_phase)
!        write(*,'(a,i6,3f6.2,e12.3)') 'kcell, fp(3), Tdiv, tleft, dVdt: ',kcell,fp(3),cp%divide_time/3600,tleft/3600,cp%dVdt
!        write(*,'(a,4e12.3,f6.3)') 'Vleft,Vphase,Vend,Vdivide0,progress: ',Vleft, Vphase, cp%V + Vleft, Vdivide0, cp%progress
!        if (kcell == 4) stop
!    endif
    if (use_Jaiswal) then   ! need to initialise CC_act
        cp%DSB = 0
        dth = (t - tswitch(2))/3600
        if (single_cell) write(nflog,'(a,3f8.3)') 'SetInitialCellCycleStatus: dth: ',t/3600, tswitch(2)/3600, dth
        call Jaiswal_update(cp,dth)
        cp%CC_act = min(cp%CC_act,0.95*CC_threshold)     ! to prevent premature mitosis
!        if (kcell_now > 0) write(nflog,'(a,i6,4f8.3)') 'SetInitialCellCycleStatus: G1 progress,initial CC_act: ',kcell,cp%progress,cp%CC_act,cp%ATM_act,cp%ATR_act
!        if (cp%CC_act > CC_threshold) write(nflog,'(a,i6,f8.1,3f8.3)') 'SetInitialCellCycleStatus: dth,T_G2,progress,CC_act: ',kcell,dth,T_G2/3600,cp%progress,cp%CC_act
        !write(nflog,*) 'stopping:'
        !write(*,*) 'stopping:'
        !stop
    endif
else    ! cell in mitosis
!    cp%fp = metab*f_CP/fg(M_phase)      ! not used
    cp%dVdt = 0
    cp%V = 2*V0
    R = par_uni(kpar)
!    cp%t_start_mitosis = -R*ccp%T_M
!    cp%t_start_mitosis = -R*cp%mitosis_duration
! Try this
    cp%t_start_mitosis = -(t - tswitch(3))
    cp%progress = (t - tswitch(3))/T_M
	ncells_mphase = ncells_mphase + 1
    cp%phase = dividing
!    write(nflog,'(a,i6,2f8.1)') 'Cell dividing: t_start_mitosis: ',kcell,cp%t_start_mitosis,cp%mitosis_duration
!    cp%progress = 1.0
!    if (kcell <= 10) write(*,'(a,2i6,2f6.3)') 'SetInitialCellCycleStatus: phase, mitosis_duration, t_start: ',kcell,cp%phase,cp%mitosis_duration/3600,cp%t_start_mitosis/3600
endif
cp%f_S_at_IR = cp%progress
if (single_cell) then
    write(*,*)
    write(*,*) 'Initial phase, progress: ',cp%phase,cp%progress
    write(nflog,*) 'Initial phase, progress: ',cp%phase,cp%progress
    write(*,*)
endif
cp%t_divide_last = -t
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine CreateMastercell
type(cell_type), pointer :: cp
type(cycle_parameters_type),pointer :: ccp

cp => master_cell
cp%ID = 0
cp%state = ALIVE
cp%generation = 1
cp%celltype = 1
ccp => cc_parameters(1)
cp%Iphase = .true.
cp%metab = phase_metabolic(1)
!cp%metab%I_rate = r_Iu	! this is just to ensure that initial growth rate is not 0
!cp%metab%C_GlnEx_prev = 0
!cp%ATP_tag = .false.
!cp%GLN_tag = .false.
cp%drug_tag = .false.
cp%radiation_tag = .false.
!cp%growth_delay = .false.
cp%G2_M = .false.
cp%p_rad_death = 0
cp%Cin(OXYGEN) = chemo(OXYGEN)%bdry_conc
cp%Cin(GLUCOSE) = chemo(GLUCOSE)%bdry_conc
#if 0
cp%Cin(LACTATE) = chemo(LACTATE)%bdry_conc
cp%Cin(GLUTAMINE) = chemo(GLUTAMINE)%bdry_conc
cp%Cin(OTHERNUTRIENT) = chemo(OTHERNUTRIENT)%bdry_conc
#endif
end subroutine


!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine checkInitialPhaseDistribution
real(8) :: Tc, Tmean, scale, b, t, R, tswitch(3), progress
integer :: k, phase, ityp = 1, kpar = 0
integer :: pcount(4), N = 10000
type(cycle_parameters_type), pointer :: ccp
type(cell_type), pointer :: cp

write(*,'(a)') 'checkInitialPhaseDistribution'
N = nlist
pcount = 0
do k = 1,N
    cp => cell_list(k)
    phase = cp%phase
    pcount(phase) = pcount(phase) + 1
enddo
write(*,'(4f8.4)') pcount/real(N)
return

ccp => cc_parameters(ityp)
tswitch(1) = ccp%T_G1 
tswitch(2) = tswitch(1) + ccp%T_S
tswitch(3) = tswitch(2) + ccp%T_G2

do k = 1,N
    Tc = DivideTime(ityp)     ! log-normally distributed
    Tmean = divide_time_mean(ityp)
    scale = Tc/Tmean
    tswitch(1) = scale*ccp%T_G1 
    tswitch(2) = tswitch(1) + scale*ccp%T_S
    tswitch(3) = tswitch(2) + scale*ccp%T_G2
    b = log(2.0)/Tc
    R = par_uni(kpar)
    t = -(1/b)*log(1 - R/2)     ! cycle progression, log-normal r.v. (t/Tc = fractional progression)
    if (t < tswitch(1)) then
        phase = G1_phase
        progress = t/ccp%T_G1
    elseif (t < tswitch(2)) then
        phase = S_phase
        progress = (t - tswitch(1))/ccp%T_S
    elseif (t < tswitch(3)) then
        phase = G2_phase
        progress = (t - tswitch(2))/ccp%T_G2
    else
        phase = M_phase
        progress = 1.0
    endif
    pcount(phase) = pcount(phase) + 1
enddo
write(*,'(4f8.4)') pcount/real(N)
end subroutine

!--------------------------------------------------------------------------------
! Not USED
!--------------------------------------------------------------------------------
subroutine setTestCell(kcell)
integer :: kcell
integer :: ityp
real(REAL_KIND) :: V0, Tdiv, Tgrowth, Tgrowth0, Tfixed, rVmax, phase_time1
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp

ityp = 1
ccp => cc_parameters(ityp)
cp => cell_list(kcell)
V0 = Vdivide0/2
Tdiv = divide_time_mean(ityp)
cp%divide_time = Tdiv
rVmax = max_growthrate(ityp)
write(nflog,'(a,3f8.0)') 'ccp%T_G1, ccp%T_S, ccp%T_G2: ',ccp%T_G1, ccp%T_S, ccp%T_G2
Tgrowth0 = ccp%T_G1 + ccp%T_S + ccp%T_G2
Tfixed = ccp%T_M + ccp%G1_mean_delay + ccp%G2_mean_delay
write(nflog,'(a,3f8.0)') 'ccp%T_M, ccp%G1_mean_delay, ccp%G2_mean_delay: ',ccp%T_M, ccp%G1_mean_delay, ccp%G2_mean_delay
Tgrowth = Tdiv - Tfixed
write(nflog,'(a,4f8.0)') 'Tdiv, Tfixed, Tgrowth0, Tgrowth: ',Tdiv, Tfixed, Tgrowth0, Tgrowth
!cp%fg = Tgrowth/Tgrowth0
!cp%divide_volume = V0 + Tgrowth*rVmax
!phase_time1 = cp%fg*ccp%T_G1
cp%phase = G1_phase
!cp%G1_time = phase_time1
cp%V = V0
cp%t_divide_last = 0
write(nflog,*) 'setTestCell'
write(nflog,*) 'phase: ',cp%phase
write(nflog,*) 'divide_time:', cp%divide_time
write(nflog,*) 'V: ',cp%V
write(nflog,*) 'divide_volume: ',cp%divide_volume
!write(nflog,*) 'fg: ',cp%fg
!write(nflog,*) 'G1_time: ',cp%G1_time
end subroutine


!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine InitConcs(ichemo)
integer :: ichemo
!integer :: i, kcell, site(3), ntvars
real(REAL_KIND) :: c0

if (istep == 0) then
	if (ichemo >= OXYGEN .and. ichemo < DRUG_A) then
		c0 = chemo(ichemo)%bdry_conc
	else	! drug or metabolite
		c0 = 0
	endif
	Caverage(ichemo) = c0
endif
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine ProcessEvent(radiation_dose)
integer :: kevent, ichemo, idrug, im, nmetab
real(REAL_KIND) :: V, C(MAX_CHEMO), radiation_dose
type(event_type) :: E
logical :: full
logical :: is_event
integer :: kcell
type(cell_type), pointer :: cp

!write(logmsg,*) 'ProcessEvent:'
!call logger(logmsg)
radiation_dose = -1
is_event = .false.
do kevent = 1,Nevents
	E = event(kevent)
	if (t_simulation >= E%time .and. .not.E%done) then
		write(nflog,'(a,i3,2f8.0,i3,2f10.4)') 'Event: ',E%etype,t_simulation,E%time,E%ichemo,E%volume,E%conc
		if (E%etype == RADIATION_EVENT) then
			radiation_dose = E%dose
			write(logmsg,'(a,i3,f8.0,f8.3)') 'RADIATION_EVENT: kevent, time, dose: ',kevent,t_simulation,E%dose
			call logger(logmsg)
!            is_radiation = (radiation_dose > 0) !.true. to cover the case of dose = 0 (control)
            is_radiation = .true.
            is_event = .true.
		elseif (E%etype == MEDIUM_EVENT) then
			write(logmsg,'(a,i6,2f8.3,f8.3,2f8.4)') 'MEDIUM_EVENT: E%time, time, volume, O2medium: ',istep,E%time/3600,t_simulation/3600,E%volume,E%O2medium
			call logger(logmsg)
            write(nflog,'(a,f8.3)') 'drug exposure time: ',(t_simulation - t_irradiation)/3600
            write(nflog,*) 'npar_uni, npar_rnor = ',npar_uni,npar_rnor
!            write(nfres,'(a,i6,f8.3)') 'washout: ',istep,t_simulation/3600
!            call washoutSF
			C = 0
			C(OXYGEN) = E%O2medium
			C(GLUCOSE) = E%glumedium
			C(DRUG_A:DRUG_A+1) = 0
            drug_conc0 = 0      ! for decay calc in grower() 
            do kcell = 1,10
                cp => cell_list(kcell)
                write(nflog,'(a,2i4,4f8.3)') 'kcell,phase,progress,DSB: ',kcell,cp%phase,cp%progress,cp%DSB(1:3,1)
            enddo
			V = E%volume
			full = E%full
			call MediumChange(V,C,full)
            is_event = .true.
		elseif (E%etype == DRUG_EVENT) then
			C = 0
			C(OXYGEN) = E%O2conc
			C(GLUCOSE) = chemo(GLUCOSE)%dose_conc
			ichemo = E%ichemo
			idrug = E%idrug
			C(ichemo) = E%conc
			V = E%volume
			write(logmsg,'(a,i4,2f8.3)') 'DRUG_EVENT: ichemo, volume, conc: ',ichemo,E%volume,E%conc
			call logger(logmsg)
			! set %present
			chemo(ichemo)%present = .true.
			chemo(ichemo)%bdry_conc = 0
			nmetab = drug(idrug)%nmetabolites
			do im = 1,nmetab
				if (chemo(ichemo + im)%used) then
					chemo(ichemo + im)%present = .true.
					chemo(ichemo + im)%bdry_conc = 0
				endif
			enddo
			full = E%full
			call MediumChange(V,C,full)
			call UpdateChemomap
!			call UpdateCbnd(0.0d0)
            is_event = .true.
            drug_time = t_simulation
		endif
		event(kevent)%done = .true.
	endif
enddo
!if (is_event) then
!    write(nflog,*) "processEvent: stopping"
!    stop
!endif
end subroutine

!----------------------------------------------------------------------------------
! Radiation treatment is stored in protocol(0)
! NOT USED NOW
!----------------------------------------------------------------------------------
subroutine Treatment(radiation_dose)
real(REAL_KIND) :: radiation_dose
integer :: i, idrug, ichemo, nmetab, im	!, ichemo_metab

radiation_dose = 0
do i = 1,protocol(0)%n
	if (t_simulation >= protocol(0)%tstart(i) .and. .not.protocol(0)%started(i)) then
		radiation_dose = protocol(0)%dose(i)
		protocol(0)%started(i) = .true.
		protocol(0)%ended(i) = .true.
		write(nflog,*) 'Radiation started: dose: ',radiation_dose
		exit
	endif
enddo
do idrug = 1,2
	ichemo = protocol(idrug)%ichemo
	if (idrug == 1) then
		nmetab = 2
	elseif (idrug == 2) then
		nmetab = 2
	endif
	do i = 1,protocol(idrug)%n
		if (i == 1 .and. t_simulation < protocol(idrug)%tstart(i)) then
			chemo(ichemo)%bdry_conc = 0
			do im = 1,nmetab
				chemo(ichemo + im)%bdry_conc = 0
			enddo
			exit
		endif
		if (t_simulation >= protocol(idrug)%tstart(i) .and. .not.protocol(idrug)%started(i)) then
			chemo(ichemo)%bdry_conc = protocol(idrug)%conc(i)
			protocol(idrug)%started(i) = .true.
			protocol(idrug)%ended(i) = .false.
			chemo(ichemo)%present = .true.
			call InitConcs(ichemo)
			call SetupMedium(ichemo)
			do im = 1,nmetab
				if (chemo(ichemo + im)%used) then
					chemo(ichemo + im)%present = .true.
					call InitConcs(ichemo + im)
					call SetupMedium(ichemo + im)
				endif
			enddo
			write(nflog,*) 'Started DRUG: ',chemo(ichemo)%name,chemo(ichemo)%bdry_conc, i
			exit
		endif
	enddo
	do i = 1,protocol(idrug)%n
		if (t_simulation >= protocol(idrug)%tend(i) .and. .not.protocol(idrug)%ended(i)) then
			chemo(ichemo)%bdry_conc = 0
			protocol(idrug)%ended(i) = .true.
			call InitConcs(ichemo)
			call SetupMedium(ichemo)
			do im = 1,nmetab
				chemo(ichemo + im)%bdry_conc = 0
				if (chemo(ichemo + im)%used) then
					call InitConcs(ichemo + im)
					call SetupMedium(ichemo + im)
				endif
			enddo
			write(nflog,*) 'Ended DRUG: ',chemo(ichemo)%name,i
			exit
		endif
	enddo
enddo	
end subroutine

!-----------------------------------------------------------------------------------------
! If the volume removed is Vr, the fraction of constituent mass that is retained
! in the medium is (Vm - Vr)/Vm.  The calculation does not apply to oxygen.
! Usually Vr = Ve.
! Revised treatment of concentrations, to allow for setting O2 in medium change
!
! Now only medium concentrations are stored, in Caverage(MAX_CHEMO+1:2*MAX_CHEMO)
! (oxygen is a special case, Caverage is actually conc at well bottom)
!-----------------------------------------------------------------------------------------
subroutine MediumChange(Ve,Ce,full)
real(REAL_KIND) :: Ve, Ce(:)
logical :: full
real(REAL_KIND) :: R, Vm, Vr, Vcells, mass(MAX_CHEMO), C
integer :: ichemo, idrug, im, iparent

write(nflog,*) 'MediumChange:'
write(nflog,'(a,f8.4)') 'Ve: ',Ve
write(nflog,'(a,13f8.4)') 'Ce: ',Ce
write(nflog,'(a,13e12.3)')'medium_M: ',chemo(OXYGEN+1:)%medium_M
if (full) then
	total_volume = Ve
	Caverage(MAX_CHEMO+1:2*MAX_CHEMO) = Ce
else
	Vcells = Ncells*Vcell_cm3   ! This is valid only if Ncells is the actual experimental cell count.
	Vm = total_volume - Vcells
	Vr = min(Vm,Ve)
	!write(nflog,'(a,4f8.4)') 'total_volume, Vcells, Vm, Vr: ',total_volume, Vcells, Vm, Vr 
	mass = (Vm - Vr)*Caverage(MAX_CHEMO+1:2*MAX_CHEMO) + Ve*Ce(:)
	total_volume = Vm - Vr + Ve + Vcells
	Caverage(MAX_CHEMO+1:2*MAX_CHEMO) = mass/(total_volume - Vcells)
endif
chemo(OXYGEN)%bdry_conc = Ce(OXYGEN)
do ichemo = GLUCOSE,DRUG_A
	C = Caverage(MAX_CHEMO+ichemo)
	Caverage(ichemo) = C
	chemo(ichemo)%Cmedium = C
	chemo(ichemo)%bdry_conc = C
enddo
do idrug = 1,1
	iparent = DRUG_A + 2*(idrug-1)
	do im = 0,1
		ichemo = iparent + im	
		chemo(ichemo)%Cmedium = Caverage(MAX_CHEMO+ichemo)
	enddo
enddo

call SetConstLevels	
write(nflog,'(a,3e12.3)') 'Const Cmedium: ',Caverage(MAX_CHEMO+1:MAX_CHEMO+2)
write(nflog,'(a,2e12.3)') 'Drug Cmedium:  ',Caverage(MAX_CHEMO+DRUG_A:MAX_CHEMO+DRUG_A+1)
write(nflog,*) 'glucose bdry_conc: ',chemo(GLUCOSE)%bdry_conc
t_lastmediumchange = istep*DELTA_T
medium_change_step = .true.
end subroutine

!-----------------------------------------------------------------------------------------
! Total cell O2 flux determines O2 at the well bottom Cex, but flux is determined by Cex.
! Use the current values for Cin and Cex to get the mass flux, which is equated to the
! flux that corresponds to the area, Kdiff and concentration gradient (Cbnd - Cex)/depth.
! Kd.A.(Cbnd - Cex)/d = Ncells.(Kin.Cex - Kout.Cin)
! => Cex = (A.Kd.Cbnd/d + Ncells.Kout.Cin)/(A.Kd/d + Ncells.Kin) 
!
! NOTE: currently WRONG for O2
!
!-----------------------------------------------------------------------------------------
subroutine SetConstLevels
integer :: ichemo, k, kcell, idrug, iparent, im
real(REAL_KIND) :: Kin, Kout, Kd, Cex, Cin, Cbnd, A, d, flux, Cin_prev, alpha, dC
real(REAL_KIND) :: tol = 1.0e-6

!ichemo = OXYGEN

do ichemo = 1,2
	if (chemo(ichemo)%constant .or. fully_mixed) then
		Cex = chemo(ichemo)%bdry_conc
		Cin = getCin(ichemo,Cex)
	else
		Kin = chemo(ichemo)%membrane_diff_in
		Kout = chemo(ichemo)%membrane_diff_out
		Cex = Caverage(MAX_CHEMO+ichemo)	! initial guess
		Cin = Caverage(ichemo)				! initial guess
		Cbnd = chemo(ichemo)%bdry_conc
		Kd = chemo(ichemo)%medium_diff_coef
		A = well_area
		d = total_volume/A
!		if (ichemo == OXYGEN) then
!			do k = 1,100
!				Cin_prev = Cin
!				Cex = (A*Kd*Cbnd/d + Ncells*Kout*Cin)/(A*Kd/d + Ncells*Kin)
!				Cin = getCin(ichemo,Cex)
!			    write(nflog,'(a,i4,3e15.6)') 'SetMediumOxygen: ',k,Cin,Cex,Cex-Cin
!				if (abs(Cin-Cin_prev)/Cin_prev < tol) exit
!			enddo
!		endif
	endif
	Caverage(ichemo) = Cin
	Caverage(MAX_CHEMO+ichemo) = Cex
	do k = 1,N1D
		alpha = real(k-1)/(N1D-1)
		chemo(ichemo)%Cmedium(k) = alpha*chemo(ichemo)%bdry_conc + (1-alpha)*Cex
	enddo
	do kcell = 1,nlist
		if (cell_list(kcell)%state == DEAD) cycle
		cell_list(kcell)%Cin(ichemo) = Cin
	!	do idrug = 1,2
	!		iparent = DRUG_A + 3*(idrug-1)
	!		do im = 0,2
	!			ichemo = iparent + im
	!			if (.not.chemo(ichemo)%present) cycle
	!			cell_list(kcell)%Cin(ichemo) = chemo(ichemo)%Cmedium(1)		! set IC conc to initial medium conc 
	!		enddo
	!	enddo
	enddo
	write(nflog,'(a,i4,2e12.3)') 'SetConstLevels: ichemo, Cex, Cin: ',ichemo,Cex,Cin
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Advance simulation through one big time step (DELTA_T)
! The concentration fields are first solved through NT_CONC subdivisions of DELTA_T,
! then the cell states are updated, which includes cell death and division.
! On either death or division cell positions are adjusted, and site concentrations
! (if necessary), and the ODE solver mappings in ODEdiff.
!-----------------------------------------------------------------------------------------
subroutine simulate_step(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: simulate_step  
use, intrinsic :: iso_c_binding
integer(c_int) :: res
integer :: kcell, site(3), hour, nthour, kpar=0
real(REAL_KIND) :: r(3), rmax, tstart, dt, dts, radiation_dose, diam_um, framp, area, diam
!integer, parameter :: NT_CONC = 6
integer :: i, ic, ichemo, ndt, iz, idrug, ityp, idiv, ndiv, Nmetabolisingcells, NpreCA, Nd, Nnew
integer :: nvars, ns, nphaseh(8), ph
real(REAL_KIND) :: dxc, ex_conc(120*CYCLE_PHASE+1)		! just for testing
real(REAL_KIND) :: DELTA_T_save, t_sim_0, SFlive
real(REAL_KIND) :: C_O2, HIF1, PDK1
type(metabolism_type), pointer :: mp
type(cell_type), pointer :: cp
integer :: phase_count(0:4), nG2
real(REAL_KIND) :: total, tadjust
real(REAL_KIND) :: fATM, fATR, fCP, dtCPdelay, dtATMdelay, dtATRdelay, ATM_DSB, DNA_rate
real(REAL_KIND) :: pATM_sum, pATR_sum, DSB_sum
real(REAL_KIND) :: SFtot, Pp, Pd, newSFtot, total_mitosis_time, V0
logical :: PEST_OK
logical :: ok = .true.
logical :: dbug

!write(nflog,*) 'istep,npar_uni,npar_rnor: ',istep,npar_uni,npar_rnor

cp => cell_list(1)
!call test_Jaiswal
!res = 1
!return

!kcell = 8
!cp=>cell_list(kcell)
!write(nflog,'(a,2i6,7f8.3)') 'cell phase,progress,DSB: ',kcell,cp%phase,cp%progress,cp%DSB(1:3,:)

if (istep < 2) call counter
mp => master_cell%metab

t_simulation = istep*DELTA_T	! seconds
!write(nfout,'(2i4,12f10.3)') istep,cp%phase,cp%progress,t_simulation/3600,cp%DSB(1:3,1),cp%Nmis

cp => cell_list(39)
!write(nfres,'(a,2i6,8f8.3)') 'cell 39: DSB: ',istep,cp%phase,cp%progress,t_simulation/3600,cp%DSB(1:2,:),cp%Nmis

!if (single_cell) write(nflog,*) 'Ncells_type: ',Ncells_type

!write(nflog,'(a,f8.2,2i6,f8.2)') 'start simulate_step: t_simulation: ',t_simulation/3600,istep,cp%phase,cp%progress
!call testmetab2
dbug = (istep < 0)
!Nmetabolisingcells = Ncells - (Ndying(1) + Ndying(2))
Nmetabolisingcells = Ncells - N_checkpoint
!if (Nmetabolisingcells == 0) then
!	call logger('# of metabolising cells = 0')
!    res = 0
!    return
!endif
if (Ncells == 0) then
	call logger('# of cells = 0')
    res = 2
    return
endif
if (limit_stop) then
	call logger('Spheroid size limit reached')
	res = 6
	return
endif

nthour = 3600/DELTA_T

if (istep == -100) then
	stop
endif

if (ngaps > 200) then
	call squeezer
endif

!if (use_metabolism) then

drug_gt_cthreshold = .false.

#if 0
!if (medium_change_step .or. (chemo(DRUG_A)%present .and..not.DRUG_A_inhibiter)) then
!if (medium_change_step .or. (chemo(DRUG_A)%present .and..not.use_inhibiter)) then
!    ndiv = 6
!else
    ndiv = 1
!endif
dt = DELTA_T/ndiv
dts = dt/NT_CONC
DELTA_T_save = DELTA_T
DELTA_T = dt
t_sim_0 = t_simulation
do idiv = 0,ndiv-1
    t_simulation = t_sim_0 + idiv*DELTA_T
    
    if (dbug) write(nflog,*) 'Solver'
    do it_solve = 1,NT_CONC
	    tstart = (it_solve-1)*dts
        !	t_simulation = (istep-1)*DELTA_T + tstart
	    t_simulation = t_simulation + tstart
	    call Solver(it_solve,t_simulation,dts,Ncells,ok)
	    if (.not.ok) then
		    res = 5
		    return
	    endif
    enddo	! end it_solve loop
    if (dbug) write(nflog,*) 'did Solver'
!	if (use_metabolism) then
!		do ityp = 1,Ncelltypes
!			HIF1 = mp%HIF1
!			C_O2 = chemo(OXYGEN)%cmedium(1)
!!			write(*,'(a,2e12.3)') 'before: HIF1,C_O2: ',HIF1,C_O2 
!			call analyticSetHIF1(C_O2,HIF1,DELTA_T)
!!			write(*,'(a,e12.3)') 'after: HIF1: ',HIF1
!			mp%HIF1 = HIF1
!			PDK1 = mp%PDK1
!			call analyticSetPDK1(HIF1,PDK1,dt)
!			mp%PDK1 = PDK1
!		enddo
!!	endif
!	!write(nflog,*) 'did Solver'
!    call GlutamineDecay
    if (Ndrugs_used > 0) then
        call CheckDrugConcs
        call CheckDrugPresence
    endif
    
!    if (t_irradiation >= 0) call GrowCells(DELTA_T,t_simulation,ok)

    if (.not.ok) then
	    res = 3
	    return
    endif
enddo	! end idiv loop
DELTA_T = DELTA_T_save
#endif

medium_change_step = .false.

radiation_dose = 0
if (use_treatment) then     ! now we use events
	call treatment(radiation_dose)
endif
if (use_events) then
	call ProcessEvent(radiation_dose)
endif
if (radiation_dose >= 0) then
	write(logmsg,'(a,f6.1)') 'Radiation dose: ',radiation_dose
	call logger(logmsg)
! could make calls to set_divide_volume and SetInitialCellCycleStatus here
! before Irradiation
!    write(nflog,*) 'before npar_uni, npar_rnor = ',npar_uni,npar_rnor
    do kcell = 1,Ncells
        cp => cell_list(kcell)
        V0 = Vdivide0/2
        call set_divide_volume(cp,V0)
        call SetInitialCellCycleStatus(kcell,cp)
    enddo
!    write(nflog,*) 'after npar_uni, npar_rnor = ',npar_uni,npar_rnor
	call Irradiation(radiation_dose, ok)
	if (.not.ok) then
		res = 3
		return
    endif
endif

if (t_irradiation >= 0) call GrowCells(DELTA_T,t_simulation,ok)

tottotDSB = 0
do kcell = 1,nlist
    cp => cell_list(kcell)
    tottotDSB = tottotDSB + sum(cp%DSB)
enddo
tottotDSB = 0
    
res = 0

call getNviable

kcell = 1
cp => cell_list(kcell)
if (use_synchronise .and. .false.) then
    call get_slowdown_factors(cp,fATM,fATR)
    fCP = fATM*fATR
    dt = DELTA_T
    dtCPdelay = dt*(1 - fCP)
    if (2-fATM-fATR == 0) then
        dtATMdelay = 0
        dtATRdelay = 0
    else
        dtATMdelay = dtCPdelay*(1 - fATM)/(2 - fATM - fATR)
        dtATRdelay = dtCPdelay*(1 - fATR)/(2 - fATM - fATR)
    endif
    if (isnan(dtATMdelay)) then
        write(*,*) 'dtATMdelay is NAN'
        write(*,*) dtCPdelay,fATM,fATR,2-fATM-fATR
        stop
    endif
    tCPdelay = tCPdelay + dtCPdelay
    tATMdelay = tATMdelay + dtATMdelay
    tATRdelay = tATRdelay + dtATRdelay 
!    write(nfphase,'(2i6,11f8.3)') istep,cp%phase,cp%progress,t_simulation/3600,cp%pATM,cp%pATR,fATM,fATR,fCP,tCPdelay/3600,tATMdelay/3600,tATRdelay/3600
    ATM_DSB = sum(cp%DSB(NHEJslow,:)) + sum(cp%DSB(HR,:))   ! complex DSB
    if (phase_log) write(nfphase,'(2i6,11f8.3)') istep,cp%phase,cp%progress,t_simulation/3600,ATM_DSB,cp%pATM
    if (cp%phase == S_phase) stop
    if (cp%phase == dividing) then
        write(*,*) 'Reached mitosis'
        res = 1
        return
    endif
endif

!if (output_DNA_rate .and. (Nirradiated > 0)) then
!    call show_S_phase_statistics()
!endif
if (compute_cycle) then
    call get_phase_distribution(phase_count)
    total = sum(phase_count)
    phase_dist = 100*phase_count/total
    tadjust = event(1)%time/3600    ! if the RADIATION event is #1
!    write(nflog,'(a,i6,5f8.2)') 'phase distribution: ',istep,phase_dist
!    write(nflog,'(a,6i6,f10.0)') 'phase_count: ',istep,phase_count,total
!    write(*,'(a,f8.3,4i8,4f8.3)') 'hour,count, phase_dist: ',real(istep)/nthour,phase_count(1:4),phase_dist(1:4)
!    write(nflog,'(a,f8.3,i8,f8.3)') 'hour, count, M%: ',real(istep)/nthour - tadjust,phase_count(M_phase),phase_dist(M_phase)
endif
if (compute_cycle .or. output_DNA_rate) then
    if (next_phase_hour > 0) then  ! check if this is a phase_hour
        if (real(istep)/nthour >= phase_hour(next_phase_hour)) then   ! record phase_dist
            write(*,*) 'Reached phase hour: ',next_phase_hour,phase_hour(next_phase_hour)
            if (next_phase_hour <= 9) then
                write(*,'(a,4i8)') 'count: ',phase_count(1:4)
                write(nflog,*) 'Reached phase hour: ',next_phase_hour,phase_hour(next_phase_hour)
                write(nflog,'(a,4i8)') 'count: ',phase_count(1:4)
                write(nflog,'(a,4f8.3)') 'dist: ',phase_dist(1:4)
!                write(nflog,'(a,4f8.3)') '%dist: ',100*phase_dist(1:4)/sum(phase_dist(1:4))
            endif
	        if (compute_cycle) then
!	            call get_phase_distribution(phase_count)
!	            total = sum(phase_count)
!	            phase_dist = 100*phase_count/total
                recorded_phase_dist(next_phase_hour,1:4) = 100*phase_dist(1:4)/sum(phase_dist(1:4))
	        endif
	        if (output_DNA_rate) then
	            call get_DNA_synthesis_rate(DNA_rate)
	            recorded_DNA_rate(next_phase_hour) = DNA_rate
                write(*,*) 'in simulate_step: DNA_rate: ',DNA_rate
	        endif
            next_phase_hour = next_phase_hour + 1
            if (next_phase_hour > nphase_hours) next_phase_hour = 0
        endif
    endif
endif

	if (use_synchronise) then
	    pATM_sum = 0
	    pATR_sum = 0
	    DSB_sum = 0
	    nG2 = 0
	    do kcell = 1,nlist
            cp => cell_list(kcell)
            if (cp%phase == G2_phase) then
                nG2 = nG2+1
                pATM_sum = pATM_sum + cp%pATM
                pATR_sum = pATR_sum + cp%pATR
                DSB_sum = DSB_sum + sum(cp%DSB)
            endif
        enddo
!        write(nfphase,'(a,f8.3,2f8.4,f8.1)') 'hour, (G2) ave pATM, pATR, DSB: ',istep*DELTA_T/3600.,pATM_sum/nG2, pATR_sum/nG2, DSB_sum/nG2
!        write(*,'(a,i4,f6.3,i4)') 'phase, progress: ',cp%phase, cp%progress, next_phase_hour
!        write(nfphase,'(2i4,2f6.1,2f6.3)') hour,cp%phase,tnow/3600,cp%G2_time/3600,cp%V/Vdivide0,cp%divide_volume/Vdivide0
    endif

if (dbug .or. mod(istep,nthour) == 0) then
    hour = istep/nthour
    nphaseh = 0
    do kcell = 1,nlist
        cp => cell_list(kcell)
        i = cp%phase
        nphaseh(i) = nphaseh(i) + 1
    enddo
    nphase(hour,:) = nphaseh
	ntphase = nphaseh + ntphase
!	write(nflog,*)
	write(nflog,'(a,i6,i4,4(a,i8))') 'istep, hour: ',istep,hour,' Nlive: ',Ncells   !, ' N reached mitosis: ',NPsurvive    ,' Napop: ',Napop    !, &
	if (.not. single_cell) write(*,'(a,i6,i4,4(a,i8))') 'istep, hour: ',istep,hour,' Nlive: ',Ncells   !, ' N reached mitosis: ',NPsurvive    ,' Napop: ',Napop    !, &
!    write(nflog,*) 'npar_uni,npar_rnor: ',npar_uni,npar_rnor
    call get_phase_distribution(phase_count)
    total = sum(phase_count(1:4))
    phase_dist = 100*phase_count/total
!    write(*,'(a,4i8,4f7.1)') 'count, phase_dist: ',phase_count(1:4),phase_dist(1:4)
!    write(nflog,'(a,4i8,4f7.1)') 'count, phase_dist: ',phase_count(1:4),phase_dist(1:4)
!	if (single_cell) call medras_compare()
!	write(nfphase,'(a,2f8.3)') 'S-phase k1, k2: ', K_ATR(2,1),K_ATR(2,2)
!	if (output_DNA_rate) then
!	    call get_DNA_synthesis_rate(DNA_rate)
!	endif
	if (use_synchronise .and. .false.) then
	    pATM_sum = 0
	    pATR_sum = 0
	    do kcell = 1,nlist
            cp => cell_list(kcell)
            pATM_sum = pATM_sum + cp%pATM
            pATR_sum = pATR_sum + cp%pATR
        enddo
        if (phase_log) write(nfphase,'(a,i4,2f8.4)') 'hour, ave pATM, pATR: ',hour,pATM_sum/nlist, pATR_sum/nlist
!        write(*,'(a,i4,f6.3,i4)') 'phase, progress: ',cp%phase, cp%progress, next_phase_hour
!        write(nfphase,'(2i4,2f6.1,2f6.3)') hour,cp%phase,tnow/3600,cp%G2_time/3600,cp%V/Vdivide0,cp%divide_volume/Vdivide0
    endif
    
    total = 0
    do kcell = 1,nlist
        cp => cell_list(kcell)
        total = total + cp%totDSB0
    enddo    
endif

! write(nflog,'(a,f8.3)') 'did simulate_step: time: ',wtime()-start_wtime
!kcell = 5388
!cp => cell_list(kcell)
!write(nflog,'(a,3i8)') 'cell, state, phase: ',kcell,cp%state,cp%phase
!write(*,'(a,3i8)') 'cell, state, phase: ',kcell,cp%state,cp%phase

istep = istep + 1
overstepped = (istep == maxhours*nthour)
if (overstepped) then
    write(*,*) 'overstepped the mark'
    call nondivided()
!    stop
    call completed
    res = 1
    return
endif

!write(*,*) 'end simulate_step: t_simulation: ',t_simulation
!call averages
! Need the phase_dist check in case NPsurvive = Nirradiated before phase_hours

! Check for completion of the run
if (compute_cycle .or. output_DNA_rate) then
    if (next_phase_hour == 0) then
        call completed
        res = 1
    endif
    return
endif
PEST_OK = .true.
if (use_PEST) then  ! check that all phase distributions have been estimated
    PEST_OK = (next_phase_hour == 0)
endif
    
!write(*,*) 'NPsurvive: ',NPsurvive, (Nirradiated - Napop)
if (is_radiation .and. (NPsurvive >= (Nirradiated - Napop - Nmitotic)) .and. PEST_OK) then  !!! needs to change
    ! getSFlive computes the average psurvive for all cells that reach mitosis,
    ! which number NPsurvive = Nirradiated - Napop.
    call getSFlive(SFlive)
    SFtot = SFlive*(Nirradiated - Napop)
    write(*,*)
    write(nflog,'(a,3i6,e12.3)') 'NPsurvive,Nirradiated,Napop,SFtot: ',NPsurvive,Nirradiated,Napop,SFtot
    write(nflog,'(a,e12.3)') 'Unadjusted SFave = SFtot/NPsurvive: ',SFtot/NPsurvive
    write(nflog,'(a,e12.3)') 'SFave including apoptosis killing: ',SFtot/Nirradiated
    if (include_daughters) then
    ! To adjust SFlive to replace a cell that reached mitosis before CA with its daughter cells.
    ! In order to compare simulated SFave with SF determined by experimental CA (clonogenic analysis),
    ! SFave needs to be calculated as the average of cells that make it to CA.
        write(nflog,*) 'Accounting for daughters: from NPsurvive: ',NPsurvive
        newSFtot = 0
        Nnew = 0
        NpreCA = 0
        Nd = 0
        total_mitosis_time = 0
        do kcell = 1,nlist
            cp => cell_list(kcell)
            if (cp%state == DEAD) cycle
            total_mitosis_time = total_mitosis_time + cp%mitosis_time
!            if (cp%psurvive > 1.0e-30) then     ! was 1.0E-8
!            write(*,'(a,2e12.3)') 'cp%mitosis_time, CA_time: ',cp%mitosis_time, CA_time
                if (cp%mitosis_time < CA_time_h*3600) then
                    NpreCA = NpreCA + 1
                    Pp = cp%psurvive
                    Pd = 1 - sqrt(1.0 - Pp)   ! Pd = psurvive for the 2 daughters: Pp = 2Pd - Pd^2
! Try this !
!    Pd = Pp/2  No good - dose = 0.01 --> SFave = 0.54
                    NPsurvive = NPsurvive + 1
                    SFtot = SFtot - Pp + 2*Pd
                    newSFtot = newSFtot + 2*Pd
                    Nnew = Nnew + 2
    !                if (kcell <= 100) then
    !                    write(*,'(a,i6,2e12.3)') 'daughters: kcell, Pp, Pd: ',kcell,Pp,Pd
    !                    write(*,'(a,2i3,2f8.3)') 'kcell,phase0,Pp,2*Pd: ',kcell,cp%phase0,Pp,2*Pd
    !                endif
!                    if (cp%phase0 < 4 .and. Pd < 0.3) then
!                        Nd = Nd + 2
!                        write(*,'(a,2i6,f8.3)') 'daughters: ', kcell, cp%phase0,Pd
!                    endif
                !else
                !    write(nflog,'(a,2f8.2)') 'mitosis_time, CA_time: ',cp%mitosis_time/3600, CA_time/3600
                !    newSFtot = newSFtot + cp%psurvive
                !    Nnew = Nnew + 1
!!                    if (cp%phase0 < 4 .and. cp%psurvive < 0.4) write(*,'(a,2i6,f8.3)') 'psurvive: ',kcell, cp%phase0,cp%psurvive
                endif
!            else
!                write(nflog,'(a,e12.3)') 'psurvive: ',cp%psurvive
!            endif
        enddo
        write(nflog,*) 'to NPsurvive: ',NPsurvive
        write(*,'(a,3i6,2f10.3)') 'NpreCA, Nd, Nnew, newSFtot, newSFave: ',NpreCA,Nd,Nnew,newSFtot,newSFtot/Nnew
        write(*,'(a,3i6,f8.4)') 'nlist,Nirradiated,NPsurvive, new SFave: ',nlist,Nirradiated,NPsurvive,newSFtot/nlist
        write(*,'(a,i8,f8.2)') 'Average mitosis_time: ',nlist,total_mitosis_time/(3600*nlist)
        write(*,'(a,f8.2,i6)') 'CA_time_h, NpreCA: ',CA_time_h,NpreCA
    endif
    if (NPsurvive > 0) then
        SFave = SFtot/(NPsurvive + Napop)
    else
        SFave = 0
    endif
    write(*,'(a,i6,2e12.3)') 'NPsurvive,SFlive,SFtot: ',NPsurvive,SFlive,SFtot
!    if (use_Napop) then
!        SFave = ((Nirradiated - Napop)/Nirradiated)*SFlive
!    else
!        SFave = SFlive
!    endif
    write(*,*)
    write(nflog,'(a,f8.2)') 'CA_time_h: ',CA_time_h
    write(nflog,'(a,2i6)') 'Npsurvive, Napop: ',Npsurvive,Napop
    write(logmsg,'(a,e12.4,f8.3)') 'SFave,log10(SFave): ',SFave,log10(SFave)
    call logger(logmsg)

    
    call completed
!    write(nflog,'(a)') 'nphase:'
!    do hour = 0,istep/nthour
!        write(nflog,'(i3,8i8)') hour,nphase(hour,:)
!    enddo
    res = 1
endif

end subroutine

!-----------------------------------------------------------------------------------------
! When the mark is overstepped, locate cells that have not reached mitosis
!-----------------------------------------------------------------------------------------
subroutine nondivided
integer :: kcell, n
type(cell_type), pointer :: cp

n = 0
do kcell = 1,nlist
    cp => cell_list(kcell)
    if (cp%state /= DEAD .and. cp%psurvive < 0) then
        n = n+1
        write(*,'(a,i6,i3,3f8.4)') 'nondivided: kcell, phase, kt2cc,ke2cc,kcc2a: ',kcell, cp%phase,cp%kt2cc,cp%ke2cc,cp%kcc2a
        write(*,*) 'fg: ',cp%fg
        write(nflog,'(a,i6,i3,3f8.4)') 'nondivided: kcell, phase, kt2cc,ke2cc,kcc2a: ',kcell, cp%phase,cp%kt2cc,cp%ke2cc,cp%kcc2a
        write(nflog,*) 'fg: ',cp%fg
    endif
enddo
write(*,*) 'Total nondivided: ',n
end subroutine

#if 0
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine medras_compare()
type(cell_type), pointer :: cp
real(8) :: totDSB, Paber, Pmit, Psurvive

cp => cell_list(1)
totDSB = sum(cp%DSB(1:3))
write(nflog,'(a,3f8.3)') 'N_DSB, totNmis, Nlethal: ',totDSB,cp%Nlethal/Klethal,cp%Nlethal
Paber = exp(-cp%Nlethal)
Pmit = exp(-mitRate*totDSB)
Psurvive = Pmit*Paber  
write(nflog,'(a,3f8.4)') 'Paber, Pmit, Psurvive: ',Paber, Pmit, Psurvive
end subroutine
#endif

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!function getSFlive result(SF)
subroutine getSFlive(SF)
real(REAL_KIND) :: SF
real(REAL_KIND) :: sfsum, totDSB, total, total0
integer :: kcell, n
type(cell_type), pointer :: cp

n = 0
sfsum = 0
total = 0
total0 = 0
do kcell = 1,nlist
    cp => cell_list(kcell)
    if (cp%state == DEAD) cycle
    if (cp%totDSB0 <= 0) then
        cycle  ! this cell was not irradiated - must be a daughter cell (can't happen)
    endif
    n = n+1
    sfsum = sfsum + cp%Psurvive
    totDSB = sum(cp%DSB)
!    if (abs(totDSB - cp%totDSB0) > 1) then
!        write(nflog,'(a,i6,2f8.0)') 'getSFlive: totDSB,totDSB0: ',kcell,totDSB,cp%totDSB0
!    endif
!    if (kcell <= 10) write(nflog,'(a,i6,2f9.1)') 'totDSB0,totDSB: ',kcell,cp%totDSB0,totDSB
    total = total + totDSB
    total0 = total0 + cp%totDSB0
    if (cp%Psurvive < 0) then
        write(*,*) 'ERROR: getSFlive: Psurvive: ',kcell,cp%Psurvive
        stop
    endif
enddo
!write(nflog,*) 'getSFlive: n, Nirradiated: ',n, Nirradiated
! n = Nirradiated - Napop
SF = sfsum/n
!write(nflog,'(a,2f10.1)') 'getSFlive: initial,final DSB totals: ',total0,total
!write(*,'(a,i6,4e12.3)') 'n,Ave totDSB,Pmit,Nlethal,Paber: ',n,totDSB/n, totPmit/n, totNlethal/n, totPaber/n
!end function
end subroutine

!-----------------------------------------------------------------------------------------
! The selection of outputs needs to be made an input parameter.
! Choices:
! SF + distribution
! SF only
! distribution only
!
! The first two can be lumped together - if distribution is not wanted it will not be read.
! I.e. need only use_SF
! For now assume always use_SF = true, because fitting without SF is no good.
! For PEST runs, use_SF corresponds to 'M' fitting
!-----------------------------------------------------------------------------------------
subroutine completed
!real(REAL_KIND) :: fract(0:4)
integer :: kcell, ph, nir(4), nmitosis,nsum, kcellmax, i, j, k, ityp
real(REAL_KIND) :: sftot_phase(4), sfsum, sfmax
integer :: tadjust
type(cycle_parameters_type), pointer :: ccp
type(cell_type), pointer :: cp
logical :: only_M_phase = .false.
logical :: PDS4 = .false.
real(REAL_KIND) :: dt, phi, PDS4_M(3) = [0.191, 0.414286, 0.732812]
real(REAL_KIND) :: normalised_phase_dist(60,0:4)   
REAL(REAL_KIND) :: ave(15), SFMave

if (overstepped) then
    SFave = 1.0E-6
    goto 99
endif
if (compute_cycle) then
    write(nflog,*) 'Completed compute cycle'
    write(*,*) 'Completed compute cycle'
    if (use_synchronise) then
        tadjust = 0
    else
        tadjust = event(1)%time/3600
    endif
!    recorded_phase_dist(1,1:NP) = phdist0      ! ?????????????
    do i = 1,nphase_hours
        write(nflog,'(f6.1,4x,4f8.3)') phase_hour(i)-tadjust,recorded_phase_dist(i,1:4)
        write(*,'(f6.1,4x,4f8.3)') phase_hour(i)-tadjust,recorded_phase_dist(i,1:4)
        if (PDS4) then
            phi = phi + (recorded_phase_dist(i,4) - PDS4_M(i))**2
        endif
    enddo
    write(*,*) 'wrote recorded_phase_dist'
    write(nflog,*) 'wrote recorded_phase_dist'
    write(*,*) 'post-mitosis G1 counts: ng11, ng12: ',ng11,ng12
    if (PDS4) then 
        write(nflog,'(a)') '    hour    expt   model   error'
        write(*,'(a)') '    hour    expt   model   error'
        do i = 1,nphase_hours
            write(nflog,'(f8.1,3f8.4)') phase_hour(i)-tadjust,PDS4_M(i),recorded_phase_dist(i,4),recorded_phase_dist(i,4) - PDS4_M(i)
            write(*,'(f8.1,3f8.4)') phase_hour(i)-tadjust,PDS4_M(i),recorded_phase_dist(i,4),recorded_phase_dist(i,4) - PDS4_M(i)
        enddo
        write(nflog,'(a,f6.3)') '    phi: ',phi
        write(*,'(a,f6.3)') '    phi: ',phi
    endif
    
    if (nphase_hours > 0) then
        if (only_M_phase) then
            write(nfres,'(20e15.6)') (recorded_phase_dist(i,4),i=1,nphase_hours)
        else
            if (normalise) then
                ityp = 1
                ccp => cc_parameters(ityp)
                control_ave(1) = 100*ccp%f_G1
                control_ave(2) = 100*ccp%f_S
                control_ave(3) = 100*ccp%f_G2
                control_ave(4) = 100*ccp%f_M

                write(*,*) 'Normalising PDs'
                write(nflog,*) 'Normalising PDs'
                write(nflog,'(a,4f8.3)') 'control: ',control_ave(1:4)
                dt = 0.5
                do i = 1,nphase_hours
                    do j = 1,4
                        normalised_phase_dist(i,j) = recorded_phase_dist(i,j)/control_ave(j)
                    enddo
                    write(nflog,'(f6.1,4f8.3)') phase_hour(i)-tadjust,normalised_phase_dist(i,1:4)
                enddo
                write(*,*) 'write PEST output'
                if (G2M_only) then
                    if (expt_tag == "CA-135") then
                        write(nfres,'(20f8.5)') (normalised_phase_dist(i,3:4),i=1,nphase_hours)
                    elseif (expt_tag == "CC-11 ") then
                        write(nfres,'(20f8.5)') (normalised_phase_dist(i,3:3),i=1,nphase_hours)
                    elseif (expt_tag == "CC-13 ") then
                        write(nfres,'(20f8.5)') (normalised_phase_dist(i,4),i=1,nphase_hours)
                    endif
                else
                    write(*,'(a,a,i6)') 'expt_tag,nphase_hours: ',expt_tag,nphase_hours
                    if (expt_tag == "CA-135") then
                        write(nfres,'(20f8.5)') (normalised_phase_dist(i,1:4),i=1,nphase_hours)
                        write(nflog,'(20f8.5)') (normalised_phase_dist(i,1:4),i=1,nphase_hours)
                    elseif (expt_tag == "CC-11 ") then
                        write(nfres,'(20f8.5)') (normalised_phase_dist(i,1:3),i=1,nphase_hours)
                        write(nflog,'(20f8.5)') (normalised_phase_dist(i,1:3),i=1,nphase_hours)
                    elseif (expt_tag == "CC-13 ") then
                        write(nfres,'(20f8.5)') (normalised_phase_dist(i,4),i=1,nphase_hours)
                        write(nflog,'(20f8.5)') (normalised_phase_dist(i,4),i=1,nphase_hours)
                    endif
                endif
                if (expt_tag == "PDSN0G") then
                    write(nfres,'(20f8.5)') (normalised_phase_dist(i,1:4),i=1,nphase_hours)
                elseif (expt_tag == "PDSN2G") then
                    write(nfres,'(20f8.5)') (normalised_phase_dist(i,4),i=1,3), &
                                            (normalised_phase_dist(i,1:4),i=4,nphase_hours)
                elseif (expt_tag == "PDSN6G") then
                    write(nfres,'(20f8.5)') (normalised_phase_dist(i,1:4),i=1,nphase_hours)
                endif
                write(nflog,'(20f10.5)') (normalised_phase_dist(i,1:4),i=1,nphase_hours)
                write(nflog,'(20f10.5)') control_ave(1:4)
            else
                write(*,*) 'Not normalising PDs'
                write(nfres,'(20e15.6)') (recorded_phase_dist(i,1:4),i=1,nphase_hours)
            endif
        endif
!        do i = 1,nphase_hours
!            write(nfres,'(f6.1,4x,4f6.1)') phase_hour(i)-tadjust,recorded_phase_dist(i,1:4)
!            write(*,'(f6.1,4x,4f6.1)') phase_hour(i)-tadjust,recorded_phase_dist(i,1:4)
!        enddo
    endif
    write(*,'(a,3f8.3)') 'Average G1, S, G2 CP delays: ',totG1delay/nG1delay,totSdelay/(3600.*nSdelay),totG2delay/nG2delay
    write(*,'(a,3f8.3)') 'Average G2 delay: pATM, pATR, total: ',G2thsum/NG2th,sum(G2thsum)/NG2th
    if (phase_log) then
        write(nfphase,'(a,3f8.3)') 'Average G2 delay: pATM, pATR, total: ',G2thsum/NG2th,sum(G2thsum)/NG2th
        write(nfphase,'(a,i6)') 'Number of cells averaged: ',NG2th
        if (npet > 0) write(*,'(a,i6,f6.3)') 'Average phase time: ',npet,phase_exit_time_sum/(3600*npet)
    endif
!    return
endif
if (output_DNA_rate) then
    write(nflog,*) 'Completed'
    write(*,*) 'Completed'
    tadjust = IR_time_h
    if (nphase_hours > 0) then
        write(*,*) 'write DNA_rate'
        write(nfres,'(20e15.6)') (recorded_DNA_rate(i),i=1,nphase_hours)
        do i = 1,nphase_hours
            write(nflog,'(f6.2,4x,4f6.3)') phase_hour(i)-tadjust,recorded_DNA_rate(i)
            write(*,'(f6.2,4x,4f6.3)') phase_hour(i)-tadjust,recorded_DNA_rate(i)
        enddo
    endif
    return
endif

! Look at average survival by IR phase
nir = 0
sftot_phase = 0
!nsum = 0
!sfsum = 0
sfmax = 0
do kcell = 1,nlist
    cp => cell_list(kcell)
    if (cp%state == DEAD) cycle
    if (cp%totDSB0 <= 0) cycle
!    nsum = nsum + 1
!    sfsum = sfsum + cp%Psurvive
!    if (cp%phase0 <= G1_checkpoint) then
!        ph = 1
!        if (cp%Psurvive > sfmax) then
!            sfmax = cp%Psurvive
!            kcellmax = kcell
!        endif
!    elseif (cp%phase0 <= S_checkpoint) then
!        ph = 2
!    elseif (cp%phase0 <= G2_checkpoint) then
!        ph = 3
!    else
!        ph = 4
!    endif
    if (cp%phase0 == G1_phase) then
        ph = 1
        if (cp%Psurvive > sfmax) then
            sfmax = cp%Psurvive
            kcellmax = kcell
        endif
    elseif (cp%phase0 == S_phase) then
        ph = 2
    elseif (cp%phase0 == G2_phase) then
        ph = 3
    else
        ph = 4
    endif
    nir(ph) = nir(ph) + 1
    sftot_phase(ph) = sftot_phase(ph) + cp%Psurvive
enddo
!nir = max(nir,1)       ! this was an error
nmitosis = sum(nir)
write(*,'(a,4i6)') 'nir: ',nir
!write(*,'(a,4f10.6)') 'Average SF by phase: ',sftot/nir
!if (use_synchronise) then
!    if (phase_log) write(nfphase,'(a,4f8.4)') 'Average SF by phase: ',sftot/nir
!endif
write(nflog,'(a,6f12.3)') 'totPmit, totPaber, tottotDSB: ',totPmit, totPaber, tottotDSB
write(*,'(a,i6,5f11.1)') 'Nmitosis, totPmit, totPaber, tottotDSB: ',int(Nmitosis),totPmit, totPaber, tottotDSB
write(*,'(a,6e12.3)') 'totPaber: ',totPaber
write(nflog,'(a,4f10.5)') 'SFtot_phase: ',SFtot_phase
write(nflog,'(a,i6,f10.5)') 'Nmitosis,SFtot: ',Nmitosis,sum(SFtot_phase)
! adjust for pre-rep doubling of misjoins
totNmisjoins(1) = 2*totNmisjoins(1)
if (use_equal_mitrates) then
    write(nflog,*) 
    write(nflog,*) '!!!! use_equal_mitrates = true !!!!'
    write(nflog,*)
endif
write(nflog,'(a,7f9.3)') 'Ave (pre, post) NDSB, Nmisjoins: ', &
    totNDSB/nmitosis,totNmisjoins/nmitosis,sum(totNmisjoins)/nmitosis
write(*,'(a,7f9.3)') 'Ave (pre, post) NDSB, Nmisjoins: ', &
    totNDSB/nmitosis,totNmisjoins/nmitosis,sum(totNmisjoins)/nmitosis

if (.false.) then   ! make this true to write BBB lines
    ! Averages
    SFMave = 0
    ave = 0
    do kcell = 1,nlist
        cp => cell_list(kcell)
        SFMave = SFMave + cp%Psurvive
        do i = 1,3
            do j = 1,2
                k = (i-1)*2 + j
                ave(k) = ave(k) + cp%DSB0(i,j)
            enddo
        enddo
        ave(7) = ave(7) + cp%t_mitosis
    enddo
    SFMave = SFMave/nlist
    ave = ave/nlist
    write(*,'(a,2f7.3,10f7.2)') 'BBB: ',log10(SFMave),log10(SFave),totNDSB/nmitosis,sum(totNmisjoins)/nmitosis,ave(1:7)
endif
!write(*,'(a,7f8.2)') 'DSB0,tIR: ',ave(1:7)
!write(logmsg,'(a,8f8.3)') 'phase_dist:      ',phase_dist(0:3)
!call logger(logmsg)
!write(*,'(a,f8.3)') 'Final average pATM: ',ATMsum/nirradiated
!write(*,'(a,f8.3)') 'Final average pATR: ',ATRsum/nirradiated

!write(nflog,'(a,2f6.2)') 'Distributions of Nlethal, totDSB: ddist_Nlethal, ddist_totDSB: ',ddist_Nlethal,ddist_totDSB
!do i = 1,NMDIST
!    write(nflog,'(i6,2(f8.3,f8.4,2x))') i,ddist_Nlethal*(i-0.5),count_Nlethal(i)/nmitosis,ddist_totDSB*(i-0.5),count_totDSB(i)/nmitosis
!enddo

99 continue
if (use_synchronise) call G2_time_distribution()
if (use_PEST) then
    if (use_SF) then
        write(nfres,'(e15.6)') log10(SFave)
    endif
!    if (nphase_hours > 0) then
!        write(*,*) 'write recorded_phase_dist'
!        write(nfres,'(20e15.6)') (recorded_phase_dist(i,1:4),i=1,nphase_hours)
!    endif
endif

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine G2_time_distribution
integer :: kcell, n, k
integer :: cnt(40)
real(REAL_KIND) :: t
type(cell_type), pointer :: cp

n = 0
cnt = 0
do kcell = 1,nlist
    cp => cell_list(kcell)
    if (cp%state == DEAD) cycle
    t = cp%G2_time/3600
    if (t == 0) cycle
    n = n+1
    k = t + 1
    cnt(k) = cnt(k) + 1
enddo
write(*,*) 'G2_time distribution (h): n: ',n
write(nflog,*) 'G2_time distribution (h): n: ',n
do k = 1,40
    if (cnt(k) > 0) then
        write(*,'(i2,a,i2,i6,f8.3)') k-1, '-', k, cnt(k), cnt(k)/real(n)
        write(nflog,'(i2,a,i2,i6,f8.3)') k-1, '-', k, cnt(k), cnt(k)/real(n)
    endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! This is for drm_monolayer_deriv.exe
!-----------------------------------------------------------------------------------------
subroutine getResults(SF, dist)
!DEC$ ATTRIBUTES DLLEXPORT :: getResults
real(REAL_KIND) :: SF
integer :: dist(:)

SF = SFave
dist = phase_dist
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine averages
real(REAL_KIND) :: ave_V, ave_dVdt, ave_fg
integer :: kcell, n
type(cell_type), pointer :: cp

ave_V = 0
ave_dVdt = 0
ave_fg = 0
n = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) continue
	n = n+1
	ave_V = ave_V + cp%V
	ave_dVdt = ave_dVdt + cp%dVdt
!	ave_fg = ave_fg + cp%fg
enddo
write(nflog,'(a,3e12.3)') 'averages: V,dVdt,fg: ',ave_V/n,ave_dVdt/n
end subroutine

!-----------------------------------------------------------------------------------------
! 
!-----------------------------------------------------------------------------------------
subroutine show_metabolism(kcell)
integer :: kcell
type(cell_type), pointer :: cp
type (metabolism_type), pointer :: mp
real(REAL_KIND) :: Cin(MAX_CHEMO)
integer :: res
!real(REAL_KIND) :: HIF1, G_rate, PP_rate, P_rate
!real(REAL_KIND) :: L_rate, A_rate, I_rate, O_rate

cp =>cell_list(kcell)
mp => cp%metab
Cin = Caverage(1:MAX_CHEMO)
!write(*,'(a,3f8.4)') 'O2, glucose, lactate: ',Cin(1:3) 
call get_metab_rates(cp, Cin, 1.0d0, res)
return

write(*,'(a,i2,3e12.3)') 'phase, V: ',cp%phase,cp%V		!I2Divide,Itotal,mp%Itotal,mp%I2Divide
write(*,'(a,3e11.3)') 'G_rate, P_rate, O_rate: ',mp%G_rate, mp%P_rate, mp%O_rate
write(*,'(a,4e11.3)') 'L_rate, A_rate, I_rate: ',mp%L_rate, mp%A_rate, mp%I_rate
write(*,'(a,4f8.4)') 'O2, glucose, lactate, H: ',Caverage(OXYGEN),Caverage(GLUCOSE),Caverage(LACTATE),mp%HIF1
if (mp%A_rate < r_Ag) then
	write(*,*) 'Not growing'
endif
write(*,*)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine showcell(kcell)
integer :: kcell
type(cell_type), pointer :: cp

cp => cell_list(kcell)
write(nflog,'(a,i6,4e12.3)') 'kcell, volume, divide_volume, dVdt, divide_time: ', &
                kcell, cp%V, cp%divide_volume, cp%dVdt, cp%divide_time
end subroutine

!-----------------------------------------------------------------------------------------
! Average volumes etc
!-----------------------------------------------------------------------------------------
subroutine showcells
integer :: kcell, n
real(REAL_KIND) :: Vsum,divVsum,dVdtsum,divtsum
!real(REAL_KIND) :: Vn   ! to normalise volumes
type(cell_type), pointer :: cp

!Vn = Vdivide0/1.6
Vsum=0
divVsum=0
dVdtsum=0
divtsum=0
n=0
do kcell = 1,nlist
    cp => cell_list(kcell)
    if (cp%state == DEAD) cycle
    n = n+1
    Vsum = Vsum + cp%V
    divVsum = divVsum + cp%divide_volume
    dVdtsum = dVdtsum + cp%dVdt
    divtsum = divtsum + cp%divide_time
enddo
write(nflog,'(a,4e12.3)') 'ave volume, divide_volume, dVdt, divide_time: ', Vsum/n,divVsum/n,dVdtsum/n,divtsum/n
!write(*,'(a,4e12.3)') 'ave volume, divide_volume, dVdt, divide_time: ', Vsum/n,divVsum/n,dVdtsum/n,divtsum/n
end subroutine


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Execute(ncpu,infile_array,inbuflen,outfile_array,outbuflen,res) BIND(C) 
!DEC$ ATTRIBUTES DLLEXPORT :: execute
use, intrinsic :: iso_c_binding
character(c_char), intent(in) :: infile_array(*), outfile_array(*)
integer(c_int) :: ncpu, inbuflen, outbuflen, res
character*(2048) :: infile, outfile, logfile
character*(13) :: fname
character*(1) :: numstr
logical :: ok, success, isopen
integer :: i

write(*,*) 'Execute 1: synch_phase, synch_fraction: ', synch_phase,synch_fraction

write(*,*) 'ncpu, inbuflen, outbuflen: ',ncpu, inbuflen, outbuflen
res = 0
use_TCP = .false.
if (use_TCP) then
	use_PEST = .false.
endif
infile = ''
do i = 1,inbuflen
	infile(i:i) = infile_array(i)
enddo
outfile = ''
do i = 1,outbuflen
	outfile(i:i) = outfile_array(i)
enddo

inquire(unit=nflog,OPENED=isopen)
if (isopen) then
	close(nflog)
endif
i = index(infile,'.')
logfile = infile(1:i)//'log'
write(*,*) 'infile: ',trim(infile)
write(*,*) 'logfile: ',trim(logfile)
!open(nflog,file='drm_monolayer.log',status='replace')
open(nflog,file=logfile,status='replace')
write(*,*) 'Execute 2: synch_phase, synch_fraction: ', synch_phase,synch_fraction

write(nflog,*) 'irun: ',res
if (phase_log) then
    write(numstr,'(i1)') res
    fname = 'phase'//numstr//'.log'
    write(nflog,'(a)') fname
    inquire(unit=nfphase,OPENED=isopen)
    if (isopen) close(nfphase)
    open(nfphase,file=fname,status='replace')
endif

res = 0
#ifdef GFORTRAN
    write(logmsg,'(a)') 'Built with GFORTRAN'
	call logger(logmsg)
#endif

logmsg = 'OS??'
#ifdef LINUX
    write(logmsg,'(a)') 'OS is Linux'
#endif
#ifdef OSX
    write(logmsg,'(a)') 'OS is OS-X'
#endif
#ifdef _WIN32
    write(logmsg,'(a)') 'OS is Windows'
#endif
#ifdef WINDOWS
    write(logmsg,'(a)') 'OS is Windows'
#endif
call logger(logmsg)

!#ifdef OPENMP
#if defined(OPENMP) || defined(_OPENMP)
    write(logmsg,'(a)') 'Executing with OpenMP'
	call logger(logmsg)
#endif

write(logmsg,*) 'inputfile:  ', trim(infile)
call logger(logmsg)
write(logmsg,*) 'outputfile: ', trim(outfile)
call logger(logmsg)

#if 0
if (use_TCP) then
	call connecter(ok)
	if (.not.ok) then
		call logger('Failed to make TCP connections')
		return
	endif
endif
#endif

DELTA_T = 600
nsteps = 100
res=0
write(*,*) 'Execute 3: synch_phase, synch_fraction: ', synch_phase,synch_fraction

call Setup(ncpu,infile,outfile,ok)
if (ok) then
!	clear_to_send = .true.
!	simulation_start = .true.
else
	call logger('=== Setup failed ===')
endif
if (ok) then
	res = 0
else
    write(nflog,*) 'Setup error'
	res = 1
endif
execute_t1 = wtime()
!call counter

!call test_Pmis
!close(nflog)
!stop
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine DisableTCP
!DEC$ ATTRIBUTES DLLEXPORT :: disableTCP
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, ALIAS:"DISABLETCP" :: disableTCP

use_TCP = .false.   ! because this is called from monolayer_main()	
end subroutine

#if 0
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Connection(awp,port,error)
TYPE(winsockport) :: awp
integer :: port, error
integer :: address = 0
!!!character*(64) :: ip_address = "127.0.0.1"C      ! need a portable way to make a null-terminated C string
character*(64) :: host_name = "localhost"

if (.not.winsock_init(1)) then
    call logger("winsock_init failed")
    stop
endif

awp%handle = 0
awp%host_name = host_name
awp%ip_port = port
awp%protocol = IPPROTO_TCP
call Set_Winsock_Port (awp,error)

if (.not.awp%is_open) then
    write(nflog,*) 'Error: connection: awp not open: ',port
else
    write(nflog,*) 'connection: awp open: ',port, error
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Connecter(ok)
logical :: ok
integer :: error

! Main connection
ok = .true.
error = 0
call Connection(awp_0,TCP_PORT_0,error)
if (awp_0%handle < 0 .or. error /= 0) then
    write(logmsg,'(a)') 'TCP connection to TCP_PORT_0 failed'
    call logger(logmsg)
    ok = .false.
    return
endif
if (.not.awp_0%is_open) then
	write(logmsg,'(a)') 'No connection to TCP_PORT_0'
    call logger(logmsg)
    ok = .false.
    return
endif
write(logmsg,'(a)') 'Connected to TCP_PORT_0  '
call logger(logmsg)

if (use_CPORT1) then
	call connection(awp_1,TCP_PORT_1,error)
	if (awp_1%handle < 0 .or. error /= 0) then
		write(logmsg,'(a)') 'TCP connection to TCP_PORT_1 failed'
		call logger(logmsg)
		ok = .false.
		return
	endif
	if (.not.awp_1%is_open) then
		write(logmsg,'(a)') 'No connection to TCP_PORT_1'
		call logger(logmsg)
		ok = .false.
		return
	endif
	write(logmsg,'(a)') 'Connected to TCP_PORT_1  '
	call logger(logmsg)
endif
! Allow time for completion of the connection
call sleeper(2)
end subroutine
#endif 

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine terminate_run(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: terminate_run 
use, intrinsic :: iso_c_binding
integer(c_int) :: res
character*(8), parameter :: quit = '__EXIT__'
integer :: error, i

call logger('terminate_run')
call Wrapup

if (res == 0) then
	call logger(' Execution successful!')
elseif (res == -1) then
	call logger(' Execution stopped')
elseif (res == 2) then
	call logger(' No more live cells')
elseif (res == 6) then
	call logger(' Spheroid size limit reached')
elseif (res == 3) then
	call logger(' === Execution failed === ERROR in GrowCells')
elseif (res == 4) then
	call logger(' === Execution failed === ERROR in diff_solver')
elseif (res == 5) then
	call logger(' === Execution failed === ERROR in Solver')
endif
write(logmsg,'(a,f10.2)') 'Execution time (min): ',(wtime() - execute_t1)/60
call logger(logmsg)

!close(nflog)

#if 0
if (use_TCP) then
	if (stopped) then
	    call winsock_close(awp_0)
	    if (use_CPORT1) call winsock_close(awp_1)
	else
	    call winsock_send(awp_0,quit,8,error)
	    call winsock_close(awp_0)
		if (use_CPORT1) then
			call winsock_send(awp_1,quit,8,error)
			call winsock_close(awp_1)
		endif
	endif
endif
#endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Wrapup
integer :: ierr, ichemo, idrug
logical :: isopen

call logger('doing wrapup ...')
ierr = 0
!if (allocated(zoffset)) deallocate(zoffset)
!if (allocated(zdomain)) deallocate(zdomain)
if (allocated(gaplist)) deallocate(gaplist,stat=ierr)
!if (allocated(occupancy)) deallocate(occupancy)
!if (allocated(cell_list)) deallocate(cell_list)
!if (allocated(allstate)) deallocate(allstate)
!if (allocated(allstatep)) deallocate(allstatep)
!if (allocated(work_rkc)) deallocate(work_rkc)
!do ichemo = 1,MAX_CHEMO
!	if (allocated(chemo(ichemo)%coef)) deallocate(chemo(ichemo)%coef)
!	if (allocated(chemo(ichemo)%conc)) deallocate(chemo(ichemo)%conc)
!	if (allocated(chemo(ichemo)%grad)) deallocate(chemo(ichemo)%grad)
!enddo
!if (allocated(ODEdiff%ivar)) deallocate(ODEdiff%ivar)
!if (allocated(ODEdiff%varsite)) deallocate(ODEdiff%varsite)
!if (allocated(ODEdiff%icoef)) deallocate(ODEdiff%icoef)
if (allocated(protocol)) then
	do idrug = 0,2	!<------  change this to a variable 
		if (allocated(protocol(idrug)%tstart)) deallocate(protocol(idrug)%tstart)
		if (allocated(protocol(idrug)%tend)) deallocate(protocol(idrug)%tend)
		if (allocated(protocol(idrug)%conc)) deallocate(protocol(idrug)%conc)
		if (allocated(protocol(idrug)%dose)) deallocate(protocol(idrug)%dose)
		if (allocated(protocol(idrug)%started)) deallocate(protocol(idrug)%started)
		if (allocated(protocol(idrug)%ended)) deallocate(protocol(idrug)%ended)
	enddo
	deallocate(protocol)
endif
call logger('deallocated all arrays')

!write(logmsg,'(a,8i8)') 'rad_count: ',rad_count
!call logger(logmsg)

! Close all open files
inquire(unit=nfout,OPENED=isopen)
if (isopen) then
	close(nfout)
	call logger('closed nfout')
endif
inquire(nfres,OPENED=isopen)
if (isopen) close(nfres)
inquire(nfphase,OPENED=isopen)
if (isopen) close(nfphase)
call logger('closed files')

if (par_zig_init) then
	call par_zigfree
endif
call logger('freed par_zig')
end subroutine

!-----------------------------------------------------------------------------------------
! To test different ways of specifying phase durations
!-----------------------------------------------------------------------------------------
subroutine checkPD
real(8) :: Tc, b, V0
integer :: kcell, counts(5), ityp = 1
logical :: ok
type(cycle_parameters_type),pointer :: ccp
type(cell_type), pointer :: cp

ccp => cc_parameters(ityp)
Tc = divide_time_mean(ityp)/3600    ! seconds -> hours
b = log(2.0)/Tc
ccp%T_G1 = -(log(1-ccp%f_G1/2))/b
ccp%T_S = -(log(1-(ccp%f_G1+ccp%f_S)/2))/b - ccp%T_G1
ccp%T_M = log(1 + ccp%f_M)/b     ! Smith & Dendy
ccp%T_G2 = Tc - ccp%T_G1 - ccp%T_S - ccp%T_M

write(*,'(a,4f8.2)') 'T_G1, T_S, T_G2, T_M (Dendy): ',ccp%T_G1, ccp%T_S, ccp%T_G2, ccp%T_M
write(nflog,'(a,4f8.2)') 'T_G1, T_S, T_G2, T_M (Dendy): ',ccp%T_G1, ccp%T_S, ccp%T_G2, ccp%T_M
! hours to seconds
ccp%T_G1 = 3600*ccp%T_G1
ccp%T_S = 3600*ccp%T_S
ccp%T_G2 = 3600*ccp%T_G2
ccp%T_M = 3600*ccp%T_M
Mnodes = 1
call ArrayInitialisation(ok)
if (.not.ok) then
    write(*,*) 'Not OK'
    stop
endif
write(*,*) 'Initialised arrays'
! Create cells and assign initial cycle status
do kcell = 1,initial_count
    write(*,*) 'kcell: ',kcell
    cp => cell_list(kcell)
    cp%ID = kcell
    cp%celltype = 1
    cp%state = ALIVE
    cp%generation = 1
    cp%rad_state = 0
    cp%mitosis_duration = get_mitosis_duration()
!    write(*,*) 'mitosis_duration: ',cp%mitosis_duration
    call set_divide_volume(cp, V0)  ! sets %divide_volume and %divide_time
!    if (kcell <= 10) write(*,'(a,i6,4f8.1)') 'kcell,T_G1,T_S,T_G2,T_M: ',kcell,ccp%T_G1,ccp%T_S,ccp%T_G2,ccp%T_M
!    write(*,*) 'divide_time, fg: ',cp%divide_time,cp%fg
!    stop
    call SetInitialCellCycleStatus(kcell,cp)
enddo
! Determine phase distribution
Ncells = initial_count
counts = 0
do kcell = 1,Ncells
    cp => cell_list(kcell)
    counts(cp%phase) = counts(cp%phase) + 1
enddo
write(*,'(a,5i6)') 'counts: ',counts
write(*,'(a,5f8.3)') 'PD: ',counts/real(Ncells)
write(*,'(a,4f8.3)') 'control: ',ccp%f_G1,ccp%f_S,ccp%f_G2,ccp%f_M
write(nflog,'(a,5i6)') 'counts: ',counts
write(nflog,'(a,5f8.3)') 'PD: ',counts/real(Ncells)
write(nflog,'(a,4f8.3)') 'control: ',ccp%f_G1,ccp%f_S,ccp%f_G2,ccp%f_M
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine counter
integer :: kcell, iph, counts(4)

counts = 0
do kcell = 1,Ncells
    iph = cell_list(kcell)%phase
    iph = min (iph,M_phase)
    counts(iph) = counts(iph) + 1
enddo
write(*,'(a,5i6)') 'counter, istep: ',istep,counts  
write(nflog,'(a,5i6)') 'counter, istep: ',istep,counts  
end subroutine

end module
