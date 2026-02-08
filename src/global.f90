! Global definitions

module global

use omp_lib
use real_kind_mod
use par_zig_mod
!use winsock
use, intrinsic :: ISO_C_BINDING

implicit none

integer, parameter :: NP = 4

integer, parameter :: TCP_PORT_0 = 5000		! main communication port (logging)
integer, parameter :: TCP_PORT_1 = 5001		! data transfer port (plotting)
integer, parameter :: NORMAL_DIST      = 1
integer, parameter :: LOGNORMAL_DIST   = 2
integer, parameter :: EXPONENTIAL_DIST = 3
integer, parameter :: CONSTANT_DIST    = 4

!integer, parameter :: DIVIDING  = 1
!integer, parameter :: QUIESCENT = 2
!integer, parameter :: DEAD      = 3
integer, parameter :: ALIVE = 1
integer, parameter :: DYING = 2
integer, parameter :: DEAD = 3
integer, parameter :: DIVIDED = 4
integer, parameter :: EVALUATED = 5

integer, parameter :: G1_phase      = 1
integer, parameter :: G1_checkpoint = 6
integer, parameter :: S_phase       = 2
integer, parameter :: S_checkpoint  = 7
integer, parameter :: G2_phase      = 3
integer, parameter :: G2_checkpoint = 8
integer, parameter :: M_phase       = 4
integer, parameter :: dividing      = 5

integer, parameter :: OUTSIDE_TAG  = -1
integer, parameter :: UNREACHABLE_TAG  = -2

integer, parameter :: DIVIDE_ALWAYS_PUSH  = 1
integer, parameter :: DIVIDE_USE_CLEAR_SITE  = 2
integer, parameter :: DIVIDE_USE_CLEAR_SITE_RANDOM  = 3

integer, parameter :: nfin=10, nfout=11, nflog=12, nfres=13, nfrun=14, nfcell=15, nftreatment=16, nfphase=17, &
					  nfpar=18, nftcp=20, nfpest = 21
integer, parameter :: neumann(3,6) = reshape((/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /), (/3,6/))

integer, parameter :: CFSE = 0
integer, parameter :: OXYGEN = 1
integer, parameter :: GLUCOSE = 2
integer, parameter :: LACTATE = 3       ! dropped
integer, parameter :: GLUTAMINE = 4     ! dropped
integer, parameter :: OTHERNUTRIENT = 5 ! dropped
integer, parameter :: DRUG_A = 3        ! was 6
integer, parameter :: TPZ_DRUG = DRUG_A
integer, parameter :: TPZ_DRUG_METAB_1 = TPZ_DRUG + 1
integer, parameter :: TPZ_DRUG_METAB_2 = TPZ_DRUG + 2
integer, parameter :: DRUG_B = DRUG_A + 3
integer, parameter :: DNB_DRUG = DRUG_B
integer, parameter :: DNB_DRUG_METAB_1 = DNB_DRUG + 1
integer, parameter :: DNB_DRUG_METAB_2 = DNB_DRUG + 2
!integer, parameter :: MAX_CHEMO = DRUG_B + 2
integer, parameter :: NUTS = DRUG_A - 1
integer, parameter :: MAX_CHEMO = NUTS + 2

integer, parameter :: GROWTH_RATE = MAX_CHEMO + 1	! (not used here, used in the GUI)
integer, parameter :: CELL_VOLUME = MAX_CHEMO + 2
integer, parameter :: O2_BY_VOL = MAX_CHEMO + 3
integer, parameter :: CYCLE_PHASE = MAX_CHEMO + 4

integer, parameter :: N_EXTRA = CYCLE_PHASE - MAX_CHEMO + 1	! = 5 = total # of variables - MAX_CHEMO
integer, parameter :: NCONST = MAX_CHEMO

integer, parameter :: TPZ_CLASS = 1
integer, parameter :: DNB_CLASS = 2
integer, parameter :: DRUG_EVENT = 1
integer, parameter :: RADIATION_EVENT = 2
integer, parameter :: MEDIUM_EVENT = 3

integer, parameter :: NTCP = 1000

integer, parameter :: DIST_NV = 20

integer, parameter :: EXTRA = 1
integer, parameter :: INTRA = 2
integer, parameter :: MAX_CELLTYPES = 2
integer, parameter :: MAX_DRUGTYPES = 2
integer, parameter :: max_nlist = 300000
integer, parameter :: NRF = 4
integer, parameter :: LIMIT_THRESHOLD = 1500

logical, parameter :: use_ODE_diffusion = .true.
logical, parameter :: compute_concentrations = .true.
logical, parameter :: use_division = .true.
logical, parameter :: use_death = .true.
logical, parameter :: use_react = .true.
logical, parameter :: use_migration = .false.	! causing an error with vacant site becoming bdry 
logical, parameter :: use_medium_flux = .true.	! flux of constituents between spheroid and medium is accounted for.
logical, parameter :: use_metabolites = .true.
logical, parameter :: use_celltype_colour = .true.

logical, parameter :: use_Cex_Cin = .true.		! assume equilibrium to derive Cin from Cex
logical, parameter :: suppress_growth = .false.

logical, parameter :: OFF_LATTICE = .false.

! Jaiswal
logical :: use_slope_threshold      ! set true if CC_threshold < 0.1, then slope_threshold = CC_threshold
real(REAL_KIND) :: slope_threshold


real(REAL_KIND), parameter :: PI = 4.0*atan(1.0)
real(REAL_KIND), parameter :: CFSE_std = 0.05
real(REAL_KIND), parameter :: small_d = 0.1e-4          ! 0.1 um -> cm

!type occupancy_type
!	integer :: indx(2)
!	real(REAL_KIND) :: C(MAX_CHEMO)
!	type(boundary_type), pointer :: bdry
!	! for FD grid weighting
!	integer :: cnr(3,8)
!	real(REAL_KIND) :: wt(8)
!end type

type metabolism_type
	real(REAL_KIND) :: HIF1
	real(REAL_KIND) :: PDK1
	real(REAL_KIND) :: I_rate_max
	real(REAL_KIND) :: G_rate
	real(REAL_KIND) :: PP_rate
	real(REAL_KIND) :: P_rate 
	real(REAL_KIND) :: L_rate
	real(REAL_KIND) :: A_rate
	real(REAL_KIND) :: I_rate
	real(REAL_KIND) :: O_rate
	real(REAL_KIND) :: Gln_rate
	real(REAL_KIND) :: ON_rate
	real(REAL_KIND) :: GA_rate
	real(REAL_KIND) :: f_G
	real(REAL_KIND) :: f_P
	real(REAL_KIND) :: f_Gln
	real(REAL_KIND) :: C_P
	real(REAL_KIND) :: C_A
	real(REAL_KIND) :: recalcable     ! if > 0, can determine rates from w = f_G/f_Gu, otherwise need full procedure
	real(REAL_KIND) :: C_GlnEx_prev
end type

type cell_type
	integer :: ID
	integer :: celltype
	integer :: site(3)
	integer :: ivin
	logical :: active
	integer :: state
	logical :: Iphase
	logical :: irradiated
	real(REAL_KIND) :: f_S_at_IR
!    integer :: nspheres             ! =1 for Iphase, =2 for Mphase
!	real(REAL_KIND) :: radius(2)	! sphere radii (um) -> cm
!	real(REAL_KIND) :: centre(3,2)  ! sphere centre positions
!	real(REAL_KIND) :: d			! centre separation distance (um) -> cm
	integer :: generation
!	real(REAL_KIND) :: conc(MAX_CHEMO)
	real(REAL_KIND) :: Cin(MAX_CHEMO)
!	real(REAL_KIND) :: Cex(MAX_CHEMO)
	real(REAL_KIND) :: dCdt(MAX_CHEMO)
	real(REAL_KIND) :: dMdt(MAX_CHEMO)      ! mumol/s
	real(REAL_KIND) :: CFSE
	real(REAL_KIND) :: dVdt             ! actual growth rate
	real(REAL_KIND) :: V			    ! actual volume cm3
	real(REAL_KIND) :: growth_rate_factor	! to introduce some random variation 
	real(REAL_KIND) :: ATP_rate_factor	! to introduce some random variation 
	real(REAL_KIND) :: divide_volume	! actual divide volume
	real(REAL_KIND) :: divide_time      ! cycle time
	real(REAL_KIND) :: fg(4)			! to make sum(T_G1, T_S, T_G2, T_M) consistent with Tdivide
	real(REAL_KIND) :: t_divide_last	! these two values are used for colony simulation
	real(REAL_KIND) :: t_divide_next
	real(REAL_KIND) :: t_S_phase
	real(REAL_KIND) :: birthtime
    real(REAL_KIND) :: t_mitosis		! time from IR to end of mitosis
!	real(REAL_KIND) :: t_anoxia
!	real(REAL_KIND) :: t_anoxia_die
!	real(REAL_KIND) :: t_aglucosia
!	real(REAL_KIND) :: t_aglucosia_die
	real(REAL_KIND) :: M
	real(REAL_KIND) :: p_rad_death
	real(REAL_KIND) :: p_drug_death(MAX_DRUGTYPES)
	real(REAL_KIND) :: t_start_mitosis
	real(REAL_KIND) :: t_start_G2
	real(REAL_KIND) :: G2_time
	real(REAL_KIND) :: mitosis
	real(REAL_KIND) :: CP_delay
!	logical :: growth_delay
!	real(REAL_KIND) :: dt_delay
!	real(REAL_KIND) :: t_growth_delay_end	! this is for suppression of growth before first division
!	integer :: N_delayed_cycles_left		! decremented by 1 at each cell division
	real(REAL_KIND) :: tag_time             ! time cell is tagged to die metabolically
	logical :: radiation_tag    !, ATP_tag, GLN_tag
	logical :: drug_tag(MAX_DRUGTYPES)
	real(REAL_KIND) :: apoptosis_delay
	integer :: rad_state    ! 0 = pre-radiation or post-mitosis, I > 0 = radiated in phase I: G1 = 1, S = 2, G2 = 3
	logical :: G2_M
!	logical :: exists
!	integer :: cnr(3,8)
!	real(REAL_KIND) :: wt(8)

	! Cell cycle 
    integer :: phase
!    logical :: G1_flag, G1S_flag, S_flag, SG2_flag, G2_flag, G2M_flag
!    real(REAL_KIND) :: G1_time, S_time, G2_time, S_duration
!    real(REAL_KIND) :: G1_V, S_V, G2_V
!    real(REAL_KIND) :: G1S_time, SG2_time, G2M_time, M_time
!    real(REAL_KIND) :: doubling_time
!    logical :: arrested ! for S-phase arrest
!    real(REAL_KIND) :: S_start_time	    ! for PI labelling
    real(REAL_KIND) :: progress, fp, mitosis_duration
!    integer :: NL1, NL2(2)
    
!    integer :: N_PL, N_IRL, N_Ch1, N_Ch2
!    logical :: irrepairable
    
    ! exponential cycle time
    real(REAL_KIND) :: G1ch_entry_time, G1ch_time, G1ch_max_delay
    real(REAL_KIND) :: Sch_entry_time, Sch_time, Sch_max_delay
    real(REAL_KIND) :: G2ch_entry_time, G2ch_time, G2ch_max_delay
    
	type(metabolism_type) :: metab
	logical :: have_rates   
	integer :: ndt
	
	! DRM section
	integer :: phase0
	real(8) :: progress0
	real(8) :: pATM, pATR, DSB(NP,2), DSB0(NP,2),totDSB0, totMis
	real(8) :: Psurvive, psurvive_nodouble, mitosis_time
	real(8) :: Nmis(2)		!, Nlethal(2)
	
	! Jaiswal section (26/09/22)
	real(REAL_KIND) :: CC_act, ATR_act, ATM_act, dCC_act_dt, kccmd, kccrd, kcc
    real(REAL_KIND) :: G2t0		! initial time in G2

    ! Greens section (01/07/2023)
    real(8) :: gconc(3)		! oxygen, glucose, drug

end type

type cycle_parameters_type
    real(REAL_KIND) :: f_G1, f_S, f_G2, f_M
    real(REAL_KIND) :: T_G1, T_S, T_G2, T_M
    real(REAL_KIND) :: G1_mean_delay, S_mean_delay, G2_mean_delay
    real(REAL_KIND) :: G2_delay_factor
    real(REAL_KIND) :: Pk_G1, Pk_S, Pk_G2
    real(REAL_KIND) :: arrest_threshold
    ! Radiation damage/repair
    real(REAL_KIND) :: eta_PL, eta_IRL
    real(REAL_KIND) :: HRR_repair_base, HRR_repair_max
!    real(REAL_KIND) :: NHEJ_repair, NHEJ_misrepair	! in mc.f90
    real(REAL_KIND) :: DIM_misrepair
    real(REAL_KIND) :: fraction_Ch1, psurvive_Ch1, psurvive_Ch2		! prob of surviving mitosis
    real(REAL_KIND) :: psurvive_PL
    real(REAL_KIND) :: aTCP, bTCP
    ! Apoptosis
    real(REAL_KIND) :: apoptosis_rate
    real(REAL_KIND) :: apoptosis_median
    real(REAL_KIND) :: apoptosis_shape
    real(REAL_KIND) :: f_apoptosis_rate_lo  ! multiplying factor for low apoptosis rate
    real(REAL_KIND) :: t_apoptosis_hi       ! duration of high rate of apoptosis
end type

type drug_type
	character*(3)   :: classname
	integer         :: drugclass
	character*(16)  :: name
	integer         :: nmetabolites
	logical         :: use_metabolites
	real(REAL_KIND) :: diff_coef(0:2)
	real(REAL_KIND) :: medium_diff_coef(0:2)
	real(REAL_KIND) :: membrane_diff_in(0:2)
	real(REAL_KIND) :: membrane_diff_out(0:2)
	real(REAL_KIND) :: halflife(0:2)
	logical         :: kills(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: Kmet0(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: C2(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: KO2(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: n_O2(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: Vmax(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: Km(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: Klesion(MAX_CELLTYPES,0:2)
	integer         :: kill_model(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: death_prob(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: Kd(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: kill_O2(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: kill_drug(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: kill_duration(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: kill_fraction(MAX_CELLTYPES,0:2)
	logical         :: sensitises(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: SER_max(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: SER_Km(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: SER_KO2(MAX_CELLTYPES,0:2)
	
	logical :: phase_dependent
	logical :: active_phase(6)
end type

!type boundary_type
!    integer :: site(3)
!    type (boundary_type), pointer :: next
!end type

type dist_type
	integer :: class
	real(REAL_KIND) :: p1, p2, p3
end type

!type, bind(C) :: field_data
!	integer(c_int) :: site(3)
!	integer(c_int) :: state
!	real(c_double) :: volume
!	real(c_double) :: conc(0:MAX_CHEMO+N_EXTRA)	! Must agree with the definition in field.h: 0 = CFSE, MAX_CHEMO+1 = dVdt, MAX_CHEMO+2 = cellvolume
!end type

type, bind(C) :: dist_data
	logical(c_bool) :: used
	real(c_double) :: dv
	real(c_double) :: v0
	real(c_double) :: prob(DIST_NV)
end type

type treatment_type
	integer :: ichemo
	integer :: n
!	character*(16) :: name
	real(REAL_KIND), allocatable :: tstart(:)
	real(REAL_KIND), allocatable :: tend(:)
	real(REAL_KIND), allocatable :: conc(:)
	real(REAL_KIND), allocatable :: dose(:)
	logical, allocatable :: started(:)
	logical, allocatable :: ended(:)
end type

type event_type
	integer :: etype
	real(REAL_KIND) :: time
	integer :: idrug				! DRUG
	integer :: ichemo				! DRUG CHEMO INDEX
	real(REAL_KIND) :: volume		! DRUG MEDIUM
	real(REAL_KIND) :: conc			! DRUG
	real(REAL_KIND) :: O2conc		! DRUG
	real(REAL_KIND) :: O2flush		! DRUG
	real(REAL_KIND) :: dose			! RADIATION
	real(REAL_KIND) :: O2medium		! MEDIUM
	real(REAL_KIND) :: glumedium	! MEDIUM
	real(REAL_KIND) :: lacmedium	! MEDIUM
	logical :: full
	logical :: done
end type	

type LQ_type
	real(REAL_KIND) :: OER_am, OER_bm
	real(REAL_KIND) :: alpha_H, beta_H
	real(REAL_KIND) :: K_ms
	real(REAL_KIND) :: death_prob
!	real(REAL_KIND) :: growth_delay_factor
!	real(REAL_KIND) :: growth_delay_N
end type

type savedata_type
    logical :: active
    character*(128) :: filebase
    real(REAL_KIND) :: dt, tstart
    integer :: nt, it
end type

type(dist_type) :: divide_dist(MAX_CELLTYPES)
!type(occupancy_type), allocatable :: occupancy(:,:,:)
type(cell_type), allocatable, target :: cell_list(:)
type(treatment_type), allocatable :: protocol(:)
type(event_type), allocatable :: event(:)
!real(REAL_KIND), allocatable, target :: Cslice(:,:,:,:)
type(cell_type), target, allocatable :: ccell_list(:)
integer, allocatable :: perm_index(:)
type(cell_type), target :: master_cell

character*(12) :: dll_version, dll_run_version
character*(12) :: gui_version, gui_run_version
integer :: initial_count

integer :: nlist, Ncells, Ncells0, ncells_mphase, lastID, Ncelltypes
integer :: nColonyMax
integer :: Ncells_type(MAX_CELLTYPES), Ndying(MAX_CELLTYPES), Nviable(MAX_CELLTYPES), Ndead(MAX_CELLTYPES)
logical :: limit_stop
!integer :: nadd_sites, Nsites, Nreuse
integer :: Ndrugs_used
integer :: Nradiation_tag(MAX_CELLTYPES), NATP_tag(MAX_CELLTYPES), NGLN_tag(MAX_CELLTYPES)
integer :: Ndrug_tag(MAX_DRUGTYPES,MAX_CELLTYPES)
integer :: Nradiation_dead(MAX_CELLTYPES), NATP_dead(MAX_CELLTYPES), NGLN_dead(MAX_CELLTYPES)
integer :: Ndrug_dead(MAX_DRUGTYPES,MAX_CELLTYPES)
!logical :: use_radiation_growth_delay_all = .true.

integer :: ndoublings
real(REAL_KIND) :: doubling_time_sum

type(cycle_parameters_type), target :: cc_parameters(MAX_CELLTYPES)

logical :: drug_gt_cthreshold(MAX_DRUGTYPES)
real(REAL_KIND) :: Cthreshold
real(REAL_KIND) :: Clabel_threshold		! for labelling drug, e.g. EDU

type(savedata_type) :: saveFACS	!, saveprofile, saveslice

integer :: istep, ndays, nsteps, it_solve, NT_CONC, NT_GUI_OUT, show_progeny, ichemo_curr, NT_DISPLAY
integer :: Mnodes, ncpu_input
integer :: Nevents
real(REAL_KIND) :: DELTA_T, DELTA_X, fluid_fraction, Vsite_cm3, Vextra_cm3, Vcell_pL, tnow, DT_DISPLAY
!real(REAL_KIND) :: dxb, dxb3, dxf, dx3
!real(REAL_KIND) :: grid_offset(3)
real(REAL_KIND) :: Vcell_cm3, medium_volume0, total_volume, well_area, t_lastmediumchange
real(REAL_KIND) :: celltype_fraction(MAX_CELLTYPES)
integer :: selected_celltype
logical :: celltype_display(MAX_CELLTYPES)
real(REAL_KIND) :: MM_THRESHOLD, anoxia_threshold, t_anoxia_limit, anoxia_death_delay, Vdivide0, dVdivide
real(REAL_KIND) :: aglucosia_threshold, t_aglucosia_limit, aglucosia_death_delay, max_growthrate(MAX_CELLTYPES)
real(REAL_KIND) :: divide_time_median(MAX_CELLTYPES), divide_time_shape(MAX_CELLTYPES), divide_time_mean(MAX_CELLTYPES)
real(REAL_KIND) :: t_simulation, execute_t1
real(REAL_KIND) :: O2cutoff(3), hypoxia_threshold
real(REAL_KIND) :: growthcutoff(3)
real(REAL_KIND) :: spcrad_value
real(REAL_KIND) :: total_dMdt
!real(REAL_KIND) :: total_flux_prev, medium_Cbnd_prev
real(REAL_KIND) :: start_wtime
real(REAL_KIND) :: IR_time_h, CA_time_h

! Metabolism parameters
real(REAL_KIND) :: f_GL     ! ratio of lactate production to glucose consumption (not used in simple metab)
real(REAL_KIND) :: f_IN     ! fraction of total intermediates rate that is N-type (from glutamine only)
real(REAL_KIND) :: f_PPu    ! normal fraction of pyruvate that is processed by the TCA (rest -> lactate)
real(REAL_KIND) :: f_Gu     ! normal fraction of glycosis (r_G) going to make intermediates
real(REAL_KIND) :: f_Glnu   ! normal fraction of glutamine (r_Gn) going to make intermediates
real(REAL_KIND) :: f_ONu    ! normal fraction of ON (r_ONn) going to make intermediates
real(REAL_KIND) :: f_Pu     ! normal fraction of pyruvates (r_P) going to make intermediates
real(REAL_KIND) :: N_GA		! number of ATP molecules generated per glucose molecule in glycosis
real(REAL_KIND) :: N_GI		! number of intermediate molecules generated per glucose molecule in glycosis
real(REAL_KIND) :: N_GP		! number of pyruvate molecules generated per glucose molecule in glycosis
real(REAL_KIND) :: N_GlnA	! number of ATP molecules generated per glutamine molecule
real(REAL_KIND) :: N_GlnI	! number of intermediate molecules generated per glutamine molecule
real(REAL_KIND) :: N_PA		! number of ATP molecules generated per pyruvate molecule in pyruvate oxidation
real(REAL_KIND) :: N_PI		! number of intermediate molecules generated per pyruvate molecule in pyruvate oxidation
real(REAL_KIND) :: N_PO		! number of O2 molecules consumed per pyruvate molecule in pyruvate oxidation
real(REAL_KIND) :: N_GlnO	! number of O2 molecules consumed per glutamine molecule in glutamine oxidation
real(REAL_KIND) :: N_ONO	! number of O2 molecules consumed per ON molecule in ON oxidation
real(REAL_KIND) :: N_ONI	! number of intermediate molecules generated per ON molecule
real(REAL_KIND) :: N_ONA	! number of ATP molecules generated per ON molecule
real(REAL_KIND) :: f_ATPg	! threshold ATP production rate fractions for cell growth
real(REAL_KIND) :: f_ATPs	! threshold ATP production rate fractions for cell survival
real(REAL_KIND) :: f_ATPramp ! multiplying factor for ramp start for reducing r_G, r_P
real(REAL_KIND) :: r_Ag		! threshold ATP production rates for cell growth
real(REAL_KIND) :: r_As		! threshold ATP production rates for cell survival
real(REAL_KIND) :: CO_H		! threshold O2 for Ofactor
real(REAL_KIND) :: CG_H		! threshold glucose for Gfactor
real(REAL_KIND) :: I_rate_max   ! = r_Iu
real(REAL_KIND) :: rIA      ! = r_Iu/r_Au
real(REAL_KIND) :: C_Gln_lo, C_Gln_hi, f_rGln_lo, f_rGln_threshold, f_rON_base
real(REAL_KIND) :: Gln_Nshare
real(REAL_KIND) :: k_glutamine_decay

! By cell
type(metabolism_type), target :: phase_metabolic(3)

type(drug_type), allocatable, target :: drug(:)

integer, allocatable :: gaplist(:)
integer :: ngaps, ndivided
integer, parameter :: max_ngaps = 200000

logical :: bdry_changed
type(LQ_type) :: LQ(MAX_CELLTYPES)
character*(2048) :: inputfile
character*(2048) :: outputfile
character*(2048) :: logmsg
character*(1024) :: header
logical :: test_case(4)

!TYPE(winsockport) :: awp_0, awp_1
logical :: use_TCP = .false.         ! turned off in para_main()
logical :: use_CPORT1 = .false.
logical :: stopped, clear_to_send
logical :: simulation_start, par_zig_init, initialized
logical :: use_radiation, use_treatment
!logical :: use_growth_suppression = .true.	! see usage in subroutine CellGrowth
logical :: use_extracellular_O2 = .false.
logical :: use_V_dependence
logical :: use_divide_time_distribution
logical :: use_exponential_cycletime
logical :: use_constant_divide_volume
!logical :: use_volume_method
logical :: use_cell_cycle
logical :: use_constant_growthrate = .false. 
logical :: use_new_drugdata = .true.
logical :: randomise_initial_volume = .true.
logical :: synchronise
logical :: is_radiation
!logical :: use_FD = .true.
logical :: use_permute
logical :: use_gaplist = .true.
!logical :: relax
logical :: medium_change_step
logical :: fully_mixed
logical :: use_parallel
logical :: simulate_colony      ! colony growth will be simulated
logical :: colony_simulation    ! colony growth is being simulated
logical :: use_metabolism
logical :: use_glutamine_decay
logical :: noSS = .false.   ! if true, SS solution is not used for C_P
logical :: dbug = .false.
logical :: bdry_debug

logical :: use_events = .true.

real(REAL_KIND) :: ysave(100000),dCreactsave(100000)

integer :: divide_option = DIVIDE_USE_CLEAR_SITE
!integer :: divide_option = DIVIDE_ALWAYS_PUSH
integer :: idbug = 0
integer :: Nbnd
integer :: seed(2)
integer :: kcell_dbug
integer :: kcell_now
integer :: kcell_test = 1

! PEST variables
logical :: use_PEST = .false.
character*(128) :: PEST_parfile, PEST_outputfile

! C_G base concentration (for checking expt. results)
real(REAL_KIND) :: C_G_base = 10.359

! DNA-PK parameters (temporary)
!logical :: DRUG_A_inhibiter = .false.   ! set true if drug is in the input file
!real(REAL_KIND) :: C_inhibiter
!real(REAL_KIND) :: a_inhibit = 0.83, b_inhibit = 0.5

integer :: rad_count(6)

! These are now set in main
logical :: use_synchronise 
integer :: synch_phase
real(REAL_KIND) :: synch_fraction
logical :: single_cell

! DRM section
logical :: DRM = .true.
logical, parameter :: use_Napop = .true.   ! use count of apoptosed cells in SFave calculation - true for consistency with CA
integer :: N_checkpoint     ! number of cells in checkpoint - not growing
integer :: ntphase(8)
integer :: NPsurvive, Nirradiated, Napop, Nmitotic, Nsecond
real(REAL_KIND), allocatable :: Psurvive(:)
!real(REAL_KIND) :: CA_time = 99*60*60   ! seconds, default value overridden by protocol
logical :: include_daughters = .true.
!logical, parameter :: phase_dist = .true.
real(REAL_KIND) :: t_irradiation, SFave, t_mitosis
logical :: use_SF = .true.
real(REAL_KIND) :: phase_hour(60)
integer :: nphase_hours, next_phase_hour
real(REAL_KIND) :: phase_dist(0:4)    ! % of cells in each phase
real(REAL_KIND) :: recorded_phase_dist(60,0:4)   ! % of cells in each phase phase_hour after IR
real(REAL_KIND) :: recorded_DNA_rate(60)         ! average S-phase DNA rate in each phase_hour after IR
integer, allocatable :: nphase(:,:)
logical, parameter :: hourly_cycle_dist = .true.
logical, parameter :: track_cell_phase = .true.
!logical, parameter :: S_phase_RR = .false.
logical :: S_phase_RR = .false.
real(REAL_KIND) :: S_phase_RR_progress
real(REAL_KIND) :: G1_delay, S_delay, G2_delay  ! test values are specified in drm_monolayer_main
real(REAL_KIND) :: totNmis = 0
real(REAL_KIND) :: G2_katm3_factor=1.0, G2_katm4_factor=1.0, G2_katr3_factor=1.0, G2_katr4_factor=1.0
!integer :: maxhours = 199
logical :: overstepped

logical, parameter :: no_S_Iliakis = .false.	! If true this suppresses Iliakis effect in S-phase
logical, parameter :: constant_S_pHR = .true.	! not used now, instead f_S_decay
real(REAL_KIND), parameter :: dose_threshold = 0	! was 1
integer :: ATR_in_S = 1		! 0 = no ATR signalling in S, 1 = signalling, no CP effect, 2 = signalling and CP effect
logical, parameter :: use_Arnould = .true.
real(REAL_KIND) :: Reffmin = 0.7, Kclus = 0.693	! for DSB clustering
integer :: expt_ID
logical :: use_equal_mitrates = .false.
logical :: use_cell_kcc_dependence = .true.

! Greens function section
logical :: greens = .false.
integer :: NgreenCells

logical :: phase_log = .false.
logical :: use_inhibiter = .false.
logical :: use_fixed_CP = .false.
logical :: compute_cycle
logical :: SFdone
logical :: use_constant_drug_conc = .true.   ! This is used in the SN39536 experiments

!real(REAL_KIND) :: mitosis_std = 0.0    ! as a fraction of mean T_M
real(REAL_KIND) :: mitosis_std = 0.1336*3600     ! Chao 2019, Erlang k=14, L = 28, hours -> seconds
logical, parameter :: drop_mitotic_cells = .false.
logical, parameter :: allow_second_mitosis = .true.
real(8), parameter :: kmit = 0.023  ! Baide et al., dose = 1 in M --> Psurvive = 0.2

real(REAL_KIND) :: phase_exit_time_sum
integer :: npet

! Drug half-life simulation
logical :: use_drug_halflife
real(REAL_KIND) :: Khalflife, drug_time, drug_conc0
real(REAL_KIND) :: rad_dose, flush_time_h

! DNA-PK inhibition parameters
real(8) :: fDNAPK, Chalf, fDNAPKmin
logical :: suppress_ATR

real(8) :: eta_Z = 4	! to set Z in eta_Arnould()

logical :: test_run = .false.	! to check Psurvive etc
LOGICAL :: use_no_random = .false.	! to turn off variation in cycle time, DSB_Gy

!DEC$ ATTRIBUTES DLLEXPORT :: nsteps, use_PEST, PEST_outputfile, DELTA_T
!DEC$ ATTRIBUTES DLLEXPORT :: use_synchronise, synch_phase, synch_fraction	!, nfphase
!!DEC$ ATTRIBUTES DLLEXPORT :: SFave, S_phase_RR, S_phase_RR_progress, G1_delay, S_delay, G2_delay, compute_cycle, use_fixed_CP
!!DEC$ ATTRIBUTES DLLEXPORT :: G2_katm3_factor, G2_katm4_factor, G2_katr3_factor, G2_katr4_factor
contains

!-----------------------------------------------------------------------------------------
! WTIME returns a reading of the wall clock time.
!-----------------------------------------------------------------------------------------
real(DP) function wtime()
!DEC$ ATTRIBUTES DLLEXPORT :: wtime
  integer :: clock_max, clock_rate, clock_reading

  call system_clock ( clock_reading, clock_rate, clock_max )
  wtime = real(clock_reading,kind=DP)/clock_rate
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine logger(msg)
character*(*) :: msg
integer :: error
logical :: logfile_isopen
character*(1) :: LF = char(94)

error = 0
inquire(unit=nflog,OPENED=logfile_isopen)
if (use_TCP) then
#if 0
    if (awp_0%is_open) then
        call winsock_send(awp_0,trim(msg)//LF,len_trim(msg)+1,error)
    elseif (logfile_isopen) then
        write(nflog,*) trim(msg)
    else
        write(99,*) trim(msg)
    endif
#endif
else
	write(*,*) trim(msg)
endif
if (logfile_isopen) then
	write(nflog,'(a)') trim(msg)
	if (error /= 0) then
	    write(nflog,'(a,i4)') 'winsock_send error: ',error
	    close(nflog)
	endif
endif
if (error /= 0) stop
end subroutine

!---------------------------------------------------------------------
! Uniform randomly generates an integer I: n1 <= I <= n2
!---------------------------------------------------------------------
integer function random_int(n1,n2,kpar)
integer :: n1,n2,kpar
integer :: k,R

if (n1 == n2) then
    random_int = n1
elseif (n1 > n2) then
    write(logmsg,*) 'ERROR: random_int: n1 > n2: ',n1,n2
    call logger(logmsg)
    stop
endif
R = par_shr3(kpar)
k = abs(R)
random_int = n1 + mod(k,(n2-n1+1))
end function

!--------------------------------------------------------------------------------
! Make a random choice of an integer from 1 - N on the basis of probabilities in
! the array p(:) (assumed to be normalized).
!--------------------------------------------------------------------------------
integer function random_choice(p,N,kpar)
integer :: N,kpar
real(REAL_KIND) :: p(:)
integer :: k
real(REAL_KIND) :: R, psum

R = par_uni(kpar)
psum = 0
do k = 1,N
    psum = psum + p(k)
    if (R <= psum) then
        random_choice = k
        return
    endif
enddo
write(logmsg,*) 'ERROR: random_choice: ',N,p
call logger(logmsg)
stop
end function

!-----------------------------------------------------------------------------------------
! Returns a unit vector with random 3D direction
!-----------------------------------------------------------------------------------------
subroutine get_random_vector3(v)
real(REAL_KIND) :: v(3)
real(REAL_KIND) :: R1, R2, s, a
integer :: kpar=0

R1 = par_uni(kpar)
R2 = par_uni(kpar)
s = sqrt(R2*(1-R2))
a = 2*PI*R1
v(1) = 2*cos(a)*s
v(2) = 2*sin(a)*s
v(3) = 1 - 2*R2
end subroutine

!--------------------------------------------------------------------------------
! Returns a permutation of the elements of a()
!--------------------------------------------------------------------------------
subroutine permute(a,n,kpar)
integer :: a(*),n,kpar
integer :: i,k,tmp

do i = 1,n
    k = random_int(1,n,kpar)
	tmp = a(i)
	a(i) = a(k)
	a(k) = tmp
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine waste_time(n,dummy)
integer :: k, n
real(REAL_KIND) :: dummy
real(REAL_KIND) :: rsum,R
integer :: kpar=0

rsum = 0
do k = 1,n
    R = par_uni(kpar)
    rsum = rsum + R
enddo
dummy = rsum
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function norm(r)
real(REAL_KIND) :: r(3)

norm = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function norm2(r)
real(REAL_KIND) :: r(3)

norm2 = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine normalize(r)
real(REAL_KIND) :: r(3)

r = r/norm(r)
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine get_vnorm(v,vnorm)
real(REAL_KIND) :: v(3), vnorm(3)
real(REAL_KIND) :: d

d = dot_product(v,v)
vnorm = v/sqrt(d)
end subroutine

!-----------------------------------------------------------------------------------------
! Changed from kcell argument to cp, because it can be a colony cell, in ccell_list.
! Assumes maximum growth rate.
! S-phase duration is fixed.  
! M-phase duration has already been assigned as cp%mitosis_duration
! Times are seconds
!-----------------------------------------------------------------------------------------
subroutine set_divide_volume(cp,V0)
type(cell_type), pointer :: cp
real(REAL_KIND) :: V0
real(REAL_KIND) :: Tdiv, fg(4), Tfixed, Tgrowth0, Tgrowth, rVmax, V
real(REAL_KIND) :: T_G1, T_S, T_G2, T_M, scale, Tsum
integer :: ityp, kpar=0
type(cycle_parameters_type), pointer :: ccp

ityp = cp%celltype
ccp => cc_parameters(ityp)

rVmax = max_growthrate(ityp)
#if 0
if (use_exponential_cycletime) then
    if (ccp%Pk_G1 > 0) then
        cp%G1ch_time = rv_exponential(ccp%Pk_G1,kpar)
    else
        cp%G1ch_time = 1.0e10
    endif
    if (ccp%Pk_S > 0) then
        cp%Sch_time = rv_exponential(ccp%Pk_S,kpar)
    else
        cp%Sch_time = 1.0e10
    endif
    if (ccp%Pk_G2 > 0) then
        cp%G2ch_time = rv_exponential(ccp%Pk_G2,kpar)
    else
        cp%G2ch_time = 1.0e10
    endif
    Tgrowth = ccp%T_G1 + ccp%T_S + ccp%T_G2
    Tfixed = cp%G1ch_time + cp%Sch_time + cp%G2ch_time + ccp%T_M
    Tdiv = Tgrowth + Tfixed
    fg = 1
else
#endif

if (test_run .OR. use_no_random) then
    cp%fg = 1.0
    cp%divide_volume = 2*V0
    cp%divide_time = divide_time_mean(1)
else
T_S = ccp%T_S
T_M = cp%mitosis_duration
Tdiv = DivideTime(ityp)     ! log-normally distributed
Tdiv = min(Tdiv,1.2*divide_time_median(ityp))
if (allow_second_mitosis) then	! makes no difference to T_G1, T_S, T_G2, T_M
	Tsum = ccp%T_G1 + ccp%T_S + ccp%T_G2 !+ ccp%T_M
	Tfixed = T_M	! was 0
	fg(M_phase) = T_M/ccp%T_M
else
	Tsum = ccp%T_G1 + ccp%T_S + ccp%T_G2
	Tfixed = T_M
	fg(M_phase) = 1.0
endif
Tgrowth = Tdiv - Tfixed
T_G1 = Tgrowth*ccp%T_G1/Tsum
T_S = Tgrowth*ccp%T_S/Tsum
T_G2 = Tgrowth*ccp%T_G2/Tsum
scale = Tgrowth/Tsum
fg(G1_phase) = T_G1/ccp%T_G1
fg(S_phase) = T_S/ccp%T_S
fg(G2_phase) = T_G2/ccp%T_G2

if (single_cell) fg = 1.0	! mean phase times
V = V0 + rVmax*(T_G1/fg(G1_phase) + T_S/fg(S_phase) + T_G2/fg(G2_phase))
cp%divide_volume = V
cp%divide_time = Tdiv   ! cycle time, varies with cell
cp%fg = fg
if (single_cell) write(nflog,'(a,4f8.4)') 'cp%fg: ',cp%fg
endif
!if (kcell_now <= 20) then
!    write(*,*)
!    write(*,*) 'set_divide_volume: ',kcell_now
!    write(*,'(a,2e12.3)') 'Total cycle volume: ',rVmax*(T_G1/fg(G1_phase) + T_S/fg(S_phase) + T_G2/fg(G2_phase)),Vdivide0/2
!    write(*,'(a,3f8.2)') 'phase times: ',T_G1/3600,T_S/3600,T_G2/3600
!    write(*,'(a,e12.3,3f8.3)') 'rVmax,fg: ',rVmax,fg(G1_phase),fg(S_phase),fg(G2_phase)
!!    write(*,'(a,6f9.2)') 'T G1, S, G2, M: ',ccp%T_G1/3600,ccp%T_S/3600,ccp%T_G2/3600,ccp%T_M/3600,(ccp%T_G1+ccp%T_S+ccp%T_G2+ccp%T_M)/3600,divide_time_mean(1)/3600
!    write(*,'(a,6f9.2)') 'T G1, S, G2, M, Tdiv: ',T_G1/3600,T_S/3600,T_G2/3600,T_M/3600,(T_G1+T_S+T_G2+T_M)/3600,Tdiv/3600
!    write(*,'(a,5e14.5)') 'set_divide_volume: rVmax,Vdivide0,cp%divide_volume: ',rVmax,Vdivide0,cp%divide_volume
!!    write(*,'(a,2e14.5)') 'growth at average growth rate: ',Vdivide0/2,(divide_time_mean(1) - ccp%T_M)*rVmax
!endif
if (single_cell) then
	write(*,'(a,4f8.3)') 'single_cell fg: ',fg
endif
end subroutine	

!--------------------------------------------------------------------------------------
! This is actually cycle time
!--------------------------------------------------------------------------------------
real(REAL_KIND) function DivideTime(ityp)
integer :: ityp
real(REAL_KIND) :: p1, p2
integer :: kpar = 0

dividetime = 0
p1 = divide_dist(ityp)%p1
p2 = divide_dist(ityp)%p2
select case (divide_dist(ityp)%class)
case (NORMAL_DIST)
	DivideTime = rv_normal(p1,p2,kpar)
case (LOGNORMAL_DIST)
	DivideTime = rv_lognormal(p1,p2,kpar)
case (CONSTANT_DIST)
	DivideTime = p1
end select
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real function DivisionTime(ityp)
integer :: ityp
integer :: kpar = 0
real(REAL_KIND), parameter :: rndfraction = 0.2

DivisionTime = rv_lognormal(divide_dist(ityp)%p1,divide_dist(ityp)%p2,kpar)
end function

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
function get_I2Divide(cp) result(I2div)
type(cell_type), pointer :: cp
real(REAL_KIND) :: I2div
integer :: ityp
type(cycle_parameters_type), pointer :: ccp
type(metabolism_type), pointer :: metabolic
	
metabolic => phase_metabolic(1)

ityp = cp%celltype
ccp => cc_parameters(ityp)
!I2div = cp%divide_time*metabolic(ityp)%I_rate_max
I2div = (ccp%T_G1 + ccp%T_S + ccp%T_G2)*I_rate_max
!		*metabolic%I_rate_max
end function

!--------------------------------------------------------------------------------------
! par_rnor is N(0,1)
! p1 = mean
! p2 = std deviation
!--------------------------------------------------------------------------------------
real(REAL_KIND) function rv_normal(p1,p2,kpar)
integer :: kpar
real(REAL_KIND) :: p1,p2
real(REAL_KIND) :: R

R = par_rnor(kpar)
rv_normal = p1+R*p2
end function

!--------------------------------------------------------------------------------------
! When Y is normal N(p1,p2) then X = exp(Y) is lognormal with
!   median = m = exp(p1)
!   shape  = s = exp(p2)
! Also median = m = mean/(s^2/2)
! kpar = parallel process number
!--------------------------------------------------------------------------------------
real(REAL_KIND) function rv_lognormal(p1,p2,kpar)
integer :: kpar
real(REAL_KIND) :: p1,p2
real(REAL_KIND) :: R,zz

R = par_rnor(kpar)
zz = p1 + R*p2
rv_lognormal = exp(zz)
end function

!--------------------------------------------------------------------------------------
! For testing.
!--------------------------------------------------------------------------------------
real(REAL_KIND) function my_rnor()
real(REAL_KIND) :: sum, R
integer :: k
integer :: kpar=0

sum = 0
do k = 1,12
    R = par_uni(kpar)
    sum = sum + R
enddo
my_rnor = sum - 6.0
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real(REAL_KIND) function rv_exponential(lambda,kpar)
real(REAL_KIND) :: lambda
integer :: kpar
real(REAL_KIND) :: R

R = par_uni(kpar)
rv_exponential = -log(1-R)/lambda
end function

!--------------------------------------------------------------------------------------
! Cumulative probability distribution for a lognormal variate with median m, shape s
! Computes Pr(X < a) where X = exp(Y) and Y is N(p1,p2), p1 = log(m), p2 = log(s)
! Pr(X < a) = Pr(Y < log(a)) = Pr(p1 + p2*R < log(a)) = Pr(R < (log(a)-p1)/p2)
! where R is N(0,1)
!--------------------------------------------------------------------------------------
real(REAL_KIND) function cum_prob_lognormal(a,p1,p2)
real(REAL_KIND) :: a, p1, p2
real(REAL_KIND) :: b, prob

b = (log(a) - p1)/p2
prob = 0.5 + 0.5*erf(b/sqrt(2.0))
cum_prob_lognormal = prob
end function

!-----------------------------------------------------------------------------------------
! Generate a random value for CFSE from a distribution with mean = average
! In the simplest case we can allow a uniform distribution about the average.
! Multiplying factor in the range (1-a, 1+a)
! Better to make it a Gaussian distribution: 
!  = average*(1+s*R)
! where R = N(0,1), s = std deviation
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function generate_CFSE(average)
real(REAL_KIND) :: average, std
integer :: kpar = 0
real(REAL_KIND) :: R

! Uniform distribution
!R = par_uni(kpar)
!generate_CFSE = (1 - a + 2*a*R)*average
! Gaussian distribution
R = par_rnor(kpar)	! N(0,1)
generate_CFSE = (1 + CFSE_std*R)*average
end function

!-----------------------------------------------------------------------------------------
! Seconds
!-----------------------------------------------------------------------------------------
function get_mitosis_duration() result(t)
real(REAL_KIND) :: t
type(cycle_parameters_type), pointer :: ccp
integer :: ityp = 1
ccp => cc_parameters(ityp)
if (single_cell .or. test_run .OR. use_no_random) then
    t = ccp%T_M
else
    t = rv_normal(ccp%T_M, mitosis_std, 0)
endif
!t = 0.43*3600
t = max(t,0.0)
!write(*,'(a,f8.3)') 'mitosis_duration: ',t/3600
end function

!-----------------------------------------------------------------------------------------
! ityp = cell type
! V0 = cell starting volume (after division) = %volume
! Two approaches:
! 1. Use Vdivide0 and dVdivide to generate a volume
! 2. Use the divide time log-normal distribution 
!    (use_V_dependence = false)
! NOT USED
!-----------------------------------------------------------------------------------------
function get_divide_volume(ityp,V0,Tdiv,fg) result(Vdiv)
integer :: ityp
real(REAL_KIND) :: V0, Tdiv, fg
real(REAL_KIND) :: Vdiv, Tfixed, Tgrowth0, Tgrowth, rVmax
real(REAL_KIND) :: b, R
integer :: kpar=0
type(cycle_parameters_type), pointer :: ccp

ccp => cc_parameters(ityp)

rVmax = max_growthrate(ityp)
Tgrowth0 = ccp%T_G1 + ccp%T_S + ccp%T_G2
Tfixed = ccp%T_M + ccp%G1_mean_delay + ccp%G2_mean_delay
if (use_divide_time_distribution) then
	Tdiv = DivideTime(ityp)
	Tgrowth = Tdiv - Tfixed
	fg = Tgrowth/Tgrowth0
	Vdiv = V0 + Tgrowth*rVmax
else
	if (use_constant_divide_volume) then
		Vdiv = Vdivide0
		Tdiv = Tgrowth0 + Tfixed
		fg = 1
	else
		R = par_uni(kpar)
		Vdiv = Vdivide0 + dVdivide*(2*R-1)
		Tgrowth = (Vdiv - V0)/rVmax
		fg = Tgrowth/Tgrowth0
		Tdiv = Tgrowth + Tfixed
	endif
endif
!write(*,'(a,4e11.3)') 'get_divide_volume: rVmax,Vdiv,Tdiv, Tmean: ',rVmax,Vdiv,Tdiv/3600,divide_time_mean(ityp)/3600

end function	

!-----------------------------------------------------------------------------------------
! ityp = cell type
! V0 = cell starting volume (after division) = %volume
! Two approaches:
! 1. Use Vdivide0 and dVdivide to generate a volume
! 2. Use the divide time log-normal distribution 
!    (a) use_V_dependence = true
!    (b) use_V_dependence = false
! NOT USED
!-----------------------------------------------------------------------------------------
function get_divide_volume1(ityp,V0,Tdiv) result(Vdiv)
integer :: ityp
real(REAL_KIND) :: V0, Tdiv
real(REAL_KIND) :: Vdiv
real(REAL_KIND) :: Tmean, b, R
integer :: kpar=0

Tmean = divide_time_mean(ityp)
if (use_divide_time_distribution) then
	Tdiv = DivideTime(ityp)
	if (use_V_dependence) then
		b = log(2.0)*(Tdiv/Tmean)
		Vdiv = V0*exp(b)
	else
		Vdiv = V0 + (Vdivide0/2)*(Tdiv/Tmean)
	endif
else
	if (use_constant_divide_volume) then
		Vdiv = Vdivide0
	else
		R = par_uni(kpar)
		Vdiv = Vdivide0 + dVdivide*(2*R-1)
	endif
	Tdiv = Tmean
endif
end function	

!--------------------------------------------------------------------------------------
! For computing total exchange of nutrients with the medium, it is necessary to sum
! the progress factor %fp over all live cells.
!--------------------------------------------------------------------------------------
function getEffectiveNcells() result(neff)
real(REAL_KIND) :: neff
type (cell_type), pointer :: cp
integer :: kcell

neff = 0
do kcell = 1,nlist
    cp => cell_list(kcell)
    if (cp%state == DYING .or. cp%state == DEAD) cycle
    neff = neff + cp%fp
enddo
end function	

!--------------------------------------------------------------------------------------
! Determine real roots r(:) of the cubic equation:
! x^3 + a.x^2 + b.x + c = 0
! If there is one real root, n=1 and the root is r(1)
! If there are three distinct real roots, n=3 and the roots are r(1), r(2), r(3)
! If there is a repeated root, n=2 and the single root is r(1), the repeated root is r(2)
!--------------------------------------------------------------------------------------
subroutine cubic_roots(a, b, c, r, n)
real(REAL_KIND) :: a, b, c, r(3)
integer :: n
real(REAL_KIND) :: QQ, RR, theta, R2, Q3, AA, BB

QQ = (a*a - 3*b)/9
RR = (2*a*a*a - 9*a*b + 27*c)/54
Q3 = QQ*QQ*QQ
R2 = RR*RR
if (R2 < Q3) then
	n = 3
	theta = acos(RR/sqrt(Q3))
	r(1) = -2*sqrt(QQ)*cos(theta/3) - a/3
	r(2) = -2*sqrt(QQ)*cos((theta+2*PI)/3) - a/3
	r(3) = -2*sqrt(QQ)*cos((theta-2*PI)/3) - a/3
else
	n = 1
	AA = -sign(1.d0,RR)*(abs(RR) + sqrt(R2 - Q3))**(1.d0/3.d0)
	if (AA == 0) then
		BB = 0
	else
		BB = QQ/AA
	endif
	r(1) = AA + BB - a/3
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Squeeze gaps out of cellist array, adjusting occupancy array.
!-----------------------------------------------------------------------------------------
subroutine squeezer()
integer :: last, kcell, site(3), indx(2), i, j, idc, n, region

!write(*,*) 'squeezer'
!call logger('squeezer')
if (ngaps == 0) return
last = nlist
kcell = 0
n = 0
do
    kcell = kcell+1
    if (cell_list(kcell)%state == DEAD) then    ! a gap
        if (kcell == last) exit
        do
            if (last == 0) then
                write(nflog,*) 'last = 0: kcell: ',kcell
                stop
            endif
            if (cell_list(last)%state == DEAD) then
                last = last-1
                n = n+1
                if (n == ngaps) exit
            else
                exit
            endif
        enddo
        if (n == ngaps) exit
        cell_list(kcell) = cell_list(last)
        last = last-1
        n = n+1
    endif
    if (n == ngaps) exit
enddo
nlist = nlist - ngaps
ngaps = 0
if (dbug) write(nflog,*) 'squeezed: ',n,nlist

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine make_perm_index(ok)
logical :: ok
integer :: npp, kcell, kpar=0

npp = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	npp = npp + 1
	perm_index(npp) = kcell
enddo
if (npp /= Ncells) then
	write(logmsg,*) 'Error: make_perm_index: npp /= Ncells: ',npp,ncells,nlist
	call logger(logmsg)
	ok = .false.
	return
endif
if (use_permute) then
	call permute(perm_index,npp,kpar)
endif
ok = .true.
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_phase_distribution(phase_count)
integer :: phase_count(0:4)
integer :: kcell, ph
type(cell_type), pointer :: cp
integer :: G1_count(2)		! 1 = pre-mitosis, 2 = post-mitosis
phase_count = 0
G1_count = 0
do kcell = 1,nlist
    cp => cell_list(kcell)
    if (cp%state == DEAD) then
        ph = 0      ! 
    elseif (cp%phase == G1_phase .or. cp%phase == G1_checkpoint) then
        ph = 1
		if (cp%rad_state == 0) then
			G1_count(2) = G1_count(2) + 1
		else
			G1_count(1) = G1_count(1) + 1
		endif
    elseif (cp%phase == S_phase .or. cp%phase == S_checkpoint) then
        ph = 2
    elseif (cp%phase == G2_phase .or. cp%phase == G2_checkpoint) then
        ph = 3
    else
        ph = 4
!		if (next_phase_hour == 3) then
!			write(nflog,'(a,4i6,2f6.2)') 'phase_count, kcell, phase, phase0, G2t0, tmit0: ',phase_count(4),kcell,cp%phase, cp%phase0, cp%G2t0, cp%t_start_mitosis/3600
!		endif
    endif
    phase_count(ph) = phase_count(ph) + 1
enddo
!write(nflog,'(a,3i8,2f8.3)') 'G1: pre-M post-M counts, total, fractions: ',G1_count, phase_count(1),real(G1_count)/max(1,sum(G1_count))
!write(*,'(a,5i6)') 'phase_count: ',phase_count(0:4)
end subroutine

!-----------------------------------------------------------------------------------------
! https://hpaulkeeler.com/simulating-poisson-random-variables-direct-method/
!-----------------------------------------------------------------------------------------
function poisson_gen(L) result(res)
real(8) :: L
integer :: k, res, kpar=0
real(8) :: p, u
k = 0
p = 1
do
	k = k+1
	u = par_uni(kpar)
	p = p*u
	if (p < L) exit
enddo
res = k-1
end function

!------------------------------------------------------------------------
! See docs\models\DNA-PK\DNA-PK.xlsx
! This is DNAPKact(C) = H(C)
! Returns a fraction between 0 and 1.
!------------------------------------------------------------------------
function logistic(C) result(y)
real(REAL_KIND) :: C, y
real(REAL_KIND) :: EC50, hill_n
!real(REAL_KIND) :: bottom, top, EC50, hillslope, percent

hill_n = 1.0
EC50 = Chalf
if (C > 0) then
    y = 1/(1 + (C/EC50)**hill_n)
else
   y = 1
endif
y = (1 - fDNAPKmin)*y + fDNAPKmin

!bottom = 0
!top = 100A
!hillslope = -1.0	!-0.6919
!EC50 = Chalf
!if (C > 0) then
!	y = ((log10(EC50) - log10(C))*hillslope)
!    percent = bottom + (top-bottom)/(1 + 10**y)
!else
!    percent = 100
!endif
!y = percent/100
!y = (1 - fDNAPKmin)*y + fDNAPKmin
!!write(nflog,'(a,4f8.3)') 'logistic: C, Chalf, fDNAPKmin, y: ',C, Chalf, fDNAPKmin, y
end function

subroutine check_logistic
integer :: i
real(8) :: C, y

do i = 1,10
	C = i*0.1
	y = logistic(C)
	write(nflog,'(2f8.4)') C, y
enddo
do i = 1,10
	C = i*2
	y = logistic(C)
	write(nflog,'(2f8.4)') C, y
enddo
do i = 1,10
	C = 20+i*3
	y = logistic(C)
	write(nflog,'(2f8.4)') C, y
enddo
end subroutine


end module
