module mc
use real_kind_mod
use global

use chemokine   ! is this OK?
use eta_module
use km10_mod

!use greens_mod

implicit none

! There are 5 pathways:
! 1 = NHEJfast fast pre-replication (simple)
! 2 = NHEJslow slow pre-replication (complex)
! 3 = HRc slow post-replication (complex)
! 4 = HRs slow post-replication (simple)
! 5 = TMEJ very slow post-replication (complex)

!! 3 = NHEJ fast post-replication (simple) --> NHEJfast
!! 4 = HR slow post-replication (complex)
!! 5 = HR slow post-replication (simple)
!! 6 = MMEJ very slow post-replication (complex)

! There are 5 pathways:
integer, parameter :: NHEJfast = 1
integer, parameter :: NHEJslow = 2
integer, parameter :: HR = 3
integer, parameter :: TMEJ = 4   ! this is really alt-EJ (was MDR)
integer, parameter :: SSA = 5
!integer, parameter :: HRc = 3
!integer, parameter :: HRs = 4
!integer, parameter :: TMEJ = 5


! Note: McMahon parameter values are from MEDRAS code, specifically medrascell.py
!       ------------------------

!real(8) :: eta = 0.0006506387875449218
real(8) :: Pcomplex = 0.4337    ! fraction of created DSBs that are complex (McMahon: complexFrac)
real(8) :: pJeggo = 0.9         ! fraction of post-replication DSBs that are fast. (NOT USED)
real(8) :: fsmin = 0.5
!real(8) :: PHR = 0.8      ! fraction of post-replication simple DSBs that are repaired by HR (slow) rather than NHEJ (fast)
real(8) :: pHRs_max, pHRc_max, pHRs_G2, pHRc_G2
logical :: use_sigmoid = .true.
real(8) :: rmin = 0.1, kdecay = 0.1, ksig = 1, csig = 8.56    ! decay function parameters
character*(2) :: phaseName(4) = ['G1','S ','G2','M ']
!real(8) :: repRate(NP) = [2.081, 0.2604, 2.081, 0.2604, 0.008462]   ! by pathway
!real(8) :: repRate(NP) = [2.081, 0.2604, 2.081, 0.2604, 0.2604, 0.008462]   ! by pathway  with HR simple
!real(8) :: repRate(NP) = [2.081, 0.2604, 0.2604, 0.2604, 0.008462]   ! by pathway  with HR simple
real(8) :: repRate(NP) = [2.081, 0.2604, 0.087, 0.008462]   ! by pathway (TMEJ was 0.008462) (McMahon: fastRepair, slowRepair, verySlowRepair)
!real(8) :: misrepRate(8) = [0.18875,0.18875,0.18247,0.18247,0.14264,0.14264,0.18875,0.18875]        ! by phase, Nlethal = 0.5*misrepRate*Nmis 
!real(8) :: misrepRate(8) = [0.18875,0.0,0.18247,0.0,0.14264,0.0,0.18875,0.0]        ! by phase, Nlethal = 0.5*misrepRate*Nmis 
!real(8) :: fidRate(NP)  = [0.98537, 0.98537, 0.98537, 1.0, 0.4393]  ! by pathway
!real(8) :: fidRate(NP)  = [0.98537, 0.98537, 0.98537, 1.0, 1.0, 0.4393]  ! by pathway  with HR simple
!real(8) :: fidRate(NP)  = [0.98537, 0.98537, 1.0, 1.0, 0.4393]  ! by pathway  with HR simple
real(8) :: fidRate(NP)  = [0.98537, 1.0, 0.4393, 0.0]  ! by pathway  with HR simple (McMahon: NHEJFidelity, NHEJFidelity, 1.0, MMEJFidelity)
!logical :: pathwayUsed(8,NP)
real(8) :: apopRate = 0.01117   ! (McMahon: apoptoticRate) (NOT USED only baseRate is used)
real(8) :: baseRate = 0.000739  ! (McMahon: baseRate)
real(8) :: mitRate(2)   !  = 0.0141    ! (McMahon: mitoticRate)
real(8) :: Msurvival = 0.05
real(8) :: Kaber = 1.0          ! now fixed, McMahon has 1.  
real(8) :: Klethal = 0.4
real(8) :: K_ATM(3,4) ! = [0.076, 0.3, 1.0, 1.0]    ! (1) and (2) are the parameters of kinase kinetics, (3) and (4) are CP slowdown parameters
real(8) :: K_ATR(3,4) ! = [0.005, 0.3, 1.0, 1.0]
real(8) :: Ztime = 0    ! hours

real(8) :: sigma_NHEJ = 0.04187
real(8) :: sigma_TMEJ = 0.08

! DNA-PK inhibition parameters
real(8) :: Chalf    ! inhibitor concentration that halves repair rate 
real(8) :: Preass   ! rate of reassignment to pathway 4, TMEJ (prob of reass/hour)
!real(8) :: KmaxInhibitRate = 0.8
!real(8) :: b_exp = 1.0
!real(8) :: b_hill = 0.5
!logical :: use_exp_inhibit = .true.

! SSA parameters
!logical :: use_SSA = .false.
!real(8) :: SSAfrac = 0.1
!real(8) :: SSA_rep_fraction = 1.0
!real(8) :: SSA_fid_fraction = 0.5
!logical :: read_SSA_parameters = .false.
!logical :: alt_EJ_suppressed = .false.

! Jaiswal formulation (26/09/22)
real(8) :: Kcc2a, Kcc2e, Kd2e, Kd2t, Ke2cc, Kt2cc, Kti2t
real(8) :: Km1, Km10, Km10t
real(8) :: kmmp, kmmd               ! ATM production, decay
real(8) :: kmrp, kmrd               ! ATR production, decay
real(8) :: kmccp, kmccrd, kmccmd    ! CC production, ATR-driven decay, ATM-driven decay
real(8) :: CC_tot, ATR_tot, ATM_tot, CC_act0, CC_threshold, norm_factor
real(8) :: km10_alfa, km10_beta     ! for G2 only
real(8) :: G1_tdelay = 4            ! delay before ATM_act updated (hours)
logical :: use_Jaiswal = .true.
logical :: vary_km10 = .true.
real(8) :: jaiswal_std = 0.6
logical :: use_ATR_S = .false.

real(8) :: ATMsum, ATRsum, Sthsum, G2thsum(2)
integer :: NSth, NG2th
real(8) :: repRateFactor(NP)

real(8) :: totPmit, totPaber, tottotDSB, totNlethal
real(8) :: totNDSB(2), totNmisjoins, totSFfactor(3)

real(8) :: tCPdelay, tATMdelay, tATRdelay
logical :: use_addATMATRtimes = .false.
real(8) :: totG1delay, totSdelay, totG2delay
integer :: nG1delay, nSdelay, nG2delay
logical :: use_G2_pATM_Nindependent = .false.
logical :: output_DNA_rate = .false.
logical :: FIX_katm1s_eq_katm1g2 = .false.  ! handle this in the input files, making katm1g2 = katm1s
logical :: negligent_G2_CP = .false.
logical :: use_DSB_CP = .false.
logical :: use_D_model = .false.

logical :: use_km10_kcc2a_dependence = .true.
real(8) :: kcc2a_ave
logical :: use_exp_slowdown = .false.
logical :: use_G1_stop = .false.    ! These flags control use of either CP delay (true) or slowdown (false)
logical :: use_S_stop = .false.
logical :: use_G1_pATM = .false.
logical :: use_S_pATM = .false.

logical :: use_G2_stop = .false.                        ! because use_Jaiswal is true
logical :: use_phase_dependent_CP_parameters = .true.   ! now always true

! Normalisation
real(8) :: control_ave(4)   ! now set equal to ccp%f_G1, ...
logical :: normalise, M_only
character(6) :: expt_tag
logical :: G2M_only = .false.

! G1 checkpoint
logical :: use_G1_CP_factor = .false.
real(8) :: G1_CP_factor = 0.0
real(8) :: G1_CP_time

! Checking TMEJ misjoining
real(8) :: misjoins(2)      ! 1 = NHEJ, 2 = TMEJ

! Checking prob of slow repair in G2
logical :: check_G2_slow = .true.
integer :: nslow_sum
real(8) :: pHR_sum, pNHEJslow_sum, fdecay_sum

! Iliakis
integer :: nIliakis
real(8) :: kIliakis     ! if kIliakis = 0, fIliakis = 1
real(8) :: fIliakis
logical :: use_Iliakis 

! Distributions of Nlethal, totDSB at mitosis
integer, parameter :: NMDIST = 100
integer :: count_Nlethal(NMDIST), count_totDSB(NMDIST)
real(8), parameter :: ddist_Nlethal = 20.0/NMDIST
real(8), parameter :: ddist_totDSB = 200.0/NMDIST


!DEC$ ATTRIBUTES DLLEXPORT :: Pcomplex, apopRate, baseRate, mitRate, Msurvival, Kaber, Klethal, K_ATM, K_ATR !, KmaxInhibitRate, b_exp, b_hill

contains

!--------------------------------------------------------------------------
! Note: decision about PEST output is based on iphase_hours. 
!--------------------------------------------------------------------------
subroutine ReadMcParameters(nfin)
integer :: nfin
integer :: iuse_baserate, iuse_exp, iphase_hours, icase, nCPparams, iph, j
real(8) :: TMEJrep, TMEJfid, SSArep, SSAfid
real(8) :: pHR_S, pfc_S, pHR_G2, pfc_G2, k3, k4

write(*,*) 'ReadMcParameters:'
read(nfin,*) iphase_hours
write(*,*) 'iphase_hours: ',iphase_hours
read(nfin,*) baseRate
write(*,*) 'baseRate: ',baseRate
read(nfin,*) mitRate(1)
read(nfin,*) mitRate(2)
read(nfin,*) Msurvival
read(nfin,*) Klethal
!read(nfin,*) nCPparams
!if (nCPparams == 1) then
!    iph = 1
!    use_phase_dependent_CP_parameters = .false.
!    read(nfin,*) K_ATM(iph,1)
!    read(nfin,*) K_ATM(iph,2)
!    read(nfin,*) K_ATM(iph,3)
!    read(nfin,*) K_ATM(iph,4)
!    read(nfin,*) K_ATR(iph,1)
!    read(nfin,*) K_ATR(iph,2)
!    read(nfin,*) K_ATR(iph,3)
!    read(nfin,*) K_ATR(iph,4)
!    do iph = 2,3        ! Note: doing this actually makes it unnecessary to set iph = 1 when .not. use_phase_dependent_CP_parameters
!        do j = 1,4
!            K_ATM(iph,j) = K_ATM(1,j)
!            K_ATR(iph,j) = K_ATR(1,j)
!        enddo
!    enddo
!elseif (nCPparams == 3) then
    write(nflog,*) 'K_ATM'
    do j = 1,4
        do iph = 1,3
            read(nfin,*) K_ATM(iph,j)
            write(*,*) j,iph,K_ATM(iph,j)
            write(nflog,*) j,iph,K_ATM(iph,j)
        enddo
     enddo
    write(nflog,*) 'K_ATR'
    do j = 1,4
        do iph = 2,3
            read(nfin,*) K_ATR(iph,j)
            write(*,*) j,iph,K_ATR(iph,j)
            write(nflog,*) j,iph,K_ATR(iph,j)
        enddo
    enddo
    
k3 = K_ATM(S_phase,3)
k4 = K_ATM(S_phase,4)
write(nflog,'(a,2f8.3)') 'k3, k4: ',k3,k4

read(nfin,*) repRate(NHEJfast)
read(nfin,*) repRate(NHEJslow)
read(nfin,*) repRate(HR)
read(nfin,*) repRate(TMEJ)
read(nfin,*) pComplex
read(nfin,*) fsmin  !pJeggo
read(nfin,*) pHRs_max
read(nfin,*) pHRc_max
read(nfin,*) rmin
read(nfin,*) ksig
read(nfin,*) csig
read(nfin,*) nIliakis
read(nfin,*) kIliakis
read(nfin,*) G1_tdelay
read(nfin,*) Chalf
read(nfin,*) Preass
read(nfin,*) dsigma_dt
read(nfin,*) sigma_NHEJ
write(nflog,*) 'sigma_NHEJ: ',sigma_NHEJ

if (use_Jaiswal) then
    read(nfin,*) Kcc2a
    read(nfin,*) Kcc2e
    read(nfin,*) Kd2e
    read(nfin,*) Kd2t
    read(nfin,*) Ke2cc
    read(nfin,*) Kt2cc
    read(nfin,*) Kti2t
    read(nfin,*) Kmccp
    read(nfin,*) Kmccmd
    read(nfin,*) Kmccrd
    read(nfin,*) Kmrp
    read(nfin,*) Kmrd
    read(nfin,*) Kmmp
    read(nfin,*) Kmmd
    read(nfin,*) CC_tot
    read(nfin,*) ATR_tot
    read(nfin,*) ATM_tot
    read(nfin,*) CC_act0
    read(nfin,*) CC_threshold
    read(nfin,*) norm_factor
    write(nflog,*) 'norm_factor: ',norm_factor
    if (CC_threshold < 0.1) then
        use_slope_threshold = .true.
        slope_threshold = CC_threshold
    else
        use_slope_threshold = .false.
        CC_threshold = CC_threshold*CC_tot
    endif
endif
call make_eta_table(sigma_NHEJ, sigma_TMEJ, fsmin)

ATMsum = 0  ! to investigate ATM dependence on parameters
ATRsum = 0  ! to investigate ATR dependence on parameters
!Sthsum = 0
NSth = 0
G2thsum = 0
NG2th = 0

use_Iliakis = (kIliakis > 0)

! For PEST runs, iphase_hours must be -1 (for M runs) or -2 (for C runs) or -3 (for MC runs)
use_SF = .false.
nphase_hours = 0
next_phase_hour = 0
phase_hour(:) = 0
output_DNA_rate = .false.
normalise = .false.
M_only = .false.
if (iphase_hours == -1) then    ! PDSN dose = 0
    expt_tag = "PDSN0G"
    compute_cycle = .true.
    normalise = .true.
    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
    nphase_hours = 2
    next_phase_hour = 1
    phase_hour(1:2) = [5.0, 11.5]   ! these are hours post irradiation, incremented when irradiation time is known (in ReadProtocol)
    ! Note: output G1, S, G2, M
elseif (iphase_hours == -2) then    ! PDSN dose = 2
    expt_tag = "PDSN2G"
    compute_cycle = .true.
    normalise = .true.
    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
    nphase_hours = 6
    next_phase_hour = 1
    phase_hour(1:6) = [1.0, 2.0, 3.0, 5.0, 8.5, 11.5]   ! these are hours post irradiation, incremented when irradiation time is known (in ReadProtocol)
elseif (iphase_hours == -3) then    ! PDSN dose = 6
    expt_tag = "PDSN6G"
    compute_cycle = .true.
    normalise = .true.
    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
    nphase_hours = 3
    next_phase_hour = 1
    phase_hour(1:3) = [5.0, 8.5, 11.5]   ! these are hours post irradiation, incremented when irradiation time is known (in ReadProtocol)
elseif (iphase_hours == 1) then
    use_SF = .true.     ! in this case SFave only is recorded
    compute_cycle = .false.

elseif (mod(iphase_hours,10) == 2) then    ! this is the compute_cycle case for CA-135
    expt_tag = "CA-135"
    compute_cycle = .true.
    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
    nphase_hours = 5
    next_phase_hour = 1
    phase_hour(1:5) = [5.0, 8.5, 11.5, 18.5, 24.5]   ! these are hours post irradiation, incremented when irradiation time is known (in ReadProtocol)
    ! Note: output G1, S, G2, M
elseif (mod(iphase_hours,10) == 9) then    ! this is the compute_cycle case for CC-13
    expt_tag = "CC-13"
    compute_cycle = .true.
    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
    nphase_hours = 8
    next_phase_hour = 1
    phase_hour(1:nphase_hours) = [1, 2, 3, 5, 8, 12, 18, 24]   
    ! Note: output M
elseif (mod(iphase_hours,10) == 5) then    ! this is the compute_cycle case for CC-11
    expt_tag = "CC-11"
    compute_cycle = .true.
    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
    nphase_hours = 5
    next_phase_hour = 1
    phase_hour(1:5) = [1.0, 1.5, 2.5, 3.5, 4.5]   ! these are hours post irradiation, incremented when irradiation time is known (in ReadProtocol)
    ! Note: output G1, S, G2

elseif (iphase_hours == 6) then    ! this is the output_DNA_rate case (EDU)
    compute_cycle = .false.
    output_DNA_rate = .true.
    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
! Generating DNA synthesis data
!    nphase_hours = 10  ! To generate DNA synthesis factor from S ATM parameters
!    phase_hour(1:10) = [0.05,0.1,0.16666,0.33333,0.5,0.75,1.25,2.25,3.25,4.25]   ! these are hours post irradiation, incremented when irradiation time is known (in ReadProtocol)
! EDU fitting
!    nphase_hours = 5    ! To fit S ATM parameters to EDU data
!    phase_hour(1:5) = [0.75,1.25,2.25,3.25,4.25]   ! these are hours post irradiation, incremented when irradiation time is known (in ReadProtocol)
    nphase_hours = 18
    do j = 1,nphase_hours
        phase_hour(j) = j*0.25
    enddo
    next_phase_hour = 1
elseif (iphase_hours == 3) then
    use_SF = .true.     ! in this case SFave is recorded and there are multiple phase distribution recording times
    nphase_hours = 4
    next_phase_hour = 1
    phase_hour(1:4) = [8, 12, 18, 24]
elseif (iphase_hours == 4) then    ! this is the synchronised IR case
    compute_cycle = .true.
    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
    nphase_hours = 1
    next_phase_hour = 1
    phase_hour(1:5) = [40, 0, 0, 0, 0]   ! these are hours post irradiation, incremented when irradiation time is known (in ReadProtocol)
elseif (mod(iphase_hours,10) == 7) then    ! this is the compute_cycle case for multiple times, no PEST
    compute_cycle = .true.
    normalise = .true.
    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
    nphase_hours = 49
    next_phase_hour = 1
    do j = 1,49
        phase_hour(j) = (j-1)*0.5
    enddo
!    phase_hour(1:25) = [0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0,24.0]   ! these are hours post irradiation, incremented when irradiation time is known (in ReadProtocol)   
elseif (mod(iphase_hours,10) == 8) then    ! this is the compute_cycle case for M%-only experiments
    compute_cycle = .true.
    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
    nphase_hours = 5
    next_phase_hour = 1
    phase_hour(1:5) = [1.0, 1.5, 2.5, 3.5, 4.5]   ! these are hours post irradiation, incremented when irradiation time is known (in ReadProtocol)
else
    if (use_PEST) then
        write(*,*) 'Error: ReadMcParameters: with PEST iphase_hours must be 1 - 8'
        write(nflog,*) 'Error: ReadMcParameters: with PEST iphase_hours must be 1 - 8'
        stop
    endif
endif
icase = iphase_hours/10
if (icase == 10) normalise = .true.
if (icase == 11) M_only = .true.
if (icase == 12) then
    normalise = .true.
    M_only = .true.
endif
write(*,*) 'iphase_hours: ',iphase_hours,'    normalise: ',normalise, '    M_only: ',M_only
write(*,*) 'nphase_hours: ',nphase_hours
write(*,*) 'phase_hour:'
write(*,'(10f6.2)') phase_hour(1:nphase_hours)
write(*,*) 'did ReadMcParameters:'
!repRate(1:2) = 0
!write(*,*) '====================================================================='
!write(*,*) 'Test!!!!!!!!!!!!!!!    SETTING repRate = 0 for NHEJ'
!write(*,*) '====================================================================='

! Initialise CP delay totals
tCPdelay = 0 
tATMdelay = 0
tATRdelay = 0
! Initialise phase CP delay totals
totG1delay = 0
totSdelay = 0
totG2delay = 0
nG1delay = 0
nSdelay = 0
nG2delay = 0

! Initialise misjoining counts
misjoins = 0

! Set use of pATM in G1
!use_G1_stop = use_G1_pATM  ! this overrides the default, which was set above (to use ATM_act in G1)
end subroutine

!!--------------------------------------------------------------------------
!! DNA-PK inhibition factor.  C is inhibitor drug concentration
!!--------------------------------------------------------------------------
!function finhibit(C) result(fin)
!real(8) :: C, fin
!
!if (use_exp_inhibit) then
!    fin = KmaxInhibit*(1 - exp(-b_exp*C))
!else
!    fin = KmaxInhibit*C/(b_hill + C)
!endif
!end function

!--------------------------------------------------------------------------
! DNA-PK inhibition rate factor.  C is inhibitor drug concentration
!--------------------------------------------------------------------------
!function inhibitRate(C) result(fin)
!real(8) :: C, fin
!
!if (use_exp_inhibit) then
!    fin = KmaxInhibitRate*(1 - exp(-b_exp*C))
!else
!    fin = KmaxInhibitRate*C/(b_hill + C)
!endif
!end function

!--------------------------------------------------------------------------
! repRate:          depends on pathway
! fidRate:          depends on phase, pathway
! misrepRate:       depends on phase (pathway?)
! pathwayUsed:      depends on phase, pathway
!--------------------------------------------------------------------------
!subroutine setupRadiation
!
!pathwayUsed = .false.
!pathwayUsed(G1_phase,1:2) = .true.
!pathwayUsed(S_phase,1:4) = .true.
!pathwayUsed(G2_phase,3:4) = .true.
!
!if (use_inhibiter) then
!    pathwayUsed(G1_phase,5) = .true.
!    pathwayUsed(S_phase,5) = .true.
!    pathwayUsed(G2_phase,5) = .true.
!endif
!end subroutine

!--------------------------------------------------------------------------
! f is the fractional progress through S-phase
! DSB0[:] is the initial number of DSBs on each repair pathway
! totDSB0 is the total initial DSB count
! If the pathway is NHEJ, a fraction of the DSBs that would be on the NHEJ
! pathway are instead switched to the MMEJ pathway when the DNA-PK inhibitor
! drug is present.
!
! Option 1: the first method, which used one pHR for S and one for G2, and
! pHRfast_max = pHR_S/2, pHRslow_max = pHR_S/2
! pHRs_G2 = pHR_G2/2, pHRc_G2 = pHR_G2/2 with pHR_G2 = 0.15 (fixed)
!
! Option 2: Bill's suggestion, making the pHRs and pHRc parameters differ,
! and preserving the total complex fraction.
!
! Note: pHRs + pHRc must be <= 1
! 
! 30/05/2023
! Add distinction between pre- and post-DNA relication DSBs.
! pre = DSB(:,1), post = DSB(:,2)
!--------------------------------------------------------------------------
subroutine cellIrradiation(cp, dose, Cdrug)
type(cell_type), pointer :: cp
real(8) :: dose, Cdrug
integer :: ityp, kpar = 0
integer :: phase
real(8) :: DSB0(NP,2)
real(8) :: totDSB0, baseDSB, fin, T_S, f_S, NG1, NNHEJ, pHRs, pHRc, pHR, pNHEJ, NHRs, NHRc
real(8) :: Pbase, Pdie, R
real(8) :: DSB_Gy = 35
real(8) :: th, Npre, Npre_s, Npre_c, Npost, Npost_s, Npost_c, Pc, x
integer, parameter :: option = 2
type(cycle_parameters_type),pointer :: ccp
logical :: use_Jeggo = .true.

cp%irradiated = .true.
ccp => cc_parameters(1)
phase = cp%phase
cp%phase0 = min(phase, M_phase)
NG1 = DSB_Gy*dose
DSB0 = 0

if (use_Jeggo) then
    If (phase == G1_phase) Then
        f_S = 0
    ElseIf (phase == S_phase) Then
        f_S = cp%progress
!        th = (istep*DELTA_T - cp%t_S_phase)/3600  ! time since start of S_phase (h)
         th = cp%progress*ccp%T_S/3600
    ElseIf (phase >= G2_phase) Then
        f_S = 1
!        th = (istep*DELTA_T - cp%t_S_phase)/3600  ! time since start of S_phase (h)
         th = (ccp%T_S + cp%progress*ccp%T_G2)/3600
        
    End If
    totDSB0 = (1 + f_S) * NG1
    DSB0(TMEJ,:) = 0
    If (phase == G1_phase) Then
        Npre = totDSB0
    Else
        Npre = NG1 * (1 - f_S)
    End If
    Npost = totDSB0 - Npre
    If (f_S > 0) Then
        pHRs = fIliakis*pHRs_max * ((1 - rmin) * fdecay(th) + rmin)
        pHRc = fIliakis*pHRc_max * ((1 - rmin) * fdecay(th) + rmin)
    else
        pHRs = 0
        pHRc = 0
    End If
    pHR = (1 - pComplex)*pHRs + pComplex*pHRc
    if (kcell_now == 1) then
        write(*,'(a,2f8.3)') 'th, fdecay(th): ',th, fdecay(th)
        write(*,'(a,3f8.3)') 'fIliakis,pHRs,pHRc: ',fIliakis,pHRs,pHRc
    endif
    if (option == 1) then
        pNHEJ = 1 - pHR
        DSB0(NHEJfast,1) = (1 - pComplex) * Npre
        DSB0(NHEJfast,2) = pJeggo * pNHEJ * Npost     ! fast
        DSB0(NHEJslow,1) = pComplex * Npre
        DSB0(NHEJslow,2) = (1 - pJeggo) * pNHEJ * Npost     ! slow
        DSB0(HR,1) = 0
        DSB0(HR,2) = pHR * Npost
        if (phase == G2_phase) then
            nslow_sum = nslow_sum + 1
            pHR_sum = pHR_sum + pHR
            pNHEJslow_sum = pNHEJslow_sum + (1 - pJeggo) * pNHEJ
            fdecay_sum = fdecay_sum + fdecay(th)
        endif
    elseif (option == 2) then
        DSB0(NHEJfast,1) = (1-pComplex)*Npre
        DSB0(NHEJfast,2) = (1-pComplex)*(1-pHRs)*Npost     ! fast
        DSB0(NHEJslow,1) = pComplex*Npre
        DSB0(NHEJslow,2) = pComplex*(1-pHRc)*Npost     ! slow
        DSB0(HR,1) = 0
        DSB0(HR,2) = pHR*Npost
        if (kcell_now == 1) then
            write(*,*) 'Npre, Npost: ',Npre,Npost
            write(*,'(2(a,3f8.1))') 'pre  DSB at IR: ',DSB0(1:3,1),'  NNHEJ: ',sum(DSB0(1:2,1))
            write(*,'(2(a,3f8.1))') 'post DSB at IR: ',DSB0(1:3,2),'  NNHEJ: ',sum(DSB0(1:2,2))
            write(*,'(a,3f8.1)')    'total DSB at IR: ',sum(DSB0(NHEJfast,:)),sum(DSB0(NHEJslow,:)),sum(DSB0(HR,:))
        endif
    endif
else
    if (phase == G1_phase) then
        f_S = 0.0
        totDSB0 = NG1
        DSB0(NHEJslow,1) = Pcomplex*totDSB0
        DSB0(NHEJfast,1) = (1 - Pcomplex)*totDSB0
    elseif (phase == S_phase) then
        th = (istep*DELTA_T - cp%t_S_phase)/3600  ! time since start of S_phase (h)
        f_S = cp%progress
        totDSB0 = (1+f_S)*NG1
        pHRs = pHRs_max*((1-rmin)*fdecay(th) +rmin)
        pHRc = pHRc_max*((1-rmin)*fdecay(th) +rmin)
        Npre = NG1*(1 - f_S)
        Npost = NG1*2*f_S
        Npre_s = Npre*(1 - pComplex)
        Npre_c = Npre*pComplex
        Npost_s = Npost*(1 - pComplex)
        Npost_c = Npost*pComplex
        NHRs = Npost_s*pHRs
        NHRc = Npost_c*pHRc
        DSB0(HR,1) = NHRs + NHRc
        DSB0(NHEJfast,1) = Npre_s + Npost_s - NHRs  ! not correct
        DSB0(NHEJslow,1) = Npre_c + Npost_c - NHRc
    elseif (phase >= G2_phase) then
        f_S = 1
        th = (istep*DELTA_T - cp%t_S_phase)/3600  ! time since start of S_phase (h)
        totDSB0 = 2*NG1
        Npost = totDSB0
        pHRs = pHRs_max*((1-rmin)*fdecay(th)+rmin)
        DSB0(HR,2) = Npost*pHRs
        DSB0(NHEJfast,2) = Npost*(1 - pHRs)
        DSB0(NHEJslow,2) = 0
    endif
endif    

cp%pATM = 0
cp%pATR = 0
if (phase == G1_phase) then
    ! Apoptosis in G1
    Pbase = exp(-baseRate*totDSB0)   ! this is 1 - Pdie
    Pdie = 1 - Pbase 
    if (single_cell) Pdie = 0
    R = par_uni(kpar)
    if (R < Pdie) then  ! cell dies of apoptosis
        cp%state = DEAD
        Napop = Napop + 1
        Ncells = Ncells - 1
        ityp = cp%celltype
        Ncells_type(ityp) = Ncells_type(ityp) - 1
        Nviable(ityp) = Ncells_type(ityp)
        Ndead(ityp) = Ndead(ityp) + 1
!        write(*,*) 'apoptotic death: ',kcell_now, phase
    endif
elseif (phase == G2_phase) then
    if (use_Jaiswal) then
        ! nothing needed to be done here
        !if (cp%CC_act > 0) write(*,*) 'Cell in G2 at IR, CC_act: ',kcell_now, cp%CC_act
    elseif (use_D_model) then
        cp%pATM = K_ATM(3,1)*totDSB0/(K_ATM(3,2) + totDSB0)
        cp%pATR = K_ATR(3,1)*totDSB0/(K_ATR(3,2) + totDSB0)
        if (cp%pATR > 0) then
            write(*,*) 'cp%pATR: ', kcell_now,K_ATR(3,1),cp%pATR
            stop
        endif
    endif
endif

cp%DSB = DSB0
cp%totDSB0 = totDSB0
cp%Nlethal = 0
!if (kcell_now <= 10) write(*,'(a,2i6,8f8.2)') 'IR: kcell, phase,DSB,f_S: ',kcell_now,phase,cp%DSB(1:4),f_S

totPmit = 0
totPaber = 0
tottotDSB = 0
totNlethal = 0

totNDSB = 0
totNmisjoins = 0
totSFfactor = 0

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
function fdecay(t) result(f)
real(8) :: t, f

if (use_sigmoid) then
    f = fsigmoid(t)
else
    f = exp(-kdecay*t)
endif
end function

!--------------------------------------------------------------------------
! From Richards's curve
! t is hours since start of S-phase, csig approx = T_S
!--------------------------------------------------------------------------
function fsigmoid(t) result(f)
real(8) :: t, f

f = 1.0/(1.0 + exp(ksig*(t - csig)))
end function


!--------------------------------------------------------------------------
! There is no need to work with pATM concentration, since the concentration
! is proportional to the mass.
! Let M = mass of pATM, r = mass rate of production, k = decay rate constant
! dM/dt = r - kM
! dM/(r - kM) = dt
! dM/(M - r/k) = -k.dt
! log(M - r/k) = -kt + c
! M(t) - r/k = A.exp(-kt)
! M(0) = r/k + A => A = M(0) - r/k
! M(t) = r/k + (M(0) - r/k)exp(-kt)
! This assumes that r is constant, therefore it is valid only for small t,
! because ATM_DSB is constantly changing.
! Note: parameters from mcradio assume time is hours, therefore the time
! step passed to the subroutine, dth, has been converted to hours.
! 11/02/22 changed K_ATM(2) & (3) to (1) & (2)
! Note that as k -> 0, M(t) -> M(0) + lim(rt*(1 - exp(-kt))/kt)
! and limit as x->0, lim(1-exp(-x)/x = lim(x - x2/2 + x^3/6 - ...)/x = 1
! therefore M(t) -> M(0) + rt as expected.
!--------------------------------------------------------------------------
subroutine updateATM(iph,pATM,ATM_DSB,dth)
integer :: iph
real(8) :: pATM, ATM_DSB, dth, pATM0
real(8) :: r, t, fz, kdecay
logical :: OK

!write(*,*) 'updateATM: iph: ',iph
OK = .true.
if (iph == G1_phase .and. .not.use_G1_pATM) OK = .false.
if (iph == S_phase .and. .not.use_S_pATM) OK = .false.
if (iph >= G2_phase) OK = .false.
if (.not.OK) then
    write(*,*) 'ERROR: updateATM: should not get here with iph, use_G1_pATM, use_S_pATM: ',iph,use_G1_pATM,use_S_pATM
    stop
endif

pATM0 = pATM
if (use_G2_pATM_Nindependent .and. iph == G2_phase) then
    r = K_ATM(iph,1)   ! rate of production of pATM
else
    r = K_ATM(iph,1)*ATM_DSB   ! rate of production of pATM
endif
t = (tnow - t_irradiation)/3600
kdecay = max(1.0e-8,K_ATM(iph,2))  ! decay rate constant
pATM = r/kdecay + (pATM - r/kdecay)*exp(-kdecay*dth)
if (isnan(pATM)) then
    write(*,*) 'NAN in updateATM: ',iph, r, kdecay, r/kdecay
    stop
endif
!if (kcell_now <= 10) write(*,'(a,i4,6f6.2)') 'updateATM: r,k,r/k,ATM_DSB,pATM: ',kcell_now,r,kdecay,r/kdecay,ATM_DSB,pATM
!write(*,'(a,i3,6f8.4)') 'updateATM: ',iph,ATM_DSB,r,k,r/k,pATM
!if (iph == 2) stop
end subroutine

!--------------------------------------------------------------------------
! Try making pATR production rate saturate, or limit pATR?
! Note: parameters from mcradio assume time is hours, therefore the time
! step passed to the subroutine, dth, has been converted to hours.
! 11/02/22 changed K_ATR(2) & (3) to (1) & (2)
!--------------------------------------------------------------------------
subroutine updateATR(iph,pATR,ATR_DSB,dth)
integer :: iph
real(8) :: pATR, ATR_DSB, dth
real(8) :: r, t, fz, k1, k2, x, xmax, km, xf

if (iph == G2_phase .and. use_Jaiswal) return
k1 = K_ATR(iph,1)
r = k1*ATR_DSB   ! rate of production of pATR
t = (tnow - t_irradiation)/3600
if (t > Ztime) then
    fz = 1
else
    fz = t/Ztime
endif
r = fz*r
k2 = K_ATR(iph,2)  ! decay rate constant
if ((iph == 3) .and. use_D_model) then
    r = 0
    k2 = max(1.0e-12,K_ATR(iph,3))  ! decay rate constant
endif
if (k2 == 0) then
    write(*,*) 'updateATR: k2 = 0: iph,K_ATR(iph,2): ',iph,K_ATR(iph,2)
    stop
endif
pATR = r/k2 + (pATR - r/k2)*exp(-k2*dth)
!if (kcell_now == 1) write(nfphase,'(a,i4,6f8.3)') 'updateATR: k1,r,k2,r/k2,ATR_DSB,pATR: ',kcell_now,k1,r,k2,r/k2,ATR_DSB,pATR
!if (pATR > 0) write(*,'(a,i4,5f8.4)') 'updateATR: r,k,r/k,ATR_DSB,pATR: ',kcell_now,r,k,r/k,ATR_DSB,pATR
end subroutine

!------------------------------------------------------------------------
! Time computed in hours, returned in secs 
! CP delay determined by ATM_act or pATM.
!------------------------------------------------------------------------
function G1_checkpoint_time(cp) result(t)
type(cell_type), pointer :: cp
real(REAL_KIND) :: t, th, atm
integer :: iph = 1

if (.not.use_G1_stop) then
    write(*,*) 'Error: G1_checkpoint_time: should not get here with use_G1_stop = ',use_G1_stop
    stop
endif
if (use_DSB_CP) then
    t = 0
    return
endif
if (use_G1_CP_factor) then
    t = G1_CP_time
    return
endif
if (use_G1_pATM) then
    atm = cp%pATM
else
    atm = cp%ATM_act
endif
th = K_ATM(iph,3)*(1 - exp(-K_ATM(iph,4)*atm))
!if (kcell_now == 1) write(*,'(a,i6,2f8.3)') 'G1_checkpoint_time (h): ',kcell_now,atm,th
if (isnan(th)) then
    write(*,*) 'NAN in G1_checkpoint_time: ',atm
    stop
endif
totG1delay = th + totG1delay
nG1delay = nG1delay + 1
t = 3600*th
end function

!------------------------------------------------------------------------
! Combined effect of ATM_act (or pATM) and ATR_act (was pATR) in S
! Time computed in hours, returned in secs
! CP delay is the sum of delays created by ATM_act and by ATR_act
! Effect of ATR_act is currently turned off.
!------------------------------------------------------------------------
function S_checkpoint_time(cp) result(t)
type(cell_type), pointer :: cp
real(REAL_KIND) :: t, th, th_ATM, th_ATR, atm
integer :: iph = 2
logical :: use_ATR = .false.

if (.not.use_S_stop) then
    write(*,*) 'Error: S_checkpoint_time: should not get here with use_S_stop = ',use_S_stop
    stop
endif
if (use_S_pATM) then
    atm = cp%pATM
else
    atm = cp%ATM_act
endif
!th_ATM = K_ATM(iph,3)*(1 - exp(-K_ATM(iph,4)*cp%pATM))
th = K_ATM(iph,3)*(1 - exp(-K_ATM(iph,4)*atm))
if (use_ATR) then
    th_ATR = K_ATR(iph,3)*(1 - exp(-K_ATR(iph,4)*cp%pATR))
else
    th_ATR = 0
endif
th = th_ATM + th_ATR
if (kcell_now <= 10) then
    write(nfphase,*)
    write(nfphase,'(a,i6,3f6.2)') 'S_checkpoint_time: ',kcell_now,th_ATM,th_ATR,th
!    write(nfphase,'(a,4f8.3)') 'pATM,katm3g2,katm4g2,1-exp(-katm4g2*pATM): ',cp%pATM,K_ATM(iph,3),K_ATM(iph,4),1-exp(-K_ATM(iph,4)*cp%pATM)
!    write(nfphase,'(a,4f8.3)') 'pATR,katr3g2,katr4g2,1-exp(-katr4g2*pATR): ',cp%pATR,K_ATR(iph,3),K_ATR(iph,4),1-exp(-K_ATR(iph,4)*cp%pATR)
endif
t = 3600*th
end function

!------------------------------------------------------------------------
! Combined effect of pATM and pATR in G2
! Time computed in hours, returned in secs
! CP delay is the sum of delays created by pATM and by pATR
!------------------------------------------------------------------------
function G2_checkpoint_time(cp) result(t)
type(cell_type), pointer :: cp
real(REAL_KIND) :: t, th, th_ATM, th_ATR, N_DSB, k3, k4
integer :: iph = 3
real(REAL_KIND) :: G2_DSB_threshold = 15

if (use_Jaiswal) then
    write(*,*) 'ERROR: G2_checkpoint_time not used with Jaiswal'
    stop
endif
if (use_D_model) then
    th_ATM = cp%pATM
!    th_ATR = cp%pATR
    th_ATR = 0  ! if effect of pATR is disabled
else
    if (use_DSB_CP) then
        N_DSB = sum(cp%DSB)
        k3 = K_ATR(iph,3)
        k4 = K_ATR(iph,4)
        th = k3*N_DSB/(k4 + N_DSB)
        if (kcell_now < 10) write(*,'(a,4f8.3)') 'N_DSB,k3,k4,th: ',N_DSB,k3,k4,th
        t = 3600*th
        return
    endif
    if (negligent_G2_CP .and. (sum(cp%DSB) < G2_DSB_threshold)) then
        th_ATM = 0
    else
        th_ATM = K_ATM(iph,3)*(1 - exp(-K_ATM(iph,4)*cp%pATM))
    endif
    th_ATR = K_ATR(iph,3)*(1 - exp(-K_ATR(iph,4)*cp%pATR))
endif
!if (kcell_now <= 100) then
!    write(nflog,'(a,i6,4f8.4)') 'cell, pATM, pATR, th_ATM, th_ATR: ',kcell_now, cp%pATM,cp%pATR,th_ATM,th_ATR
!    write(*,'(a,i6,4f8.4)') 'cell, pATM, pATR, th_ATM, th_ATR: ',kcell_now, cp%pATM,cp%pATR,th_ATM,th_ATR
!endif
th = th_ATM + th_ATR
totG2delay = th + totG2delay
nG2delay = nG2delay + 1
if (kcell_now <= 10) then
    write(nfphase,*)
    write(nfphase,'(a,i6,3f6.2)') 'G2_checkpoint_time: ',kcell_now,th_ATM,th_ATR,th
    write(nfphase,'(a,4f8.3)') 'pATM,katm3g2,katm4g2,1-exp(-katm4g2*pATM): ',cp%pATM,K_ATM(iph,3),K_ATM(iph,4),1-exp(-K_ATM(iph,4)*cp%pATM)
    write(nfphase,'(a,4f8.3)') 'pATR,katr3g2,katr4g2,1-exp(-katr4g2*pATR): ',cp%pATR,K_ATR(iph,3),K_ATR(iph,4),1-exp(-K_ATR(iph,4)*cp%pATR)
endif
t = 3600*th
if (is_radiation) then
    G2thsum(1) = G2thsum(1) + th_ATM
    G2thsum(2) = G2thsum(2) + th_ATR
    NG2th = NG2th + 1
endif
end function

!------------------------------------------------------------------------
! Only G2 is affected by pATR
! Best to turn off pATR in S-phase by setting katr1s = katr3s = 0
! Now slowdown is acting in G1-phase and S-phase, and only affected by pATM
! fslow in computed for G2, but not used because G2 is controlled by Jaiswal.
! Slowdown must occur only for cells that were irradiated.  
! This means for cells that were present at IR.  In other words:
! Initially, every cell has cp%irradiated = false
! At IR, every cell has cp%irradiated = true
! At cell division, each daughter cell should have cp%irradiated = false
! Fixed: 30/03/2023
!------------------------------------------------------------------------
subroutine get_slowdown_factors(cp,fATM,fATR)
type(cell_type), pointer :: cp
integer :: iph
real(REAL_KIND) :: fATM, fATR
real(REAL_KIND) :: pATM, pATR, k3, k4, N_DSB, atm, atr
logical :: OK

if (.not.cp%irradiated) then
    fATM = 1
    fATR = 1
    return
endif
iph = cp%phase
if (iph > G2_phase) then
    write(*,*) 'Error: get_slowdown_factors called with iph = ',iph
    stop
endif
OK = .true.
if ((iph == G1_phase) .and. use_G1_stop) OK = .false.
if ((iph == S_phase) .and. use_S_stop) OK = .false.
if (.not.OK) then
    write(*,*) 'Error: get_slowdown_factors: iph, use_G1_stop, use_S_stop: ',iph, use_G1_stop, use_S_stop
    stop
endif

atr = 0
fATR = 1
if (use_DSB_CP) then
    N_DSB = sum(cp%DSB)
    k3 = K_ATM(iph,3)
    k4 = K_ATM(iph,4)
!    if (kcell_now == 1) write(*,'(a,i6,4f8.3)') 'k3,k4,N_DSB: ',kcell_now,k3,k4,N_DSB,k3*N_DSB/(k4 + N_DSB)
    fATM = max(0.01,1 - k3*N_DSB/(k4 + N_DSB))
    return
endif
if (iph == G1_phase) then
!    atm = 0
!    if (kcell_now == 1) write(*,*) '**** G1 atm = 0 in get_slowdown_factors ****'
    if (use_G1_pATM ) then
        atm = cp%pATM
    else
        atm = cp%ATM_act
    endif
elseif (iph == S_phase) then
!    atm = 0
!    if (kcell_now == 1) write(*,*) '**** S atm = 0 in get_slowdown_factors ****'
    if (use_S_pATM ) then
        atm = cp%pATM
    else
        atm = cp%ATM_act
        if (use_ATR_S) atr = cp%ATR_act
    endif
endif
k3 = K_ATM(iph,3)
k4 = K_ATM(iph,4)
if (single_cell) then
    write(*,'(a,3e12.3)') 'k3,k4,ATM_act: ',k3,k4,cp%ATM_act
    write(nflog,'(a,3e12.3)') 'k3,k4,ATM_act: ',k3,k4,cp%ATM_act
endif
!if (single_cell) write(nflog,'(a,i6,3e12.4)') 'get_slowdown_factors: kcell,k3,k4,atm: ',kcell_now,k3,k4,atm
if (use_exp_slowdown) then
    fATM = exp(-k4*atm)
else
    if ((k4 + atm) > 0) then
        fATM = max(0.01,1 - k3*atm/(k4 + atm))
    else
        fATM = 1.0
    endif
endif
if (iph == S_phase .and. use_ATR_S) then
    k3 = K_ATR(iph,3)   !*G2_katr3_factor
    k4 = K_ATR(iph,4)   !*G2_katr4_factor
    fATR = max(0.01,1 - k3*atr/(k4 + atr))
endif
!write(*,'(a,i3,4f8.4)') 'iph,fATR,fATM: ',iph,atr,atm,fATR,fATM
!if (iph == 2 .and. kcell_now <= 10) then
!    write(*,'(a,i3,5f8.4)') 'iph,fATR,fATM,fslow,: ',iph,atr,atm,fATR,fATM,fATR*fATM
!    write(nflog,'(a,i3,5f8.4)') 'iph,fATR,fATM,fslow,: ',iph,atr,atm,fATR,fATM,fATR*fATM
!endif
end subroutine

!------------------------------------------------------------------------
! Whether a phase uses CP delay or slowdown is handled here.
!------------------------------------------------------------------------
function slowdown(cp) result(fslow)
type(cell_type), pointer :: cp
real(REAL_KIND) :: fslow
type(cycle_parameters_type),pointer :: ccp
integer :: iph
real(REAL_KIND) :: fATM, fATR, dt, dtCPdelay

ccp => cc_parameters(1)
if (use_fixed_CP) then
    ! Testing with fixed checkpoint delays, specified in drm_monolayer_many
    if (cp%phase == G1_phase) then
        fslow = ccp%T_G1/(ccp%T_G1 + G1_delay)
    elseif (cp%phase == S_phase) then
        fslow = ccp%T_S/(ccp%T_S + S_delay)
    elseif (cp%phase == G2_phase) then
        fslow = ccp%T_G2/(ccp%T_G2 + G2_delay)
    endif
else
    if (use_phase_dependent_CP_parameters) then
        iph = cp%phase
    else
        iph = 1
    endif
    if ((iph == 1 .and. use_G1_stop) .or. &
        (iph == 2 .and. use_S_stop) .or. &
        (iph == 3 .and. (use_G2_stop .or. use_jaiswal))) then
        fslow = 1.0
    else
        dt = DELTA_T
        call get_slowdown_factors(cp,fATM, fATR)
        if (use_addATMATRtimes) then
            dt = DELTA_T
            fslow = max(0.0,fATM + fATR - 1)
        else
            fslow = fATM*fATR
            if (single_cell) write(nflog,'(a,3f8.3)') 'fATM, fATR, fslow: ',fATM, fATR, fslow
            if (kcell_now == -3) then
                write(*,'(a,i6,i4,3f6.3)') 'fslow: ',kcell_now,iph,fATM,fATR,fslow
                write(nflog,'(a,i6,i4,3f6.3)') 'fslow: ',kcell_now,iph,fATM,fATR,fslow
            endif
        endif
        dtCPdelay = dt*(1 - fslow)
        totSdelay = totSdelay + dtCPdelay
!        if (kcell_now == 10) write(*,'(a,i4,4f8.3)') 'slowdown: ',kcell_now,fATM,fATR,fslow,cp%progress
    endif
endif
!if (kcell_now == 1) write(*,'(a,2i6,3f6.3)') 'kcell, phase, fslow: ',kcell_now, cp%phase, fslow, fATM, fATR
end function

!------------------------------------------------------------------------
!------------------------------------------------------------------------
subroutine get_CP_delay(cp)
type(cell_type), pointer :: cp

if (cp%phase == G1_phase) then
    cp%CP_delay = G1_checkpoint_time(cp)
!    if (kcell_now <= 100) write(*,'(a,i6,f8.3)') 'G1 CP_delay: ',kcell_now,cp%CP_delay/3600
elseif (cp%phase == S_phase) then
    cp%CP_delay = S_checkpoint_time(cp)
elseif (cp%phase == G2_phase .and. .not.use_Jaiswal) then
    cp%CP_delay = G2_checkpoint_time(cp)
!    write(*,'(a,i6,f8.3)') 'G2 CP_delay: ',kcell_now,cp%CP_delay/3600
!    write(nflog,'(a,i6,f8.3)') 'G2 CP_delay: ',kcell_now,cp%CP_delay/3600
endif
end subroutine

!------------------------------------------------------------------------
! Note: DNA repair is happening in updateRepair (pathwayRepair), and since 
! misrepair is computed at the same time, best to leave that as it is.  This
! means that DSB(:) is being updated with the full time step, not using the
! reduced time step employed for the Jaiswal equations.  Therefore damage D
! is assumed to be fixed for this updating.
!------------------------------------------------------------------------
subroutine Jaiswal_update(cp, dth)
type(cell_type), pointer :: cp
real(8) :: dth
real(8) :: dt = 0.001
real(8) :: D_ATR, D_ATM, CC_act, ATR_act, ATM_act, CC_inact, ATR_inact, ATM_inact, CC_act0
real(8) :: dCC_act_dt, dATR_act_dt, dATM_act_dt, t, t_G2, Kkcc2a, DSB(NP)
integer :: iph, it, Nt
type(cycle_parameters_type),pointer :: ccp
logical :: dbug

iph = cp%phase
if (iph >= 7) iph = iph - 6     ! checkpoint phase numbers --> phase number, to continue ATM and ATR processes through checkpoints
if (iph > G2_phase) then
    return
endif
Nt = int(dth/dt + 0.5)
do it = 1,NP
    DSB(it) = sum(cp%DSB(it,:))
enddo
dbug = (iph == 2 .and. (kcell_now <= 0))
if (iph == G1_phase) then
    D_ATM = DSB(NHEJslow)*norm_factor
    ATM_act = cp%ATM_act
elseif (iph == S_phase) then
    D_ATM = (DSB(HR) + DSB(NHEJslow))*norm_factor
    ATM_act = cp%ATM_act
    if (use_ATR_S) then
        D_ATR = DSB(HR)*norm_factor
        ATR_act = cp%ATR_act
    endif
elseif (iph == G2_phase) then
    D_ATR = DSB(HR)*norm_factor
    D_ATM = (DSB(HR) + DSB(NHEJslow))*norm_factor
    CC_act = cp%CC_act
    CC_act0 = CC_act
    ATR_act = cp%ATR_act
    ATM_act = cp%ATM_act
    t_G2 = (tnow - cp%t_start_G2)/3600
    if (t_G2 > 40) then
        ATR_act = 0
    elseif (t_G2 > 30) then
        ATR_act = ATR_act*(t_G2 - 30)/(40 - 30)
    endif
!    if (kcell_now == 1) write(*,*) '**** G2 D_ATM = D_ATR = 0 in jaiswal_update ****'
!    D_ATR = 0
!    D_ATM = 0
else
    return
endif
if (single_cell) write(nflog,*) 'CC_act: ',CC_act
if (use_km10_kcc2a_dependence) then
    Kkcc2a = cp%Kcc2a
else
    Kkcc2a = Kcc2a
endif
do it = 1,Nt
    ATM_inact = ATM_tot - ATM_act
    ATR_inact = ATR_tot - ATR_act
    if (iph == G2_phase) then
        CC_inact = CC_tot - CC_act
        ATR_inact = ATR_tot - ATR_act
        dCC_act_dt = (Kkcc2a + CC_act) * CC_inact / (Kmccp + CC_inact) - cp%Kt2cc * ATM_act * CC_act / (Kmccmd + CC_act) - cp%Ke2cc * ATR_act * CC_act / (Kmccrd + CC_act)
        dATR_act_dt = Kd2e * D_ATR * ATR_inact / (Kmrp + ATR_inact) - Kcc2e * ATR_act * CC_act / (Kmrd + CC_act)
        CC_act = CC_act + dt * dCC_act_dt
        CC_act = max(CC_act, 0.0)
        ATR_act = ATR_act + dt * dATR_act_dt
        ATR_act = min(ATR_act, ATR_tot)
    elseif (use_ATR_S) then
        dATR_act_dt = Kd2e * D_ATR * ATR_inact / (Kmrp + ATR_inact) - Kcc2e * ATR_act * CC_act / (Kmrd + CC_act)
        ATR_act = ATR_act + dt * dATR_act_dt
        ATR_act = min(ATR_act, ATR_tot)
    endif
    dATM_act_dt = Kd2t * D_ATM * ATM_inact / (Kmmp + ATM_inact) - Kti2t * ATM_act / (Kmmd + ATM_act)    
    ATM_act = ATM_act + dt*dATM_act_dt
    ATM_act = min(ATM_act, ATM_tot)

    t = it*dt
enddo
if (dbug) then
!    write(*,'(a,2i4,4f8.4)') 'kcell,iph,ATR,ATM: ',kcell_now,iph,ATR_act,ATM_act
!    write(nflog,'(a,2i4,4f8.4)') 'kcell,iph,ATR,ATM: ',kcell_now,iph,ATR_act,ATM_act
endif
cp%ATM_act = ATM_act
if (iph == G2_phase) then
    cp%CC_act = CC_act
    cp%ATR_act = ATR_act
    cp%dCC_act_dt = dCC_act_dt
elseif (iph == S_phase .and. use_ATR_S) then
    cp%ATR_act = ATR_act
endif
t = t_simulation/3600.
end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
subroutine test_Jaiswal
real(8) :: t, dth
integer :: it, Nt
type(cell_type), pointer :: cp

dth = 0.1
!Nt = int(20./dth)
Nt = 1
cp => cell_list(1)
cp%CC_act = 0
cp%DSB = 0
cp%ATM_act = 0
cp%ATR_act = 0
t = 0
write(nflog,*) 'test_Jaiswal: Nt: ',Nt
do it = 1,Nt
    call Jaiswal_update(cp, dth)
    t = t + dth
    write(nflog,'(2f8.4)') t,cp%CC_act
enddo
end subroutine

!--------------------------------------------------------------------------
! To determine repair on a given pathway
! Using parameters from mcradio, dth in hours
!--------------------------------------------------------------------------
subroutine pathwayRepair(path, dth, N0, N)
integer :: path
real(8) :: dth, N0, N

if (dth >= 0) then
    N = N0*exp(-repRate(path)*dth*repRateFactor(path))
else
    N = 0
endif
end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
function misrepairRate(initialBreaks, finalBreaks, eta) result(Pmis)
real(8) :: initialBreaks, finalBreaks, eta, Pmis
real(8) :: repairedBreaks, atanNum, atanDen

repairedBreaks = initialBreaks-finalBreaks
if (repairedBreaks < 1E-10) then
    Pmis = 0
    return
endif
atanNum = sqrt(3.0)*eta*repairedBreaks
atanDen = 2 + eta*(2*initialBreaks*finalBreaks*eta + initialBreaks + finalBreaks) 
Pmis = 1 - 2 * atan(atanNum/atanDen) / (repairedBreaks*sqrt(3.0)*eta)

end function

!------------------------------------------------------------------------
! To use parameters from mcradio, need to convert time from secs to hours
!------------------------------------------------------------------------
subroutine updateRepair(cp,dt)
type(cell_type), pointer :: cp
real(8) :: dt
integer :: phase
real(8) :: DSB(NP,2), dNlethal
real(8) :: DSB0(NP,2)
real(8) :: totDSB0, totDSB, Pmis, Nmis, dNmis(NP), totDSBinfid0, totDSBinfid, ATR_DSB, ATM_DSB, dth, binMisProb
real(8) :: Cdrug, inhibrate, Nreassign
real(8) :: f_S, eta_NHEJ, eta_TMEJ, tIR, eta
real(8) :: th_since_IR
logical :: pathUsed(NP)
integer :: k, iph, jpp  ! jpp = 1 for pre, = 2 for post
logical :: use_DSBinfid = .true.
real(8) :: DSB_min = 1.0e-3
logical :: use_totMis = .true.      ! was false!
logical :: use_ATM != .false.
logical :: dbug
logical :: use_G1_tdelay = .false.
logical :: do_G1_Jaiswal

dth = dt/3600   ! hours
use_ATM = .not.use_fixed_CP
phase = cp%phase
!if (single_cell .and. phase > G2_phase) then
!    write(nflog,*) 'updateRepair: phase: ',phase
!    stop
!endif
iph = phase
DSB = cp%DSB
!write(*,'(a,i4,4f8.1)') 'phase,DSB: ',phase,DSB(1:4)
repRateFactor = 1.0
! DNA-PK inhibition, DSB reassignment
if (use_constant_drug_conc) then
    Cdrug = Caverage(MAX_CHEMO + DRUG_A)
else
    Cdrug = cp%Cin(DRUG_A)
endif
!inhibrate = inhibitRate(Cdrug)
!if (inhibrate > 0) then
!    do k = 1,NP-1
!        Nreassign = DSB(k)*inhibrate*dt
!        DSB(k) = DSB(k) - Nreassign
!        DSB(TMEJ) = DSB(TMEJ) + Nreassign
!    enddo
!endif
! Now the effect of inhibition is to reduce the repair rate.
! Chalf is the drug concentration that reduces repair rate by 1/2.  fRR = exp(-k*C/Chalf)
! fRR = 0.5 when C = Chalf, 0.5 = exp(-k), k = -log(0.5) = 0.693
if (Chalf == 0) then
    write(nflog,*) 'ERROR: Chalf = 0'
    stop
endif
repRateFactor(1:2) = exp(-0.693*Cdrug/Chalf) 
!if (kcell_now == 1) write(nflog,'(a,4f10.4)') 'Cdrug IC,EC,Chalf,fRR: ',Cdrug,Caverage(MAX_CHEMO + DRUG_A),Chalf,repRateFactor(1)
! Reassignment to pathway 4 is (tentatively) constant
! Preass is an input parameter = prob of reassignment per hour
if (Preass > 0 .and. phase >= S_phase) then
    do jpp = 1,2
    do k = 1,2  ! only NHEJ pathways
        Nreassign = DSB(k,jpp)*Preass*dth
        DSB(k,jpp) = DSB(k,jpp) - Nreassign
!        if (use_SSA) then
!            DSB(TMEJ) = DSB(TMEJ) + Nreassign*(1 - SSAfrac)
!            DSB(SSA) = DSB(SSA) + Nreassign*SSAfrac
!        else
            DSB(TMEJ,jpp) = DSB(TMEJ,jpp) + Nreassign
!        endif
    enddo
    enddo
endif
DSB0 = DSB     ! initial DSBs for this time step

! Revised approach with no fidelity, just misrejoining

totDSB0 = sum(DSB0)
!if (totDSB0 == 0) then
!    write(nflog,*) 'totDSB0 = 0'
!    return
!endif
!if (kcell_now == 1) write(*,'(a,2i6,3f8.2)') 'updateRepair: kcell, phase, DSB0: ',kcell_now,cp%phase,cp%DSB(1:3)

ATM_DSB = sum(DSB(NHEJslow,:)) + sum(DSB(HR,:))   ! complex DSB
ATR_DSB = sum(DSB(HR,:))

!if (iph >= 7) iph = iph - 6     ! checkpoint phase numbers --> phase number, to continue pATM and pATR processes through checkpoints
!if (iph <= G2_phase) then      ! not for 4 (M_phase) or 5 (dividing)
!    if (iph == G2_phase .and. use_Jaiswal) then
!        call Jaiswal_update(cp, dth)
!    else
!        call updateATM(iph,cp%pATM,ATM_DSB,dth)     ! updates the pATM mass through time step = dth (if use_Jaiswal, G1 and S)
!        if (iph > G1_phase) call updateATR(iph,cp%pATR,ATR_DSB,dth)     ! updates the pATR mass through time step = dth (if use_Jaiswal, S only)
!    endif
!endif
!if (use_G1_TDELAY) then
!    ! We want to update Jaiswal for a cell in G1 only if G1_tdelay has elapsed since IR
!    ! Current time is tnow (sec).  IR time is t_irradiation
!    ! Time since IR is th_since_IR(hours)
!    th_since_IR = (tnow - t_irradiation)/3600
!if (.not.(iph == G1_phase .and. th_since_IR < G1_tdelay)) then
!    call Jaiswal_update(cp, dth)
!else
!    if (G1_tdelay == 0) write(*,*) 'Error: no Jaiswal_update with G1_tdelay = 0'
!endif

if (iph == G1_phase) then
    if (use_G1_tdelay) then
        ! We want to update Jaiswal for a cell in G1 only if G1_tdelay has elapsed since IR
        ! Current time is tnow (sec).  IR time is t_irradiation
        ! Time since IR is th_since_IR(hours)
        th_since_IR = (tnow - t_irradiation)/3600
        do_G1_Jaiswal = (th_since_IR > G1_tdelay).and.(G1_tdelay > 0)
    else
        ! We want to update Jaiswal for a cell in G1 only if the cell has undergone mitosis post-IR
        do_G1_Jaiswal = (cp%birthtime > t_irradiation)      ! (cp%generation > 1)
    endif
endif
if (((iph == G1_phase).and.do_G1_Jaiswal).or.(iph == S_phase)) then
    call Jaiswal_update(cp,dth)
endif

if ((iph == 1 .and. use_G1_pATM) .or. (iph == 2 .and. use_S_pATM)) then 
    call updateATM(iph,cp%pATM,ATM_DSB,dth)     ! updates the pATM mass through time step = dth
endif

DSB = 0
do k = 1,NP
do jpp = 1,2
!   if (dbug .and. DSB0(k) > 0) write(*,*) 'pathwayRepair: k,DSB0(k): ',kcell_now,k,DSB0(k)
    call pathwayRepair(k, dth, DSB0(k,jpp), DSB(k,jpp))
    if (DSB(k,jpp) < DSB_min) DSB(k,jpp) = 0
enddo
enddo
! DSB0(k) is the count before repair, DSB(k) is the count after repair

! The following commented out code follows MEDRAS
#if 0
totDSB = sum(DSB)
totDSBinfid0 = 0
totDSBinfid = 0
do k = 1,NP
    if (DSB0(k) > 0 .and. fidRate(k) < 1.0) then
        totDSBinfid0 = totDSBinfid0 + DSB0(k)
        totDSBinfid = totDSBinfid + DSB(k)
    endif
enddo

! Testing eta dependence
!eta_G1 = eta_G2 ! -> Pmis = 0.01 - 0.07
!eta_G2 = eta_G1 ! -> Pmis = 0.1 - 0.127
! Therefore small eta -> small Pmis -> smaller Nmis (Pmis is ~0.05 for eta_G2, ~0.11 for eta_G1) -> bigger SF

if (phase == G1_phase) then
    eta = eta_G1
elseif (phase == S_phase) then
    eta = eta_G1 + cp%progress*(eta_G2 - eta_G1)
elseif (phase == G2_phase) then
    eta = eta_G2
endif
if (use_DSBinfid) then
    Pmis = misrepairRate(totDSBinfid0, totDSBinfid, eta)
else
    Pmis = misrepairRate(totDSB0, totDSB, eta)
endif
cp%totMis = cp%totMis + Pmis*(totDSB0 - totDSB)
binMisProb = cp%totMis/(cp%totDSB0 - totDSB)
if (single_cell) write(nflog,'(a,2e12.3)') 'Pmis, binMisProb: ',Pmis,binMisProb
Nmis = 0
dNmis = 0
do k = 1,NP
!    if (pathUsed(k) .and. fidRate(k) < 1.0) then
    if (fidRate(k) < 1.0) then
        if (use_totMis) then
            dNmis(k) = (DSB0(k) - DSB(k))*(1 - fidRate(k)*(1-binMisProb))
            Nmis = Nmis + dNmis(k)
        else
            dNmis(k) = (DSB0(k) - DSB(k))*(1 - fidRate(k)*(1-Pmis))
            ! Approx = (number of repairs)*(1 - fidRate) = NDSB*repRate*(1 - fidRate)
!            if (k == 2 .or. k == 4) write(nfphase,'(a,i4,4f8.4)') 'k: ',k,(DSB0(k) - DSB(k)),fidRate(k),Pmis,fidRate(k)*(1-Pmis)
            Nmis = Nmis + dNmis(k)
!            if (DSB0(k) > DSB(k)) then
!                write(nflog,'(a,i3,3f8.3,e12.3)') 'dNmis: ',k,(DSB0(k) - DSB(k)),Pmis,(1 - fidRate(k)*(1-Pmis)),dNmis
!            endif
        endif
    endif
enddo
write(nflog,'(a,6f8.3)') 'dNmis, Nmis: ',dNmis,Nmis
if (isnan(Nmis)) then
    write(*,*) 'updateRepair: Nmis isnan'
    stop
endif
#endif

if (phase == G1_phase) then
    f_S = 0.0
elseif (phase == S_phase) then
    f_S = cp%progress
elseif (phase >= G2_phase) then
    f_S = 1.0
endif
tIR = (t_simulation - t_irradiation)/3600   ! time since IR, in hours
tIR = max(tIR,0.0)
eta_NHEJ = eta_lookup(phase, NHEJfast, f_S, tIR) 
eta_TMEJ = eta_lookup(phase, TMEJ, f_S, tIR) 

Nmis = 0
! For NHEJ pathways
totDSB0 = sum(DSB0(NHEJfast,:)) + sum(DSB0(NHEJslow,:))
totDSB = sum(DSB(NHEJfast,:)) + sum(DSB(NHEJslow,:))
Pmis = misrepairRate(totDSB0, totDSB, eta_NHEJ)
Nmis = Nmis + Pmis*(totDSB0 - totDSB)
misjoins(1) = misjoins(1) + Pmis*(totDSB0 - totDSB)
!write(nflog,'(a,4f10.4)') 'totDSB0, totDSB, eta_NHEJ, Pmis: ',totDSB0, totDSB, eta_NHEJ, Pmis
!if (tIR == 2.0) then
!    write(nflog,*) 'tIR: ',tIR
!    write(nflog,*) 'Varying eta to find how Pmis depends on eta'
!    do k = 1,11
!        eta = eta_NHEJ*(0.5 + (k-1)*0.1)
!        Pmis = misrepairRate(totDSB0, totDSB, eta)
!        write(nflog,'(a,i4,2f8.4)') 'k, eta, Pmis: ',k,eta,Pmis
!    enddo
!endif
if (sum(DSB0(TMEJ,:)) > 0) then
    ! For TMEJ pathway
    Pmis = misrepairRate(sum(DSB0(TMEJ,:)), sum(DSB(TMEJ,:)), eta_TMEJ)
    Nmis = Nmis + Pmis*(sum(DSB0(TMEJ,:)) - sum(DSB(TMEJ,:)))
    misjoins(2) = misjoins(2) + Pmis*(totDSB0 - totDSB)
endif
cp%DSB = DSB
!dNlethal = Klethal*misrepRate(phase)*Nmis   ! (was 0.5)  1.65 needs to be another parameter -> Klethal
dNlethal = Klethal*Nmis   ! Here Klethal ~ 2.1*0.19 = 0.4, i.e. using a single average misreprate
cp%Nlethal = cp%Nlethal + dNlethal

end subroutine

!------------------------------------------------------------------------
! Clonogenic survival probability at first mitosis (McMahon, mcradio)
! cp%phase0 is the cell's phase at IR
! Note that this assumes that cells died of apoptosis in G1 at baseRate
! (see cellIrradiation())
!------------------------------------------------------------------------
subroutine survivalProbability(cp)
type(cell_type), pointer :: cp
real(8) :: DSB(NP,2), totDSB(2), Nlethal,Nmisjoins, Paber, Pbase, Papop, Pmit, Psurv
integer :: k, jpp

DSB = cp%DSB
do jpp = 1,2
    totDSB(jpp) = sum(DSB(:,jpp))
    totNDSB(jpp) = totNDSB(jpp) + totDSB(jpp)
enddo
Nlethal = cp%Nlethal
if (klethal > 0) then
    Nmisjoins = Nlethal/Klethal
else
    Nmisjoins = 0
endif
totNmisjoins = totNmisjoins + Nmisjoins

Pmit = exp(-(mitRate(1)*totDSB(1) + mitRate(2)*totDSB(2)))
do k = 1,2
    totSFfactor(k) = totSFfactor(k) + exp(-(mitRate(k)*totDSB(k)))
enddo
if (cp%phase0 < M_phase) then   ! G1, S, G2
    Paber = exp(-Nlethal)
!    Pmit = exp(-mitRate*totDSB)
    cp%Psurvive = Pmit*Paber  
!    if (kcell_now == 1) write(*,'(a,2f8.4,2e12.3)') 'totDSB,Nlethal,Pmit,Paber: ',totDSB,Nlethal,Pmit,Paber  
else    ! M_phase or dividing
    Paber = 1
!    Pmit = exp(-mitRate*totDSB)
    cp%Psurvive = Pmit*Msurvival
!    write(nflog,'(a,i6,3f10.4)') 'IR in mitosis: ',kcell_now,totDSB,cp%totDSB0,Pmit
endif
totSFfactor(3) = totSFfactor(3) + Paber
if (kcell_now <= 10) then
!    write(*,*)
    write(*,'(a,2i3,3f8.3)') 'cell,phase0,totDSB,psurvive: ',kcell_now,cp%phase0,totDSB,cp%Psurvive
!    write(*,'(a,i3,2f8.2)') 'phase0,Nmisjoins,totDSB: ',cp%phase0,Nmisjoins,totDSB
!    write(*,'(a,3e12.3)') 'Paber,Pmit,Psurvive: ',Paber,Pmit,cp%Psurvive
!    write(*,*)
!    write(nflog,'(a,2f10.2)') 'Nlethal, totDSB: ',Nlethal,totDSB
endif
NPsurvive = NPsurvive + 1   ! this is the count of cells for which Psurvive has been computed
cp%mitosis_time = tnow      ! needed to determine if mitosis occurs before or after CA
!Psurvive(NPsurvive) = cp%Psurvive

ATMsum = ATMsum + cp%pATM
ATRsum = ATRsum + cp%pATR

totPmit = totPmit + Pmit
totPaber = totPaber + Paber
tottotDSB = tottotDSB + sum(totDSB)
totNlethal = totNlethal + Nlethal

if (Nlethal > NMDIST*ddist_Nlethal) then
    count_Nlethal(NMDIST) = count_Nlethal(NMDIST) + 1
else
    k = Nlethal/ddist_Nlethal + 1
    count_Nlethal(k) = count_Nlethal(k) + 1
endif
if (sum(totDSB) > NMDIST*ddist_totDSB) then
    count_totDSB(NMDIST) = count_totDSB(NMDIST) + 1
else
    k = sum(totDSB)/ddist_totDSB + 1
    count_totDSB(k) = count_totDSB(k) + 1
endif

end subroutine

!------------------------------------------------------------------------
! Get average DNA growth factor for cells in S-phase
!------------------------------------------------------------------------
subroutine get_DNA_synthesis_rate(DNA_rate)
real(8) :: DNA_rate
type(cell_type), pointer :: cp
integer :: kcell, iph, cnt
real(8) :: atm, k3m, k4m, fATM, atr, k3r, k4r, fATR, rate, rate_sum, atm_ave

!write(*,'(a)') 'get_DNA_synthesis_rate'
k3m = K_ATM(S_phase,3)
k4m = K_ATM(S_phase,4)
k3r = K_ATR(S_phase,3)
k4r = K_ATR(S_phase,4)
fATR = 1.0
cnt = 0
rate_sum = 0
atm_ave = 0
do kcell = 1,nlist
    cp => cell_list(kcell)
    iph = cp%phase
    if (iph == S_phase) then
        cnt = cnt + 1
        if (use_S_pATM) then
            atm = cp%pATM
        else
            atm = cp%ATM_act
        endif
        atm_ave = atm_ave + atm
!        fATM = max(0.01,1 - k3*atm/(k4 + atm))  ! 0.01 ??
        if (use_exp_slowdown) then
            fATM = exp(-k4m*atm)
        else
            fATM = max(0.01,1 - k3m*atm/(k4m + atm))
        endif
        if (use_ATR_S) then
            atr = cp%ATR_act
            fATR = max(0.01,1 - k3r*atr/(k4r + atr))
        endif
        rate = fATM*fATR
!        if (kcell<= 10) write(*,'(a,i4,5f8.4)') 'DNA_rate: ',kcell, atr, atm, fATR, fATM, rate
        rate_sum = rate_sum + rate
    endif
enddo
DNA_rate = rate_sum/cnt
atm_ave = atm_ave/cnt
!write(*,'(a,3f8.3,e12.3)') 'DNA growth rate factor: ',DNA_rate,k3,k4,atm_ave
end subroutine

#if 0
!------------------------------------------------------------------------
! Write average S-phase ATM_DSB, pATM, fATM
!------------------------------------------------------------------------
subroutine show_S_phase_statistics()
real(8) :: hour
type(cell_type), pointer :: cp
integer :: kcell, iph, cnt, nthour
real(8) :: pATM, k3, k4, fATM, fATM_ave, pATM_ave, ATM_DSB, ATM_DSB_ave

!write(nflog,'(a)') 'S_phase_statistics'
nthour = 3600/DELTA_T
hour = real(istep)/nthour
!k3 = K_ATM(S_phase,3)
!k4 = K_ATM(S_phase,4)
cnt = 0
ATM_DSB_ave = 0
fATM_ave = 0
pATM_ave = 0
do kcell = 1,nlist
    cp => cell_list(kcell)
    iph = cp%phase
    if (iph == S_phase) then
        cnt = cnt + 1
        ATM_DSB = cp%DSB(NHEJslow) + cp%DSB(HR)   ! complex DSB
        ATM_DSB_ave = ATM_DSB_ave + ATM_DSB
        pATM = cp%ATM_act
        pATM_ave = pATM_ave + pATM
!        fATM = max(0.01,1 - k3*pATM/(k4 + pATM))
!        fATM_ave = fATM_ave + fATM
    endif
enddo
ATM_DSB_ave = ATM_DSB_ave/cnt
!fATM_ave = fATM_ave/cnt
pATM_ave = pATM_ave/cnt
write(nflog,'(a,i6,f8.2,i6,3f8.3)') 'S stats: ',istep, hour,cnt, ATM_DSB_ave,pATM_ave

end subroutine
#endif

end module
