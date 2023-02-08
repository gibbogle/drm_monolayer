module mc
use real_kind_mod
use global

use chemokine   ! is this OK?
use eta_module

implicit none

! There are 5 pathways:
! 1 = NHEJs fast pre-replication (simple)
! 2 = NHEJc slow pre-replication (complex)
! 3 = HRc slow post-replication (complex)
! 4 = HRs slow post-replication (simple)
! 5 = TMEJ very slow post-replication (complex)

!! 3 = NHEJ fast post-replication (simple) --> NHEJs
!! 4 = HR slow post-replication (complex)
!! 5 = HR slow post-replication (simple)
!! 6 = MMEJ very slow post-replication (complex)

! There are 5 pathways:
integer, parameter :: NHEJs = 1
integer, parameter :: NHEJc = 2
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
!real(8) :: PHR = 0.8      ! fraction of post-replication simple DSBs that are repaired by HR (slow) rather than NHEJ (fast)
real(8) :: pHRs_S, pHRc_S, pHRs_G2, pHRc_G2
logical :: use_sigmoid = .true.
real(8) :: rmin = 0.1, kdecay = 0.1, ksig = 1, csig = 8.56    ! decay function parameters
character*(2) :: phaseName(4) = ['G1','S ','G2','M ']
!real(8) :: repRate(NP) = [2.081, 0.2604, 2.081, 0.2604, 0.008462]   ! by pathway
!real(8) :: repRate(NP) = [2.081, 0.2604, 2.081, 0.2604, 0.2604, 0.008462]   ! by pathway  with HR simple
!real(8) :: repRate(NP) = [2.081, 0.2604, 0.2604, 0.2604, 0.008462]   ! by pathway  with HR simple
real(8) :: repRate(NP) = [2.081, 0.2604, 0.2604, 0.008462, 0.0]   ! by pathway (TMEJ was 0.008462) (McMahon: fastRepair, slowRepair, verySlowRepair)
!real(8) :: misrepRate(8) = [0.18875,0.18875,0.18247,0.18247,0.14264,0.14264,0.18875,0.18875]        ! by phase, Nlethal = 0.5*misrepRate*Nmis 
!real(8) :: misrepRate(8) = [0.18875,0.0,0.18247,0.0,0.14264,0.0,0.18875,0.0]        ! by phase, Nlethal = 0.5*misrepRate*Nmis 
!real(8) :: fidRate(NP)  = [0.98537, 0.98537, 0.98537, 1.0, 0.4393]  ! by pathway
!real(8) :: fidRate(NP)  = [0.98537, 0.98537, 0.98537, 1.0, 1.0, 0.4393]  ! by pathway  with HR simple
!real(8) :: fidRate(NP)  = [0.98537, 0.98537, 1.0, 1.0, 0.4393]  ! by pathway  with HR simple
real(8) :: fidRate(NP)  = [0.98537, 1.0, 0.4393, 0.0, 0.0]  ! by pathway  with HR simple (McMahon: NHEJFidelity, NHEJFidelity, 1.0, MMEJFidelity)
!logical :: pathwayUsed(8,NP)
real(8) :: apopRate = 0.01117   ! (McMahon: apoptoticRate) (NOT USED only baseRate is used)
real(8) :: baseRate = 0.000739  ! (McMahon: baseRate)
real(8) :: mitRate  = 0.0141    ! (McMahon: mitoticRate)
real(8) :: Msurvival = 0.05
real(8) :: Kaber = 1.0          ! now fixed, McMahon has 1.  
real(8) :: Klethal = 0.4
!real(8) :: K_ATM(4) = [1.0, 0.076, 0.3, 10.6]   ! Now fix K_ATM(1), K_ATR(1) at 1.0
!real(8) :: K_ATR(4) = [1.0, 0.005, 0.3, 3.4]
logical :: use_phase_dependent_CP_parameters
real(8) :: K_ATM(3,4) ! = [0.076, 0.3, 1.0, 1.0]    ! (1) and (2) are the parameters of kinase kinetics, (3) and (4) are CP slowdown parameters
real(8) :: K_ATR(3,4) ! = [0.005, 0.3, 1.0, 1.0]
real(8) :: Ztime = 0    ! hours
logical :: use_exp_slowdown = .true.

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
logical :: use_SSA = .false.
real(8) :: SSAfrac = 0.1
!real(8) :: SSA_rep_fraction = 1.0
!real(8) :: SSA_fid_fraction = 0.5
!logical :: read_SSA_parameters = .false.
!logical :: alt_EJ_suppressed = .false.

! Jaiswal formulation (26/09/22)
real(8) :: Kcc2a, Kcc2e, Kd2e, Kd2t, Ke2cc, Km1, Km2, Km10, Kt2cc, Kti2t, Km10t
real(8) :: CC_tot, ATR_tot, ATM_tot, CC_act0, CC_threshold, norm_factor
logical :: use_Jaiswal = .true.
logical :: vary_km10 = .true.

real(8) :: ATMsum, ATRsum, Sthsum, G2thsum(2)
integer :: NSth, NG2th
real(8) :: repRateFactor(NP)

real(8) :: totPmit, totPaber, tottotDSB, totNlethal

real(8) :: tCPdelay, tATMdelay, tATRdelay
logical :: use_addATMATRtimes = .false.
logical :: use_G1_stop = .true.
logical :: use_S_stop = .false.
logical :: use_G2_stop = .false.
real(8) :: totG1delay, totSdelay, totG2delay
integer :: nG1delay, nSdelay, nG2delay
logical :: use_G2_pATM_Nindependent = .false.
logical :: output_DNA_rate = .false.
logical :: FIX_katm1s_eq_katm1g2 = .false.  ! handle this in the input files, making katm1g2 = katm1s
logical :: negligent_G2_CP = .false.
logical :: use_DSB_CP = .false.
logical :: use_D_model = .true.

! Normalisation
!real(8) :: control_ave(4) = [35.746, 48.898, 13.735, 1.621]
real(8) :: control_ave(4)   ! now set equal to ccp%f_G1, ...
logical :: normalise, M_only

! G1 checkpoint
logical :: use_G1_CP_factor = .false.
real(8) :: G1_CP_factor = 0.0
real(8) :: G1_CP_time

! Checking TMEJ misjoining
real(8) :: misjoins(2)      ! 1 = NHEJ, 2 = TMEJ


!DEC$ ATTRIBUTES DLLEXPORT :: Pcomplex, apopRate, baseRate, mitRate, Msurvival, Kaber, Klethal, K_ATM, K_ATR !, KmaxInhibitRate, b_exp, b_hill

contains

!--------------------------------------------------------------------------
! Note: decision about PEST output is based on iphase_hours. 
!--------------------------------------------------------------------------
subroutine ReadMcParameters(nfin)
integer :: nfin
integer :: iuse_baserate, iuse_exp, iphase_hours, icase, nCPparams, iph, j
real(8) :: TMEJrep, TMEJfid, SSArep, SSAfid
real(8) :: pHR_S, pfc_S, pHR_G2, pfc_G2

write(*,*) 'ReadMcParameters:'
read(nfin,*) iphase_hours
write(*,*) 'iphase_hours: ',iphase_hours
read(nfin,*) baseRate
write(*,*) 'baseRate: ',baseRate
read(nfin,*) mitRate
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
    use_phase_dependent_CP_parameters = .true.      ! now always true
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
!else
!    write(*,*) 'ERROR: wrong nCPparams: ',nCPparams
!    stop
!endif

read(nfin,*) Pcomplex
read(nfin,*) pHRs_S
read(nfin,*) pHRc_S
read(nfin,*) rmin
read(nfin,*) ksig
read(nfin,*) Chalf
read(nfin,*) Preass
read(nfin,*) TMEJrep
read(nfin,*) sigma_TMEJ
read(nfin,*) SSArep
read(nfin,*) SSAfrac
repRate(4) = TMEJrep
!repRate(5) = SSArep
dsigma_dt = SSArep      ! to use to investigate sigma_NHEJ growing with time since IR

if (use_Jaiswal) then
    read(nfin,*) Kcc2a
    read(nfin,*) Kcc2e
    read(nfin,*) Kd2e
    read(nfin,*) Kd2t
    read(nfin,*) Ke2cc
    read(nfin,*) Km1
    read(nfin,*) Km10
    read(nfin,*) Kt2cc
    read(nfin,*) Kti2t
    read(nfin,*) Km10t
    read(nfin,*) CC_tot
    read(nfin,*) ATR_tot
    read(nfin,*) ATM_tot
    read(nfin,*) CC_act0
    read(nfin,*) CC_threshold
    read(nfin,*) norm_factor
    if (CC_threshold < 0.1) then
        use_slope_threshold = .true.
        slope_threshold = CC_threshold
    else
        use_slope_threshold = .false.
        CC_threshold = CC_threshold*CC_tot
    endif
    Km2 = Km1
endif

call make_eta_table(sigma_NHEJ, sigma_TMEJ)

ATMsum = 0  ! to investigate ATM dependence on parameters
ATRsum = 0  ! to investigate ATR dependence on parameters
!Sthsum = 0
NSth = 0
G2thsum = 0
NG2th = 0

! For PEST runs, iphase_hours must be -1 (for M runs) or -2 (for C runs) or -3 (for MC runs)
use_SF = .false.
nphase_hours = 0
next_phase_hour = 0
phase_hour(:) = 0
output_DNA_rate = .false.
normalise = .false.
M_only = .false.
if (iphase_hours == 1) then
    use_SF = .true.     ! in this case SFave only is recorded
    compute_cycle = .false.
!elseif (iphase_hours == -2 .or. iphase_hours == -102 .or. iphase_hours == -112) then    ! this is the compute_cycle case for CA-135
elseif (mod(iphase_hours,10) == 2) then    ! this is the compute_cycle case for CA-135
    compute_cycle = .true.
    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
    nphase_hours = 3
    next_phase_hour = 1
    phase_hour(1:5) = [5.0, 8.5, 11.5, 0.0, 0.0]   ! these are hours post irradiation, incremented when irradiation time is known (in ReadProtocol)
elseif (mod(iphase_hours,10) == 5) then    ! this is the compute_cycle case for CC-11
    CC11 = .true.
    compute_cycle = .true.
    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
    nphase_hours = 5
    next_phase_hour = 1
    phase_hour(1:5) = [1.0, 1.5, 2.5, 3.5, 4.5]   ! these are hours post irradiation, incremented when irradiation time is known (in ReadProtocol)
elseif (iphase_hours == 6) then    ! this is the output_DNA_rate case
    CC11 = .true.
    compute_cycle = .false.
    output_DNA_rate = .true.
    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
! Generating DNA synthesis data
!    nphase_hours = 10  ! To generate DNA synthesis factor from S ATM parameters
!    phase_hour(1:10) = [0.05,0.1,0.16666,0.33333,0.5,0.75,1.25,2.25,3.25,4.25]   ! these are hours post irradiation, incremented when irradiation time is known (in ReadProtocol)
! EDU fitting
    nphase_hours = 5    ! To fit S ATM parameters to EDU data
    phase_hour(1:5) = [0.75,1.25,2.25,3.25,4.25]   ! these are hours post irradiation, incremented when irradiation time is known (in ReadProtocol)
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
    CC11 = .true.
    compute_cycle = .true.
    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
    nphase_hours = 25
    next_phase_hour = 1
    phase_hour(1:25) = [0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0,24.0]   ! these are hours post irradiation, incremented when irradiation time is known (in ReadProtocol)   
elseif (mod(iphase_hours,10) == 8) then    ! this is the compute_cycle case for M%-only experiments
    CC11 = .true.
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

! Test changes to repRate(HR)
!write(nflog,*) '!!!!!!!!!!!!!!!!!!!! changing repRate(HR) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!repRate(HR) = 0.4*repRate(HR)
!write(nflog,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
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
! pHRs_S = pHR_S/2, pHRc_S = pHR_S/2
! pHRs_G2 = pHR_G2/2, pHRc_G2 = pHR_G2/2 with pHR_G2 = 0.15 (fixed)
!
! Option 2: Bill's suggestion, making the pHRs and pHRc parameters differ,
! and preserving the total complex fraction.
!
! Note: pHRs + pHRc must be <= 1
!--------------------------------------------------------------------------
subroutine cellIrradiation(cp, dose, Cdrug)
type(cell_type), pointer :: cp
real(8) :: dose, Cdrug
integer :: ityp, kpar = 0
integer :: phase
real(8) :: DSB0(NP)
real(8) :: totDSB0, baseDSB, fin, T_S, f_S, NG1, NNHEJ, pHRs, pHRc, NHRs, NHRc
real(8) :: Pbase, Pdie, R
real(8) :: DSB_Gy = 35
real(8) :: th, Npre, Npre_s, Npre_c, Npost, Npost_s, Npost_c
integer, parameter :: option = 2

phase = cp%phase
cp%phase0 = min(phase, M_phase)
NG1 = DSB_Gy*dose
DSB0 = 0

#if 0
if (phase == G1_phase) then
    f_S = 0
    pHRs = 0
    pHRc = 0
else
    if (phase >= G2_phase) then
        f_S = 1.0
        pHRs = pHRs_G2
        pHRc = pHRc_G2
    else
        f_S = cp%progress
        pHRs = pHRs_S
        pHRc = pHRc_S
    endif    
endif

totDSB0 = (1+f_S)*NG1
NHRs = f_S*pHRs*2*NG1
NHRc = f_S*pHRc*2*NG1
DSB0(HR) = NHRs + NHRc
if (option == 1) then
    DSB0(NHEJc) = Pcomplex*(totDSB0 - DSB0(HR)) 
    DSB0(NHEJs) = (1 - Pcomplex)*(totDSB0 - DSB0(HR)) 
elseif (option == 2) then
    DSB0(NHEJc) = Pcomplex*totDSB0 - NHRc
    DSB0(NHEJs) = (1 - Pcomplex)*totDSB0 - NHRs
endif
#endif

if (phase == G1_phase) then
    f_S = 0.0
    totDSB0 = NG1
    DSB0(NHEJc) = Pcomplex*totDSB0
    DSB0(NHEJs) = (1 - Pcomplex)*totDSB0
else
    th = (istep*DELTA_T - cp%t_S_phase)/3600  ! time since start of S_phase (h)
    if (phase >= G2_phase) then
        f_S = 1.0
    else
        f_S = cp%progress
    endif 
    totDSB0 = (1+f_S)*NG1
    pHRs = pHRs_S*((1-rmin)*fdecay(th) +rmin)
    pHRc = pHRc_S*((1-rmin)*fdecay(th) +rmin)
    Npre = NG1*(1 - f_S)
    Npost = NG1*2*f_S
    Npre_s = Npre*(1 - pComplex)
    Npre_c = Npre*pComplex
    Npost_s = Npost*(1 - pComplex)
    Npost_c = Npost*pComplex
    NHRs = Npost_s*pHRs
    NHRc = Npost_c*pHRc
    DSB0(HR) = NHRs + NHRc
    DSB0(NHEJs) = Npre_s + Npost_s - NHRs
    DSB0(NHEJc) = Npre_c + Npost_c - NHRc       
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
if (kcell_now <= 10) write(*,'(a,2i6,8f8.2)') 'IR: kcell, phase,DSB,f_S: ',kcell_now,phase,cp%DSB(1:4),f_S

totPmit = 0
totPaber = 0
tottotDSB = 0
totNlethal = 0

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

pATM0 = pATM
if (use_G2_pATM_Nindependent .and. iph == G2_phase) then
    r = K_ATM(iph,1)   ! rate of production of pATM
else
    r = K_ATM(iph,1)*ATM_DSB   ! rate of production of pATM
endif
t = (tnow - t_irradiation)/3600
if (Ztime > 0) then
    if (t > Ztime) then
        fz = 0
    else
        fz = 1 - t/Ztime
    endif
else    ! no reduction in pATM production
    fz = 1
endif
r = fz*r
kdecay = max(1.0e-8,K_ATM(iph,2))  ! decay rate constant
if ((iph == 3) .and. use_D_model) then
    r = 0
    kdecay = max(1.0e-8,K_ATM(iph,3))  ! decay rate constant
endif
pATM = r/kdecay + (pATM - r/kdecay)*exp(-kdecay*dth)
if (isnan(pATM)) then
    write(*,*) 'NAN in updateATM: ',iph, r, kdecay, r/kdecay
    stop
endif
!if (kcell_now == 1) write(*,'(a,i4,6f6.2)') 'updateATM: fz,r,k,r/k,ATM_DSB,pATM: ',kcell_now,fz,r,kdecay,r/kdecay,ATM_DSB,pATM
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
! Effect of pATM
! Time computed in hours, returned in secs 
! CP delay determined by pATM only.
!------------------------------------------------------------------------
function G1_checkpoint_time(cp) result(t)
type(cell_type), pointer :: cp
real(REAL_KIND) :: t, th
integer :: iph = 1

if (use_DSB_CP) then
    t = 0
    return
endif
if (use_G1_CP_factor) then
    t = G1_CP_time
    return
endif
th = K_ATM(iph,3)*(1 - exp(-K_ATM(iph,4)*cp%pATM))
!if (kcell_now == 1) write(*,'(a,i6,2f8.3)') 'G1_checkpoint_time (h): ',kcell_now,cp%pATM,th
if (isnan(th)) then
    write(*,*) 'NAN in G1_checkpoint_time: ',cp%pATM
    stop
endif
totG1delay = th + totG1delay
nG1delay = nG1delay + 1
t = 3600*th
!if (is_radiation) then
!    Sthsum = Sthsum + th
!    NSth = NSth + 1
!endif
end function

!------------------------------------------------------------------------
! Combined effect of pATM and pATR
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
! Now slowdown is acting only in S-phase, and only affected by pATM
!------------------------------------------------------------------------
subroutine get_slowdown_factors(cp,iph,fATM,fATR)
type(cell_type), pointer :: cp
integer :: iph
real(REAL_KIND) :: fATM, fATR
real(REAL_KIND) :: pATM, pATR, k3, k4, N_DSB

fATR = 1
if (iph /= S_phase) then
    fATM = 1
    return
endif
if (use_DSB_CP) then
    N_DSB = sum(cp%DSB)
    k3 = K_ATM(iph,3)
    k4 = K_ATM(iph,4)
    fATR = 1
!    if (kcell_now == 1) write(*,'(a,i6,4f8.3)') 'k3,k4,N_DSB: ',kcell_now,k3,k4,N_DSB,k3*N_DSB/(k4 + N_DSB)
    fATM = max(0.01,1 - k3*N_DSB/(k4 + N_DSB))
    return
endif
pATM = cp%pATM
k3 = K_ATM(iph,3)
k4 = K_ATM(iph,4)
if (use_exp_slowdown) then
    fATM = exp(-k4*pATM)
else
    fATM = max(0.01,1 - k3*pATM/(k4 + pATM))
endif
!write(*,'(a,i3,4f8.4)') 'fATM: ',iph,k3,k4,pATM,fATM
!if (iph > G1_phase) then
!    pATR = cp%pATR
!    k3 = K_ATR(iph,3)   !*G2_katr3_factor
!    k4 = K_ATR(iph,4)   !*G2_katr4_factor
!    fATR = max(0.01,1 - k3*pATR/(k4 + pATR))
!else
!    fATR = 1.0
!endif
!write(*,'(a,i3,4f8.4)') 'fATR: ',iph,k3,k4,pATR,fATR
end subroutine

!------------------------------------------------------------------------
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
        (iph == 3 .and. use_G2_stop)) then
        fslow = 1.0
    else
        dt = DELTA_T
        call get_slowdown_factors(cp,iph,fATM, fATR)
        if (use_addATMATRtimes) then
            dt = DELTA_T
            fslow = max(0.0,fATM + fATR - 1)
        else
            fslow = fATM*fATR
!            if (kcell_now == 1) write(*,'(a,i6,i4,f6.3)') 'fslow: ',kcell_now,iph,fslow
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
subroutine G2_Jaiswal_update(cp, dth)
type(cell_type), pointer :: cp
real(8) :: dth
real(8) :: dt = 0.001
real(8) :: D, CC_act, ATR_act, ATM_act, CC_inact, ATR_inact, ATM_inact
real(8) :: dCC_act_dt, dATR_act_dt, dATM_act_dt, t, T_G2, kkm10
integer :: it, Nt
type(cycle_parameters_type),pointer :: ccp

if (vary_km10 .and. .not. is_radiation) then
    ccp => cc_parameters(1)
    T_G2 = ccp%T_G2*cp%fg(3)/3600
    kkm10 = 5 + 3*(T_G2 - 2.5)
!    write(*,'(a,3f8.3)') 'T_G2, fg(3), kkm10: ',T_G2,cp%fg(3),kkm10
else
    kkm10 = Km10
endif
Nt = int(dth/dt + 0.5)
!if (single_cell) write(nflog,*) 'G2_Jaiswal_update: Nt: ',Nt
D = sum(cp%DSB(1:4))*norm_factor
CC_act = cp%CC_act
ATR_act = cp%ATR_act
ATM_act = cp%ATM_act
do it = 1,Nt
    CC_inact = CC_tot - CC_act
    ATR_inact = ATR_tot - ATR_act
    ATM_inact = ATM_tot - ATM_act
    dCC_act_dt = (Kcc2a + CC_act) * CC_inact / (kkm10 + CC_inact) - Kt2cc * ATM_act * CC_act / (Km10t + CC_act) - Ke2cc * ATR_act * CC_act / (kkm10 + CC_act)
    dATR_act_dt = Kd2e * D * ATR_inact / (kkm10 + ATR_inact) - Kcc2e * ATR_act * CC_act / (kkm10 + CC_act)
    dATM_act_dt = Kd2t * D * ATM_inact / (Km1 + ATM_inact) - Kti2t * ATM_act / (Km1 + ATM_act)
    
    CC_act = CC_act + dt * dCC_act_dt
    ATR_act = ATR_act + dt * dATR_act_dt
    ATM_act = ATM_act + dt * dATM_act_dt
    
    t = it*dt
!    write(nflog,'(i6,f8.4,2f10.6)') it,t,CC_act,dCC_act_dt
enddo
if (kcell_now == 1) write(*,'(a,4f8.4)') 'ATR_act, Kd2e,D,ATR_inact: ',ATR_act, Kd2e,D,ATR_inact
if (kcell_now <= 10) write(*,'(a,i8,4f8.4)') 'CC_act: ',kcell_now,ATR_act,ATM_act,CC_act,dCC_act_dt
cp%CC_act = CC_act
cp%ATR_act = ATR_act
cp%ATM_act = ATM_act
cp%dCC_act_dt = dCC_act_dt
t = t_simulation/3600.
!cp%progress = (cp%CC_act - CC_act0)/(CC_threshold - CC_act0)
if (single_cell) write(*,'(a,f6.2,4e12.3)') 'G2_J: t, vars: ',t,cp%CC_act,cp%ATR_act,cp%ATM_act,D
!if (kcell_now == 3) write(*,'(a,f6.2,4e12.3)') 'G2_J: t, vars: ',t,cp%CC_act,dCC_act_dt
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
    call G2_Jaiswal_update(cp, dth)
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
! 
!------------------------------------------------------------------------
subroutine updateRepair(cp,dt)
type(cell_type), pointer :: cp
real(8) :: dt
integer :: phase
real(8) :: DSB(NP), dNlethal
real(8) :: DSB0(NP)
real(8) :: totDSB0, totDSB, Pmis, Nmis, dNmis(NP), totDSBinfid0, totDSBinfid, ATR_DSB, ATM_DSB, dth, binMisProb
real(8) :: Cdrug, inhibrate, Nreassign
real(8) :: f_S, eta_NHEJ, eta_TMEJ, tIR
logical :: pathUsed(NP)
integer :: k, iph
logical :: use_DSBinfid = .true.
real(8) :: DSB_min = 1.0e-3
logical :: use_totMis = .true.      ! was false!
logical :: use_ATM != .false.
logical :: dbug

dth = dt/3600   ! hours
use_ATM = .not.use_fixed_CP
phase = cp%phase
!if (use_phase_dependent_CP_parameters) then    ! always now
    iph = phase
!else
!    iph = 1
!endif
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
if (Preass > 0) then
    do k = 1,2  ! only NHEJ pathways
        Nreassign = DSB(k)*Preass*dth
        DSB(k) = DSB(k) - Nreassign
        if (use_SSA) then
            DSB(TMEJ) = DSB(TMEJ) + Nreassign*(1 - SSAfrac)
            DSB(SSA) = DSB(SSA) + Nreassign*SSAfrac
        else
            DSB(TMEJ) = DSB(TMEJ) + Nreassign
        endif
    enddo
endif
DSB0 = DSB     ! initial DSBs for this time step

! Revised approach with no fidelity, just misrejoining

totDSB0 = sum(DSB0)
if (totDSB0 == 0) return

!if (kcell_now == 1) write(*,'(a,2i6,3f8.2)') 'updateRepair: kcell, phase, DSB0: ',kcell_now,cp%phase,cp%DSB(1:3)

ATM_DSB = DSB(NHEJc) + DSB(HR)   ! complex DSB
ATR_DSB = DSB(HR)
if (iph >= 7) iph = iph - 6     ! checkpoint phase numbers --> phase number, to continue pATM and pATR processes through checkpoints
if (iph <= G2_phase) then      ! not for 4 (M_phase) or 5 (dividing)
    if (iph == G2_phase .and. use_Jaiswal) then
        call G2_Jaiswal_update(cp, dth)
    else
        call updateATM(iph,cp%pATM,ATM_DSB,dth)     ! updates the pATM mass through time step = dth
        if (iph > G1_phase) call updateATR(iph,cp%pATR,ATR_DSB,dth)     ! updates the pATR mass through time step = dth
    endif
endif

DSB = 0
do k = 1,NP
!   if (dbug .and. DSB0(k) > 0) write(*,*) 'pathwayRepair: k,DSB0(k): ',kcell_now,k,DSB0(k)
    call pathwayRepair(k, dth, DSB0(k), DSB(k))
    if (DSB(k) < DSB_min) DSB(k) = 0
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
eta_NHEJ = eta_lookup(phase, NHEJs, f_S, tIR) 
eta_TMEJ = eta_lookup(phase, TMEJ, f_S, tIR) 

Nmis = 0
! For NHEJ pathways
totDSB0 = DSB0(NHEJs) + DSB0(NHEJc)
totDSB = DSB(NHEJs) + DSB(NHEJc)
Pmis = misrepairRate(totDSB0, totDSB, eta_NHEJ)
Nmis = Nmis + Pmis*(totDSB0 - totDSB)
misjoins(1) = misjoins(1) + Pmis*(totDSB0 - totDSB)
! For TMEJ pathway
Pmis = misrepairRate(DSB0(TMEJ), DSB(TMEJ), eta_TMEJ)
Nmis = Nmis + Pmis*(DSB0(TMEJ) - DSB(TMEJ))
misjoins(2) = misjoins(2) + Pmis*(totDSB0 - totDSB)
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
real(8) :: DSB(NP), totDSB, Nlethal,Paber, Pbase, Papop, Pmit, Psurv

DSB = cp%DSB
totDSB = sum(DSB)
Nlethal = cp%Nlethal
if (kcell_now == 1) write(nflog,'(a,3i4,2e12.3)') 'kcell, phase0, phase, totDSB, Nlethal: ',kcell_now, cp%phase0, cp%phase, totDSB, Nlethal
if (cp%phase0 == G1_phase) then     ! same as S, G2
!    Pbase = exp(-baseRate*cp%totDSB0)   ! this is handled in cellIrradiation() = 1 - Pdie
    Paber = exp(-Kaber*Nlethal) ! Kaber = 1, not needed because Klethal is a parameter
    Pmit = exp(-mitRate*totDSB)
    cp%Psurvive = Paber*Pmit
!    if (kcell_now < 40) write(*,'(a,i3,2f5.1,3e11.3)') 'G1: Nlethal,totDSB,Paber,Pmit: ',kcell_now,Nlethal,totDSB,Paber,Pmit,Paber*Pmit
!       cp%Psurvive = Paber*Pbase*Papop*Pmit
!    if (kcell_now == 35) write(*,*) 'cell #35 Psurvive: ',cp%Psurvive
elseif (cp%phase0 < M_phase) then   ! G1, S, G2
    Paber = exp(-Nlethal)
    Pmit = exp(-mitRate*totDSB)
    cp%Psurvive = Pmit*Paber  
    if (kcell_now == 1) write(*,'(a,2f8.4,2e12.3)') 'totDSB,Nlethal,Pmit,Paber: ',totDSB,Nlethal,Pmit,Paber  
else    ! M_phase or dividing
    Paber = 1
    Pmit = exp(-mitRate*totDSB)
    cp%Psurvive = Pmit*Msurvival
endif
NPsurvive = NPsurvive + 1   ! this is the count of cells for which Psurvive has been computed
cp%mitosis_time = tnow      ! needed to determine if mitosis occurs before or after CA
if (istep > 1278) write(*,*) 'reached M: phase0: ',istep,cp%phase0
!Psurvive(NPsurvive) = cp%Psurvive

ATMsum = ATMsum + cp%pATM
ATRsum = ATRsum + cp%pATR

totPmit = totPmit + Pmit
totPaber = totPaber + Paber
tottotDSB = tottotDSB + totDSB
totNlethal = totNlethal + Nlethal

end subroutine

!------------------------------------------------------------------------
! Get average DNA growth factor for cells in S-phase
!------------------------------------------------------------------------
subroutine get_DNA_synthesis_rate(DNA_rate)
real(8) :: DNA_rate
type(cell_type), pointer :: cp
integer :: kcell, iph, cnt
real(8) :: pATM, k1, k2, fATM, rate_sum, pATM_ave

!write(*,'(a)') 'get_DNA_synthesis_rate'
k1 = K_ATM(S_phase,3)
k2 = K_ATM(S_phase,4)
cnt = 0
rate_sum = 0
pATM_ave = 0
do kcell = 1,nlist
    cp => cell_list(kcell)
    iph = cp%phase
    if (iph == S_phase) then
        cnt = cnt + 1
        pATM = cp%pATM
        pATM_ave = pATM_ave + pATM
        fATM = max(0.01,1 - k1*pATM/(k2 + pATM))
        rate_sum = rate_sum + fATM
    endif
enddo
DNA_rate = rate_sum/cnt
pATM_ave = pATM_ave/cnt
!write(*,'(a,3f8.3,e12.3)') 'DNA growth rate factor: ',DNA_rate,k1,k2,pATM_ave
end subroutine

!------------------------------------------------------------------------
! Write average S-phase ATM_DSB, pATM, fATM
!------------------------------------------------------------------------
subroutine show_S_phase_statistics()
real(8) :: hour
type(cell_type), pointer :: cp
integer :: kcell, iph, cnt, nthour
real(8) :: pATM, k1, k2, fATM, fATM_ave, pATM_ave, ATM_DSB, ATM_DSB_ave

!write(nflog,'(a)') 'S_phase_statistics'
nthour = 3600/DELTA_T
hour = real(istep)/nthour
k1 = K_ATM(S_phase,3)
k2 = K_ATM(S_phase,4)
cnt = 0
ATM_DSB_ave = 0
fATM_ave = 0
pATM_ave = 0
do kcell = 1,nlist
    cp => cell_list(kcell)
    iph = cp%phase
    if (iph == S_phase) then
        cnt = cnt + 1
        ATM_DSB = cp%DSB(NHEJc) + cp%DSB(HR)   ! complex DSB
        ATM_DSB_ave = ATM_DSB_ave + ATM_DSB
        pATM = cp%pATM
        pATM_ave = pATM_ave + pATM
        fATM = max(0.01,1 - k1*pATM/(k2 + pATM))
        fATM_ave = fATM_ave + fATM
    endif
enddo
ATM_DSB_ave = ATM_DSB_ave/cnt
fATM_ave = fATM_ave/cnt
pATM_ave = pATM_ave/cnt
write(nflog,'(i6,f8.2,i6,3f8.3)') istep, hour,cnt, ATM_DSB_ave,pATM_ave,fATM_ave

end subroutine

end module
