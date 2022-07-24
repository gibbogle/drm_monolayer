module mc
use real_kind_mod
use global

implicit none

! There are 6 pathways:
! 1 = NHEJs fast pre-replication (simple)
! 2 = NHEJc slow pre-replication (complex)
! 3 = HRc slow post-replication (complex)
! 4 = HRs slow post-replication (simple)
! 5 = MDR very slow post-replication (complex)

!! 3 = NHEJ fast post-replication (simple) --> NHEJs
!! 4 = HR slow post-replication (complex)
!! 5 = HR slow post-replication (simple)
!! 6 = MMEJ very slow post-replication (complex)

! There are 4 pathways:
integer, parameter :: NHEJs = 1
integer, parameter :: NHEJc = 2
integer, parameter :: HR = 3
integer, parameter :: MDR = 4   ! this is really alt-EJ
integer, parameter :: SSA = 5
!integer, parameter :: HRc = 3
!integer, parameter :: HRs = 4
!integer, parameter :: MDR = 5

real(8) :: eta = 0.0006506387875449218
real(8) :: Pcomplex = 0.4337    ! fraction of created DSBs that are complex
real(8) :: PHRsimple = 0.8      ! fraction of post-replication simple DSBs that are repaired by HR (slow) rather than NHEJ (fast)
character*(2) :: phaseName(4) = ['G1','S ','G2','M ']
!real(8) :: repRate(NP) = [2.081, 0.2604, 2.081, 0.2604, 0.008462]   ! by pathway
!real(8) :: repRate(NP) = [2.081, 0.2604, 2.081, 0.2604, 0.2604, 0.008462]   ! by pathway  with HR simple
!real(8) :: repRate(NP) = [2.081, 0.2604, 0.2604, 0.2604, 0.008462]   ! by pathway  with HR simple
real(8) :: repRate(NP) = [2.081, 0.2604, 0.2604, 0.008462, 0.0]   ! by pathway (MDR was 0.008462)
!real(8) :: misrepRate(8) = [0.18875,0.18875,0.18247,0.18247,0.14264,0.14264,0.18875,0.18875]        ! by phase, Nlethal = 0.5*misrepRate*Nmis 
!real(8) :: misrepRate(8) = [0.18875,0.0,0.18247,0.0,0.14264,0.0,0.18875,0.0]        ! by phase, Nlethal = 0.5*misrepRate*Nmis 
!real(8) :: fidRate(NP)  = [0.98537, 0.98537, 0.98537, 1.0, 0.4393]  ! by pathway
!real(8) :: fidRate(NP)  = [0.98537, 0.98537, 0.98537, 1.0, 1.0, 0.4393]  ! by pathway  with HR simple
!real(8) :: fidRate(NP)  = [0.98537, 0.98537, 1.0, 1.0, 0.4393]  ! by pathway  with HR simple
real(8) :: fidRate(NP)  = [0.98537, 0.98537, 1.0, 0.4393, 0.0]  ! by pathway  with HR simple ! MDR was 0.4393
!logical :: pathwayUsed(8,NP)
real(8) :: apopRate = 0.01117   ! NOT USED
real(8) :: baseRate = 0.000739
real(8) :: mitRate  = 0.0141
real(8) :: Msurvival = 0.05
real(8) :: Kaber = 1.0          ! now fixed, McMahon has 1.  
real(8) :: Klethal = 0.4
real(8) :: MDRrep
real(8) :: MDRfid
!real(8) :: K_ATM(4) = [1.0, 0.076, 0.3, 10.6]   ! Now fix K_ATM(1), K_ATR(1) at 1.0
!real(8) :: K_ATR(4) = [1.0, 0.005, 0.3, 3.4]
logical :: use_phase_dependent_CP_parameters
real(8) :: K_ATM(3,4) ! = [0.076, 0.3, 1.0, 1.0]    ! (1) and (2) are the parameters of kinase kinetics, (3) and (4) are CP slowdown parameters
real(8) :: K_ATR(3,4) ! = [0.005, 0.3, 1.0, 1.0]

! DNA-PK inhibition parameters
real(8) :: Chalf    ! inhibitor concentration that halves repair rate 
real(8) :: Preass   ! rate of reassignment to pathway 4 (prob of reass/hour)
!real(8) :: KmaxInhibitRate = 0.8
!real(8) :: b_exp = 1.0
!real(8) :: b_hill = 0.5
!logical :: use_exp_inhibit = .true.

! SSA parameters
logical :: use_SSA = .true.
real(8) :: SSA_fraction = 0.1
real(8) :: SSA_rep_fraction = 1.0
real(8) :: SSA_fid_fraction = 0.5
logical :: read_SSA_parameters = .false.
logical :: alt_EJ_suppressed = .false.

real(8) :: ATMsum, ATRsum, Sthsum, G2thsum
integer :: NSth, NG2th
real(8) :: repRateFactor(NP)

real(8) :: totPmit, totPaber, tottotDSB, totNlethal

real(8) :: tCPdelay, tATMdelay, tATRdelay
logical :: use_addATMATRtimes = .false.
logical :: use_stops = .true.
logical :: use_G1_stop = .true.
logical :: use_S_stop = .false.
logical :: use_G2_stop = .true.
real(8) :: totG1delay, totSdelay, totG2delay
integer :: nG1delay, nSdelay, nG2delay
logical :: use_G2_pATM_Nindependent = .true.
logical :: output_DNA_rate = .false.

!DEC$ ATTRIBUTES DLLEXPORT :: Pcomplex, PHRsimple, apopRate, baseRate, mitRate, Msurvival, Kaber, Klethal, K_ATM, K_ATR !, KmaxInhibitRate, b_exp, b_hill

contains

!--------------------------------------------------------------------------
! Note: decision about PEST output is based on iphase_hours.
!--------------------------------------------------------------------------
subroutine ReadMcParameters(nfin)
integer :: nfin
integer :: iuse_baserate, iuse_exp, iphase_hours, nCPparams, iph, j
write(*,*) 'ReadMcParameters:'
read(nfin,*) iphase_hours
write(*,*) 'iphase_hours: ',iphase_hours
read(nfin,*) baseRate
write(*,*) 'baseRate: ',baseRate
read(nfin,*) mitRate
read(nfin,*) Msurvival
read(nfin,*) Klethal
read(nfin,*) nCPparams
if (nCPparams == 1) then
    iph = 1
    use_phase_dependent_CP_parameters = .false.
    read(nfin,*) K_ATM(iph,1)
    read(nfin,*) K_ATM(iph,2)
    read(nfin,*) K_ATM(iph,3)
    read(nfin,*) K_ATM(iph,4)
    read(nfin,*) K_ATR(iph,1)
    read(nfin,*) K_ATR(iph,2)
    read(nfin,*) K_ATR(iph,3)
    read(nfin,*) K_ATR(iph,4)
    do iph = 2,3        ! Note: doing this actually makes it unnecessary to set iph = 1 when .not. use_phase_dependent_CP_parameters
        do j = 1,4
            K_ATM(iph,j) = K_ATM(1,j)
            K_ATR(iph,j) = K_ATR(1,j)
        enddo
    enddo
elseif (nCPparams == 3) then
    use_phase_dependent_CP_parameters = .true.
    do j = 1,4
        do iph = 1,3
            read(nfin,*) K_ATM(iph,j)
            write(*,*) j,iph,K_ATM(iph,j)
        enddo
     enddo
    do j = 1,4
        do iph = 2,3
            read(nfin,*) K_ATR(iph,j)
            write(*,*) j,iph,K_ATR(iph,j)
        enddo
    enddo
else
    write(*,*) 'ERROR: wrong nCPparams: ',nCPparams
    stop
endif
read(nfin,*) Pcomplex
read(nfin,*) PHRsimple
read(nfin,*) Chalf
read(nfin,*) Preass
read(nfin,*) MDRrep
read(nfin,*) MDRfid
if (use_SSA .and. read_SSA_parameters) then
    read(nfin,*) SSA_fraction
    read(nfin,*) SSA_rep_fraction
    read(nfin,*) SSA_fid_fraction
endif
repRate(4) = MDRrep
fidRate(4) = MDRfid
!read(nfin,*) KmaxInhibitRate
!read(nfin,*) b_exp
!read(nfin,*) b_hill
!read(nfin,*) iuse_exp
!use_exp_inhibit = (iuse_exp == 1)
!misrepRate(1) = Kaber*misrepRate(1)

if (use_SSA) then
    repRate(SSA) = SSA_rep_fraction*repRate(MDR)
    fidRate(SSA) = SSA_fid_fraction*fidRate(MDR)
    if (alt_EJ_suppressed) then
        SSA_fraction = 1.0
    endif
endif

ATMsum = 0  ! to investigate ATM dependence on parameters
ATRsum = 0  ! to investigate ATR dependence on parameters
!Sthsum = 0
NSth = 0
!G2thsum = 0
NG2th = 0

! For PEST runs, iphase_hours must be -1 (for M runs) or -2 (for C runs) or -3 (for MC runs)
use_SF = .false.
nphase_hours = 0
next_phase_hour = 0
phase_hour(:) = 0
output_DNA_rate = .false.
if (iphase_hours == -1) then
    use_SF = .true.     ! in this case SFave only is recorded
    compute_cycle = .false.
elseif (iphase_hours == -2) then    ! this is the compute_cycle case for CA-135
    compute_cycle = .true.
    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
    nphase_hours = 3
    next_phase_hour = 1
    phase_hour(1:5) = [4.5, 8.0, 11.0, 0.0, 0.0]   ! these are hours post irradiation, incremented when irradiation time is known (in ReadProtocol)
elseif (iphase_hours == -5) then    ! this is the compute_cycle case for CC-11
    CC11 = .true.
    compute_cycle = .true.
    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
    nphase_hours = 5
    next_phase_hour = 1
    phase_hour(1:5) = [0.5, 1.0, 2.0, 3.0, 4.0]   ! these are hours post irradiation, incremented when irradiation time is known (in ReadProtocol)
elseif (iphase_hours == -6) then    ! this is the output_DNA_rate case
    CC11 = .true.
    compute_cycle = .false.
    output_DNA_rate = .true.
    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
    nphase_hours = 5
    next_phase_hour = 1
    phase_hour(1:5) = [0.5, 1.0, 2.0, 3.0, 4.0]   ! these are hours post irradiation, incremented when irradiation time is known (in ReadProtocol)
elseif (iphase_hours == -3) then
    use_SF = .true.     ! in this case SFave is recorded and there are multiple phase distribution recording times
    nphase_hours = 4
    next_phase_hour = 1
    phase_hour(1:4) = [8, 12, 18, 24]
elseif (iphase_hours == -4) then    ! this is the synchronised IR case
    compute_cycle = .true.
    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
    nphase_hours = 1
    next_phase_hour = 1
    phase_hour(1:5) = [40, 0, 0, 0, 0]   ! these are hours post irradiation, incremented when irradiation time is known (in ReadProtocol)
else
    if (use_PEST) then
        write(*,*) 'Error: ReadMcParameters: with PEST iphase_hours must be -1,-2,-3'
        write(nflog,*) 'Error: ReadMcParameters: with PEST iphase_hours must be -1,-2,-3'
        stop
    endif
endif
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
! Bill 
! So the number assigned to HR, NHR = PHR (1 + fS) NG1. ….(1)
! Similarly, NNHEJ = (1-PHR)(1 + fS)NG1.
! But NNHEJ = NNHEJ,S + NNHEJ,C
! And Pcomplex = NNHEJ, C /NNHEJ
! Thus Pcomplex = NNHEJ,C /(1-PHR)(1 + fS)NG1.
! So NNHEJ,C = Pcomplex (1-PHR)(1 + fS)NG1 ….(2)
! But NNHEJ,S = NNHEJ – NNHEJ,C
! Therefore NNHEJ,S = (1-PHR)(1 + fS)NG1 - Pcomplex (1-PHR)(1 + fS)NG1
!                   = (1-Pcomplex)(1-PHR)(1 + fS)NG1 ….(3)
!--------------------------------------------------------------------------
subroutine cellIrradiation(cp, dose, Cdrug)
type(cell_type), pointer :: cp
real(8) :: dose, Cdrug
integer :: ityp, kpar = 0
integer :: phase
real(8) :: DSB0(NP)
real(8) :: totDSB0, baseDSB, fin, T_S, f_S, PHR, NG1, NNHEJ
real(8) :: Pbase, Pdie, R
real(8) :: DSB_Gy = 35

phase = cp%phase
cp%phase0 = phase
NG1 = DSB_Gy*dose
DSB0 = 0
if (phase == G1_phase) then
    f_S = 0
    PHR = 0
elseif (phase == G2_phase) then
    f_S = 1.0
    PHR = PHRsimple
elseif (phase == S_phase) then
    f_S = cp%progress
    PHR = PHRsimple
endif

totDSB0 = (1+f_S)*NG1
NNHEJ = NG1*(1-f_S) + 2*NG1*f_S*(1-PHR)
DSB0(HR) = totDSB0 - NNHEJ
DSB0(NHEJc) = Pcomplex*NNHEJ
DSB0(NHEJs) = (1 - Pcomplex)*NNHEJ

!if (kcell_now == 76) write(*,'(a,i4,3f8.4)') 'cellIrradiation: phase, f_S, totDSB0, NNHEJ: ',phase, f_S, totDSB0, NNHEJ

!DSB0(HR) = PHR*(1 + f_S)*NG1
!DSB0(NHEJc) = Pcomplex*(1-PHR)*(1 + f_S)*NG1
!DSB0(NHEJs) = (1-Pcomplex)*(1-PHR)*(1 + f_S)*NG1

!baseDSB = DSB_Gy*dose

if (phase == G1_phase) then
!    totDSB0 = baseDSB
!    DSB0(NHEJs) = (1 - Pcomplex)*totDSB0
!    DSB0(NHEJc) = Pcomplex*totDSB0
    ! Apoptosis in G1
    Pbase = exp(-baseRate*totDSB0)   ! this is 1 - Pdie
    Pdie = 1 - Pbase 
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
endif
!if (kcell_now == 1) then
!if (phase == G1_phase) then
!    write(*,'(a,i4,f8.3,5f8.2)') 'IR: phase,f_S,DSB0: ',phase, f_S, DSB0(1:3), totDSB0,sum(DSB0(1:2))
!    stop
!endif
!elseif (phase >= G2_phase) then
!    totDSB0 = 2*baseDSB
!    DSB0(NHEJs) = (1 - Pcomplex)*(1 - PHRsimple)*totDSB0
!    DSB0(HRs) = ((Pcomplex + PHRsimple) - Pcomplex*PHRsimple)*totDSB0
!
!else    ! S_phase
!    if (kcell_now == 1) write(nflog,*) 'in S_phase'
!    DSB0(NHEJs) = (1 - Pcomplex)*((1-f_S) + 2*f_S*(1 - PHRsimple))*baseDSB
!    DSB0(NHEJc) = Pcomplex*(1-f_S)*baseDSB
!    DSB0(HRc) = Pcomplex*2*f_S*baseDSB
!    DSB0(HRs) = (1 - Pcomplex)*2*f_S*PHRsimple*baseDSB  ! This pathway is not used if PHRsimple = 0
!    totDSB0 = (1+f_S)*baseDSB
!endif

!fin = finhibit(Cdrug)
!DSB0(6) = fin*sum(DSB0(1:3))
!DSB0(1:3) = (1 - fin)*DSB0(1:3)
cp%DSB = DSB0
cp%totDSB0 = totDSB0
cp%Nlethal = 0
if (kcell_now <= 10) write(*,'(a,2i6,8f8.2)') 'IR: kcell, phase,DSB,f_S: ',kcell_now,phase,cp%DSB(1:3),f_S

totPmit = 0
totPaber = 0
tottotDSB = 0
totNlethal = 0

end subroutine

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
! Note: parameters from mcradio assume time is hours, therefore the time
! step passed to the subroutine, dth, has been converted to hours.
! 11/02/22 changed K_ATM(2) & (3) to (1) & (2)
!--------------------------------------------------------------------------
subroutine updateATM(iph,pATM,ATM_DSB,dth)
integer :: iph
real(8) :: pATM, ATM_DSB, dth, pATM0
real(8) :: r, kdecay

pATM0 = pATM
if (use_G2_pATM_Nindependent .and. iph == G2_phase) then
    r = K_ATM(iph,1)   ! rate of production of pATM
else
    r = K_ATM(iph,1)*ATM_DSB   ! rate of production of pATM
endif
kdecay = K_ATM(iph,2)  ! decay rate constant
pATM = r/kdecay + (pATM - r/kdecay)*exp(-kdecay*dth)
if (isnan(pATM)) then
    write(*,*) 'NAN in updateATM: ',iph, r, kdecay, r/kdecay
    stop
endif
!if (kcell_now == 1) write(*,'(a,i4,5f8.2)') 'updateATM: r,k,r/k,ATM_DSB,pATM: ',kcell_now,r,k,r/k,ATM_DSB,pATM
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
real(8) :: r, k, x, xmax, km, xf

r = K_ATR(iph,1)*ATR_DSB   ! rate of production of pATR
k = K_ATR(iph,2)  ! decay rate constant
if (k == 0) then
    write(*,*) 'updateATR: k = 0: iph,K_ATR(iph,2): ',iph,K_ATR(iph,2)
    stop
endif
pATR = r/k + (pATR - r/k)*exp(-k*dth)
!if (kcell_now == 76) write(*,'(a,i4,5f8.4)') 'updateATR: r,k,r/k,ATR_DSB,pATR: ',kcell_now,r,k,r/k,ATR_DSB,pATR
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
real(REAL_KIND) :: t, th, th_ATM, th_ATR
integer :: iph = 3

th_ATM = K_ATM(iph,3)*(1 - exp(-K_ATM(iph,4)*cp%pATM))
th_ATR = K_ATR(iph,3)*(1 - exp(-K_ATR(iph,4)*cp%pATR))
th = th_ATM + th_ATR
totG2delay = th + totG2delay
nG2delay = nG2delay + 1
!if (kcell_now == 3) write(*,'(a,i6,2f8.3)') 'G2_checkpoint_time (ATM, ATR) (h): ',kcell_now,th_ATM,th_ATR 
t = 3600*th
!if (is_radiation) then
!    G2thsum = G2thsum + th
!    NG2th = NG2th + 1
!endif
end function

!------------------------------------------------------------------------
! Only G2 is affected by pATR
! Best to turn off pATR in S-phase by setting katr1s = katr3s = 0
!------------------------------------------------------------------------
subroutine get_slowdown_factors(cp,iph,fATM,fATR)
type(cell_type), pointer :: cp
integer :: iph
real(REAL_KIND) :: fATM, fATR
real(REAL_KIND) :: pATM, pATR, k3, k4

pATM = cp%pATM
k3 = K_ATM(iph,3)
k4 = K_ATM(iph,4)
fATM = max(0.01,1 - k3*pATM/(k4 + pATM))
!write(*,'(a,i3,4f8.4)') 'fATM: ',iph,k1,k2,pATM,fATM
if (iph > G1_phase) then
    pATR = cp%pATR
    k3 = K_ATR(iph,3)   !*G2_katr3_factor
    k4 = K_ATR(iph,4)   !*G2_katr4_factor
    fATR = max(0.01,1 - k3*pATR/(k4 + pATR))
else
    fATR = 1.0
endif
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
!    write(*,'(a,i6,f8.3)') 'G1 CP_delay: ',kcell_now,cp%CP_delay/3600
elseif (cp%phase == G2_phase) then
    cp%CP_delay = G2_checkpoint_time(cp)
endif
end subroutine

!--------------------------------------------------------------------------
! To determine repair on a given pathway
! Using parameters from mcradio, dth in hours
!--------------------------------------------------------------------------
subroutine pathwayRepair(path, dth, N0, N)
integer :: path
real(8) :: dth, N0, N
!real(8) :: Kreduction = 1.0

if (dth >= 0) then
    N = N0*exp(-repRate(path)*dth*repRateFactor(path))  !!!! TESTING reducing repair rate
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
logical :: pathUsed(NP)
integer :: k, iph
logical :: use_DSBinfid = .true.
real(8) :: DSB_min = 1.0e-3
real(8) :: eta_G1 = 0.0006506
real(8) :: eta_G2 = 0.0003326
logical :: use_totMis = .false.
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
repRateFactor = 1.0
! DNA-PK inhibition, DSB reassignment
Cdrug = cp%Cin(DRUG_A)
!inhibrate = inhibitRate(Cdrug)
!if (inhibrate > 0) then
!    do k = 1,NP-1
!        Nreassign = DSB(k)*inhibrate*dt
!        DSB(k) = DSB(k) - Nreassign
!        DSB(MDR) = DSB(MDR) + Nreassign
!    enddo
!endif
! Now the effect of inhibition is to reduce the repair rate.
! Chalf is the drug concentration that reduces repair rate by 1/2.  fRR = exp(-k*C/Chalf)
! fRR = 0.5 when C = Chalf, 0.5 = exp(-k), k = -log(0.5) = 0.693
repRateFactor(1:2) = exp(-0.693*Cdrug/Chalf) 
!if (kcell_now == 1) write(*,'(a,3f10.4)') 'Cdrug,Chalf,fRR: ',Cdrug,Chalf,repRateFactor(1)
! Reassignment to pathway 4 is (tentatively) constant
! Preass is an input parameter = prob of reassignment per hour
if (Preass > 0) then
    do k = 1,2  ! only NHEJ pathways
        Nreassign = DSB(k)*Preass*dth
        DSB(k) = DSB(k) - Nreassign
        if (use_SSA) then
            DSB(MDR) = DSB(MDR) + Nreassign*(1 - SSA_fraction)
            DSB(SSA) = DSB(SSA) + Nreassign*SSA_fraction
        else
            DSB(MDR) = DSB(MDR) + Nreassign
        endif
    enddo
endif

DSB0 = DSB     ! initial DSBs for this time step
totDSB0 = sum(DSB0)
if (totDSB0 == 0) return

!if (kcell_now == 1) write(*,'(a,2i6,3f8.2)') 'updateRepair: kcell, phase, DSB0: ',kcell_now,cp%phase,cp%DSB(1:3)

!if (use_ATM) then
!ATM_DSB = DSB(2) + DSB(4)   ! complex DSB
!ATR_DSB = DSB(4) + DSB(5)
ATM_DSB = DSB(NHEJc) + DSB(HR)   ! complex DSB
ATR_DSB = DSB(HR)
!if (phase == G1_phase) then
!    call updateATR(cp%pATR,ATR_DSB,dth)     ! updates the pATR mass through time step = dth
!elseif (phase == S_phase) then
!    call updateATM(cp%pATM,ATM_DSB,dth)     ! updates the pATM mass through time step = dth
!    call updateATR(cp%pATR,ATR_DSB,dth)     ! updates the pATR mass through time step = dth
!else
    if (iph >= 7) iph = iph - 6     ! checkpoint phase numbers --> phase number, to continue pATM and pATR processes through checkpoints
    if (iph <= 3) then      ! not for 4 (M_phase) or 5 (dividing)
        call updateATM(iph,cp%pATM,ATM_DSB,dth)     ! updates the pATM mass through time step = dth
        if (iph > 1) call updateATR(iph,cp%pATR,ATR_DSB,dth)     ! updates the pATR mass through time step = dth
    endif
!endif

!pathUsed = pathwayUsed(phase,:)
dbug = (kcell_now == -3) .and. (cp%phase0 == G2_phase)
DSB = 0
do k = 1,NP
!    if (pathUsed(k)) then
        if (dbug .and. DSB0(k) > 0) write(*,*) 'pathwayRepair: k,DSB0(k): ',kcell_now,k,DSB0(k)
        call pathwayRepair(k, dth, DSB0(k), DSB(k))
!    endif
    if (DSB(k) < DSB_min) DSB(k) = 0
enddo
totDSB = sum(DSB)
totDSBinfid0 = 0
totDSBinfid = 0
do k = 1,NP
    if (DSB0(k) > 0 .and. fidRate(k) < 1.0) then
        totDSBinfid0 = totDSBinfid0 + DSB0(k)
        totDSBinfid = totDSBinfid + DSB(k)
    endif
enddo
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
!if (Pmis > 0.2) write(*,'(a,i6,e12.3)') 'Pmis: ',kcell_now,Pmis
cp%totMis = cp%totMis + Pmis*(totDSB0 - totDSB)
binMisProb = cp%totMis/(cp%totDSB0 - totDSB)
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
!write(nfphase,'(a,5f8.3)') 'dNmis: ',dNmis,Nmis
if (isnan(Nmis)) then
    write(*,*) 'updateRepair: Nmis isnan'
    stop
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
real(8) :: DSB(NP), totDSB, Nlethal,Paber, Pbase, Papop, Pmit, Psurv

DSB = cp%DSB
totDSB = sum(DSB)
Nlethal = cp%Nlethal
if (kcell_now == 1) write(nflog,'(a,3i4,2e12.3)') 'kcell, phase0, phase, totDSB, Nlethal: ',kcell_now, cp%phase0, cp%phase, totDSB, Nlethal
if (cp%phase0 == G1_phase) then
!    Pbase = exp(-baseRate*cp%totDSB0)   ! this is handled in cellIrradiation() = 1 - Pdie
    Paber = exp(-Kaber*Nlethal) ! Kaber = 1, not needed because Klethal is a parameter
    Pmit = exp(-mitRate*totDSB)
    cp%Psurvive = Paber*Pmit
!    if (kcell_now < 40) write(*,'(a,i3,2f5.1,3e11.3)') 'G1: Nlethal,totDSB,Paber,Pmit: ',kcell_now,Nlethal,totDSB,Paber,Pmit,Paber*Pmit
!       cp%Psurvive = Paber*Pbase*Papop*Pmit
!    if (kcell_now == 35) write(*,*) 'cell #35 Psurvive: ',cp%Psurvive
elseif (cp%phase0 < M_phase) then   ! S, G2
    Paber = exp(-Kaber*Nlethal)
    Pmit = exp(-mitRate*totDSB)
    cp%Psurvive = Pmit*Paber  
    if (kcell_now == 1) write(*,'(a,2f8.4,2e12.3)') 'totDSB,Nlethal,Pmit,Paber: ',totDSB,Nlethal,Pmit,Paber   
else    ! M_phase or dividing
    Paber = 1
    Pmit = exp(-mitRate*totDSB)
    cp%Psurvive = Pmit*Msurvival
endif
NPsurvive = NPsurvive + 1   ! this is the count of cells for which Psurvive has been computed
Psurvive(NPsurvive) = cp%Psurvive

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

write(*,'(a)') 'get_DNA_synthesis_rate'
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
write(*,'(a,3f8.3,e12.3)') 'DNA growth rate factor: ',DNA_rate,k1,k2,pATM_ave
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
