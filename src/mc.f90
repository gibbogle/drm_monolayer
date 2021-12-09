module mc
use real_kind_mod
use global

implicit none

! There are 6 pathways:
! 1 = NHEJ fast pre-replication (simple)
! 2 = NHEJ slow pre-replication (complex)
! 3 = NHEJ fast post-replication (simple)
! 4 = HR slow post-replication (complex)
! 5 = HR slow post-replication (simple)
! 6 = MMEJ very slow post-replication (complex)

real(8) :: eta = 0.0006506387875449218
real(8) :: Pcomplex = 0.4337    ! fraction of created DSBs that are complex
real(8) :: PHRsimple = 0.9      ! fraction of post-replication simple DSBs that are repaired by HR (slow) rather than NHEJ (fast)
!character*(2) :: phaseName(8) = ['G1','','S ','','G2','','M ','']
!real(8) :: repRate(NP) = [2.081, 0.2604, 2.081, 0.2604, 0.008462]   ! by pathway
real(8) :: repRate(NP) = [2.081, 0.2604, 2.081, 0.2604, 0.2604, 0.008462]   ! by pathway  with HR simple
real(8) :: misrepRate(8) = [0.18875,0.0,0.18247,0.0,0.14264,0.0,0.18875,0.0]        ! by phase, Nlethal = 0.5*misrepRate*Nmis 
!real(8) :: fidRate(NP)  = [0.98537, 0.98537, 0.98537, 1.0, 0.4393]  ! by pathway
real(8) :: fidRate(NP)  = [0.98537, 0.98537, 0.98537, 1.0, 1.0, 0.4393]  ! by pathway  with HR simple
logical :: pathwayUsed(8,NP)
real(8) :: apopRate = 0.01117
real(8) :: baseRate = 0.000739
real(8) :: mitRate  = 0.0141
real(8) :: Msurvival = 0.5
real(8) :: Kaber = 1.0          ! added, McMahon has 1 
real(8) :: Klethal = 1.65
real(8) :: K_ATM(4) = [1.0, 0.01, 1.0, 10.0]
real(8) :: K_ATR(4) = [1.0, 0.01, 0.1, 10.0]

! DNA-PK inhibition parameters
real(8) :: KmaxInhibit = 0.8
real(8) :: b_exp = 1.0
real(8) :: b_hill = 0.5
logical :: use_exp_inhibit = .true.
logical :: use_base_apop_rate = .true.
contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine ReadMcParameters(nfin)
integer :: nfin
integer :: iuse_baserate, iuse_exp

read(nfin,*) baseRate
read(nfin,*) mitRate
read(nfin,*) Kaber
read(nfin,*) Klethal
read(nfin,*) K_ATM(1)
read(nfin,*) K_ATM(2)
read(nfin,*) K_ATM(3)
read(nfin,*) K_ATM(4)
read(nfin,*) K_ATR(1)
read(nfin,*) K_ATR(2)
read(nfin,*) K_ATR(3)
read(nfin,*) K_ATR(4)
read(nfin,*) Pcomplex
read(nfin,*) PHRsimple
read(nfin,*) KmaxInhibit
read(nfin,*) b_exp
read(nfin,*) b_hill
read(nfin,*) iuse_exp
read(nfin,*) iuse_baserate
use_base_apop_rate = (iuse_baserate == 1)
use_exp_inhibit = (iuse_exp == 1)
end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
function finhibit(C) result(fin)
real(8) :: C, fin

if (use_exp_inhibit) then
    fin = KmaxInhibit*(1 - exp(-b_exp*C))
else
    fin = KmaxInhibit*C/(b_hill + C)
endif
end function

!--------------------------------------------------------------------------
! repRate:          depends on pathway
! fidRate:          depends on phase, pathway
! misrepRate:       depends on phase (pathway?)
! pathwayUsed:      depends on phase, pathway
!--------------------------------------------------------------------------
subroutine setupRadiation

pathwayUsed = .false.
pathwayUsed(G1_phase,1:2) = .true.
pathwayUsed(S_phase,1:5) = .true.
pathwayUsed(G2_phase,3:5) = .true.

if (use_inhibiter) then
    pathwayUsed(G1_phase,6) = .true.
    pathwayUsed(S_phase,6) = .true.
    pathwayUsed(G2_phase,6) = .true.
endif
end subroutine

!--------------------------------------------------------------------------
! f is the fractional progress through S-phase
! DSB0[:] is the initial number of DSBs on each repair pathway
! totDSB0 is the total initial DSB count
! If the pathway is NHEJ, a fraction of the DSBs that would be on the NHEJ
! pathway are instead switched to the MMEJ pathway when the DNA-PK inhibitor
! drug is present.
! Let the inhibitor concentration = Cdrug. The fraction going to ALTEJ is finhibit(Cdrug).
!--------------------------------------------------------------------------
subroutine cellIrradiation(cp, dose, Cdrug)
type(cell_type), pointer :: cp
real(8) :: dose, Cdrug
integer :: ityp, kpar = 0
integer :: phase
real(8) :: DSB0(NP)
real(8) :: totDSB0, baseDSB, fin, T_S, f
real(8) :: Pbase, Pdie, R
real(8) :: DSB_Gy = 35

phase = cp%phase
cp%phase0 = phase
f = 0
if (phase == S_checkpoint) then
    f = 1.0
elseif (phase == S_phase) then
    T_S = cp%S_time - cp%G1S_time     ! duration of S_phase for this cell, which entered S at G1S_time
    f = (tnow - cp%G1S_time)/T_S
endif

baseDSB = DSB_Gy*dose
DSB0 = 0
if (phase <= G1_checkpoint) then
    totDSB0 = baseDSB
    DSB0(1) = (1 - Pcomplex)*totDSB0
    DSB0(2) = Pcomplex*totDSB0
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
elseif (phase >= G2_phase) then
    totDSB0 = 2*baseDSB
    DSB0(3) = (1 - Pcomplex)*totDSB0
    DSB0(4) = Pcomplex*totDSB0
else    ! S_phase
    DSB0(1) = (1 - Pcomplex)*(1-f)*baseDSB
    DSB0(2) = Pcomplex*(1-f)*baseDSB
    DSB0(3) = (1 - Pcomplex)*2*f*baseDSB*(1 - PHRsimple)
    DSB0(4) = Pcomplex*2*f*baseDSB
    DSB0(5) = (1 - Pcomplex)*2*f*baseDSB*PHRsimple
    totDSB0 = (1+f)*baseDSB
endif
fin = finhibit(Cdrug)
DSB0(6) = fin*sum(DSB0(1:3))
DSB0(1:3) = (1 - fin)*DSB0(1:3)
cp%DSB = DSB0
cp%totDSB0 = totDSB0
cp%Nlethal = 0
if (kcell_now == -4) write(*,'(a,2i6,8f6.1)') 'IR: kcell, phase,DSB: ',kcell_now,phase,cp%DSB

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
!--------------------------------------------------------------------------
subroutine updateATM(pATM,ATM_DSB,dth)
real(8) :: pATM, ATM_DSB, dth
real(8) :: r, k

r = K_ATM(2)*ATM_DSB   ! rate of production of pATM
k = K_ATM(3)  ! decay rate constant
pATM = r/k + (pATM - r/k)*exp(-k*dth)
if (kcell_now == -4) write(*,'(a,i4,5f8.4)') 'updateATM: r,k,r/k,ATM_DSB,pATM: ',kcell_now,r,k,r/k,ATM_DSB,pATM
end subroutine

!--------------------------------------------------------------------------
! Try making pATR production rate saturate, or limit pATR?
! Note: parameters from mcradio assume time is hours, therefore the time
! step passed to the subroutine, dth, has been converted to hours.
!--------------------------------------------------------------------------
subroutine updateATR(pATR,ATR_DSB,dth)
real(8) :: pATR, ATR_DSB, dth
real(8) :: r, k, x, xmax, km, xf

r = K_ATR(2)*ATR_DSB   ! rate of production of pATR
k = K_ATR(3)  ! decay rate constant
pATR = r/k + (pATR - r/k)*exp(-k*dth)
!if (kcell_now == 3) write(*,'(a,i4,5f8.4)') 'updateATR: r,k,r/k,ATR_DSB,pATR: ',kcell_now,r,k,r/k,ATR_DSB,pATR
!if (pATR > 0) write(*,'(a,i4,5f8.4)') 'updateATR: r,k,r/k,ATR_DSB,pATR: ',kcell_now,r,k,r/k,ATR_DSB,pATR
end subroutine

!------------------------------------------------------------------------
! Effect of pATM
! Using parameters from mcradio, time computed in hours, returned in secs
!------------------------------------------------------------------------
function S_checkpoint_time(cp) result(t)
type(cell_type), pointer :: cp
real(REAL_KIND) :: t, th

th = K_ATM(4)*(1 - exp(-K_ATM(1)*cp%pATM))
if (kcell_now == -4) write(*,'(a,i6,f8.3)') 'S_checkpoint_time (h): ',kcell_now,th
t = 3600*th
end function

!------------------------------------------------------------------------
! Effect of pATM and pATR
! Using parameters from mcradio, time computed in hours, returned in secs
!------------------------------------------------------------------------
function G2_checkpoint_time(cp) result(t)
type(cell_type), pointer :: cp
real(REAL_KIND) :: t, th, th_ATM, th_ATR

th_ATM = K_ATM(4)*(1 - exp(-K_ATM(1)*cp%pATM))
th_ATR = K_ATR(4)*(1 - exp(-K_ATR(1)*cp%pATR))
th = th_ATM + th_ATR
!if (kcell_now == 3) write(*,'(a,i6,2f8.3)') 'G2_checkpoint_time (ATM, ATR) (h): ',kcell_now,th_ATM,th_ATR
t = 3600*th
end function

!--------------------------------------------------------------------------
! To determine repair on a given pathway
! Using parameters from mcradio, dth in hours
!--------------------------------------------------------------------------
subroutine pathwayRepair(path, dth, N0, N)
integer :: path
real(8) :: dth, N0, N
real(8) :: Kreduction = 1.0

if (dth >= 0) then
    N = N0*exp(-repRate(path)*dth*Kreduction)  !!!! TESTING reducing repair rate
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
real(8) :: DSB(NP), dNlethal
real(8) :: DSB0(NP)
real(8) :: totDSB0, totDSB, Pmis, Nmis, totDSBinfid0, totDSBinfid, ATR_DSB, ATM_DSB, dth
logical :: pathUsed(NP)
integer :: k
logical :: use_DSBinfid = .true.

dth = dt/3600   ! hours
phase = cp%phase
DSB = cp%DSB
DSB0 = DSB     ! initial DSBs for this time step
totDSB0 = sum(DSB0)
if (totDSB0 == 0) return

!if (kcell_now == 1) write(*,'(a,2i6,8f6.2)') 'updateRepair: kcell, phase, DSB0, Nlethal: ',kcell_now,cp%phase,cp%DSB,cp%Nlethal
ATM_DSB = DSB(2) + DSB(4)   ! complex DSB
ATR_DSB = DSB(4) + DSB(5)
if (phase == G1_phase) then
    call updateATR(cp%pATR,ATR_DSB,dth)     ! updates the pATR mass through time step = dth
elseif (phase == S_phase) then
    call updateATM(cp%pATM,ATM_DSB,dth)     ! updates the pATM mass through time step = dth
    call updateATR(cp%pATR,ATR_DSB,dth)     ! updates the pATR mass through time step = dth
elseif (phase == G2_phase) then
    call updateATM(cp%pATM,ATM_DSB,dth)     ! updates the pATM mass through time step = dth
    call updateATR(cp%pATR,ATR_DSB,dth)     ! updates the pATR mass through time step = dth
endif
if (kcell_now == -4) write(*,'(a,2i6,2e12.3)') 'pATM, pATR: ',kcell_now,phase,cp%pATM,cp%pATR

pathUsed = pathwayUsed(phase,:)
DSB = 0
do k = 1,NP
    if (pathUsed(k)) then
        call pathwayRepair(k, dth, DSB0(k), DSB(k))
    endif
enddo
totDSB = sum(DSB)
totDSBinfid0 = 0
totDSBinfid = 0
do k = 1,NP
    if (pathUsed(k) .and. fidRate(k) < 1.0) then
        totDSBinfid0 = totDSBinfid0 + DSB0(k)
        totDSBinfid = totDSBinfid + DSB(k)
    endif
enddo
if (use_DSBinfid) then
    Pmis = misrepairRate(totDSBinfid0, totDSBinfid, eta)
else
    Pmis = misrepairRate(totDSB0, totDSB, eta)
endif
Nmis = 0
do k = 1,NP
    if (pathUsed(k) .and. fidRate(k) < 1.0) then
        Nmis = Nmis + (DSB0(k) - DSB(k))*(1 - fidRate(k)*(1 - Pmis))
    endif
enddo
cp%DSB = DSB
dNlethal = Klethal*misrepRate(phase)*Nmis   ! (was 0.5)  1.65 needs to be another parameter -> Klethal
cp%Nlethal = cp%Nlethal + dNlethal
end subroutine

!------------------------------------------------------------------------
! Clonogenic survival probability at first mitosis (McMahon, mcradio)
! cp%phase0 is the cell's phase at IR
!------------------------------------------------------------------------
subroutine survivalProbability(cp)
type(cell_type), pointer :: cp
real(8) :: DSB(NP), totDSB, Nlethal,Paber, Pbase, Papop, Pmit, Psurv

DSB = cp%DSB
totDSB = sum(DSB)
Nlethal = cp%Nlethal
if (cp%phase0 == G1_phase) then
!    Pbase = exp(-baseRate*cp%totDSB0)   ! this is 1 - Pdie
    Paber = exp(-Kaber*Nlethal)
    Pmit = exp(-mitRate*totDSB)
!       cp%Psurvive = Paber*Pbase*Papop*Pmit
    cp%Psurvive = Paber*Pmit
elseif (cp%phase0 < M_phase) then
    Paber = exp(-Kaber*Nlethal)
    Pmit = exp(-mitRate*totDSB)
    cp%Psurvive = Pmit*Paber        
else    ! M_phase or dividing
    Pmit = exp(-mitRate*totDSB)
    cp%Psurvive = Pmit*Msurvival
endif
NPsurvive = NPsurvive + 1
Psurvive(NPsurvive) = cp%Psurvive

end subroutine

end module
