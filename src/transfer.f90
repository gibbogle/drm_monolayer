! Transferring results to the GUI

module transfer

use global
use chemokine
use metabolism
use cycle_mod
use, intrinsic :: iso_c_binding

#include "../src/version.h"

implicit none

type, bind(C) :: celldata_type
	integer(c_int) :: tag
	real(c_double) :: radius
	real(c_double) :: centre(3)
	integer(c_int) :: celltype
	integer(c_int) :: status
end type

type, bind(C) :: fielddata_type
    integer(c_int) :: NX, NY, NZ, NCONST
    real(c_double) :: DX
    type(c_ptr) :: Conc_ptr   ! Cslice(NX,NY,NZ,NCONST)
    integer(c_int) :: ncells
    type(c_ptr) :: cell_ptr
end type

real(REAL_KIND), allocatable :: SFlookup(:)
real(REAL_KIND) :: dp_survive

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_DLL_build_version(version_array,array_len) BIND(C) 
!DEC$ ATTRIBUTES DLLEXPORT :: get_dll_build_version
use, intrinsic :: iso_c_binding
character(c_char) :: version_array(12)
integer(c_int) :: array_len
integer :: k

dll_version = DLL_BUILD_VERSION
gui_version = GUI_BUILD_VERSION
!write(nflog,*) 'get_DLL_build_version: ',dll_version
do k = 1,12
	version_array(k) = dll_version(k:k)
!	write(nflog,'(i2,a,a)') k,' ',version_array(k)
	if (version_array(k) == ' ') then
		version_array(k) = char(0)
		array_len = k
		exit
	endif
enddo
!write(nflog,*) 'array_len: ',array_len
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine get_dimensions(nsteps_dim, deltat, ntdisplay, maxchemo, nextra, cused) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dimensions
use, intrinsic :: iso_c_binding
integer(c_int) :: nsteps_dim, ntdisplay, maxchemo, nextra
real(c_double) :: deltat
logical(c_bool) :: cused(*)
integer :: ichemo

write(nflog,*) 'get_dimensions: MAX_CHEMO: ',MAX_CHEMO
nsteps_dim = nsteps
deltat = DELTA_T
ntdisplay = NT_DISPLAY
maxchemo = MAX_CHEMO
nextra = N_EXTRA
do ichemo = 1,MAX_CHEMO
	cused(ichemo+1) = chemo(ichemo)%used
enddo
cused(1) = .true.			! CFSE
cused(MAX_CHEMO+2) = .true.	! Growth rate
end subroutine


!-----------------------------------------------------------------------------------------
! Rendered cell colour may depend on stage, state, receptor expression level.
! col(:) = (r,g,b)
!-----------------------------------------------------------------------------------------
subroutine CellColour(kcell,highlight,col)
integer :: kcell, highlight, col(3)
integer :: stage, status
integer, parameter :: WHITE(3) = (/255,255,255/)
integer, parameter :: RED(3) = (/255,0,0/)
integer, parameter :: GREEN(3) = (/0,255,0/)
integer, parameter :: BLUE(3) = (/0,0,255/)
integer, parameter :: DEEPRED(3) = (/200,0,0/)
integer, parameter :: DEEPBLUE(3) = (/30,20,255/)
integer, parameter :: DEEPGREEN(3) = (/0,150,0/)
integer, parameter :: LIGHTRED(3) = (/255,70,90/)
integer, parameter :: LIGHTBLUE(3) = (/0,200,255/)
integer, parameter :: LIGHTGREEN(3) = (/50,255,150/)
integer, parameter :: DEEPORANGE(3) = (/240,70,0/)
integer, parameter :: LIGHTORANGE(3) = (/255,130,0/)
integer, parameter :: YELLOW(3) = (/255,255,0/)
integer, parameter :: DEEPPURPLE(3) = (/180,180,30/)
integer, parameter :: LIGHTPURPLE(3) = (/230,230,100/)
integer, parameter :: DEEPBROWN(3) = (/130,70,0/)
integer, parameter :: LIGHTBROWN(3) = (/200,100,0/)
integer, parameter :: GRAY(3) = (/128,128,128/)

integer, parameter :: Qt_white = 3
integer, parameter :: Qt_black = 2
integer, parameter :: Qt_red = 7
integer, parameter :: Qt_darkRed = 13
integer, parameter :: Qt_green = 8
integer, parameter :: Qt_darkGreen = 14
integer, parameter :: Qt_blue = 9
integer, parameter :: Qt_darkBlue = 15
integer, parameter :: Qt_cyan = 10
integer, parameter :: Qt_darkCyan = 16
integer, parameter :: Qt_magenta = 11
integer, parameter :: Qt_darkMagenta = 17
integer, parameter :: Qt_yellow = 12
integer, parameter :: Qt_darkYellow = 18
integer, parameter :: Qt_gray = 5
integer, parameter :: Qt_darkGray = 4
integer, parameter :: Qt_lightGray = 6

if (highlight == 0) then
    col = LIGHTORANGE
else
    col = LIGHTRED
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Pack the colours (r,g,b) into an integer.
!-----------------------------------------------------------------------------------------
integer function rgb(col)
integer :: col(3)

rgb = ishft(col(1),16) + ishft(col(2),8) + col(3)
end function

!-----------------------------------------------------------------------------------------
! EC(:) is concentrations at the bottom, next to the cell layer
! cmedium(:) is average medium concentrations
!-----------------------------------------------------------------------------------------
subroutine getMediumConc(EC,cmedium)
real(REAL_KIND) :: EC(:), cmedium(:)

EC = Caverage(MAX_CHEMO+1:2*MAX_CHEMO)
cmedium = Cmediumave(:)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine getNecroticFraction(necrotic_fraction,totvol_cm3)
real(REAL_KIND) :: necrotic_fraction, totvol_cm3	! vol_cm3 not used here, needed in scell
real(REAL_KIND) :: cellvol_cm3, dvol
!necrotic_fraction = (Nsites-Ncells)/real(Nsites)
cellvol_cm3 = Ncells*DELTA_X**3
dvol = totvol_cm3-cellvol_cm3
necrotic_fraction = dvol/totvol_cm3
if (necrotic_fraction < 0.005) necrotic_fraction = 0
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
function rint(i) result(r)
integer :: i
real(REAL_KIND) :: r
r = i
end function

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine getHypoxicCount(nhypoxic)
integer :: nhypoxic(3)
integer :: kcell, i

nhypoxic = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	do i = 1,3
		if (cell_list(kcell)%Cin(OXYGEN) < O2cutoff(i)) nhypoxic(i) = nhypoxic(i) + 1
	enddo
enddo
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine getClonoHypoxicCount(nclonohypoxic)
integer :: nclonohypoxic(3)
integer :: kcell, i, idrug
logical :: tagged

nclonohypoxic = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	if (cell_list(kcell)%state == DYING) cycle
	if (cell_list(kcell)%radiation_tag) cycle
	tagged = .false.
	do idrug = 1,MAX_DRUGTYPES
	    if (cell_list(kcell)%drug_tag(idrug)) tagged = .true.
	enddo
	if (tagged) cycle
	do i = 1,3
		if (cell_list(kcell)%Cin(OXYGEN) < O2cutoff(i)) then
			nclonohypoxic(i) = nclonohypoxic(i) + 1
		endif
	enddo
enddo
end subroutine

!--------------------------------------------------------------------------------
! Need to compare growth rate with a fraction of average growth rate
! Should ngrowth(i) include cells at checkpoints?
! nogrow is the count of viable cells that are not growing
!--------------------------------------------------------------------------------
subroutine getGrowthCount(ngrowth, nogrow, nphase,nmutations, nclono)
integer :: ngrowth(3), nogrow(:), nphase(:), nmutations, nclono
integer :: kcell, i, ityp, iphase, nch1, nch2, ip
real(REAL_KIND) :: r_mean(2), ps1, ps2, ps, SF, rclono, alfa
real(REAL_KIND) :: N50 = 5.64	! number of divisions to get 50 cells NOT CORRECT
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp

if (use_cell_cycle) then
    r_mean = max_growthrate
else
    r_mean = Vdivide0/(2*divide_time_mean)
endif
ngrowth = 0
nogrow = 0
nphase = 0
nmutations = 0
rclono = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	ityp = cp%celltype
	if (cp%state == DEAD) cycle
!	if (.not.(cp%state == DYING) .and. cp%metab%A_rate < r_Ag) then
!		nogrow(ityp) = nogrow(ityp) + 1
!	endif
	do i = 1,3
		if (cp%dVdt < growthcutoff(i)*r_mean(ityp) .or. &
		    cp%phase == G1_checkpoint .or. cp%phase == S_checkpoint .or. cp%phase == G2_checkpoint) then
		    ngrowth(i) = ngrowth(i) + 1
		endif
	enddo
	iphase = min(cp%phase,M_phase)
	nphase(iphase) = nphase(iphase) + 1
!	nmutations = nmutations + cp%N_Ch1 + cp%N_Ch2
	ccp => cc_parameters(ityp)
#if 0
	if (.not.cp%state == DYING) then
		nch1 = cp%N_Ch1
		if (nch1 > 0) then
			ps1 = ccp%psurvive_Ch1**nch1
		else
			ps1 = 1
		endif
		nch2 = cp%N_Ch2
		if (nch2 > 0) then
			ps2 = ccp%psurvive_Ch2**nch2
		else
			ps2 = 1
		endif
		ps = ps1*ps2
		if (ps == 1) then
			rclono = rclono + 1
		elseif (ps > 0) then
!			rclono = rclono + ps**N50
!		    call clonogenic(ps,50,SF)
            ip = ps/dp_survive
            alfa = ps - ip*dp_survive
            SF = (1-alfa)*SFlookup(ip) + alfa*SFlookup(ip+1)
			rclono = rclono + SF
		endif
	endif
#endif
enddo
nclono = rclono + 0.5
end subroutine

!--------------------------------------------------------------------------------
! From the probability of survival of division, estimate the probability of
! growing a colony with more than Nlim cells.
!--------------------------------------------------------------------------------
subroutine clonogenic(ityp, p_survive, Nlim, SF)
real(REAL_KIND) :: p_survive, SF
integer :: ityp, Nlim
real(REAL_KIND) :: R, p_d
integer :: nc, irun, nruns, idiv, ndiv, n, ntemp, i, kpar = 0
logical :: ended, dbug
integer :: n_colony_days = 10

p_d = 1 - p_survive
dbug = .false.
nruns = 10000
ndiv = 24*3600*n_colony_days/divide_time_mean(ityp) + 1
nc = 0
do irun = 1,nruns
	n = 1
	ended = .false.
	do idiv = 1,ndiv
		ntemp = n
		do i = 1,ntemp
!			call random_number(R)
            R = par_uni(kpar)
			if (R < p_d) then
				n = n - 1
				if (n == 0) then
					ended = .true.
					exit
				endif
			else
				n = n + 1
			endif
	        if (n >= Nlim) then
	            nc = nc + 1
	            ended = .true.
	            exit
	        endif
		enddo
	    if (dbug) then
	        write(*,*) 'irun, n: ',irun,n
	    endif
		if (ended) exit
	enddo
!	if (n == 0) then
!		kbox = 0
!	else
!		if (boxsize == 1) then
!			kbox = n
!		else
!			kbox = n/boxsize + 1
!		endif
!	endif
!	if (kbox > nboxes) then
!		write(*,*) 'Error: irun, n, kbox, nboxes: ',irun, n, kbox, nboxes
!		stop
!	endif
!	ncol(kbox) = ncol(kbox) + 1
enddo
SF = nc/real(nruns)
end subroutine

!--------------------------------------------------------------------------------
! Generate and store a lookup table for SF as a function of p (division survival prob).
!--------------------------------------------------------------------------------
subroutine GenerateSFlookup(ityp)
integer :: ityp
integer :: np_survive
integer :: Nlim = 50
integer :: i
real(REAL_KIND) :: p_survive

write(nflog,*)
write(nflog,*) 'GenerateSFlookup: ityp: ',ityp
dp_survive = 0.01
np_survive = 100
if (allocated(SFlookup)) deallocate(SFlookup)
allocate(SFlookup(0:np_survive))
do i = 0,np_survive
    p_survive = i*dp_survive
    if (i == 0) then
        SFlookup(i) = 0
    elseif (i == np_survive) then
        SFlookup(i) = 1
    else
        call clonogenic(ityp, p_survive, Nlim, SFlookup(i))
    endif
!    write(nflog,'(i4,2f8.5)') i,p_survive,SFlookup(i)
enddo
write(nflog,*) 'did GenerateSFlookup: ityp: ',ityp
end subroutine

!--------------------------------------------------------------------------------
! Compute total uptake rate for a constituent
! NOT USED
!--------------------------------------------------------------------------------
subroutine sum_dMdt(ichemo)
integer :: ichemo
integer :: kcell, Nh, Nc
real(REAL_KIND) :: C, metab, dMdt, asum, msum, Csum

if (ichemo > GLUCOSE) then
	write(*,*) 'Error: sum_dMdt: only for oxygen and glucose'
	stop
endif
Nh = chemo(ichemo)%Hill_N
asum = 0
Csum = 0
msum = 0
Nc = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	Nc = Nc + 1
	C = cell_list(kcell)%Cin(ichemo)
	Csum = Csum + C
	metab = C**Nh/(chemo(ichemo)%MM_C0**Nh + C**Nh)
	msum = msum + metab
	dMdt = metab*chemo(ichemo)%max_cell_rate 
	asum = asum + dMdt
enddo
total_dMdt = total_dMdt + asum
!write(*,'(a,2i6,2e12.3)') 'sum_dMdt: ',ichemo,Nc,asum,total_dMdt*3600 
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine getNviable	!(Nviable, Nlive)
!integer :: Nviable(:) 
integer :: Nlive(MAX_CELLTYPES)
integer :: kcell, ityp, idrug, nd
logical :: tag
type(cell_type), pointer :: cp

Nviable = 0
Nlive = 0
nd = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) then
	    nd = nd+1
	    cycle
	endif
    ityp = cp%celltype
    Nlive(ityp) = Nlive(ityp) + 1
!    if (cp%state /= DYING) write(*,*) 'getNviable: ',kcell,cp%state
	if (cp%state == DYING) cycle
	Nviable(ityp) = Nviable(ityp) + 1
enddo
write(nflog,'(a,6i6)') 'getNviable: ',Nviable(1),Ncells_type(1),nd,Napop,Ndying(1),Nmitotic
if (Nlive(1) /= Ncells_type(1)) then
	write(*,'(a,5i8)') 'Error: getNviable: Nlive /= Ncells_type(1), nd, Napop, Nmitotic: ',Nlive(1),Ncells_type(1),nd,Napop, Nmitotic
	write(nflog,'(a,5i8)') 'Error: getNviable: Nlive /= Ncells_type(1), nd, Napop, Nmitotic: ',Nlive(1),Ncells_type(1),nd,Napop, Nmitotic
	stop
endif
if (Nviable(1) /= Ncells_type(1) - Ndying(1)) then
	write(*,'(a,4i8)') 'Error: getNviable: Nviable /= Ncells_type(1) - Ndying, Nmitotic: ',Nviable(1),Ncells_type(1),Ndying(1),Nmitotic
	write(nflog,'(a,4i8)') 'Error: getNviable: Nviable /= Ncells_type(1) - Ndying, Nmitotic: ',Nviable(1),Ncells_type(1),Ndying(1),Nmitotic
	stop
endif

end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
function getstatus(cp) result(status)
type(cell_type), pointer :: cp
integer :: status

!if (cp%anoxia_tag) then
!	status = 2	! tagged to die of anoxia
!elseif (cp%aglucosia_tag) then
!	status = 4	! tagged to die of aglucosia
!if (cp%ATP_tag) then
!	status = 2	! tagged to die of low ATP
!else
if (cp%radiation_tag) then
	status = 10
elseif (cp%drug_tag(1)) then
	status = 11
elseif (cp%drug_tag(1)) then
	status = 12
elseif (cp%Cin(OXYGEN) < hypoxia_threshold) then
	status = 1	! radiobiological hypoxia
!elseif (cp%mitosis > 0) then
elseif (cp%V > 0.9*cp%divide_volume) then  ! just a surrogate for mitosis
	status = 3	! in mitosis
else
	status = 0
endif
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_constituents(nvars,cvar_index,nvarlen,name_array,narraylen) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_constituents
use, intrinsic :: iso_c_binding
character(c_char) :: name_array(0:*)
integer(c_int) :: nvars, cvar_index(0:*), nvarlen, narraylen
integer :: ivar, k, ichemo
character*(24) :: name
character(c_char) :: c

write(nflog,*) 'get_constituents'
nvarlen = 24
ivar = 0
k = ivar*nvarlen
cvar_index(ivar) = 0	! CFSE
name = 'CFSE'
call copyname(name,name_array(k),nvarlen)
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	ivar = ivar + 1
	k = ivar*nvarlen
	cvar_index(ivar) = ichemo
	name = chemo(ichemo)%name
	write(nflog,*) 'get_constituents: ',ichemo,name
	call copyname(name,name_array(k),nvarlen)
enddo
ivar = ivar + 1
k = ivar*nvarlen
cvar_index(ivar) = GROWTH_RATE
name = 'Growth rate'
call copyname(name,name_array(k),nvarlen)
ivar = ivar + 1
k = ivar*nvarlen
cvar_index(ivar) = CELL_VOLUME
name = 'Cell volume'
call copyname(name,name_array(k),nvarlen)
ivar = ivar + 1
k = ivar*nvarlen
cvar_index(ivar) = O2_BY_VOL
name = 'Cell O2xVol'
call copyname(name,name_array(k),nvarlen)
ivar = ivar + 1
k = ivar*nvarlen
cvar_index(ivar) = CYCLE_PHASE
name = 'Cycle phase'
call copyname(name,name_array(k),nvarlen)
nvars = ivar + 1
write(nflog,*) 'did get_constituents'
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine copyname(name,name_array,n)
character*(*) :: name
character :: name_array(*)
integer :: n
integer :: k

do k = 1,n
	name_array(k) = name(k:k)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Get number of live cells
!-----------------------------------------------------------------------------------------
subroutine get_nFACS(n) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_nfacs
use, intrinsic :: iso_c_binding
integer(c_int) :: n
integer :: k, kcell

!call logger('get_nFACS')
n = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	n = n+1
enddo
write(nflog,*) 'get_nFACS: n: ',n
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_FACS(facs_data,vmin,vmax,vmin_log,vmax_log, scale_drug_by_vol) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_facs 
use, intrinsic :: iso_c_binding
real(c_double) :: facs_data(*)
real(c_double) :: vmin(*), vmax(*)
real(c_double) :: vmin_log(*), vmax_log(*)
logical(c_bool),value :: scale_drug_by_vol
integer :: k, kcell, iextra, ichemo, ivar, nvars, var_index(32)
real(REAL_KIND) :: cfse_min, val, val_log
type(cell_type), pointer :: cp

call logger('get_FACS')
nvars = 1	! CFSE
var_index(nvars) = 0
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	nvars = nvars + 1
	var_index(nvars) = ichemo
!	write(*,*) 'ichemo: ',ichemo,' nvars: ',nvars,'  ',chemo(ichemo)%name
enddo
do iextra = 1,N_EXTRA-1
	nvars = nvars + 1
	var_index(nvars) = MAX_CHEMO + iextra
enddo
!write(nflog,*) 'nvars: ',nvars 
vmin(1:nvars) = 1.0e10
vmax(1:nvars) = -1.0e10
!vmin_log(1:nvars) = 1.0e10
!vmax_log(1:nvars) = -1.0e10
cfse_min = 1.0e20
k = 0
do kcell = 1,nlist
    cp =>cell_list(kcell)
	if (cp%state == DEAD) cycle
	do ivar = 1,nvars
		ichemo = var_index(ivar)
		if (ichemo == 0) then
			val = cp%CFSE
			cfse_min = min(val,cfse_min)
		elseif (ichemo <= MAX_CHEMO) then
			val = cp%Cin(ichemo)
			if (scale_drug_by_vol) then
			    val = val*cp%V/1.0e-9
			endif
		elseif (ichemo == GROWTH_RATE) then
			val = cp%dVdt/1.0e-9		! -> pL
		elseif (ichemo == CELL_VOLUME) then
			val = cp%V/1.0e-9
		elseif (ichemo == O2_BY_VOL) then
			val = cp%V*cp%Cin(OXYGEN)/1.0e-9
		elseif (ichemo == CYCLE_PHASE) then
			val = cp%phase
		endif
		k = k+1
		facs_data(k) = val
		vmin(ivar) = min(val,vmin(ivar))
		vmax(ivar) = max(val,vmax(ivar))
!		if (val <= 1.0e-8) then
!			val_log = -8
!		else
!			val_log = log10(val)
!		endif
!		vmin_log(ivar) = min(val_log,vmin_log(ivar))
!		vmax_log(ivar) = max(val_log,vmax_log(ivar))
	enddo
enddo
call logger('Did get_FACS')
end subroutine

!-----------------------------------------------------------------------------------------
! nhisto is the number of histogram boxes
! vmax(ivar) is the maximum value for variable ivar
! Probably need to adjust vmax() to a roundish value
!
! Compute 3 distributions: 1 = both cell types
!                          2 = type 1
!                          3 = type 2
! Stack three cases in vmax() and histo_data()
!-----------------------------------------------------------------------------------------
subroutine get_histo(nhisto, histo_data, vmin, vmax, histo_data_log, vmin_log, vmax_log, scale_drug_by_vol) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_histo
use, intrinsic :: iso_c_binding
integer(c_int),value :: nhisto
real(c_double) :: vmin(*), vmax(*), histo_data(*)
real(c_double) :: vmin_log(*), vmax_log(*), histo_data_log(*)
logical(c_bool),value :: scale_drug_by_vol
real(REAL_KIND) :: val, val_log
integer :: n(3), i, ih, k, kcell, ict, ichemo, ivar, nvars, var_index(32)
integer,allocatable :: cnt(:,:,:)
real(REAL_KIND),allocatable :: dv(:,:), valmin(:,:), valmax(:,:)
integer,allocatable :: cnt_log(:,:,:)
integer :: phase_cnt(3,7), ivar_phase
real(REAL_KIND),allocatable :: dv_log(:,:), valmin_log(:,:), valmax_log(:,:)
type(cell_type), pointer :: cp
!real(REAL_KIND) :: vmin_log(100), vmax_log(100)
!real(REAL_KIND),allocatable :: histo_data_log(:)

!call logger('get_histo')
nvars = 1	! CFSE
var_index(nvars) = 0
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	nvars = nvars + 1
	var_index(nvars) = ichemo
enddo
nvars = nvars + 1
var_index(nvars) = GROWTH_RATE
nvars = nvars + 1
var_index(nvars) = CELL_VOLUME
nvars = nvars + 1
var_index(nvars) = O2_BY_VOL
nvars = nvars + 1
var_index(nvars) = CYCLE_PHASE
ivar_phase = nvars

allocate(cnt(3,nvars,nhisto))
allocate(dv(3,nvars))
allocate(valmin(3,nvars))
allocate(valmax(3,nvars))
allocate(cnt_log(3,nvars,nhisto))
allocate(dv_log(3,nvars))
allocate(valmin_log(3,nvars))
allocate(valmax_log(3,nvars))
!allocate(histo_data_log(10000))
cnt = 0
valmin = 1.0e10
valmax = -1.0e10
cnt_log = 0
valmin_log = 1.0e10
valmax_log = -1.0e10
n = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	ict = cp%celltype
	do ivar = 1,nvars
		ichemo = var_index(ivar)
		if (ichemo == 0) then
			val = cp%CFSE
		elseif (ichemo <= MAX_CHEMO) then
			val = cp%Cin(ichemo)
			if (scale_drug_by_vol .and. ichemo >= DRUG_A) then
			    val = val*cp%V/1.0e-9
			endif
		elseif (ichemo == GROWTH_RATE) then
			val = cp%dVdt/1.0e-9		! -> pL
		elseif (ichemo == CELL_VOLUME) then
			val = Vcell_pL*cp%V/1.0e-9
		elseif (ichemo == O2_BY_VOL) then
			val = cp%Cin(OXYGEN)*Vcell_pL*cp%V/1.0e-9
		elseif (ichemo == CYCLE_PHASE) then
			val = cp%phase
		endif
		valmax(ict+1,ivar) = max(valmax(ict+1,ivar),val)	! cell type 1 or 2
		valmax(1,ivar) = max(valmax(1,ivar),val)			! both
		if (val <= 1.0e-8) then
			val_log = -8
		else
			val_log = log10(val)
		endif
		valmin_log(ict+1,ivar) = min(valmin_log(ict+1,ivar),val_log)	! cell type 1 or 2
		valmin_log(1,ivar) = min(valmin_log(1,ivar),val_log)			! both
		valmax_log(ict+1,ivar) = max(valmax_log(ict+1,ivar),val_log)	! cell type 1 or 2
		valmax_log(1,ivar) = max(valmax_log(1,ivar),val_log)			! both
	enddo
	n(ict+1) = n(ict+1) + 1
	n(1) = n(1) + 1
enddo
do ivar = 1,nvars
	ichemo = var_index(ivar)
	if (ichemo == CELL_VOLUME) then
		valmin(:,ivar) = Vcell_pL*0.8
		valmin_log(:,ivar) = log10(Vcell_pL*0.8)
	else
		valmin(:,ivar) = 0
	endif
enddo

dv = (valmax - valmin)/nhisto
!write(nflog,*) 'dv'
!write(nflog,'(e12.3)') dv 
dv_log = (valmax_log - valmin_log)/nhisto
!write(nflog,*) 'dv_log'
!write(nflog,'(e12.3)') dv_log
phase_cnt = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	ict = cp%celltype
	do ivar = 1,nvars
		ichemo = var_index(ivar)
		if (ichemo == 0) then
			val = cp%CFSE
		elseif (ichemo <= MAX_CHEMO) then
			val = cp%Cin(ichemo)
			if (scale_drug_by_vol .and. ichemo >= DRUG_A) then
			    val = val*cp%V/1.0e-9
			endif
		elseif (ichemo == GROWTH_RATE) then
			val = cp%dVdt/1.0e-9		! -> pL
		elseif (ichemo == CELL_VOLUME) then
			val = Vcell_pL*cp%V/1.0e-9
		elseif (ichemo == O2_BY_VOL) then
			val = cp%Cin(OXYGEN)*Vcell_pL*cp%V/1.0e-9
		elseif (ichemo == CYCLE_PHASE) then
			val = cell_list(kcell)%phase
			phase_cnt(1,cp%phase) = phase_cnt(1,cp%phase) + 1
			phase_cnt(ict+1,cp%phase) = phase_cnt(ict+1,cp%phase) + 1
		endif
		k = (val-valmin(1,ivar))/dv(1,ivar) + 1
		k = min(k,nhisto)
		k = max(k,1)
		cnt(1,ivar,k) = cnt(1,ivar,k) + 1
		k = (val-valmin(ict+1,ivar))/dv(ict+1,ivar) + 1
		k = min(k,nhisto)
		k = max(k,1)
		cnt(ict+1,ivar,k) = cnt(ict+1,ivar,k) + 1
		if (val <= 1.0e-8) then
			val_log = -8
		else
			val_log = log10(val)
		endif
		k = (val_log-valmin_log(1,ivar))/dv_log(1,ivar) + 1
		k = min(k,nhisto)
		k = max(k,1)
		cnt_log(1,ivar,k) = cnt_log(1,ivar,k) + 1
		k = (val_log-valmin_log(ict+1,ivar))/dv_log(ict+1,ivar) + 1
		k = min(k,nhisto)
		k = max(k,1)
		cnt_log(ict+1,ivar,k) = cnt_log(ict+1,ivar,k) + 1
	enddo
enddo

do i = 1,3
	if (n(i) == 0) then
		vmin((i-1)*nvars+1:i*nvars) = 0
		vmax((i-1)*nvars+1:i*nvars) = 0
		histo_data((i-1)*nvars*nhisto+1:i*nhisto*nvars) = 0
		vmin_log((i-1)*nvars+1:i*nvars) = 0
		vmax_log((i-1)*nvars+1:i*nvars) = 0
		histo_data_log((i-1)*nvars*nhisto+1:i*nhisto*nvars) = 0
	else
		do ivar = 1,nvars
			vmin((i-1)*nvars+ivar) = valmin(i,ivar)
			vmax((i-1)*nvars+ivar) = valmax(i,ivar)
			do ih = 1,nhisto
				k = (i-1)*nvars*nhisto + (ivar-1)*nhisto + ih
				histo_data(k) = (100.*cnt(i,ivar,ih))/n(i)
			enddo
			vmin_log((i-1)*nvars+ivar) = valmin_log(i,ivar)
			vmax_log((i-1)*nvars+ivar) = valmax_log(i,ivar)
			do ih = 1,nhisto
				k = (i-1)*nvars*nhisto + (ivar-1)*nhisto + ih
				histo_data_log(k) = (100.*cnt_log(i,ivar,ih))/n(i)
			enddo
		enddo
	endif
	! Now overwrite CYCLE_PHASE data
	ivar = ivar_phase
	do ih = 1,nhisto
		k = (i-1)*nvars*nhisto + (ivar-1)*nhisto + ih
		if (ih <= 7) then
			histo_data(k) = (100.*phase_cnt(i,ih))/n(i)
		else
			histo_data(k) = 0
		endif
	enddo
enddo
deallocate(cnt)
deallocate(dv)
deallocate(valmin)
deallocate(valmax)
deallocate(cnt_log)
deallocate(dv_log)
deallocate(valmin_log)
deallocate(valmax_log)
call logger('did get_histo')
end subroutine

!--------------------------------------------------------------------------------
! Always the z axis
!--------------------------------------------------------------------------------
subroutine get_concdata(nvars, ns, dxc, ex_conc) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_concdata
use, intrinsic :: iso_c_binding
integer(c_int) :: nvars, ns
real(c_double) :: dxc, ex_conc(0:*)
real(REAL_KIND) :: A, d, cave(MAX_CHEMO)
integer :: k, ks, ichemo, offset

!call logger('get_concdata')
nvars = MAX_CHEMO
ns = N1D
A = well_area
d = total_volume/A
dxc = d/(ns-1)
cave = 0
do ks = 1,ns
	do ichemo = 1, nvars
		offset = ichemo*ns
		k = offset - 1 + ks
		ex_conc(k) = chemo(ichemo)%Cmedium(ks)
		cave(ichemo) = cave(ichemo) + ex_conc(k)
	enddo
enddo
write(nflog,'(a,4e12.3)') 'cave: ',cave/ns
end subroutine

!-----------------------------------------------------------------------------------------
! Note:
! DRUG_A = TPZ_DRUG e.g. SN30000
! DRUG_B = DNB_DRUG e.g. PR104A 
!
! IC = average intracellular concentration
! EC = average medium concentration <<<<< NOTE
!-----------------------------------------------------------------------------------------
subroutine get_values(nvars,varID,ysim)
!DEC$ ATTRIBUTES DLLEXPORT :: get_values
integer :: nvars
character*(24) :: varID(nvars)
real(REAL_KIND) :: ysim(nvars)
integer :: ivar, ityp
!integer :: Nviable(MAX_CELLTYPES)
integer :: Nlive(MAX_CELLTYPES)
real(REAL_KIND) :: plate_eff(MAX_CELLTYPES)

do ivar = 1,nvars
	if (varID(ivar) == 'OXYGEN_EC') then
		ysim(ivar) = Cmediumave(OXYGEN)
	elseif (varID(ivar) == 'GLUCOSE_EC') then
		ysim(ivar) = Cmediumave(GLUCOSE)
		write(*,*) 'get_values: GLUCOSE_EC: ',ivar,ysim(ivar)
!	elseif (varID(ivar) == 'LACTATE_EC') then
!		ysim(ivar) = Cmediumave(LACTATE)
!	elseif (varID(ivar) == 'DRUG_A_EC') then
!		ysim(ivar) = Cmediumave(DRUG_A)
!	elseif (varID(ivar) == 'DRUG_A_METAB1_EC') then
!		ysim(ivar) = Cmediumave(DRUG_A+1)
!	elseif (varID(ivar) == 'DRUG_A_METAB2_EC') then
!		ysim(ivar) = Cmediumave(DRUG_A+2)
!	elseif (varID(ivar) == 'DRUG_B_EC') then
!		ysim(ivar) = Cmediumave(DRUG_B)
!	elseif (varID(ivar) == 'DRUG_B_METAB1_EC') then
!		ysim(ivar) = Cmediumave(DRUG_B+1)
!	elseif (varID(ivar) == 'DRUG_B_METAB2_EC') then
!		ysim(ivar) = Cmediumave(DRUG_B+2)
	elseif (varID(ivar) == 'OXYGEN_IC') then
		ysim(ivar) = Caverage(OXYGEN)
	elseif (varID(ivar) == 'GLUCOSE_IC') then
		ysim(ivar) = Caverage(GLUCOSE)
!	elseif (varID(ivar) == 'LACTATE_IC') then
!		ysim(ivar) = Caverage(LACTATE)
	elseif (varID(ivar) == 'DRUG_A_IC') then
		ysim(ivar) = Caverage(DRUG_A)
	elseif (varID(ivar) == 'DRUG_A_METAB1_IC') then
		ysim(ivar) = Caverage(DRUG_A+1)
!	elseif (varID(ivar) == 'DRUG_A_METAB2_IC') then
!		ysim(ivar) = Caverage(DRUG_A+2)
!	elseif (varID(ivar) == 'DRUG_B_IC') then
!		ysim(ivar) = Caverage(DRUG_B)
!	elseif (varID(ivar) == 'DRUG_B_METAB1_IC') then
!		ysim(ivar) = Caverage(DRUG_B+1)
!	elseif (varID(ivar) == 'DRUG_B_METAB2_IC') then
!		ysim(ivar) = Caverage(DRUG_B+2)
	elseif (varID(ivar) == 'NCELLS') then
		ysim(ivar) = Ncells			! for now, total live cells
	elseif (varID(ivar) == 'PE') then
!		call getNviable(Nviable, Nlive)
		Nlive = Ncells_type
		do ityp = 1,Ncelltypes
			if (Nlive(ityp) > 0) then
				plate_eff(ityp) = real(Nviable(ityp))/Nlive(ityp)
			else
				plate_eff(ityp) = 0
			endif
		enddo
		ysim(ivar) = plate_eff(1)	! for now, just type 1 cells
	elseif (varID(ivar) == 'RADIATION') then
		ysim(ivar) = -1
	else
		write(*,*) 'varID is not in the list of possible IDs: ',varID(ivar)
		stop
	endif
enddo
end subroutine



!-----------------------------------------------------------------------------------------
! Computes the fluorescence distribution that would result if cells were killed and stained
! with propidium iodide (PI).
! What should be done with cells that are already dead?
!-----------------------------------------------------------------------------------------
subroutine get_PI_dist(nbins, fract, max_fract, max_fluor) bind(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_pi_dist
use, intrinsic :: iso_c_binding
integer(c_int), value :: nbins
real(c_double) :: fract(*), max_fract, max_fluor
integer :: kcell, i
real :: fluor, dfluor, maxf
type(cell_type), pointer :: cp

maxf = 0
max_fluor = 4
dfluor = max_fluor/nbins
fract(1:nbins) = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	if (cp%phase <= G1_checkpoint) then
		fluor = 1
	elseif (cp%phase == S_phase) then
!		fluor = 1 + min((tnow - cp%G1S_time)/(cp%S_time - cp%G1S_time), 1.0)
!DRM        fluor = 1 + min(cp%S_time/cp%S_duration, 1.0)
	else
		fluor = 2
	endif
	fluor = fluor*cp%V/Vcell_cm3
	maxf = max(maxf,fluor)
	i = min(int(fluor/dfluor + 1),nbins)
	fract(i) = fract(i) + 1
enddo
fract(1:nbins) = fract(1:nbins)/sum(fract(1:nbins))
max_fract = maxval(fract(1:nbins))
write(nflog,*) 'get_PI_dist: max_fract, maxf: ',max_fract, maxf
!write(nflog,'(10f7.4)') fract(1:nbins)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_string(bufptr) bind(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_string
use, intrinsic :: iso_c_binding
type(c_ptr) :: bufptr
character(c_char) :: buf(1024)
integer :: buflen
character*(1024), save, target :: string

string = 'A test string'
buflen = len(trim(string))
!write(*,*) 'buflen: ',buflen
string(buflen+1:buflen+1) = char(0)
bufptr = c_loc(string)
end subroutine


end module

