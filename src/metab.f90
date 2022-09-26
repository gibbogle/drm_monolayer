! To test simple cell metabolism model
! Concentration of ATP varies by < 10%  https://en.wikipedia.org/wiki/Glycolysis#Intermediates_for_other_pathways

! Units:
!     time				s = seconds
!     distance			cm 
!     volume			cm^3
!     mass				micromole = 10^-6 mol = mumol
!     flux				mumol/s
!     concentration		mumol/cm^3 = mM
!
! Need to comment out 'use chemokine' when used in the test program metab.exe
!
! Question: How do the results of this model translate into cell rate of volume growth? 
!
! This version includes lumped amino acids - AA
!----------------------------------------------------------------------------------------------------------------
module FCN_mod
use real_kind_mod
implicit none

real(REAL_KIND) :: FCN_a1, FCN_b1, FCN_c1
real(REAL_KIND) :: FCN_r_Pm_base, FCN_fPn, FCN_Km_P, FCN_fGn, FCN_r_G, FCN_N_GA, FCN_V, FCN_K1, FCN_K2, FCN_C_L

end module

!----------------------------------------------------------------------------------------------------------------
module metabolism
use real_kind_mod
use global
#if .not. EXCEL
use chemokine
#endif
implicit none

integer, parameter :: FGP_SOLVER_MAXATP_TANDEM = 1
integer, parameter :: FGP_SOLVER_MAXATP_STAGED = 2
integer, parameter :: FGP_SOLVER_SURVIVAL_STAGED = 3

logical :: use_glutamine, from_excel, solved 

real(REAL_KIND), parameter :: f_I_threshold = 0.5	! NOT USED
real(REAL_KIND) :: I_threshold
! From spheroid-abm, the max rates of consumption of oxygen and glucose are: 
!   oxygen:  6.25e-17 moles/cell/s
!   glucose: 6.80e-17 moles/cell/s
! We work with mumol/cell/sec, and convert these values by scaling by 1.0e6, to give
!   oxygen:  6.25e-11 mumol/cell/s
!   glucose: 6.80e-11 mumol/cell/s

real(REAL_KIND) :: Hill_Km_O2
real(REAL_KIND) :: Hill_N_O2
real(REAL_KIND) :: Hill_Km_G	! Hill Km for dependence of glycolysis rate on glucose
real(REAL_KIND) :: Hill_N_G		! Hill N for dependence of glycolysis rate on glucose 
real(REAL_KIND) :: G_maxrate, O2_maxrate

#if EXCEL
real(REAL_KIND), parameter :: PI = 4*atan(1.0)
#endif
integer :: knt

contains

!--------------------------------------------------------------------------
! Currently this sets up parameters for type 1 cells only.
! For tandem case we need to incorporate glutamine.
!--------------------------------------------------------------------------
subroutine SetupMetabolism(mp,ok)
type(metabolism_type), pointer :: mp
logical :: ok

ok = .true.
Hill_Km_O2 = chemo(OXYGEN)%MM_C0
Hill_N_O2 = chemo(OXYGEN)%Hill_N
O2_maxrate = chemo(OXYGEN)%max_cell_rate
Hill_Km_G = chemo(GLUCOSE)%MM_C0
Hill_N_G = chemo(GLUCOSE)%Hill_N
G_maxrate = chemo(GLUCOSE)%max_cell_rate
write(nflog,'(a,2e12.3)') 'O2_maxrate, G_maxrate: ',O2_maxrate,G_maxrate

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine get_metab_rates(cp, Cin, C_GlnEx, res)
integer :: res
type(cell_type), pointer :: cp
real(REAL_KIND) :: Cin(:), C_GlnEx

call f_metab(cp,Cin(:), C_GlnEx, res)
end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine f_metab(cp, Cin, C_GlnEx, res)
integer :: res
real(REAL_KIND) :: Cin(:), C_GlnEx
type(cell_type), pointer :: cp
type(metabolism_type), pointer :: mp
real(REAL_KIND) :: C_O2, C_G, metab_O2, metab_G

!cp => cell_list(kcell_now)     ! BAD can't use a global variable that changes and use OMP 
mp => cp%metab
res = 0
C_O2 = max(0.0,Cin(OXYGEN))
C_G = max(0.0,Cin(GLUCOSE))
metab_O2 = O2_metab(C_O2)
metab_G = glucose_metab(C_G)
mp%O_rate = metab_O2*O2_maxrate		! consumption
mp%G_rate = metab_G*G_maxrate		! consumption

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
function f_MM(C,Km,N) result(v)
real(REAL_KIND) :: C, Km, v
integer :: N

v = C**N/(Km**N + C**N)
end function

!----------------------------------------------------------------------------------
! Computes metabolism rate as a fraction of the maximum cell rate
! Use the "soft landing" option for Hill_N = 1 if MM_threshold = 0
!----------------------------------------------------------------------------------
function O2_metab(C) result(metab)
integer :: ichemo
real(REAL_KIND) :: C
real(REAL_KIND) :: metab

ichemo = OXYGEN
if (ichemo == OXYGEN) then
	if (chemo(ichemo)%Hill_N == 2) then
		if (C > 0) then
			metab = C*C/(chemo(ichemo)%MM_C0*chemo(ichemo)%MM_C0 + C*C)
		else
			metab = 0
		endif
	else
		if (C > 0) then
			metab = C/(chemo(ichemo)%MM_C0 + C)
		else
			metab = 0
		endif
	endif
endif
end function

!----------------------------------------------------------------------------------
! Computes metabolism rate as a fraction of the maximum cell rate
!----------------------------------------------------------------------------------
function glucose_metab(C) result(metab)
real(REAL_KIND) :: C, metab
real(REAL_KIND) :: Kmin, Kmax, Kmm1
real(REAL_KIND) :: Vmax1, Vmax2, Kmm2, n1, n2
real(REAL_KIND) :: fV = 0.6
real(REAL_KIND) :: fK = 0.08
real(REAL_KIND) :: fboost = 2
real(REAL_KIND) :: Cboost = 0.1
logical :: variable_Km = .false.
logical :: double_Km = .false.
logical :: use_boost = .false.

if (C == 0) then
	metab = 0
	return
endif
if (use_boost) then

elseif (double_Km) then
	Kmm1 = Hill_Km_G
	Kmm2 = fK*Kmm1
	n1 = Hill_N_G
	n2 = 1
	metab = fV*C**n1/(Kmm1**n1 + C**n1) + (1 - fV)*C**n2/(Kmm2**n2 + C**n2)
elseif (variable_Km) then
	Kmax = Hill_Km_G	! These are completely arbitrary values
	Kmin = Kmax/15
	Kmm1 = 1*Kmin
	metab = C*(Kmm1 + C)/(Kmin*Kmm1 + Kmax*C + C*(Kmm1 + C))
else
	metab = C**Hill_N_G /(C**Hill_N_G + Hill_Km_G**Hill_N_G)
endif
end function

#if 0
!#if EXCEL
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
#endif


end module


