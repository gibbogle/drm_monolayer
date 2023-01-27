module eta_module

use global

implicit none

integer, parameter :: Neta = 11
real(8) :: eta_table(Neta,5)

contains
!--------------------------------------------------------------------------
! From Eq 15 in McMahon2021
! S is sigma
!--------------------------------------------------------------------------
function theta15(R, S) result(thetaVal)
real(8) :: R, S
real(8) :: R2, S2
real(8) :: termOne, termTwo, termThree, scaling, thetaVal

R2 = R*R
S2 = S*S
termOne   = sqrt(2*pi)*R2*S*erf(R*sqrt(2.0)/S)
termTwo   = -exp(-4*R2/(2*S2))*(S2 -R2*S2)
termThree = S2 - 3*R2*S2
scaling   = 2*pi*S2/(R*R2)
thetaVal = (termOne + termTwo + termThree)*scaling
end function	

!--------------------------------------------------------------------------
! From Eq 14 in McMahon2021
! S is sigma
!--------------------------------------------------------------------------
function etafun(R,S) result(eta)
real(8) :: R, S, eta

eta = (6/(4*pi*R**3))*theta15(R,S)
end function	

!--------------------------------------------------------------------------
! S_NHEJ, S_TMEJ are sigma_NHEJ, sigma_TMEJ
!--------------------------------------------------------------------------
subroutine make_eta_table(S_NHEJ, S_TMEJ)
real(8) :: S_NHEJ, S_TMEJ
real(8) :: V, R, eta
integer :: k

do k = 1,Neta
	V = 1 + (k-1.0)/(Neta-1.0)
	R = V**(1./3.)
	eta = etafun(R,S_NHEJ)
	eta_table(k,1:2) = eta
	eta = etafun(R,S_TMEJ)
	eta_table(k,3:5) = eta
!	write(*,'(a,f8.4,2f6.2,f10.6)') 'sigma, V, R, eta: ', sigma, V, R, eta
enddo

end subroutine

!--------------------------------------------------------------------------
! f_S = fractional progress through S-phase. =0 in G1, =1 in G2
! Two sigma values have been defined: sigma_NHEJ and sigma_TMEJ
! For each pathway, lookup tables for eta as a function of f_S have been 
! precomputed, requiring only interpolation.
! The number of table entries is Neta (= 11), and 
! eta_table(1,path) = eta_G1(path), eta_table(Neta,path) = eta_G2(path)
!--------------------------------------------------------------------------
function eta_lookup(phase,path,f_S) result(eta)
integer :: phase, path
real(8) :: f_S, eta
real(8) :: f
integer :: k1, k2

if (phase == G1_phase) then
    eta = eta_table(1,path)    ! eta_G1(path)
elseif (phase > S_phase) then
    eta = eta_table(Neta,path) ! eta_G2(path)
else
    k1 = f_S*(Neta-1) + 1
    if (k1 < Neta) then
        k2 = k1+1
        f = f_S - real(k1-1)/(Neta-1)
        eta = (1 - f)*eta_table(k1,path) + f*eta_table(k2,path)
    else
        eta = eta_table(Neta,path)
    endif
endif
end function

end module