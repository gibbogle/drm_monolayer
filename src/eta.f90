module eta_module

use global

implicit none

integer, parameter :: Neta = 11
integer, parameter :: NtIR = 9
real(8) :: dsigma_dt
real(8) :: eta_table(Neta,NtIR,5)

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
real(8) :: V, R, S, eta
integer :: k, it

write(nflog,*) 'make_eta_table: ',S_NHEJ, S_TMEJ
do k = 1,Neta
do it = 1,NtIR
	V = 1 + (k-1.0)/(Neta-1.0)
	R = V**(1./3.)
	S = S_NHEJ + (it-1)*dsigma_dt
	eta = etafun(R,S)
	eta_table(k,it,1:2) = eta
	eta = etafun(R,S_TMEJ)
	eta_table(k,it,3:5) = eta
	write(nflog,'(a,2i4,2f10.6)') 'k, it, S, eta_NHEJ: ', k, it, S, eta_table(k,it,1)
enddo
enddo
end subroutine

!--------------------------------------------------------------------------
! f_S = fractional progress through S-phase. =0 in G1, =1 in G2
! Two sigma values have been defined: sigma_NHEJ and sigma_TMEJ
! For each pathway, lookup tables for eta as a function of f_S, tIR have been 
! precomputed, requiring only interpolation.
! The number of table entries is Neta*NtIR (= 11*9)
!--------------------------------------------------------------------------
function eta_lookup(phase,path,f_S,tIR) result(eta)
integer :: phase, path
real(8) :: f_S, tIR, eta
real(8) :: fk, fi
integer :: k1, k2, it1, it2

!if (phase == G1_phase) then
!    eta = eta_table(1,path)
!elseif (phase > S_phase) then
!    eta = eta_table(Neta,path)
!else
!    k1 = f_S*(Neta-1) + 1
!    if (k1 < Neta) then
!        k2 = k1+1
!        f = f_S - real(k1-1)/(Neta-1)
!        eta = (1 - f)*eta_table(k1,path) + f*eta_table(k2,path)
!    else
!        eta = eta_table(Neta,path)
!    endif
!endif
if (tIR >= NtIR - 1) then
    it1 = NtIR-1
    it2 = NtIR
    fi = 1
else
    it1 = tIR + 1   ! rounded down
    it2 = it1 + 1
    fi = tIR - real(it1-1)/(NtIR-1)
endif
if (phase == G1_phase) then
    k1 = 1
    k2 = 2
    fk = 0
elseif (phase > S_phase) then
    k1 = Neta - 1
    k2 = Neta
    fk = 1.0
else
    k1 = f_S*(Neta-1) + 1
    if (k1 < Neta) then
        k2 = k1+1
        fk = f_S - real(k1-1)/(Neta-1)
    else
        k1 = Neta - 1
        k2 = Neta
        fk = 1.0 
    endif
endif
eta = (1-fk)*(1-fi)*eta_table(k1,it1,path) + fk*(1-fi)*eta_table(k2,it1,path) &
    + (1-fk)*fi*eta_table(k1,it2,path) + fk*fi*eta_table(k2,it2,path)
!write(*,'(a,2f7.3,4i3,e12.3)') 'eta_lookup: f_S,tIR,k1,k2,it1,it2,eta: ',f_S,tIR,k1,k2,it1,it2,eta
end function

end module