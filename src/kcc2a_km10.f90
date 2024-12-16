module km10_mod
implicit none

real(8) :: CC_factor

contains
!-------------------------------------------------------------------------------------------
! Find x (kcc) such that f(x) = tmitosis(CC_tot,x,kmccp) - Tph = 0
! x(n+1) = x(n) - f(x(n))/f'(x(n))
!-------------------------------------------------------------------------------------------
subroutine newton(x0,t0,CC_tot,kmccp,Tph)
real(8) :: x0, t0, CC_tot, kmccp, Tph
real(8) :: f0, f1, dfdx, dx, x1, x,t
integer :: n

dx = 0.01
n = 0
x = x0
t = t0
do
    n = n+1
    t = tmitosis(CC_tot,x,kmccp)
    f0 = t - Tph
    dfdx = (tmitosis(CC_tot,x+dx,kmccp) - t)/dx
    if (dfdx == 0) then
        write(*,*) 'dfdx = 0'
        exit
    endif
    x1 = x - f0/dfdx
    if (abs(x-x1)/x < 0.002) exit
!    write(*,'(i4,4f8.4)') n,dfdx,x0,x1,f0
    x = x1
    t = tmitosis(CC_tot,x,kmccp)
    if (n > 20) then
        write(*,*) 'newton, n > 20'
        write(*,'(a,3f8.3)') 'Initial x0, t0, Tph: ',x0,t0,Tph
        write(*,'(a,2f8.3)') 'Final x, t: ',x,t
        stop
    endif
enddo
x0 = x
t0 = t
end subroutine

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
function tmitosis(CC_tot,kcc,km) result(t)
real(8) :: CC_tot, kcc, km
real(8) :: t, y, dt, dydt
integer :: it, Nt = 5000

dt = 0.002
t = 0
y = 0
do it = 1,Nt
    dydt = (Kcc + y)*(CC_tot-y)/(km + (CC_tot-y))
    dydt = max(dydt,0.0)
    y = y + dydt*dt
    t = t + dt
    if (y > CC_factor*CC_tot) then
        exit
    endif
enddo
end function

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
function get_Kcc(kmccp, CC_tot, CC_threshold_factor, T_G2) result(kcc)
real(8) :: kmccp, CC_tot, CC_threshold_factor, T_G2, kcc
real(8) :: x0, t0
real(8),parameter :: alfa = -1.0, beta = 0.45

CC_factor = CC_threshold_factor
!x0 = alfa + beta*kmccp      ! initial guess
x0 = 0.5
t0 = tmitosis(CC_tot,x0,kmccp)
!write(*,*) 'x0, t0, T_G2: ',x0,t0,T_G2
call newton(x0,t0,CC_tot,kmccp,1.0*T_G2)    ! find x0 = kcc such that tmitosis = T_G2
kcc = x0
!write(*,*) 'kcc: ',kcc
end function

end module