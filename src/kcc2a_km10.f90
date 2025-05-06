module km10_mod
implicit none

real(8) :: CC_factor
logical :: ddbug

contains
!-------------------------------------------------------------------------------------------
! Find x (kcc) such that f(x) = tmitosis(CC_tot,x,kmccp) - Tph = 0
! x(n+1) = x(n) - f(x(n))/f'(x(n))
!-------------------------------------------------------------------------------------------
subroutine newton(x0,t0,CC_tot,kmccp,Tph)
real(8) :: x0, t0, CC_tot, kmccp, Tph
real(8) :: f0, f1, dfdx, dx, x1, x, t, t1
integer :: n
integer :: Nits = 30

dx = 0.02   ! was 0.01
n = 0
x = x0
t = t0
do
    n = n+1
    t = tmitosis(CC_tot,x,kmccp)
    f0 = t - Tph
    t1 = tmitosis(CC_tot,x+dx,kmccp)
    dfdx = (t1 - t)/dx
    if (ddbug) write(*,'(a,i4,5e12.3)') 'n,x,t,t1,f0,dfdx: ',n,x,t,t1,f0,dfdx
    if (dfdx == 0) then
        write(*,*) 'dfdx = 0'
        stop
    endif
    x1 = x - f0/dfdx
    x1 = max(x1,0.01)   ! prevent x going < 0
    if (abs(x-x1)/x < 0.002) exit
    if (ddbug) write(*,'(i4,4f8.4)') n,dfdx,x0,x1,f0
    x = x1
    t = tmitosis(CC_tot,x,kmccp)
    if (abs(Tph - t) < 0.005) exit
    if (n > Nits) then
        write(*,*) 'newton, n > Nits: ',Nits
        write(*,'(a,3f10.5)') 'Initial x0, t0, Tph: ',x0,t0,Tph
        write(*,'(a,2f10.5)') 'Final x, t: ',x,t
        write(*,'(a,3f10.5)') 'CC_tot, CC_factor,kmccp: ',CC_tot,CC_factor,kmccp
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
integer :: it, Nt = 10000

dt = 0.002
t = 0
y = 0
do it = 1,Nt
    dydt = (Kcc + y)*(CC_tot-y)/(km + (CC_tot-y))
!    dydt = max(dydt,0.0)
    if (dydt < 0) then
        write(*,*) 'dydt: ',dydt
        write(*,*) 'Kcc,Km,y,CC_tot-y: ',Kcc,Km,y,CC_tot-y
        stop
    endif
    y = y + dydt*dt
    if (ddbug) write(*,'(a,3e12.3)') 'y, dydt, dy: ',y,dydt,dydt*dt
    t = t + dt
    if (y > CC_factor*CC_tot) then
        exit
    endif
enddo
if (ddbug) write(*,*) 'it = ',it
end function


!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
function get_Kcc(kmccp, CC_tot, CC_threshold_factor, T_G2) result(kcc)
real(8) :: kmccp, CC_tot, CC_threshold_factor, T_G2, kcc
real(8) :: x0, t0
!real(8),parameter :: alfa = -1.0, beta = 0.45
real(8),parameter :: a = 0.015, b = 0.1737, c = 0.1075 

ddbug = (CC_threshold_factor < 0)
!ddbug = .true.
CC_factor = abs(CC_threshold_factor)
!x0 = alfa + beta*kmccp      ! initial guess
x0 = 0.5
if (kmccp < 1.0) x0 = 0.1
if (kmccp > 5.0) x0 = 2.0
!x0 = a*kmccp**2 + b*kmccp + c
t0 = tmitosis(CC_tot,x0,kmccp)
if (ddbug) write(*,*) 'x0, t0, T_G2: ',x0,t0,T_G2
call newton(x0,t0,CC_tot,kmccp,1.0*T_G2)    ! find x0 = kcc such that tmitosis = T_G2
kcc = x0
!write(*,*) 'kmccp, kcc: ',kmccp,kcc
end function

end module