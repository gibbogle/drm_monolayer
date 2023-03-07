!-------------------------------------------------------------------------------------------
! Find x (km10) such that f(x) = tmitosis(CC_tot,kcc2a,x) - Tph = 0
! x(n+1) = x(n) - f(x(n))/f'(x(n))
!-------------------------------------------------------------------------------------------
subroutine newton(x0,t0,CC_tot,kcc2a,Tph)
real(8) :: x0, t0, CC_tot, kcc2a, Tph
real(8) :: f0, f1, dfdx, dx
integer :: n

dx = 0.01
n = 0
do
    n = n+1
    t0 = tmitosis(CC_tot,kcc2a,x0)
    f0 = t0 - Tph
    dfdx = (tmitosis(CC_tot,kcc2a,x0+dx) - t0)/dx
    if (dfdx == 0) exit
    x1 = x0 - f0/dfdx
    if (abs(x0-x1)/x0 < 0.001) exit
!    write(*,'(i4,4f8.4)') n,dfdx,x0,x1,f0
    x0 = x1
    t0 = tmitosis(CC_tot,kcc2a,x0)
    if (n > 10) stop
enddo
end subroutine

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
function tmitosis(CC_tot,kcc2a,km10) result(t)
real(8) :: CC_tot, kcc2a, km10
real(8) :: t, y, dt, dydt
integer :: it, Nt = 5000

dt = 0.002
t = 0
y = 0
do it = 1,Nt
    dydt = (Kcc2a + y)*(CC_tot-y)/(Km10 + (CC_tot-y))
    dydt = max(dydt,0.0)
    y = y + dydt*dt
    t = t + dt
    if (y > 0.9*CC_tot) then
        exit
    endif
enddo
end function

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
subroutine fit(n,x,y,alfa,beta)
real(8) :: x(n), y(n), alfa, beta
integer :: n
real(8) :: xave, yave, num, den
integer :: i

xave = sum(x)/n
yave = sum(y)/n

!write(*,*) 'xave,yave: ',xave,yave
num = 0
den = 0
do i = 1,n
    num = num + (x(i)-xave)*(y(i)-yave)
    den = den + (x(i)-xave)**2
enddo
!write(*,*) 'num, den: ',num,den
beta = num/den
alfa = yave - beta*xave

end subroutine

!-------------------------------------------------------------------------------------------
! To determine the relationship between kcc2a and km10 for given CC_tot in order to
! generate a specified time to CC threshold with no radiation dose.
!-------------------------------------------------------------------------------------------
subroutine kcc2a_km10(CC_tot, alfa, beta)
implicit none
real(8) :: CC_tot, alfa(3), beta(3)
integer, parameter :: nka = 8
real(8) :: Tphase(3) = [5.636, 8.556, 3.951]    ! G1, S, G2
real(8) :: km10, dt, t, y, dydt, dkm, t0, x0
real(8) :: kcc2a(nka), km10sol(nka)
real(8) :: tmitosis
integer :: iph, ika, ikm, it, Nt=1000

write(*,*) 'kcc2a_km10'
dt = 0.01
dkm = 1
!CC_tot = 5
do iph = 1,3
    do ika = 1,nka
        kcc2a(ika) = 1 + ika*0.5
        do ikm = 1,30
            km10 = ikm*dkm
            t = tmitosis(CC_tot,kcc2a(ika),km10)
!            write(*,'(i4,3f8.3)') iph, kcc2a(ika), km10, t
            if (t > Tphase(iph)) then
                ! now zero in on km10 to solve t = Tphase(iph)
                t0 = t
                x0 = km10
                call newton(x0,t0,CC_tot,kcc2a(ika),Tphase(iph))    ! find km10 = x0 s.t. t = Tphase(iph)
                km10sol(ika) = x0
!                write(*,'(a,i4,2f8.3)') 'Solved km10: ',iph,kcc2a(ika),km10sol(ika)
                exit
            endif
        enddo
    enddo
!    write(*,'(a,i4,f8.1)') 'For iph, CC_tot: ',iph,CC_tot
!    write(*,'(a,10f8.3)') 'kcc2a: ',kcc2a
!    write(*,'(a,10f8.3)') 'km10 : ',km10sol
    ! linear fit
    call fit(nka,kcc2a,km10sol,alfa(iph),beta(iph))   ! km10 = alfa + beta*kcc2a
!    write(*,'(a,2f8.4)') 'alfa, beta: ',alfa(iph),beta(iph)
!    write(*,*)
enddo
end subroutine