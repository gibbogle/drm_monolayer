module greens_mod

use global

implicit none

contains

subroutine set_greens(Ngreen) bind(C)
!DEC$ ATTRIBUTES DLLEXPORT :: set_greens  
use, intrinsic :: iso_c_binding
integer(c_int), value :: Ngreen

if (Ngreen == 0) then
    greens = .false.
else
    greens = .true.
    NgreenCells = Ngreen
endif

end subroutine

subroutine get_times(dt_drm) bind(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_times  
use, intrinsic :: iso_c_binding
real(c_float) :: dt_drm
dt_drm = DELTA_T/60     ! time step in minutes
end subroutine

subroutine cell_step(kcell, conc) bind(C)
!DEC$ ATTRIBUTES DLLEXPORT :: cell_step  
use, intrinsic :: iso_c_binding
integer(c_int), value :: kcell
real(c_double) :: conc(*)

!write(*,'(a,3f6.3)') 'conc: ',conc(1:3)
cell_list(kcell)%gconc = conc(1:3)

end subroutine

subroutine write_binary(a, n) bind(C)
!DEC$ ATTRIBUTES DLLEXPORT :: write_binary 
use, intrinsic :: iso_c_binding
integer, parameter :: nf = 22
integer(c_int), value :: n
real(c_float) :: a(n)

write(*,*) 'do write'
open(nf, file = 'binary.dat', status = 'replace', form = 'unformatted', access = 'stream')
write(nf) a
close(nf)
write(*,*) 'did write'

a = 0
write(*,*) 'do read'
open(nf, file = 'binary.dat', status = 'old', form = 'unformatted', access = 'stream')
read(nf) a
write(*,*) 'did read'
write(*,*) a(1:10)

end subroutine

end module