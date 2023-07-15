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

subroutine cell_step(kcell, conc) bind(C)
!DEC$ ATTRIBUTES DLLEXPORT :: cell_step  
use, intrinsic :: iso_c_binding
integer(c_int), value :: kcell
real(c_double) :: conc(*)

!write(*,'(a,3f6.3)') 'conc: ',conc(1:3)
cell_list(kcell)%gconc = conc(1:3)

end subroutine


end module