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

end module