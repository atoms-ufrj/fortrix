!    This file is part of Fortrix (version 0.1)
!
!    Fortrix is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    Fortrix is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.!
!
!    You should have received a copy of the GNU General Public License
!    along with Fortrix.  If not, see <http://www.gnu.org/licenses/>.

module scalar_module

use base_matrix_module

implicit none

type, extends(base_matrix) :: scalar
  contains
    procedure, pass :: print => scalar_print
end type scalar

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function new_scalar( a ) result( b )
    real(rb), intent(in) :: a
    type(scalar)         :: b
    b%rows = 1
    b%cols = 1
    b%size = 1
    allocate( b%ptr )
    call mproc % allocate( b%ptr, 1 )
    call mproc % upload( 1, [a], b%ptr )
  end function new_scalar

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine scalar_print( x, fmt )
    class(scalar), intent(in) :: x
    character(*),  intent(in) :: fmt
    real(rb) :: a(1)
    call mproc % download( 1, x % ptr, a )
    call mproc % print( 1, 1, [a], fmt )
  end subroutine scalar_print

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end module scalar_module
