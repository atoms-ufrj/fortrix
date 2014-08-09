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

module dmatrix_module

use base_matrix_module

implicit none

type, extends(base_matrix) :: dmatrix
  contains
    procedure, pass :: print => dmatrix_print
end type dmatrix

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine dmatrix_print( x, fmt )
    class(dmatrix), intent(in) :: x
    character(*),    intent(in) :: fmt
    real(rb), allocatable :: a(:), b(:,:)
    integer :: i
    allocate( a(x%size), b(x%size,x%size) )
    call mproc % download( x%size, x % ptr, a )
    b = 0.0_rb
    forall (i=1:x%rows) b(i,i) = a(i)
    call mproc % print( x%size, x%size, b, fmt )
  end subroutine dmatrix_print

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end module dmatrix_module
