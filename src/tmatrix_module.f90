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

module tmatrix_module

use base_matrix_module

implicit none

type, extends(base_matrix) :: tmatrix
  contains
    procedure, pass :: print => tmatrix_print
end type tmatrix

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine tmatrix_print( x, fmt )
    class(tmatrix), intent(in) :: x
    character(*),  intent(in)  :: fmt
    real(rb), allocatable :: a(:,:)
    type(c_ptr) :: y
    call mproc % allocate( y, x%size )
    call mproc % transpose( x%rows, x%cols, x%ptr, y )
    allocate( a(x%cols,x%rows) )
    call mproc % download( x%size, y, a )
    call mproc % deallocate( y )
    call mproc % print( x%cols, x%rows, a, fmt )
  end subroutine tmatrix_print

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end module tmatrix_module
