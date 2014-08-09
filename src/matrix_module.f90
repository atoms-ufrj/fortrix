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

module matrix_module

use base_matrix_module

implicit none

type, extends(base_matrix) :: matrix
  contains
    procedure, pass :: print => matrix_print
end type matrix

! Define a polymorphic constructor for type(matrix):
interface new_matrix
  module procedure :: new_matrix_0d
  module procedure :: new_matrix_1d
  module procedure :: new_matrix_2d
end interface new_matrix

interface ones
  module procedure :: ones_1d
  module procedure :: ones_2d
end interface ones

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function new_matrix_0d( a ) result( b )
    real(rb), intent(in) :: a
    type(matrix)         :: b
    b % rows = 1
    b % cols = 1
    b % size = 1
    allocate( b % ptr )
    call mproc % allocate( b % ptr, 1 )
    call mproc % upload( 1, [a], b % ptr )
  end function new_matrix_0d

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function new_matrix_1d( a ) result( b )
    real(rb), intent(in) :: a(:)
    type(matrix)         :: b
    integer(c_int) :: m
    m = size(a)
    b % rows = m
    b % cols = 1
    b % size = m
    allocate( b % ptr )
    call mproc % allocate( b % ptr, m )
    call mproc % upload( m, a, b % ptr )
  end function new_matrix_1d

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function new_matrix_2d( a ) result( b )
    real(rb), intent(in) :: a(:,:)
    type(matrix)         :: b
    integer(c_int) :: m, n
    m = size(a,1)
    n = size(a,2)
    b % rows = m
    b % cols = n
    b % size = m*n
    allocate( b%ptr )
    call mproc % allocate( b % ptr, b % size )
    call mproc % upload( b % size, a, b % ptr )
  end function new_matrix_2d

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function ones_1d( m ) result( b )
    integer, intent(in) :: m
    type(matrix)        :: b
    b % rows = m
    b % cols = 1
    b % size = m
    allocate( b % ptr )
    call mproc % allocate( b % ptr, m )
    call mproc % replicate( m, mproc % one, b % ptr )
  end function ones_1d

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function ones_2d( m, n ) result( b )
    integer, intent(in) :: m, n
    type(matrix)        :: b
    b % rows = m
    b % cols = n
    b % size = m*n
    allocate( b%ptr )
    call mproc % allocate( b % ptr, b % size )
    call mproc % replicate( b % size, mproc % one, b % ptr )
  end function ones_2d

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine matrix_print( x, fmt )
    class(matrix), intent(in) :: x
    character(*),  intent(in) :: fmt
    real(rb), allocatable :: a(:,:)
    allocate( a(x%rows,x%cols) )
    call mproc % download( x%size, x % ptr, a )
    call mproc % print( x%rows, x%cols, a, fmt )
  end subroutine matrix_print

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end module matrix_module
