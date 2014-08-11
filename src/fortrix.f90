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

module Fortrix

! TO DO LIST:
!
! - Test results of operations
! - Implement and test linear system solution

use iso_c_binding
use matrix_module
use dmatrix_module
use tmatrix_module
use scalar_module

implicit none

interface assignment(=)
  module procedure assign_matrix
  module procedure assign_tmatrix
  module procedure assign_dmatrix
  module procedure assign_scalar
  module procedure assign_matrix_from_tmatrix
  module procedure assign_matrix_from_scalar
end interface

interface operator (.d.)
  module procedure :: diag_matrix
  module procedure :: diag_tmatrix
end interface

interface Dinv
  module procedure :: inv_diag_matrix
  module procedure :: inv_diag_tmatrix
end interface

interface operator (.t.)
  module procedure :: transpose_matrix
  module procedure :: transpose_tmatrix
  module procedure :: transpose_dmatrix
end interface

! ALL REQUIRE TESTING:
interface operator (+)
  module procedure :: sum_matrix_matrix
  module procedure :: sum_tmatrix_tmatrix
  module procedure :: sum_dmatrix_dmatrix
  module procedure :: sum_scalar_scalar
  module procedure :: sum_matrix_tmatrix
  module procedure :: sum_tmatrix_matrix
  module procedure :: sum_matrix_dmatrix
  module procedure :: sum_dmatrix_matrix
  module procedure :: sum_tmatrix_dmatrix
  module procedure :: sum_dmatrix_tmatrix
end interface

! ALL REQUIRE TESTING:
interface operator (-)
  module procedure :: reciprocal_matrix
  module procedure :: reciprocal_scalar
  module procedure :: subtraction_matrix_matrix
  module procedure :: subtraction_tmatrix_tmatrix
  module procedure :: subtraction_dmatrix_dmatrix
  module procedure :: subtraction_scalar_scalar
  module procedure :: subtraction_matrix_tmatrix
  module procedure :: subtraction_tmatrix_matrix
  module procedure :: subtraction_matrix_dmatrix
  module procedure :: subtraction_dmatrix_matrix
  module procedure :: subtraction_tmatrix_dmatrix
  module procedure :: subtraction_dmatrix_tmatrix
end interface

! Hadamard products:
interface operator (.o.)
  module procedure :: hadamard_matrix_matrix
  module procedure :: hadamard_tmatrix_tmatrix
  module procedure :: hadamard_matrix_tmatrix
  module procedure :: hadamard_tmatrix_matrix
end interface

interface operator (*)
  ! Scaling:
  module procedure :: product_scalar_matrix
  module procedure :: product_matrix_scalar
  module procedure :: product_scalar_tmatrix
  module procedure :: product_tmatrix_scalar
  module procedure :: product_scalar_dmatrix
  module procedure :: product_dmatrix_scalar
  module procedure :: product_scalar_scalar
  ! Matrix product:
  module procedure :: product_matrix_matrix
  module procedure :: product_matrix_tmatrix
  module procedure :: product_tmatrix_matrix
  module procedure :: product_tmatrix_tmatrix
  module procedure :: product_matrix_dmatrix
  module procedure :: product_dmatrix_matrix
  module procedure :: product_tmatrix_dmatrix
  module procedure :: product_dmatrix_tmatrix
  module procedure :: product_dmatrix_dmatrix
end interface

interface operator (**)
  module procedure :: pow_matrix
  module procedure :: pow_tmatrix
  module procedure :: pow_scalar
end interface

interface exp
  module procedure :: exp_matrix
end interface

interface log
  module procedure :: log_matrix
end interface

interface Diag
  module procedure :: diag_matrix
  module procedure :: diag_tmatrix
end interface

interface inv
  module procedure :: inv_matrix
  module procedure :: inv_tmatrix
  module procedure :: inv_dmatrix
  module procedure :: inv_scalar
end interface

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Initialization and finalization:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine fortrix_startup( thread )
    integer, intent(in), optional :: thread
    if (present(thread)) then
      call mproc % initialize( thread )
    else
      call mproc % initialize( 1 )
    end if
  end subroutine fortrix_startup

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine fortrix_shutdown( )
    call mproc % finalize( )
  end subroutine fortrix_shutdown

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Printing:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine print_matrix( a, fmt )
    class(base_matrix), intent(in)           :: a
    character(*),       intent(in), optional :: fmt
    if (present(fmt)) then
      call a % print( fmt )
    else
      call a % print( "ES11.3" )
    end if
    if (a%volatile) call a % destroy()
  end subroutine print_matrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Assignment:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine assign_matrix( b, a )
    type(matrix),  intent(inout) :: b
    class(matrix), intent(in)    :: a
    call b % assign_from( a )
  end subroutine assign_matrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine assign_tmatrix( b, a )
    type(tmatrix),  intent(inout) :: b
    class(tmatrix), intent(in)    :: a
    call b % assign_from( a )
  end subroutine assign_tmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine assign_dmatrix( b, a )
    type(dmatrix),  intent(inout) :: b
    class(dmatrix), intent(in)    :: a
    call b % assign_from( a )
  end subroutine assign_dmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine assign_scalar( b, a )
    type(scalar),       intent(inout) :: b
    class(base_matrix), intent(in)    :: a
    if (a%size /= 1) call mproc % error("a scalar can only be assigned from a 1-by-1 matrix")
    call b % assign_from( a )
  end subroutine assign_scalar

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine assign_matrix_from_tmatrix( b, a )
    type(matrix),   intent(inout) :: b
    class(tmatrix), intent(in)    :: a
    if (b%volatile) then
      allocate( b%ptr )
      call mproc % allocate( b%ptr, a%size )
      b%volatile = .false.
    else if (b%size /= a%size) then
      call mproc % deallocate( b%ptr )
      call mproc % allocate( b%ptr, a%size )
    end if
    call mproc % transpose( a%rows, a%cols, a%ptr, b%ptr )
    if (a%volatile) call a % destroy()
    b%rows = a%cols
    b%cols = a%rows
    b%size = a%size
  end subroutine assign_matrix_from_tmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine assign_matrix_from_scalar( b, a )
    type(matrix),  intent(inout) :: b
    class(scalar), intent(in)    :: a
    call b % assign_from( a )
  end subroutine assign_matrix_from_scalar

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Diagonalization:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function diag_matrix( a ) result( b )
    class(matrix), intent(in) :: a
    type(dmatrix)             :: b
    if (a%cols /= 1_c_int) &
      call mproc % error("trying to apply the diagonal function to a non-vector matrix")
    call b % clone( a )
  end function diag_matrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function diag_tmatrix( a ) result( b )
    class(tmatrix), intent(in) :: a
    type(dmatrix)              :: b
    if (a%rows /= 1_c_int) &
      call mproc%error("trying to apply diagonal function to a non-vector matrix")
    call b % clone( a )
    b%rows = b%cols
    b%cols = 1_c_int
  end function diag_tmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function inv_diag_matrix( a ) result( b )
    class(matrix),  intent(in) :: a
    type(dmatrix)              :: b
    if (a%cols /= 1_c_int) &
      call mproc % error("trying to apply the diagonal function to a non-vector matrix")
    call b % clone( a )
    call mproc % inv( b%size, b%ptr )
  end function inv_diag_matrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function inv_diag_tmatrix( a ) result( b )
    class(tmatrix), intent(in) :: a
    type(dmatrix)              :: b
    if (a%rows /= 1_c_int) &
      call mproc%error("trying to apply diagonal function to a non-vector matrix")
    call b % clone( a )
    call mproc % inv( b%size, b%ptr )
    b%rows = b%cols
    b%cols = 1_c_int
  end function inv_diag_tmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Transposition:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function transpose_matrix( a ) result( b )
    class(matrix), intent(in) :: a
    type(tmatrix)             :: b
    call b % clone( a )
  end function transpose_matrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function transpose_tmatrix( a ) result( b )
    class(tmatrix), intent(in) :: a
    type(matrix)               :: b
    call b % clone( a )
  end function transpose_tmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function transpose_dmatrix( a ) result( b )
    class(dmatrix), intent(in) :: a
    type(dmatrix)              :: b
    call b % clone( a )
  end function transpose_dmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Matrix addition:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function sum_matrix_matrix( a, b ) result( c )
    class(matrix), intent(in) :: a
    class(matrix), intent(in) :: b
    type(matrix)              :: c
    call c % sum_mm( a, b )
  end function sum_matrix_matrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function sum_tmatrix_tmatrix( a, b ) result( c )
    class(tmatrix), intent(in) :: a
    class(tmatrix), intent(in) :: b
    type(tmatrix)              :: c
    call c % sum_mm( a, b )
  end function sum_tmatrix_tmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function sum_dmatrix_dmatrix( a, b ) result( c )
    class(dmatrix), intent(in) :: a
    class(dmatrix), intent(in) :: b
    type(dmatrix)              :: c
    call c % sum_mm( a, b )
  end function sum_dmatrix_dmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function sum_scalar_scalar( a, b ) result( c )
    class(scalar), intent(in) :: a
    class(scalar), intent(in) :: b
    type(scalar)              :: c
    call c % sum_mm( a, b )
  end function sum_scalar_scalar

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function sum_matrix_tmatrix( a, b ) result( c )
    class(matrix),  intent(in) :: a
    class(tmatrix), intent(in) :: b
    type(matrix)               :: c
    call c % sum_mt( a, b )
  end function sum_matrix_tmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function sum_tmatrix_matrix( a, b ) result( c )
    class(tmatrix), intent(in) :: a
    class(matrix),  intent(in) :: b
    type(matrix)               :: c
    call c % sum_mt( b, a )
  end function sum_tmatrix_matrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function sum_matrix_dmatrix( a, b ) result( c )
    class(matrix),  intent(in) :: a
    class(dmatrix), intent(in) :: b
    type(matrix)               :: c
    call c % sum_md( a, b )
  end function sum_matrix_dmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function sum_dmatrix_matrix( a, b ) result( c )
    class(dmatrix), intent(in) :: a
    class(matrix),  intent(in) :: b
    type(matrix)               :: c
    call c % sum_md( b, a )
  end function sum_dmatrix_matrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function sum_tmatrix_dmatrix( a, b ) result( c )
    class(tmatrix), intent(in) :: a
    class(dmatrix), intent(in) :: b
    type(tmatrix)              :: c
    call c % sum_md( a, b )
  end function sum_tmatrix_dmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function sum_dmatrix_tmatrix( a, b ) result( c )
    class(dmatrix), intent(in) :: a
    class(tmatrix), intent(in) :: b
    type(tmatrix)              :: c
    call c % sum_md( b, a )
  end function sum_dmatrix_tmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Reciprocal:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function reciprocal_matrix( a ) result( b )
    class(matrix), intent(in)  :: a
    type(matrix)               :: b
    call b % clone( a )
    call mproc % scale( b%size, mproc%minus_one, b%ptr )
  end function reciprocal_matrix

  function reciprocal_scalar( a ) result( b )
    class(scalar), intent(in)  :: a
    type(scalar)               :: b
    call b % clone( a )
    call mproc % scale( b%size, mproc%minus_one, b%ptr )
  end function reciprocal_scalar

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Matrix subtraction:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function subtraction_matrix_matrix( a, b ) result( c )
    class(matrix), intent(in) :: a
    class(matrix), intent(in) :: b
    type(matrix)              :: c
    call c % subtract_mm( a, b )
  end function subtraction_matrix_matrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function subtraction_tmatrix_tmatrix( a, b ) result( c )
    class(tmatrix), intent(in) :: a
    class(tmatrix), intent(in) :: b
    type(tmatrix)              :: c
    call c % subtract_mm( a, b )
  end function subtraction_tmatrix_tmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function subtraction_dmatrix_dmatrix( a, b ) result( c )
    class(dmatrix), intent(in) :: a
    class(dmatrix), intent(in) :: b
    type(dmatrix)              :: c
    call c % subtract_mm( a, b )
  end function subtraction_dmatrix_dmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function subtraction_scalar_scalar( a, b ) result( c )
    class(scalar), intent(in) :: a
    class(scalar), intent(in) :: b
    type(scalar)              :: c
    call c % subtract_mm( a, b )
  end function subtraction_scalar_scalar

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function subtraction_matrix_tmatrix( a, b ) result( c )
    class(matrix),  intent(in) :: a
    class(tmatrix), intent(in) :: b
    type(matrix)               :: c
    call c % subtract_mt( a, b )
  end function subtraction_matrix_tmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function subtraction_tmatrix_matrix( a, b ) result( c )
    class(tmatrix), intent(in) :: a
    class(matrix),  intent(in) :: b
    type(matrix)               :: c
    call c % subtract_mt( b, a )
  end function subtraction_tmatrix_matrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function subtraction_matrix_dmatrix( a, b ) result( c )
    class(matrix),  intent(in) :: a
    class(dmatrix), intent(in) :: b
    type(matrix)               :: c
    call c % subtract_md( a, b )
  end function subtraction_matrix_dmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function subtraction_dmatrix_matrix( a, b ) result( c )
    class(dmatrix), intent(in) :: a
    class(matrix),  intent(in) :: b
    type(matrix)               :: c
    call c % subtract_md( b, a )
  end function subtraction_dmatrix_matrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function subtraction_tmatrix_dmatrix( a, b ) result( c )
    class(tmatrix), intent(in) :: a
    class(dmatrix), intent(in) :: b
    type(tmatrix)              :: c
    call c % subtract_md( a, b )
  end function subtraction_tmatrix_dmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function subtraction_dmatrix_tmatrix( a, b ) result( c )
    class(dmatrix), intent(in) :: a
    class(tmatrix), intent(in) :: b
    type(tmatrix)              :: c
    call c % subtract_md( b, a )
  end function subtraction_dmatrix_tmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Hadamard product:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function hadamard_matrix_matrix( a, b ) result( c )
    class(matrix), intent(in) :: a
    class(matrix), intent(in) :: b
    type(matrix)              :: c
    if ((a%rows /= b%rows).or.(a%cols /= b%cols)) &
      call mproc % error("trying to carry out hadamard product with unconformable matrices")
    call c % hadamard( a, b )
  end function hadamard_matrix_matrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function hadamard_tmatrix_tmatrix( a, b ) result( c )
    class(tmatrix), intent(in) :: a
    class(tmatrix), intent(in) :: b
    type(tmatrix)              :: c
    if ((a%rows /= b%rows).or.(a%cols /= b%cols)) &
      call mproc % error("trying to carry out hadamard product with unconformable matrices")
    call c % hadamard( a, b )
  end function hadamard_tmatrix_tmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function hadamard_matrix_tmatrix( a, b ) result( c )
    class(matrix),  intent(in) :: a
    class(tmatrix), intent(in) :: b
    type(matrix)               :: c
    if ((a%rows /= b%cols).or.(a%cols /= b%rows)) &
      call mproc % error("trying to carry out hadamard product with unconformable matrices")
    call assign_matrix_from_tmatrix( c, b )
    call mproc % multiply( c%size, a%ptr, c%ptr )
    if (a%volatile) call a % destroy()
    c%volatile = .true.
  end function hadamard_matrix_tmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function hadamard_tmatrix_matrix( a, b ) result( c )
    class(tmatrix), intent(in) :: a
    class(matrix),  intent(in) :: b
    type(matrix)               :: c
    if ((a%cols /= b%rows).or.(a%rows /= b%cols)) &
      call mproc % error("trying to carry out hadamard product with unconformable matrices")
    call assign_matrix_from_tmatrix( c, a )
    call mproc % multiply( c%size, b%ptr, c%ptr )
    if (b%volatile) call b % destroy()
    c%volatile = .true.
  end function hadamard_tmatrix_matrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Scaling:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function product_scalar_matrix( a, b ) result( c )
    class(scalar), intent(in) :: a
    class(matrix), intent(in) :: b
    type(matrix)              :: c
    call c % scaling( a, b )
  end function product_scalar_matrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function product_matrix_scalar( a, b ) result( c )
    class(matrix), intent(in) :: a
    class(scalar), intent(in) :: b
    type(matrix)              :: c
    call c % scaling( b, a )
  end function product_matrix_scalar

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function product_scalar_tmatrix( a, b ) result( c )
    class(scalar),  intent(in) :: a
    class(tmatrix), intent(in) :: b
    type(tmatrix)              :: c
    call c % scaling( a, b )
  end function product_scalar_tmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function product_tmatrix_scalar( a, b ) result( c )
    class(tmatrix), intent(in) :: a
    class(scalar),  intent(in) :: b
    type(tmatrix)              :: c
    call c % scaling( b, a )
  end function product_tmatrix_scalar

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function product_scalar_dmatrix( a, b ) result( c )
    class(scalar),  intent(in) :: a
    class(dmatrix), intent(in) :: b
    type(dmatrix)              :: c
    call c % scaling( a, b )
  end function product_scalar_dmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function product_dmatrix_scalar( a, b ) result( c )
    class(dmatrix), intent(in) :: a
    class(scalar),  intent(in) :: b
    type(dmatrix)              :: c
    call c % scaling( b, a )
  end function product_dmatrix_scalar

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function product_scalar_scalar( a, b ) result( c )
    class(scalar), intent(in) :: a
    class(scalar), intent(in) :: b
    type(scalar)              :: c
    call c % scaling( a, b )
  end function product_scalar_scalar

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Matrix product:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function product_matrix_matrix( a, b ) result( c )
    class(matrix), intent(in) :: a
    class(matrix), intent(in) :: b
    type(matrix)              :: c
    call c % product_mm( a, b )
  end function product_matrix_matrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function product_matrix_tmatrix( a, b ) result( c )
    class(matrix),  intent(in) :: a
    class(tmatrix), intent(in) :: b
    type(matrix)               :: c
    call c % product_mt( a, b )
  end function product_matrix_tmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function product_tmatrix_matrix( a, b ) result( c )
    class(tmatrix), intent(in) :: a
    class(matrix),  intent(in) :: b
    type(matrix)               :: c
    call c % product_tm( a, b )
  end function product_tmatrix_matrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function product_tmatrix_tmatrix( a, b ) result( c )
    class(tmatrix), intent(in) :: a
    class(tmatrix), intent(in) :: b
    type(tmatrix)              :: c
    call c % product_mm( b, a )
  end function product_tmatrix_tmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function product_matrix_dmatrix( a, b ) result( c )
    class(matrix),  intent(in) :: a
    class(dmatrix), intent(in) :: b
    type(matrix)               :: c
    call c % product_md( a, b )
  end function product_matrix_dmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function product_dmatrix_matrix( a, b ) result( c )
    class(dmatrix), intent(in) :: a
    class(matrix),  intent(in) :: b
    type(matrix)               :: c
    call c % product_dm( a, b )
  end function product_dmatrix_matrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function product_tmatrix_dmatrix( a, b ) result( c )
    class(tmatrix), intent(in) :: a
    class(dmatrix), intent(in) :: b
    type(tmatrix)              :: c
    call c % product_dm( b, a )
  end function product_tmatrix_dmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function product_dmatrix_tmatrix( a, b ) result( c )
    class(dmatrix), intent(in) :: a
    class(tmatrix), intent(in) :: b
    type(tmatrix)              :: c
    call c % product_md( b, a )
  end function product_dmatrix_tmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function product_dmatrix_dmatrix( a, b ) result( c )
    class(dmatrix), intent(in) :: a
    class(dmatrix), intent(in) :: b
    type(dmatrix)              :: c
    if (a%rows /= b%rows) call mproc % error("trying to multiply unconformable matrices")
    call c % hadamard( a, b )
  end function product_dmatrix_dmatrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Exponential:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function exp_matrix( a ) result( b )
    class(matrix), intent(in) :: a
    type(matrix)              :: b
    call b % clone( a )
    call mproc % exp( b%size, b%ptr )
  end function exp_matrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Natural logarithm:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function log_matrix( a ) result( b )
    class(matrix), intent(in) :: a
    type(matrix)              :: b
    call b % clone( a )
    call mproc % log( b%size, b%ptr )
  end function log_matrix

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Power:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function pow_matrix( a, b ) result( c )
    class(matrix), intent(in)  :: a
    class(scalar), intent(in)  :: b
    type(matrix)               :: c
    call c % clone( a )
    call mproc % pow( c%size, b%ptr, c%ptr )
    if (b%volatile) call b % destroy()
  end function pow_matrix

  function pow_tmatrix( a, b ) result( c )
    class(tmatrix), intent(in)  :: a
    class(scalar),  intent(in)  :: b
    type(tmatrix)               :: c
    call c % clone( a )
    call mproc % pow( c%size, b%ptr, c%ptr )
    if (b%volatile) call b % destroy()
  end function pow_tmatrix

  function pow_scalar( a, b ) result( c )
    class(scalar), intent(in)  :: a
    class(scalar), intent(in)  :: b
    type(scalar)               :: c
    call c % clone( a )
    call mproc % pow( c%size, b%ptr, c%ptr )
    if (b%volatile) call b % destroy()
  end function pow_scalar

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Entrywise Inverse:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function inv_matrix( a ) result( b )
    class(matrix), intent(in) :: a
    type(matrix)              :: b
    call b % clone( a )
    call mproc % inv( b%size, b%ptr )
  end function inv_matrix

  function inv_tmatrix( a ) result( b )
    class(tmatrix), intent(in) :: a
    type(tmatrix)              :: b
    call b % clone( a )
    call mproc % inv( b%size, b%ptr )
  end function inv_tmatrix

  function inv_dmatrix( a ) result( b )
    class(dmatrix), intent(in) :: a
    type(dmatrix)              :: b
    call mproc % error("cannot compute the entrywise inverse of a diagonal matrix")
  end function inv_dmatrix

  function inv_scalar( a ) result( b )
    class(scalar), intent(in) :: a
    type(scalar)              :: b
    call b % clone( a )
    call mproc % inv( b%size, b%ptr )
  end function inv_scalar

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end module fortrix
