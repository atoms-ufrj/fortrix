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

module base_matrix_module

use iso_c_binding
use matrix_handler_cuda_module

type, abstract :: base_matrix

  integer(c_int) :: rows = 0              !< Number of rows
  integer(c_int) :: cols = 0              !< Number of columns
  integer(c_int) :: size = 0              !< Number of stored entries
  logical :: volatile = .true.            !< 1/0 if stored data is volatile/permanent
  type(c_ptr), pointer :: ptr => null()   !< Pointer to the stored data

  contains

    procedure, pass :: deallocate  => base_matrix_deallocate
    procedure, pass :: destroy     => base_matrix_destroy
    procedure, pass :: clone       => base_matrix_clone
    procedure, pass :: assign_from => base_matrix_assign_from
    procedure, pass :: sum_mm      => base_matrix_sum_mm
    procedure, pass :: sum_mt      => base_matrix_sum_mt
    procedure, pass :: sum_md      => base_matrix_sum_md
    procedure, pass :: subtract_mm => base_matrix_subtract_mm
    procedure, pass :: subtract_mt => base_matrix_subtract_mt
    procedure, pass :: subtract_md => base_matrix_subtract_md
    procedure, pass :: hadamard    => base_matrix_hadamard
    procedure, pass :: scaling     => base_matrix_scaling
    procedure, pass :: product_mm  => base_matrix_product_mm
    procedure, pass :: product_mt  => base_matrix_product_mt
    procedure, pass :: product_tm  => base_matrix_product_tm
    procedure, pass :: product_md  => base_matrix_product_md
    procedure, pass :: product_dm  => base_matrix_product_dm

    procedure(base_matrix_print), deferred, pass :: print

end type base_matrix

abstract interface

  subroutine base_matrix_print( x, fmt )
    import :: base_matrix
    class(base_matrix), intent(in) :: x
    character(*),       intent(in) :: fmt
  end subroutine base_matrix_print

end interface

type(matrix_handler_cuda) :: mproc

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine base_matrix_deallocate( this )
    class(base_matrix), intent(inout) :: this
    if (this%size == 0) call mproc % error( "trying to deallocate a matrix that is not allocated" )
    call mproc % deallocate( this%ptr )
    deallocate( this%ptr )
    this%rows = 0
    this%cols = 0
    this%size = 0
    this%volatile = .true.
  end subroutine base_matrix_deallocate

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine base_matrix_destroy( this )
    class(base_matrix), intent(in) :: this
    type(c_ptr), pointer :: aux
    call mproc % deallocate( this%ptr )
    aux => this%ptr
    deallocate( aux )
  end subroutine base_matrix_destroy

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine base_matrix_clone( b, a )
    class(base_matrix), intent(out) :: b
    class(base_matrix), intent(in)  :: a
    b%rows = a%rows
    b%cols = a%cols
    b%size = a%size
    if (a%volatile) then
      b%ptr => a%ptr
    else
      allocate(b%ptr)
      call mproc % allocate( b%ptr, b%size )
      call mproc % copy( b%size, a%ptr, b%ptr )
    end if
  end subroutine base_matrix_clone

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine base_matrix_assign_from( b, a )
    class(base_matrix), intent(inout) :: b
    class(base_matrix), intent(in)    :: a
    if (a%volatile) then
      if (b%volatile) then
        b%volatile = .false.
      else
        call mproc % deallocate( b%ptr )
        deallocate( b%ptr )
      end if
      b%ptr => a%ptr
    else
      if (b%volatile) then
        allocate( b%ptr )
        call mproc % allocate( b%ptr, a%size )
        b%volatile = .false.
      else if (b%size /= a%size) then
        call mproc % deallocate( b%ptr )
        call mproc % allocate( b%ptr, a%size )
      end if
      call mproc % copy( a%size, a%ptr, b%ptr )
    end if
    b%rows = a%rows
    b%cols = a%cols
    b%size = a%size
  end subroutine base_matrix_assign_from

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine base_matrix_sum_mm( c, a, b )
    class(base_matrix), intent(inout) :: c
    class(base_matrix), intent(in)    :: a, b
    if ((a%rows /= b%rows).or.(a%cols /= b%cols)) &
      call mproc % error("trying to sum unconformable matrices")
    c%rows = a%rows
    c%cols = a%cols
    c%size = a%size
    if (a%volatile) then
      c%ptr => a%ptr
      call mproc % add( c%size, b%ptr, c%ptr )
      if (b%volatile) call base_matrix_destroy( b )
    else
      if (b%volatile) then
        c%ptr => b%ptr
      else
        allocate( c%ptr )
        call mproc % allocate( c%ptr, c%size )
        call mproc % copy( c%size, b%ptr, c%ptr )
      end if
      call mproc % add( c%size, a%ptr, c%ptr )
    end if
  end subroutine base_matrix_sum_mm

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine base_matrix_sum_mt( c, a, b )
    class(base_matrix), intent(inout) :: c
    class(base_matrix), intent(in)    :: a, b
    if ((a%rows /= b%cols).or.(a%cols /= b%rows)) &
      call mproc % error("trying to sum unconformable matrices")
    c%rows = a%rows
    c%cols = a%cols
    c%size = a%size
    allocate( c%ptr )
    call mproc % allocate( c%ptr, c%size )
    call mproc % add_transp( c%rows, c%cols, a%ptr, b%ptr, c%ptr )
    if (a%volatile) call base_matrix_destroy( a )
    if (b%volatile) call base_matrix_destroy( b )
  end subroutine base_matrix_sum_mt

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine base_matrix_sum_md( c, a, b )
    class(base_matrix), intent(inout) :: c
    class(base_matrix), intent(in)    :: a, b
    if ((a%rows /= b%rows).or.(a%cols /= b%rows)) &
      call mproc % error("trying to sum unconformable matrices")
    c%rows = a%rows
    c%cols = a%cols
    c%size = a%size
    if (a%volatile) then
      c%ptr => a%ptr
    else
      allocate( c%ptr )
      call mproc % allocate( c%ptr, c%size )
      call mproc % copy( c%size, a%ptr, c%ptr )
    end if
    call mproc % add_diag( c%rows, b%ptr, c%ptr )
  end subroutine base_matrix_sum_md

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine base_matrix_subtract_mm( c, a, b )
    class(base_matrix), intent(inout) :: c
    class(base_matrix), intent(in)    :: a, b
    if ((a%rows /= b%rows).or.(a%cols /= b%cols)) &
      call mproc % error("trying to subtract unconformable matrices")
    c%rows = a%rows
    c%cols = a%cols
    c%size = a%size
    if (a%volatile) then
      c%ptr => a%ptr
    else
      allocate( c%ptr )
      call mproc % allocate( c%ptr, c%size )
      call mproc % copy( c%size, a%ptr, c%ptr )
    end if
    call mproc % subtract( c%size, b%ptr, c%ptr )
    if (b%volatile) call base_matrix_destroy( b )
  end subroutine base_matrix_subtract_mm

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine base_matrix_subtract_mt( c, a, b )
    class(base_matrix), intent(inout) :: c
    class(base_matrix), intent(in)    :: a, b
    if ((a%rows /= b%cols).or.(a%cols /= b%rows)) &
      call mproc % error("trying to subtract unconformable matrices")
    c%rows = a%rows
    c%cols = a%cols
    c%size = a%size
    allocate( c%ptr )
    call mproc % allocate( c%ptr, c%size )
    call mproc % subt_transp( c%rows, c%cols, a%ptr, b%ptr, c%ptr )
    if (a%volatile) call base_matrix_destroy( a )
    if (b%volatile) call base_matrix_destroy( b )
  end subroutine base_matrix_subtract_mt

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine base_matrix_subtract_md( c, a, b )
    class(base_matrix), intent(inout) :: c
    class(base_matrix), intent(in)    :: a, b
    if ((a%rows /= b%rows).or.(a%cols /= b%rows)) &
      call mproc % error("trying to subtract unconformable matrices")
    c%rows = a%rows
    c%cols = a%cols
    c%size = a%size
    if (a%volatile) then
      c%ptr => a%ptr
    else
      allocate( c%ptr )
      call mproc % allocate( c%ptr, c%size )
      call mproc % copy( c%size, a%ptr, c%ptr )
    end if
    call mproc % subt_diag( c%rows, b%ptr, c%ptr )
  end subroutine base_matrix_subtract_md

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine base_matrix_hadamard( c, a, b )
    class(base_matrix), intent(inout) :: c
    class(base_matrix), intent(in)    :: a, b
    c%rows = a%rows
    c%cols = a%cols
    c%size = a%size
    if (a%volatile) then
      c%ptr => a%ptr
      call mproc % multiply( c%size, b%ptr, c%ptr )
      if (b%volatile) call base_matrix_destroy( b )
    else
      if (b%volatile) then
        c%ptr => b%ptr
      else
        allocate( c%ptr )
        call mproc % allocate( c%ptr, c%size )
        call mproc % copy( c%size, b%ptr, c%ptr )
      end if
      call mproc % multiply( c%size, a%ptr, c%ptr )
    end if
  end subroutine base_matrix_hadamard

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine base_matrix_scaling( c, a, b )
    class(base_matrix), intent(inout) :: c
    class(base_matrix), intent(in)    :: a, b
    c%rows = b%rows
    c%cols = b%cols
    c%size = b%size
    if (b%volatile) then
      c%ptr => b%ptr
    else
      allocate( c%ptr )
      call mproc % allocate( c%ptr, c%size )
      call mproc % copy( b%size, b%ptr, c%ptr )
    end if
    call mproc % scale( b%size, a%ptr, c%ptr )
    if (a%volatile) call base_matrix_destroy( a )
  end subroutine base_matrix_scaling

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine base_matrix_product_mm( c, a, b )
    class(base_matrix), intent(inout) :: c
    class(base_matrix), intent(in)    :: a, b
    if (a%cols /= b%rows) call mproc % error("trying to multiply unconformable matrices")
    c%rows = a%rows
    c%cols = b%cols
    c%size = c%rows * c%cols
    allocate( c%ptr )
    call mproc % allocate( c%ptr, c%size )
    call mproc % product_mm( a%rows, a%cols, b%cols, a%ptr, b%ptr, c%ptr )
    if (a%volatile) call base_matrix_destroy( a )
    if (b%volatile) call base_matrix_destroy( b )
  end subroutine base_matrix_product_mm

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine base_matrix_product_mt( c, a, b )
    class(base_matrix), intent(inout) :: c
    class(base_matrix), intent(in)    :: a, b
    if (a%cols /= b%cols) call mproc % error("trying to multiply unconformable matrices")
    c%rows = a%rows
    c%cols = b%rows
    c%size = c%rows * c%cols
    allocate( c%ptr )
    call mproc % allocate( c%ptr, c%size )
    call mproc % product_mt( a%rows, a%cols, b%rows, a%ptr, b%ptr, c%ptr )
    if (a%volatile) call base_matrix_destroy( a )
    if (b%volatile) call base_matrix_destroy( b )
  end subroutine base_matrix_product_mt

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine base_matrix_product_tm( c, a, b )
    class(base_matrix), intent(inout) :: c
    class(base_matrix), intent(in)    :: a, b
    if (a%rows /= b%rows) call mproc % error("trying to multiply unconformable matrices")
    c%rows = a%cols
    c%cols = b%cols
    c%size = c%rows * c%cols
    allocate( c%ptr )
    call mproc % allocate( c%ptr, c%size )
    call mproc % product_tm( a%cols, a%rows, b%cols, a%ptr, b%ptr, c%ptr )
    if (a%volatile) call base_matrix_destroy( a )
    if (b%volatile) call base_matrix_destroy( b )
  end subroutine base_matrix_product_tm

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine base_matrix_product_md( c, a, b )
    class(base_matrix), intent(inout) :: c
    class(base_matrix), intent(in)    :: a, b
    if (a%cols /= b%rows) call mproc % error("trying to multiply unconformable matrices")
    if (a%rows == 1_c_int) then
      call base_matrix_hadamard( c, a, b )
    else
      c%rows = a%rows
      c%cols = a%cols
      c%size = a%size
      allocate( c%ptr )
      call mproc % allocate( c%ptr, c%size )
      call mproc % product_md( a%rows, a%cols, a%ptr, b%ptr, c%ptr )
      if (a%volatile) call base_matrix_destroy( a )
      if (b%volatile) call base_matrix_destroy( b )
    end if
  end subroutine base_matrix_product_md

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine base_matrix_product_dm( c, a, b )
    class(base_matrix), intent(inout) :: c
    class(base_matrix), intent(in)    :: a, b
    if (a%rows /= b%rows) call mproc % error("trying to multiply unconformable matrices")
    if (b%cols == 1_c_int) then
      call base_matrix_hadamard( c, a, b )
    else
      c%rows = b%rows
      c%cols = b%cols
      c%size = b%size
      allocate( c%ptr )
      call mproc % allocate( c%ptr, c%size )
      call mproc % product_dm( b%rows, b%cols, a%ptr, b%ptr, c%ptr )
      if (a%volatile) call base_matrix_destroy( a )
      if (b%volatile) call base_matrix_destroy( b )
    end if
  end subroutine base_matrix_product_dm

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end module base_matrix_module
