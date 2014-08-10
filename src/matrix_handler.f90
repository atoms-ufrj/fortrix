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

module matrix_handler_module

use iso_c_binding

implicit none

#if defined(double)
  integer, parameter :: rb = 8
#else
  integer, parameter :: rb = 4
#endif

type, abstract :: matrix_handler
  type(c_ptr) :: handle
  type(c_ptr) :: one, zero, minus_one
  contains
    procedure(matrix_handler_initialize),  deferred, pass   :: initialize
    procedure(matrix_handler_finalize),    deferred, pass   :: finalize
    procedure(matrix_handler_allocate),    deferred, nopass :: allocate
    procedure(matrix_handler_deallocate),  deferred, nopass :: deallocate
    procedure(matrix_handler_upload),      deferred, nopass :: upload
    procedure(matrix_handler_download),    deferred, nopass :: download
    procedure(matrix_handler_copy),        deferred, pass   :: copy
    procedure(matrix_handler_scale),       deferred, pass   :: scale
    procedure(matrix_handler_replicate),   deferred, nopass :: replicate
    procedure(matrix_handler_add),         deferred, pass   :: add
    procedure(matrix_handler_subtract),    deferred, pass   :: subtract
    procedure(matrix_handler_add_diag),    deferred, pass   :: add_diag
    procedure(matrix_handler_subt_diag),   deferred, pass   :: subt_diag
    procedure(matrix_handler_multiply),    deferred, nopass :: multiply
    procedure(matrix_handler_exp),         deferred, nopass :: exp
    procedure(matrix_handler_log),         deferred, nopass :: log
    procedure(matrix_handler_inv),         deferred, nopass :: inv
    procedure(matrix_handler_pow),         deferred, nopass :: pow
    procedure(matrix_handler_transpose),   deferred, pass   :: transpose
    procedure(matrix_handler_add_transp),  deferred, pass   :: add_transp
    procedure(matrix_handler_subt_transp), deferred, pass   :: subt_transp
    procedure(matrix_handler_product_mm),  deferred, pass   :: product_mm
    procedure(matrix_handler_product_mt),  deferred, pass   :: product_mt
    procedure(matrix_handler_product_md),  deferred, pass   :: product_md
    procedure(matrix_handler_product_dm),  deferred, pass   :: product_dm

    procedure, nopass :: print => matrix_handler_print
    procedure, nopass :: error => matrix_handler_error

end type matrix_handler

abstract interface

  subroutine matrix_handler_initialize( this, thread )
    import :: matrix_handler
    class(matrix_handler), intent(inout) :: this
    integer,                    intent(in)    :: thread
  end subroutine matrix_handler_initialize

  subroutine matrix_handler_finalize( this )
    import :: matrix_handler
    class(matrix_handler), intent(inout) :: this
  end subroutine matrix_handler_finalize

  subroutine matrix_handler_allocate( x, n )
    import :: c_ptr
    type(c_ptr),           intent(inout) :: x
    integer,               intent(in)    :: n
  end subroutine matrix_handler_allocate

  subroutine matrix_handler_deallocate( x )
    import :: c_ptr
    type(c_ptr),           intent(inout) :: x
  end subroutine matrix_handler_deallocate

  subroutine matrix_handler_upload( n, x, y )
    import :: rb, c_int, c_ptr
    integer(c_int), intent(in)         :: n
    real(rb),       intent(in), target :: x(*)
    type(c_ptr),    intent(in)         :: y
  end subroutine matrix_handler_upload

  subroutine matrix_handler_download( n, x, y )
    import :: rb, c_int, c_ptr
    integer(c_int), intent(in)         :: n
    type(c_ptr),    intent(in)         :: x
    real(rb),       intent(in), target :: y(*)
  end subroutine matrix_handler_download

  subroutine matrix_handler_copy( this, n, x, y )
    import :: matrix_handler, c_int, c_ptr
    class(matrix_handler), intent(inout) :: this
    integer(c_int),        intent(in) :: n
    type(c_ptr),           intent(in) :: x
    type(c_ptr),           intent(in) :: y
  end subroutine matrix_handler_copy

  subroutine matrix_handler_scale( this, n, x, y )
    import :: matrix_handler, c_int, c_ptr
    class(matrix_handler), intent(inout) :: this
    integer(c_int),        intent(in) :: n
    type(c_ptr),           intent(in) :: x
    type(c_ptr),           intent(in) :: y
  end subroutine matrix_handler_scale

  subroutine matrix_handler_replicate( n, x, y )
    import :: c_int, c_ptr
    integer(c_int),             intent(in) :: n
    type(c_ptr),                intent(in) :: x
    type(c_ptr),                intent(in) :: y
  end subroutine matrix_handler_replicate

  subroutine matrix_handler_add( this, n, x, y )
    import :: matrix_handler, c_int, c_ptr
    class(matrix_handler), intent(inout) :: this
    integer(c_int),        intent(in) :: n
    type(c_ptr),           intent(in) :: x
    type(c_ptr),           intent(in) :: y
  end subroutine matrix_handler_add

  subroutine matrix_handler_subtract( this, n, x, y )
    import :: matrix_handler, c_int, c_ptr
    class(matrix_handler), intent(inout) :: this
    integer(c_int),        intent(in) :: n
    type(c_ptr),           intent(in) :: x
    type(c_ptr),           intent(in) :: y
  end subroutine matrix_handler_subtract

  subroutine matrix_handler_add_diag( this, n, x, y )
    import :: matrix_handler, c_int, c_ptr
    class(matrix_handler), intent(inout) :: this
    integer(c_int),        intent(in) :: n
    type(c_ptr),           intent(in) :: x
    type(c_ptr),           intent(in) :: y
  end subroutine matrix_handler_add_diag

  subroutine matrix_handler_subt_diag( this, n, x, y )
    import :: matrix_handler, c_int, c_ptr
    class(matrix_handler), intent(inout) :: this
    integer(c_int),        intent(in) :: n
    type(c_ptr),           intent(in) :: x
    type(c_ptr),           intent(in) :: y
  end subroutine matrix_handler_subt_diag

  subroutine matrix_handler_multiply( n, x, y )
    import :: c_int, c_ptr
    integer(c_int), intent(in) :: n
    type(c_ptr),    intent(in) :: x
    type(c_ptr),    intent(in) :: y
  end subroutine matrix_handler_multiply

  subroutine matrix_handler_exp( n, x )
    import :: c_int, c_ptr
    integer(c_int), intent(in) :: n
    type(c_ptr),    intent(in) :: x
  end subroutine matrix_handler_exp

  subroutine matrix_handler_log( n, x )
    import :: c_int, c_ptr
    integer(c_int), intent(in) :: n
    type(c_ptr),    intent(in) :: x
  end subroutine matrix_handler_log

  subroutine matrix_handler_inv( n, x )
    import :: c_int, c_ptr
    integer(c_int), intent(in) :: n
    type(c_ptr),    intent(in) :: x
  end subroutine matrix_handler_inv

  subroutine matrix_handler_pow( n, x, y )
    import :: c_int, c_ptr
    integer(c_int), intent(in) :: n
    type(c_ptr),    intent(in) :: x, y
  end subroutine matrix_handler_pow

  subroutine matrix_handler_transpose( this, m, n, x, y )
    import :: matrix_handler, c_int, c_ptr
    class(matrix_handler), intent(inout) :: this
    integer(c_int),        intent(in)    :: m, n
    type(c_ptr),           intent(in)    :: x, y
  end subroutine matrix_handler_transpose

  subroutine matrix_handler_add_transp( this, m, n, x, y, z )
    import :: matrix_handler, c_int, c_ptr
    class(matrix_handler), intent(inout) :: this
    integer(c_int),        intent(in)    :: m, n
    type(c_ptr),           intent(in)    :: x, y, z
  end subroutine matrix_handler_add_transp

  subroutine matrix_handler_subt_transp( this, m, n, x, y, z )
    import :: matrix_handler, c_int, c_ptr
    class(matrix_handler), intent(inout) :: this
    integer(c_int),        intent(in)    :: m, n
    type(c_ptr),           intent(in)    :: x, y, z
  end subroutine matrix_handler_subt_transp

  subroutine matrix_handler_product_mm( this, m, k, n, x, y, z )
    import :: matrix_handler, c_int, c_ptr
    class(matrix_handler), intent(inout) :: this
    integer(c_int),        intent(in)    :: m, k, n
    type(c_ptr),           intent(in)    :: x, y, z
  end subroutine matrix_handler_product_mm

  subroutine matrix_handler_product_mt( this, m, k, n, x, y, z )
    import :: matrix_handler, c_int, c_ptr
    class(matrix_handler), intent(inout) :: this
    integer(c_int),        intent(in)    :: m, k, n
    type(c_ptr),           intent(in)    :: x, y, z
  end subroutine matrix_handler_product_mt

  subroutine matrix_handler_product_md( this, m, n, x, y, z )
    import :: matrix_handler, c_int, c_ptr
    class(matrix_handler), intent(inout) :: this
    integer(c_int),        intent(in)    :: m, n
    type(c_ptr),           intent(in)    :: x, y, z
  end subroutine matrix_handler_product_md

  subroutine matrix_handler_product_dm( this, m, n, x, y, z )
    import :: matrix_handler, c_int, c_ptr
    class(matrix_handler), intent(inout) :: this
    integer(c_int),        intent(in)    :: m, n
    type(c_ptr),           intent(in)    :: x, y, z
  end subroutine matrix_handler_product_dm

end interface

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine matrix_handler_print( m, n, array, fmt )
    integer(c_int), intent(in) :: m, n
    real(rb),       intent(in) :: array(*)
    character(*),   intent(in) :: fmt
    integer :: i, j
    character(3) :: C
    write(C,'(I3)') m*n
    do i = 1, m
      write(*,fmt='('//C//fmt//')') (array((j-1)*m+i),j=1,n)
    end do
  end subroutine matrix_handler_print

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine matrix_handler_error( msg )
    character(*), intent(in) :: msg
    write(*,'("Error: ",A,".")') msg
    stop
  end subroutine matrix_handler_error

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end module matrix_handler_module
