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

module matrix_handler_blas_module

!#if defined(double)
!#  define cudaRexp    cudaDexp
!#  define cudaRlog    cudaDlog
!#  define cudaRinv    cudaDinv
!#  define cudaRpow    cudaDpow
!#  define cudaRprod   cudaDprod
!#  define cudaRreplic cudaDreplic
!#  define cublasRcopy cublasDcopy
!#  define cublasRaxpy cublasDaxpy
!#  define cublasRscal cublasDscal
!#  define cublasRgbmv cublasDgbmv
!#  define cublasRgemv cublasDgemv
!#  define cublasRger  cublasDger
!#  define cublasRgemm cublasDgemm
!#  define cublasRgeam cublasDgeam
!#  define cublasRdgmm cublasDdgmm
!#else
!#  define cudaRexp    cudaSexp
!#  define cudaRlog    cudaSlog
!#  define cudaRinv    cudaSinv
!#  define cudaRpow    cudaSpow
!#  define cudaRprod   cudaSprod
!#  define cudaRreplic cudaSreplic
!#  define cublasRcopy cublasScopy
!#  define cublasRaxpy cublasSaxpy
!#  define cublasRscal cublasSscal
!#  define cublasRgbmv cublasSgbmv
!#  define cublasRgemv cublasSgemv
!#  define cublasRger  cublasSger
!#  define cublasRgemm cublasSgemm
!#  define cublasRgeam cublasSgeam
!#  define cublasRdgmm cublasSdgmm
!#endif

use iso_c_binding
use matrix_handler_module
!use cuda_interface

implicit none

integer(c_int), parameter, private :: zero = 0_c_int, &
                                      one  = 1_c_int, &
                                      crb  = int(rb,c_int)

type, extends(matrix_handler) :: matrix_handler_blas
  contains
    procedure, pass   :: initialize  => matrix_handler_blas_initialize
    procedure, pass   :: finalize    => matrix_handler_blas_finalize
    procedure, nopass :: allocate    => matrix_handler_blas_allocate
    procedure, nopass :: deallocate  => matrix_handler_blas_deallocate
    procedure, nopass :: upload      => matrix_handler_blas_upload
    procedure, nopass :: download    => matrix_handler_blas_download
    procedure, pass   :: copy        => matrix_handler_blas_copy
    procedure, pass   :: scale       => matrix_handler_blas_scale
    procedure, nopass :: replicate   => matrix_handler_blas_replicate
    procedure, pass   :: add         => matrix_handler_blas_add
    procedure, pass   :: subtract    => matrix_handler_blas_subtract
    procedure, pass   :: add_diag    => matrix_handler_blas_add_diag
    procedure, pass   :: subt_diag   => matrix_handler_blas_subt_diag
    procedure, nopass :: multiply    => matrix_handler_blas_multiply
    procedure, nopass :: exp         => matrix_handler_blas_exp
    procedure, nopass :: log         => matrix_handler_blas_log
    procedure, nopass :: inv         => matrix_handler_blas_inv
    procedure, nopass :: pow         => matrix_handler_blas_pow
    procedure, pass   :: transpose   => matrix_handler_blas_transpose
    procedure, pass   :: add_transp  => matrix_handler_blas_add_transp
    procedure, pass   :: subt_transp => matrix_handler_blas_subt_transp
    procedure, pass   :: product_mm  => matrix_handler_blas_product_mm
    procedure, pass   :: product_mt  => matrix_handler_blas_product_mt
    procedure, pass   :: product_tm  => matrix_handler_blas_product_tm
    procedure, pass   :: product_md  => matrix_handler_blas_product_md
    procedure, pass   :: product_dm  => matrix_handler_blas_product_dm
end type matrix_handler_blas

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine matrix_handler_blas_initialize( this, thread )
    class(matrix_handler_blas), intent(inout) :: this
    integer,                    intent(in)    :: thread
!    integer(c_int) :: ndevs, device, stat
!    stat = cudaGetDeviceCount(ndevs)
!    if (stat /= 0) stop "Error using CUDA."
!    write(*,'("Number of CUDA-enabled devices available is ",I2,".")') ndevs
!    device = mod(thread,ndevs)
!    stat = cudaSetDevice( device )
!    if (stat /= 0) stop "Error using CUDA."
!    write(*,'("Calculations performed in device ",I2,".")') device
!    stat = cublasCreate( this % handle )
!    if (stat /= 0) stop "Error using CUBLAS."
!    write(*,'("CUBLAS library context successfully created.")')
!    stat = cublasSetPointerMode( this % handle, CUBLAS_POINTER_MODE_DEVICE )
!    if (stat /= 0) stop "Error using CUBLAS."
!    write(*,'("CUBLAS pointer mode set to ""device"". ")')
!    call this % allocate( this % zero, 1 )
!    call this % upload( 1, [0.0_rb], this % zero )
!    call this % allocate( this % one, 1 )
!    call this % upload( 1, [1.0_rb], this % one )
!    call this % allocate( this % minus_one, 1 )
!    call this % upload( 1, [-1.0_rb], this % minus_one )
  end subroutine matrix_handler_blas_initialize

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine matrix_handler_blas_finalize( this )
    class(matrix_handler_blas), intent(inout) :: this
!    integer(c_int) :: stat
!    call this % deallocate( this % zero )
!    call this % deallocate( this % one )
!    call this % deallocate( this % minus_one )
!    stat = cublasDestroy( this % handle )
!    if (stat /= 0) stop "Error using CUBLAS."
!    write(*,'("CUBLAS library context successfully destroyed. ")')
  end subroutine matrix_handler_blas_finalize

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Allocates device memory for n entries and points x to its leading address.
  subroutine matrix_handler_blas_allocate( x, n )
    type(c_ptr),    intent(inout) :: x
    integer(c_int), intent(in)    :: n
    
!    integer(c_int) :: stat
!    stat = cudaMalloc( x, int(n*rb,c_size_t) )
!    if (stat /= 0) stop "Error using CUDA."
  end subroutine matrix_handler_blas_allocate

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Frees the device memory space pointed by x.
  subroutine matrix_handler_blas_deallocate( x )
    type(c_ptr), intent(inout) :: x
!    integer(c_int) :: stat
!    stat = cudaFree( x )
!    if (stat /= 0) stop "Error using CUDA."
  end subroutine matrix_handler_blas_deallocate

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Uploads n entries of array x to pointer y.
  subroutine matrix_handler_blas_upload( n, x, y )
    integer(c_int), intent(in)         :: n
    real(rb),       intent(in), target :: x(*)
    type(c_ptr),    intent(in)         :: y
!    integer(c_int) :: stat
!    type(c_ptr)    :: xp
!    xp = c_loc(x(1))
!    stat = cublasSetVector( n, crb, xp, one, y, one )
!    if (stat /= 0) stop "Error using CUBLAS."
  end subroutine matrix_handler_blas_upload

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Downloads n entries from pointer x to array y.
  subroutine matrix_handler_blas_download( n, x, y )
    integer(c_int), intent(in)         :: n
    type(c_ptr),    intent(in)         :: x
    real(rb),       intent(in), target :: y(*)
!    integer(c_int) :: stat
!    type(c_ptr)    :: yp
!    yp = c_loc(y(1))
!    stat = cublasGetVector( n, crb, x, one, yp, one )
!    if (stat /= 0) stop "Error using CUBLAS."
  end subroutine matrix_handler_blas_download

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Copies n entries from pointer x to pointer y.
  subroutine matrix_handler_blas_copy( this, n, x, y )
    class(matrix_handler_blas), intent(inout) :: this
    integer(c_int),             intent(in) :: n
    type(c_ptr),                intent(in) :: x
    type(c_ptr),                intent(in) :: y
!    integer(c_int) :: stat
!    stat = cublasRcopy( this % handle, n, x, one, y, one )
!    if (stat /= 0) stop "Error using CUBLAS."
  end subroutine matrix_handler_blas_copy

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Scale n entries of pointer y by the first entry of pointer x.
  subroutine matrix_handler_blas_scale( this, n, x, y )
    class(matrix_handler_blas), intent(inout) :: this
    integer(c_int),             intent(in) :: n
    type(c_ptr),                intent(in) :: x
    type(c_ptr),                intent(in) :: y
!    integer(c_int) :: stat
!    stat = cublasRscal( this % handle, n, x, y, one )
!    if (stat /= 0) stop "Error using CUBLAS."
  end subroutine matrix_handler_blas_scale

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Replicate the first entry of pointer x as the n first entries of pointer y.
  subroutine matrix_handler_blas_replicate( n, x, y )
    integer(c_int), intent(in) :: n
    type(c_ptr),    intent(in) :: x
    type(c_ptr),    intent(in) :: y
!    call cudaRreplic( n, x, y )
  end subroutine matrix_handler_blas_replicate

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Adds n entries of pointer x to the corresponding entries of pointer y.
  subroutine matrix_handler_blas_add( this, n, x, y )
    class(matrix_handler_blas), intent(inout) :: this
    integer(c_int),             intent(in) :: n
    type(c_ptr),                intent(in) :: x
    type(c_ptr),                intent(in) :: y
!    integer(c_int) :: stat
!    stat = cublasRaxpy( this % handle, n, this % one, x, one, y, one )
!    if (stat /= 0) stop "Error using CUBLAS."
  end subroutine matrix_handler_blas_add

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Subtracts n entries of pointer x from the corresponding entries of pointer y.
  subroutine matrix_handler_blas_subtract( this, n, x, y )
    class(matrix_handler_blas), intent(inout) :: this
    integer(c_int),             intent(in) :: n
    type(c_ptr),                intent(in) :: x
    type(c_ptr),                intent(in) :: y
!    integer(c_int) :: stat
!    stat = cublasRaxpy( this % handle, n, this % minus_one, x, one, y, one )
!    if (stat /= 0) stop "Error using CUBLAS."
  end subroutine matrix_handler_blas_subtract

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Adds n entries of pointer x to the diagonal of an n-by-n matrix pointed by y.
  subroutine matrix_handler_blas_add_diag( this, n, x, y )
    class(matrix_handler_blas), intent(inout) :: this
    integer(c_int),             intent(in) :: n
    type(c_ptr),                intent(in) :: x
    type(c_ptr),                intent(in) :: y
!    integer(c_int) :: stat
!    stat = cublasRaxpy( this % handle, n, this % one, x, one, y, n )
!    if (stat /= 0) stop "Error using CUBLAS."
  end subroutine matrix_handler_blas_add_diag

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Subtracts n entries of pointer x from the diagonal of an n-by-n matrix pointed by y.
  subroutine matrix_handler_blas_subt_diag( this, n, x, y )
    class(matrix_handler_blas), intent(inout) :: this
    integer(c_int),             intent(in) :: n
    type(c_ptr),                intent(in) :: x
    type(c_ptr),                intent(in) :: y
!    integer(c_int) :: stat
!    stat = cublasRaxpy( this % handle, n, this % minus_one, x, one, y, n )
!    if (stat /= 0) stop "Error using CUBLAS."
  end subroutine matrix_handler_blas_subt_diag

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Multiply n entries of pointer x to the corresponding entries of pointer y.
  subroutine matrix_handler_blas_multiply( n, x, y )
    integer(c_int), intent(in) :: n
    type(c_ptr),    intent(in) :: x
    type(c_ptr),    intent(in) :: y
!    call cudaRprod( n, x, y )
  end subroutine matrix_handler_blas_multiply

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Applies the expoential function to n entries of pointer x.
  subroutine matrix_handler_blas_exp( n, x )
    integer(c_int), intent(in) :: n
    type(c_ptr),    intent(in) :: x
!    call cudaRexp( n, x )
  end subroutine matrix_handler_blas_exp

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Applies the natural logarithm function to n entries of pointer x.
  subroutine matrix_handler_blas_log( n, x )
    integer(c_int), intent(in) :: n
    type(c_ptr),    intent(in) :: x
!    call cudaRlog( n, x )
  end subroutine matrix_handler_blas_log

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Calculates the inverse of each of n entries of pointer x.
  subroutine matrix_handler_blas_inv( n, x )
    integer(c_int), intent(in) :: n
    type(c_ptr),    intent(in) :: x
!    call cudaRinv( n, x )
  end subroutine matrix_handler_blas_inv

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Elevates n entries of pointer y to the power of the first entry of pointer x.
  subroutine matrix_handler_blas_pow( n, x, y )
    integer(c_int), intent(in) :: n
    type(c_ptr),    intent(in) :: x, y
!    call cudaRpow( n, y, x )
  end subroutine matrix_handler_blas_pow

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Transposes an m-by-n matrix pointed by x and saves the result in y.
  subroutine matrix_handler_blas_transpose( this, m, n, x, y )
    class(matrix_handler_blas), intent(inout) :: this
    integer(c_int),             intent(in)    :: m, n
    type(c_ptr),                intent(in)    :: x, y
!    integer(c_int) :: stat
!    stat = cublasRgeam( this%handle, one, one, n, m, this%one, x, m, this%zero, x, m, y, n )
!    if (stat /= 0) stop "Error using CUBLAS."
  end subroutine matrix_handler_blas_transpose

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Adds the transpose of an n-by-m matrix y to the m-by-n matrix x and saves the result in z.
  subroutine matrix_handler_blas_add_transp( this, m, n, x, y, z )
    class(matrix_handler_blas), intent(inout) :: this
    integer(c_int),             intent(in)    :: m, n
    type(c_ptr),                intent(in)    :: x, y, z
!    integer(c_int) :: stat
!    stat = cublasRgeam( this%handle, zero, one, m, n, this%one, x, m, this%one, y, n, z, m )
!    if (stat /= 0) stop "Error using CUBLAS."
  end subroutine matrix_handler_blas_add_transp

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Subtracts the transpose of an n-by-m matrix y from the m-by-n matrix x and saves the result in z.
  subroutine matrix_handler_blas_subt_transp( this, m, n, x, y, z )
    class(matrix_handler_blas), intent(inout) :: this
    integer(c_int),             intent(in)    :: m, n
    type(c_ptr),                intent(in)    :: x, y, z
!    integer(c_int) :: stat
!    stat = cublasRgeam( this%handle, zero, one, m, n, this%one, x, m, this%minus_one, y, n, z, m )
!    if (stat /= 0) stop "Error using CUBLAS."
  end subroutine matrix_handler_blas_subt_transp

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ! Performs the product z = x*y, where x is an m-by-k matrix and y is a k-by-n matrix.
  subroutine matrix_handler_blas_product_mm( this, m, k, n, x, y, z )
    class(matrix_handler_blas), intent(inout) :: this
    integer(c_int),             intent(in)    :: m, k, n
    type(c_ptr),                intent(in)    :: x, y, z
!    integer(c_int) :: stat
!    stat = cublasRgemm( this%handle, zero, zero, m, n, k, this%one, x, m, y, k, this%zero, z, m )
!    if (stat /= 0) stop "Error using CUBLAS."
  end subroutine matrix_handler_blas_product_mm

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ! Performs the product z = x*transp(y), where x is an m-by-k matrix and y is a n-by-k matrix.
  subroutine matrix_handler_blas_product_mt( this, m, k, n, x, y, z )
    class(matrix_handler_blas), intent(inout) :: this
    integer(c_int),             intent(in)    :: m, k, n
    type(c_ptr),                intent(in)    :: x, y, z
!    integer(c_int) :: stat
!    stat = cublasRgemm( this%handle, zero, one, m, n, k, this%one, x, m, y, n, this%zero, z, m )
!    if (stat /= 0) stop "Error using CUBLAS."
  end subroutine matrix_handler_blas_product_mt

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ! Performs the product z = transp(x)*y, where x is an k-by-m matrix and y is a k-by-n matrix.
  subroutine matrix_handler_blas_product_tm( this, m, k, n, x, y, z )
    class(matrix_handler_blas), intent(inout) :: this
    integer(c_int),             intent(in)    :: m, k, n
    type(c_ptr),                intent(in)    :: x, y, z
!    integer(c_int) :: stat
!    stat = cublasRgemm( this%handle, one, zero, m, n, k, this%one, x, k, y, k, this%zero, z, m )
!    if (stat /= 0) stop "Error using CUBLAS."
  end subroutine matrix_handler_blas_product_tm

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Performs the product z = x*diag(y), where x is an m-by-n matrix and y is an n-dimensional vector.
  subroutine matrix_handler_blas_product_md( this, m, n, x, y, z )
    class(matrix_handler_blas), intent(inout) :: this
    integer(c_int),             intent(in)    :: m, n
    type(c_ptr),                intent(in)    :: x, y, z
!    integer(c_int) :: stat
!    stat = cublasRdgmm( this%handle, one, m, n, x, m, y, one, z, m )
!    if (stat /= 0) stop "Error using CUBLAS."
  end subroutine matrix_handler_blas_product_md

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Performs the product z = diag(x)*y, where x is an m-dimensional vector and y is an m-by-n matrix.
  subroutine matrix_handler_blas_product_dm( this, m, n, x, y, z )
    class(matrix_handler_blas), intent(inout) :: this
    integer(c_int),             intent(in)    :: m, n
    type(c_ptr),                intent(in)    :: x, y, z
!    integer(c_int) :: stat
!    stat = cublasRdgmm( this%handle, zero, m, n, y, m, x, one, z, m )
!    if (stat /= 0) stop "Error using CUBLAS."
  end subroutine matrix_handler_blas_product_dm

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end module matrix_handler_blas_module

