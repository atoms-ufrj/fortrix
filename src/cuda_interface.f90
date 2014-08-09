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

module cuda_interface

use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_size_t

integer(c_int), parameter :: CUBLAS_STATUS_SUCCESS          = 0,  &
                             CUBLAS_STATUS_NOT_INITIALIZED  = 1,  &
                             CUBLAS_STATUS_ALLOC_FAILED     = 3,  &
                             CUBLAS_STATUS_INVALID_VALUE    = 7,  &
                             CUBLAS_STATUS_ARCH_MISMATCH    = 8,  &
                             CUBLAS_STATUS_MAPPING_ERROR    = 11, &
                             CUBLAS_STATUS_EXECUTION_FAILED = 13, &
                             CUBLAS_STATUS_INTERNAL_ERROR   = 14

integer(c_int), parameter :: CUBLAS_FILL_MODE_LOWER = 0, &
                             CUBLAS_FILL_MODE_UPPER = 1

integer(c_int), parameter :: CUBLAS_DIAG_NON_UNIT = 0, &
                             CUBLAS_DIAG_UNIT     = 1

integer(c_int), parameter :: CUBLAS_SIDE_LEFT  = 0, &
                             CUBLAS_SIDE_RIGHT = 1

integer(c_int), parameter :: CUBLAS_OP_N = 0, &
                             CUBLAS_OP_T = 1, &
                             CUBLAS_OP_C = 2

integer(c_int), parameter :: CUBLAS_POINTER_MODE_HOST   = 0, &
                             CUBLAS_POINTER_MODE_DEVICE = 1

integer(c_int), parameter :: CUBLAS_ATOMICS_NOT_ALLOWED   = 0, &
                             CUBLAS_ATOMICS_ALLOWED       = 1

interface

  ! CUDA RUNTIME FUNCTIONS:

  function cudaGetDeviceCount( deviceCount ) result( flag ) &
    bind(C,name='cudaGetDeviceCount')
    import :: c_int
    integer(c_int) :: deviceCount
    integer(c_int) :: flag
  end function

  function cudaGetDevice( device ) result( flag ) &
    bind(C,name='cudaGetDevice')
    import :: c_int
    integer(c_int) :: device
    integer(c_int) :: flag
  end function

  function cudaSetDevice( device ) result( flag ) &
    bind(C,name='cudaSetDevice')
    import :: c_int
    integer(c_int), value :: device
    integer(c_int)        :: flag
  end function

  function cudaFree( devPtr ) result( flag ) &
    bind(C,name='cudaFree')
    import :: c_ptr, c_int
    type(c_ptr), value :: devPtr
    integer(c_int)     :: flag
  end function

  function cudaMalloc( devPtr, size ) result( flag ) &
    bind(C,name='cudaMalloc')
    import :: c_ptr, c_size_t, c_int
    type(c_ptr)              :: devPtr
    integer(c_size_t), value :: size
    integer(c_int)           :: flag
  end function

  ! UDER-DEFINED MATRIX OPERATION FUNCTIONS:

  subroutine cudaSexp( n, a ) bind(C,name='cudaSexp')
    import :: c_ptr, c_int
    integer(c_int), value :: n
    type(c_ptr),    value :: a
  end subroutine

  subroutine cudaDexp( n, a ) bind(C,name='cudaDexp')
    import :: c_ptr, c_int
    integer(c_int), value :: n
    type(c_ptr),    value :: a
  end subroutine

  subroutine cudaSlog( n, a ) bind(C,name='cudaSlog')
    import :: c_ptr, c_int
    integer(c_int), value :: n
    type(c_ptr),    value :: a
  end subroutine

  subroutine cudaDlog( n, a ) bind(C,name='cudaDlog')
    import :: c_ptr, c_int
    integer(c_int), value :: n
    type(c_ptr),    value :: a
  end subroutine

  subroutine cudaSinv( n, a ) bind(C,name='cudaSinv')
    import :: c_ptr, c_int
    integer(c_int), value :: n
    type(c_ptr),    value :: a
  end subroutine

  subroutine cudaDinv( n, a ) bind(C,name='cudaDinv')
    import :: c_ptr, c_int
    integer(c_int), value :: n
    type(c_ptr),    value :: a
  end subroutine

  subroutine cudaSpow( n, a, b ) bind(C,name='cudaSpow')
    import :: c_ptr, c_int
    integer(c_int), value :: n
    type(c_ptr),    value :: a, b
  end subroutine

  subroutine cudaDpow( n, a, b ) bind(C,name='cudaDpow')
    import :: c_ptr, c_int
    integer(c_int), value :: n
    type(c_ptr),    value :: a, b
  end subroutine

  subroutine cudaSprod( n, a, b) bind(C,name='cudaSprod')
    import :: c_ptr, c_int
    integer(c_int), value :: n
    type(c_ptr),    value :: a, b
  end subroutine

  subroutine cudaDprod( n, a, b) bind(C,name='cudaDprod')
    import :: c_ptr, c_int
    integer(c_int), value :: n
    type(c_ptr),    value :: a, b
  end subroutine

  subroutine cudaSreplic( n, a, b ) bind(C,name='cudaSreplic')
    import :: c_ptr, c_int
    integer(c_int), value :: n
    type(c_ptr),    value :: a, b
  end subroutine

  subroutine cudaDreplic( n, a, b ) bind(C,name='cudaDreplic')
    import :: c_ptr, c_int
    integer(c_int), value :: n
    type(c_ptr),    value :: a, b
  end subroutine

  ! CUBLAS HELPER FUNCTIONS:

  function cublasCreate( handle ) result( stat ) &
    bind(C,name='cublasCreate_f')
    import :: c_ptr, c_int
    type(c_ptr)    :: handle
    integer(c_int) :: stat
  end function

  function cublasDestroy( handle ) result( stat ) &
    bind(C,name='cublasDestroy_f')
    import :: c_ptr, c_int
    type(c_ptr), value :: handle
    integer(c_int)     :: stat
  end function

  function cublasGetVersion( handle, version ) result( stat ) &
    bind(C,name='cublasGetVersion_f')
    import :: c_ptr, c_int
    type(c_ptr), value :: handle
    integer(c_int)     :: version
    integer(c_int)     :: stat
  end function

  function cublasSetPointerMode( handle, mode ) result( stat ) &
    bind(C,name='cublasSetPointerMode_f')
    import :: c_ptr, c_int
    type(c_ptr),    value :: handle
    integer(c_int), value :: mode
    integer(c_int)        :: stat
  end function

  function cublasSetVector( n, elemSize, x, incx, y, incy ) result( stat ) &
    bind(C,name='cublasSetVector')
    import :: c_ptr, c_int
    integer(c_int), value :: n, elemSize, incx, incy
    type(c_ptr),    value :: x, y
    integer(c_int)        :: stat
  end function

  function cublasGetVector( n, elemSize, x, incx, y, incy ) result( stat ) &
    bind(C,name='cublasGetVector')
    import :: c_ptr, c_int
    integer(c_int), value :: n, elemSize, incx, incy
    type(c_ptr),    value :: x, y
    integer(c_int)        :: stat
  end function

  ! LEVEL-1 CUBLAS FUNCTIONS:

  function cublasScopy( handle, n, x, incx, y, incy ) result( stat ) &
    bind(C,name='cublasScopy_f')
    import :: c_ptr, c_int
    type(c_ptr),    value :: handle
    integer(c_int), value :: n, incx, incy
    type(c_ptr),    value :: x, y
    integer(c_int)        :: stat
  end function

  function cublasDcopy( handle, n, x, incx, y, incy ) result( stat ) &
    bind(C,name='cublasDcopy_f')
    import :: c_ptr, c_int
    type(c_ptr),    value :: handle
    integer(c_int), value :: n, incx, incy
    type(c_ptr),    value :: x, y
    integer(c_int)        :: stat
  end function

  function cublasSaxpy( handle, n, alpha, x, incx, y, incy ) result( stat ) &
    bind(C,name='cublasSaxpy_f')
    import :: c_ptr, c_int
    type(c_ptr),    value :: handle
    integer(c_int), value :: n, incx, incy
    type(c_ptr),    value :: alpha, x, y
    integer(c_int)        :: stat
  end function

  function cublasDaxpy( handle, n, alpha, x, incx, y, incy ) result( stat ) &
    bind(C,name='cublasDaxpy_f')
    import :: c_ptr, c_int
    type(c_ptr),    value :: handle
    integer(c_int), value :: n, incx, incy
    type(c_ptr),    value :: alpha, x, y
    integer(c_int)        :: stat
  end function

  function cublasSscal( handle, n, alpha, x, incx ) result( stat ) &
    bind(C,name='cublasSscal_f')
    import :: c_ptr, c_int
    type(c_ptr),    value :: handle
    integer(c_int), value :: n, incx
    type(c_ptr),    value :: alpha, x
    integer(c_int)        :: stat
  end function

  function cublasDscal( handle, n, alpha, x, incx ) result( stat ) &
    bind(C,name='cublasDscal_f')
    import :: c_ptr, c_int
    type(c_ptr),    value :: handle
    integer(c_int), value :: n, incx
    type(c_ptr),    value :: alpha, x
    integer(c_int)        :: stat
  end function

  ! LEVEL-2 CUBLAS FUNCTIONS:

  function cublasSgbmv( handle, trans, m, n, kl, ku, alpha, A, lda, &
                        x, incx, beta, y, incy ) result( stat ) &
    bind(C,name='cublasSgbmv_f')
    import :: c_ptr, c_int
    type(c_ptr),    value :: handle
    integer(c_int), value :: trans, m, n, kl, ku, lda, incx, incy
    type(c_ptr),    value :: alpha, A, x, beta, y
    integer(c_int)        :: stat
  end function

  function cublasDgbmv( handle, trans, m, n, kl, ku, alpha, A, lda, &
                        x, incx, beta, y, incy ) result( stat ) &
    bind(C,name='cublasDgbmv_f')
    import :: c_ptr, c_int
    type(c_ptr),    value :: handle
    integer(c_int), value :: trans, m, n, kl, ku, lda, incx, incy
    type(c_ptr),    value :: alpha, A, x, beta, y
    integer(c_int)        :: stat
  end function

  function cublasSgemv(handle, trans, m, n, alpha, A, lda, x, incx, &
                       beta, y, incy ) result( stat ) &
    bind(C,name='cublasSgemv_f')
    import :: c_ptr, c_int
    type(c_ptr),    value :: handle
    integer(c_int), value :: trans, m, n, lda, incx, incy
    type(c_ptr),    value :: alpha, A, x, beta, y
    integer(c_int)        :: stat
  end function

  function cublasDgemv( handle, trans, m, n, alpha, A, lda, x, incx, &
                        beta, y, incy ) result( stat ) &
    bind(C,name='cublasDgemv_f')
    import :: c_ptr, c_int
    type(c_ptr),    value :: handle
    integer(c_int), value :: trans, m, n, lda, incx, incy
    type(c_ptr),    value :: alpha, A, x, beta, y
    integer(c_int)        :: stat
  end function

  function cublasSger( handle, m, n, alpha, x, incx, &
                       y, incy, A, lda ) result( stat ) &
    bind(C,name='cublasSger_f')
    import :: c_ptr, c_int
    type(c_ptr),    value :: handle
    integer(c_int), value :: m, n, incx, incy, lda
    type(c_ptr),    value :: alpha, x, y, A
    integer(c_int)        :: stat
  end function

  function cublasDger( handle, m, n, alpha, x, incx, &
                       y, incy, A, lda ) result( stat ) &
    bind(C,name='cublasDger_f')
    import :: c_ptr, c_int
    type(c_ptr),    value :: handle
    integer(c_int), value :: m, n, incx, incy, lda
    type(c_ptr),    value :: alpha, x, y, A
    integer(c_int)        :: stat
  end function

  ! LEVEL-3 CUBLAS FUNCTIONS:

  function cublasSgemm( handle, transa, transb, m, n, k, alpha, &
                        A, lda, B, ldb, beta, C, ldc ) result( stat ) &
    bind(C,name='cublasSgemm_f')
    import :: c_ptr, c_int
    type(c_ptr),    value :: handle
    integer(c_int), value :: transa, transb, m, n, k, lda, ldb, ldc
    type(c_ptr),    value :: alpha, A, B, beta, C
    integer(c_int)        :: stat
  end function

  function cublasDgemm( handle, transa, transb, m, n, k, alpha, &
                        A, lda, B, ldb, beta, C, ldc ) result( stat ) &
    bind(C,name='cublasDgemm_f')
    import :: c_ptr, c_int
    type(c_ptr),    value :: handle
    integer(c_int), value :: transa, transb, m, n, k, lda, ldb, ldc
    type(c_ptr),    value :: alpha, A, B, beta, C
    integer(c_int)        :: stat
  end function

  ! CUBLAS BLAS-LIKE EXTENSIONS:

  function cublasSgeam( handle, transa, transb, m, n, alpha, A, lda, &
                        beta, B, ldb, C, ldc ) result ( stat )       &
    bind(C,name='cublasSgeam_f')
    import :: c_ptr, c_int
    type(c_ptr),    value :: handle
    integer(c_int), value :: transa, transb, m, n, lda, ldb, ldc
    type(c_ptr),    value :: alpha, beta, A, B, C
    integer(c_int)        :: stat
  end function

  function cublasDgeam( handle, transa, transb, m, n, alpha, A, lda, &
                        beta, B, ldb, C, ldc ) result ( stat )       &
    bind(C,name='cublasDgeam_f')
    import :: c_ptr, c_int
    type(c_ptr),    value :: handle
    integer(c_int), value :: transa, transb, m, n, lda, ldb, ldc
    type(c_ptr),    value :: alpha, beta, A, B, C
    integer(c_int)        :: stat
  end function

  function cublasSdgmm( handle, mode, m, n, A, lda,       &
                        x, incx, C, ldc ) result ( stat ) &
    bind(C,name='cublasSdgmm_f')
    import :: c_ptr, c_int
    type(c_ptr),    value :: handle
    integer(c_int), value :: mode, m, n, lda, incx, ldc
    type(c_ptr),    value :: A, x, C
    integer(c_int)        :: stat
  end function

  function cublasDdgmm( handle, mode, m, n, A, lda,       &
                        x, incx, C, ldc ) result ( stat ) &
    bind(C,name='cublasDdgmm_f')
    import :: c_ptr, c_int
    type(c_ptr),    value :: handle
    integer(c_int), value :: mode, m, n, lda, incx, ldc
    type(c_ptr),    value :: A, x, C
    integer(c_int)        :: stat
  end function

end interface

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Error handling for CUBLAS functions:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine cublas( stat )
    integer(c_int) :: stat
    character(21) :: error
    if (stat /= 0_c_int) then
      select case (stat)
        case (1);  error = "not initialized"
        case (3);  error = "allocation failed"
        case (7);  error = "invalid value"
        case (8);  error = "architecture mismatch"
        case (11); error = "mapping error"
        case (13); error = "execution failed"
        case (14); error = "internal error"
      end select
      write(*,'("Error using CUBLAS: ",A,".")') trim(error)
      stop
    end if
  end subroutine cublas

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Error handling for CUDA API functions:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine cuda( stat )
    integer(c_int) :: stat
    if (stat /= 0) then
      write(*,'("Error using CUDA.")')
      stop
    end if
    ! Utilizar função:
    !char *cudaGetErrorString (cudaError_t error)
  end subroutine cuda

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end module cuda_interface
