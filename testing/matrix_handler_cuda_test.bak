program matrix_handler_cuda_test

#if defined(double)
#  define rb 8
#else
#  define rb 4
#endif

use iso_c_binding
use matrix_handler_cuda_module

implicit none

type(matrix_handler_cuda) :: handle
type(c_ptr) :: x, y, z
integer(c_int), parameter :: n = 10
real(rb) :: x_array(10) = real([2,3,4,5,6,7,8,9,0,1],rb)

call handle % initialize( thread = 1 )
call handle % allocate( x, n )
call handle % upload( 10, x_array, x )
call print_matrix( 5, 2, x )

call handle % allocate( y, n )
call handle % copy( 10, x, y )
call print_matrix( 5, 2, y )

call handle % scale( 10, x, y )
call print_matrix( 5, 2, y )

call handle % replicate( 10, x, y )
call print_matrix( 5, 2, y )

call handle % add( 10, x, y )
call print_matrix( 5, 2, y )

call handle % multiply( 10, x, y )
call print_matrix( 5, 2, y )

call handle % copy( 10, x, y )
call print_matrix( 5, 2, y )

call handle % allocate( z, 25 )
call handle % product_mt( 5, 2, 5, x, y, z )
call print_matrix( 5, 5, z )

call handle % transpose( 5, 2, x, y )
call print_matrix( 2, 5, y )

call handle % product_mm( 5, 2, 5, x, y, z )
call print_matrix( 5, 5, z )

call print_matrix( 2, 1, y )

call handle % product_md( 5, 2, x, y, z )
call print_matrix( 5, 2, z )

call print_matrix( 5, 1, y )

call handle % product_dm( 5, 2, y, x, z )
call print_matrix( 5, 2, z )

call handle % deallocate( x )
call handle % deallocate( y )

call handle % finalize()

contains

  subroutine print_matrix( m, n, x )
    integer,     intent(in) :: m, n
    type(c_ptr), intent(in) :: x
    integer :: i
    real(8) :: x_array(m,n)
    call handle % download( m*n, x, x_array )
    print*
    do i = 1, m
      print*, x_array(i,:)
    end do
  end subroutine print_matrix

end program matrix_handler_cuda_test
