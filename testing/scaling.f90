program test

use Fortrix

type(matrix) :: A, B
type(scalar) :: S
type(dmatrix) :: D, E

call Fortrix_Startup()

call print_matrix( new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8,9 ], [3,3] ) ), "F6.1" )

S = new_scalar(2.0_rb)
B = new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8,9 ], [3,3] ) )
D = .d.new_matrix( [real(rb) :: 1,2,3,4 ] )

do while (.true.)

  write(*,'(/,"Testing scaling of a temporary matrix with a temporary scalar...")')
  A = new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8,9 ], [3,3] ) ) * new_scalar(2.0_rb)
  print*	
  call print_matrix( A, "F6.1" )

  write(*,'(/,"Testing scaling of a temporary matrix with a permanent scalar...")')
  A = new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8,9 ], [3,3] ) ) * S
  print*	
  call print_matrix( A, "F6.1" )

  write(*,'(/,"Testing scaling of a permanent matrix with a temporary scalar...")')
  A = B * new_scalar(2.0_rb)
  print*	
  call print_matrix( A, "F6.1" )

  write(*,'(/,"Testing scaling of a permanent matrix with a permanent scalar...")')
  A = B * S
  print*	
  call print_matrix( A, "F6.1" )

  write(*,'(/,"Testing switched scaling of a permanent matrix with a permanent scalar...")')
  A = S * B
  print*	
  call print_matrix( A, "F6.1" )

  write(*,'(/,"Testing scaling of a diagonal matrix with a temporary scalar...")')
  E = D * new_scalar(2.0_rb)
  print*	
  call print_matrix( E, "F6.1" )

end do

call Fortrix_Shutdown()

end program test
