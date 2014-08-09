program test

use Fortrix

type(matrix) :: A, B
type(scalar) :: S, T
type(dmatrix) :: D, E

call Fortrix_Startup()

call print_matrix( new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8,9 ], [3,3] ) ), "F6.1" )

write(*,'(/,"Testing first matrix assigment from a temporary matrix...")')
B = new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8,9 ], [3,3] ) )
print*
call print_matrix( B, "F6.1" )

write(*,'(/,"Testing second matrix assigment from a temporary matrix...")')
B = new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8,9 ], [3,3] ) )
print*
call print_matrix( B, "F6.1" )

do while (.true.)

  write(*,'(/,"Testing first matrix assigment from a temporary matrix...")')
  A = new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8,9 ], [3,3] ) )
  print*	
  call print_matrix( A, "F6.1" )

  write(*,'(/,"Testing second matrix assigment from a temporary matrix of equal size...")')
  A = new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8,9 ], [3,3] ) )
  print*	
  call print_matrix( A, "F6.1" )

  write(*,'(/,"Testing second matrix assigment from a temporary matrix of different size...")')
  A = new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [2,4] ) )
  print*	
  call print_matrix( A, "F6.1" )

  call A%deallocate()

  write(*,'(/,"Testing first matrix assigment from a permanent matrix...")')
  A = B
  print*
  call print_matrix( A, "F6.1" )

  write(*,'(/,"Testing second matrix assigment from a permanent matrix...")')
  A = B
  print*
  call print_matrix( A, "F6.1" )

  write(*,'(/,"Testing matrix assigment from itself...")')
  A = A
  print*
  call print_matrix( A, "F6.1" )

  call A%deallocate()

  write(*,'(/,"Matrix to be transposed:")')
  print*
  call print_matrix( new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [2,4] ) ), "F6.1" )

  write(*,'(/,"Testing first matrix assigment from a transposed matrix...")')
  A = .t.new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [2,4] ) )
  print*
  call print_matrix( A, "F6.1" )

  write(*,'(/,"Testing second matrix assigment from a transposed matrix...")')
  A = .t.new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [2,4] ) )
  print*
  call print_matrix( A, "F6.1" )

  call A%deallocate()

  write(*,'(/,"Testing first matrix assigment from a diagonal matrix...")')
  D = .d.new_matrix([real(rb) :: 1,2,3,4 ])
  print*
  call print_matrix( D, "F6.1" )

  write(*,'(/,"Testing second matrix assigment from a diagonal matrix...")')
  D = .d.new_matrix([real(rb) :: 1,2,3,4 ])
  print*
  call print_matrix( D, "F6.1" )

  write(*,'(/,"Testing matrix assigment from a permanent diagonal matrix...")')
  E = D
  print*
  call print_matrix( E, "F6.1" )

  call D%deallocate()
  call E%deallocate()

!  write(*,'(/,"Testing first matrix assigment from a temporary scalar...")')
!  A = new_scalar(1.0_rb)
!  print*
!  call print_matrix( A, "F6.1" )

!  write(*,'(/,"Testing second matrix assigment from a temporary scalar...")')
!  A = new_scalar(2.0_rb)
!  print*
!  call print_matrix( A, "F6.1" )

  write(*,'(/,"Testing scalar assigment from a temporary scalar...")')
  S = new_scalar(2.0_rb)
  print*
  call print_matrix( S, "F6.1" )

  write(*,'(/,"Testing scalar assigment from a permanent scalar...")')
  T = S
  print*
  call print_matrix( T, "F6.1" )

!  write(*,'(/,"Testing scalar assigment from a permanent matrix...")')
!  T = A
!  print*
!  call print_matrix( T, "F6.1" )

!  call A % deallocate()

end do

call Fortrix_Shutdown()

end program test
