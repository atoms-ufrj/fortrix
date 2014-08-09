program test

use Fortrix

type(matrix) :: A, B, C

call Fortrix_Startup()

A = new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [2,4] ) )
B = new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [4,2] ) )

do while (.true.)

  write(*,'(/,"Testing product of a temporary matrix with a temporary matrix...")')
  C = new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [2,4] ) ) .o. &
      new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [2,4] ) )
  print*	
  call print_matrix( C, "F6.1" )

  write(*,'(/,"Testing product of a temporary matrix with a permanent matrix...")')
  C = new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [2,4] ) ) .o. A
  print*	
  call print_matrix( C, "F6.1" )

  write(*,'(/,"Testing product of a permanent matrix with a temporary matrix...")')
  C = B .o. new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [4,2] ) )
  print*	
  call print_matrix( C, "F6.1" )

  write(*,'(/,"Testing product of a permanent matrix with a permanent matrix...")')
  C = A .o. A
  print*	
  call print_matrix( C, "F6.1" )

  print*, "A = "
  call print_matrix( A, "F6.1" )

  print*, "B = "
  call print_matrix( B, "F6.1" )

  write(*,'(/,"Testing product of a normal matrix with a transposed matrix...")')
  C = A .o. (.t.B)
  print*
  call print_matrix( C, "F6.1" )

  write(*,'(/,"Testing product of a transposed matrix with a normal matrix...")')
  C = (.t.A) .o. B
  print*
  call print_matrix( C, "F6.1" )

end do

call Fortrix_Shutdown()

end program test
