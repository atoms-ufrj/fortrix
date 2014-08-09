program test

use Fortrix

type(matrix)  :: A, B, C
type(dmatrix) :: D, E, F

call Fortrix_Startup()

A = new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [2,4] ) )
B = new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [4,2] ) )

D = .d.new_matrix( [real(rb) :: 1,2,3,4 ] )

do while (.true.)

  write(*,'(/,"Testing product of a temporary dmatrix with a temporary dmatrix...")')
  E = (.d.new_matrix( [real(rb) :: 1,2,3,4 ] ))*(.d.new_matrix( [real(rb) :: 1,2,3,4 ] ))
  print*
  call print_matrix( E, "F6.1" )

  write(*,'(/,"Testing product of a temporary dmatrix with a permanent dmatrix...")')
  E = (.d.new_matrix( [real(rb) :: 1,2,3,4 ] ))*D
  print*
  call print_matrix( E, "F6.1" )

  write(*,'(/,"Testing product of a permanent dmatrix with a temporary dmatrix...")')
  E = D*(.d.new_matrix( [real(rb) :: 1,2,3,4 ] ))
  print*
  call print_matrix( E, "F6.1" )

  write(*,'(/,"Testing product of a permanent dmatrix with a permanent dmatrix...")')
  F = D*E
  print*
  call print_matrix( F, "F6.1" )

  write(*,'(/,"Testing product of a permanent dmatrix with itself...")')
  F = D*D
  print*
  call print_matrix( F, "F6.1" )

  write(*,'(/,"Testing product of a temporary matrix with a temporary dmatrix...")')
  C = new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [2,4] ) )*(.d.new_matrix( [real(rb) :: 1,2,3,4 ] ))
  print*
  call print_matrix( C, "F6.1" )

  write(*,'(/,"Testing product of a temporary matrix with a permanent dmatrix...")')
  C = new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [2,4] ) )*D
  print*
  call print_matrix( C, "F6.1" )

  write(*,'(/,"Testing product of a permanent matrix with a permanent dmatrix...")')
  C = A*D
  print*
  call print_matrix( C, "F6.1" )

  write(*,'(/,"Testing product of a temporary dmatrix with a temporary matrix...")')
  C = (.d.new_matrix( [real(rb) :: 1,2,3,4 ] ))*new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [4,2] ) )
  print*
  call print_matrix( C, "F6.1" )

  write(*,'(/,"Testing product of a permanent dmatrix with a temporary matrix...")')
  C = D*new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [4,2] ) )
  print*
  call print_matrix( C, "F6.1" )

  write(*,'(/,"Testing product of a temporary dmatrix with a permanent matrix...")')
  C = (.d.new_matrix( [real(rb) :: 1,2,3,4 ] ))*B
  print*
  call print_matrix( C, "F6.1" )

  write(*,'(/,"Testing product of a permanent dmatrix with a permanent matrix...")')
  C = D*B
  print*
  call print_matrix( C, "F6.1" )

end do

call Fortrix_Shutdown()

end program test
