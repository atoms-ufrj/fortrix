type(matrix) :: A, B, Bt, C

call Fortrix_Startup()

A = new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [2,4] ) )
B = new_matrix( reshape( [real(rb) :: 2,3,4,5,6,7,8,9 ], [2,4] ) )
Bt = .t.B

do while (.true.)

  write(*,'(60("="))')

  print*, "A = "
  call print_matrix( A, "F6.1" )

  print*, "B = "
  call print_matrix( B, "F6.1" )

  print*, "Bt = "
  call print_matrix( Bt, "F6.1" )

  write(*,'(/,"Testing '//trim(operation)//' of a temporary matrix with a temporary matrix...")')
  C = new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [2,4] ) ) op &
      new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [2,4] ) )
  print*	
  call print_matrix( C, "F6.1" )

  write(*,'(/,"Testing '//trim(operation)//' of a temporary matrix with a permanent matrix...")')
  C = new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [2,4] ) ) op A
  print*	
  call print_matrix( C, "F6.1" )

  write(*,'(/,"Testing '//trim(operation)//' of a permanent matrix with a temporary matrix...")')
  C = A op new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [2,4] ) )
  print*	
  call print_matrix( C, "F6.1" )

  write(*,'(/,"Testing '//trim(operation)//' of a permanent matrix with a permanent matrix...")')
  C = A op B
  print*
  call print_matrix( C, "F6.1" )

  write(*,'(/,"Testing '//trim(operation)//' of a normal matrix with a transposed matrix...")')
  C = A op (.t.Bt)
  print*
  call print_matrix( C, "F6.1" )

  write(*,'(/,"Testing '//trim(operation)//' of a transposed matrix with a normal matrix...")')
  C = (.t.A) op Bt
  print*
  call print_matrix( C, "F6.1" )

end do

call Fortrix_Shutdown()
