type(matrix)  :: A, B, Bt, C

logical       :: exec, loop, output
character(10) :: arg
loop = iargc() > 0
if (loop) then
  call getarg( 1, arg )
  loop = trim(arg) == "loop"
end if
output = .not.loop
exec = .true.

call Fortrix_Startup()

A = new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [2,4] ) )
B = new_matrix( reshape( [real(rb) :: 2,3,4,5,6,7,8,9 ], [2,4] ) )
Bt = .t.B

do while (exec)

  if (output) write(*,'(60("="))')

  if (output) print*, "A = "
  if (output) call print_matrix( A, "F6.1" )

  if (output) print*, "B = "
  if (output) call print_matrix( B, "F6.1" )

  if (output) print*, "Bt = "
  if (output) call print_matrix( Bt, "F6.1" )

  if (output) write(*,'(/,"Testing '//trim(operation)//' of a temporary matrix with a temporary matrix...",/)')
  C = new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [2,4] ) ) op &
      new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [2,4] ) )
  if (output) call print_matrix( C, "F6.1" )

  if (output) write(*,'(/,"Testing '//trim(operation)//' of a temporary matrix with a permanent matrix...",/)')
  C = new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [2,4] ) ) op A
  if (output) call print_matrix( C, "F6.1" )

  if (output) write(*,'(/,"Testing '//trim(operation)//' of a permanent matrix with a temporary matrix...",/)')
  C = A op new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [2,4] ) )
  if (output) call print_matrix( C, "F6.1" )

  if (output) write(*,'(/,"Testing '//trim(operation)//' of a permanent matrix with a permanent matrix...",/)')
  C = A op B
  if (output) call print_matrix( C, "F6.1" )

  if (output) write(*,'(/,"Testing '//trim(operation)//' of a normal matrix with a transposed matrix...",/)')
  C = A op (.t.Bt)
  if (output) call print_matrix( C, "F6.1" )

  if (output) write(*,'(/,"Testing '//trim(operation)//' of a transposed matrix with a normal matrix...",/)')
  C = (.t.A) op Bt
  if (output) call print_matrix( C, "F6.1" )

  exec = loop

end do

call Fortrix_Shutdown()
