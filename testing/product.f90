program test

use Fortrix

type(matrix) :: A, B, C

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
B = new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [4,2] ) )

do while (exec)

  if (output) then
    write(*,'(/, "A = ")')
    call print_matrix( A, "F6.1" )
    write(*,'(/, "B = ")')
    call print_matrix( B, "F6.1" )
  end if

  if (output) write(*,'(/,"Testing product of a temporary matrix with a temporary matrix...",/)')
  C = new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [2,4] ) ) * &
      new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [4,2] ) )
  if (output) call print_matrix( C, "F6.1" )

  if (output) write(*,'(/,"Testing product of a temporary matrix with a permanent matrix...",/)')
  C = new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [2,4] ) ) * B
  if (output) call print_matrix( C, "F6.1" )

  if (output) write(*,'(/,"Testing product of a permanent matrix with a temporary matrix...",/)')
  C = A * new_matrix( reshape( [real(rb) :: 1,2,3,4,5,6,7,8 ], [4,2] ) )
  if (output) call print_matrix( C, "F6.1" )

  if (output) write(*,'(/,"Testing product of a permanent matrix with a permanent matrix...",/)')
  C = A * B
  if (output) call print_matrix( C, "F6.1" )

  if (output) write(*,'(/,"Testing product of a normal matrix with a transposed matrix...",/)')
  C = A * (.t.A)
  if (output) call print_matrix( C, "F6.1" )

  if (output) write(*,'(/,"Testing product of a transposed matrix with a normal matrix...",/)')
  C = (.t.B) * B
  if (output) call print_matrix( C, "F6.1" )

  if (output) write(*,'(/,"Testing product of a transposed matrix with a normal matrix...",/)')
  C = (.t.B) * B
  if (output) call print_matrix( C, "F6.1" )

  if (output) write(*,'(/,"Testing product of a transposed matrix with a transposed matrix...",/)')
  C = (.t.B) * (.t.A)
  if (output) call print_matrix( C, "F6.1" )

  if (output) write(*,'(/,"Testing direct product of a transposed matrix with a transposed matrix...",/)')
  if (output) call print_matrix( (.t.B) * (.t.A), "F6.1" )

  exec = loop

end do

call Fortrix_Shutdown()

end program test
