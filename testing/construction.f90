program construction

use Fortrix

implicit none

call fortrix_startup()

do while (.true.)

call print_matrix( new_scalar(1.0_rb), "F6.1" )

print*
call print_matrix( new_matrix(2.0_rb), "F6.1" )

print*
call print_matrix( new_matrix(real([1,2,3,4],rb)), "F6.1" )

print*
call print_matrix( new_matrix(real(reshape([1,2,3,4],[2,2]),rb)), "F6.1" )

print*
call print_matrix( .d.new_matrix(real([1,2,3,4],rb)), "F6.1" )

print*
call print_matrix( .t.(.t.new_matrix(real(reshape([1,2,3,4],[2,2]),rb))), "F6.1" )

end do

call fortrix_shutdown()

end program construction
