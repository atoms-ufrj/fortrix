program unifac

use fortrix
use scilab_style

implicit none

!character(100) :: input_file

type(scilab_matrix_list) :: list

call list % read( 5 )
call list % print( 6 )

end program unifac
