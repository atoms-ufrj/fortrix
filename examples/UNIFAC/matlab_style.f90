module matlab_style

use mGlobal
use mString

implicit none

type matlab_matrix
  character(sl) :: name
  real(rb), allocatable :: value(:,:)
  type(matlab_matrix), pointer :: next => null()
  contains
    procedure :: read => matlab_matrix_read
end type

type matlab_matrix_list
  type(matlab_matrix), pointer :: first => null()
  contains
    procedure :: read => matlab_matrix_list_read
    procedure :: print => matlab_matrix_list_print
    procedure :: value  => matlab_matrix_list_value
end type

contains

!-------------------------------------------------------------------------------

  subroutine matlab_matrix_read( this, unit, name )
    class(matlab_matrix), intent(inout) :: this
    integer,              intent(in)    :: unit
    character(*),         intent(in)    :: name
    integer :: narg, row, ncols, i
    character(sl) :: arg(30)
    real(8), allocatable :: aux(:,:)
    this % name = name
    call next_command( unit, narg, arg )
    row = 0
    ncols = narg
    allocate( this % value(row,ncols) )
    do while (trim(arg(1)(1:1)) /= "]")
      row = row + 1
      if (narg /= ncols) call error( "wrong number of columns in row", &
                                     int2str(row), "of matrix", name )
      allocate( aux(row,ncols) )
      aux(1:row-1,:) = this % value
      do i = 1, ncols
        aux(row,i) = str2real( arg(i) )
      end do
      deallocate( this % value )
      call move_alloc( aux, this % value )
      call next_command( unit, narg, arg )
    end do
  end subroutine matlab_matrix_read

!-------------------------------------------------------------------------------

  subroutine matlab_matrix_list_read( this, unit )
    class(matlab_matrix_list), intent(inout) :: this
    integer,                   intent(in)    :: unit
    integer :: narg
    character(sl) :: arg(3)
    type(matlab_matrix), pointer :: matrix
    call next_command( unit, narg, arg )
    do while (narg > 0)
      if ( (narg > 1).and.(arg(2)(1:1) == "=").and.(arg(3)(1:1) == "[") ) then
        allocate( matrix )
        call matrix % read( unit, trim(arg(1)) )
        matrix % next => this % first
        this % first => matrix
      end if
      call next_command( unit, narg, arg )
    end do
  end subroutine matlab_matrix_list_read

!-------------------------------------------------------------------------------

  subroutine matlab_matrix_list_print( this )
    class(matlab_matrix_list), intent(inout) :: this
    type(matlab_matrix), pointer :: current
    integer :: i
    character(3) :: C
    current => this % first
    do while (associated(current))
      write(*,'(A," = [")') trim(current % name)
      write(C,'(I3)') size(current % value, 2)
      do i = 1, size(current % value,1)
        write(*,'('//C//'F6.1)') current % value(i,:)
      end do
      write(*,'("];")')
      current => current % next
    end do
  end subroutine matlab_matrix_list_print

!-------------------------------------------------------------------------------

  function matlab_matrix_list_value( this, name ) result( value )
    class(matlab_matrix_list), intent(inout) :: this
    character(*),              intent(in)    :: name
    real(rb),                  pointer       :: value(:,:)
    type(matlab_matrix), pointer :: current
    logical :: found
    current => this % first
    found = .false.
    do while (associated(current).and.(.not.found))
      found = trim(current % name) == name
      if (found) then
        allocate( value(size(current % value,1),size(current % value,2)) )
        value = current % value
      end if
      current => current % next
    end do
  end function matlab_matrix_list_value

!-------------------------------------------------------------------------------

end module matlab_style


program matlab_style_test
use matlab_style
type(matlab_matrix_list) :: list
call list % read( 5 )
call list % print
print*
print*, list % value( "A" )
end program matlab_style_test
