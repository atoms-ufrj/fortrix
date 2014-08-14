module scilab_style

use mString

implicit none

type scilab_matrix
  character(sl) :: name
  real(rb), allocatable :: value(:,:)
  type(scilab_matrix), pointer :: next => null()
  contains
    procedure :: read => scilab_matrix_read
    procedure :: print => scilab_matrix_print
end type

type scilab_matrix_list
  type(scilab_matrix), pointer :: first => null()
  type(scilab_matrix), pointer :: last => null()
  contains
    procedure :: read => scilab_matrix_list_read
    procedure :: scilab_matrix_list_add_vector
    procedure :: scilab_matrix_list_add_matrix
    generic   :: add => scilab_matrix_list_add_vector, &
                        scilab_matrix_list_add_matrix
    procedure :: print => scilab_matrix_list_print
    procedure :: value  => scilab_matrix_list_value
end type

contains

!-------------------------------------------------------------------------------

  subroutine scilab_matrix_read( this, unit, name )
    class(scilab_matrix), intent(inout) :: this
    integer,              intent(in)    :: unit
    character(*),         intent(in)    :: name
    integer :: narg, row, ncols
    character(sl) :: arg(30)
    real(8), allocatable :: aux(:,:)
    call next_command( unit, narg, arg )
    row = 0
    ncols = narg
    allocate( this % value(row,ncols) )
    do while (arg(1)(1:1) /= "]")
      row = row + 1
      if (narg /= ncols) call error( "wrong number of columns in row", &
                                     int2str(row), "of matrix", name )
      allocate( aux(row,ncols) )
      aux(1:row-1,:) = this % value
      aux(row,:) = str2real( arg(1:narg) )
      deallocate( this % value )
      call move_alloc( aux, this % value )
      call next_command( unit, narg, arg )
    end do
  end subroutine scilab_matrix_read

!-------------------------------------------------------------------------------

  subroutine scilab_matrix_print( this, unit )
    class(scilab_matrix), intent(in) :: this
    integer,              intent(in) :: unit
    integer :: i, nrows, ncols
    nrows = size(this % value,1)
    ncols = size(this % value,2)
    if ( (nrows == 1).and.(ncols == 1) ) then
      write(unit,'(A," = ",A,";")') &
        trim(this % name), trim(real2str(this % value(1,1)))
    else
      write(unit,'(A," = [")') trim(this % name)
      do i = 1, nrows
        write(unit,'(A)') trim(join(real2str(this % value(i,:))))
      end do
      write(unit,'("];")')
    end if
  end subroutine scilab_matrix_print

!-------------------------------------------------------------------------------

  subroutine scilab_matrix_list_read( this, unit )
    class(scilab_matrix_list), intent(inout) :: this
    integer,                   intent(in)    :: unit
    integer :: narg
    character(sl) :: arg(3)
    type(scilab_matrix), pointer :: matrix
    call next_command( unit, narg, arg )
    do while (narg > 0)
      if ((narg > 1).and.(arg(2) == "=")) then
        allocate( matrix )
        matrix % name = trim(arg(1))
        if (arg(3) == "[") then
          call matrix % read( unit, matrix % name )
        else
          allocate( matrix % value(1,1) )
          matrix % value = str2real(arg(3)(1:len_trim(arg(3))-1));
        end if
      end if
      if (associated(this%last)) then
        this % last % next => matrix
      else
        this % first => matrix    
      end if
      this % last => matrix
      call next_command( unit, narg, arg )
    end do
  end subroutine scilab_matrix_list_read

!-------------------------------------------------------------------------------

  subroutine scilab_matrix_list_add_matrix( this, name, value )
    class(scilab_matrix_list), intent(inout) :: this
    character(*),              intent(in)    :: name
    real(rb),                  intent(in)    :: value(:,:)
    type(scilab_matrix), pointer :: matrix
    allocate( matrix )
    matrix % name = trim(name)
    allocate( matrix % value(size(value,1),size(value,2)) )
    matrix % value = value
    if (associated(this%last)) then
      this % last % next => matrix
    else
      this % first => matrix    
    end if
    this % last => matrix
  end subroutine scilab_matrix_list_add_matrix

!-------------------------------------------------------------------------------

  subroutine scilab_matrix_list_add_vector( this, name, value )
    class(scilab_matrix_list), intent(inout) :: this
    character(*),              intent(in)    :: name
    real(rb),                  intent(in)    :: value(:)
    real(rb), allocatable :: matrix(:,:)
    allocate( matrix(size(value),1) )
    matrix(:,1) = value
    call scilab_matrix_list_add_matrix( this, name, matrix )
  end subroutine scilab_matrix_list_add_vector

!-------------------------------------------------------------------------------

  subroutine scilab_matrix_list_print( this, unit )
    class(scilab_matrix_list), intent(inout) :: this
    integer,                   intent(in)    :: unit
    type(scilab_matrix), pointer :: current
    current => this % first
    do while (associated(current))
      call current % print( unit )
      write(unit,*)
      current => current % next
    end do
  end subroutine scilab_matrix_list_print

!-------------------------------------------------------------------------------

  function scilab_matrix_list_value( this, name ) result( value )
    class(scilab_matrix_list), intent(inout) :: this
    character(*),              intent(in)    :: name
    real(rb),                  pointer       :: value(:,:)
    type(scilab_matrix), pointer :: current
    logical :: found
    found = .false.
    current => this % first
    do while (associated(current).and.(.not.found))
      found = trim(current % name) == trim(name)
      if (found) then
        allocate( value(size(current % value,1),size(current % value,2)) )
        value = current % value
      end if
      current => current % next
    end do
  end function scilab_matrix_list_value

!-------------------------------------------------------------------------------

end module scilab_style
