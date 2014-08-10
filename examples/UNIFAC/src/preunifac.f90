program preunifac

implicit none

integer :: inp, io, i, j, ntotal, nc, nsg, iMain, nmg
real(8) :: iR, iQ, Aij, Aji, Bij, Bji, Cij, Cji
character(100) :: name
character(256) :: compound_file, group_file, parameter_file
type tcompound
  character(100) :: name
  integer :: nsg
  integer, allocatable :: index(:), nrep(:)
end type tcompound
type(tcompound), allocatable :: compound(:), aux(:)
integer, allocatable :: iaux(:), index(:), nu(:,:), Main(:), M(:,:)
real(8), allocatable :: R(:), Q(:), A(:,:,:)

! Read file names:
write(*,'("Enter names of the compound, group, and parameter files:")')
write(*,'(">")',advance="no"); read(*,'(A)') compound_file
write(*,'(">")',advance="no"); read(*,'(A)') group_file
write(*,'(">")',advance="no"); read(*,'(A)') parameter_file

! Read all compounds from the compound file:
nc = 0
allocate( compound(nc) )
open( newunit = inp, file = compound_file, status = "old" )
read(inp,'(A)',iostat=io) name
write(*,'(/, "Compounds available:")')
do while ((io == 0).and.(trim(name) /= ""))
  nc = nc + 1
  allocate( aux(nc) )
  aux(1:nc-1) = compound
  aux(nc) % name = trim(name)
  read(inp,*) aux(nc) % nsg
  allocate( aux(nc) % index(aux(nc) % nsg) )
  allocate( aux(nc) % nrep(aux(nc) % nsg) )
  do i = 1, aux(nc) % nsg
    read(inp,*) aux(nc) % nrep(i), aux(nc) % index(i)
  end do
  call move_alloc( aux,compound )
  write(*,'(I2,") ",A)') i, trim(name)
  read(inp,'(A)',iostat=io) name
end do
close( inp )

! Select the mixture components among the available compounds:
ntotal = nc
nc = 0
allocate( index(nc) )
write(*,'("Enter mixture components (empty to finish):")')
write(*,'("> ")',advance="no")
read(*,'(A)',iostat=io) name
do while ((io == 0).and.(trim(name) /= ""))
  i = 1
  do while ((i <= ntotal).and.(trim(name) /= trim(compound(i) % name)))
    i = i + 1
  end do
  if (i <= ntotal) then
    nc = nc + 1
    allocate( iaux(nc) )
    iaux = [index,i]
    call move_alloc( iaux, index )
  else
    write(*,'("error: unavailable compound ",A)') trim(name)
    stop
  end if
  write(*,'("> ")',advance="no")
  read(*,'(A)',iostat=io) name
end do
allocate( aux(nc) )
aux = compound(index)
call move_alloc( aux, compound )
deallocate( index )
write(*,'(/,"The mixture components are:")')
do i = 1, nc
  write(*,'(I2,") ",A)') i, trim(compound(i) % name)
end do

! Identify involved subgroups:
ntotal = 0
do i = 1, nc
  ntotal = max(ntotal,maxval(compound(i) % index))
end do
allocate( index(ntotal) )
index = 0
do i = 1, nc
  index(compound(i) % index) = 1
end do
nsg = 0
do i = 1, ntotal
  if (index(i) /= 0) then
    nsg = nsg + 1
    index(i) = nsg
  end if
end do

! Build the incidence matrix:
allocate( nu(nc,nsg) )
nu = 0
do i = 1, nc
  do j = 1, compound(i) % nsg
    nu(i,index(compound(i) % index(j))) = compound(i) % nrep(j)
  end do
end do

! Retrieve parameters of involved subgroups:
allocate( Main(nsg), R(nsg), Q(nsg) )
open( newunit = inp, file = group_file, status = "old" )
read(inp,'(/)')
read(inp,'(I5,A15,I5,A17,2F12.0)',iostat=io) i, name, iMain, name, iR, iQ
do while ((io == 0).and.(i > 0).and.(i <= ntotal))
  j = index(i)
  if (j /= 0) then
    Main(j) = iMain
    R(j) = iR
    Q(j) = iQ
  end if
  read(inp,'(I5,A15,I5,A17,2F12.0)',iostat=io) i, name, iMain, name, iR, iQ
end do
close(inp)

! Build the selection matrix:
ntotal = maxval(Main)
index = 0
do i = 1, nsg
  index(Main(i)) = 1
end do
nmg = 0
do i = 1, ntotal
  if (index(i) /= 0) then
    nmg = nmg + 1
    index(i) = nmg
  end if
end do
allocate( M(nsg,nmg) )
M = 0
forall(i=1:nsg) M(i,index(Main(i))) = 1

! Retrieve interaction parameters:
allocate( A(3,nmg,nmg) )
A = 0.0_8
open( newunit = inp, file = parameter_file, status = "old" )
read(inp,'(//)')
read(inp,'(I7,I8,6F13.0)',iostat=io) i, j, Aij, Aji, Bij, Bji, Cij, Cji
do while ((io == 0).and.(i > 0).and.(i <= ntotal))
  if (j <= ntotal) then
    if ((index(i) /= 0).and.(index(j) /= 0)) then
      A(:,index(i),index(j)) = [Aij, Bij, Cij]
      A(:,index(j),index(i)) = [Aji, Bji, Cji]
    end if
  end if
  read(inp,'(I7,I8,6F13.0)',iostat=io) i, j, Aij, Aji, Bij, Bji, Cij, Cji
end do
close(inp)

call print_all( 6 )

contains

  subroutine print_all( unit )
    integer, intent(in) :: unit

    character(3) :: cn, mat = "ABC"

    write(unit,'("nc  = ",I3,";")') nc
    write(unit,'("nsg = ",I3,";")') nsg
    write(unit,'("nmg = ",I3,";")') nmg

    ! Print the incidence matrix:
    write(unit,'(/,"nu = [")')
    write(cn,'(I3)') nsg
    do i = 1, nc
      write(unit,'('//cn//'I3)') nu(i,:)
    end do
    write(unit,'("];")')

    ! Print volume and area parameters:
    write(unit,'(/,"R = [")')
    do i = 1, nsg
      write(unit,'(F7.4)') R(i)
    end do
    write(unit,'("];")')

    write(unit,'(/,"Q = [")')
    do i = 1, nsg
      write(unit,'(F7.4)') Q(i)
    end do
    write(unit,'("];")')

    ! Print the selection matrix:
    write(unit,'(/,"M = [")')
    write(cn,'(I3)') nmg
    do i = 1, nsg
      write(unit,'('//cn//'I3)') M(i,:)
    end do
    write(unit,'("];")')

    ! Print interaction parameter matrices:
    mat = "ABC"
    do i = 1, 3
      write(unit,'(/,A1," = [")') mat(i:i)
      do j = 1, nmg
        write(unit,'('//cn//'F12.4)') A(i,j,:)
      end do
      write(unit,'("];")')
    end do

  end subroutine print_all

end program preunifac
