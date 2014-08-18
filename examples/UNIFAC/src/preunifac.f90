program preunifac

use scilab_style

implicit none

integer :: inp, dat, out, io, i, j, ntotal, nsubs, nc, nsg, iMain, nmg
real(8) :: iR, iQ, Aij, Aji, Bij, Bji, Cij, Cji
character(100) :: name, model
character(256) :: compound_file, group_file, parameter_file, input_file, mixture
type tcompound
  character(100) :: name
  integer :: nsg
  integer, allocatable :: index(:), nrep(:)
end type tcompound
type(tcompound), allocatable :: compound(:), aux(:), component(:)
integer, allocatable :: iaux(:), index(:), nu(:,:), Main(:), M(:,:)
real(8), allocatable :: R(:), Q(:), A(:,:,:)
type(scilab_matrix_list) :: list

if (iargc() < 2) then
  write(*,'("Usage: ./preunifac model input-file")')
  stop
end if

call getarg( 1, model )
call getarg( 2, input_file )

compound_file = "parameters/" // trim(model) // ".compounds"
group_file = "parameters/" // trim(model) // ".groups"
parameter_file = "parameters/" // trim(model) // ".parameters"

! Read all compounds from the compound file:
nsubs = 0
allocate( compound(nsubs) )
open( newunit = dat, file = compound_file, status = "old" )
read(dat,'(A)',iostat=io) name
do while ((io == 0).and.(trim(name) /= ""))
  nsubs = nsubs + 1
  allocate( aux(nsubs) )
  aux(1:nsubs-1) = compound
  aux(nsubs) % name = trim(name)
  read(dat,*) aux(nsubs) % nsg
  allocate( aux(nsubs) % index(aux(nsubs) % nsg) )
  allocate( aux(nsubs) % nrep(aux(nsubs) % nsg) )
  do i = 1, aux(nsubs) % nsg
    read(dat,*) aux(nsubs) % nrep(i), aux(nsubs) % index(i)
  end do
  call move_alloc( aux, compound )
  read(dat,'(A)',iostat=io) name
end do
close( dat )

open( newunit = inp, file = input_file, status = "old" )

read(inp,'(A)',iostat=io) name
do while ((io == 0).and.(trim(name) /= ""))

  ! Select the mixture components among the available compounds:
  nc = 0
  allocate( index(nc) )
  do while ((io == 0).and.(trim(name) /= ""))
    i = 1
    do while ((i <= nsubs).and.(trim(name) /= trim(compound(i) % name)))
      i = i + 1
    end do
    if (i <= nsubs) then
      nc = nc + 1
      allocate( iaux(nc) )
      iaux = [index,i]
      call move_alloc( iaux, index )
    else
      write(*,'("error: unavailable compound ",A)') trim(name)
      stop
    end if
    read(inp,'(A)',iostat=io) name
  end do
  allocate( component(nc) )
  component = compound(index)
  deallocate( index )

  ! Identify involved subgroups:
  ntotal = 0
  do i = 1, nc
    ntotal = max(ntotal,maxval(component(i) % index))
  end do
  allocate( index(ntotal) )
  index = 0
  do i = 1, nc
    index(component(i) % index) = 1
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
    do j = 1, component(i) % nsg
      nu(i,index(component(i) % index(j))) = component(i) % nrep(j)
    end do
  end do
  call list % add( "nu", real(nu,rb) )

  ! Retrieve parameters of involved subgroups:
  allocate( Main(nsg), R(nsg), Q(nsg) )
  open( newunit = dat, file = group_file, status = "old" )
  read(dat,'(/)')
  read(dat,'(I5,A15,I5,A17,2F12.0)',iostat=io) i, name, iMain, name, iR, iQ
  do while ((io == 0).and.(i > 0).and.(i <= ntotal))
    j = index(i)
    if (j /= 0) then
      Main(j) = iMain
      R(j) = iR
      Q(j) = iQ
    end if
    read(dat,'(I5,A15,I5,A17,2F12.0)',iostat=io) i, name, iMain, name, iR, iQ
  end do
  close( dat )
  call list % add( "R", R )
  call list % add( "Q", Q )

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
  call list % add( "M", real(M,rb) )

  ! Retrieve interaction parameters:
  allocate( A(3,nmg,nmg) )
  A = 0.0_8
  open( newunit = dat, file = parameter_file, status = "old" )
  read(dat,'(//)')
  read(dat,'(I7,I8,6F13.0)',iostat=io) i, j, Aij, Aji, Bij, Bji, Cij, Cji
  do while ((io == 0).and.(i > 0).and.(i <= ntotal))
    if (j <= ntotal) then
      if ((index(i) /= 0).and.(index(j) /= 0)) then
        A(:,index(i),index(j)) = [Aij, Bij, Cij]
        A(:,index(j),index(i)) = [Aji, Bji, Cji]
      end if
    end if
    read(dat,'(I7,I8,6F13.0)',iostat=io) i, j, Aij, Aji, Bij, Bji, Cij, Cji
  end do
  close( dat )
  call list % add( "A", A(1,:,:) )
  call list % add( "B", A(2,:,:) )
  call list % add( "C", A(3,:,:) )

  mixture = trim(component(1) % name)
  do i = 2, nc
    mixture = trim(mixture)//"_"//trim(component(i) % name)
  end do

  open( newunit = out, file = trim(mixture)//"."//trim(model), status = "replace" )
  do i = 1, nc
    write(out,'("// ",A)') trim(component(i) % name)
  end do
  write(out,'(/,"// model: ",A,/)') trim(model)
  call list % print( out )
  close(out)
  write(*,'("Output file: ",A)') trim(mixture)//"."//trim(model)

  deallocate( index, component, nu, Main, R, Q, M, A )

  read(inp,'(A)',iostat=io) name

end do

close( inp )
end program preunifac
