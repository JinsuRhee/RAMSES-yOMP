!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_log(message, rname)
  use amr_commons
  implicit none
  character(LEN=*)::message
  character(LEN=*)::rname
  

  write(*,*) '!!----- SUBSUB log'
  write(*,*) '   in ', trim(rname)
  write(*,*) '      ', trim(message)
  write(*,*) '   myid = ', myid
  write(*,*) '!!-----'

end subroutine subsub_log
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_findcell(x, y, z, ind_cell, ind_grid, ind_level, flag)
  !!-----
  !! This routine finds a cell which (x, y, z) is located at
  !! (x,y,z) should be in a code unit
  !!-----
  use amr_commons
  implicit none
  real(dp), intent(in) :: x, y, z
  integer, intent(out) :: ind_cell, ind_grid, ind_level
  logical :: flag

  ! Local variables
  integer :: ilevel, ind, ix, iy, iz, igrid, i
  integer :: ncache, ngrid, found
  integer :: index_cell, index_grid, iskip
  real(dp),dimension(1:twotondim,1:3):: xc
  real(dp),dimension(1:3)::skip_loc
  real(dp) :: cx, cy, cz, dr_cell
  real(dp) :: dx
  logical :: ng, okay

  ! compute cell-centered position first relative to grid center
  do ind=1,twotondim
    iz=(ind-1)/4
    iy=(ind-1-4*iz)/2
    ix=(ind-1-2*iy-4*iz)
    xc(ind,1)=(dble(ix)-0.5D0)
    xc(ind,2)=(dble(iy)-0.5D0)
    xc(ind,3)=(dble(iz)-0.5D0)
  end do
  
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)


  ! Find a coarse cell containing (x,y,z)
  dx=0.5D0**levelmin
  ncache=active(levelmin)%ngrid

  ng = .false.
  found = 0
  ind_cell = -1
  ind_grid = -1

  !$omp parallel do default(shared) &
  !$omp & private(igrid, ngrid, i, index_grid, index_cell, &
  !$omp & ind, iskip, cx, cy, cz, dr_cell) &
  !$omp & schedule(dynamic)
  do igrid=1,ncache,nvector

    if(found /= 0) cycle

    ngrid=MIN(nvector,ncache-igrid+1)
    do i=1,ngrid
      index_grid=active(levelmin)%igrid(igrid+i-1)

      do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        index_cell=iskip+index_grid
        cx=(xg(index_grid,1)+xc(ind,1)*dx-skip_loc(1))
        cy=(xg(index_grid,2)+xc(ind,2)*dx-skip_loc(2))
        cz=(xg(index_grid,3)+xc(ind,3)*dx-skip_loc(3))

        dr_cell=MAX(ABS(cx-x),ABS(cy-y),ABS(cz-z))
        if(dr_cell.le.dx/2.0) then
          !$omp critical
          if(ind_cell .ge. 0) ng=.true. ! debugger
          found = 1
          ind_cell=index_cell
          ind_grid=index_grid
          !$omp end critical
          exit
        endif
      enddo

      if(found .gt. 0) exit
    end do
  enddo
  !$omp end parallel do

  if(ng)then
    call subsub_log('two more than coarse cell found', 'subsub_findcell_1')
    stop
  endif

  ! No corase cell in this id
  if(found .eq. 0) then
    flag=.false.
    return
  else
    flag=.true.
  endif

  ! Find the leaf cell containing (x,y,z) from the selected coarse cell
  ind_level = levelmin
  ng = .false.
  do
    if(son(ind_cell)==0) then
      exit
    endif
    
    ind_level = ind_level + 1
    dx=0.5D0**ind_level

    index_grid = son(ind_cell)
    do ind=1, twotondim
      iskip = ncoarse + (ind-1)*ngridmax

      index_cell = iskip + index_grid
      cx=(xg(index_grid,1)+xc(ind,1)*dx-skip_loc(1))
      cy=(xg(index_grid,2)+xc(ind,2)*dx-skip_loc(2))
      cz=(xg(index_grid,3)+xc(ind,3)*dx-skip_loc(3))
      dr_cell=MAX(ABS(cx-x),ABS(cy-y),ABS(cz-z))
      if(dr_cell.le.dx/2.0) then
        ind_cell=index_cell
        ind_grid=index_grid
        exit
      endif
    enddo

    if(ind_level .gt. nlevelmax) then
      ng = .true.
      exit
    endif
  enddo

  if(ng)then
    call subsub_log('leaf cell is not found', 'subsub_findcell_2')
    stop
  endif
end subroutine subsub_findcell
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_update
  !!-----
  !! This routine creates/removes a new subsub obj based on the update sink list
  !!-----
  use pm_commons
  use amr_commons
  use subsub_commons
  use mpi_mod
  implicit none

  !! Local variables
  integer :: i, j, subsub_sinkinmyid, info, subsub_nsink_old
  integer :: ind, i0, id0, iend
  integer, dimension(1:nsink) :: sink_hash, sink_hash_next, sink_ismatch
  integer, dimension(:), allocatable :: subsub_idmatch

  subsub_nsink_old = subsub_nsink

  sink_hash = -1
  sink_hash_next = -1
  sink_ismatch = -1

  !! Make Hash table
  do i=1, nsink
    ind = mod(abs(idsink(i)), nsink) + 1
    if(sink_hash(ind) .lt. 0) then
      sink_hash(ind) = i
    else
      i0 = sink_hash(ind)
      do
        if(sink_hash_next(i0) .lt. 0) then
          sink_hash_next(i0) = i
          exit
        else
          i0 = sink_hash_next(i0)
        endif
      enddo
    endif
  enddo

  !! Check the status based on the new sink idlist
  if(subsub_end .gt. 0) then
    allocate(subsub_idmatch(1:subsub_end))

    subsub_idmatch = -1

    do i=1, subsub_end
      id0 = subsub_obj(i)%sink_id

      !! match by hash
      ind = mod(abs(id0), nsink) + 1
      i0 = sink_hash(ind)

      if(i0 .lt. 0) cycle

      do
        if(id0 .eq. idsink(i0)) then !! survive
          subsub_idmatch(i) = i0
          sink_ismatch(i0) = 1
          exit
        else
          if(sink_hash_next(i0) .ge. 1) then
            i0 = sink_hash_next(i0)
          else
            exit
          endif
        endif
      enddo
    enddo

    !! remove
    iend = 1
    do i=1, subsub_end
      if(subsub_idmatch(i) .lt. 0) cycle
      subsub_obj(iend) = subsub_obj(i)  !! do not need to make a copy (i >= iend)
      subsub_obj(iend)%sink_ind = subsub_idmatch(i)
      iend = iend + 1
    enddo

    iend = iend - 1

    do i=iend+1, subsub_end
      call subsub_deallocate(subsub_obj(i))
    enddo

    subsub_end = iend
  endif

  !! Input newly formed sinks
  do i=1, nsink
    if(sink_ismatch(i).gt.0) cycle
    call subsub_create(i)
  enddo




  !!----- Update subsub_nsink
#ifndef WITHOUTMPI
  subsub_sinkinmyid = subsub_end
  subsub_nsink = 0 !! initialize to 0
  call MPI_ALLREDUCE(subsub_sinkinmyid,subsub_nsink,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)

  if(subsub_nsink .ne. nsink) then
    if(myid.eq.1) then
      call subsub_log('some sinks miss in domains', 'subsub_update')
      write(*,*) ' nsink = ', nsink
      write(*,*) ' updated = ', subsub_nsink
    endif
    call clean_stop
  endif
#else
  subsub_nsink = nsink
#endif





do i=1, ncpu
  if(myid .eq. i .and. subsub_end .gt. 0) then
    do j=1, subsub_end
      if(idsink(subsub_obj(j)%sink_ind) .ne. subsub_obj(j)%sink_id) then
        write(*,*) 'wrong indexing', i, j
      endif
    enddo
  endif
  do j=1, 10000000
  enddo
  call MPI_BARRIER(MPI_COMM_WORLD,info)
enddo



end subroutine subsub_update
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_backup(filename)
  use subsub_commons
  use amr_commons
  use pm_commons
  use mpi_mod
  implicit none
  character(len=80) :: filename
  type(subsub_type), dimension(1:nsink) :: subsub_dump

  ! Local variables
  integer :: i, j, unit_out

  if (verbose) write(*,*) 'Entering backup_subsub'

  !!----- Gather all subsub_obj
#ifndef WITHOUTMPI
  call subsub_gather(subsub_dump)
#else
  do i=1, nsink
    subsub_dump(i) = subsub_obj(i)
  enddo
#endif

  !!----- Write
  if(myid .eq. 1) then
    open(newunit=unit_out, file=TRIM(filename), form='unformatted')

    rewind(unit_out)
    ! Write header
    write(unit_out) subsub_nsink
    write(unit_out) subsub_level
    write(unit_out) subsub_nhydro
    write(unit_out) ndim

    do i=1, subsub_nsink
      write(unit_out) subsub_dump(i)%sink_ind          !! sink index
      write(unit_out) subsub_dump(i)%sink_id           !! sink id
      write(unit_out) subsub_dump(i)%mass_tot          !! total mass
      write(unit_out) subsub_dump(i)%vxc, subsub_dump(i)%vyc, subsub_dump(i)%vzc !! momentum x, y, z
      write(unit_out) subsub_dump(i)%clevel            !! cell level

      do j=1, ndim
        write(unit_out) subsub_dump(i)%vg(:,j)           !! vx, vy, vz
      enddo

      do j=1, subsub_nhydro
        write(unit_out) subsub_dump(i)%hydro(:,j)
      enddo
    enddo
  endif 

  !!----- Deallocate
  if(myid .eq. 1) then
    do i=1, nsink
      call subsub_deallocate(subsub_dump(i))
    enddo
  endif
end subroutine subsub_backup
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_gather(subsub_dump)
  use subsub_commons
  use amr_commons
  use pm_commons
  use mpi_mod
  implicit none

  type(subsub_type), dimension(1:nsink) :: subsub_dump

  ! Local Variables
  integer :: i, j, info, iend
  integer, dimension(1:ncpu) :: subsub_dumpN
  integer :: subsub_nn1, subsub_nn2

  integer,parameter::tag1=6421,tag2=6422,tag3=6423
  integer,parameter::tag4=6424,tag5=6425,tag6=6426
  integer,parameter::tag7=6427,tag8=6428,tag9=6429


  !!----- Gather each size of obj
  if(myid .eq. 1) then
    subsub_dumpN = 0
    subsub_dumpN(1) = subsub_end

    do i=2, ncpu
      call MPI_RECV(subsub_dumpN(i), 1, MPI_INTEGER, i-1, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
    enddo
  else
    call MPI_SEND(subsub_end, 1, MPI_INTEGER, 0, tag1, MPI_COMM_WORLD, info)
  endif

  !!----- Debug
  !do i=1, ncpu
  !  if(myid .eq. i) write(*,*) i, subsub_end
  !  do j=1, 10000000
  !  enddo
  !  call MPI_BARRIER(MPI_COMM_WORLD,info)
  !enddo


  if(myid .eq. 1) then 
    if(sum(subsub_dumpN) .ne. nsink)then
      if(myid.eq.1) call subsub_log('nsink is not compatible with the sum of subsub_obj in each cpu', 'subsub_gather')

      write(*,*) ' nsink = ', nsink
      write(*,*) ' tot obj = ', sum(subsub_dumpN)
      do i=1, ncpu
        write(*,*) subsub_dumpN(i)
      enddo
      stop
    endif
  endif

  !!----- Allocate first
  if(myid .eq. 1) then
    do i=1, nsink
      call subsub_allocate(subsub_dump(i))
    enddo

    do i=1, subsub_end
      subsub_dump(i) = subsub_obj(i)
    enddo
  endif

  !!----- Receive
  subsub_nn1 = subsub_level**3 * ndim
  subsub_nn2 = subsub_level**3 * subsub_nhydro
  iend = subsub_end+1
  if(myid .eq. 1) then
    do i=2, ncpu
      if(subsub_dumpN(i).eq.0) cycle

      do j=1, subsub_dumpN(i)

        call MPI_RECV(subsub_dump(iend)%sink_ind, 1, MPI_INTEGER,  i-1, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_dump(iend)%sink_id,  1, MPI_INTEGER,  i-1, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_dump(iend)%mass_tot, 1, subsub_mpidp, i-1, tag3, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_dump(iend)%vxc, 1,      subsub_mpidp, i-1, tag4, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_dump(iend)%vyc, 1,      subsub_mpidp, i-1, tag5, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_dump(iend)%vzc, 1,      subsub_mpidp, i-1, tag6, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_dump(iend)%clevel, 1, MPI_INTEGER,    i-1, tag7, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)

        call MPI_RECV(subsub_dump(iend)%vg, subsub_nn1,    subsub_mpidp, i-1, tag8, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_dump(iend)%hydro, subsub_nn2, subsub_mpidp, i-1, tag9, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        iend = iend + 1
      enddo
    enddo
  else
    if(subsub_end .gt. 0) then
      do i=1, subsub_end
        call MPI_SEND(subsub_obj(i)%sink_ind, 1, MPI_INTEGER,  0, tag1, MPI_COMM_WORLD,info)
        call MPI_SEND(subsub_obj(i)%sink_id,  1, MPI_INTEGER,  0, tag2, MPI_COMM_WORLD,info)
        call MPI_SEND(subsub_obj(i)%mass_tot, 1, subsub_mpidp, 0, tag3, MPI_COMM_WORLD,info)
        call MPI_SEND(subsub_obj(i)%vxc, 1,      subsub_mpidp, 0, tag4, MPI_COMM_WORLD,info)
        call MPI_SEND(subsub_obj(i)%vyc, 1,      subsub_mpidp, 0, tag5, MPI_COMM_WORLD,info)
        call MPI_SEND(subsub_obj(i)%vzc, 1,      subsub_mpidp, 0, tag6, MPI_COMM_WORLD,info)
        call MPI_SEND(subsub_obj(i)%clevel, 1, MPI_INTEGER,    0, tag7, MPI_COMM_WORLD,info)

        call MPI_SEND(subsub_obj(i)%vg, subsub_nn1,    subsub_mpidp, 0, tag8, MPI_COMM_WORLD,info)
        call MPI_SEND(subsub_obj(i)%hydro, subsub_nn1, subsub_mpidp, 0, tag9, MPI_COMM_WORLD,info)
      enddo
    endif
  endif
end subroutine subsub_gather
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_allocate(subsub_dummy)
  use subsub_commons
  use amr_commons
  implicit none
  type(subsub_type) :: subsub_dummy

  allocate(subsub_dummy%vg(1:subsub_level**3, 1:ndim))
  allocate(subsub_dummy%hydro(1:subsub_level**3, 1:subsub_nhydro))  ! 1 for density
end subroutine subsub_allocate
subroutine subsub_deallocate(subsub_dummy)
  use subsub_commons
  implicit none
  type(subsub_type) :: subsub_dummy

  subsub_dummy%sink_ind = 0
  subsub_dummy%sink_id = 0
  subsub_dummy%mass_tot = 0d0
  subsub_dummy%vxc = 0d0
  subsub_dummy%vyc = 0d0
  subsub_dummy%vzc = 0d0
  subsub_dummy%clevel = 0
     
  deallocate(subsub_dummy%vg)
  deallocate(subsub_dummy%hydro)
end subroutine subsub_deallocate
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_printsink
use amr_commons
use subsub_commons
use mpi_mod
implicit none
integer i, j, info
do i=1, ncpu
  if(myid .eq. i .and. subsub_end .gt. 0) write(*,*) i, subsub_end, subsub_nsink
  do j=1, 10000000
  enddo
  call MPI_BARRIER(MPI_COMM_WORLD,info)
enddo
end subroutine subsub_printsink

!subroutine subsub_cellinsphere(x, y, z, r, ind_cell, ind_level)
! !!-----
!  !! This routine finds cells that touch the (x, y, z, r) sphere
!  !! (x,y,z) should be in a code unit
!  use amr_commons
!  implicit none
!  real(dp), intent(in) :: x, y, z, r
!  integer, dimension(:), allocatable, intent(out) :: ind_cell, ind_level
!  
!
!  ! Local variables
!  integer :: ilevel, ind, ix, iy, iz, igrid, i
!  integer :: ncache, ngrid, found, lind
!  integer :: index_cell, index_grid, iskip
!  real(dp),dimension(1:twotondim,1:3):: xc
!  real(dp),dimension(1:3)::skip_loc
!  real(dp) :: cx, cy, cz, dr_cell
!  real(dp) :: dx
!  logical :: ng, okay
!
!  ! compute cell-centered position first relative to grid center
!  do ind=1,twotondim
!    iz=(ind-1)/4
!    iy=(ind-1-4*iz)/2
!    ix=(ind-1-2*iy-4*iz)
!    xc(ind,1)=(dble(ix)-0.5D0)
!    xc(ind,2)=(dble(iy)-0.5D0)
!    xc(ind,3)=(dble(iz)-0.5D0)
!  end do
!  
!  skip_loc=(/0.0d0,0.0d0,0.0d0/)
!  skip_loc(1)=dble(icoarse_min)
!  skip_loc(2)=dble(jcoarse_min)
!  skip_loc(3)=dble(kcoarse_min)
!
!
!  ! Find a coarse cell containing (x,y,z)
!  dx=0.5D0**levelmin
!  ncache=active(levelmin)%ngrid
!
!  ng = .false.
!  found = 0
!  ind_cell = -1
!  ind_level = -1
!  lind = 0
!
!  allocate(ind_cell(1:ngridmax*twotondim))
!  allocate(ind_level(1:ngridmax*twotondim))
!
!  !$omp parallel do default(shared) &
!  !$omp & private(igrid, ngrid, i, index_grid, index_cell, &
!  !$omp & ind, iskip, cx, cy, cz, dr_cell) &
!  !$omp & schedule(dynamic)
!  do igrid=1,ncache,nvector
!   ngrid=MIN(nvector,ncache-igrid+1)
!    do i=1,ngrid
!      index_grid=active(ilevel)%igrid(igrid+i-1)
!
!
!
!      
!      do ind=1,twotondim
!        iskip=ncoarse+(ind-1)*ngridmax
!        index_cell=iskip+index_grid
!        cx=(xg(index_grid,1)+xc(ind,1)*dx-skip_loc(1))
!        cy=(xg(index_grid,2)+xc(ind,2)*dx-skip_loc(2))
!        cz=(xg(index_grid,3)+xc(ind,3)*dx-skip_loc(3))
!
!        dr_cell=MAX(ABS(cx-x),ABS(cy-y),ABS(cz-z))
!        if(dr_cell.le.dx/2.0) then
!          !$omp critical
!          if(ind_cell .ge. 0) ng=.true. ! debugger
!          found = 1
!          ind_cell=index_cell
!          ind_grid=index_grid
!          !$omp end critical
!          exit
!        endif
!      enddo
!
!      if(found .gt. 0) exit
!    end do
