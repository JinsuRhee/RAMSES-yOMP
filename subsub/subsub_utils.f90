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
subroutine subsub_findcell(x, y, z, ind_cell, ind_grid, ind_level, subsub_flag)
  !!-----
  !! This routine finds a cell which (x, y, z) is located at
  !! (x,y,z) should be in a code unit
  !!-----
  use amr_commons
  implicit none
  real(dp), intent(in) :: x, y, z
  integer, intent(out) :: ind_cell, ind_grid, ind_level
  logical :: subsub_flag

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
    subsub_flag=.false.
    return
  else
    subsub_flag=.true.
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
subroutine subsub_finddomain(x, y, z, mpinum)
  !!-----
  !! This routine finds a domain containing (x, y, z)
  !! Should be changed to more claver one
  !!-----
  use amr_commons
  implicit none
  real(dp), intent(in) :: x, y, z
  integer :: mpinum

  ! Local variables
  integer :: ilevel, ind, ix, iy, iz, igrid, i
  integer :: ncache, ngrid, found
  integer :: index_cell, index_grid, iskip
  integer :: ind_grid, ind_cell
  real(dp),dimension(1:twotondim,1:3):: xc
  real(dp),dimension(1:3)::skip_loc
  real(dp) :: cx, cy, cz, dr_cell
  real(dp) :: dx
  logical :: ng, okay

  mpinum = -1
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
    call subsub_log('two more than coarse cell found', 'subsub_finddomain')
    stop
  endif

  ! No corase cell in this id
  if(found .eq. 0) then
    !flag=.false.
    mpinum = -1
  else
    !flag=.true.
    mpinum = myid
  endif

  return
end subroutine subsub_finddomain
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
  use subsub_parameters
  use subsub_commons
  use mpi_mod
  implicit none

  !! Local variables
  integer :: i, j, subsub_sinkinmyid, info, subsub_nsink_old
  integer :: ind, i0, id0, iend, ismy, isokay
  integer, dimension(1:nsink) :: sink_hash, sink_hash_next, sink_ismatch, sink_ismatch_all
  integer, dimension(1:nsink) :: ismysink, ismysink_dump
  integer, dimension(1:nsink) :: nowmysink, nowmysink_dump
  integer, dimension(:), allocatable :: subsub_idmatch
  integer :: nx_loc
  real(dp):: scale
  character(LEN=5) :: nchar
  character(LEN=80) :: filename, dirname
  real(dp) :: ttsta, ttend  

  if(myid .eq. 1) then
    ttsta=MPI_WTIME()
  endif

  subsub_nsink_old = subsub_nsink

  sink_hash = -1
  sink_hash_next = -1
  sink_ismatch = 0

  !!----- Make Hash table
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

    !! ----- RHEE -----
    !! Maybe, do not need to implement OMP here at the moment?
    !! ----------------
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


  !!----- Make newly formed sinks
#ifndef WITHOUTMPI
  sink_ismatch_all = 0
  call MPI_ALLREDUCE(sink_ismatch,sink_ismatch_all,nsink,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
  sink_ismatch_all = sink_ismatch
#endif

  do i=1, nsink
    if(sink_ismatch_all(i).gt.0) cycle
    call subsub_create(i)
  enddo

  !!----- Update Domain
#ifndef WITHOUTMPI
  nowmysink(:) = 0
  nowmysink_dump(:) = 0

  if(subsub_end .gt. 0) then
    do i=1, subsub_end
      nowmysink_dump(subsub_obj(i)%sink_ind) = myid
    enddo
  endif
  call MPI_ALLREDUCE(nowmysink_dump,nowmysink,nsink,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)

  ismysink(:) = 0
  ismysink_dump(:) = 0

  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  
  do i=1, nsink
    ismy = -1
    call subsub_finddomain(xsink(i,1)/scale, xsink(i,2)/scale, xsink(i,3)/scale, ismy)
    if(ismy .gt. 0) ismysink_dump(i) = myid
  enddo
  call MPI_ALLREDUCE(ismysink_dump,ismysink,nsink,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)


  do i=1, nsink
    if(nowmysink(i) .ne. ismysink(i)) then
      call subsub_sendtoanother(i, nowmysink(i), ismysink(i))
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, info)
  enddo
#endif


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

  !!----- Update Property
  call subsub_compute


  !!----- Save coarse timestep output
  if(subsub_savecoarse) then
    dirname = 'SUBSUBPROPS'
    call create_output_dirs_nobar(dirname)
    call title(nstep_coarse,nchar)
    filename = TRIM(dirname)//'/subsub_'//TRIM(nchar)//'.dat'
    call subsub_backup(filename)
  endif


!! simple debugging
!do i=1, ncpu
!  if(myid .eq. i .and. subsub_end .gt. 0) then
!    do j=1, subsub_end
!      if(idsink(subsub_obj(j)%sink_ind) .ne. subsub_obj(j)%sink_id) then
!        write(*,*) 'wrong indexing', i, j
!      endif
!    enddo
!  endif
!  do j=1, 10000000
!  enddo
!  call MPI_BARRIER(MPI_COMM_WORLD,info)
!enddo
#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD, info)
#endif

  if(myid .eq. 1) then
    ttend=MPI_WTIME()
    write(*,*) ' Time elapsed in SUBSUB_compute [sec]', sngl(ttend-ttsta)
  endif

  

end subroutine subsub_update
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_compute
  use amr_commons
  use pm_commons
  use hydro_commons
  use subsub_commons
  use subsub_parameters
  use mpi_mod
  implicit none

  !! Local variables
  integer :: i, j
  real(dp):: scale, mtot_old, mtot_new, mbh_old, mbh_new
  integer :: ind_cell, ind_grid, ind_level, nx_loc, sinkind
  logical :: subflag
  real(dp):: dx, tff, tffnew, threepi2, fourpi, vffnew
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v


  if(subsub_end .eq. 0) return

  !! constants
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)

  !! update sink by sink
  do i=1, subsub_end

    sinkind = subsub_obj(i)%sink_ind

    call subsub_findcell(xsink(sinkind,1)/scale, xsink(sinkind,2)/scale, xsink(sinkind,3)/scale, &
       ind_cell, ind_grid, ind_level, subflag)

    !! debugger
    if(.not. subflag) then
      call subsub_log('sink is not found in this domain', 'subsub_compute')
      write(*,*) 'myid = ', myid
      write(*,*) 'sinkid = ', idsink(sinkind)
      stop
    endif
    !! update mass_tot
    dx=0.5D0**ind_level

    mtot_old = subsub_obj(i)%mass_tot
    mbh_old = subsub_obj(i)%sink_mass
    mtot_new = (subsub_boxlen**ndim) * max(uold(ind_cell,1), subsub_smallr)
    mbh_new = msink(sinkind)

   !! update density by the mass change
   !!----- RHEE ------
   !! - Delta M = mass_inflow - mass_outflow
   !! - mass_in:
   !! 		spherical inflow / rotating inflow / anisotropic inflow
   !! - mass_out:
   !! 		by jet or outflow?
   !!		by tidal mass loss?
   !!
   !! - At the moment,  if Delta M > 0: all by spherical inflow (v_in || vec(r) )
   !! 					if Delta M < 0: all by spherical outflow (v_out || vec(r) )
   !! BC as the mass flux
   !! Mass flux as the uniform mass in/outflow
   !! Infall velocity of v_ff
!   if(idsink(sinkind).eq.829) then
!     subsub_obj(i)%hydro(:,1) = mtot_old / subsub_boxlen**3
!     !write(*,*) '%456456 rho_min = ', minval(subsub_obj(i)%hydro(:,1))
!     !write(*,*) '%456456 rho_max = ', maxval(subsub_obj(i)%hydro(:,1))
!     !write(*,*) '%456456 new rho = ', uold(ind_cell,1)
!     !write(*,*) '%456456 BH mass = ', msink(subsub_obj(i)%sink_ind) 
!     !write(*,*) '%456456 BH mass = ', msink(subsub_obj(i)%sink_ind) * scale_d*scale_l**3 / 2d33 / 1e6
!
!     call subsub_iget3ind(j, 2, 2, 2)
!     write(*,*) '%456456 den(2,2,2)', subsub_obj(i)%hydro(j,1)
!     write(*,*) '%456456 den(1,1,1)', subsub_obj(i)%hydro(1,1)
!   endif
   call subsub_computefine(i, mtot_new-mtot_old, mbh_new-mbh_old)

!if(mtot_new-mtot_old .gt. 0) write(*,*) 'good sink = ', idsink(sinkind)
!   if(idsink(sinkind) .eq. 829) then
!     write(*,*) '%456456 mass_old = ', mtot_old
!     write(*,*) '%456456 mass_new = ', mtot_new
!     write(*,*) '%456456 dm =', mtot_new - mtot_old
!     write(*,*) '%456456 dmbh =', mbh_new - mbh_old
!     write(*,*) '%456546 dmtot =', mtot_new - mtot_old + mbh_new - mbh_old
!     write(*,*) '%456456 by ct = ', sum(subsub_obj(i)%hydro(:,1)) * (subsub_boxlen/subsub_ngrid)**ndim
!     call subsub_iget3ind(j, 2, 2, 2)
!     write(*,*) '%456456 den(2,2,2)', subsub_obj(i)%hydro(j,1)
!     write(*,*) '%456456 den(1,1,1)', subsub_obj(i)%hydro(1,1)
!     !stop
!   endif
   !rho_flux = ((mtot_new - mtot_old) / dx**2) / 6.0D0 !! Mass/area
   !rho_flux = rho_flux / (dx / real(subsub_level,kind=dp))**3.0D0 !! from ghost cells (rho/area)

   subsub_obj(i)%mass_tot = mtot_new
   subsub_obj(i)%sink_mass = msink(sinkind)
  enddo
  

end subroutine subsub_compute
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
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_input(subsub_dummy)
  use amr_commons
  use subsub_commons
  use mpi_mod
  implicit none
  type(subsub_type) :: subsub_dummy


  subsub_end = subsub_end + 1
  if(subsub_end .gt. size(subsub_obj)) then
   if(myid .eq. 1) call subsub_log('subsub_object exceeds nsinkmax', 'subsub_input')
   call clean_stop
  endif

  subsub_obj(subsub_end) = subsub_dummy


end subroutine subsub_input
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_get3ind(i, ix, iy, iz)
  use subsub_commons
  implicit none

  integer :: i, ix, iy, iz

  iz = mod((i-1)/subsub_ngrid**2, subsub_ngrid)
  iy = mod(((i-1)-(subsub_ngrid**2)*iz)/subsub_ngrid, subsub_ngrid)
  ix = (i-1) - (subsub_ngrid**2) * iz - subsub_ngrid*iy

  ix = ix + 1
  iy = iy + 1
  iz = iz + 1

end subroutine subsub_get3ind
subroutine subsub_iget3ind(i, ix, iy, iz)
  use subsub_commons
  implicit none

  integer :: i, ix, iy, iz

  i = (iz-1)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix
end subroutine subsub_iget3ind
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_sortbyind(subsub_dummy, nn)
  use subsub_commons
  implicit none

  integer :: nn
  type(subsub_type), dimension(1:nn) :: subsub_dummy

  !! Local variables
  integer :: i, j
  type(subsub_type) :: temp

  do i=2, nn
    temp = subsub_dummy(i)
    j = i-1

    do while (j >= 1 .and. subsub_dummy(j)%sink_ind > temp%sink_ind)
      subsub_dummy(j+1) = subsub_dummy(j)
      j = j - 1
    enddo

    subsub_dummy(j+1) = temp
  enddo
end subroutine subsub_sortbyind
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_getdist2(ix, iy, iz, dd2)
  use subsub_commons
  implicit none

  integer :: ix, iy, iz
  real(dp) :: dd2, rx, ry, rz
  real(dp) :: dx
  real(dp), dimension(3) :: center

  dx = subsub_boxlen / dble(subsub_ngrid)
  center = (/subsub_boxlen/2.0D0, subsub_boxlen/2.0D0, subsub_boxlen/2.0D0/)

  rx = (dble(ix) - 0.5D0) * dx
  ry = (dble(iy) - 0.5D0) * dx
  rz = (dble(iz) - 0.5D0) * dx
  
  dd2 = (rx-center(1))**2 + (ry-center(2))**2 + (rz-center(2))**2
  
end subroutine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_finddt(subsub_dt, varr)
  use subsub_commons
  use amr_commons
  implicit none

  real(dp) :: subsub_dt
  real(dp), dimension(1:10) :: varr

  !! Local variables
  real(dp) :: subsub_maxv
  real(dp) :: subsub_dx
  real(dp) :: subsub_dtmax

  subsub_dx = subsub_boxlen / dble(subsub_ngrid)  
  subsub_maxv = maxval(varr)

  !! including BC
  subsub_dtmax = subsub_cfl * (subsub_dx / subsub_maxv)

  do
    if(subsub_dt .lt. subsub_dtmax) exit
    subsub_dt = subsub_dt * 0.5D0
  enddo
end subroutine subsub_finddt
!################################################################
!################################################################
!################################################################
!################################################################  
!! TO DO SOMEDAYS

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
