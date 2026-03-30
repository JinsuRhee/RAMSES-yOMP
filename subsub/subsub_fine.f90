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
  real(dp) :: ttsta, ttend, tcheck(10)

  if(myid .eq. 1) then
    ttsta=MPI_WTIME()
  endif

  if(nsink .eq. 0) return

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


    deallocate(subsub_idmatch)
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
    if(mod(abs(idsink(i)),ncpu)+1 .ne. myid) cycle
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

  !!----- Update Domain
  !! This should be ahead of subsub_compute
  call subsub_updatedomain

  !!----- Update Property
  call subsub_compute


  !!----- Save coarse timestep output
  if(subsub_savecoarse) then
#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD, info)
#endif
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

  call MPI_ALLREDUCE(subsub_tcheck_cg(1), subsub_tcheck_cg_global(1), 1, subsub_mpidp, MPI_MAX, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(subsub_tcheck_cg(2), subsub_tcheck_cg_global(2), 1, subsub_mpidp, MPI_MAX, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(subsub_tcheck_cg(3), subsub_tcheck_cg_global(3), 1, subsub_mpidp, MPI_MAX, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(subsub_tcheck_cg(4), subsub_tcheck_cg_global(4), 1, subsub_mpidp, MPI_MAX, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(subsub_tcheck_cg(5), subsub_tcheck_cg_global(5), 1, subsub_mpidp, MPI_MAX, MPI_COMM_WORLD, info)

  call MPI_ALLREDUCE(subsub_ncheck_cg(1), subsub_ncheck_cg_global(1), 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(subsub_ncheck_cg(2), subsub_ncheck_cg_global(2), 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(subsub_ncheck_cg(3), subsub_ncheck_cg_global(3), 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, info)

if(myid.eq.1) then
  write(*,*)'     poisson      = ', subsub_tcheck_cg_global(1)
  write(*,*)'     force        = ', subsub_tcheck_cg_global(2)
  write(*,*)'     dt           = ', subsub_tcheck_cg_global(3)
  write(*,*)'     kick         = ', subsub_tcheck_cg_global(4)
  write(*,*)'     drift        = ', subsub_tcheck_cg_global(5)
  write(*,*)'     poisson_nall = ', subsub_ncheck_cg_global(1)
  write(*,*)'     poisson_n1   = ', subsub_ncheck_cg_global(2)
  write(*,*)'     poisson_n2   = ', subsub_ncheck_cg_global(3)
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
  !real(dp):: dx, tff, tffnew, threepi2, fourpi, vffnew
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v


  if(subsub_end .eq. 0) return

  !! constants
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)

  !! update sink by sink
  do i=1, subsub_end

    sinkind = subsub_obj(i)%sink_ind

    !call subsub_findcell(xsink(sinkind,1)/scale, xsink(sinkind,2)/scale, xsink(sinkind,3)/scale, &
    !   ind_cell, ind_grid, ind_level, subflag)

    !! debugger
    !if(.not. subflag) then
    !  call subsub_log('sink is not found in this domain', 'subsub_compute')
    !  write(*,*) 'myid = ', myid
    !  write(*,*) 'sinkid = ', idsink(sinkind)
    !  stop
    !endif


    !! update mass_tot
    mtot_old = subsub_obj(i)%mass_tot
    mbh_old = subsub_obj(i)%sink_mass
    mtot_new = (subsub_boxlen**ndim) * max(subsub_obj(i)%uold(1,1), subsub_smallr)
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

!! neglect mass conversion (starts from initial mass)
   !subsub_obj(i)%mass_tot = mtot_new
   subsub_obj(i)%sink_mass = msink(sinkind)
  enddo
  

end subroutine subsub_compute