!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_subsub
  use amr_commons
  use pm_commons
  use subsub_commons
  use mpi_mod
  implicit none
  !!-----
  !! This routine works to initialize the subsub variable to sink particles in myid
  !!-----

  !! Local Varaibles
  integer :: i, j
  integer,dimension(1:nsink)::itemp
  integer :: subsub_sinkinmyid, info
  integer :: ilun
  character(LEN=5)::nchar


  call subsub_precision_mpi()
  
  !!-----
  !! Check
  !!-----
  if(subsub_on .and. .not. sink) then
    if(myid .eq. 1) then
      call subsub_log('sink should be .true.', 'init_subsub')
    endif
    call clean_stop
  endif

  if(.not. hydro)return

  if(nsink.eq.0)return

  !!-----
  !! Retrieve from the restart
  !!-----
 

  !!-----
  !! Create and input obj
  !!-----

  !! TODO Get by IO
  if(verbose)write(*,*)'Entering init_subsub'

  allocate(subsub_obj(1:nsinkmax))
  subsub_end = 0

  do i=1, nsink
    call subsub_create(i)
  enddo

#ifndef WITHOUTMPI
  subsub_sinkinmyid = subsub_end
  subsub_nsink = 0 !! initialize to 0
  call MPI_ALLREDUCE(subsub_sinkinmyid,subsub_nsink,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)

  if(subsub_nsink .ne. nsink) then
    if(myid.eq.1) call subsub_log('some sinks miss in domains', 'init_subsub_2')
    call clean_stop
  endif
#else
  subsub_nsink = nsink
#endif

end subroutine init_subsub
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_create(sinkid)
  use subsub_commons
  use amr_commons
  use hydro_commons
  use pm_commons
  use mpi_mod
  implicit none
  !!-----
  !! This routine create subsub_obj and input the information of cell that the sink inherited
  !!-----
  integer, intent(in) ::sinkid

  

  !! Local Varaibles
  type(subsub_type) :: subsub_dummy
  integer :: i, ilevel, ind, ix, iy, iz, ncache, igrid, ngrid, nx_loc
  real(dp):: scale, scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp):: dx, dx_loc, vol_loc

  !integer,dimension(1:nvector)::ind_grid,ind_cell
  integer ::ind_grid, ind_cell, ind_level

  real(dp):: x,y,z, dxx, dyy, dzz, drr
  integer :: dlev
  logical :: okay, flag
  real(dp),dimension(1:3)::xbound,skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc

  !!-----
  !! Initialize
  !!-----
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)

  !!-----
  !! Retrieve related properties from the inherited cell
  !!-----
  flag = .false.
  call subsub_findcell(xsink(sinkid,1)/scale, xsink(sinkid,2)/scale, xsink(sinkid,3)/scale, &
    ind_cell, ind_grid, ind_level, flag)

  if(.not. flag) return !! no matched cell in this domain


  !! if not restart
  dx=0.5D0**ind_level
  subsub_dummy%sink_ind = sinkid
  subsub_dummy%sink_id  = idsink(sinkid)
  subsub_dummy%mass_tot=(dx**ndim) * max(uold(ind_cell,1), smallr)
  subsub_dummy%vxc=uold(ind_cell,2)
  subsub_dummy%vyc=uold(ind_cell,3)
  subsub_dummy%vzc=uold(ind_cell,4)
  subsub_dummy%clevel=ind_level


  !!-----
  !! DEBUG
  !!-----
!  skip_loc(1)=dble(icoarse_min)
!  skip_loc(2)=dble(jcoarse_min)
!  skip_loc(3)=dble(kcoarse_min)
!
!  do ilevel=levelmin,nlevelmax-1
!    if(active(ilevel)%ngrid == 0) cycle
!    if(ind_grid .ge. active(ilevel)%igrid(1)) dlev=ilevel
!  enddo
!
!  okay=.false.
!  do ind=1, twotondim
!    x=(xg(ind_grid,1)+xc(ind,1)-skip_loc(1))*scale
!    y=(xg(ind_grid,2)+xc(ind,2)-skip_loc(2))*scale
!    z=(xg(ind_grid,3)+xc(ind,3)-skip_loc(3))*scale
!    dxx=x-xsink(sinkid,1)
!    dyy=y-xsink(sinkid,2)
!    dzz=z-xsink(sinkid,3)
!    drr=MAX(ABS(dxx), ABS(dyy), ABS(dzz))
!
!    if(drr .lt. 0.5d0**dlev) okay=.true.
!  enddo
!
!  if(.not. okay)then
!    write(*,*) 'no good'
!  else
!    write(*,*) 'good'
!  endif

                       
  !!-----
  !! Initialize internal grid
  !!-----
  call subsub_allocate(subsub_dummy)

  subsub_dummy%vg(:,1) = 0.  !! zero velocity IC for test
  subsub_dummy%vg(:,2) = 0.
  subsub_dummy%vg(:,3) = 0.

  subsub_dummy%hydro(:,1) =  max(uold(ind_cell,1), smallr)  !! uniform density IC

  !!-----
  !! Input to subsub_obj array
  !!-----
  call subsub_input(subsub_dummy)

  !!-----
  !! Deallocate
  !!-----
  call subsub_deallocate(subsub_dummy)

end subroutine subsub_create
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