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
  integer :: ismy
  character(LEN=5)::nchar

  integer :: ii, ix, iy, iz
  real(dp) :: rx, ry, rz
  
  !! TODO Get by IO
  if(verbose)write(*,*)'Entering init_subsub'

  !!-----
  call subsub_precision_mpi()
  subsub_boxlen = boxlen/(2D0**nlevelmax)
  subsub_nn = subsub_ngrid**ndim
  subsub_dx = subsub_boxlen/dble(subsub_ngrid)

  !!-----
  !! allocate
  !!-----
  if(allocated(subsub_phi)) deallocate(subsub_phi)
  allocate(subsub_phi(1:subsub_nn))

  if(allocated(subsub_fg)) deallocate(subsub_fg)
  allocate(subsub_fg(1:subsub_nn, 1:ndim))


  !! if memory overhead becomes worse, these can be allocated in the subroutine
  if(allocated(subsub_rho_old)) deallocate(subsub_rho_old)
  allocate(subsub_rho_old(1:subsub_nn))

  !! Density update
  if(allocated(subsub_vg_old)) deallocate(subsub_vg_old)
  allocate(subsub_vg_old(1:subsub_nn, 1:ndim))

  !if(allocated(subsub_vg_up)) deallocate(subsub_vg_up)
  !allocate(subsub_vg_up(1:subsub_nn, 1:ndim))

  !if(allocated(subsub_vg_down)) deallocate(subsub_vg_down)
  !allocate(subsub_vg_down(1:subsub_nn, 1:ndim))

  !if(allocated(subsub_flux_up)) deallocate(subsub_flux_up)
  !allocate(subsub_flux_up(1:subsub_nn, 1:ndim))

  !if(allocated(subsub_flux_down)) deallocate(subsub_flux_down)
  !allocate(subsub_flux_down(1:subsub_nn, 1:ndim))

  !! CG related
  if(allocated(subsub_cgrhs)) deallocate(subsub_cgrhs)
  allocate(subsub_cgrhs(1:subsub_nn))

  if(allocated(subsub_cgLphi)) deallocate(subsub_cgLphi)
  allocate(subsub_cgLphi(1:subsub_nn))

  if(allocated(subsub_cgRes)) deallocate(subsub_cgRes)
  allocate(subsub_cgRes(1:subsub_nn))

  if(allocated(subsub_cgp)) deallocate(subsub_cgp)
  allocate(subsub_cgp(1:subsub_nn))

  if(allocated(subsub_cgLp)) deallocate(subsub_cgLp)
  allocate(subsub_cgLp(1:subsub_nn))

  if(allocated(subsub_dd2)) deallocate(subsub_dd2)
  allocate(subsub_dd2(1:subsub_nn))

  !$omp parallel do collapse(3) private(ii, ix, iy, iz, rx, ry, rz)
  do ix=1, subsub_ngrid
  do iy=1, subsub_ngrid
  do iz=1, subsub_ngrid
    ii = (iz-1)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix
    rx = (dble(ix) - 0.5D0) * subsub_dx
    ry = (dble(iy) - 0.5D0) * subsub_dx
    rz = (dble(iz) - 0.5D0) * subsub_dx

    subsub_dd2(ii) = (rx-subsub_boxlen/2.0D0)**2 + (ry-subsub_boxlen/2.0D0)**2 + (rz-subsub_boxlen/2.0D0)**2
  enddo
  enddo
  enddo
  !$omp end parallel do

  !!-----
  !! Check
  !!-----
  if(subsub_on .and. .not. sink) then
    if(myid .eq. 1) then
      call subsub_log('sink should be .true.', 'init_subsub')
    endif
    call clean_stop
  endif

  allocate(subsub_obj(1:nsinkmax))

  !! Debugger tool (memory clear)
  do i=1, nsinkmax
    !call subsub_allocate(subsub_obj(i))
    call subsub_deallocate(subsub_obj(i))
  enddo

  subsub_end = 0
  subsub_nsink = nsink

  if(.not. hydro)return
  if(nsink.eq.0)return

  !!-----
  !! Input arrays
  !!-----
  if(nrestart.gt.0)then
    !!-----
    !! Retrieve from the restart
    !!-----
    call subsub_readdump()
  else
    !!-----
    !! Create and input obj
    !!-----
    do i=1, nsink
      ismy = mod(abs(idsink(i)), ncpu) + 1
      if(myid.eq.ismy) then
        call subsub_create(i)
      endif
    enddo
  endif


  call subsub_updatedomain




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