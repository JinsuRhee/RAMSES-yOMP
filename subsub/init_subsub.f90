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

  integer :: ii, ix, iy, iz, ii2
  real(dp) :: rx, ry, rz
  
  !! TODO Get by IO
  if(verbose)write(*,*)'Entering init_subsub'

  !!-----
  call subsub_precision_mpi()
  subsub_boxlen = boxlen/(2D0**nlevelmax)
  subsub_nn = subsub_ngrid**ndim
  subsub_dx = subsub_boxlen/dble(subsub_ngrid)

  subsub_nnface = subsub_nn - (subsub_ngrid-2)**ndim

  subsub_ngrid2 = subsub_ngrid * subsub_ngrid

  subsub_poisson_softening = subsub_dx*0.5D0

  subsub_howlong = 0.0D0
  
  !!-----
  !! allocate
  !!-----
  if(allocated(subsub_dd)) deallocate(subsub_dd)
  allocate(subsub_dd(1:subsub_nn))

  if(allocated(subsub_faceind)) deallocate(subsub_faceind)
  allocate(subsub_faceind(1:subsub_nnface))

  if(allocated(subsub_faceindx)) deallocate(subsub_faceindx)
  allocate(subsub_faceindx(1:subsub_nnface))

  if(allocated(subsub_faceindy)) deallocate(subsub_faceindy)
  allocate(subsub_faceindy(1:subsub_nnface))

  if(allocated(subsub_faceindz)) deallocate(subsub_faceindz)
  allocate(subsub_faceindz(1:subsub_nnface))

  if(allocated(subsub_hydrobc)) deallocate(subsub_hydrobc)
  allocate(subsub_hydrobc(1:2, 1:ndim, 1:subsub_ngrid, 1:subsub_ngrid, 1:subsub_nhydro))

  !subsub_ncloudmax = (2*ir_cloud+1)**ndim
  !allocate(subsub_clouds(1:subsub_ncloudmax, 1:subsub_nhydro+ndim))
  !allocate(subsub_clouds_ind(1:subsub_ncloudmax))
  !allocate(subsub_edgeBC(1:nsink))


  ii2 = 1
  subsub_poisson_softening = subsub_dx
  !$omp parallel do collapse(3) private(ii, ix, iy, iz, rx, ry, rz)
  do ix=1, subsub_ngrid
  do iy=1, subsub_ngrid
  do iz=1, subsub_ngrid
    ii = (iz-1)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix
    rx = (dble(ix) - 0.5D0) * subsub_dx
    ry = (dble(iy) - 0.5D0) * subsub_dx
    rz = (dble(iz) - 0.5D0) * subsub_dx

    subsub_dd(ii) = (rx-subsub_boxlen/2.0D0)**2 + (ry-subsub_boxlen/2.0D0)**2 + (rz-subsub_boxlen/2.0D0)**2 + subsub_poisson_softening**2
    subsub_dd(ii) = sqrt(subsub_dd(ii))

    !! Zero BC
    if(ix.eq. 1 .or. ix .eq. subsub_ngrid .or. &
      iy.eq. 1 .or. iy .eq. subsub_ngrid .or. &
      iz.eq. 1 .or. iz .eq. subsub_ngrid) then

      !$omp critical
      subsub_faceind(ii2) = ii
      subsub_faceindx(ii2) = ix
      subsub_faceindy(ii2) = iy
      subsub_faceindz(ii2) = iz
      ii2 = ii2 + 1
      !$omp end critical
    endif
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