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
  !! TODO Get by IO
  if(verbose)write(*,*)'Entering init_subsub'

  call subsub_precision_mpi()
  subsub_boxlen = boxlen/(2D0**nlevelmax)

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


  allocate(subsub_obj(1:nsinkmax))

  !! Debugger tool (memory clear)
  do i=1, nsinkmax
    call subsub_allocate(subsub_obj(i))
    call subsub_deallocate(subsub_obj(i))
  enddo

  subsub_end = 0

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
      call subsub_create(i)
    enddo
  endif

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