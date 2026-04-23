SUBROUTINE subsub_coarse
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

  real(dp) :: scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  if(myid .eq. 1) then
    write(*,*) ' Time elapsed in SUBSUB_compute [sec]', sngl(subsub_howlong)
    subsub_howlong = 0.0D0
  endif

  if(nsink .eq. 0) return

 


  subsub_v2sink = 0.0D0
  subsub_cs2sink = 0.0D0

  !!----- Update Sink Object
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
#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD, info)
#endif



END SUBROUTINE subsub_coarse