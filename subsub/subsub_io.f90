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
    call subsub_allocate(subsub_dump(i))
    subsub_dump(i) = subsub_obj(i)
  enddo
#endif

  !!----- Write
  if(myid .eq. 1) then
    open(newunit=unit_out, file=TRIM(filename), form='unformatted')

    rewind(unit_out)
    ! Write header
    write(unit_out) subsub_nsink
    write(unit_out) subsub_ngrid
    write(unit_out) subsub_nhydro
    write(unit_out) ndim

    do i=1, subsub_nsink
      write(unit_out) subsub_dump(i)%sink_ind          !! sink index
      write(unit_out) subsub_dump(i)%sink_id           !! sink id
      write(unit_out) subsub_dump(i)%sink_mass         !! sink id
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
subroutine subsub_readdump
  use subsub_commons
  use amr_commons
  use pm_commons
  use mpi_mod
  implicit none

  !! Local variables
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::filename,filedir
  integer :: i, j, unit_out, info
  integer :: dumint
  real(dp) :: dumdp, dumdp2(3)
  real(dp), dimension(1:subsub_ngrid**ndim) :: dumdparr
  type(subsub_type), dimension(1:nsink) :: subsub_dummy

  call title(nrestart,nchar)
  if(IOGROUPSIZEREP>0)call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)

  if(IOGROUPSIZEREP>0) then
    filedir='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/'
  else
    filedir='output_'//TRIM(nchar)//'/'
  endif

  !! Read by the host node
  !! ----- RHEE -----
  !! Let's leave reading the entire objects by myid=1 and spreading afterward
  !! If the memory usage becomes worse, read and spread each object one by one
  !! ----------------
  if(myid .eq. 1)then
    filename=TRIM(filedir)//'subsub_'//TRIM(nchar)//'.out'
    open(newunit=unit_out, file=TRIM(filename), form='unformatted')

    read(unit_out) dumint
    subsub_nsink = dumint

    read(unit_out) dumint
    if(dumint .ne. subsub_ngrid) then
      call subsub_log('subsub_ngrid mismatched', 'subsub_readdump')
      write(*,*) 'in namelist: ', subsub_ngrid
      write(*,*) 'in the previous run: ', dumint
      stop
    endif
    !subsub_level = dumint

    read(unit_out) dumint
    if(dumint .ne. subsub_nhydro) then
      call subsub_log('subsub_nhydro mismatched', 'subsub_readdump')
      write(*,*) 'in namelist: ', subsub_nhydro
      write(*,*) 'in the previous run: ', dumint
      stop
    endif

    read(unit_out) dumint
    if(dumint .ne. ndim) then
      call subsub_log('ndim mismatched', 'subsub_readdump')
      write(*,*) 'in namelist: ', ndim
      write(*,*) 'in the previous run: ', dumint
      stop
    endif

    !allocate(subsub_dummy(1:nsink))
    
    do i=1, subsub_nsink
      call subsub_allocate(subsub_dummy(i))

      read(unit_out) dumint
      subsub_dummy(i)%sink_ind = dumint

      read(unit_out) dumint
      subsub_dummy(i)%sink_id = dumint

      read(unit_out) dumdp
      subsub_dummy(i)%sink_mass = dumdp

      read(unit_out) dumdp
      subsub_dummy(i)%mass_tot = dumdp

      read(unit_out) dumdp2
      subsub_dummy(i)%vxc = dumdp2(1)
      subsub_dummy(i)%vyc = dumdp2(2)
      subsub_dummy(i)%vzc = dumdp2(3)

      read(unit_out) dumint
      subsub_dummy(i)%clevel = dumint

      do j=1, ndim
        read(unit_out) dumdparr
        subsub_dummy(i)%vg(:,j) = dumdparr
      enddo

      do j=1, subsub_nhydro
        read(unit_out) dumdparr
        subsub_dummy(i)%hydro(:,j) = dumdparr
      enddo

      !! debugger
      if(idsink(subsub_dummy(i)%sink_ind) .ne. subsub_dummy(i)%sink_id) then
        write(*,*) 'order of sink particles changed?'
        ! if so, spread part should be changed
        stop
      endif
    enddo
    close(unit_out)
  endif

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD, info)
#endif
 
  !!-----
  !! Sort by sink_ind
  !!-----
  call subsub_sortbyind(subsub_dummy, nsink)

  !!-----
  !! Spread
  !!----- 
#ifndef WITHOUTMPI
  call subsub_spread(subsub_dummy)
#else
  do i=1, subsub_nsink
    subsub_obj(i) = subsub_dummy(i)
  enddo
#endif

  !!-----
  !! Deallocate
  !!-----
  if(myid.eq.1) then
    do i=1, subsub_nsink
      call subsub_deallocate(subsub_dummy(i))
    enddo
  endif
  !deallocate(subsub_dummy)
end subroutine
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
  integer :: subsub_nn0, subsub_nn1, subsub_nn2

  integer,parameter::tag1=6421,tag2=6422,tag3=6423
  integer,parameter::tag4=6424,tag5=6425,tag6=6426
  integer,parameter::tag7=6427,tag8=6428,tag9=6429
  integer,parameter::tag10=6430, tag11=6431, tag12=6432


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
  subsub_nn0 = subsub_ngrid**ndim
  subsub_nn1 = (subsub_ngrid**ndim) * ndim
  subsub_nn2 = (subsub_ngrid**ndim) * subsub_nhydro
  iend = subsub_end+1
  if(myid .eq. 1) then
    do i=2, ncpu
      if(subsub_dumpN(i).eq.0) cycle

      do j=1, subsub_dumpN(i)

        call MPI_RECV(subsub_dump(iend)%sink_ind, 1, MPI_INTEGER,  i-1, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_dump(iend)%sink_id,  1, MPI_INTEGER,  i-1, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_dump(iend)%sink_mass,1, subsub_mpidp, i-1, tag3, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_dump(iend)%mass_tot, 1, subsub_mpidp, i-1, tag4, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_dump(iend)%vxc,      1, subsub_mpidp, i-1, tag5, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_dump(iend)%vyc,      1, subsub_mpidp, i-1, tag6, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_dump(iend)%vzc,      1, subsub_mpidp, i-1, tag7, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_dump(iend)%clevel,   1, MPI_INTEGER,  i-1, tag8, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)

        call MPI_RECV(subsub_dump(iend)%vg(1,1), subsub_nn1,    subsub_mpidp, i-1, tag9, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_dump(iend)%hydro(1,1), subsub_nn2, subsub_mpidp, i-1, tag10, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_dump(iend)%phi(1), subsub_nn0,     subsub_mpidp, i-1, tag11, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_dump(iend)%fg(1,1), subsub_nn1,    subsub_mpidp, i-1, tag12, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        iend = iend + 1
      enddo
    enddo
  else
    if(subsub_end .gt. 0) then
      do i=1, subsub_end
        call MPI_SEND(subsub_obj(i)%sink_ind, 1, MPI_INTEGER,  0, tag1, MPI_COMM_WORLD,info)
        call MPI_SEND(subsub_obj(i)%sink_id,  1, MPI_INTEGER,  0, tag2, MPI_COMM_WORLD,info)
        call MPI_SEND(subsub_obj(i)%sink_mass,1, subsub_mpidp, 0, tag3, MPI_COMM_WORLD,info)
        call MPI_SEND(subsub_obj(i)%mass_tot, 1, subsub_mpidp, 0, tag4, MPI_COMM_WORLD,info)
        call MPI_SEND(subsub_obj(i)%vxc,      1, subsub_mpidp, 0, tag5, MPI_COMM_WORLD,info)
        call MPI_SEND(subsub_obj(i)%vyc,      1, subsub_mpidp, 0, tag6, MPI_COMM_WORLD,info)
        call MPI_SEND(subsub_obj(i)%vzc,      1, subsub_mpidp, 0, tag7, MPI_COMM_WORLD,info)
        call MPI_SEND(subsub_obj(i)%clevel,   1, MPI_INTEGER,  0, tag8, MPI_COMM_WORLD,info)

        call MPI_SEND(subsub_obj(i)%vg(1,1), subsub_nn1,    subsub_mpidp, 0, tag9, MPI_COMM_WORLD,info)
        call MPI_SEND(subsub_obj(i)%hydro(1,1), subsub_nn2, subsub_mpidp, 0, tag10, MPI_COMM_WORLD,info)
        call MPI_SEND(subsub_obj(i)%phi(1), subsub_nn0,     subsub_mpidp, 0, tag11, MPI_COMM_WORLD,info)
        call MPI_SEND(subsub_obj(i)%fg(1,1), subsub_nn1,      subsub_mpidp, 0, tag12, MPI_COMM_WORLD,info)
      enddo
    endif
  endif
end subroutine subsub_gather
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_spread(subsub_dummy)
  use subsub_commons
  use amr_commons
  use pm_commons
  use mpi_mod
  implicit none

  type(subsub_type), dimension(1:nsink) :: subsub_dummy

  !! Local variables
  integer :: i, j, ismy, info
  integer :: subsub_nn0, subsub_nn1, subsub_nn2, nx_loc

  integer,parameter::tag1=7421,tag2=7422,tag3=7423
  integer,parameter::tag4=7424,tag5=7425,tag6=7426
  integer,parameter::tag7=7427,tag8=7428,tag9=7429
  integer,parameter::tag10=7430, tag11=7431, tag12=7432

  integer, dimension(1:nsink) :: ismysink_dump, ismysink
  integer :: subsub_sinkinmyid
  real(dp):: scale
  !real(dp):: dx, dx_max

  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)

  ismysink = 0
  ismysink_dump = 0

  !!-----
  !! Check ownwership
  !!-----
  do i=1, nsink
    ismy = -1
    call subsub_finddomain(xsink(i,1)/scale, xsink(i,2)/scale, xsink(i,3)/scale, ismy)
    if(ismy .gt. 0) ismysink_dump(i) = myid
  enddo

  call MPI_ALLREDUCE(ismysink_dump,ismysink,nsink,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)

  !!-----
  !! Spread
  !!-----
  subsub_end = 0
  subsub_nn0 = subsub_ngrid**ndim
  subsub_nn1 = (subsub_ngrid**ndim) * ndim
  subsub_nn2 = (subsub_ngrid**ndim) * subsub_nhydro

  do i=1, nsink
    if(myid .eq. 1) then
      if(ismysink(i).eq.1)then
        subsub_end = subsub_end + 1
        call subsub_allocate(subsub_obj(subsub_end))
        subsub_obj(subsub_end) = subsub_dummy(i)
      else
        call MPI_SEND(subsub_dummy(i)%sink_ind, 1, MPI_INTEGER,  ismysink(i)-1, tag1, MPI_COMM_WORLD,info)
        call MPI_SEND(subsub_dummy(i)%sink_id,  1, MPI_INTEGER,  ismysink(i)-1, tag2, MPI_COMM_WORLD,info)
        call MPI_SEND(subsub_dummy(i)%sink_mass,1, subsub_mpidp, ismysink(i)-1, tag3, MPI_COMM_WORLD,info)
        call MPI_SEND(subsub_dummy(i)%mass_tot, 1, subsub_mpidp, ismysink(i)-1, tag4, MPI_COMM_WORLD,info)
        call MPI_SEND(subsub_dummy(i)%vxc,      1, subsub_mpidp, ismysink(i)-1, tag5, MPI_COMM_WORLD,info)
        call MPI_SEND(subsub_dummy(i)%vyc,      1, subsub_mpidp, ismysink(i)-1, tag6, MPI_COMM_WORLD,info)
        call MPI_SEND(subsub_dummy(i)%vzc,      1, subsub_mpidp, ismysink(i)-1, tag7, MPI_COMM_WORLD,info)
        call MPI_SEND(subsub_dummy(i)%clevel,   1, MPI_INTEGER,  ismysink(i)-1, tag8, MPI_COMM_WORLD,info)
        
        call MPI_SEND(subsub_dummy(i)%vg(1,1), subsub_nn1,    subsub_mpidp, ismysink(i)-1, tag9, MPI_COMM_WORLD,info)
        call MPI_SEND(subsub_dummy(i)%hydro(1,1), subsub_nn2, subsub_mpidp, ismysink(i)-1, tag10, MPI_COMM_WORLD,info)
        call MPI_SEND(subsub_dummy(i)%phi(1), subsub_nn0,     subsub_mpidp, ismysink(i)-1, tag11, MPI_COMM_WORLD,info)
        call MPI_SEND(subsub_dummy(i)%fg(1,1), subsub_nn1,    subsub_mpidp, ismysink(i)-1, tag12, MPI_COMM_WORLD,info)
      endif
    else
      if(ismysink(i).eq.myid) then
        subsub_end = subsub_end + 1
        call subsub_allocate(subsub_obj(subsub_end))

        call MPI_RECV(subsub_obj(subsub_end)%sink_ind, 1, MPI_INTEGER,  0, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_obj(subsub_end)%sink_id,  1, MPI_INTEGER,  0, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_obj(subsub_end)%sink_mass,1, subsub_mpidp, 0, tag3, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_obj(subsub_end)%mass_tot, 1, subsub_mpidp, 0, tag4, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_obj(subsub_end)%vxc,      1, subsub_mpidp, 0, tag5, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_obj(subsub_end)%vyc,      1, subsub_mpidp, 0, tag6, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_obj(subsub_end)%vzc,      1, subsub_mpidp, 0, tag7, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_obj(subsub_end)%clevel,   1, MPI_INTEGER,  0, tag8, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)

        call MPI_RECV(subsub_obj(subsub_end)%vg(1,1), subsub_nn1,    subsub_mpidp, 0, tag9, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_obj(subsub_end)%hydro(1,1), subsub_nn2, subsub_mpidp, 0, tag10, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_obj(subsub_end)%phi(1), subsub_nn0,     subsub_mpidp, 0, tag11, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_obj(subsub_end)%fg(1,1), subsub_nn1,    subsub_mpidp, 0, tag12, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
      endif
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, info)
  enddo

  !!----- Update subsub_nsink
  subsub_sinkinmyid = subsub_end
  subsub_nsink = 0 !! initialize to 0
  call MPI_ALLREDUCE(subsub_sinkinmyid,subsub_nsink,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)

  if(subsub_nsink .ne. nsink) then
    if(myid.eq.1) then
      call subsub_log('some sinks miss in domains', 'subsub_spread')
      write(*,*) ' nsink = ', nsink
      write(*,*) ' updated = ', subsub_nsink
    endif
    call clean_stop
  endif

end subroutine subsub_spread
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_sendtoanother(ind, from, to)
  use subsub_commons
  use amr_commons
  use mpi_mod
  implicit none

  integer :: ind, from, to

  !! Local variables
  integer :: i, j
  integer :: ind_send, info
  integer :: subsub_nn0, subsub_nn1, subsub_nn2

  type(subsub_type) :: subsub_dummy

  integer,parameter::tag1=8421,tag2=8422,tag3=8423
  integer,parameter::tag4=8424,tag5=8425,tag6=8426
  integer,parameter::tag7=8427,tag8=8428,tag9=8429
  integer,parameter::tag10=8430, tag11=8431, tag12=8432

  if(myid .ne. from .and. myid .ne. to) return

  call subsub_allocate(subsub_dummy)

  if(myid .eq. from) then
    do i=1, subsub_end
      if(subsub_obj(i)%sink_ind .eq. ind) then
        ind_send = i
        exit
      endif
    enddo
    subsub_dummy = subsub_obj(ind_send)
  endif

  subsub_nn0 = subsub_ngrid**ndim
  subsub_nn1 = (subsub_ngrid**ndim) * ndim
  subsub_nn2 = (subsub_ngrid**ndim) * subsub_nhydro
  if(myid .eq. from) then
    call MPI_SEND(subsub_dummy%sink_ind, 1, MPI_INTEGER,  to-1, tag1, MPI_COMM_WORLD,info)
    call MPI_SEND(subsub_dummy%sink_id,  1, MPI_INTEGER,  to-1, tag2, MPI_COMM_WORLD,info)
    call MPI_SEND(subsub_dummy%sink_mass,1, subsub_mpidp, to-1, tag3, MPI_COMM_WORLD,info)
    call MPI_SEND(subsub_dummy%mass_tot, 1, subsub_mpidp, to-1, tag4, MPI_COMM_WORLD,info)
    call MPI_SEND(subsub_dummy%vxc,      1, subsub_mpidp, to-1, tag5, MPI_COMM_WORLD,info)
    call MPI_SEND(subsub_dummy%vyc,      1, subsub_mpidp, to-1, tag6, MPI_COMM_WORLD,info)
    call MPI_SEND(subsub_dummy%vzc,      1, subsub_mpidp, to-1, tag7, MPI_COMM_WORLD,info)
    call MPI_SEND(subsub_dummy%clevel,   1, MPI_INTEGER,  to-1, tag8, MPI_COMM_WORLD,info)
        
    call MPI_SEND(subsub_dummy%vg(1,1), subsub_nn1,    subsub_mpidp, to-1, tag9, MPI_COMM_WORLD,info)
    call MPI_SEND(subsub_dummy%hydro(1,1), subsub_nn2, subsub_mpidp, to-1, tag10, MPI_COMM_WORLD,info)
    call MPI_SEND(subsub_dummy%phi(1), subsub_nn0,     subsub_mpidp, to-1, tag11, MPI_COMM_WORLD,info)
    call MPI_SEND(subsub_dummy%fg(1,1), subsub_nn1,    subsub_mpidp, to-1, tag12, MPI_COMM_WORLD,info)
  else
    call MPI_RECV(subsub_dummy%sink_ind, 1, MPI_INTEGER,  from-1, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
    call MPI_RECV(subsub_dummy%sink_id,  1, MPI_INTEGER,  from-1, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
    call MPI_RECV(subsub_dummy%sink_mass,1, subsub_mpidp, from-1, tag3, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
    call MPI_RECV(subsub_dummy%mass_tot, 1, subsub_mpidp, from-1, tag4, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
    call MPI_RECV(subsub_dummy%vxc,      1, subsub_mpidp, from-1, tag5, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
    call MPI_RECV(subsub_dummy%vyc,      1, subsub_mpidp, from-1, tag6, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
    call MPI_RECV(subsub_dummy%vzc,      1, subsub_mpidp, from-1, tag7, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
    call MPI_RECV(subsub_dummy%clevel,   1, MPI_INTEGER,  from-1, tag8, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)

    call MPI_RECV(subsub_dummy%vg(1,1), subsub_nn1,    subsub_mpidp, from-1, tag9, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
    call MPI_RECV(subsub_dummy%hydro(1,1), subsub_nn2, subsub_mpidp, from-1, tag10, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
    call MPI_RECV(subsub_dummy%phi(1), subsub_nn0,     subsub_mpidp, from-1, tag11, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
    call MPI_RECV(subsub_dummy%fg(1,1), subsub_nn1,    subsub_mpidp, from-1, tag12, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
  endif
  
  if(myid .eq. from) then
    if(subsub_end .gt. 1)then
      do i=ind_send, subsub_end-1
        subsub_obj(i) = subsub_obj(i+1)
      enddo
    endif

    call subsub_deallocate(subsub_obj(subsub_end))
    subsub_end = subsub_end - 1
  else if(myid .eq. to) then
    subsub_end = subsub_end + 1

    subsub_obj(subsub_end) = subsub_dummy

    call subsub_sortbyind(subsub_obj(1:subsub_end), subsub_end)
  endif
  
  call subsub_deallocate(subsub_dummy)
end subroutine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_allocate(subsub_dummy)
  use subsub_commons
  use amr_commons
  implicit none
  type(subsub_type) :: subsub_dummy

  !!----- RHEE -----
  !! since in some cases, subsub_dummy is weridly allocated. I don't know why
  !!----------------
  call subsub_deallocate(subsub_dummy)


  allocate(subsub_dummy%vg(1:subsub_ngrid**ndim, 1:ndim))
  allocate(subsub_dummy%hydro(1:subsub_ngrid**ndim, 1:subsub_nhydro))  ! 1 for density
  allocate(subsub_dummy%phi(1:subsub_ngrid**ndim))
  allocate(subsub_dummy%fg(1:subsub_ngrid**ndim, 1:ndim))

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
     
  if(allocated(subsub_dummy%vg)) deallocate(subsub_dummy%vg)
  if(allocated(subsub_dummy%hydro)) deallocate(subsub_dummy%hydro)
  if(allocated(subsub_dummy%phi)) deallocate(subsub_dummy%phi)
  if(allocated(subsub_dummy%fg)) deallocate(subsub_dummy%fg)

end subroutine subsub_deallocate
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_create(sinkind)
  use subsub_commons
  use amr_commons
  use hydro_commons
  use pm_commons
  use mpi_mod
  implicit none
  !!-----
  !! This routine create subsub_obj and input the information of cell that the sink inherited
  !!-----
  integer, intent(in) ::sinkind
  
  

  !! Local Varaibles
  type(subsub_type) :: subsub_dummy
  integer :: i, ilevel, ind, ix, iy, iz, ncache, igrid, ngrid, nx_loc
  real(dp):: scale
  real(dp):: dx

  !integer,dimension(1:nvector)::ind_grid,ind_cell
  integer ::ind_grid, ind_cell, ind_level

  real(dp):: x,y,z, dxx, dyy, dzz, drr
  integer :: dlev
  logical :: okay
  logical :: subsub_flag
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
  subsub_flag = .false.
  call subsub_findcell(xsink(sinkind,1)/scale, xsink(sinkind,2)/scale, xsink(sinkind,3)/scale, &
    ind_cell, ind_grid, ind_level, subsub_flag)

  if(.not. subsub_flag) return !! no matched cell in this domain


  !! if not restart

  !!-----
  !! Initialize object
  !!-----
  call subsub_allocate(subsub_dummy)

  dx=0.5D0**ind_level
  subsub_dummy%sink_ind = sinkind
  subsub_dummy%sink_id  = idsink(sinkind)
  subsub_dummy%sink_mass = msink(sinkind)
  subsub_dummy%mass_tot = (subsub_boxlen**ndim) * max(uold(ind_cell,1), smallr)

  subsub_dummy%vxc=uold(ind_cell,2)
  subsub_dummy%vyc=uold(ind_cell,3)
  subsub_dummy%vzc=uold(ind_cell,4)
  subsub_dummy%clevel=ind_level
  
  subsub_dummy%vg(:,1) = 0.  !! zero velocity IC for test
  subsub_dummy%vg(:,2) = 0.
  subsub_dummy%vg(:,3) = 0.

  subsub_dummy%hydro(:,1) =  max(uold(ind_cell,1), smallr)  !! uniform density IC
  subsub_dummy%phi(:) = 0.

  subsub_dummy%fg(:,1) = 0.
  subsub_dummy%fg(:,2) = 0.
  subsub_dummy%fg(:,3) = 0.

  !!-----
  !! Input to subsub_obj array
  !!-----
  call subsub_input(subsub_dummy)

  !!-----
  !! Deallocate
  !!-----
  call subsub_deallocate(subsub_dummy)


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
!    dxx=x-xsink(sinkind,1)
!    dyy=y-xsink(sinkind,2)
!    dzz=z-xsink(sinkind,3)
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


end subroutine subsub_create