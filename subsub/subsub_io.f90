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
  integer :: i, j, unit_out, icell
  real(dp) :: scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v


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
  !! Constant
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  if(myid .eq. 1) then
    open(newunit=unit_out, file=TRIM(filename), form='unformatted')

    rewind(unit_out)
    ! Write header
    write(unit_out) subsub_nsink
    write(unit_out) subsub_ngrid
    write(unit_out) subsub_nhydro
    write(unit_out) subsub_boxlen
    write(unit_out) ndim
    write(unit_out) aexp
    write(unit_out) scale_l
    write(unit_out) scale_t
    write(unit_out) scale_d
    write(unit_out) scale_v
    write(unit_out) scale_nH
    write(unit_out) scale_T2

    if(subsub_nsink .eq. 0) then
      close(unit_out)
      return
    endif

    do i=1, subsub_nsink
!write(*,*) maxval(subsub_dump(i)%hydro(:,2)), minval(subsub_dump(i)%hydro(:,2))
      write(unit_out) subsub_dump(i)%sink_ind          !! sink index
      write(unit_out) subsub_dump(i)%sink_id           !! sink id
      write(unit_out) subsub_dump(i)%sink_mass         !! sink id
      write(unit_out) subsub_dump(i)%mass_tot          !! total mass
      write(unit_out) subsub_dump(i)%mass_cell         !! cell mass
      do icell=0, twondim
        write(unit_out) subsub_dump(i)%uold(icell,:)              !! uold of the hosting cell
      enddo
      !write(unit_out) subsub_dump(i)%vxc, subsub_dump(i)%vyc, subsub_dump(i)%vzc !! momentum x, y, z
      write(unit_out) subsub_dump(i)%clevel            !! cell level
      write(unit_out) subsub_dump(i)%domain            !! cell level

      !do j=1, ndim
      !  write(unit_out) subsub_dump(i)%vg(:,j)           !! vx, vy, vz
      !enddo

      do j=1, subsub_nhydro
        write(unit_out) subsub_dump(i)%hydro(:,j)
      enddo

      write(unit_out) subsub_dump(i)%phi
      write(unit_out) subsub_dump(i)%phi_bh
    enddo
    

    close(unit_out)
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
  integer :: i, j, unit_out, info, ivar, icell
  integer :: dumint
  real(dp) :: dumdp, dumdp2(subsub_nhydro)
  real(dp), dimension(1:subsub_nn) :: dumdparr
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

    read(unit_out) ! boxlen

    read(unit_out) dumint
    if(dumint .ne. ndim) then
      call subsub_log('ndim mismatched', 'subsub_readdump')
      write(*,*) 'in namelist: ', ndim
      write(*,*) 'in the previous run: ', dumint
      stop
    endif

    read(unit_out) !aexp
    read(unit_out) !scale_l
    read(unit_out) !scale_t
    read(unit_out) !scale_d
    read(unit_out) !scale_v
    read(unit_out) !scale_nH
    read(unit_out) !scale_T2

    !allocate(subsub_dummy(1:nsink))
    if(subsub_nsink .eq. 0) then
      close(unit_out)
      return
    endif

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

      read(unit_out) dumdp
      subsub_dummy(i)%mass_cell = dumdp
  
      do icell=0, twondim
        read(unit_out) dumdp2
        do ivar=1, subsub_nhydro
          subsub_dummy(i)%uold(icell,ivar) = dumdp2(ivar)
        enddo
      enddo

      read(unit_out) dumint
      subsub_dummy(i)%clevel = dumint

      read(unit_out) dumint
      subsub_dummy(i)%domain = dumint
  
      !do j=1, ndim
      !  read(unit_out) dumdparr
      !  subsub_dummy(i)%vg(:,j) = dumdparr
      !enddo
  
      do j=1, subsub_nhydro
        read(unit_out) dumdparr
        subsub_dummy(i)%hydro(:,j) = dumdparr(:)
      enddo

      read(unit_out) dumdparr
      subsub_dummy(i)%phi(:) = dumdparr(:)

      read(unit_out) dumdparr
      subsub_dummy(i)%phi_bh(:) = dumdparr(:)
  
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
  if(myid .eq. 1)then
    do i=1, nsink
      write(*,*) subsub_dummy(i)%sink_id
    enddo
    call subsub_sortbyind(subsub_dummy, nsink)
  endif

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
  use subsub_parameters
  use amr_commons
  use pm_commons
  use mpi_mod
  implicit none

  type(subsub_type), dimension(1:nsink) :: subsub_dump

  ! Local Variables
  integer :: i, j, info, iend, myobj, ivar, icell, icell0
  integer, dimension(1:ncpu) :: subsub_dumpN
  integer :: subsub_nn0, subsub_nn1, subsub_nn2

  integer,parameter::tag1=6421,tag2=6422,tag3=6423
  integer,parameter::tag4=6424,tag5=6425,tag6=6426
  integer,parameter::tag7=6427,tag8=6428,tag9=6429
  integer,parameter::tag10=6430, tag11=6431, tag12=6432

  real(dp), dimension(1:subsub_mpidblpren+subsub_nhydro*(twondim+1)) :: dblarr
  integer,  dimension(1:subsub_mpidblpren+subsub_nhydro*(twondim+1)) :: intarr
  real(dp) :: tcheck(20)

  if(nsink.eq.0) return

  !!----- Allocate first
  if(myid .eq. 1) then
    do i=1, nsink
      call subsub_allocate(subsub_dump(i))
    enddo
  endif

  !!----- Receive
  subsub_nn0 = subsub_nn
  subsub_nn1 = subsub_mpidblpren+subsub_nhydro*(twondim+1)
  subsub_nn2 = (subsub_nn) * subsub_nhydro

  iend = 0
  do i=1, nsink
    myobj = mod(abs(idsink(i)), ncpu) + 1

    if(myobj .eq. 1 .and. myid .eq. 1) then
      iend = iend + 1
      subsub_dump(i) = subsub_obj(iend)
    else
      if(myid .eq. 1) then
        call MPI_RECV(intarr(1), subsub_nn1, MPI_INTEGER,  myobj-1, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(dblarr(1), subsub_nn1, subsub_mpidp, myobj-1, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)

        call MPI_RECV(subsub_dump(i)%hydro(1,1), subsub_nn2, subsub_mpidp, myobj-1, tag3, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_dump(i)%phi(1), subsub_nn0,     subsub_mpidp, myobj-1, tag4, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_dump(i)%phi_bh(1), subsub_nn0,  subsub_mpidp, myobj-1, tag5, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)

        subsub_dump(i)%sink_ind = intarr(1)
        subsub_dump(i)%sink_id = intarr(2)
        subsub_dump(i)%clevel = intarr(3)
        subsub_dump(i)%domain = intarr(4)

        subsub_dump(i)%sink_mass = dblarr(1)
        subsub_dump(i)%mass_tot = dblarr(2)
        subsub_dump(i)%mass_cell = dblarr(3)

        icell0 = subsub_mpidblpren+1
        do icell=0, twondim
          do ivar=1, subsub_nhydro
            subsub_dump(i)%uold(icell,ivar) = dblarr(icell0)
            icell0 = icell0 + 1
          enddo
        enddo

      else if(myid .eq. myobj) then
        iend = iend + 1
        intarr(1) = subsub_obj(iend)%sink_ind
        intarr(2) = subsub_obj(iend)%sink_id
        intarr(3) = subsub_obj(iend)%clevel
        intarr(4) = subsub_obj(iend)%domain

        dblarr(1)   = subsub_obj(iend)%sink_mass
        dblarr(2)   = subsub_obj(iend)%mass_tot
        dblarr(3)   = subsub_obj(iend)%mass_cell

        icell0 = subsub_mpidblpren+1
        do icell=0, twondim
          do ivar=1, subsub_nhydro
            dblarr(icell0) = subsub_obj(iend)%uold(icell,ivar)
            icell0 = icell0 + 1
          enddo
        enddo

        call MPI_SEND(intarr(1), subsub_nn1, MPI_INTEGER,  0, tag1, MPI_COMM_WORLD, info)
        call MPI_SEND(dblarr(1), subsub_nn1, subsub_mpidp, 0, tag2, MPI_COMM_WORLD, info)
        
        call MPI_SEND(subsub_obj(iend)%hydro(1,1), subsub_nn2, subsub_mpidp, 0, tag3, MPI_COMM_WORLD, info)
        call MPI_SEND(subsub_obj(iend)%phi(1), subsub_nn0,     subsub_mpidp, 0, tag4, MPI_COMM_WORLD, info)
        call MPI_SEND(subsub_obj(iend)%phi_bh(1), subsub_nn0,  subsub_mpidp, 0, tag5, MPI_COMM_WORLD, info)
      endif
    endif

    if(iend .gt. subsub_end)then
      call subsub_log('iend exceeds my subsub_end', 'subsub_gather')
      write(*,*) iend, subsub_end, i, idsink(i)
      stop
    endif

    !call MPI_BARRIER(MPI_COMM_WORLD, info)
  enddo
end subroutine subsub_gather
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_spread(subsub_dummy)
  use subsub_commons
  use subsub_parameters
  use amr_commons
  use pm_commons
  use mpi_mod
  implicit none

  type(subsub_type), dimension(1:nsink) :: subsub_dummy

  !! Local variables
  integer :: i, j, ismy, info, icell, icell0, ivar
  integer :: subsub_nn0, subsub_nn1, subsub_nn2, nx_loc

  integer,parameter::tag1=7421,tag2=7422,tag3=7423
  integer,parameter::tag4=7424,tag5=7425,tag6=7426
  integer,parameter::tag7=7427,tag8=7428,tag9=7429
  integer,parameter::tag10=7430, tag11=7431, tag12=7432

  integer, dimension(1:nsink) :: ismysink_dump, ismysink
  integer :: subsub_sinkinmyid
  real(dp):: scale
  real(dp), dimension(1:subsub_mpidblpren+subsub_nhydro*(twondim+1)) :: dblarr
  integer,  dimension(1:subsub_mpidblpren+subsub_nhydro*(twondim+1)) :: intarr
  !real(dp):: dx, dx_max

  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)

  ismysink = 0
  ismysink_dump = 0

  !!-----
  !! Check ownwership
  !!-----
  
  do i=1, nsink
    !ismy = -1
    !call subsub_finddomain(xsink(i,1)/scale, xsink(i,2)/scale, xsink(i,3)/scale, ismy)
    !if(ismy .gt. 0) ismysink_dump(i) = myid
    !ismysink_dump(i) = mod(abs(idsink(i)),ncpu) + 1
    ismysink(i) = mod(abs(idsink(i)),ncpu) + 1
  enddo
  
  !call MPI_ALLREDUCE(ismysink_dump,ismysink,nsink,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)

  !!-----
  !! Spread
  !!-----
  subsub_end = 0
  subsub_nn0 = subsub_nn
  subsub_nn1 = subsub_mpidblpren+subsub_nhydro*(twondim+1)
  subsub_nn2 = (subsub_nn) * subsub_nhydro

  do i=1, nsink
    if(myid .eq. 1) then
      if(ismysink(i).eq.1)then
        subsub_end = subsub_end + 1
        call subsub_allocate(subsub_obj(subsub_end))
        subsub_obj(subsub_end) = subsub_dummy(i)
      else

        intarr(1) = subsub_dummy(i)%sink_ind
        intarr(2) = subsub_dummy(i)%sink_id
        intarr(3) = subsub_dummy(i)%clevel
        intarr(4) = subsub_dummy(i)%domain

        dblarr(1) = subsub_dummy(i)%sink_mass
        dblarr(2) = subsub_dummy(i)%mass_tot
        dblarr(3) = subsub_dummy(i)%mass_cell

        icell0 = subsub_mpidblpren+1
        do icell=0, twondim
          do ivar=1, subsub_nhydro
            dblarr(icell0) = subsub_dummy(i)%uold(icell,ivar)
            icell0 = icell0 + 1
          enddo
        enddo

        call MPI_SEND(intarr(1), subsub_nn1, MPI_INTEGER,  ismysink(i)-1, tag1, MPI_COMM_WORLD, info)
        call MPI_SEND(dblarr(1), subsub_nn1, subsub_mpidp, ismysink(i)-1, tag2, MPI_COMM_WORLD, info)
        
        call MPI_SEND(subsub_dummy(i)%hydro(1,1), subsub_nn2, subsub_mpidp, ismysink(i)-1, tag3, MPI_COMM_WORLD, info)
        call MPI_SEND(subsub_dummy(i)%phi(1), subsub_nn0,     subsub_mpidp, ismysink(i)-1, tag4, MPI_COMM_WORLD, info)
        call MPI_SEND(subsub_dummy(i)%phi_bh(1), subsub_nn0,  subsub_mpidp, ismysink(i)-1, tag5, MPI_COMM_WORLD, info)
        
      endif
    else
      if(ismysink(i).eq.myid) then
        subsub_end = subsub_end + 1
        call subsub_allocate(subsub_obj(subsub_end))

        call MPI_RECV(intarr(1), subsub_nn1, MPI_INTEGER,  0, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(dblarr(1), subsub_nn1, subsub_mpidp, 0, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)

        call MPI_RECV(subsub_obj(subsub_end)%hydro(1,1), subsub_nn2, subsub_mpidp, 0, tag3, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_obj(subsub_end)%phi(1), subsub_nn0,     subsub_mpidp, 0, tag4, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(subsub_obj(subsub_end)%phi_bh(1), subsub_nn0,  subsub_mpidp, 0, tag5, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)

        subsub_obj(subsub_end)%sink_ind = intarr(1)
        subsub_obj(subsub_end)%sink_id = intarr(2)
        subsub_obj(subsub_end)%clevel = intarr(3)
        subsub_obj(subsub_end)%domain = intarr(4)

        subsub_obj(subsub_end)%sink_mass = dblarr(1)
        subsub_obj(subsub_end)%mass_tot = dblarr(2)
        subsub_obj(subsub_end)%mass_cell = dblarr(3)

        icell0 = subsub_mpidblpren+1
        do icell=0, twondim
          do ivar=1, subsub_nhydro
            subsub_obj(subsub_end)%uold(icell,ivar) = dblarr(icell0)
            icell0 = icell0 + 1
          enddo
        enddo        
      endif
    endif
    !call MPI_BARRIER(MPI_COMM_WORLD, info)
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
!subroutine subsub_sendtoanother(ind, from, to)
!  use subsub_commons
!  use amr_commons
!  use mpi_mod
!  implicit none
!
!  integer :: ind, from, to
!
!  !! Local variables
!  integer :: i, j
!  integer :: ind_send, info
!  integer :: subsub_nn0, subsub_nn1, subsub_nn2
!
!  type(subsub_type) :: subsub_dummy
!
!  integer,parameter::tag1=8421,tag2=8422,tag3=8423
!  integer,parameter::tag4=8424,tag5=8425,tag6=8426
!  integer,parameter::tag7=8427,tag8=8428,tag9=8429
!  integer,parameter::tag10=8430, tag11=8431, tag12=8432
!
!  if(myid .ne. from .and. myid .ne. to) return
!
!  call subsub_allocate(subsub_dummy)
!
!  if(myid .eq. from) then
!    do i=1, subsub_end
!      if(subsub_obj(i)%sink_ind .eq. ind) then
!        ind_send = i
!        exit
!      endif
!    enddo
!    subsub_dummy = subsub_obj(ind_send)
!  endif
!
!  subsub_nn0 = subsub_nn
!  subsub_nn1 = (subsub_nn) * ndim
!  subsub_nn2 = (subsub_nn) * subsub_nhydro
!  if(myid .eq. from) then
!    call MPI_SEND(subsub_dummy%sink_ind, 1, MPI_INTEGER,  to-1, tag1, MPI_COMM_WORLD,info)
!    call MPI_SEND(subsub_dummy%sink_id,  1, MPI_INTEGER,  to-1, tag2, MPI_COMM_WORLD,info)
!    call MPI_SEND(subsub_dummy%sink_mass,1, subsub_mpidp, to-1, tag3, MPI_COMM_WORLD,info)
!    call MPI_SEND(subsub_dummy%mass_tot, 1, subsub_mpidp, to-1, tag4, MPI_COMM_WORLD,info)
!    call MPI_SEND(subsub_dummy%vxc,      1, subsub_mpidp, to-1, tag5, MPI_COMM_WORLD,info)
!    call MPI_SEND(subsub_dummy%vyc,      1, subsub_mpidp, to-1, tag6, MPI_COMM_WORLD,info)
!    call MPI_SEND(subsub_dummy%vzc,      1, subsub_mpidp, to-1, tag7, MPI_COMM_WORLD,info)
!    call MPI_SEND(subsub_dummy%clevel,   1, MPI_INTEGER,  to-1, tag8, MPI_COMM_WORLD,info)
!        
!    !call MPI_SEND(subsub_dummy%vg(1,1), subsub_nn1,    subsub_mpidp, to-1, tag9, MPI_COMM_WORLD,info)
!    call MPI_SEND(subsub_dummy%hydro(1,1), subsub_nn2, subsub_mpidp, to-1, tag10, MPI_COMM_WORLD,info)
!    call MPI_SEND(subsub_dummy%phi(1), subsub_nn0,     subsub_mpidp, to-1, tag11, MPI_COMM_WORLD,info)
!    !call MPI_SEND(subsub_dummy%fg(1,1), subsub_nn1,    subsub_mpidp, to-1, tag12, MPI_COMM_WORLD,info)
!  else
!    call MPI_RECV(subsub_dummy%sink_ind, 1, MPI_INTEGER,  from-1, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
!    call MPI_RECV(subsub_dummy%sink_id,  1, MPI_INTEGER,  from-1, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
!    call MPI_RECV(subsub_dummy%sink_mass,1, subsub_mpidp, from-1, tag3, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
!    call MPI_RECV(subsub_dummy%mass_tot, 1, subsub_mpidp, from-1, tag4, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
!    call MPI_RECV(subsub_dummy%vxc,      1, subsub_mpidp, from-1, tag5, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
!    call MPI_RECV(subsub_dummy%vyc,      1, subsub_mpidp, from-1, tag6, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
!    call MPI_RECV(subsub_dummy%vzc,      1, subsub_mpidp, from-1, tag7, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
!    call MPI_RECV(subsub_dummy%clevel,   1, MPI_INTEGER,  from-1, tag8, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
!
!    !call MPI_RECV(subsub_dummy%vg(1,1), subsub_nn1,    subsub_mpidp, from-1, tag9, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
!    call MPI_RECV(subsub_dummy%hydro(1,1), subsub_nn2, subsub_mpidp, from-1, tag10, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
!    call MPI_RECV(subsub_dummy%phi(1), subsub_nn0,     subsub_mpidp, from-1, tag11, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
!    !call MPI_RECV(subsub_dummy%fg(1,1), subsub_nn1,    subsub_mpidp, from-1, tag12, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
!  endif
!  
!  if(myid .eq. from) then
!    if(subsub_end .gt. 1)then
!      do i=ind_send, subsub_end-1
!        subsub_obj(i) = subsub_obj(i+1)
!      enddo
!    endif
!
!    call subsub_deallocate(subsub_obj(subsub_end))
!    subsub_end = subsub_end - 1
!  else if(myid .eq. to) then
!    subsub_end = subsub_end + 1
!
!    subsub_obj(subsub_end) = subsub_dummy
!
!    call subsub_sortbyind(subsub_obj(1:subsub_end), subsub_end)
!  endif
!  
!  call subsub_deallocate(subsub_dummy)
!end subroutine
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


  !allocate(subsub_dummy%vg(1:subsub_ngrid**ndim, 1:ndim))
  allocate(subsub_dummy%hydro(1:subsub_ngrid**ndim, 1:subsub_nhydro))  ! 1 for density
  allocate(subsub_dummy%phi(1:subsub_ngrid**ndim))
  allocate(subsub_dummy%phi_bh(1:subsub_ngrid**ndim))
  allocate(subsub_dummy%uold(0:twondim, 1:subsub_nhydro))
  allocate(subsub_dummy%edgeBC(1:twotondim,1:subsub_nhydro))
  allocate(subsub_dummy%faceBC(1:2,1:ndim,1:subsub_ngrid,1:subsub_ngrid,1:subsub_nhydro)) ! face direction, dimension, index1, index2, nvar
  

  subsub_dummy%new = .true.
  subsub_dummy%hydro(:,:) = 0.0D0
  subsub_dummy%phi(:) = 0.0D0
  subsub_dummy%phi_bh(:) = 0.0D0
  subsub_dummy%uold(:,:) = 0.0D0
  subsub_dummy%edgeBC(:,:) = 0.0D0
  subsub_dummy%faceBC(:,:,:,:,:) = 0.0D0
  subsub_dummy%v2 = -1.0D0 !! initially negative (working as a flag)
  subsub_dummy%cs2 = -1.0D0

end subroutine subsub_allocate
subroutine subsub_deallocate(subsub_dummy)
  use subsub_commons
  implicit none
  type(subsub_type) :: subsub_dummy

  subsub_dummy%sink_ind = 0
  subsub_dummy%sink_id = 0
  subsub_dummy%sink_mass = 0.0D0
  subsub_dummy%mass_tot = 0.0D0
  subsub_dummy%mass_cell = 0.0D0
  subsub_dummy%clevel = 0
  subsub_dummy%domain = 0
  subsub_dummy%v2 = 0.0D0
  subsub_dummy%cs2 = 0.0D0
     
  if(allocated(subsub_dummy%hydro)) deallocate(subsub_dummy%hydro)
  if(allocated(subsub_dummy%phi)) deallocate(subsub_dummy%phi)
  if(allocated(subsub_dummy%phi_bh)) deallocate(subsub_dummy%phi_bh)
  if(allocated(subsub_dummy%uold)) deallocate(subsub_dummy%uold)
  if(allocated(subsub_dummy%edgeBC)) deallocate(subsub_dummy%edgeBC)
  if(allocated(subsub_dummy%faceBC)) deallocate(subsub_dummy%faceBC)

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
  integer, intent(in) :: sinkind
  
  

  !! Local Varaibles
  !type(subsub_type) :: subsub_dummy
  !integer :: i, ilevel, ind, ix, iy, iz, ncache, igrid, ngrid, nx_loc
  !real(dp):: scale
  !real(dp):: dx

  !integer,dimension(1:nvector)::ind_grid,ind_cell
  !integer ::ind_grid, ind_cell, ind_level

  !real(dp):: x,y,z, dxx, dyy, dzz, drr
  !integer :: dlev
  !logical :: okay
  !logical :: subsub_flag
  real(dp),dimension(1:3)::xbound,skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc

  !!-----
  !! Initialize
  !!-----

  !!-----
  !! Retrieve related properties from the inherited cell
  !!-----
  !subsub_flag = .false.
  !call subsub_findcell(xsink(sinkind,1)/scale, xsink(sinkind,2)/scale, xsink(sinkind,3)/scale, &
  !  ind_cell, ind_grid, ind_level, subsub_flag)
  !if(.not. subsub_flag) return !! no matched cell in this domain

  !!-----
  !! Initialize object
  !!-----
  subsub_end = subsub_end + 1

  if(subsub_end .gt. size(subsub_obj)) then
   if(myid .eq. 1) call subsub_log('subsub_object exceeds nsinkmax', 'subsub_create')
   call clean_stop
  endif

  call subsub_allocate(subsub_obj(subsub_end))

  subsub_obj(subsub_end)%sink_ind = sinkind
  subsub_obj(subsub_end)%sink_id  = idsink(sinkind)
  subsub_obj(subsub_end)%sink_mass = msink(sinkind)

end subroutine subsub_create
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_updatedomain(ilevel)
  !!-----
  !! This routine updates the domain number for the sink particles and the uold of the cell hosting sink particles
  use amr_commons
  use subsub_commons
  use subsub_parameters
  use pm_commons
  use mpi_mod
  use hydro_commons
  implicit none

  integer :: ilevel

  !! Local variables
  integer :: i, j, ismy, info, myind, myobj, nx_loc, ivar
  integer :: icell, icell0
  integer :: ind_cell, ind_grid, ind_level
  integer, dimension(1:nsink) :: mysink, mysink_dump
  real(dp) :: scale
  logical :: subsub_flag
  integer :: tag1=9421, tag2=9422

  real(dp), dimension(0:twondim, 1:subsub_nhydro) :: hdummy
  real(dp), dimension(1:subsub_nhydro*(twondim+1)) :: dblarr
  real(dp), dimension(1:twotondim, 1:subsub_nhydro) :: edgeBC
  integer, dimension(1:1) :: intarr

  integer :: ix, iy, iz, ii
  real(dp) :: r, u, v, w
  real(dp) :: xc, yc, zc, rvx, rvy, rvz, rvv, pold, ekin, v2, cs2, ekin_new, rhotot

  !real(dp), allocatable(1:nsink) :: v2_local, cs2_local


  if(nsink .eq. 0) return

  !!
!  allocate(v2_local(1:nsink))
!  allocate(cs2_local(1:nsink))
!
!  v2_local(:) = 0.0D0
!  cs2_local(:) = 0.0D0

  !! constants
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)

  !subsub_ngrid2 = subsub_ngrid * subsub_ngrid

  !! find ownership
#ifndef WITHOUTMPI
  mysink(:) = 0
  mysink_dump(:) = 0
  
  do i=1, nsink
    ismy = -1
    call subsub_finddomain(xsink(i,1)/scale, xsink(i,2)/scale, xsink(i,3)/scale, ismy)
    if(ismy .gt. 0) mysink_dump(i) = myid
  enddo

  call MPI_ALLREDUCE(mysink_dump,mysink,nsink,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
  mysink(:) = myid
#endif

  !! send/recv uold
  myind = 0
  do i=1, nsink

#ifndef WITHOUTMPI
    myobj = mod(abs(idsink(i)), ncpu) + 1
#else
    myobj = 1
#endif

    !! If this sink is in my domain, (mysink.eq.myid) retrieve cell information and send it to the MPI who has this sink object (myobj)
    !! If this sink is mine (mysink.eq.myid) and allocated to me (myobj .eq. myid) do not send/recv
    if(myobj .eq. myid) then
      myind = myind + 1

      if(mysink(i) .eq. myid) then
        subsub_flag = .false.
        call subsub_findcell(xsink(i,1)/scale, xsink(i,2)/scale, xsink(i,3)/scale, &
          ind_cell, ind_grid, ind_level, subsub_flag)

        if(.not. subsub_flag)then
          call subsub_log('cannot find the cell', 'subsub_updatedomain')
          stop
        endif

        !call subsub_computeBC(i, edgeBC)

        subsub_obj(myind)%uold(0,:) = uold(ind_cell,1:subsub_nhydro)
        !do icell=1, twondim
        !  do ivar=1, subsub_nhydro
        !    subsub_obj(myind)%edgeBC(icell,ivar) = edgeBC(icell, ivar)
        !  enddo
        !enddo
        !! initial mass by the total cell mass

        subsub_obj(myind)%mass_tot = (subsub_boxlen**ndim) * max(uold(ind_cell,1), subsub_dfloor)
        subsub_obj(myind)%mass_cell = (subsub_boxlen**ndim) * max(uold(ind_cell,1), subsub_dfloor)
        subsub_obj(myind)%clevel = ind_level
        subsub_obj(myind)%domain = myid


      else
        call MPI_RECV(dblarr(1), subsub_nhydro, subsub_mpidp, mysink(i)-1, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)
        call MPI_RECV(intarr(1), 1, MPI_INTEGER,  mysink(i)-1, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, info)

        !icell0 = 1
        !do icell=1, twondim
        !  do ivar=1, subsub_nhydro
        !    subsub_obj(myind)%edgeBC(icell,ivar) = dblarr(icell0)
        !    icell0 = icell0 + 1
        !  enddo
        !enddo

        do ivar=1, subsub_nhydro
          subsub_obj(myind)%uold(0,ivar) = dblarr(ivar)
        enddo

        
        subsub_obj(myind)%mass_tot = (subsub_boxlen**ndim) * max(subsub_obj(myind)%uold(0,1), subsub_dfloor)
        subsub_obj(myind)%mass_cell = (subsub_boxlen**ndim) * max(subsub_obj(myind)%uold(0,1), subsub_dfloor)
        subsub_obj(myind)%clevel = intarr(1)
        subsub_obj(myind)%domain = mysink(i)

      endif
    else
      if(mysink(i) .eq. myid) then
        subsub_flag = .false.
        call subsub_findcell(xsink(i,1)/scale, xsink(i,2)/scale, xsink(i,3)/scale, &
          ind_cell, ind_grid, ind_level, subsub_flag)

        if(.not. subsub_flag)then
          call subsub_log('cannot find the cell', 'subsub_updatedomain2')
          stop
        endif
             
        !call subsub_nbcellinput(ind_cell, ind_level, hdummy)
        !call subsub_computeBC(i, edgeBC)
        do ivar=1, subsub_nhydro
          dblarr(ivar) = uold(ind_cell,ivar)
        enddo

        !icell0 = 1
        !do icell=1, twondim
        !  do ivar=1, subsub_nhydro
        !    dblarr(icell0) = edgeBC(icell, ivar)
        !    icell0 = icell0 + 1
        !  enddo
        !enddo
        intarr(1) = ind_level


        call MPI_SEND(dblarr(1), subsub_nhydro, subsub_mpidp, myobj-1, tag1, MPI_COMM_WORLD,info)
        call MPI_SEND(intarr(1), 1, MPI_INTEGER,  myobj-1, tag2, MPI_COMM_WORLD,info)
      endif
    endif 
  enddo

  !!----- Set Cells
  if(subsub_end .gt. 0) then
    do myind=1, subsub_end

      !!----- Initialize for newly allocated
      if(subsub_obj(myind)%new) then 

        !v2 = 0.0D0
        !cs2 = 0.0D0
        !rhotot = 0.0D0
        !$omp parallel do private(pold, ekin, ekin_new)
        do j=1, subsub_nn
          subsub_obj(myind)%hydro(j,1) = max(subsub_obj(myind)%uold(0,1), subsub_dfloor)
          subsub_obj(myind)%hydro(j,2) = 0.0D0!dblarr(2)
          subsub_obj(myind)%hydro(j,3) = 0.0D0!dblarr(3)
          subsub_obj(myind)%hydro(j,4) = 0.0D0!dblarr(4)

          ekin = 0.5D0 * (subsub_obj(myind)%uold(0,2)**2 + subsub_obj(myind)%uold(0,3)**2 + subsub_obj(myind)%uold(0,4)**2) / subsub_obj(myind)%uold(0,1)
          ekin_new = 0.5D0 * (subsub_obj(myind)%hydro(j,2)**2 + subsub_obj(myind)%hydro(j,3)**2 + subsub_obj(myind)%hydro(j,4)**2) / subsub_obj(myind)%hydro(j,1)

          pold = max((subsub_obj(myind)%uold(0,5) - ekin) * (gamma - 1.0D0), subsub_pfloor)
          subsub_obj(myind)%hydro(j,5) = pold/(gamma-1.0D0) + ekin_new
            

          !! Mass-weighted v2 and cs2
          !rhotot = rhotot + subsub_obj(myind)%hydro(j,1)
          !v2 = v2 + ekin_new
          !cs2 = cs2 + pold / subsub_obj(myind)%hydro(j,1) * subsub_obj(myind)%hydro(j,1)
        enddo
        !$omp end parallel do


        !subsub_obj(myind)%v2 = v2 / rhotot
        !subsub_obj(myind)%cs2 = cs2 / rhotot
        !v2_local(subsub_dummy%sink_ind) = subsub_obj(myind)%v2
        !cs2_local(subsub_dummy%sink_ind) = subsub_obj(myind)%cs2


  
        subsub_obj(myind)%mass_tot = (subsub_boxlen**ndim) * max(dblarr(1), subsub_dfloor)
        subsub_obj(myind)%mass_cell = (subsub_boxlen**ndim) * max(dblarr(1), subsub_dfloor)
        subsub_obj(myind)%sink_mass = msink(subsub_obj(myind)%sink_ind)
  
        !!----- Some Debugger
  
        !! DEBUG MODE FOR SELF-GRAVITY TEST
        !! )) DEBUGG GRAV((         <- this is for grep
        if(subsub_dev_gravonly .eq. 1)then
          !$omp parallel do private(ix,iy,iz, xc, yc, zc, rvx, rvy, rvz, rvv)
          do j=1, subsub_nn
            subsub_obj(myind)%hydro(j,2) = 0.0D0
            subsub_obj(myind)%hydro(j,3) = 0.0D0
            subsub_obj(myind)%hydro(j,4) = 0.0D0
            subsub_obj(myind)%hydro(j,5) = 0.0D0
          enddo
          !$omp end parallel do
        endif
  
        !! DEBUG MODE FOR SELF-GRAVITY TEST
        !! )) DEBUGG HYDRO((         <- this is for grep
        if(subsub_dev_hydroonly .eq. 1)then
          !$omp parallel do private(ix,iy,iz, xc, yc, zc, rvx, rvy, rvz, rvv)
          do j=1, subsub_nn
            call subsub_get3ind(j, ix, iy, iz)
  
            if(iy.ge.subsub_ngrid/2)then
              subsub_obj(myind)%hydro(j,1) = 1.0D-1
              subsub_obj(myind)%hydro(j,2) = 1.0D-3*subsub_obj(myind)%hydro(j,1)
              subsub_obj(myind)%hydro(j,3) = 0.
              subsub_obj(myind)%hydro(j,4) = 0.
              subsub_obj(myind)%hydro(j,5) = 1.0D-10
            else
              subsub_obj(myind)%hydro(j,1) = 5.0D-1
              subsub_obj(myind)%hydro(j,2) = 1.0D-4*subsub_obj(myind)%hydro(j,1)
              subsub_obj(myind)%hydro(j,3) = 0.
              subsub_obj(myind)%hydro(j,4) = 0.
              subsub_obj(myind)%hydro(j,5) = 1.0D-10
            endif
          enddo
          !$omp end parallel do
        endif
  
        !! DEBUG MODE FOR SELF-GRAVITY TEST
        !! )) DEBUGG GRAV+HYDRO((         <- this is for grep
        if(subsub_dev_gravhydro .eq. 1)then
          !$omp parallel do private(ix,iy,iz, xc, yc, zc, rvx, rvy, rvz, rvv)
          do j=1, subsub_nn
            subsub_obj(myind)%hydro(j,1) = dblarr(1)
            subsub_obj(myind)%hydro(j,2) = 0.0D0
            subsub_obj(myind)%hydro(j,3) = 0.0D0
            subsub_obj(myind)%hydro(j,4) = 0.0D0
            subsub_obj(myind)%hydro(j,5) = subsub_pfloor/(gamma-1.0D0)
          enddo
          !$omp end parallel do
        endif
  
        !! DEBUG MODE FOR SELF-GRAVITY TEST
        !! )) DEBUGG GRAV+HYDRO+BH((         <- this is for grep 
        if(subsub_dev_gravhydrobh .eq. 1)then
          !$omp parallel do private(ix,iy,iz, xc, yc, zc, rvx, rvy, rvz, rvv)
          do j=1, subsub_nn
            call subsub_get3ind(j, ix, iy, iz)
  
            subsub_obj(myind)%hydro(j,1) = subsub_obj(myind)%uold(0,1)
  
            xc = (dble(ix)-0.5D0)*subsub_dx - subsub_boxlen/2.0D0
            yc = (dble(iy)-0.5D0)*subsub_dx - subsub_boxlen/2.0D0
            zc = (dble(iz)-0.5D0)*subsub_dx - subsub_boxlen/2.0D0
            rvv = 1.0D-6 * MAX((subsub_boxlen/2.0D0*sqrt(3.D0) - subsub_dd(j)),0.0D0)/subsub_boxlen/2.0D0
  
            rvx = sqrt(rvv / (1.0D0 + (xc/yc)**2))
            if(yc.lt.0) rvx = -rvx
            rvy = rvx*(-xc/yc)
            rvz = 0.0D0
  
            !write(*,*) vsink(subsub_obj(myind)%sink_ind,1)
            !write(*,*) subsub_obj(myind)%uold(0,2)/subsub_obj(myind)%uold(0,1)
            !write(*,*) vsink(subsub_obj(myind)%sink_ind,1)-subsub_obj(myind)%uold(0,2)/subsub_obj(myind)%uold(0,1)
  
            subsub_obj(myind)%hydro(j,2) = rvx * subsub_obj(myind)%hydro(j,1)
            subsub_obj(myind)%hydro(j,3) = rvy * subsub_obj(myind)%hydro(j,1)
            subsub_obj(myind)%hydro(j,4) = rvz * subsub_obj(myind)%hydro(j,1)
            subsub_obj(myind)%hydro(j,5) = subsub_pfloor/(gamma-1.0D0)
          enddo
          !$omp end parallel do
        endif
  
  
        !! DEBUG MODE FOR SELF-GRAVITY TEST
        !! )) DEBUGG SPHERICAL IC((         <- this is for grep 
        if(subsub_dev_icsphere .eq. 1) then
          !!----- RHEE -----
          !! For spherical IC test
          !!----------------
          !$omp parallel do collapse(2) private(iz, ii)
          do ix=1, subsub_ngrid
          do iy=1, subsub_ngrid
          do iz=1, subsub_ngrid
            ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
            if(ix**2 + iy**2 + iz**2 .gt. subsub_ngrid**2) subsub_obj(myind)%hydro(ii,1) = subsub_dfloor
          enddo
          enddo
          enddo
          !$omp end parallel do
        endif
      endif

      subsub_obj(myind)%new = .false.
    enddo
  endif

  !! Communicate sink arrays
!#ifndef WITHOUTMPI
!  call MPI_ALLREDUCE(v2_local(1), subsub_v2sink(1), nsink, subsub_mpidp, MPI_SUM, MPI_COMM_WORLD, info)
!  call MPI_ALLREDUCE(cs2_local(1), subsub_cs2sink(1), nsink, subsub_mpidp, MPI_SUM, MPI_COMM_WORLD, info)
!#else
!  subsub_v2sink(:) = v2_local(:)
!  subsub_cs2sink(:) = subsub_cs2sink(:)
!#endif

!  deallocate(v2_local)
!  deallocate(cs2_local)
end subroutine subsub_updatedomain
!################################################################
!################################################################
!################################################################
!################################################################