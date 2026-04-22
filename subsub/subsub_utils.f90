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
recursive subroutine subsub_walkgrid(xc, skip_loc, sinkid, ind0, ilevel0, x, y, z, r, clouds, n, nmax)
  use amr_commons
  use pm_commons
  use hydro_commons
  use subsub_parameters
  use hydro_parameters, ONLY: gamma
  implicit none

  integer :: ind0, ilevel0, sinkid, nmax
  real(dp), dimension(1:nmax,1:subsub_nhydro+ndim) :: clouds
  integer :: n
  real(dp) :: x, y, z, r
  real(dp),dimension(1:twotondim,1:3) :: xc
  real(dp),dimension(1:3)::skip_loc

  !! Local variables
  integer :: ind, ind2, idim, ilevelc
  integer :: ind_level, index_grid, index_cell, iskip
  integer :: index_grid2
  real(dp) :: dx, dr_cell
  integer :: ix, iy, iz, ipart, jpart

  
  real(dp),dimension(1:twotondim, 1:3) :: center

  integer :: npart1, next_part
  logical :: okay, okay2, okay3
  real(dp) :: rho, u, v, w, p, ekin

  ilevelc = ilevel0 + 1

  if(ilevelc .gt. nlevelmax-nlevelsheld) return
  index_grid = ind0


if(myid.eq.26) then
  write(*,*) ind0, ilevelc
endif

  !! Get cell center positons
  dx=0.5D0**ilevelc
  do ind=1, twotondim
    do idim=1, ndim
      center(ind, idim) = xg(index_grid,idim)+xc(ind,idim)*dx-skip_loc(idim)
    enddo
  enddo

  okay2 = .false.
  do ind=1, twotondim

    iskip = ncoarse + (ind-1)*ngridmax
    index_cell = iskip + index_grid

    okay = .true.
    if (center(ind,1)+dx/2.d0 < x-r .or. center(ind,1)-dx/2.d0 > x+r) okay = .false.
    if (center(ind,2)+dx/2.d0 < y-r .or. center(ind,2)-dx/2.d0 > y+r) okay = .false.
    if (center(ind,3)+dx/2.d0 < z-r .or. center(ind,3)-dx/2.d0 > z+r) okay = .false.

if(myid.eq.26) then
  write(*,*) ilevelc, ind, index_grid
  if(okay)then
    write(*,*) ' --- okay here', son(index_cell)
    write(*,*)
  endif
endif

    !! One of the edges is located in this cell
    if(okay) then
      index_grid2 = son(index_cell)

      if(index_grid2 .eq. 0) then !! This is a leaf cell. Collect coulds.
        okay2 = .true.
        exit
      else
        call subsub_walkgrid(xc, skip_loc, sinkid, index_grid2, ilevelc, x, y, z, r, clouds, n, nmax)

      endif
    endif
  enddo

if(myid.eq.26) then
  write(*,*) 'found ??'
  if(okay2) write(*,*) ' -> yes'
  if(okay2) write(*,*) index_grid, numbp(index_grid), ilevelc
endif

  !! This grid has leaf cells contained to the box
  if(okay2) then
    npart1=numbp(index_grid)
    if(npart1 .le. 0) return

    ipart = headp(index_grid)
    do jpart=1, npart1
      next_part=nextp(ipart)

      !! cloud for this sink
      if(is_cloud(typep(ipart)) .and. idp(ipart) .eq. -sinkid) then

        n = n + 1
        
        !! Contained to which cell?
       do ind2=1, twotondim
          if(xp(ipart,1) .lt. center(ind2,1)-dx/2.0D0 .or. xp(ipart,1) .gt. center(ind2,1)+dx/2.0D0) cycle
          if(xp(ipart,2) .lt. center(ind2,2)-dx/2.0D0 .or. xp(ipart,2) .gt. center(ind2,2)+dx/2.0D0) cycle
          if(xp(ipart,3) .lt. center(ind2,3)-dx/2.0D0 .or. xp(ipart,3) .gt. center(ind2,3)+dx/2.0D0) cycle
         

          iskip = ncoarse + (ind2-1)*ngridmax
          index_cell = iskip + index_grid

          if(n.gt.nmax) then
            write(*,*) 'more cloud particles than maximum?'
            write(*,*) 'n_cloud_max = ', nmax
            write(*,*) 'myid = ', myid
            write(*,*) 'sinkid = ', sinkid
            stop
          endif

          !! Save as primitive
          rho = MAX(uold(index_cell,1), subsub_dfloor)
          u = uold(index_cell,2)/rho
          v = uold(index_cell,3)/rho
          w = uold(index_cell,4)/rho
          ekin = 0.5D0 * (u**2 + v**2 + w**2) * rho
          p = (uold(index_cell,5) - ekin)*(gamma-1.0D0)

          clouds(n,1) = rho
          clouds(n,2) = u - vp(ipart,1)
          clouds(n,3) = v - vp(ipart,2)
          clouds(n,4) = w - vp(ipart,3)
          clouds(n,5) = p

          clouds(n,6) = xp(ipart,1) - x
          clouds(n,7) = xp(ipart,2) - y
          clouds(n,8) = xp(ipart,3) - z
          exit
        enddo


      endif
      ipart=next_part  ! Go to next particle
    enddo
  endif
end
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_test
  use subsub_parameters
  use subsub_commons
  use amr_commons
  use pm_commons
  use mpi_mod
  implicit none
  logical :: okay2
  integer :: nnn, ind_cell, ncache, igrid, ngrid, index_grid, npart1, ipart, jpart, next_part
  real(dp) :: rrx, rry, rrz
  integer :: i,j,k, ilevel, l , m , n, ii, jj, kk, pp
  integer :: ll2, mm2, nn2, ii2, jj2, kk2, pp2

  if(subsub_end .gt. 0) then
    do i=1, subsub_end
      if(subsub_obj(i)%sink_id .eq. 3) write(*,*) 'sink 3 in ', myid
    enddo
  endif
 
  if(myid.eq.1) then
    do i=1, nsink
      write(*,*) ' sink ind & id ' , idsink(i), i
    enddo
  endif
  call MPI_BARRIER(MPI_COMM_WORLD, nnn)

  okay2=.false.
  nnn = 0
  rrx = 0.
  rry = 0.
  rrz = 0.
  j = 0
  k = 0
  pp = 0


ipart = 0
do i = 1, npartmax
  if (levelp(i) > 0) then
    ipart = ipart+1

    if(is_cloud(typep(i))) k = k + 1

    if(is_cloud(typep(i)) .and. idp(i) .eq. -3) pp = pp + 1
  endif
enddo

l = ipart
m = k

k = 0


do ilevel=levelmin, nlevelmax
  igrid = headl(myid, ilevel)
  if(numbl(myid, ilevel).le.0)cycle
  do i=1, numbl(myid, ilevel)
    npart1 = numbp(igrid)
    j = j + npart1
    if(npart1 .eq. 0) cycle

    ipart = headp(igrid)
    do jpart=1, npart1
      if(is_cloud(typep(ipart))) k = k + 1
      if(is_cloud(typep(ipart)) .and. idp(ipart) .eq. -3) then
        nnn = nnn + 1
        rrx = rrx + xp(ipart,1)
        rry = rry + xp(ipart,2)
        rrz = rrz + xp(ipart,3)
        ind_cell = igrid
      endif
      next_part = nextp(ipart)
      ipart = next_part
    enddo

    igrid = next(igrid)
  enddo
enddo

ii = 0
jj = 0
kk = 0
do ilevel=levelmin, nlevelmax
  ncache = active(ilevel)%ngrid
  do igrid=1, ncache
    index_grid = active(ilevel)%igrid(igrid)

    npart1 = numbp(index_grid)
    if(npart1 .le. 0) cycle
    ii = ii + npart1

    ipart = headp(index_grid)
    do jpart=1, npart1
      if(is_cloud(typep(ipart))) jj = jj + 1
      if(is_cloud(typep(ipart)) .and. idp(ipart) .eq. -3) kk = kk + 1
      ipart = nextp(ipart)
    enddo
  enddo
enddo
!  ind_cell = 0
!  ncache = active(levelmin)%ngrid
!  do igrid=1, ncache, nvector
!    ngrid = MIN(nvector,ncache-igrid+1)
!    do i=1, ngrid
!      index_grid=active(levelmin)%igrid(igrid+i-1)
!
!      npart1 = numbp(index_grid)
!      j = j + npart1
!      if(npart1 .eq. 0) cycle
!      ipart = headp(index_grid)
!      do jpart=1, npart1
!        if(is_cloud(typep(ipart))) k = k + 1
!        if(is_cloud(typep(ipart)) .and. idp(ipart) .eq. -3) then
!          nnn = nnn + 1
!          rrx = rrx + xp(ipart,1)
!          rry = rry + xp(ipart,2)
!          rrz = rrz + xp(ipart,3)
!          ind_cell = index_grid
!        endif
!        next_part = nextp(ipart)
!        ipart=next_part
!      enddo
!    enddo
!  enddo
!!----- Update subsub_nsink

!!----- Update subsub_nsink
  call MPI_ALLREDUCE(l,ll2,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,i)
  call MPI_ALLREDUCE(ii,ii2,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,i)
  call MPI_ALLREDUCE(m,mm2,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,i)
  call MPI_ALLREDUCE(jj,jj2,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,i)
  call MPI_ALLREDUCE(kk,kk2,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,i)
  call MPI_ALLREDUCE(pp,pp2,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,i)

  if(myid.eq.1) then
    write(*,*) 'total particles ', ll2, ii2
    write(*,*) 'total clouds ', mm2, jj2
    write(*,*) ' 2109 ?', pp2, kk2
  endif


call MPI_BARRIER(MPI_COMM_WORLD, nnn)
do i=1, ncpu
  if(myid .eq. i) then
    write(*,*) 'statistics', myid
    write(*,*) ' by part ', l, ' by grid ', j, ' by_grid2 ', ii
    write(*,*) 'C by part ', m, ' by grid', k, ' by_grid2 ', jj
    write(*,*) '  '
  endif
  call MPI_BARRIER(MPI_COMM_WORLD, nnn)
enddo
  !write(*,*) 'statistics : ', myid, j,k
  !if(nnn .gt. 0)then
  !write(*,*) '-- ', ind_cell
  !write(*,*) '-- ', nnn, myid
  !write(*,*) '-- ', rrx/nnn, rry/nnn, rrz/nnn
  !endif

  call MPI_BARRIER(MPI_COMM_WORLD, nnn)

  call clean_stop
end subroutine
subroutine subsub_getcloud(sinkid, x, y, z, clouds, ncloud, nmax)

  use amr_commons
  use pm_commons
  use pm_parameters
  use subsub_parameters
  use mpi_mod
  implicit none

  integer :: sinkid, nmax, ncloud
  real(dp) :: x, y, z
  !real(dp), dimension(1:twotondim, 1:ndim), intent(in) :: edge
  real(dp), dimension(1:nmax,1:subsub_nhydro+ndim) :: clouds 


  !! Local variables
  integer :: i, j, k, idim, igrid
  integer :: ind, ix, iy, iz
  real(dp) :: dx_min, rmax, dx, dr_cell
  integer :: ncache
  integer :: found, ind_cell, ind_grid, ngrid, index_grid, index_cell, iskip
  real(dp),dimension(1:twotondim,1:3):: xc
  real(dp),dimension(1:3)::skip_loc, center, tcenter

  logical :: okay, okay2

integer :: npart1, ipart, jpart, next_part, ilevel, nnn
real(dp) :: rrx, rry, rrz

  !!----- Pre
  dx_min=0.5d0**(nlevelmax-nlevelsheld)/aexp
  rmax=dble(ir_cloud)*dx_min

  tcenter = (/x, y, z/)
  !!----- Get Coarse First
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

  found = 0
  ind_cell = -1
  ind_grid = -1

  !$omp parallel do default(shared) &
  !$omp & private(ngrid, i, idim, index_grid, index_cell, &
  !$omp & ind, iskip, center, dr_cell) &
  !$omp & schedule(dynamic)
  do igrid=1,ncache,nvector

    if(found /= 0) cycle

    ngrid=MIN(nvector,ncache-igrid+1)
    do i=1,ngrid
      index_grid=active(levelmin)%igrid(igrid+i-1)

      do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        index_cell=iskip+index_grid

        dr_cell = 0.0d0
        do idim=1, ndim
          center(idim) = xg(index_grid,idim)+xc(ind,idim)*dx-skip_loc(idim)
          dr_cell = MAX(dr_cell, ABS(center(idim)-tcenter(idim)))
        enddo
        if(dr_cell.le.dx/2.0) then
          !$omp critical
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


if(found .eq. 0) then
  write(*,*) 'not matched !?'
  stop
endif

!FROM HERE

  !!----- Recursively walking through AMR cells and retrieve clouds
  ncloud = 0
if(myid.eq.26) write(*,*) ind_grid

if(myid.eq.26)then
  okay2=.false.
  nnn = 0
  rrx = 0.
  rry = 0.
  rrz = 0.

  ind_cell = 0
  ncache = active(levelmin)%ngrid
  do igrid=1, ncache, nvector
    ngrid = MIN(nvector,ncache-igrid+1)
    do i=1, ngrid
      index_grid=active(levelmin)%igrid(igrid+i-1)

      npart1 = numbp(index_grid)
      if(npart1 .eq. 0) cycle
      ipart = headp(index_grid)
      do jpart=1, npart1
        
        if(is_cloud(typep(ipart)) .and. idp(ipart) .eq. -sinkid) then
          nnn = nnn + 1
          rrx = rrx + xp(ipart,1)
          rry = rry + xp(ipart,2)
          rrz = rrz + xp(ipart,3)
          ind_cell = index_grid
        endif
        next_part = nextp(ipart)
        ipart=next_part
      enddo
    enddo
  enddo
  write(*,*) ind_grid , ' =?= ', ind_cell, sinkid
  write(*,*) nnn
  write(*,*) rrx/nnn, rry/nnn, rrz/nnn
  stop

  do ilevel=levelmin, nlevelmax
    ncache=active(ilevel)%ngrid

    do igrid=1, ncache, nvector
      ngrid=MIN(nvector,ncache-igrid+1)

      do i=1, ngrid
        index_grid=active(ilevel)%igrid(igrid+i-1)

        okay=.false.
        do ind=1, 8
          iskip=ncoarse+(ind-1)*ngridmax
          ind_cell = iskip+index_grid
    
          okay = (son(ind_cell)==0)

          if(okay) exit
        enddo

        if(okay) then
          npart1 = numbp(index_grid)
          if(npart1 .le. 0) cycle

          ipart = headp(index_grid)

          do jpart=1, npart1
            next_part = nextp(ipart)

            if(is_cloud(typep(ipart)) .and. idp(ipart) .eq. -sinkid) then
              nnn = nnn + 1
              rrx = rrx + xp(ipart,1)
              rry = rry + xp(ipart,2)
              rrz = rrz + xp(ipart,3)
            endif
              !okay2=.true.
            !if(idp(ipart) .eq. sinkid) okay2=.true.

            if(okay2) exit
          enddo
        endif

        if(okay2) exit
      enddo
      if(okay2) exit
    enddo
    if(okay2) exit
  enddo
    write(*,*) 'found here === ', index_grid, ilevel
    !write(*,*) xp(ipart,1), xp(ipart,2), xp(ipart,3)
    write(*,*) x, y, z

    write(*,*) nnn, sinkid
    write(*,*) rrx/nnn, rry/nnn, rrz/nnn
    !stop

    dx=0.5D0**ilevel

    do ind=1, twotondim
      write(*,*) ' edge = ', ind
      write(*,*) 'x0 = ', xg(index_grid,1)+xc(ind,1)*dx-skip_loc(1) - dx/2.0D0
      write(*,*) 'x1 = ', xg(index_grid,1)+xc(ind,1)*dx-skip_loc(1) + dx/2.0D0
      write(*,*) 'y0 = ', xg(index_grid,2)+xc(ind,2)*dx-skip_loc(2) - dx/2.0D0
      write(*,*) 'y1 = ', xg(index_grid,2)+xc(ind,2)*dx-skip_loc(2) + dx/2.0D0
      write(*,*) 'z0 = ', xg(index_grid,3)+xc(ind,3)*dx-skip_loc(3) - dx/2.0D0
      write(*,*) 'z1 = ', xg(index_grid,3)+xc(ind,3)*dx-skip_loc(3) + dx/2.0D0
      write(*,*) ' --- '

    enddo
  
    write(*,*) x-rmax, x+rmax
    write(*,*) y-rmax, y+rmax
    write(*,*) z-rmax, z+rmax

    stop
endif

  

  call subsub_walkgrid(xc, skip_loc, sinkid, son(ind_cell), levelmin, x, y, z, rmax, clouds, ncloud, nmax)

if(myid.eq.26)  write(*,*) 'cloud found well? : ', myid, ncloud

end subroutine
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
      write(*,*) i, j, nn
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
  
  dd2 = (rx-center(1))**2 + (ry-center(2))**2 + (rz-center(3))**2
  
end subroutine subsub_getdist2
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
  real(dp) :: subsub_dtmax

  subsub_maxv = maxval(varr)

  if(subsub_maxv .le. 0.0D0) return
  
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
subroutine subsub_getnbor(ind_cell, ind_nbcell, ind_nbgrid, ind_pos)
  use amr_commons
  implicit none
  integer::ind_cell
  integer, dimension(0:twondim)::ind_nbcell, ind_nbgrid, ind_pos
  !-----------------------------------------------------------------
  ! This subroutine determines the 2*ndim neighboring cells
  ! cells of the input cell (ind_cell).
  ! If for some reasons they don't exist, the routine returns
  ! 0.
  !
  ! copied from star_formation.f90
  !-----------------------------------------------------------------
  integer::i,j,iskip, ig, ih
  integer::pos,ind_grid_father
  integer,dimension(1:8,1:6)::ggg,hhh


  ggg(1:8,1)=(/1,0,1,0,1,0,1,0/); hhh(1:8,1)=(/2,1,4,3,6,5,8,7/)
  ggg(1:8,2)=(/0,2,0,2,0,2,0,2/); hhh(1:8,2)=(/2,1,4,3,6,5,8,7/)
  ggg(1:8,3)=(/3,3,0,0,3,3,0,0/); hhh(1:8,3)=(/3,4,1,2,7,8,5,6/)
  ggg(1:8,4)=(/0,0,4,4,0,0,4,4/); hhh(1:8,4)=(/3,4,1,2,7,8,5,6/)
  ggg(1:8,5)=(/5,5,5,5,0,0,0,0/); hhh(1:8,5)=(/5,6,7,8,1,2,3,4/)
  ggg(1:8,6)=(/0,0,0,0,6,6,6,6/); hhh(1:8,6)=(/5,6,7,8,1,2,3,4/)

  ! Get father cell position in the grid
  pos = (ind_cell-ncoarse-1)/ngridmax+1

  ! Get father grid
  ind_grid_father=ind_cell-ncoarse-(pos-1)*ngridmax
  
  ! Get neighboring father grids
  ind_nbgrid(0) = ind_grid_father
  do j=1, twondim
    ind_nbgrid(j)=son( nbor(ind_grid_father,j) )
  enddo

  ind_nbcell(0:twondim)=0
  ind_pos(0:twondim)=0
  ind_pos(0) = pos
  ind_nbcell(0) = ind_cell
  do j=1, twondim
    ig=ggg(pos,j)
    ih=hhh(pos,j)
    iskip=ncoarse+(ih-1)*ngridmax
    if(ind_nbgrid(ig)>0)then
      ind_nbcell(j)=iskip+ind_nbgrid(ig)
      ind_pos(j)=ih-1
    endif
  enddo

end subroutine subsub_getnbor
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_nbcellinput(ind_cell, ind_level, hvar)
  use amr_commons
  use subsub_commons
  use hydro_commons
  use hydro_parameters, ONLY: gamma
  implicit none

  integer :: ind_cell, ind_level
  real(dp), dimension(0:twondim, 1:subsub_nhydro) :: hvar

  !! Local variables
  integer :: i, j, ivar
  integer, dimension(0:twondim) :: ind_nbcell, ind_nbgrid, ind_pos
  real(dp), dimension(1:subsub_nhydro) :: hdum

  integer :: ix, iy, iz, ind
  real(dp),dimension(1:twotondim,1:3):: xc
  real(dp),dimension(1:3)::skip_loc
  real(dp) :: cx, cy, cz, dx, cx2, cy2, cz2

  real(dp) :: r, u, v, w, r0, u0, v0, w0, e0

  call subsub_getnbor(ind_cell, ind_nbcell, ind_nbgrid, ind_pos)

  !do ind=1,twotondim
  !  iz=(ind-1)/4
  !  iy=(ind-1-4*iz)/2
  !  ix=(ind-1-2*iy-4*iz)
  !  xc(ind,1)=(dble(ix)-0.5D0)
  !  xc(ind,2)=(dble(iy)-0.5D0)
  !  xc(ind,3)=(dble(iz)-0.5D0)
  !end do
  !skip_loc=(/0.0d0,0.0d0,0.0d0/)
  !skip_loc(1)=dble(icoarse_min)
  !skip_loc(2)=dble(jcoarse_min)
  !skip_loc(3)=dble(kcoarse_min)
  

  !!-----
  !! Host Cell
  !!-----
  do ivar=1, subsub_nhydro
    hvar(0,ivar) = uold(ind_cell, ivar)
  enddo

  return

  !!-----
  !! Shift to the object frame
  !!-----
  r0 = hvar(0,1)
  u0 = hvar(0,2)/r0
  v0 = hvar(0,3)/r0
  w0 = hvar(0,4)/r0
  e0 = hvar(0,5) - 0.5D0*r0*(u0**2 + v0**2 + w0**2)


  hvar(0,2:4) = 0.0D0
  hvar(0,5) = e0

  j = 0

  hdum(1:subsub_nhydro) = 0.0D0
  do i=1, twondim
    if(ind_nbcell(i).gt.0)then
      j = j+1

      do ivar=1, subsub_nhydro
        hvar(i,ivar) = uold(ind_nbcell(i),ivar)
      enddo

      
      !! Shift to the object frame
      r = hvar(i,1)
      u = hvar(i,2) / r
      v = hvar(i,3) / r
      w = hvar(i,4) / r

      hvar(i,2) = hvar(i,2) - u0*r
      hvar(i,3) = hvar(i,3) - v0*r
      hvar(i,4) = hvar(i,4) - w0*r

      
      hvar(i,5) = hvar(i,5) - (u0*u*r + v0*v*r + w0*w*r) + 0.5D0*r*(u0**2 + v0**2 + w0**2)

      hdum(1) = hdum(1) + hvar(i,1)
      hdum(2) = hdum(2) + hvar(i,2)/hvar(i,1)
      hdum(3) = hdum(3) + hvar(i,3)/hvar(i,1)
      hdum(4) = hdum(4) + hvar(i,4)/hvar(i,1)
      hdum(5) = hdum(5) + (gamma - 1.0D0)*(hvar(i,5) - 0.5D0*(hvar(i,2)**2 + hvar(i,3)**2 + hvar(i,4)**2)/hvar(i,1))
    endif
  enddo

  if(j.lt.3) then
    write(*,*) 'how it could be with j < 3?'
    stop 
  endif

  if(j.ne.6) then
    hdum(1) = hdum(1) / dble(j)
    hdum(2) = hdum(2) / dble(j)
    hdum(3) = hdum(3) / dble(j)
    hdum(4) = hdum(4) / dble(j)
    hdum(5) = hdum(5) / dble(j)
    !! Give the average density and velocity and pressure to empty arrays
    do i=1, twondim
      if(ind_nbcell(i).eq.0)then
        hvar(i,1) = hdum(1)
        hvar(i,2) = hdum(2) * hdum(1)
        hvar(i,3) = hdum(3) * hdum(1)
        hvar(i,4) = hdum(4) * hdum(1)
        hvar(i,5) = hdum(5) / (gamma - 1.0D0) + 0.5D0 * hdum(1) * (hdum(2)**2 + hdum(3)**2 + hdum(4)**2)
      endif
    enddo
  endif

  !! Fill empty neighbors to the mean values
  !do i=1, twondim
  !  if(ind_nbcell(i).eq.0)then
  !    do ivar=1, subsub_nhydro
  !      hvar(i,ivar) = hdum(ivar) / dble(j)
  !    enddo
  !  endif
  !enddo

end subroutine subsub_nbcellinput

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
