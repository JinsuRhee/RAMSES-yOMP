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
subroutine subsub_findcell(x, y, z, ind_cell, ind_grid, ind_level, flag)
  !!-----
  !! This routine finds a cell which (x, y, z) is located at
  !! (x,y,z) should be in a code unit
  use amr_commons
  implicit none
  real(dp), intent(in) :: x, y, z
  integer, intent(out) :: ind_cell, ind_grid, ind_level
  logical :: flag

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
    flag=.false.
    return
  else
    flag=.true.
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
