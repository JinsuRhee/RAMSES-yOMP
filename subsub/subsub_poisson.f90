!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_computefine(ind, delMgas, delMBH)
!!-----
!! This routine updates density field in a subsub object
!! 1) By prescribed boundary face flux
!! 2) Upwind flux

  use subsub_commons
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  integer :: ind
  real(dp) :: delMgas, delMBH

  !!----- Local Variables
  integer :: i
  integer :: subsub_nn, subsub_step
  
  real(dp) :: subsub_dt, subsub_t0, subsub_dtmax, subsub_dx, subsub_maxv
  real(dp), dimension(1:2,1:ndim) :: bc_vv, bc_flux
  real(dp), dimension(1:10) :: varr
  real(dp) :: m_in, m_out
  real(dp) :: rho_ave, tff, vff
  real(dp) :: delMtot

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::threepi2, fourpi
  !! Constant
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  threepi2=3.0d0*ACOS(-1.0d0)**2
  fourpi=4d0*acos(-1d0)

  subsub_nn = subsub_ngrid**ndim

  !! Check inflow/outflow type
  delMtot = delMgas + delMBH
  if(delMgas .gt. 0) then ! this part will be specified with different outflow assumptions
    m_in = delMgas
    m_out = 0.0D0
  else
    m_in = 0.0D0
    m_out = abs(delMgas)
  endif


  !! BC for inflow
  if(m_in .gt. 0) then
    if(subsub_inflowtype .eq. 1)then !! spherical inflow by free-fall
      rho_ave = subsub_obj(ind)%mass_tot + msink(subsub_obj(ind)%sink_ind)
      rho_ave = rho_ave / subsub_boxlen**ndim
      tff = sqrt(threepi2/8./fourpi/rho_ave/(6.67d-8*scale_d*scale_t**2))
      vff = subsub_boxlen/tff

      bc_vv(1,:) = vff
      bc_vv(2,:) = -vff
      bc_flux(1,:) = (m_in / 6.0D0 / subsub_boxlen**2) / (dtold(levelmin))
      bc_flux(2,:) = -(m_in / 6.0D0 / subsub_boxlen**2) / (dtold(levelmin))
    else
    endif
  endif
  if(m_out .gt. 0) then
    !!---- RHEE -----
    !! Outflow model will be added here
    !!---------------
  endif


  !!----- Find dt based on CFL condition
  subsub_t0 = 0.0D0

  !!----- Start Subcycle
  subsub_dt = dtold(levelmin)
  do
    !!----- Find dt by CFL
    varr(:) = 0.0D0
    varr(1) = maxval(subsub_obj(ind)%vg)
    varr(2) = maxval(bc_vv)
   
    call subsub_finddt(subsub_dt, dtold(levelmin), varr)

    !!----- Compute Phi
    call subsub_poissionnew(ind)

    !!----- Compute Force
    call subsub_forcenew(ind)

    !!----- Compute Velocity 
    call subsub_vnew(ind, subsub_dt)

    !!----- Compute Density (Poission + Continuity)
    call subsub_dnew(ind, subsub_dt, subsub_dx, bc_vv, bc_flux)


    subsub_t0 = subsub_t0 + subsub_dt 

    !!----- Total mass update
    subsub_obj(ind)%mass_tot = sum(subsub_obj(ind)%hydro(:,1)) * (subsub_boxlen/subsub_ngrid)**ndim
    if(subsub_t0 .ge. dtold(levelmin)) exit
  enddo

end subroutine subsub_computefine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_dnew(objind, subsub_dt, subsub_dx, bc_vv, bc_flux)
  use subsub_commons
  use subsub_parameters
  use amr_commons
  implicit none
  integer :: objind
  real(dp) :: subsub_dt, subsub_dx
  real(dp), dimension(1:2,1:ndim) :: bc_vv, bc_flux

  !! Local variables
  integer :: i, ii, ix, iy, iz, subsub_nn
  real(dp) :: bc_vx_up, bc_vx_down
  real(dp) :: bc_vy_up, bc_vy_down
  real(dp) :: bc_vz_up, bc_vz_down
  real(dp) :: bc_fx_up, bc_fx_down
  real(dp) :: bc_fy_up, bc_fy_down
  real(dp) :: bc_fz_up, bc_fz_down
  real(dp) :: subsub_delta

  real(dp), dimension(:), allocatable :: rho_old
  real(dp), dimension(:,:), allocatable :: vg_old
  real(dp), dimension(:,:), allocatable :: vg_up, vg_down
  real(dp), dimension(:,:), allocatable :: flux_up, flux_down

  subsub_nn = subsub_ngrid**ndim

  !!----- Allocate
  allocate(rho_old(1:subsub_nn))
  allocate(vg_old(1:subsub_nn,1:ndim))
  allocate(vg_up(1:subsub_nn,1:ndim))
  allocate(vg_down(1:subsub_nn,1:ndim))
  allocate(flux_up(1:subsub_nn,1:ndim))
  allocate(flux_down(1:subsub_nn,1:ndim))

  !!----- Save Old
  rho_old(:) = subsub_obj(objind)%hydro(:,1)
  do i=1, ndim
    vg_old(:,i) = subsub_obj(objind)%vg(:,i)
  enddo

  bc_vx_up = bc_vv(2,1)
  bc_vy_up = bc_vv(2,2)
  bc_vz_up = bc_vv(2,3)
  bc_vx_down = bc_vv(1,1)
  bc_vy_down = bc_vv(1,2)
  bc_vz_down = bc_vv(1,3)

  bc_fx_up = bc_flux(2,1)
  bc_fy_up = bc_flux(2,2)
  bc_fz_up = bc_flux(2,3)
  bc_fx_down = bc_flux(1,1)
  bc_fy_down = bc_flux(1,2)
  bc_fz_down = bc_flux(1,3)

  !!-----
  !! 1/2 Velocity
  !!-----
  
  !$omp parallel default(shared) private(i, ix, iy, iz, ii, subsub_delta)

  !$omp do
  do i=1, subsub_nn
    call subsub_get3ind(i, ix, iy, iz)


    !! X velocity
    if(ix.ne.subsub_ngrid) then
      call subsub_iget3ind(ii, ix+1, iy, iz)
      vg_up(i,1) = (vg_old(ii,1) + vg_old(i,1))/2.0D0
    else
      vg_up(i,1) = bc_vx_up
    endif

    if(ix.ne.1) then
      call subsub_iget3ind(ii, ix-1, iy, iz)
      vg_down(i,1) = (vg_old(ii,1) + vg_old(i,1))/2.0D0
    else
      vg_down(i,1) = bc_vx_down
    endif

    !! Y velocity
    if(iy.ne.subsub_ngrid) then
      call subsub_iget3ind(ii, ix, iy+1, iz)
      vg_up(i,2) = (vg_old(ii,2) + vg_old(i,2))/2.0D0
    else
      vg_up(i,2) = bc_vy_up
    endif

    if(iy.ne.1) then
      call subsub_iget3ind(ii, ix, iy-1, iz)
      vg_down(i,2) = (vg_old(ii,2) + vg_old(i,2))/2.0D0
    else
      vg_down(i,2) = bc_vy_down
    endif

    !! Z velocity
    if(iz.ne.subsub_ngrid) then
      call subsub_iget3ind(ii, ix, iy, iz+1)
      vg_up(i,3) = (vg_old(ii,3) + vg_old(i,3))/2.0D0
    else
      vg_up(i,3) = bc_vz_up
    endif

    if(iz.ne.1) then
      call subsub_iget3ind(ii, ix, iy, iz-1)
      vg_down(i,3) = (vg_old(ii,3) + vg_old(i,3))/2.0D0
    else
      vg_down(i,3) = bc_vz_down
    endif
  enddo
  !$omp end do

  !!-----
  !! 1/2 Flux (by UPWIND)
  !!-----
  !$omp do
  do i=1, subsub_nn
    call subsub_get3ind(i, ix, iy, iz)

    !! X UP
    if(vg_up(i,1) .GT. 0) then
      flux_up(i,1) = rho_old(i)*vg_up(i,1)
    else
      if(ix.ne.subsub_ngrid)then
        call subsub_iget3ind(ii, ix+1, iy, iz)
        flux_up(i,1) = rho_old(ii)*vg_up(i,1)
      else
        flux_up(i,1) = bc_fx_up
      endif
    endif

    !! X DOWN
    if(vg_down(i,1) .LT. 0) then
      flux_down(i,1) = rho_old(i)*vg_down(i,1)
    else
      if(ix.ne.1)then
        call subsub_iget3ind(ii, ix-1, iy, iz)
        flux_down(i,1) = rho_old(ii)*vg_down(i,1)
      else
        flux_down(i,1) = bc_fx_down
      endif
    endif

    !! Y UP
    if(vg_up(i,2) .GT. 0) then
      flux_up(i,2) = rho_old(i)*vg_up(i,2)
    else
      if(iy.ne.subsub_ngrid)then
        call subsub_iget3ind(ii, ix, iy+1, iz)
        flux_up(i,2) = rho_old(ii)*vg_up(i,2)
      else
        flux_up(i,2) = bc_fy_up
      endif
    endif

    !! Y DOWN
    if(vg_down(i,2) .LT. 0) then
      flux_down(i,2) = rho_old(i)*vg_down(i,2)
    else
      if(iy.ne.1)then
        call subsub_iget3ind(ii, ix, iy-1, iz)
        flux_down(i,2) = rho_old(ii)*vg_down(i,2)
      else
        flux_down(i,2) = bc_fy_down
      endif
    endif

    !! Z UP
    if(vg_up(i,3) .GT. 0) then
      flux_up(i,3) = rho_old(i)*vg_up(i,3)
    else
      if(iz.ne.subsub_ngrid)then
        call subsub_iget3ind(ii, ix, iy, iz+1)
        flux_up(i,3) = rho_old(ii)*vg_up(i,3)
      else
        flux_up(i,3) = bc_fz_up
      endif
    endif

    !! Z DOWN
    if(vg_down(i,3) .LT. 0) then
      flux_down(i,3) = rho_old(i)*vg_down(i,3)
    else
      if(iz.ne.1)then
        call subsub_iget3ind(ii, ix, iy, iz-1)
        flux_down(i,3) = rho_old(ii)*vg_down(i,3)
      else
        flux_down(i,3) = bc_fz_down
      endif
    endif
  enddo
  !$omp end do

  !!-----
  !! Update density
  !!-----
  !$omp do
  do i=1, subsub_nn
    call subsub_get3ind(i, ix, iy, iz)

    subsub_delta = flux_up(i,1) - flux_down(i,1) + &
      flux_up(i,2) - flux_down(i,2) + &
      flux_up(i,3) - flux_down(i,3)

    subsub_delta = subsub_delta * subsub_dt / subsub_dx

    subsub_obj(objind)%hydro(i,1) = max(rho_old(i) - subsub_delta, subsub_densityfloor)
  enddo
  !$omp end do

  !$omp end parallel


  !!----- Deallocate
  deallocate(rho_old)
  deallocate(vg_old)
  deallocate(vg_up)
  deallocate(vg_down)
  deallocate(flux_up)
  deallocate(flux_down)

end subroutine subsub_dnew
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_poissionnew(objind)
  use subsub_commons
  use subsub_parameters
  use pm_commons
  use amr_commons
  use cooling_module, ONLY:twopi

  implicit none
  integer :: objind

  !! Local variables
  integer :: i, ix, iy, iz, ii, niter, iter
  integer :: subsub_nn
  integer :: xu, xd, yu, yd, zu, zd
  real(dp) :: mtot, gconst, dd2, subsub_errmax, subsub_error0, fourpiG, dx2
  real(dp), dimension(:), allocatable :: subsub_phi
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  
  !!----- Constants
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  mtot = subsub_obj(objind)%mass_tot! + subsub_obj(objind)%sink_mass
  gconst = 6.67d-8*scale_d*scale_t**2
  fourpiG = gconst * 2d0*twopi
  subsub_nn = subsub_ngrid**ndim

  niter = subsub_poission_niter
  dx2 = (subsub_boxlen / dble(subsub_ngrid)) ** 2
  !!----- Initialize
  allocate(subsub_phi(1:subsub_nn))
  subsub_phi(:) = subsub_obj(objind)%phi(:)
  
  !! BC
  !$omp parallel do default(shared) private(ix, iy, iz, dd2)
  do i=1, subsub_nn
    call subsub_get3ind(i, ix, iy, iz)

    if(ix.eq.1 .or. ix .eq. subsub_ngrid .or. &
       iy.eq.1 .or. iy .eq. subsub_ngrid .or. &
       iz.eq.1 .or. iz .eq. subsub_ngrid) then

       call subsub_getdist2(ix, iy, iz, dd2)

       subsub_phi(i) = -gconst * mtot / sqrt(dd2 + subsub_softening**2)

    endif
  enddo
  !$omp end parallel do

  subsub_obj(objind)%phi(:) = subsub_phi(:)

  !! Relexation  
  subsub_errmax = 0.0D0
  do iter=1, niter
    
    subsub_error0 = 0.0D0

    !$omp parallel do default(shared) private(ii, ix, iy, iz) &
    !$omp & private(xu, xd, yu, yd, zu, zd) collapse(3) &
    !$omp & reduction(max:subsub_error0)
    do iz=2, subsub_ngrid-1
    do iy=2, subsub_ngrid-1
    do ix=2, subsub_ngrid-1
      call subsub_iget3ind(ii, ix, iy, iz)

      call subsub_iget3ind(xu, ix+1, iy, iz)
      call subsub_iget3ind(xd, ix-1, iy, iz)
      call subsub_iget3ind(yu, ix, iy+1, iz)
      call subsub_iget3ind(yd, ix, iy-1, iz)
      call subsub_iget3ind(zu, ix, iy, iz+1)
      call subsub_iget3ind(zd, ix, iy, iz-1)

      subsub_phi(ii) = (&
        subsub_obj(objind)%phi(xu) + subsub_obj(objind)%phi(xd) + &
        subsub_obj(objind)%phi(yu) + subsub_obj(objind)%phi(yd) + &
        subsub_obj(objind)%phi(zu) + subsub_obj(objind)%phi(zd) - &
        dx2 * subsub_obj(objind)%hydro(ii,1) * fourpiG ) / 6.0D0
      
      subsub_error0 = max(abs(subsub_phi(ii) - subsub_obj(objind)%phi(ii)), subsub_smallr)
      
    enddo
    enddo
    enddo
    !$omp end parallel do
    
    subsub_errmax = subsub_error0 / maxval( abs(subsub_obj(objind)%phi) ) 

    subsub_obj(objind)%phi(:) = subsub_phi(:)

    if(subsub_errmax .lt. subsub_poissiontolerance) exit
  enddo

  deallocate(subsub_phi)
end subroutine subsub_poissionnew
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_forcenew(objind)
  use subsub_commons
  use pm_commons
  use amr_commons

  implicit none
  integer :: objind

  !! Local variables
  integer :: i, ii, ix, iy, iz
  integer :: xu, xd, yu, yd, zu, zd
  real(dp) :: dx

  !!----- Constants
  dx = subsub_boxlen / dble(subsub_ngrid)

  !!----- initialize
  subsub_obj(objind)%fg(:,:) = 0.0D0

  !!----- Force from Phi
  !$omp parallel do default(shared) private(ii, ix, iy, iz) &
  !$omp & private(xu, xd, yu, yd, zu, zd) collapse(3)
  do iz=2, subsub_ngrid-1
  do iy=2, subsub_ngrid-1
  do ix=2, subsub_ngrid-1
      call subsub_iget3ind(ii, ix, iy, iz)

      call subsub_iget3ind(xu, ix+1, iy, iz)
      call subsub_iget3ind(xd, ix-1, iy, iz)
      call subsub_iget3ind(yu, ix, iy+1, iz)
      call subsub_iget3ind(yd, ix, iy-1, iz)
      call subsub_iget3ind(zu, ix, iy, iz+1)
      call subsub_iget3ind(zd, ix, iy, iz-1)

      subsub_obj(objind)%fg(ii,1) = - (subsub_obj(objind)%phi(xu) - subsub_obj(objind)%phi(xd)) / (2.0D0*dx)
      subsub_obj(objind)%fg(ii,2) = - (subsub_obj(objind)%phi(yu) - subsub_obj(objind)%phi(yd)) / (2.0D0*dx)
      subsub_obj(objind)%fg(ii,3) = - (subsub_obj(objind)%phi(zu) - subsub_obj(objind)%phi(zd)) / (2.0D0*dx)
  enddo
  enddo
  enddo
  !$omp end parallel do
end subroutine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_vnew(objind, dt)
  use subsub_commons
  use pm_commons
  use amr_commons

  implicit none
  integer :: objind
  real(dp) :: dt

  !! Local variables
  integer :: i, idim, subsub_nn

  subsub_nn = subsub_ngrid**ndim

  !$omp parallel do default(shared) private(idim)
  do i=1, subsub_nn
    do idim=1, ndim
      subsub_obj(objind)%vg(i, idim) = subsub_obj(objind)%vg(i, idim) + subsub_obj(objind)%fg(i, idim) * dt
    enddo
  enddo
  !$omp end parallel do
end subroutine subsub_vnew