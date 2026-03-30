subroutine subsub_dnew(objind, subsub_dt, bc_vv, bc_flux)
  use subsub_commons
  use subsub_parameters
  use amr_commons
  implicit none
  integer :: objind
  real(dp) :: subsub_dt
  real(dp), dimension(1:2,1:ndim) :: bc_vv, bc_flux

  !! Local variables
  integer :: i, ii, ix, iy, iz
  real(dp) :: bc_vx_up, bc_vx_down
  real(dp) :: bc_vy_up, bc_vy_down
  real(dp) :: bc_vz_up, bc_vz_down
  real(dp) :: bc_fx_up, bc_fx_down
  real(dp) :: bc_fy_up, bc_fy_down
  real(dp) :: bc_fz_up, bc_fz_down
  real(dp) :: subsub_delta

  !!----- Save Old
  subsub_rho_old(:) = subsub_obj(objind)%hydro(:,1)
  do i=1, ndim
    subsub_vg_old(:,i) = subsub_obj(objind)%hydro(:,i+1) / subsub_rho_old(:)
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
      !call subsub_iget3ind(ii, ix+1, iy, iz)
      ii = (iz-1)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix+1
      
      subsub_vg_up(i,1) = (subsub_vg_old(ii,1) + subsub_vg_old(i,1))/2.0D0
    else
      subsub_vg_up(i,1) = bc_vx_up
    endif

    if(ix.ne.1) then
      !call subsub_iget3ind(ii, ix-1, iy, iz)
      ii = (iz-1)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix-1

      subsub_vg_down(i,1) = (subsub_vg_old(ii,1) + subsub_vg_old(i,1))/2.0D0
    else
      subsub_vg_down(i,1) = bc_vx_down
    endif

    !! Y velocity
    if(iy.ne.subsub_ngrid) then
      !call subsub_iget3ind(ii, ix, iy+1, iz)
      ii = (iz-1)*subsub_ngrid**2 + (iy)*subsub_ngrid + ix

      subsub_vg_up(i,2) = (subsub_vg_old(ii,2) + subsub_vg_old(i,2))/2.0D0
    else
      subsub_vg_up(i,2) = bc_vy_up
    endif

    if(iy.ne.1) then
      !call subsub_iget3ind(ii, ix, iy-1, iz)
      ii = (iz-1)*subsub_ngrid**2 + (iy-2)*subsub_ngrid + ix

      subsub_vg_down(i,2) = (subsub_vg_old(ii,2) + subsub_vg_old(i,2))/2.0D0
    else
      subsub_vg_down(i,2) = bc_vy_down
    endif

    !! Z velocity
    if(iz.ne.subsub_ngrid) then
      !call subsub_iget3ind(ii, ix, iy, iz+1)
      ii = (iz)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix

      subsub_vg_up(i,3) = (subsub_vg_old(ii,3) + subsub_vg_old(i,3))/2.0D0
    else
      subsub_vg_up(i,3) = bc_vz_up
    endif

    if(iz.ne.1) then
      !call subsub_iget3ind(ii, ix, iy, iz-1)
      ii = (iz-2)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix

      subsub_vg_down(i,3) = (subsub_vg_old(ii,3) + subsub_vg_old(i,3))/2.0D0
    else
      subsub_vg_down(i,3) = bc_vz_down
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
    if(subsub_vg_up(i,1) .GT. 0) then
      subsub_flux_up(i,1) = subsub_rho_old(i)*subsub_vg_up(i,1)
    else
      if(ix.ne.subsub_ngrid)then
        !call subsub_iget3ind(ii, ix+1, iy, iz)
        ii = (iz-1)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix+1

        subsub_flux_up(i,1) = subsub_rho_old(ii)*subsub_vg_up(i,1)
      else
        subsub_flux_up(i,1) = bc_fx_up
      endif
    endif

    !! X DOWN
    if(subsub_vg_down(i,1) .LT. 0) then
      subsub_flux_down(i,1) = subsub_rho_old(i)*subsub_vg_down(i,1)
    else
      if(ix.ne.1)then
        !call subsub_iget3ind(ii, ix-1, iy, iz)
        ii = (iz-1)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix-1

        subsub_flux_down(i,1) = subsub_rho_old(ii)*subsub_vg_down(i,1)
      else
        subsub_flux_down(i,1) = bc_fx_down
      endif
    endif

    !! Y UP
    if(subsub_vg_up(i,2) .GT. 0) then
      subsub_flux_up(i,2) = subsub_rho_old(i)*subsub_vg_up(i,2)
    else
      if(iy.ne.subsub_ngrid)then
        !call subsub_iget3ind(ii, ix, iy+1, iz)
        ii = (iz-1)*subsub_ngrid**2 + (iy)*subsub_ngrid + ix

        subsub_flux_up(i,2) = subsub_rho_old(ii)*subsub_vg_up(i,2)
      else
        subsub_flux_up(i,2) = bc_fy_up
      endif
    endif

    !! Y DOWN
    if(subsub_vg_down(i,2) .LT. 0) then
      subsub_flux_down(i,2) = subsub_rho_old(i)*subsub_vg_down(i,2)
    else
      if(iy.ne.1)then
        !call subsub_iget3ind(ii, ix, iy-1, iz)
        ii = (iz-1)*subsub_ngrid**2 + (iy-2)*subsub_ngrid + ix

        subsub_flux_down(i,2) = subsub_rho_old(ii)*subsub_vg_down(i,2)
      else
        subsub_flux_down(i,2) = bc_fy_down
      endif
    endif

    !! Z UP
    if(subsub_vg_up(i,3) .GT. 0) then
      subsub_flux_up(i,3) = subsub_rho_old(i)*subsub_vg_up(i,3)
    else
      if(iz.ne.subsub_ngrid)then
        !call subsub_iget3ind(ii, ix, iy, iz+1)
        ii = (iz)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix

        subsub_flux_up(i,3) = subsub_rho_old(ii)*subsub_vg_up(i,3)
      else
        subsub_flux_up(i,3) = bc_fz_up
      endif
    endif

    !! Z DOWN
    if(subsub_vg_down(i,3) .LT. 0) then
      subsub_flux_down(i,3) = subsub_rho_old(i)*subsub_vg_down(i,3)
    else
      if(iz.ne.1)then
        !call subsub_iget3ind(ii, ix, iy, iz-1)
        ii = (iz-2)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix

        subsub_flux_down(i,3) = subsub_rho_old(ii)*subsub_vg_down(i,3)
      else
        subsub_flux_down(i,3) = bc_fz_down
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

    subsub_delta = subsub_flux_up(i,1) - subsub_flux_down(i,1) + &
      subsub_flux_up(i,2) - subsub_flux_down(i,2) + &
      subsub_flux_up(i,3) - subsub_flux_down(i,3)

    subsub_delta = subsub_delta * subsub_dt / subsub_dx

    subsub_obj(objind)%hydro(i,1) = max(subsub_rho_old(i) - subsub_delta, subsub_densityfloor)
  enddo
  !$omp end do

  !$omp end parallel


end subroutine subsub_dnew
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_poissonnew(objind)
  use subsub_commons
  use subsub_parameters
  use pm_commons
  use amr_commons
  use cooling_module, ONLY:twopi
  use mpi_mod
  implicit none
  integer :: objind

  !! Local variables
  
  
  
  !!----- Initial Guess
  subsub_phi(:) = subsub_obj(objind)%phi(:)
subsub_tcheck(20) = MPI_WTIME()  
  !!----- BC
  call subsub_poissonnew_BC(objind, subsub_phi, subsub_nn)
subsub_tcheck(21) = MPI_WTIME()
  !!----- Copy to the Old Poisson
  subsub_obj(objind)%phi(:) = subsub_phi(:)

  !!----- Poisson Solver
  if(subsub_poisson_type .eq. 1) then
    !! By 6-point Jaccobian relaxation
    call subsub_poissonnew_6point(objind)
  else if(subsub_poisson_type .eq. 2) then
    !! By Conjugate-Gradient
    call subsub_poissonnew_cg(objind)
  else
    call subsub_log('No Possion Solver found', 'subsub_poissonnew')
    stop
  endif
subsub_tcheck(22) = MPI_WTIME()
  !!----- Update
  subsub_obj(objind)%phi(:) = subsub_phi(:)
  
end subroutine subsub_poissonnew
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_poissonnew_BC(objind, inphi, nn)
  use subsub_commons
  use amr_commons
  use cooling_module, ONLY:twopi

  integer :: objind, nn
  real(dp), dimension(1:nn) :: inphi

  real(dp) ::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  integer :: i, ix, iy, iz
  real(dp) :: mtot, gconst, dd2, fourpiG, dx2

  !!----- Constants
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  mtot = subsub_obj(objind)%mass_tot! + subsub_obj(objind)%sink_mass
  gconst = 6.67d-8*scale_d*scale_t**2
  fourpiG = gconst * 2d0*twopi

  !$omp parallel do default(shared) private(ix, iy, iz, dd2)
  do i=1, nn
    call subsub_get3ind(i, ix, iy, iz)

    if(ix.eq.1 .or. ix .eq. subsub_ngrid .or. &
       iy.eq.1 .or. iy .eq. subsub_ngrid .or. &
       iz.eq.1 .or. iz .eq. subsub_ngrid) then

       call subsub_getdist2(ix, iy, iz, dd2)

       inphi(i) = -gconst * mtot / sqrt(dd2 + subsub_poisson_softening**2)

    endif
  enddo
  !$omp end parallel do
end subroutine subsub_poissonnew_BC
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_poissonnew_zeroBC(A)
  use subsub_commons

  real(dp), dimension(1:subsub_nn) :: A

  !! Local variables
  integer :: i, ix, iy, iz
  !$omp parallel do default(shared) private(ix, iy, iz)
  do i=1, subsub_nn
    call subsub_get3ind(i, ix, iy, iz)

    if(ix.eq.1 .or. ix .eq. subsub_ngrid .or. &
       iy.eq.1 .or. iy .eq. subsub_ngrid .or. &
       iz.eq.1 .or. iz .eq. subsub_ngrid) then

       A(i) = 0.0D0

    endif
  enddo
  !$omp end parallel do
end subroutine subsub_poissonnew_zeroBC
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_poissonnew_6point(objind)
  use subsub_commons
  use subsub_parameters
  use amr_commons
  use cooling_module, ONLY:twopi
  implicit none

  integer :: objind

  !! Local variables
  integer :: iter, ix, iy, iz, ii
  integer :: niter
  real(dp) :: subsub_errmax, subsub_error0, subsub_error1, dx2, fourpiG, gconst
  integer :: xu, xd, yu, yd, zu, zd
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  
  !!----- Constants
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  gconst = 6.67d-8*scale_d*scale_t**2
  fourpiG = gconst * 2d0*twopi

  !! Relexation
  niter = subsub_poisson_niter
  dx2 = (subsub_boxlen / dble(subsub_ngrid)) ** 2
  subsub_errmax = 0.0D0


  do iter=1, niter
    
    subsub_error0 = 0.0D0
    subsub_error1 = 0.0D0

    !$omp parallel do default(shared) private(ii, ix, iy, iz) &
    !$omp & private(xu, xd, yu, yd, zu, zd) collapse(3) &
    !$omp & reduction(+:subsub_error0, subsub_error1)
    do iz=2, subsub_ngrid-1
    do iy=2, subsub_ngrid-1
    do ix=2, subsub_ngrid-1
      !call subsub_iget3ind(ii, ix, iy, iz)

      !call subsub_iget3ind(xu, ix+1, iy, iz)
      !call subsub_iget3ind(xd, ix-1, iy, iz)
      !call subsub_iget3ind(yu, ix, iy+1, iz)
      !call subsub_iget3ind(yd, ix, iy-1, iz)
      !call subsub_iget3ind(zu, ix, iy, iz+1)
      !call subsub_iget3ind(zd, ix, iy, iz-1)

      ii = (iz-1)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix

      xu = (iz-1)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix+1
      xd = (iz-1)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix-1
      yu = (iz-1)*subsub_ngrid**2 + (iy)*subsub_ngrid + ix
      yd = (iz-1)*subsub_ngrid**2 + (iy-2)*subsub_ngrid + ix
      zu = (iz)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix
      zd = (iz-2)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix


      subsub_phi(ii) = (&
        subsub_obj(objind)%phi(xu) + subsub_obj(objind)%phi(xd) + &
        subsub_obj(objind)%phi(yu) + subsub_obj(objind)%phi(yd) + &
        subsub_obj(objind)%phi(zu) + subsub_obj(objind)%phi(zd) - &
        dx2 * subsub_obj(objind)%hydro(ii,1) * fourpiG ) / 6.0D0


      !! Error by L^2
      subsub_error0 = subsub_error0 + (subsub_phi(ii) - subsub_obj(objind)%phi(ii))**2
      subsub_error1 = subsub_error1 + subsub_obj(objind)%phi(ii)**2
    enddo
    enddo
    enddo
    !$omp end parallel do

    subsub_obj(objind)%phi(:) = subsub_phi(:)

    subsub_errmax = subsub_error0 / max(subsub_error1, subsub_smallr)
    if(subsub_errmax .lt. subsub_poisson_tolerance) exit
  enddo

end subroutine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_poissonnew_cg(objind)
  use subsub_commons
  use subsub_parameters
  use amr_commons
  use cooling_module, ONLY:twopi
  implicit none

  integer :: objind

  !! Local variables
  integer :: i, ix, iy, iz, ii, iter
  integer :: niter
  real(dp) :: subsub_errmax, subsub_error0, dx2, fourpiG, gconst
  integer :: xu, xd, yu, yd, zu, zd
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

  real(dp) :: rr_old, rr_new, pLp, alpha, beta, relres, rhs2, dx2inv, mtot, dd2
  logical :: iterflag
  real(dp), dimension(1:3) :: subsub_center
  real(dp) :: rx, ry, rz

  !!----- Constants
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  gconst = 6.67d-8*scale_d*scale_t**2
  fourpiG = gconst * 2d0*twopi

  mtot = subsub_obj(objind)%mass_tot

  niter = subsub_poisson_niter
  dx2inv = 1.0D0 / (subsub_dx**2)

  subsub_center = (/subsub_boxlen/2.0D0, subsub_boxlen/2.0D0, subsub_boxlen/2.0D0/)

  !!----- Initialize
  subsub_cgrhs(:) = 0.0D0
  subsub_cgLphi(:) = 0.0D0
  subsub_cgRes(:) = 0.0D0
  subsub_cgp(:) = 0.0D0
  subsub_cgLp(:) = 0.0D0

  subsub_debugn = 0
  !!-----
  !! Ax = b
  !! r_n = Ax_n - b
  !! alpha = [ (r_n)^T r_n ] / [ (x_n)^T Ax_n]
  !! x_n+1 = x_n + alpha * r_n
  !! r_n+1 = r_n - alpha * Ap_n
  !! beta = [ (r_n+1)^T r_n+1 ] / [ (r_n)^T r_n]
  !! p_n+1 = r_n+1 + beta * p_n
  !!
  !!
  !! x_0 = subsub_obj(objind)%phi (initial guess)
  !! x_n = subsub_phi
  !! r_n = subsub_cgRes
  !! p_n = subsub_cgP


  !!----- Build RHS & com rhs2 & Build Lphi & cgRes , cgP & rr_old
  rhs2 = 0.0D0
  rr_old = 0.0D0 
  !$omp parallel default(shared) private(i, ix, iy, iz) &
  !$omp & private(xu, xd, yu, yd, zu, zd)

  !$omp do collapse(3) reduction(+:rhs2, rr_old)
  do ix=2, subsub_ngrid-1
  do iy=2, subsub_ngrid-1
  do iz=2, subsub_ngrid-1
    !call subsub_iget3ind(i, ix, iy, iz)
    i = (iz-1)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix

    subsub_cgrhs(i) = -fourpiG * subsub_obj(objind)%hydro(i,1)
    rhs2 = rhs2 + subsub_cgrhs(i)*subsub_cgrhs(i)

    xu = (iz-1)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix+1
    xd = (iz-1)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix-1
    yu = (iz-1)*subsub_ngrid**2 + (iy)*subsub_ngrid + ix
    yd = (iz-1)*subsub_ngrid**2 + (iy-2)*subsub_ngrid + ix
    zu = (iz)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix
    zd = (iz-2)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix

    subsub_cgLphi(i) = (&
      6.0D0*subsub_phi(i) - subsub_phi(xu) - subsub_phi(xd) &
      - subsub_phi(yu) - subsub_phi(yd) - subsub_phi(zu) - subsub_phi(zd) )  * dx2inv

    subsub_cgRes(i) = subsub_cgrhs(i) - subsub_cgLphi(i)
    subsub_cgp(i) = subsub_cgRes(i)
    rr_old = rr_old + subsub_cgRes(i)*subsub_cgRes(i)
  enddo
  enddo
  enddo
  !$omp end do

  !$omp end parallel

  if(sqrt(rr_old / max(rhs2, subsub_smallr)) .le. subsub_poisson_tolerance) then
    return
  endif

  !!----- Iteration
  iterflag = .false.

  !$omp parallel default(shared) &
  !$omp & private(ix, iy, iz, ii, xu, xd, yu, yd, zu, zd) &
  !$omp & private(dd2, rx, ry, rz)
  do iter=1, niter

    !$omp single
    pLp = 0.0D0
    rr_new = 0.0D0
    iterflag = .false.
    !$omp end single

    !!----- compute Lp & pLp
    !call subsub_poissonnew_cg_Laplacian(subsub_cgp, subsub_cgLp, subsub_nn)

    !$omp do collapse(3) reduction(+:pLp)
    do iz=1, subsub_ngrid
    do iy=1, subsub_ngrid
    do ix=1, subsub_ngrid
      ii = (iz-1)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix
      
      !! zero BC
      if(ix.eq.1 .or. ix.eq.subsub_ngrid)then
        subsub_cgLp(ii) = 0
        cycle
      endif
      if(iy.eq.1 .or. iy.eq.subsub_ngrid)then
        subsub_cgLp(ii) = 0
        cycle
      endif
      if(iz.eq.1 .or. iz.eq.subsub_ngrid)then
        subsub_cgLp(ii) = 0
        cycle
      endif

      xu = (iz-1)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix+1
      xd = (iz-1)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix-1
      yu = (iz-1)*subsub_ngrid**2 + (iy)*subsub_ngrid + ix
      yd = (iz-1)*subsub_ngrid**2 + (iy-2)*subsub_ngrid + ix
      zu = (iz)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix
      zd = (iz-2)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix

      subsub_cgLp(ii) = (&
        6.0D0*subsub_cgp(ii) - subsub_cgp(xu) - subsub_cgp(xd) &
        - subsub_cgp(yu) - subsub_cgp(yd) - subsub_cgp(zu) - subsub_cgp(zd) )  * dx2inv

      pLp = pLp + subsub_cgp(ii) * subsub_cgLp(ii)
    enddo
    enddo
    enddo
    !$omp end do

    !$omp single
    if(abs(pLp) .lt. subsub_smallr) then
      iterflag=.true.
      alpha = 0.0D0
    else
      alpha = rr_old / pLp
    endif
    !$omp end single

    !!----- update phi & set BC & update Res & compute new rr
    
    !$omp do collapse(3) reduction(+:rr_new)
    do iz=1, subsub_ngrid
    do iy=1, subsub_ngrid
    do ix=1, subsub_ngrid
      ii = (iz-1)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix
 
      subsub_cgRes(ii) = subsub_cgRes(ii) - alpha*subsub_cgLp(ii)

      !! BC
      if(ix.eq.1 .or. ix.eq.subsub_ngrid)then
        rx = (dble(ix) - 0.5D0) * subsub_dx
        ry = (dble(iy) - 0.5D0) * subsub_dx
        rz = (dble(iz) - 0.5D0) * subsub_dx
        dd2 = (rx-subsub_center(1))**2 + (ry-subsub_center(2))**2 + (rz-subsub_center(3))**2

        subsub_phi(ii) = -gconst * mtot / sqrt(dd2 + subsub_poisson_softening**2)
        cycle
      endif
      if(iy.eq.1 .or. iy.eq.subsub_ngrid)then
        rx = (dble(ix) - 0.5D0) * subsub_dx
        ry = (dble(iy) - 0.5D0) * subsub_dx
        rz = (dble(iz) - 0.5D0) * subsub_dx
        dd2 = (rx-subsub_center(1))**2 + (ry-subsub_center(2))**2 + (rz-subsub_center(3))**2

        subsub_phi(ii) = -gconst * mtot / sqrt(dd2 + subsub_poisson_softening**2)
        cycle
      endif
      if(iz.eq.1 .or. iz.eq.subsub_ngrid)then
        rx = (dble(ix) - 0.5D0) * subsub_dx
        ry = (dble(iy) - 0.5D0) * subsub_dx
        rz = (dble(iz) - 0.5D0) * subsub_dx
        dd2 = (rx-subsub_center(1))**2 + (ry-subsub_center(2))**2 + (rz-subsub_center(3))**2
        
        subsub_phi(ii) = -gconst * mtot / sqrt(dd2 + subsub_poisson_softening**2)
        cycle
      endif

      rr_new = rr_new + subsub_cgRes(ii)*subsub_cgRes(ii)
      subsub_phi(ii) = subsub_phi(ii) + alpha * subsub_cgp(ii)
    enddo
    enddo
    enddo
    !$omp end do

    !$omp single
    relres = sqrt(rr_new / max(rhs2, subsub_smallr))
    if(relres .lt. subsub_poisson_tolerance) iterflag=.true.
    beta = rr_new / max(rr_old, subsub_smallr)
    !$omp end single

    !!----- Update cgP & set BC
    !$omp do collapse(3)
    do iz=1, subsub_ngrid
    do iy=1, subsub_ngrid
    do ix=1, subsub_ngrid
      ii = (iz-1)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix
      
      !! zero BC
      if(ix.eq.1 .or. ix.eq.subsub_ngrid)then
        subsub_cgp(ii) = 0
        cycle
      endif
      if(iy.eq.1 .or. iy.eq.subsub_ngrid)then
        subsub_cgp(ii) = 0
        cycle
      endif
      if(iz.eq.1 .or. iz.eq.subsub_ngrid)then
        subsub_cgp(ii) = 0
        cycle
      endif

      subsub_cgp(ii) = subsub_cgRes(ii) + beta * subsub_cgp(ii)
    enddo
	enddo
	enddo
    !$omp end do

    !$omp single
    rr_old = rr_new
    !$omp end single

    if(iterflag) exit
  enddo
  !$omp end parallel

  subsub_debugn = iter
  !!----- Deallocations
  !if(allocated(subsub_cgrhs)) deallocate(subsub_cgrhs)
  !if(allocated(subsub_cgLphi)) deallocate(subsub_cgLphi)
  !if(allocated(subsub_cgRes)) deallocate(subsub_cgRes)
  !if(allocated(subsub_cgp)) deallocate(subsub_cgp)
  !if(allocated(subsub_cgLp)) deallocate(subsub_cgLp)

end subroutine subsub_poissonnew_cg
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_poissonnew_cg_Laplacian(phi, Lphi, nn)
  use subsub_commons

  implicit none
  integer :: nn
  real(dp), dimension(1:nn) :: phi, Lphi

  !! Local variables
  integer :: i, ix, iy, iz, ii
  integer :: xu, xd, yu, yd, zu, zd
  real(dp) :: dx2inv

  !! Constants
  dx2inv = 1.0D0 / (subsub_dx**2)

  !! Initialize
  Lphi(:) = 0.0D0

  !$omp parallel do default(shared) collapse(3) &
  !$omp & private(ix, iy, iz, ii, xu, xd, yu, yd, zu, zd)
  do iz=2, subsub_ngrid-1
  do iy=2, subsub_ngrid-1
  do ix=2, subsub_ngrid-1
      !call subsub_iget3ind(ii, ix, iy, iz)
      !call subsub_iget3ind(xu, ix+1, iy, iz)
      !call subsub_iget3ind(xd, ix-1, iy, iz)
      !call subsub_iget3ind(yu, ix, iy+1, iz)
      !call subsub_iget3ind(yd, ix, iy-1, iz)
      !call subsub_iget3ind(zu, ix, iy, iz+1)
      !call subsub_iget3ind(zd, ix, iy, iz-1)

      ii = (iz-1)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix

      xu = (iz-1)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix+1
      xd = (iz-1)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix-1
      yu = (iz-1)*subsub_ngrid**2 + (iy)*subsub_ngrid + ix
      yd = (iz-1)*subsub_ngrid**2 + (iy-2)*subsub_ngrid + ix
      zu = (iz)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix
      zd = (iz-2)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix

      Lphi(ii) = (&
        6.0D0*phi(ii) - phi(xu) - phi(xd) &
        - phi(yu) - phi(yd) - phi(zu) - phi(zd) )  * dx2inv
  enddo
  enddo
  enddo
  !$omp end parallel do
end subroutine subsub_poissonnew_cg_Laplacian
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_poissonnew_cg_interiordot(A, B, idot)
  use subsub_commons
  implicit none

  real(dp), dimension(1:subsub_nn) :: A, B
  real(dp) :: idot

  !! Local variables
  integer :: i, ix, iy, iz, ii

  idot = 0.0D0

  !$omp parallel do default(shared) collapse(3) &
  !$omp & private(ix, iy, iz, ii) reduction(+:idot)
  do iz=2, subsub_ngrid-1
  do iy=2, subsub_ngrid-1
  do ix=2, subsub_ngrid-1
      !call subsub_iget3ind(ii, ix, iy, iz)

      ii = (iz-1)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix

      idot = idot + A(ii) * B(ii)
  enddo
  enddo
  enddo
  !$omp end parallel do

end subroutine
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
  subsub_fg(:,:) = 0.0D0

  !!----- Force from Phi
  !$omp parallel do default(shared) private(ii, ix, iy, iz) &
  !$omp & private(xu, xd, yu, yd, zu, zd) collapse(3)
  do iz=1, subsub_ngrid
  do iy=1, subsub_ngrid
  do ix=1, subsub_ngrid

      ii = (iz-1)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix
      if(ix.eq.1 .or. ix.eq.subsub_ngrid)then
        xu = (iz-1)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix+1
        xd = (iz-1)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix-1
        subsub_fg(ii,1) = - (subsub_obj(objind)%phi(xu) - subsub_obj(ii)) / dx

here when two more boundaries overlapped
      !call subsub_iget3ind(ii, ix, iy, iz)

      !call subsub_iget3ind(xu, ix+1, iy, iz)
      !call subsub_iget3ind(xd, ix-1, iy, iz)
      !call subsub_iget3ind(yu, ix, iy+1, iz)
      !call subsub_iget3ind(yd, ix, iy-1, iz)
      !call subsub_iget3ind(zu, ix, iy, iz+1)
      !call subsub_iget3ind(zd, ix, iy, iz-1)

      

      xu = (iz-1)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix+1
      xd = (iz-1)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix-1
      yu = (iz-1)*subsub_ngrid**2 + (iy)*subsub_ngrid + ix
      yd = (iz-1)*subsub_ngrid**2 + (iy-2)*subsub_ngrid + ix
      zu = (iz)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix
      zd = (iz-2)*subsub_ngrid**2 + (iy-1)*subsub_ngrid + ix

      subsub_fg(ii,1) = - (subsub_obj(objind)%phi(xu) - subsub_obj(objind)%phi(xd)) / (2.0D0*dx)
      subsub_fg(ii,2) = - (subsub_obj(objind)%phi(yu) - subsub_obj(objind)%phi(yd)) / (2.0D0*dx)
      subsub_fg(ii,3) = - (subsub_obj(objind)%phi(zu) - subsub_obj(objind)%phi(zd)) / (2.0D0*dx)
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
  integer :: i, idim

  !$omp parallel do default(shared) private(idim)
  do i=1, subsub_nn
    do idim=1, ndim
      subsub_obj(objind)%hydro(i, idim+1) = subsub_obj(objind)%hydro(i, idim+1) + &
        subsub_fg(i, idim) * dt * subsub_obj(objind)%hydro(i, 1)
    enddo
  enddo
  !$omp end parallel do
end subroutine subsub_vnew