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
  use mpi_mod
  implicit none
  integer :: ind
  real(dp) :: delMgas, delMBH

  !!----- Local Variables
  integer :: i
  integer :: subsub_step, subsub_nstep
  
  real(dp) :: subsub_dt, subsub_t0, subsub_dtmax, subsub_maxv
  real(dp), dimension(1:2,1:ndim) :: bc_vv, bc_flux
  real(dp), dimension(1:10) :: varr
  real(dp) :: m_in, m_out
  real(dp) :: rho_ave, tff, vff
  real(dp) :: delMtot

  real(dp) :: scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp) :: threepi2, fourpi

  real(dp) :: tcheck(1:10)

  !! Constant
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  threepi2=3.0d0*ACOS(-1.0d0)**2
  fourpi=4d0*acos(-1d0)

  !! Check inflow/outflow type
  delMtot = delMgas + delMBH
  if(delMgas .gt. 0) then ! this part will be specified with different outflow assumptions
    m_in = delMgas
    m_out = 0.0D0
  else
    m_in = 0.0D0
    m_out = abs(delMgas)
  endif


  !!-----
  !! Set BC
  !!-----
  bc_vv(:,:) = 0.0D0
  bc_flux(:,:) = 0.0D0

  !!-----for inflow
  if(m_in .gt. 0) then
    if(subsub_inflowtype .eq. 1)then !! spherical inflow by free-fall
      rho_ave = subsub_obj(ind)%mass_tot! + msink(subsub_obj(ind)%sink_ind)
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


  !!----- 
  subsub_t0 = 0.0D0
  subsub_nstep = 0

  !!----- Start Subcycle
  subsub_dt = dtold(levelmin)


  !!----- By Conjugrate-Gradient & KDK
  call subsub_cgkdk(ind, bc_vv, bc_flux)


end subroutine subsub_computefine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_cgkdk(objind, bc_vv, bc_flux)
  !!----- CG
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
  !!------

  !!----- RHEE -----
  !! Do we need to update mass_tot after the first kick?
  !!----------------
  use subsub_commons
  use amr_commons
  use cooling_module, ONLY:twopi
  use mpi_mod
  implicit none
  integer :: objind
  real(dp), dimension(1:2,1:ndim) :: bc_vv, bc_flux

  !! Local variables
  integer:: i, iter
  integer :: ii, ix, iy, iz
  integer :: xu, xd, yu, yd, zu, zd
  real(dp) :: mpinow
  real(dp) :: maxv, maxv_local
  logical :: isexit

  real(dp) :: subsub_dt, subsub_t0, subsub_dtmax
  integer :: subsub_nstep, subsub_ngrid2
  logical :: okay, skip_cg
  real(dp) :: rx, ry, rz

  real(dp) :: gconst, mtot, fourpiG, fact1
  real(dp) ::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

  !!CG
  real(dp) :: alpha, beta
  real(dp) :: rhs2, rr_old, rr_new, dx2inv, pLp, relres

  !!UPWIND
  real(dp), dimension(1:2, 1:ndim) :: vgdummy, fluxdummy
  real(dp) :: deltarho

  !!-----
  !! constants
  !!-----
  subsub_dt = dtold(levelmin)
  subsub_t0 = 0.0D0
  subsub_nstep = 0

  subsub_ncheck_cg(:) = 0
  subsub_tcheck_cg(:) = 0.0D0
    ! (1) - Poisson
    ! (2) - Force
    ! (3) - DT
    ! (4) - Kick
    ! (5) - Drift
  subsub_ngrid2 = subsub_ngrid**2
  dx2inv = 1.0D0 / (subsub_dx**2)

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  gconst = 6.67d-8*scale_d*scale_t**2

  mtot = subsub_obj(objind)%mass_tot
  fact1 = -gconst * mtot
  fourpiG = gconst * 2.0D0 * twopi



  isexit = .false.
  !!----- Main Iteration loop for the time loop
  !$omp parallel &
  !$omp private(i, ix, iy, iz, ii, rx, ry, rz, okay) &
  !$omp private(xu, xd, yu, yd, zu, zd, iter) &
  !$omp private(vgdummy, fluxdummy, deltarho)

  !! Main Time loop
  
  do
    if(isexit) exit
#ifndef WITHOUTMPI
    !$omp single
    mpinow = MPI_WTIME()
    !$omp end single
#endif

    !!!! Initial guess from the old potential
    !$omp single
    subsub_phi(:) = subsub_obj(objind)%phi(:)
    !$omp end single

    !!-----
    !! Compute Phi by CG
    !!-----

    !!!! Fix the boundary potential
    !$omp do collapse(3)
    do ix=1, subsub_ngrid
    do iy=1, subsub_ngrid
    do iz=1, subsub_ngrid
      if(ix.eq.1 .or. ix .eq. subsub_ngrid .or. &
         iy.eq.1 .or. iy .eq. subsub_ngrid .or. &
         iz.eq.1 .or. iz .eq. subsub_ngrid) then
        ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix

        subsub_phi(ii) = fact1 / sqrt(subsub_dd2(ii) + subsub_poisson_softening**2)

        !! copy to the old phi array
        !subsub_obj(objind)%phi(ii) = subsub_phi(ii)
      endif
    enddo
    enddo
    enddo
    !$omp end do

    !!!! Initialize CG arrays
    !$omp single
    !subsub_cgrhs(:) = 0.0D0
    !subsub_cgLphi(:) = 0.0D0
    !subsub_cgRes(:) = 0.0D0
    !subsub_cgp(:) = 0.0D0
    !subsub_cgLp(:) = 0.0D0
    skip_cg = .false.
    !subsub_fg(:,:) = 0.0D0
    rhs2 = 0.0D0
    rr_old = 0.0D0
    !$omp end single

    !!!! Build Initial CG arrays
    !$omp do collapse(3) reduction(+:rhs2, rr_old)
    do ix=1, subsub_ngrid
    do iy=1, subsub_ngrid
    do iz=1, subsub_ngrid
      ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix

      !! Zero BC
      if(ix.eq. 1 .or. ix .eq. subsub_ngrid .or. &
         iy.eq. 1 .or. iy .eq. subsub_ngrid .or. &
         iz.eq. 1 .or. iz .eq. subsub_ngrid) then
        subsub_cgrhs(ii) = 0.0D0
        subsub_cgLphi(ii) = 0.0D0
        subsub_cgRes(ii) = 0.0D0
        subsub_cgp(ii) = 0.0D0
        cycle
      endif

      subsub_cgrhs(ii) = -fourpiG * subsub_obj(objind)%hydro(ii,1)

      rhs2 = rhs2 + subsub_cgrhs(ii)*subsub_cgrhs(ii)

      xu = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix+1
      xd = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix-1
      yu = (iz-1)*subsub_ngrid2 + (iy  )*subsub_ngrid + ix
      yd = (iz-1)*subsub_ngrid2 + (iy-2)*subsub_ngrid + ix
      zu = (iz  )*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
      zd = (iz-2)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix

      subsub_cgLphi(ii) = (&
        6.0D0*subsub_phi(ii) - subsub_phi(xu) - subsub_phi(xd) &
        - subsub_phi(yu) - subsub_phi(yd) - subsub_phi(zu) - subsub_phi(zd) ) * dx2inv

      subsub_cgRes(ii) = subsub_cgrhs(ii) - subsub_cgLphi(ii)
      subsub_cgp(ii) = subsub_cgRes(ii)
      rr_old = rr_old + subsub_cgRes(ii)*subsub_cgRes(ii)
    enddo
    enddo
    enddo
    !$omp end do

    !$omp single
    if(sqrt(rr_old / max(rhs2, subsub_smallr)) .le. subsub_poisson_tolerance) then
      skip_cg = .true.
    endif
    !$omp end single


    !!!! CG Loop
    do iter=1, subsub_poisson_niter
      if(skip_cg) exit

      !$omp single
      pLp = 0.0D0
      rr_new = 0.0D0
      subsub_ncheck_cg(2) = subsub_ncheck_cg(2) + 1
      !$omp end single

      !!!!!! Compute Lp & pLp
      !$omp do collapse(3) reduction(+:pLp)
      do ix=1, subsub_ngrid
      do iy=1, subsub_ngrid
      do iz=1, subsub_ngrid
        ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix

        !! Zero BC
        if(ix.eq. 1 .or. ix .eq. subsub_ngrid .or. &
           iy.eq. 1 .or. iy .eq. subsub_ngrid .or. &
           iz.eq. 1 .or. iz .eq. subsub_ngrid) then
          subsub_cgLp(ii) = 0.0D0
          cycle
        endif

        xu = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix+1
        xd = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix-1
        yu = (iz-1)*subsub_ngrid2 + (iy  )*subsub_ngrid + ix
        yd = (iz-1)*subsub_ngrid2 + (iy-2)*subsub_ngrid + ix
        zu = (iz  )*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
        zd = (iz-2)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix

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
        skip_cg = .true.
        alpha = 0.0D0
      else
        alpha = rr_old / pLp
      endif
      !$omp end single

      !!!!!! Update
      !$omp do collapse(3) reduction(+:rr_new)
      do ix=1, subsub_ngrid
      do iy=1, subsub_ngrid
      do iz=1, subsub_ngrid
        ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix

        subsub_cgRes(ii) = subsub_cgRes(ii) - alpha*subsub_cgLp(ii)

        !! BC
        if(ix.eq. 1 .or. ix .eq. subsub_ngrid .or. &
           iy.eq. 1 .or. iy .eq. subsub_ngrid .or. &
           iz.eq. 1 .or. iz .eq. subsub_ngrid) then

          subsub_phi(ii) = fact1 / sqrt(subsub_dd2(ii) + subsub_poisson_softening**2)
          cycle
        endif


        rr_new = rr_new + subsub_cgRes(ii) * subsub_cgRes(ii)
        subsub_phi(ii) = subsub_phi(ii) + alpha * subsub_cgp(ii)
      enddo
      enddo
      enddo
      !$omp end do

      !$omp single
      relres = sqrt(rr_new / max(rhs2, subsub_smallr))
      if(relres .lt. subsub_poisson_tolerance) skip_cg = .true.
      beta = rr_new / max(rr_old, subsub_smallr)
      !$omp end single

      !! Update2
      !$omp do collapse(3)
      do ix=1, subsub_ngrid
      do iy=1, subsub_ngrid
      do iz=1, subsub_ngrid
        ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix

        !! Zero BC
        if(ix.eq. 1 .or. ix .eq. subsub_ngrid .or. &
           iy.eq. 1 .or. iy .eq. subsub_ngrid .or. &
           iz.eq. 1 .or. iz .eq. subsub_ngrid) then
          subsub_cgp(ii) = 0.0D0
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
    enddo !! CG Iteration loop

    !$omp barrier

#ifndef WITHOUTMPI
    !$omp single
    subsub_tcheck_cg(1) = subsub_tcheck_cg(1) + MPI_WTIME() - mpinow
    mpinow = MPI_WTIME()
    !$omp end single
#endif


    !!-----
    !! Compute g from phi
    !!-----
    !$omp do collapse(3)
    do ix=1, subsub_ngrid
    do iy=1, subsub_ngrid
    do iz=1, subsub_ngrid
      ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix

      xu = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix+1
      xd = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix-1
      yu = (iz-1)*subsub_ngrid2 + (iy  )*subsub_ngrid + ix
      yd = (iz-1)*subsub_ngrid2 + (iy-2)*subsub_ngrid + ix
      zu = (iz  )*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
      zd = (iz-2)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix

      if(ix.eq.1) then
        subsub_fg(ii,1) = - (subsub_phi(xu) - subsub_phi(ii)) / subsub_dx
      else if(ix .eq. subsub_ngrid) then
        subsub_fg(ii,1) = - (subsub_phi(ii) - subsub_phi(xd)) / subsub_dx
      else
        subsub_fg(ii,1) = - (subsub_phi(xu) - subsub_phi(xd)) / (2.0D0*subsub_dx)
      endif

      if(iy.eq.1) then
        subsub_fg(ii,2) = - (subsub_phi(yu) - subsub_phi(ii)) / subsub_dx
      else if(iy .eq. subsub_ngrid) then
        subsub_fg(ii,2) = - (subsub_phi(ii) - subsub_phi(yd)) / subsub_dx
      else
        subsub_fg(ii,2) = - (subsub_phi(yu) - subsub_phi(yd)) / (2.0D0*subsub_dx)
      endif

      if(iz.eq.1) then
        subsub_fg(ii,3) = - (subsub_phi(zu) - subsub_phi(ii)) / subsub_dx
      else if(iz .eq. subsub_ngrid) then
        subsub_fg(ii,3) = - (subsub_phi(ii) - subsub_phi(zd)) / subsub_dx
      else
        subsub_fg(ii,3) = - (subsub_phi(zu) - subsub_phi(zd)) / (2.0D0*subsub_dx)
      endif
    enddo
    enddo
    enddo
    !$omp end do

    !$omp barrier

#ifndef WITHOUTMPI
    !$omp single
    subsub_tcheck_cg(2) = subsub_tcheck_cg(2) + MPI_WTIME() - mpinow
    mpinow = MPI_WTIME()
    !$omp end single
#endif


    !!-----
    !! Compute dt based on velocity
    !!-----
    !$omp do
    do i=1, subsub_nn
      subsub_vg_old(i,1) = subsub_obj(objind)%hydro(i,2) / subsub_obj(objind)%hydro(i,1)
      subsub_vg_old(i,2) = subsub_obj(objind)%hydro(i,3) / subsub_obj(objind)%hydro(i,1)
      subsub_vg_old(i,3) = subsub_obj(objind)%hydro(i,4) / subsub_obj(objind)%hydro(i,1)
    enddo
    !$omp end do


    !$omp single
    maxv = maxval(abs(bc_vv))
    maxv_local = 0.0D0
    !$omp end single
  
  
    !$omp do reduction(max:maxv_local)
    do i=1, subsub_nn
      if(abs(subsub_vg_old(i,1)) .ge. maxv_local) then
        maxv_local = abs(subsub_vg_old(i,1))
      endif
  
      if(abs(subsub_vg_old(i,2)) .ge. maxv_local) then
        maxv_local = abs(subsub_vg_old(i,2))
      endif
  
      if(abs(subsub_vg_old(i,3)) .ge. maxv_local) then
        maxv_local = abs(subsub_vg_old(i,3))
      endif
    enddo
    !$omp end do
  
    !$omp single
    maxv = max(maxv, maxv_local)

    if(maxv .lt. subsub_smallr) then
      subsub_dtmax = dtold(levelmin) - subsub_t0
    else
      subsub_dtmax = subsub_cfl * (subsub_dx / maxv)
    endif
  
    do
      if(subsub_dt .lt. subsub_dtmax) exit
      subsub_dt = subsub_dt * 0.5D0
    enddo
    subsub_dt = min(subsub_dt, dtold(levelmin) - subsub_t0)
    !$omp end single

#ifndef WITHOUTMPI
    !$omp single
    subsub_tcheck_cg(3) = subsub_tcheck_cg(3) + MPI_WTIME() - mpinow
    mpinow = MPI_WTIME()
    !$omp end single
#endif

    !$omp barrier

    !!-----
    !! Kick: Velocity update by 1/2
    !!-----
    !$omp do
    do i=1, subsub_nn
      subsub_vg_old(i,1) = subsub_vg_old(i,1) + subsub_fg(i,1) * (subsub_dt/2.0D0)
      subsub_vg_old(i,2) = subsub_vg_old(i,2) + subsub_fg(i,2) * (subsub_dt/2.0D0)
      subsub_vg_old(i,3) = subsub_vg_old(i,3) + subsub_fg(i,3) * (subsub_dt/2.0D0)
      !subsub_obj(objind)%hydro(i,2) = subsub_obj(objind)%hydro(i,2) + &
      !  subsub_fg(i,1) * (subsub_dt/2.0D0) * subsub_obj(objind)%hydro(i,1)
      !subsub_obj(objind)%hydro(i,3) = subsub_obj(objind)%hydro(i,3) + &
      !  subsub_fg(i,2) * (subsub_dt/2.0D0) * subsub_obj(objind)%hydro(i,1)
      !subsub_obj(objind)%hydro(i,4) = subsub_obj(objind)%hydro(i,4) + &
      !  subsub_fg(i,3) * (subsub_dt/2.0D0) * subsub_obj(objind)%hydro(i,1)
    enddo
    !$omp end do

    !$omp barrier

#ifndef WITHOUTMPI
    !$omp single
    subsub_tcheck_cg(4) = subsub_tcheck_cg(4) + MPI_WTIME() - mpinow
    mpinow = MPI_WTIME()
    !$omp end single
#endif


    !!-----
    !! Drift: Rho update by 1
    !!-----

    !! vgdummy, fluxdummy : (1:2,:) down , up
    !!                      (:, 1:3) x, y, z
    !$omp single
    subsub_rho_old(:) = subsub_obj(objind)%hydro(:,1)
    !$omp end single

    !$omp do collapse(3)
    do ix=1, subsub_ngrid
    do iy=1, subsub_ngrid
    do iz=1, subsub_ngrid
      ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix

      xu = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix+1
      xd = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix-1
      yu = (iz-1)*subsub_ngrid2 + (iy  )*subsub_ngrid + ix
      yd = (iz-1)*subsub_ngrid2 + (iy-2)*subsub_ngrid + ix
      zu = (iz  )*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
      zd = (iz-2)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix

      !! X
      if(ix .ne. subsub_ngrid) then
        vgdummy(2,1) =  (subsub_vg_old(xu,1) + subsub_vg_old(ii,1)) / 2.0D0
      else
        vgdummy(2,1) = bc_vv(2,1)
      endif

      if(ix .ne. 1) then
        vgdummy(1,1) = (subsub_vg_old(xd,1) + subsub_vg_old(ii,1)) / 2.0D0
      else
        vgdummy(1,1) = bc_vv(1,1)
      endif

      !! Y velocity
      if(iy .ne. subsub_ngrid) then
        vgdummy(2,2) =  (subsub_vg_old(yu,2) + subsub_vg_old(ii,2)) / 2.0D0
      else
        vgdummy(2,2) = bc_vv(2,2)
      endif

      if(iy .ne. 1) then
        vgdummy(1,2) = (subsub_vg_old(yd,2) + subsub_vg_old(ii,2)) / 2.0D0
      else
        vgdummy(1,2) = bc_vv(1,2)
      endif

      !! Z velocity
      if(iz .ne. subsub_ngrid) then
        vgdummy(2,3) =  (subsub_vg_old(zu,3) + subsub_vg_old(ii,3)) / 2.0D0
      else
        vgdummy(2,3) = bc_vv(2,3)
      endif

      if(iz .ne. 1) then
        vgdummy(1,3) = (subsub_vg_old(zd,3) + subsub_vg_old(ii,3)) / 2.0D0
      else
        vgdummy(1,3) = bc_vv(1,3)
      endif

      !! Flux update (X-up)
      if(vgdummy(2,1) .gt. 0) then
        fluxdummy(2,1) = subsub_obj(objind)%hydro(ii,1)*vgdummy(2,1)
      else
        if(ix .ne. subsub_ngrid) then
          fluxdummy(2,1) = subsub_obj(objind)%hydro(xu,1)*vgdummy(2,1)
        else
          fluxdummy(2,1) = bc_flux(2,1)
        endif
      endif

      !! (X-down)
      if(vgdummy(1,1) .lt. 0) then
        fluxdummy(1,1) = subsub_obj(objind)%hydro(ii,1)*vgdummy(1,1)
      else
        if(ix .ne. 1) then
          fluxdummy(1,1) = subsub_obj(objind)%hydro(xd,1)*vgdummy(1,1)
        else
          fluxdummy(1,1) = bc_flux(1,1)
        endif
      endif

      !! Flux update (Y-up)
      if(vgdummy(2,2) .gt. 0) then
        fluxdummy(2,2) = subsub_obj(objind)%hydro(ii,1)*vgdummy(2,2)
      else
        if(iy .ne. subsub_ngrid) then
          fluxdummy(2,2) = subsub_obj(objind)%hydro(yu,1)*vgdummy(2,2)
        else
          fluxdummy(2,2) = bc_flux(2,2)
        endif
      endif

      !! (Y-down)
      if(vgdummy(1,2) .lt. 0) then
        fluxdummy(1,2) = subsub_obj(objind)%hydro(ii,1)*vgdummy(1,2)
      else
        if(iy .ne. 1) then
          fluxdummy(1,2) = subsub_obj(objind)%hydro(yd,1)*vgdummy(1,2)
        else
          fluxdummy(1,2) = bc_flux(1,2)
        endif
      endif

      !! Flux update (Z-up)
      if(vgdummy(2,3) .gt. 0) then
        fluxdummy(2,3) = subsub_obj(objind)%hydro(ii,1)*vgdummy(2,3)
      else
        if(iz .ne. subsub_ngrid) then
          fluxdummy(2,3) = subsub_obj(objind)%hydro(zu,1)*vgdummy(2,3)
        else
          fluxdummy(2,3) = bc_flux(2,3)
        endif
      endif

      !! (Z-down)
      if(vgdummy(1,3) .lt. 0) then
        fluxdummy(1,3) = subsub_obj(objind)%hydro(ii,1)*vgdummy(1,3)
      else
        if(iz .ne. 1) then
          fluxdummy(1,3) = subsub_obj(objind)%hydro(zd,1)*vgdummy(1,3)
        else
          fluxdummy(1,3) = bc_flux(1,3)
        endif
      endif

      !! (update density)
      deltarho = fluxdummy(2,1) - fluxdummy(1,1) + &
        fluxdummy(2,2) - fluxdummy(1,2) + &
        fluxdummy(2,3) - fluxdummy(1,3)

      deltarho = deltarho * subsub_dt / subsub_dx

      subsub_rho_old(ii) = max(subsub_rho_old(ii) - deltarho, subsub_densityfloor)
    enddo
    enddo
    enddo
    !$omp end do

    !$omp barrier

    !$omp single
    subsub_obj(objind)%hydro(:,1) = subsub_rho_old(:)
 
#ifndef WITHOUTMPI
    subsub_tcheck_cg(5) = subsub_tcheck_cg(5) + MPI_WTIME() - mpinow
    mpinow = MPI_WTIME()
#endif
    !$omp end single



    !!-----
    !! Poisson by the updated density
    !!-----
    !!!! Fix the boundary potential
    
    !$omp do collapse(3)
    do ix=1, subsub_ngrid
    do iy=1, subsub_ngrid
    do iz=1, subsub_ngrid
      if(ix.eq.1 .or. ix .eq. subsub_ngrid .or. &
         iy.eq.1 .or. iy .eq. subsub_ngrid .or. &
         iz.eq.1 .or. iz .eq. subsub_ngrid) then
        ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix

        subsub_phi(ii) = fact1 / sqrt(subsub_dd2(ii) + subsub_poisson_softening**2)

        !! copy to the old phi array
        !subsub_obj(objind)%phi(ii) = subsub_phi(ii)
      endif
    enddo
    enddo
    enddo
    !$omp end do

    !!!! Initialize CG arrays
    !$omp single
    !subsub_cgrhs(:) = 0.0D0
    !subsub_cgLphi(:) = 0.0D0
    !subsub_cgRes(:) = 0.0D0
    !subsub_cgp(:) = 0.0D0
    !subsub_cgLp(:) = 0.0D0
    !subsub_fg(:,:) = 0.0D0
    skip_cg = .false.
    rhs2 = 0.0D0
    rr_old = 0.0D0
    !$omp end single

    !!!! Build Initial CG arrays
    !$omp do collapse(3) reduction(+:rhs2, rr_old)
    do ix=1, subsub_ngrid
    do iy=1, subsub_ngrid
    do iz=1, subsub_ngrid
      ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix

      !! Zero BC
      if(ix.eq. 1 .or. ix .eq. subsub_ngrid .or. &
         iy.eq. 1 .or. iy .eq. subsub_ngrid .or. &
         iz.eq. 1 .or. iz .eq. subsub_ngrid) then
        subsub_cgrhs(ii) = 0.0D0
        subsub_cgLphi(ii) = 0.0D0
        subsub_cgRes(ii) = 0.0D0
        subsub_cgp(ii) = 0.0D0
        cycle
      endif

      subsub_cgrhs(ii) = -fourpiG * subsub_obj(objind)%hydro(ii,1)

      rhs2 = rhs2 + subsub_cgrhs(ii)*subsub_cgrhs(ii)

      xu = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix+1
      xd = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix-1
      yu = (iz-1)*subsub_ngrid2 + (iy  )*subsub_ngrid + ix
      yd = (iz-1)*subsub_ngrid2 + (iy-2)*subsub_ngrid + ix
      zu = (iz  )*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
      zd = (iz-2)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix

      subsub_cgLphi(ii) = (&
        6.0D0*subsub_phi(ii) - subsub_phi(xu) - subsub_phi(xd) &
        - subsub_phi(yu) - subsub_phi(yd) - subsub_phi(zu) - subsub_phi(zd) ) * dx2inv

      subsub_cgRes(ii) = subsub_cgrhs(ii) - subsub_cgLphi(ii)
      subsub_cgp(ii) = subsub_cgRes(ii)
      rr_old = rr_old + subsub_cgRes(ii)*subsub_cgRes(ii)
    enddo
    enddo
    enddo
    !$omp end do

    !$omp single
    if(sqrt(rr_old / max(rhs2, subsub_smallr)) .le. subsub_poisson_tolerance) then
      skip_cg = .true.
    endif
    !$omp end single


    !!!! CG Loop
    do iter=1, subsub_poisson_niter
      if(skip_cg) exit

      !$omp single
      pLp = 0.0D0
      rr_new = 0.0D0
      subsub_ncheck_cg(3) = subsub_ncheck_cg(3) + 1
      !$omp end single

      !!!!!! Compute Lp & pLp
      !$omp do collapse(3) reduction(+:pLp)
      do ix=1, subsub_ngrid
      do iy=1, subsub_ngrid
      do iz=1, subsub_ngrid
        ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix

        !! Zero BC
        if(ix.eq. 1 .or. ix .eq. subsub_ngrid .or. &
           iy.eq. 1 .or. iy .eq. subsub_ngrid .or. &
           iz.eq. 1 .or. iz .eq. subsub_ngrid) then
          subsub_cgLp(ii) = 0.0D0
          cycle
        endif

        xu = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix+1
        xd = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix-1
        yu = (iz-1)*subsub_ngrid2 + (iy  )*subsub_ngrid + ix
        yd = (iz-1)*subsub_ngrid2 + (iy-2)*subsub_ngrid + ix
        zu = (iz  )*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
        zd = (iz-2)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix

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
        skip_cg = .true.
        alpha = 0.0D0
      else
        alpha = rr_old / pLp
      endif
      !$omp end single

      !!!!!! Update
      !$omp do collapse(3) reduction(+:rr_new)
      do ix=1, subsub_ngrid
      do iy=1, subsub_ngrid
      do iz=1, subsub_ngrid
        ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix

        subsub_cgRes(ii) = subsub_cgRes(ii) - alpha*subsub_cgLp(ii)

        !! BC
        if(ix.eq. 1 .or. ix .eq. subsub_ngrid .or. &
           iy.eq. 1 .or. iy .eq. subsub_ngrid .or. &
           iz.eq. 1 .or. iz .eq. subsub_ngrid) then

          subsub_phi(ii) = fact1 / sqrt(subsub_dd2(ii) + subsub_poisson_softening**2)
          cycle
        endif


        rr_new = rr_new + subsub_cgRes(ii) * subsub_cgRes(ii)
        subsub_phi(ii) = subsub_phi(ii) + alpha * subsub_cgp(ii)
      enddo
      enddo
      enddo
      !$omp end do

      !$omp single
      relres = sqrt(rr_new / max(rhs2, subsub_smallr))
      if(relres .lt. subsub_poisson_tolerance) skip_cg = .true.
      beta = rr_new / max(rr_old, subsub_smallr)
      !$omp end single

      !! Update2
      !$omp do collapse(3)
      do ix=1, subsub_ngrid
      do iy=1, subsub_ngrid
      do iz=1, subsub_ngrid
        ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix

        !! Zero BC
        if(ix.eq. 1 .or. ix .eq. subsub_ngrid .or. &
           iy.eq. 1 .or. iy .eq. subsub_ngrid .or. &
           iz.eq. 1 .or. iz .eq. subsub_ngrid) then
          subsub_cgp(ii) = 0.0D0
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
    enddo !! CG Iteration loop

    !$omp barrier

    !$omp single
    subsub_obj(objind)%phi(:) = subsub_phi(:)
#ifndef WITHOUTMPI
    subsub_tcheck_cg(1) = subsub_tcheck_cg(1) + MPI_WTIME() - mpinow
    mpinow = MPI_WTIME()
#endif
    !$omp end single

    !!-----
    !! Compute g from phi
    !!-----
    !$omp do collapse(3)
    do ix=1, subsub_ngrid
    do iy=1, subsub_ngrid
    do iz=1, subsub_ngrid
      ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix

      xu = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix+1
      xd = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix-1
      yu = (iz-1)*subsub_ngrid2 + (iy  )*subsub_ngrid + ix
      yd = (iz-1)*subsub_ngrid2 + (iy-2)*subsub_ngrid + ix
      zu = (iz  )*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
      zd = (iz-2)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix

      if(ix.eq.1) then
        subsub_fg(ii,1) = - (subsub_phi(xu) - subsub_phi(ii)) / subsub_dx
      else if(ix .eq. subsub_ngrid) then
        subsub_fg(ii,1) = - (subsub_phi(ii) - subsub_phi(xd)) / subsub_dx
      else
        subsub_fg(ii,1) = - (subsub_phi(xu) - subsub_phi(xd)) / (2.0D0*subsub_dx)
      endif

      if(iy.eq.1) then
        subsub_fg(ii,2) = - (subsub_phi(yu) - subsub_phi(ii)) / subsub_dx
      else if(iy .eq. subsub_ngrid) then
        subsub_fg(ii,2) = - (subsub_phi(ii) - subsub_phi(yd)) / subsub_dx
      else
        subsub_fg(ii,2) = - (subsub_phi(yu) - subsub_phi(yd)) / (2.0D0*subsub_dx)
      endif

      if(iz.eq.1) then
        subsub_fg(ii,3) = - (subsub_phi(zu) - subsub_phi(ii)) / subsub_dx
      else if(iz .eq. subsub_ngrid) then
        subsub_fg(ii,3) = - (subsub_phi(ii) - subsub_phi(zd)) / subsub_dx
      else
        subsub_fg(ii,3) = - (subsub_phi(zu) - subsub_phi(zd)) / (2.0D0*subsub_dx)
      endif
    enddo
    enddo
    enddo
    !$omp end do

    !$omp barrier

#ifndef WITHOUTMPI
    !$omp single
    subsub_tcheck_cg(2) = subsub_tcheck_cg(2) + MPI_WTIME() - mpinow
    mpinow = MPI_WTIME()
    !$omp end single
#endif

    !!-----
    !! Second Kick: Velocity update by 1/2
    !!-----
    !$omp do
    do i=1, subsub_nn
      subsub_obj(objind)%hydro(i,2) = (subsub_vg_old(i,1) + subsub_fg(i,1) * (subsub_dt/2.0D0)) * subsub_obj(objind)%hydro(i,1)
      subsub_obj(objind)%hydro(i,3) = (subsub_vg_old(i,2) + subsub_fg(i,2) * (subsub_dt/2.0D0)) * subsub_obj(objind)%hydro(i,1)
      subsub_obj(objind)%hydro(i,4) = (subsub_vg_old(i,3) + subsub_fg(i,3) * (subsub_dt/2.0D0)) * subsub_obj(objind)%hydro(i,1)
    enddo
    !$omp end do

    !$omp barrier

#ifndef WITHOUTMPI 
    !$omp single
    subsub_tcheck_cg(4) = subsub_tcheck_cg(4) + MPI_WTIME() - mpinow
    mpinow = MPI_WTIME()
    !$omp end single
#endif

    !! loop-control 
    !$omp single
    subsub_t0 = subsub_t0 + subsub_dt
    subsub_nstep = subsub_nstep + 1
    
    if(subsub_t0 .ge. dtold(levelmin))then
      isexit = .true.
    endif
    !$omp end single
  enddo !! Main Time loop
  !$omp end parallel


  subsub_ncheck_cg(1) = subsub_nstep
end subroutine subsub_cgkdk
!################################################################
!################################################################
!################################################################
!################################################################