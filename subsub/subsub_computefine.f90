!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_computefine(objind)
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
  integer :: objind

  !!----- Local Variables
  !integer :: i
  integer :: subsub_nstep
  
  real(dp) :: subsub_dt, subsub_t0
  !real(dp), dimension(1:2,1:ndim) :: bc_vv, bc_flux
  !real(dp), dimension(1:2,1:ndim) :: bc_rho, bc_p
  !real(dp), dimension(1:10) :: varr
  !real(dp) :: m_in, m_out
  !real(dp) :: rho_ave, tff, vff, mdot_face
  !real(dp) :: delMtot, rho_in

  !real(dp) :: scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  !real(dp) :: threepi2, fourpi

  !real(dp) :: tcheck(1:10)

  !!----- 
  subsub_t0 = 0.0D0
  subsub_nstep = 0

  !!----- Start Subcycle
  subsub_dt = dtold(levelmin)

  !!----- By Conjugrate-Gradient & KDK & Godunov
  call subsub_cgkdk(objind)


end subroutine subsub_computefine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_cgkdk(objind)
  !!----- RHEE -----
  !! Do we need to update mass_tot after the first kick?
  !!----------------
  use amr_commons
  use subsub_commons
  use subsub_parameters
  use pm_commons
  use cooling_module, ONLY:twopi
  use hydro_commons, ONLY: gamma
  use mpi_mod
  implicit none
  integer :: objind
  

  !! Local variables
  integer:: i, j, iter, sinkind, ivar
  integer :: ii, ix, iy, iz
  real(dp) :: mpinow
  real(dp) :: maxv, maxv_local
  real(dp) :: newMgas, newMBH
  logical :: isexit

  real(dp) :: subsub_dt, subsub_t0, subsub_dtmax
  integer :: subsub_nstep
  real(dp) :: gconst, mtot_now, mtot_old, mtot_new, mtot_cell, fourpiG, fact1, fact2, pi
  !real(dp) ::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  
  real(dp) :: delMass, delMBH, delMgas
  real(dp) :: mass_predicted, rho_ave, tff, vff, threepi2, rho_needed
  !! HYDRO
  real(dp) :: rho, u, v, w, r, e, vx, vy, vz, ekin, p, cs, vxn, vyn, vzn, ekin_floor
  real(dp), dimension(1:subsub_nhydro) :: Utmp
  real(dp) :: debugt

  real(dp) :: subsub_hydrobc_old(1:2, 1:ndim, 1:subsub_nhydro)

  debugt=MPI_WTIME()
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

  maxv = 0.0D0
  

  !call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  !gconst = 6.67d-8*scale_d*scale_t**2

  gconst=1d0
  pi=twopi/2d0
  threepi2=3.0d0*ACOS(-1.0d0)**2
  if(cosmo)gconst=3d0/8d0/pi*omega_m*aexp

  mtot_old = subsub_obj(objind)%mass_tot
  mtot_now = mtot_old
  mtot_new = 0.0D0

  mtot_cell = subsub_obj(objind)%mass_cell

  fact1 = -gconst * mtot_old
  fact2 = -gconst * subsub_obj(objind)%sink_mass
  fourpiG = gconst * 2.0D0 * twopi



  !!-----
  !! Initial BC Set
  !!-----
  sinkind = subsub_obj(objind)%sink_ind

  delMbh = msink(sinkind) - subsub_obj(objind)%sink_mass
  delMgas= subsub_obj(objind)%mass_cell - subsub_obj(objind)%mass_tot !! cell mass is already updated
  delMass = delMgas + delMbh

  subsub_hydrobc(1,1,:) = subsub_obj(objind)%uold(1,:) !! X-
  subsub_hydrobc(2,1,:) = subsub_obj(objind)%uold(2,:) !! X+
  subsub_hydrobc(1,2,:) = subsub_obj(objind)%uold(3,:) !! Y-
  subsub_hydrobc(2,2,:) = subsub_obj(objind)%uold(4,:) !! Y+
  subsub_hydrobc(1,3,:) = subsub_obj(objind)%uold(5,:) !! Z-
  subsub_hydrobc(2,3,:) = subsub_obj(objind)%uold(6,:) !! Z+

  
  !! DEBUG MODE FOR SELF-GRAVITY TEST
  !! )) DEBUGG GRAV((         <- this is for grep
  if(subsub_dev_gravonly .eq. 1)then
    do i=1, 2
      do j=1,3
        subsub_hydrobc(i,j,1) = subsub_dfloor
        subsub_hydrobc(i,j,2) = 0.0D0
        subsub_hydrobc(i,j,3) = 0.0D0
        subsub_hydrobc(i,j,4) = 0.0D0
        subsub_hydrobc(i,j,5) = 0.0D0
      enddo
    enddo

    fact2 = 0.0D0
  endif

  !! DEBUG MODE FOR SELF-GRAVITY TEST
  !! )) DEBUGG HYDRO((         <- this is for grep
  if(subsub_dev_hydroonly .eq. 1)then
    fact2 = 0.0D0
  endif
  
  !! DEBUG MODE FOR SELF-GRAVITY TEST
  !! )) DEBUGG GRAV+HYDRO((         <- this is for grep
  if(subsub_dev_gravhydro .eq. 1)then
    do i=1, 2
      do j=1,3
        subsub_hydrobc(i,j,1) = subsub_dfloor
        subsub_hydrobc(i,j,2) = 0.0D0
        subsub_hydrobc(i,j,3) = 0.0D0
        subsub_hydrobc(i,j,4) = 0.0D0
        subsub_hydrobc(i,j,5) = 0.0D0
      enddo
    enddo

    fact2 = 0.0D0
    subsub_debugtag = -1
  endif

  !! DEBUG MODE FOR SELF-GRAVITY TEST
  !! )) DEBUGG GRAV+HYDRO+BH((         <- this is for grep
  if(subsub_dev_gravhydrobh .eq. 1)then
    do i=1, 2
      do j=1,3
        subsub_hydrobc(i,j,1) = subsub_dfloor
        subsub_hydrobc(i,j,2) = 0.0D0
        subsub_hydrobc(i,j,3) = 0.0D0
        subsub_hydrobc(i,j,4) = 0.0D0
        subsub_hydrobc(i,j,5) = 0.0D0
      enddo
    enddo
  endif
  

  do i=1,2
    do j=1,3
      call subsub_enforce_floors(subsub_hydrobc(i,j,:))

      !! save the original bc
      subsub_hydrobc_old(i,j,:) = subsub_hydrobc(i,j,:)
    enddo
  enddo

  !!----- Main Iteration loop for the time loop
  isexit = .false.
  !$omp parallel &
  !$omp private(i, ix, iy, iz, ii) &
  !$omp private(iter) &
  !$omp private(rho, u, v, w, e, vx, vy, vz, ekin, p, cs, vxn, vyn, vzn, ekin_floor)
  !! Main Time loop
  do
    !$omp barrier
    if(isexit) exit
#ifndef WITHOUTMPI
    !$omp single
    mpinow = MPI_WTIME()
    !$omp end single
#endif

!!$omp single
!if(subsub_obj(objind)%sink_id .eq. 2) then
!  write(*,*) myid
!  write(*,*) subsub_obj(objind)%mass_cell, mtot_now
!  write(*,*)'     niter        = ', subsub_nstep
!endif
!!$omp end single

!if(mpinow-debugt .gt. 120.)then
!  !$omp single
!  write(*,*)'     myid         = ', myid
!  write(*,*)'     poisson      = ', subsub_tcheck_cg(1)
!  write(*,*)'     force        = ', subsub_tcheck_cg(2)
!  write(*,*)'     dt           = ', subsub_tcheck_cg(3)
!  write(*,*)'     kick         = ', subsub_tcheck_cg(4)
!  write(*,*)'     drift        = ', subsub_tcheck_cg(5)
!  write(*,*)'     poisson_nall = ', subsub_ncheck_cg(1)
!  write(*,*)'     poisson_n1   = ', subsub_ncheck_cg(2)
!  write(*,*)'     poisson_n2   = ', subsub_ncheck_cg(3)
!  write(*,*)'     poisson error= ', subsub_relres
!  write(*,*)'     dt/dt        = ', dtold(levelmin)/subsub_dt
!  write(*,*)'     maxv         = ', maxv
!  write(*,*)'     maxrho       = ', maxval(subsub_hydro(:,1))
!  write(*,*)'     minrho       = ', minval(subsub_hydro(:,1))
!  write(*,*)'     niter        = ', subsub_nstep
!  write(*,*)' '
!  !$omp end single
!endif

!if(subsub_nstep .gt. 10000) then
!if(mod(subsub_nstep,1000) .eq. 0) then
!if(myid .eq. 25) then
!!$omp single
!  write(*,*)'  ==================', myid
!  write(*,*)'     poisson      = ', subsub_tcheck_cg(1)
!  write(*,*)'     force        = ', subsub_tcheck_cg(2)
!  write(*,*)'     dt           = ', subsub_tcheck_cg(3)
!  write(*,*)'     kick         = ', subsub_tcheck_cg(4)
!  write(*,*)'     drift        = ', subsub_tcheck_cg(5)
!  write(*,*)'     poisson_nall = ', subsub_ncheck_cg(1)
!  write(*,*)'     poisson_n1   = ', subsub_ncheck_cg(2)
!  write(*,*)'     poisson_n2   = ', subsub_ncheck_cg(3)
!  write(*,*)'     poisson error= ', subsub_relres
!  write(*,*)'     dt/dt        = ', dtold(levelmin)/subsub_dt
!  write(*,*)'     dt           = ', subsub_dt
!  write(*,*)'     maxv         = ', maxv
!  if(subsub_nstep .gt. 1) then
!  write(*,*)'     maxrho       = ', maxval(subsub_hydro(:,1))
!  write(*,*)'     minrho       = ', minval(subsub_hydro(:,1))
!  write(*,*)'     maxPx       = ', maxval(subsub_hydro(:,2))
!  write(*,*)'     maxE       = ', maxval(subsub_hydro(:,5))
!  endif
!  write(*,*)'     niter        = ', subsub_nstep
!  write(*,*)' '
!!$omp end single
!endif
!endif

!!-----------------------------------------------------------------
!! Save old arrays
!!-----------------------------------------------------------------

    !$omp do private(ivar, Utmp)
    do i=1, subsub_nn
      subsub_phi(i) = subsub_obj(objind)%phi(i)

      do ivar=1,subsub_nhydro
        Utmp(ivar) = subsub_obj(objind)%hydro(i,ivar)
      enddo
      call subsub_enforce_floors(Utmp)

      do ivar=1,subsub_nhydro
        subsub_hydro(i,ivar) = Utmp(ivar)
      enddo
    enddo
    !$omp end do


!$omp single
if(subsub_obj(objind)%sink_id .eq. 3)then
  subsub_debugtag = 1
  write(*,*)'     niter        = ', subsub_nstep
  write(*,*)'     maxrho       = ', maxval(subsub_hydro(:,1))
  write(*,*)'     minrho       = ', minval(subsub_hydro(:,1))

  write(*,*)'     maxE       = ', maxval(subsub_hydro(:,5))
  write(*,*)'     minE       = ', minval(subsub_hydro(:,5))

  write(*,*)'     maxPx       = ', maxval(subsub_hydro(:,2))
  write(*,*)'     minPx       = ', minval(subsub_hydro(:,2))

  write(*,*)'     maxPy       = ', maxval(subsub_hydro(:,3))
  write(*,*)'     minPy       = ', minval(subsub_hydro(:,3))

  write(*,*)'     maxPz       = ', maxval(subsub_hydro(:,4))
  write(*,*)'     minPz       = ', minval(subsub_hydro(:,4))
endif
!$omp end single

!!-----------------------------------------------------------------
!! Compute Phi by CG (rho_n -> phi_n)
!!  if subsub_nstep>1, use the previous subsub_phi
!!-----------------------------------------------------------------
    if(subsub_nstep .eq. 0)then
      call subsub_poisson_cg(mtot_now)
    endif


    
#ifndef WITHOUTMPI
    !$omp single
    subsub_tcheck_cg(1) = subsub_tcheck_cg(1) + MPI_WTIME() - mpinow
    mpinow = MPI_WTIME()
    !$omp end single
#endif

!!-----------------------------------------------------------------
!! Compute g_n+1
!!-----------------------------------------------------------------
    call subsub_poisson_g(fact2)

#ifndef WITHOUTMPI
    !$omp single
    subsub_tcheck_cg(2) = subsub_tcheck_cg(2) + MPI_WTIME() - mpinow
    mpinow = MPI_WTIME()
    !$omp end single
#endif


!!-----------------------------------------------------------------
!! Determine dt by CFL (updated velocity + sound speed)
!!-----------------------------------------------------------------
    !$omp single
    maxv = 0.0D0
    maxv = max(maxv, maxval( abs(subsub_hydrobc(1,1,2:4)) / max(subsub_hydrobc(1,1,1), subsub_dfloor) ))
    maxv = max(maxv, maxval( abs(subsub_hydrobc(2,1,2:4)) / max(subsub_hydrobc(2,1,1), subsub_dfloor) ))

    maxv = max(maxv, maxval( abs(subsub_hydrobc(1,2,2:4)) / max(subsub_hydrobc(1,2,1), subsub_dfloor) ))
    maxv = max(maxv, maxval( abs(subsub_hydrobc(2,2,2:4)) / max(subsub_hydrobc(2,2,1), subsub_dfloor) ))

    maxv = max(maxv, maxval( abs(subsub_hydrobc(1,3,2:4)) / max(subsub_hydrobc(1,3,1), subsub_dfloor) ))
    maxv = max(maxv, maxval( abs(subsub_hydrobc(2,3,2:4)) / max(subsub_hydrobc(2,3,1), subsub_dfloor) ))

    maxv_local = 0.0D0
    !$omp end single
  
  
    !$omp do reduction(max:maxv_local)
    do i=1, subsub_nn
      rho = max(subsub_hydro(i,1), subsub_dfloor)
      if(rho .le. subsub_dfloor) cycle

      u = subsub_hydro(i,2)
      v = subsub_hydro(i,3)
      w = subsub_hydro(i,4)
      e = subsub_hydro(i,5)

      vx = u/rho
      vy = v/rho
      vz = w/rho

      ekin = 0.5D0 * (u**2 + v**2 + w**2) / rho
      p = max( (gamma-1.0D0)*(e - ekin), subsub_pfloor )

      if(rho .le. subsub_dfloor) then
        p = min(p,subsub_pfloor)
      endif

      cs = sqrt(gamma * p / rho)

      maxv_local = max(maxv_local, abs(vx) + cs)
      maxv_local = max(maxv_local, abs(vy) + cs)
      maxv_local = max(maxv_local, abs(vz) + cs)
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

    do 
      if(subsub_dt * 2.0D0 .gt. subsub_dtmax) exit
      subsub_dt = 2.0D0 * subsub_dt
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

!!-----------------------------------------------------------------
!! Adjust BC To match mass conservation
!!-----------------------------------------------------------------
    !$omp single
    if(subsub_nstep .lt. 0) then
!!----- RHEE -----
!! turn off at the moment
!!----------------
      mass_predicted = (mtot_cell - mtot_old) / dtold(levelmin) * subsub_t0 + mtot_old

      
      
      if(mass_predicted .gt. subsub_massinfac * mtot_now)then
        !! now the total cell mass is less than the predicted mass
        !! hence, turn on the mass inflwo with the free-fall normal

        rho_ave = mtot_now / (subsub_boxlen**3)
        tff = sqrt(threepi2/8./fourpiG/rho_ave)
        vff = subsub_boxlen/tff

        rho_needed = 6.0D0 * (mass_predicted - mtot_now) / vff / subsub_dt / subsub_boxlen**2

        if(rho_needed .gt. subsub_dfloor) then
          r = subsub_hydrobc(1,1,1)
          u = subsub_hydrobc(1,1,2) / r
          v = subsub_hydrobc(1,1,3) / r
          w = subsub_hydrobc(1,1,4) / r
          p = max((gamma - 1.0D0)*(subsub_hydrobc(1,1,5) - 0.5D0*r*(u*u + v*v + w*w)), subsub_pfloor)
  
          subsub_hydrobc(1,1,1) = rho_needed
          subsub_hydrobc(1,1,2) = subsub_hydrobc(1,1,1)*vff
          subsub_hydrobc(1,1,3) = 0.0D0
          subsub_hydrobc(1,1,4) = 0.0D0
          subsub_hydrobc(1,1,5) = p/(gamma-1.0D0) * 0.5D0*rho_needed*vff*vff
  
  
          r = subsub_hydrobc(2,1,1)
          u = subsub_hydrobc(2,1,2) / r
          v = subsub_hydrobc(2,1,3) / r
          w = subsub_hydrobc(2,1,4) / r
          p = max((gamma - 1.0D0)*(subsub_hydrobc(2,1,5) - 0.5D0*r*(u*u + v*v + w*w)), subsub_pfloor)
  
          subsub_hydrobc(2,1,1) = rho_needed
          subsub_hydrobc(2,1,2) = -subsub_hydrobc(2,1,1)*vff
          subsub_hydrobc(2,1,3) = 0.0D0
          subsub_hydrobc(2,1,4) = 0.0D0
          subsub_hydrobc(2,1,5) = p/(gamma-1.0D0) * 0.5D0*rho_needed*vff*vff
  
  
          r = subsub_hydrobc(1,2,1)
          u = subsub_hydrobc(1,2,2) / r
          v = subsub_hydrobc(1,2,3) / r
          w = subsub_hydrobc(1,2,4) / r
          p = max((gamma - 1.0D0)*(subsub_hydrobc(1,2,5) - 0.5D0*r*(u*u + v*v + w*w)), subsub_pfloor)
  
          subsub_hydrobc(1,2,1) = rho_needed
          subsub_hydrobc(1,2,2) = 0.0D0
          subsub_hydrobc(1,2,3) = subsub_hydrobc(1,2,1)*vff
          subsub_hydrobc(1,2,4) = 0.0D0
          subsub_hydrobc(1,2,5) = p/(gamma-1.0D0) * 0.5D0*rho_needed*vff*vff
  
  
  
          r = subsub_hydrobc(2,2,1)
          u = subsub_hydrobc(2,2,2) / r
          v = subsub_hydrobc(2,2,3) / r
          w = subsub_hydrobc(2,2,4) / r
          p = max((gamma - 1.0D0)*(subsub_hydrobc(2,2,5) - 0.5D0*r*(u*u + v*v + w*w)), subsub_pfloor)
  
          subsub_hydrobc(2,2,1) = rho_needed
          subsub_hydrobc(2,2,2) = 0.0D0
          subsub_hydrobc(2,2,3) = -subsub_hydrobc(2,2,1)*vff
          subsub_hydrobc(2,2,4) = 0.0D0
          subsub_hydrobc(2,2,5) = p/(gamma-1.0D0) * 0.5D0*rho_needed*vff*vff
  
  
  
          r = subsub_hydrobc(1,3,1)
          u = subsub_hydrobc(1,3,2) / r
          v = subsub_hydrobc(1,3,3) / r
          w = subsub_hydrobc(1,3,4) / r
          p = max((gamma - 1.0D0)*(subsub_hydrobc(1,3,5) - 0.5D0*r*(u*u + v*v + w*w)), subsub_pfloor)
  
          subsub_hydrobc(1,3,1) = rho_needed
          subsub_hydrobc(1,3,2) = 0.0D0
          subsub_hydrobc(1,3,3) = 0.0D0
          subsub_hydrobc(1,3,4) = subsub_hydrobc(1,3,1)*vff
          subsub_hydrobc(1,3,5) = p/(gamma-1.0D0) * 0.5D0*rho_needed*vff*vff
  
  
  
          r = subsub_hydrobc(2,3,1)
          u = subsub_hydrobc(2,3,2) / r
          v = subsub_hydrobc(2,3,3) / r
          w = subsub_hydrobc(2,3,4) / r
          p = max((gamma - 1.0D0)*(subsub_hydrobc(2,3,5) - 0.5D0*r*(u*u + v*v + w*w)), subsub_pfloor)
  
          subsub_hydrobc(2,3,1) = rho_needed
          subsub_hydrobc(2,3,2) = 0.0D0
          subsub_hydrobc(2,3,3) = 0.0D0
          subsub_hydrobc(2,3,4) = -subsub_hydrobc(2,3,1)*vff
          subsub_hydrobc(2,3,5) = p/(gamma-1.0D0) * 0.5D0*rho_needed*vff*vff
        endif
      else 
        subsub_hydrobc(:,:,:) = subsub_hydrobc_old(:,:,:)
      endif
    endif
    !$omp end single

!
!
!      if(mass_predicted .gt. mass_atthis)then
!        !! Adjust conservative variables to match the mass inflow by free-fall
!
!        frho = (mass_predicted - mass_atthis) / subsub_dt / 6.0D0 / (subsub_boxlen**2)
!        
!        rho_ave = mtot / (subsub_boxlen**3)
!        tff = sqrt(threepi2/8./fourpi/rho_ave/(6.67d-8*scale_d*scale_t**2))
!        vff = subsub_boxlen/tff
!
!        rho_bc = frho / max(abs(vff), subsub_smallr)
!
!        subsub_hydrobc(1,1,2) = subsub_hydrobc(1,1,1)*vff
!        subsub_hydrobc(1,1,3) = 0.0D0
!        subsub_hydrobc(1,1,4) = 0.0D0
!
!        subsub_hydrobc(2,1,2) = -subsub_hydrobc(2,1,1)*vff
!        subsub_hydrobc(2,1,3) = 0.0D0
!        subsub_hydrobc(2,1,4) = 0.0D0
!
!        subsub_hydrobc(1,2,2) = 0.0D0
!        subsub_hydrobc(1,2,3) = subsub_hydrobc(1,2,1)*vff
!        subsub_hydrobc(1,2,4) = 0.0D0
!
!        subsub_hydrobc(2,2,2) = 0.0D0
!        subsub_hydrobc(2,2,3) = -subsub_hydrobc(2,2,1)*vff
!        subsub_hydrobc(2,2,4) = 0.0D0
!
!        subsub_hydrobc(1,3,2) = 0.0D0
!        subsub_hydrobc(1,3,3) = 0.0D0
!        subsub_hydrobc(1,3,4) = subsub_hydrobc(1,3,1)*vff
!
!        subsub_hydrobc(2,3,2) = 0.0D0
!        subsub_hydrobc(2,3,3) = 0.0D0
!        subsub_hydrobc(2,3,4) = -subsub_hydrobc(2,3,1)*vff
!
!
!      else !! cease mass inflow
!        ! zero normal velocity
!        subsub_hydrobc(:,1,2) = 0.0D0
!        subsub_hydrobc(:,2,3) = 0.0D0
!        subsub_hydrobc(:,3,4) = 0.0D0
!      endif
!    endif
!    !$omp end single

!!-----------------------------------------------------------------
!! First Kick (U_n -> U_n+1/2)
!!-----------------------------------------------------------------
    call subsub_kick(subsub_dt/2.0D0)

!$omp single
if(subsub_obj(objind)%sink_id .eq. 3)then
  write(*,*)'         first kick'
  write(*,*)'     niter        = ', subsub_nstep
  write(*,*)'     maxrho       = ', maxval(subsub_hydro(:,1))
  write(*,*)'     minrho       = ', minval(subsub_hydro(:,1))

  write(*,*)'     maxE       = ', maxval(subsub_hydro(:,5))
  write(*,*)'     minE       = ', minval(subsub_hydro(:,5))

  write(*,*)'     maxPx       = ', maxval(subsub_hydro(:,2))
  write(*,*)'     minPx       = ', minval(subsub_hydro(:,2))

  write(*,*)'     maxPy       = ', maxval(subsub_hydro(:,3))
  write(*,*)'     minPy       = ', minval(subsub_hydro(:,3))

  write(*,*)'     maxPz       = ', maxval(subsub_hydro(:,4))
  write(*,*)'     minPz       = ', minval(subsub_hydro(:,4))
endif
!$omp end single

#ifndef WITHOUTMPI
    !$omp single
    subsub_tcheck_cg(4) = subsub_tcheck_cg(4) + MPI_WTIME() - mpinow
    mpinow = MPI_WTIME()
    !$omp end single
#endif


!!-----------------------------------------------------------------
!! Drift (U_n -> U_n+1)
!!-----------------------------------------------------------------
    call subsub_drift(subsub_dt)

!$omp single
if(subsub_obj(objind)%sink_id .eq. 3)then
  write(*,*)'            after drift'
  write(*,*)'     niter        = ', subsub_nstep
  write(*,*)'     maxrho       = ', maxval(subsub_hydro(:,1))
  write(*,*)'     minrho       = ', minval(subsub_hydro(:,1))

  write(*,*)'     maxE       = ', maxval(subsub_hydro(:,5))
  write(*,*)'     minE       = ', minval(subsub_hydro(:,5))

  write(*,*)'     maxPx       = ', maxval(subsub_hydro(:,2))
  write(*,*)'     minPx       = ', minval(subsub_hydro(:,2))

  write(*,*)'     maxPy       = ', maxval(subsub_hydro(:,3))
  write(*,*)'     minPy       = ', minval(subsub_hydro(:,3))

  write(*,*)'     maxPz       = ', maxval(subsub_hydro(:,4))
  write(*,*)'     minPz       = ', minval(subsub_hydro(:,4))
endif
!$omp end single

    !$omp single
#ifndef WITHOUTMPI
    subsub_tcheck_cg(5) = subsub_tcheck_cg(5) + MPI_WTIME() - mpinow
    mpinow = MPI_WTIME()
#endif
    !$omp end single

!!-----------------------------------------------------------------
!! Update Mtot
!!-----------------------------------------------------------------
    !$omp single
    mtot_now = 0.0D0
    !$omp end single

    !$omp do reduction(+:mtot_now)
    do i=1, subsub_nn
      mtot_now = mtot_now + subsub_hydro(i,1)
    enddo
    !$omp end do

    !$omp single
    mtot_now = mtot_now * subsub_dx**ndim
    !mtot = mtot_new
    !$omp end single
    
!!-----------------------------------------------------------------
!! New Phi based on the updated density (rho_n+1 -> phi_n+1)
!!-----------------------------------------------------------------
    call subsub_poisson_cg(mtot_now)
#ifndef WITHOUTMPI
    !$omp single
    subsub_tcheck_cg(1) = subsub_tcheck_cg(1) + MPI_WTIME() - mpinow
    mpinow = MPI_WTIME()
    !$omp end single
#endif

!!-----------------------------------------------------------------
!! Compute g_n+1
!!-----------------------------------------------------------------
    call subsub_poisson_g(fact2)
 
#ifndef WITHOUTMPI
    !$omp single
    subsub_tcheck_cg(2) = subsub_tcheck_cg(2) + MPI_WTIME() - mpinow
    mpinow = MPI_WTIME()
    !$omp end single
#endif

!!-----------------------------------------------------------------
!! Second Kick
!!-----------------------------------------------------------------
    call subsub_kick(subsub_dt/2.0D0)
    
#ifndef WITHOUTMPI 
    !$omp single
    subsub_tcheck_cg(4) = subsub_tcheck_cg(4) + MPI_WTIME() - mpinow
    mpinow = MPI_WTIME()
    !$omp end single
#endif

!!-----------------------------------------------------------------
!! Update old arrays
!!-----------------------------------------------------------------
    !$omp do
    do i=1, subsub_nn
      subsub_obj(objind)%phi(i) = subsub_phi(i)
      subsub_obj(objind)%hydro(i,1) = subsub_hydro(i,1)
      subsub_obj(objind)%hydro(i,2) = subsub_hydro(i,2)
      subsub_obj(objind)%hydro(i,3) = subsub_hydro(i,3)
      subsub_obj(objind)%hydro(i,4) = subsub_hydro(i,4)
      subsub_obj(objind)%hydro(i,5) = subsub_hydro(i,5)
    enddo
    !$omp end do

    !$omp single
    subsub_obj(objind)%mass_tot = mtot_now
    !$omp end single

!$omp single
if(subsub_obj(objind)%sink_id .eq. 3)then
  write(*,*)'              niter        = ', subsub_nstep
  write(*,*)'              maxrho       = ', maxval(subsub_hydro(:,1))
  write(*,*)'              minrho       = ', minval(subsub_hydro(:,1))
endif
!$omp end single
!!-----------------------------------------------------------------
!! Main Loop Control
!!-----------------------------------------------------------------
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
subroutine subsub_poisson_cg(mtot)
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
  use subsub_commons
  use subsub_parameters
  use amr_commons
  use cooling_module, ONLY:twopi
  implicit none

  real(dp) :: mtot

  !! Local variables
  integer :: i, ix, iy, iz, ii, iter


  integer :: xu, xd, yu, yd, zu, zd
  !real(dp) ::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

  !!CG
  real(dp) :: gconst, fact1, fourpiG, pi
  real(dp) :: dx2inv

  !!-----
  !! constants
  !!-----

  !! DEBUG MODE FOR SELF-GRAVITY TEST
  !! )) DEBUGG HYDRO((         <- this is for grep
  if(subsub_dev_hydroonly .eq. 1)then
    return
  endif

  dx2inv = 1.0D0 / (subsub_dx**2)

  !call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  !gconst = 6.67d-8*scale_d*scale_t**2
  pi=twopi/2d0
  gconst = 1d0
  if(cosmo)gconst=3d0/8d0/pi*omega_m*aexp

  fact1 = -gconst * mtot
  fourpiG = gconst * 2.0D0 * twopi



  

  !!!! Fix the boundary potential
  !$omp do private(ii)
  do i=1, subsub_nnface
    ii = subsub_faceind(i)
    subsub_phi(ii) = fact1 / subsub_dd(ii)
  enddo
  !$omp end do

  !!!! Initialize CG arrays
  !$omp single
  subsub_skipcg = .false.
  subsub_rhs2 = 0.0D0
  subsub_rrold = 0.0D0
  !$omp end single

  !!!! Build Initial CG arrays
  !$omp do private(ii)
  do i=1, subsub_nnface
    ii = subsub_faceind(i)
    subsub_cgrhs(ii) = 0.0D0
    subsub_cgLphi(ii) = 0.0D0
    subsub_cgRes(ii) = 0.0D0
    subsub_cgp(ii) = 0.0D0
  enddo
  !$omp end do

  !$omp do collapse(2) reduction(+:subsub_rhs2, subsub_rrold) &
  !$omp & private(ix, iy, iz, ii) &
  !$omp & private(xu, xd, yu, yd, zu, zd)
  do iz=2, subsub_ngrid-1
  do iy=2, subsub_ngrid-1
  do ix=2, subsub_ngrid-1
    ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix

    subsub_cgrhs(ii) = -fourpiG * subsub_hydro(ii,1)

    subsub_rhs2 = subsub_rhs2 + subsub_cgrhs(ii)*subsub_cgrhs(ii)

    xu = ii + 1
    xd = ii - 1
    yu = ii + subsub_ngrid
    yd = ii - subsub_ngrid
    zu = ii + subsub_ngrid2
    zd = ii - subsub_ngrid2

    subsub_cgLphi(ii) = (&
      6.0D0*subsub_phi(ii) - subsub_phi(xu) - subsub_phi(xd) &
      - subsub_phi(yu) - subsub_phi(yd) - subsub_phi(zu) - subsub_phi(zd) ) * dx2inv

    subsub_cgRes(ii) = subsub_cgrhs(ii) - subsub_cgLphi(ii)
    subsub_cgp(ii) = subsub_cgRes(ii)
    subsub_rrold = subsub_rrold + subsub_cgRes(ii)*subsub_cgRes(ii)

  enddo
  enddo
  enddo
  !$omp end do

  !$omp single
  if(sqrt(subsub_rrold / max(subsub_rhs2, subsub_smallr)) .le. subsub_poisson_tolerance) then
    subsub_skipcg = .true.
  endif
  !$omp end single


  !!!! CG Loop
  do iter=1, subsub_poisson_niter
    !$omp barrier
    if(subsub_skipcg) exit

    !$omp single
    subsub_pLp = 0.0D0
    subsub_rrnew = 0.0D0
    subsub_ncheck_cg(2) = subsub_ncheck_cg(2) + 1
    !$omp end single

    !!!!!! Compute Lp & pLp
    !$omp do private(ii)
    do i=1, subsub_nnface
      ii = subsub_faceind(i)
      subsub_cgLp(ii) = 0.0D0
    enddo
    !$omp end do


    !$omp do collapse(2) reduction(+:subsub_pLp) &
    !$omp & private(ix, iy, iz, ii) &
    !$omp & private(xu, xd, yu, yd, zu, zd)
    do iz=2, subsub_ngrid-1
    do iy=2, subsub_ngrid-1
    do ix=2, subsub_ngrid-1
      ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix

      xu = ii + 1
      xd = ii - 1
      yu = ii + subsub_ngrid
      yd = ii - subsub_ngrid
      zu = ii + subsub_ngrid2
      zd = ii - subsub_ngrid2

      subsub_cgLp(ii) = (&
        6.0D0*subsub_cgp(ii) - subsub_cgp(xu) - subsub_cgp(xd) &
        - subsub_cgp(yu) - subsub_cgp(yd) - subsub_cgp(zu) - subsub_cgp(zd) )  * dx2inv
      subsub_pLp = subsub_pLp + subsub_cgp(ii) * subsub_cgLp(ii)
    enddo
    enddo
    enddo
    !$omp end do

    !$omp single
    if(abs(subsub_pLp) .lt. subsub_smallr) then
      subsub_skipcg = .true.
      subsub_alpha = 0.0D0
    else
      subsub_alpha = subsub_rrold / subsub_pLp
    endif
    !$omp end single

    !!!!!! Update
    !$omp do private(ii)
    do i=1, subsub_nnface
      ii = subsub_faceind(i)
      subsub_phi(ii) = fact1 / subsub_dd(ii)
    enddo
    !$omp end do


    !$omp do collapse(2) reduction(+:subsub_rrnew) private(ix, iy, iz, ii)
    do iz=2, subsub_ngrid-1
    do iy=2, subsub_ngrid-1
    do ix=2, subsub_ngrid-1
      ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix

      subsub_cgRes(ii) = subsub_cgRes(ii) - subsub_alpha*subsub_cgLp(ii)
      subsub_rrnew = subsub_rrnew + subsub_cgRes(ii) * subsub_cgRes(ii)
      subsub_phi(ii) = subsub_phi(ii) + subsub_alpha * subsub_cgp(ii)
    enddo
    enddo
    enddo
    !$omp end do


    !$omp single
    subsub_relres = sqrt(subsub_rrnew / max(subsub_rhs2, subsub_smallr))
    if(subsub_relres .lt. subsub_poisson_tolerance) subsub_skipcg = .true.
    subsub_beta = subsub_rrnew / max(subsub_rrold, subsub_smallr)
    !$omp end single

    !! Update2
    !$omp do private(ii)
    do i=1, subsub_nnface
      ii = subsub_faceind(i)
      subsub_cgp(ii) = 0.0D0
    enddo
    !$omp end do

    !$omp do collapse(2) private(ix, iy, iz, ii)
    do iz=2, subsub_ngrid-1
    do iy=2, subsub_ngrid-1
    do ix=2, subsub_ngrid-1
      ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix

      subsub_cgp(ii) = subsub_cgRes(ii) + subsub_beta * subsub_cgp(ii)
    enddo
    enddo
    enddo
    !$omp end do

    !$omp single
    subsub_rrold = subsub_rrnew
    !$omp end single
  enddo !! CG Iteration loop

  !$omp barrier

end subroutine subsub_poisson_cg
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_poisson_g(fact)
  use subsub_commons
  use subsub_parameters

  implicit none
  real(dp) :: fact
  !! Local variables
  integer :: i, ix, iy, iz, ii, xu, xd, yu, yd, zu, zd
  !real(dp) ::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp) :: subsub_twodx

  !! DEBUG MODE FOR SELF-GRAVITY TEST
  !! )) DEBUGG HYDRO((         <- this is for grep
  if(subsub_dev_hydroonly .eq. 1)then
    return
  endif

  !!-----
  !! Add Additional Potential from the sink particle
  !!------
  subsub_twodx = subsub_dx * 2.0D0

  !$omp do
  do i=1, subsub_nn
    subsub_phi(i) = subsub_phi(i) - fact / subsub_dd(i)
  enddo
  !$omp end do


  !$omp do collapse(2) &
  !$omp & private(ix, iy, iz, ii, xu, xd, yu, yd, zu, zd)
  do iz=2, subsub_ngrid-1
  do iy=2, subsub_ngrid-1
  do ix=2, subsub_ngrid-1
    ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix

    xu = ii + 1
    xd = ii - 1
    yu = ii + subsub_ngrid
    yd = ii - subsub_ngrid
    zu = ii + subsub_ngrid2
    zd = ii - subsub_ngrid2

    subsub_fg(ii,1) = - (subsub_phi(xu) - subsub_phi(xd)) / subsub_twodx
    subsub_fg(ii,2) = - (subsub_phi(yu) - subsub_phi(yd)) / subsub_twodx
    subsub_fg(ii,3) = - (subsub_phi(zu) - subsub_phi(zd)) / subsub_twodx

    !! Additional Potential from the sink particle
    subsub_fg(ii,1) = subsub_fg(ii,1) - (fact / subsub_dd(xu) - fact / subsub_dd(xd)) / subsub_twodx
    subsub_fg(ii,2) = subsub_fg(ii,2) - (fact / subsub_dd(yu) - fact / subsub_dd(yd)) / subsub_twodx
    subsub_fg(ii,3) = subsub_fg(ii,3) - (fact / subsub_dd(zu) - fact / subsub_dd(zd)) / subsub_twodx
  enddo
  enddo
  enddo
  !$omp end do

  !$omp do private(ii, ix, iy, iz)
  do i=1, subsub_nnface
    ii = subsub_faceind(i)
    ix = subsub_faceindx(i)
    iy = subsub_faceindy(i)
    iz = subsub_faceindz(i)

    xu = ii + 1
    xd = ii - 1
    yu = ii + subsub_ngrid
    yd = ii - subsub_ngrid
    zu = ii + subsub_ngrid2
    zd = ii - subsub_ngrid2

    if(ix.eq.1) then
      subsub_fg(ii,1) = - (subsub_phi(xu) - subsub_phi(ii)) / subsub_dx

      subsub_fg(ii,1) = subsub_fg(ii,1) - (fact / subsub_dd(xu) - fact / subsub_dd(ii)) / subsub_dx
    else if(ix .eq. subsub_ngrid) then
      subsub_fg(ii,1) = - (subsub_phi(ii) - subsub_phi(xd)) / subsub_dx

      subsub_fg(ii,1) = subsub_fg(ii,1) - (fact / subsub_dd(ii) - fact / subsub_dd(xd)) / subsub_dx
    else
      subsub_fg(ii,1) = - (subsub_phi(xu) - subsub_phi(xd)) / (2.0D0*subsub_dx)

      subsub_fg(ii,1) = subsub_fg(ii,1) - (fact / subsub_dd(xu) - fact / subsub_dd(xd)) / subsub_twodx
    endif

    if(iy.eq.1) then
      subsub_fg(ii,2) = - (subsub_phi(yu) - subsub_phi(ii)) / subsub_dx

      subsub_fg(ii,2) = subsub_fg(ii,2) - (fact / subsub_dd(yu) - fact / subsub_dd(ii)) / subsub_dx
    else if(iy .eq. subsub_ngrid) then
      subsub_fg(ii,2) = - (subsub_phi(ii) - subsub_phi(yd)) / subsub_dx

      subsub_fg(ii,2) = subsub_fg(ii,2) - (fact / subsub_dd(ii) - fact / subsub_dd(yd)) / subsub_dx
    else
      subsub_fg(ii,2) = - (subsub_phi(yu) - subsub_phi(yd)) / (2.0D0*subsub_dx)

      subsub_fg(ii,2) = subsub_fg(ii,2) - (fact / subsub_dd(yu) - fact / subsub_dd(yd)) / subsub_twodx
    endif

    if(iz.eq.1) then
      subsub_fg(ii,3) = - (subsub_phi(zu) - subsub_phi(ii)) / subsub_dx

      subsub_fg(ii,3) = subsub_fg(ii,3) - (fact / subsub_dd(zu) - fact / subsub_dd(ii)) / subsub_dx
    else if(iz .eq. subsub_ngrid) then
      subsub_fg(ii,3) = - (subsub_phi(ii) - subsub_phi(zd)) / subsub_dx

      subsub_fg(ii,3) = subsub_fg(ii,3) - (fact / subsub_dd(ii) - fact / subsub_dd(zd)) / subsub_dx
    else
      subsub_fg(ii,3) = - (subsub_phi(zu) - subsub_phi(zd)) / (2.0D0*subsub_dx)

      subsub_fg(ii,3) = subsub_fg(ii,3) - (fact / subsub_dd(zu) - fact / subsub_dd(zd)) / subsub_twodx
    endif
  enddo
  !$omp end do
end subroutine subsub_poisson_g
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_kick(dt)
  use subsub_commons
  use subsub_parameters
  use hydro_parameters, ONLY: gamma
  implicit none

  real(dp) :: dt

  !! Local variables
  integer :: i
  real(dp) :: vx, vy, vz, vxn, vyn, vzn, rho, ekin_floor

  !! DEBUG MODE FOR SELF-GRAVITY TEST
  !! )) DEBUGG HYDRO((         <- this is for grep
  if(subsub_dev_hydroonly .eq. 1) then
    return
  endif

  !$omp do private(vx, vy, vz, vxn, vyn, vzn, rho, ekin_floor)
  do i=1, subsub_nn
    rho = max(subsub_hydro(i,1), subsub_dfloor)

    vx = subsub_hydro(i,2) / rho
    vy = subsub_hydro(i,3) / rho
    vz = subsub_hydro(i,4) / rho

    vxn = vx + subsub_fg(i,1) * dt
    vyn = vy + subsub_fg(i,2) * dt
    vzn = vz + subsub_fg(i,3) * dt

    !! Momentum kick
    subsub_hydro(i,2) = rho * vxn
    subsub_hydro(i,3) = rho * vyn
    subsub_hydro(i,4) = rho * vzn

    !! Energy : rho * (v dot g)
    subsub_hydro(i,5) = subsub_hydro(i,5) + 0.5D0 * rho * &
      ( (vx + vxn)*subsub_fg(i,1) + (vy + vyn)*subsub_fg(i,2) + (vz + vzn)*subsub_fg(i,3) ) * dt

    ekin_floor = 0.5D0 * rho * (vxn*vxn + vyn*vyn + vzn*vzn)
    subsub_hydro(i,5) = max(subsub_hydro(i,5), ekin_floor + subsub_pfloor / (gamma - 1.0D0))
  enddo
  !$omp end do
end subroutine subsub_kick
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_drift(dt)
  use amr_commons
  use subsub_commons
  use subsub_parameters
  use hydro_parameters, ONLY:gamma
  use, intrinsic :: ieee_arithmetic
  implicit none

  real(dp) :: dt

  !! Local variables
  integer :: i, ivar
  integer :: ii, ix, iy, iz
  real(dp) :: ekin
  real(dp), dimension(1:subsub_nhydro) :: Utmp
  !real(dp), dimension(1:2, 1:ndim, 1:subsub_nhydro, 1:5) :: hbc


  !! copy to old arrays & save primitive variables
  !$omp do collapse(2) private(ekin, ii, ix, iy, iz, ivar, Utmp)
  do iz=1, subsub_ngrid
  do iy=1, subsub_ngrid
  do ix=1, subsub_ngrid
    ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix

    do ivar=1, subsub_nhydro
      Utmp(ivar) = subsub_hydro(ii,ivar)
    enddo
    call subsub_enforce_floors(Utmp)
    do ivar=1, subsub_nhydro
      subsub_hydro(ii,ivar) = Utmp(ivar)
    enddo

    subsub_hdummy(ix,iy,iz,1,1) = subsub_hydro(ii,1)
    subsub_hdummy(ix,iy,iz,2,1) = subsub_hydro(ii,2)
    subsub_hdummy(ix,iy,iz,3,1) = subsub_hydro(ii,3)
    subsub_hdummy(ix,iy,iz,4,1) = subsub_hydro(ii,4)
    subsub_hdummy(ix,iy,iz,5,1) = subsub_hydro(ii,5)

    subsub_hdummy(ix,iy,iz,1,2) = subsub_hydro(ii,1)
    subsub_hdummy(ix,iy,iz,2,2) = subsub_hydro(ii,2) / subsub_hydro(ii,1)
    subsub_hdummy(ix,iy,iz,3,2) = subsub_hydro(ii,3) / subsub_hydro(ii,1)
    subsub_hdummy(ix,iy,iz,4,2) = subsub_hydro(ii,4) / subsub_hydro(ii,1)

    ekin = 0.5D0 * subsub_hydro(ii,1) * (subsub_hdummy(ix,iy,iz,2,2)**2 + subsub_hdummy(ix,iy,iz,3,2)**2 + subsub_hdummy(ix,iy,iz,4,2)**2)
    subsub_hdummy(ix,iy,iz,5,2) = max((gamma - 1.0D0) * (subsub_hdummy(ix,iy,iz,5,1) - ekin), subsub_pfloor)

    !! DEBUG MODE FOR SELF-GRAVITY TEST
    !! )) DEBUGG GRAV((         <- this is for grep
    if(subsub_dev_gravonly .eq. 1)then
      subsub_hdummy(ix,iy,iz,5,2) = 0.0D0
    endif

  enddo
  enddo
  enddo
  !$omp end do


  select case(subsub_RiemannType)
  case(1)
!!-----------------------------------------------------------------
!! Rusanov Solver
!!-----------------------------------------------------------------

  call subsub_hydroRiemann_Rusanov(dt)

  end select


  !$omp do private(Utmp, ivar)
  do i=1, subsub_nn

    do ivar=1, subsub_nhydro
      Utmp(ivar) = subsub_hydro(i,ivar)
    enddo
    call subsub_enforce_floors(Utmp)
    do ivar=1, subsub_nhydro
      subsub_hydro(i,ivar) = Utmp(ivar)
    enddo
  enddo
  !$omp end do

end subroutine subsub_drift
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_hydroRiemann_Rusanov(dt)!, hbc)
  use amr_commons
  use subsub_commons
  use subsub_parameters
  use hydro_parameters, ONLY:gamma
  implicit none

  real(dp) :: dt
  real(dp), dimension(1:2, 1:ndim, 1:subsub_nhydro, 1:5) :: hbc
  !! Local variables
  integer :: ivar, i, j
  integer :: ix, iy, iz
  integer :: ii, xu, xd, yu, yd, zu, zd
  integer :: ixu, ixd, iyu, iyd, izu, izd
  real(dp) :: amaxL, amaxR, lam, ekin, bc_vv, bc_cs, drho, dfac, pnew
  
  real(dp), dimension(1:subsub_nhydro) :: FL, FR, Utmp, bc_flux, bc_cons
  real(dp), dimension(1:2, 1:ndim) :: csarr_bc
  logical :: isodd, okay


  !1 conservative old
  !2 primitive old
  !3 Fx
  !4 Fy
  !5 Fz


  lam = dt / subsub_dx
  isodd = mod(subsub_ngrid,2) .eq. 1

  !----- Set boundary
  !! Conservative
  hbc(1,1,:,1) = subsub_hydrobc(1,1,:)
  hbc(2,1,:,1) = subsub_hydrobc(2,1,:)
  hbc(1,2,:,1) = subsub_hydrobc(1,2,:)
  hbc(2,2,:,1) = subsub_hydrobc(2,2,:)
  hbc(1,3,:,1) = subsub_hydrobc(1,3,:)
  hbc(2,3,:,1) = subsub_hydrobc(2,3,:)

  !! Primitive
  !(X-)
  hbc(1,1,1,2) = hbc(1,1,1,1)
  hbc(1,1,2,2) = hbc(1,1,2,1) / hbc(1,1,1,1)
  hbc(1,1,3,2) = hbc(1,1,3,1) / hbc(1,1,1,1)
  hbc(1,1,4,2) = hbc(1,1,4,1) / hbc(1,1,1,1)

  ekin = 0.5D0 * hbc(1,1,1,1) * (hbc(1,1,2,2)**2 + hbc(1,1,3,2)**2 + hbc(1,1,4,2)**2)
  hbc(1,1,5,2) = max((gamma - 1.0D0) * (hbc(1,1,5,1) - ekin),subsub_pfloor)

  !(X+)  
  hbc(2,1,1,2) = hbc(2,1,1,1)
  hbc(2,1,2,2) = hbc(2,1,2,1) / hbc(2,1,1,1)
  hbc(2,1,3,2) = hbc(2,1,3,1) / hbc(2,1,1,1)
  hbc(2,1,4,2) = hbc(2,1,4,1) / hbc(2,1,1,1)

  ekin = 0.5D0 * hbc(2,1,1,1) * (hbc(2,1,2,2)**2 + hbc(2,1,3,2)**2 + hbc(2,1,4,2)**2)
  hbc(2,1,5,2) = max((gamma - 1.0D0) * (hbc(2,1,5,1) - ekin),subsub_pfloor)

  !(Y-) 
  hbc(1,2,1,2) = hbc(1,2,1,1)
  hbc(1,2,2,2) = hbc(1,2,2,1) / hbc(1,2,1,1)
  hbc(1,2,3,2) = hbc(1,2,3,1) / hbc(1,2,1,1)
  hbc(1,2,4,2) = hbc(1,2,4,1) / hbc(1,2,1,1)

  ekin = 0.5D0 * hbc(1,2,1,1) * (hbc(1,2,2,2)**2 + hbc(1,2,3,2)**2 + hbc(1,2,4,2)**2)
  hbc(1,2,5,2) = max((gamma - 1.0D0) * (hbc(1,2,5,1) - ekin),subsub_pfloor)

  !(Y+) 
  hbc(2,2,1,2) = hbc(2,2,1,1)
  hbc(2,2,2,2) = hbc(2,2,2,1) / hbc(2,2,1,1)
  hbc(2,2,3,2) = hbc(2,2,3,1) / hbc(2,2,1,1)
  hbc(2,2,4,2) = hbc(2,2,4,1) / hbc(2,2,1,1)

  ekin = 0.5D0 * hbc(2,2,1,1) * (hbc(2,2,2,2)**2 + hbc(2,2,3,2)**2 + hbc(2,2,4,2)**2)
  hbc(2,2,5,2) = max((gamma - 1.0D0) * (hbc(2,2,5,1) - ekin),subsub_pfloor)

  !(Z-) 
  hbc(1,3,1,2) = hbc(1,3,1,1)
  hbc(1,3,2,2) = hbc(1,3,2,1) / hbc(1,3,1,1)
  hbc(1,3,3,2) = hbc(1,3,3,1) / hbc(1,3,1,1)
  hbc(1,3,4,2) = hbc(1,3,4,1) / hbc(1,3,1,1)

  ekin = 0.5D0 * hbc(1,3,1,1) * (hbc(1,3,2,2)**2 + hbc(1,3,3,2)**2 + hbc(1,3,4,2)**2)
  hbc(1,3,5,2) = max((gamma - 1.0D0) * (hbc(1,3,5,1) - ekin),subsub_pfloor)

  !(Z+) 
  hbc(2,3,1,2) = hbc(2,3,1,1)
  hbc(2,3,2,2) = hbc(2,3,2,1) / hbc(2,3,1,1)
  hbc(2,3,3,2) = hbc(2,3,3,1) / hbc(2,3,1,1)
  hbc(2,3,4,2) = hbc(2,3,4,1) / hbc(2,3,1,1)

  ekin = 0.5D0 * hbc(2,3,1,1) * (hbc(2,3,2,2)**2 + hbc(2,3,3,2)**2 + hbc(2,3,4,2)**2)
  hbc(2,3,5,2) = max((gamma - 1.0D0) * (hbc(2,3,5,1) - ekin),subsub_pfloor)

  !! DEBUG MODE FOR SELF-GRAVITY TEST
  !! )) DEBUGG GRAV((         <- this is for grep
  if(subsub_dev_gravonly .eq. 1)then
    do i=1,2
      do j=1,3
        hbc(i,j,5,1) = 0.0D0
        hbc(i,j,5,2) = 0.0D0
      enddo
    enddo
  endif


  !! X_FX-
  hbc(1,1,1,3) = hbc(1,1,2,1)
  hbc(1,1,2,3) = hbc(1,1,2,1)*hbc(1,1,2,2) + hbc(1,1,5,2)
  hbc(1,1,3,3) = hbc(1,1,3,1)*hbc(1,1,2,2)
  hbc(1,1,4,3) = hbc(1,1,4,1)*hbc(1,1,2,2)
  hbc(1,1,5,3) = (hbc(1,1,5,1) + hbc(1,1,5,2)) * hbc(1,1,2,2)

  !! X_FX+
  hbc(2,1,1,3) = hbc(2,1,2,1)
  hbc(2,1,2,3) = hbc(2,1,2,1)*hbc(2,1,2,2) + hbc(2,1,5,2)
  hbc(2,1,3,3) = hbc(2,1,3,1)*hbc(2,1,2,2)
  hbc(2,1,4,3) = hbc(2,1,4,1)*hbc(2,1,2,2)
  hbc(2,1,5,3) = (hbc(2,1,5,1) + hbc(2,1,5,2)) * hbc(2,1,2,2)

  !! Y_FY-
  hbc(1,2,1,4) = hbc(1,2,3,1)
  hbc(1,2,2,4) = hbc(1,2,2,1)*hbc(1,2,3,2)
  hbc(1,2,3,4) = hbc(1,2,3,1)*hbc(1,2,3,2) + hbc(1,2,5,2)
  hbc(1,2,4,4) = hbc(1,2,4,1)*hbc(1,2,3,2)
  hbc(1,2,5,4) = (hbc(1,2,5,1) + hbc(1,2,5,2)) * hbc(1,2,3,2)


  !! Y_FY+
  hbc(2,2,1,4) = hbc(2,2,3,1)
  hbc(2,2,2,4) = hbc(2,2,2,1)*hbc(2,2,3,2)
  hbc(2,2,3,4) = hbc(2,2,3,1)*hbc(2,2,3,2) + hbc(2,2,5,2)
  hbc(2,2,4,4) = hbc(2,2,4,1)*hbc(2,2,3,2)
  hbc(2,2,5,4) = (hbc(2,2,5,1) + hbc(2,2,5,2)) * hbc(2,2,3,2)


  !! Z_FZ-
  hbc(1,3,1,5) = hbc(1,3,4,1)
  hbc(1,3,2,5) = hbc(1,3,2,1)*hbc(1,3,4,2)
  hbc(1,3,3,5) = hbc(1,3,3,1)*hbc(1,3,4,2)
  hbc(1,3,4,5) = hbc(1,3,4,1)*hbc(1,3,4,2) + hbc(1,3,5,2)
  hbc(1,3,5,5) = (hbc(1,3,5,1) + hbc(1,3,5,2)) * hbc(1,3,4,2)

  !! Z_FZ+
  hbc(2,3,1,5) = hbc(2,3,4,1)
  hbc(2,3,2,5) = hbc(2,3,2,1)*hbc(2,3,4,2)
  hbc(2,3,3,5) = hbc(2,3,3,1)*hbc(2,3,4,2)
  hbc(2,3,4,5) = hbc(2,3,4,1)*hbc(2,3,4,2) + hbc(2,3,5,2)
  hbc(2,3,5,5) = (hbc(2,3,5,1) + hbc(2,3,5,2)) * hbc(2,3,4,2)


  

  !!----- Compute sound speed first
  !$omp do collapse(2) private(ix, iy, iz)
  do iz=1, subsub_ngrid
  do iy=1, subsub_ngrid
  do ix=1, subsub_ngrid
    if(subsub_hdummy(ix,iy,iz,1,2) .gt. subsub_dfloor) then
      subsub_csarr(ix,iy,iz) = sqrt(gamma * subsub_hdummy(ix,iy,iz,5,2) / subsub_hdummy(ix,iy,iz,1,2))
    else
      subsub_csarr(ix,iy,iz) = sqrt(gamma * subsub_pfloor / subsub_dfloor)
    endif
  enddo
  enddo
  enddo
  !$omp end do

  !! DEBUG MODE FOR SELF-GRAVITY TEST
  !! )) DEBUGG GRAV ((         <- this is for grep
  if(subsub_dev_gravonly .eq. 1)then
    !$omp do collapse(2) private(ix, iy, iz)
    do iz=1, subsub_ngrid
    do iy=1, subsub_ngrid
    do ix=1, subsub_ngrid
      subsub_csarr(ix,iy,iz) = 0.0D0
    enddo
    enddo
    enddo
  endif



  do i=1, 2
    do j=1, ndim
      if(hbc(i,j,1,2) .gt. subsub_dfloor) then
        csarr_bc(i,j) = sqrt(gamma * hbc(i,j,5,2) / hbc(i,j,1,2))
      else
        csarr_bc(i,j) = sqrt(gamma * subsub_pfloor / subsub_dfloor)
      endif
    enddo
  enddo

  !! DEBUG MODE FOR SELF-GRAVITY TEST
  !! )) DEBUGG GRAV((         <- this is for grep
  if(subsub_dev_gravonly .eq. 1)then
    do i=1,2
    do j=1, ndim
      csarr_bc(i,j) = 0.0D0
    enddo
    enddo
  endif

  !!----- Compute Flux
  !$omp do collapse(2) private(ix, iy, iz)
  do iz=1, subsub_ngrid
  do iy=1, subsub_ngrid
  do ix=1, subsub_ngrid
    subsub_hdummy(ix,iy,iz,1,3) = subsub_hdummy(ix,iy,iz,2,1)
    subsub_hdummy(ix,iy,iz,2,3) = subsub_hdummy(ix,iy,iz,2,1)*subsub_hdummy(ix,iy,iz,2,2) + subsub_hdummy(ix,iy,iz,5,2)
    subsub_hdummy(ix,iy,iz,3,3) = subsub_hdummy(ix,iy,iz,3,1)*subsub_hdummy(ix,iy,iz,2,2)
    subsub_hdummy(ix,iy,iz,4,3) = subsub_hdummy(ix,iy,iz,4,1)*subsub_hdummy(ix,iy,iz,2,2)
    subsub_hdummy(ix,iy,iz,5,3) = (subsub_hdummy(ix,iy,iz,5,1) + subsub_hdummy(ix,iy,iz,5,2)) * subsub_hdummy(ix,iy,iz,2,2)

    subsub_hdummy(ix,iy,iz,1,4) = subsub_hdummy(ix,iy,iz,3,1)
    subsub_hdummy(ix,iy,iz,2,4) = subsub_hdummy(ix,iy,iz,2,1)*subsub_hdummy(ix,iy,iz,3,2)
    subsub_hdummy(ix,iy,iz,3,4) = subsub_hdummy(ix,iy,iz,3,1)*subsub_hdummy(ix,iy,iz,3,2) + subsub_hdummy(ix,iy,iz,5,2)
    subsub_hdummy(ix,iy,iz,4,4) = subsub_hdummy(ix,iy,iz,4,1)*subsub_hdummy(ix,iy,iz,3,2)
    subsub_hdummy(ix,iy,iz,5,4) = (subsub_hdummy(ix,iy,iz,5,1) + subsub_hdummy(ix,iy,iz,5,2)) * subsub_hdummy(ix,iy,iz,3,2)

    subsub_hdummy(ix,iy,iz,1,5) = subsub_hdummy(ix,iy,iz,4,1)
    subsub_hdummy(ix,iy,iz,2,5) = subsub_hdummy(ix,iy,iz,2,1)*subsub_hdummy(ix,iy,iz,4,2)
    subsub_hdummy(ix,iy,iz,3,5) = subsub_hdummy(ix,iy,iz,3,1)*subsub_hdummy(ix,iy,iz,4,2)
    subsub_hdummy(ix,iy,iz,4,5) = subsub_hdummy(ix,iy,iz,4,1)*subsub_hdummy(ix,iy,iz,4,2) + subsub_hdummy(ix,iy,iz,5,2)
    subsub_hdummy(ix,iy,iz,5,5) = (subsub_hdummy(ix,iy,iz,5,1) + subsub_hdummy(ix,iy,iz,5,2)) * subsub_hdummy(ix,iy,iz,4,2)
  enddo
  enddo
  enddo
  !$omp end do

  !! Solver
!!----- RHEE -----
!! Simpler version
!!----------------
!! Solver

!!------------------------------------------------------
!! X-axis
!!------------------------------------------------------

  !(left boundary)
    !$omp do collapse(2) private(ii, ivar, ixu, ixd, iy, iz, Utmp, okay, pnew, dfac, bc_vv, bc_cs, bc_flux, bc_cons, amaxL, amaxR, FL, FR)
    do iz=1, subsub_ngrid
    do iy=1, subsub_ngrid
  
      ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + 1
      ixu = 1+1
      ixd = 1-1

      bc_vv = hbc(1,1,2,2)
      bc_cs  = csarr_bc(1,1)
      bc_flux(:) = hbc(1,1,:,3)
      bc_cons(:) = hbc(1,1,:,1)

      !! DEBUG MODE FOR SELF-GRAVITY TEST
      !! )) DEBUGG HYDRO ((         <- this is for grep
      if(subsub_dev_hydroonly .eq. 1)then
        ixd = subsub_ngrid
        bc_vv = subsub_hdummy(ixd,iy,iz,2,2)
        bc_cs = subsub_csarr(ixd,iy,iz)
        bc_flux(:) = subsub_hdummy(ixd,iy,iz,:,3)
        bc_cons(:) = subsub_hdummy(ixd,iy,iz,:,1)
      endif


      amaxL = max(abs(subsub_hdummy(1,iy,iz,2,2)) + subsub_csarr(1,iy,iz), abs(bc_vv) + bc_cs)
      amaxR = max(abs(subsub_hdummy(1,iy,iz,2,2)) + subsub_csarr(1,iy,iz), abs(subsub_hdummy(ixu,iy,iz,2,2)) + subsub_csarr(ixu,iy,iz))

      do ivar=1,subsub_nhydro
        FL(ivar) = 0.5D0*(bc_flux(ivar) + subsub_hdummy(1,iy,iz,ivar,3)) - &
          0.5D0*amaxL*(subsub_hdummy(1,iy,iz,ivar,1) - bc_cons(ivar))

        FR(ivar) = 0.5D0*(subsub_hdummy(1,iy,iz,ivar,3) + subsub_hdummy(ixu,iy,iz,ivar,3)) - &
          0.5D0*amaxR*(subsub_hdummy(ixu,iy,iz,ivar,1) - subsub_hdummy(1,iy,iz,ivar,1))
      enddo



      !! Before update, check the negative density & pressure
      do ivar=1,subsub_nhydro
        Utmp(ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
      enddo

      okay = .true.
      if(Utmp(1) .lt. subsub_dfloor) okay=.false.
      pnew = (Utmp(5) - 0.5D0/max(Utmp(1),subsub_dfloor)*(Utmp(2)**2 + Utmp(3)**2 + Utmp(4)**2))*(gamma - 1.0D0)
      if(pnew .lt. subsub_pfloor) okay=.false.

      if(.not. okay)then
        !!----- RHEE -----
        !! Not physical, but limit negative pressure and density
        dfac = 0.5D0

        do
          do ivar=1,subsub_nhydro
            Utmp(ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))*dfac
          enddo

          okay = .true.
          if(Utmp(1) .lt. subsub_dfloor) okay=.false.
          pnew = (Utmp(5) - 0.5D0/max(Utmp(1),subsub_dfloor)*(Utmp(2)**2 + Utmp(3)**2 + Utmp(4)**2))*(gamma - 1.0D0)
          if(pnew .lt. subsub_pfloor) okay=.false.

          if(.not.okay) dfac = dfac * 0.5D0

          if(okay) exit
          if(dfac .lt. 1.0D-12)then !! fail to update this cell
            do ivar=1,subsub_nhydro
              Utmp(ivar) = subsub_hydro(ii,ivar)
            enddo
            exit
          endif
        enddo
      endif

      !! Update
      do ivar=1,subsub_nhydro
        subsub_hydro(ii,ivar) = Utmp(ivar)
      enddo

    enddo
    enddo
    !$omp end do

  !(interior)
  do ix=2, subsub_ngrid-1
    !$omp do collapse(2) private(ii, ivar, ixu, ixd, iy, iz, Utmp, okay, pnew, dfac, amaxL, amaxR, FL, FR)
    do iy=1, subsub_ngrid
    do iz=1, subsub_ngrid
      ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
      ixu = ix+1
      ixd = ix-1

      amaxL = max(abs(subsub_hdummy(ix,iy,iz,2,2)) + subsub_csarr(ix,iy,iz), abs(subsub_hdummy(ixd,iy,iz,2,2)) + subsub_csarr(ixd,iy,iz))
      amaxR = max(abs(subsub_hdummy(ix,iy,iz,2,2)) + subsub_csarr(ix,iy,iz), abs(subsub_hdummy(ixu,iy,iz,2,2)) + subsub_csarr(ixu,iy,iz))

      do ivar=1,subsub_nhydro
        FL(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,3) + subsub_hdummy(ixd,iy,iz,ivar,3)) - &
          0.5D0*amaxL*(subsub_hdummy(ix,iy,iz,ivar,1) - subsub_hdummy(ixd,iy,iz,ivar,1))
        FR(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,3) + subsub_hdummy(ixu,iy,iz,ivar,3)) - &
          0.5D0*amaxR*(subsub_hdummy(ixu,iy,iz,ivar,1) - subsub_hdummy(ix,iy,iz,ivar,1))
      enddo

      !! Before update, check the negative density & pressure
      do ivar=1,subsub_nhydro
        Utmp(ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
      enddo

      okay = .true.
      if(Utmp(1) .lt. subsub_dfloor) okay=.false.
      pnew = (Utmp(5) - 0.5D0/max(Utmp(1),subsub_dfloor)*(Utmp(2)**2 + Utmp(3)**2 + Utmp(4)**2))*(gamma - 1.0D0)
      if(pnew .lt. subsub_pfloor) okay=.false.

      if(.not. okay)then
        !!----- RHEE -----
        !! Not physical, but limit negative pressure and density
        dfac = 0.5D0

        do
          do ivar=1,subsub_nhydro
            Utmp(ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))*dfac
          enddo

          okay = .true.
          if(Utmp(1) .lt. subsub_dfloor) okay=.false.
          pnew = (Utmp(5) - 0.5D0/max(Utmp(1),subsub_dfloor)*(Utmp(2)**2 + Utmp(3)**2 + Utmp(4)**2))*(gamma - 1.0D0)
          if(pnew .lt. subsub_pfloor) okay=.false.

          if(.not.okay) dfac = dfac * 0.5D0

          if(okay) exit
          if(dfac .lt. 1.0D-12)then !! fail to update this cell
            do ivar=1,subsub_nhydro
              Utmp(ivar) = subsub_hydro(ii,ivar)
            enddo
            exit
          endif
        enddo
      endif

      !! Update
      do ivar=1,subsub_nhydro
        subsub_hydro(ii,ivar) = Utmp(ivar)
      enddo


    enddo
    enddo
    !$omp end do
  enddo

  !(right boundary)
    !$omp do collapse(2) private(ii, ivar, ixu, ixd, iy, iz, Utmp, okay, pnew, dfac, bc_vv, bc_cs, bc_flux, bc_cons, amaxL, amaxR, FL, FR)
    do iz=1, subsub_ngrid
    do iy=1, subsub_ngrid
  
      ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + subsub_ngrid
      ixu = subsub_ngrid+1
      ixd = subsub_ngrid-1

      bc_vv = hbc(2,1,2,2)
      bc_cs  = csarr_bc(2,1)
      bc_flux(:) = hbc(2,1,:,3)
      bc_cons(:) = hbc(2,1,:,1)

      !! DEBUG MODE FOR SELF-GRAVITY TEST
      !! )) DEBUGG HYDRO ((         <- this is for grep
      if(subsub_dev_hydroonly .eq. 1)then
        ixu = 1
        bc_vv = subsub_hdummy(ixu,iy,iz,2,2)
        bc_cs = subsub_csarr(ixu,iy,iz)
        bc_flux(:) = subsub_hdummy(ixu,iy,iz,:,3)
        bc_cons(:) = subsub_hdummy(ixu,iy,iz,:,1)
      endif

      amaxL = max(abs(subsub_hdummy(subsub_ngrid,iy,iz,2,2)) + subsub_csarr(subsub_ngrid,iy,iz), abs(subsub_hdummy(ixd,iy,iz,2,2)) + subsub_csarr(ixd,iy,iz))
      amaxR = max(abs(subsub_hdummy(subsub_ngrid,iy,iz,2,2)) + subsub_csarr(subsub_ngrid,iy,iz), abs(bc_vv) + bc_cs)

      do ivar=1,subsub_nhydro
        FL(ivar) = 0.5D0*(subsub_hdummy(subsub_ngrid,iy,iz,ivar,3) + subsub_hdummy(ixd,iy,iz,ivar,3)) - &
          0.5D0*amaxL*(subsub_hdummy(subsub_ngrid,iy,iz,ivar,1) - subsub_hdummy(ixd,iy,iz,ivar,1))
        FR(ivar) = 0.5D0*(bc_flux(ivar) + subsub_hdummy(subsub_ngrid,iy,iz,ivar,3)) - &
          0.5D0*amaxR*(bc_cons(ivar) - subsub_hdummy(subsub_ngrid,iy,iz,ivar,1))
      enddo


      !! Before update, check the negative density & pressure
      do ivar=1,subsub_nhydro
        Utmp(ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
      enddo

      okay = .true.
      if(Utmp(1) .lt. subsub_dfloor) okay=.false.
      pnew = (Utmp(5) - 0.5D0/max(Utmp(1),subsub_dfloor)*(Utmp(2)**2 + Utmp(3)**2 + Utmp(4)**2))*(gamma - 1.0D0)
      if(pnew .lt. subsub_pfloor) okay=.false.

      if(.not. okay)then
        !!----- RHEE -----
        !! Not physical, but limit negative pressure and density
        dfac = 0.5D0

        do
          do ivar=1,subsub_nhydro
            Utmp(ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))*dfac
          enddo

          okay = .true.
          if(Utmp(1) .lt. subsub_dfloor) okay=.false.
          pnew = (Utmp(5) - 0.5D0/max(Utmp(1),subsub_dfloor)*(Utmp(2)**2 + Utmp(3)**2 + Utmp(4)**2))*(gamma - 1.0D0)
          if(pnew .lt. subsub_pfloor) okay=.false.

          if(.not.okay) dfac = dfac * 0.5D0

          if(okay) exit
          if(dfac .lt. 1.0D-12)then !! fail to update this cell
            do ivar=1,subsub_nhydro
              Utmp(ivar) = subsub_hydro(ii,ivar)
            enddo
            exit
          endif
        enddo
      endif

      !! Update
      do ivar=1,subsub_nhydro
        subsub_hydro(ii,ivar) = Utmp(ivar)
      enddo

    enddo
    enddo
    !$omp end do


!!------------------------------------------------------
!! Y-axis
!!------------------------------------------------------
  !(left boundary)
    !$omp do collapse(2) private(ii, ivar, iyu, iyd, ix, iz, Utmp, okay, pnew, dfac, bc_vv, bc_cs, bc_flux, bc_cons, amaxL, amaxR, FL, FR)
    do iz=1, subsub_ngrid
    do ix=1, subsub_ngrid
  
      ii = (iz-1)*subsub_ngrid2 + ix
      iyu = 1+1
      iyd = 1-1

      bc_vv = hbc(1,2,3,2)
      bc_cs  = csarr_bc(1,2)
      bc_flux(:) = hbc(1,2,:,4)
      bc_cons(:) = hbc(1,2,:,1)

      !! DEBUG MODE FOR SELF-GRAVITY TEST
      !! )) DEBUGG HYDRO ((         <- this is for grep
      if(subsub_dev_hydroonly .eq. 1)then
        iyd = subsub_ngrid
        bc_vv = subsub_hdummy(ix,iyd,iz,3,2)
        bc_cs = subsub_csarr(ix,iyd,iz)
        bc_flux(:) = subsub_hdummy(ix,iyd,iz,:,4)
        bc_cons(:) = subsub_hdummy(ix,iyd,iz,:,1)
      endif


      amaxL = max(abs(subsub_hdummy(ix,1,iz,3,2)) + subsub_csarr(ix,1,iz), abs(bc_vv) + bc_cs)
      amaxR = max(abs(subsub_hdummy(ix,1,iz,3,2)) + subsub_csarr(ix,1,iz), abs(subsub_hdummy(ix,iyu,iz,3,2)) + subsub_csarr(ix,iyu,iz))

      do ivar=1,subsub_nhydro
        FL(ivar) = 0.5D0*(bc_flux(ivar) + subsub_hdummy(ix,1,iz,ivar,4)) - &
          0.5D0*amaxL*(subsub_hdummy(ix,1,iz,ivar,1) - bc_cons(ivar))

        FR(ivar) = 0.5D0*(subsub_hdummy(ix,1,iz,ivar,4) + subsub_hdummy(ix,iyu,iz,ivar,4)) - &
          0.5D0*amaxR*(subsub_hdummy(ix,iyu,iz,ivar,1) - subsub_hdummy(ix,1,iz,ivar,1))
      enddo


      !! Avoid Negative Density

      !! Before update, check the negative density & pressure
      do ivar=1,subsub_nhydro
        Utmp(ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
      enddo

      okay = .true.
      if(Utmp(1) .lt. subsub_dfloor) okay=.false.
      pnew = (Utmp(5) - 0.5D0/max(Utmp(1),subsub_dfloor)*(Utmp(2)**2 + Utmp(3)**2 + Utmp(4)**2))*(gamma - 1.0D0)
      if(pnew .lt. subsub_pfloor) okay=.false.

      if(.not. okay)then
        !!----- RHEE -----
        !! Not physical, but limit negative pressure and density
        dfac = 0.5D0

        do
          do ivar=1,subsub_nhydro
            Utmp(ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))*dfac
          enddo

          okay = .true.
          if(Utmp(1) .lt. subsub_dfloor) okay=.false.
          pnew = (Utmp(5) - 0.5D0/max(Utmp(1),subsub_dfloor)*(Utmp(2)**2 + Utmp(3)**2 + Utmp(4)**2))*(gamma - 1.0D0)
          if(pnew .lt. subsub_pfloor) okay=.false.

          if(.not.okay) dfac = dfac * 0.5D0

          if(okay) exit
          if(dfac .lt. 1.0D-12)then !! fail to update this cell
            do ivar=1,subsub_nhydro
              Utmp(ivar) = subsub_hydro(ii,ivar)
            enddo
            exit
          endif
        enddo
      endif

      !! Update
      do ivar=1,subsub_nhydro
        subsub_hydro(ii,ivar) = Utmp(ivar)
      enddo

    enddo
    enddo
    !$omp end do


  !(interior)
  do iy=2, subsub_ngrid-1
    !$omp do collapse(2) private(ii, ivar, iyu, iyd, ix, iz, Utmp, okay, pnew, dfac, amaxL, amaxR, FL, FR)
    do iz=1, subsub_ngrid
    do ix=1, subsub_ngrid
      ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
      iyu = iy+1
      iyd = iy-1

      amaxL = max(abs(subsub_hdummy(ix,iy,iz,3,2)) + subsub_csarr(ix,iy,iz), abs(subsub_hdummy(ix,iyd,iz,3,2)) + subsub_csarr(ix,iyd,iz))
      amaxR = max(abs(subsub_hdummy(ix,iy,iz,3,2)) + subsub_csarr(ix,iy,iz), abs(subsub_hdummy(ix,iyu,iz,3,2)) + subsub_csarr(ix,iyu,iz))

      do ivar=1,subsub_nhydro
        FL(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,4) + subsub_hdummy(ix,iyd,iz,ivar,4)) - &
          0.5D0*amaxL*(subsub_hdummy(ix,iy,iz,ivar,1) - subsub_hdummy(ix,iyd,iz,ivar,1))
        FR(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,4) + subsub_hdummy(ix,iyu,iz,ivar,4)) - &
          0.5D0*amaxR*(subsub_hdummy(ix,iyu,iz,ivar,1) - subsub_hdummy(ix,iy,iz,ivar,1))
      enddo

      !! Before update, check the negative density & pressure
      do ivar=1,subsub_nhydro
        Utmp(ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
      enddo

      okay = .true.
      if(Utmp(1) .lt. subsub_dfloor) okay=.false.
      pnew = (Utmp(5) - 0.5D0/max(Utmp(1),subsub_dfloor)*(Utmp(2)**2 + Utmp(3)**2 + Utmp(4)**2))*(gamma - 1.0D0)
      if(pnew .lt. subsub_pfloor) okay=.false.

      if(.not. okay)then
        !!----- RHEE -----
        !! Not physical, but limit negative pressure and density
        dfac = 0.5D0

        do
          do ivar=1,subsub_nhydro
            Utmp(ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))*dfac
          enddo

          okay = .true.
          if(Utmp(1) .lt. subsub_dfloor) okay=.false.
          pnew = (Utmp(5) - 0.5D0/max(Utmp(1),subsub_dfloor)*(Utmp(2)**2 + Utmp(3)**2 + Utmp(4)**2))*(gamma - 1.0D0)
          if(pnew .lt. subsub_pfloor) okay=.false.

          if(.not.okay) dfac = dfac * 0.5D0

          if(okay) exit
          if(dfac .lt. 1.0D-12)then !! fail to update this cell
            do ivar=1,subsub_nhydro
              Utmp(ivar) = subsub_hydro(ii,ivar)
            enddo
            exit
          endif
        enddo
      endif

      !! Update
      do ivar=1,subsub_nhydro
        subsub_hydro(ii,ivar) = Utmp(ivar)
      enddo

    enddo
    enddo
    !$omp end do
  enddo


  !(right boundary)
    !$omp do collapse(2) private(ii, ivar, iyu, iyd, ix, iz, Utmp, okay, pnew, dfac, bc_vv, bc_cs, bc_flux, bc_cons, amaxL, amaxR, FL, FR)
    do iz=1, subsub_ngrid
    do ix=1, subsub_ngrid
  
      ii = (iz-1)*subsub_ngrid2 + (subsub_ngrid-1)*subsub_ngrid + ix
      iyu = subsub_ngrid+1
      iyd = subsub_ngrid-1

      bc_vv = hbc(2,2,3,2)
      bc_cs  = csarr_bc(2,2)
      bc_flux(:) = hbc(2,2,:,4)
      bc_cons(:) = hbc(2,2,:,1)

      !! DEBUG MODE FOR SELF-GRAVITY TEST
      !! )) DEBUGG HYDRO ((         <- this is for grep
      if(subsub_dev_hydroonly .eq. 1)then
        iyu = 1
        bc_vv = subsub_hdummy(ix,iyu,iz,3,2)
        bc_cs = subsub_csarr(ix,iyu,iz)
        bc_flux(:) = subsub_hdummy(ix,iyu,iz,:,4)
        bc_cons(:) = subsub_hdummy(ix,iyu,iz,:,1)
      endif


      amaxL = max(abs(subsub_hdummy(ix,subsub_ngrid,iz,3,2)) + subsub_csarr(ix,subsub_ngrid,iz), abs(subsub_hdummy(ix,iyd,iz,3,2)) + subsub_csarr(ix,iyd,iz))
      amaxR = max(abs(subsub_hdummy(ix,subsub_ngrid,iz,3,2)) + subsub_csarr(ix,subsub_ngrid,iz), abs(bc_vv) + bc_cs)

      do ivar=1,subsub_nhydro
        FL(ivar) = 0.5D0*(subsub_hdummy(ix,subsub_ngrid,iz,ivar,4) + subsub_hdummy(ix,iyd,iz,ivar,4)) - &
          0.5D0*amaxL*(subsub_hdummy(ix,subsub_ngrid,iz,ivar,1) - subsub_hdummy(ix,iyd,iz,ivar,1))
        FR(ivar) = 0.5D0*(bc_flux(ivar) + subsub_hdummy(ix,subsub_ngrid,iz,ivar,4)) - &
          0.5D0*amaxR*(bc_cons(ivar) - subsub_hdummy(ix,subsub_ngrid,iz,ivar,1))
      enddo


      !! Before update, check the negative density & pressure
      do ivar=1,subsub_nhydro
        Utmp(ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
      enddo

      okay = .true.
      if(Utmp(1) .lt. subsub_dfloor) okay=.false.
      pnew = (Utmp(5) - 0.5D0/max(Utmp(1),subsub_dfloor)*(Utmp(2)**2 + Utmp(3)**2 + Utmp(4)**2))*(gamma - 1.0D0)
      if(pnew .lt. subsub_pfloor) okay=.false.

      if(.not. okay)then
        !!----- RHEE -----
        !! Not physical, but limit negative pressure and density
        dfac = 0.5D0

        do
          do ivar=1,subsub_nhydro
            Utmp(ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))*dfac
          enddo

          okay = .true.
          if(Utmp(1) .lt. subsub_dfloor) okay=.false.
          pnew = (Utmp(5) - 0.5D0/max(Utmp(1),subsub_dfloor)*(Utmp(2)**2 + Utmp(3)**2 + Utmp(4)**2))*(gamma - 1.0D0)
          if(pnew .lt. subsub_pfloor) okay=.false.

          if(.not.okay) dfac = dfac * 0.5D0

          if(okay) exit
          if(dfac .lt. 1.0D-12)then !! fail to update this cell
            do ivar=1,subsub_nhydro
              Utmp(ivar) = subsub_hydro(ii,ivar)
            enddo
            exit
          endif
        enddo
      endif

      !! Update
      do ivar=1,subsub_nhydro
        subsub_hydro(ii,ivar) = Utmp(ivar)
      enddo

    enddo
    enddo
    !$omp end do

!!------------------------------------------------------
!! Z-axis
!!------------------------------------------------------
  !(left boundary)
    !$omp do collapse(2) private(ii, ivar, izu, izd, ix, iy, Utmp, okay, pnew, dfac, bc_vv, bc_cs, bc_flux, bc_cons, amaxL, amaxR, FL, FR)
    do iy=1, subsub_ngrid
    do ix=1, subsub_ngrid
  
      ii = (iy-1)*subsub_ngrid + ix
      izu = 1+1
      izd = 1-1

      bc_vv = hbc(1,3,4,2)
      bc_cs  = csarr_bc(1,3)
      bc_flux(:) = hbc(1,3,:,5)
      bc_cons(:) = hbc(1,3,:,1)

      !! DEBUG MODE FOR SELF-GRAVITY TEST
      !! )) DEBUGG HYDRO ((         <- this is for grep
      if(subsub_dev_hydroonly .eq. 1)then
        izd = subsub_ngrid
        bc_vv = subsub_hdummy(ix,iy,izd,4,2)
        bc_cs = subsub_csarr(ix,iy,izd)
        bc_flux(:) = subsub_hdummy(ix,iy,izd,:,5)
        bc_cons(:) = subsub_hdummy(ix,iy,izd,:,1)
      endif

      amaxL = max(abs(subsub_hdummy(ix,iy,1,4,2)) + subsub_csarr(ix,iy,1), abs(bc_vv) + bc_cs)
      amaxR = max(abs(subsub_hdummy(ix,iy,1,4,2)) + subsub_csarr(ix,iy,1), abs(subsub_hdummy(ix,iy,izu,4,2)) + subsub_csarr(ix,iy,izu))

      do ivar=1,subsub_nhydro
        FL(ivar) = 0.5D0*(bc_flux(ivar) + subsub_hdummy(ix,iy,1,ivar,5)) - &
          0.5D0*amaxL*(subsub_hdummy(ix,iy,1,ivar,1) - bc_cons(ivar))

        FR(ivar) = 0.5D0*(subsub_hdummy(ix,iy,1,ivar,5) + subsub_hdummy(ix,iy,izu,ivar,5)) - &
          0.5D0*amaxR*(subsub_hdummy(ix,iy,izu,ivar,1) - subsub_hdummy(ix,iy,1,ivar,1))
      enddo

      !! Before update, check the negative density & pressure
      do ivar=1,subsub_nhydro
        Utmp(ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
      enddo

      okay = .true.
      if(Utmp(1) .lt. subsub_dfloor) okay=.false.
      pnew = (Utmp(5) - 0.5D0/max(Utmp(1),subsub_dfloor)*(Utmp(2)**2 + Utmp(3)**2 + Utmp(4)**2))*(gamma - 1.0D0)
      if(pnew .lt. subsub_pfloor) okay=.false.

      if(.not. okay)then
        !!----- RHEE -----
        !! Not physical, but limit negative pressure and density
        dfac = 0.5D0

        do
          do ivar=1,subsub_nhydro
            Utmp(ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))*dfac
          enddo

          okay = .true.
          if(Utmp(1) .lt. subsub_dfloor) okay=.false.
          pnew = (Utmp(5) - 0.5D0/max(Utmp(1),subsub_dfloor)*(Utmp(2)**2 + Utmp(3)**2 + Utmp(4)**2))*(gamma - 1.0D0)
          if(pnew .lt. subsub_pfloor) okay=.false.

          if(.not.okay) dfac = dfac * 0.5D0

          if(okay) exit
          if(dfac .lt. 1.0D-12)then !! fail to update this cell
            do ivar=1,subsub_nhydro
              Utmp(ivar) = subsub_hydro(ii,ivar)
            enddo
            exit
          endif
        enddo
      endif

      !! Update
      do ivar=1,subsub_nhydro
        subsub_hydro(ii,ivar) = Utmp(ivar)
      enddo

    enddo
    enddo
    !$omp end do


  !(interior)
  do iz=2, subsub_ngrid-1
    !$omp do collapse(2) private(ii, ivar, izu, izd, ix, iy, Utmp, okay, pnew, dfac, amaxL, amaxR, FL, FR)
    do iy=1, subsub_ngrid
    do ix=1, subsub_ngrid
      ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
      izu = iz+1
      izd = iz-1

      amaxL = max(abs(subsub_hdummy(ix,iy,iz,4,2)) + subsub_csarr(ix,iy,iz), abs(subsub_hdummy(ix,iy,izd,4,2)) + subsub_csarr(ix,iy,izd))
      amaxR = max(abs(subsub_hdummy(ix,iy,iz,4,2)) + subsub_csarr(ix,iy,iz), abs(subsub_hdummy(ix,iy,izu,4,2)) + subsub_csarr(ix,iy,izu))

      do ivar=1,subsub_nhydro
        FL(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,5) + subsub_hdummy(ix,iy,izd,ivar,5)) - &
          0.5D0*amaxL*(subsub_hdummy(ix,iy,iz,ivar,1) - subsub_hdummy(ix,iy,izd,ivar,1))
        FR(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,5) + subsub_hdummy(ix,iy,izu,ivar,5)) - &
          0.5D0*amaxR*(subsub_hdummy(ix,iy,izu,ivar,1) - subsub_hdummy(ix,iy,iz,ivar,1))
      enddo

      !! Before update, check the negative density & pressure
      do ivar=1,subsub_nhydro
        Utmp(ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
      enddo

      okay = .true.
      if(Utmp(1) .lt. subsub_dfloor) okay=.false.
      pnew = (Utmp(5) - 0.5D0/max(Utmp(1),subsub_dfloor)*(Utmp(2)**2 + Utmp(3)**2 + Utmp(4)**2))*(gamma - 1.0D0)
      if(pnew .lt. subsub_pfloor) okay=.false.

      if(.not. okay)then
        !!----- RHEE -----
        !! Not physical, but limit negative pressure and density
        dfac = 0.5D0

        do
          do ivar=1,subsub_nhydro
            Utmp(ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))*dfac
          enddo

          okay = .true.
          if(Utmp(1) .lt. subsub_dfloor) okay=.false.
          pnew = (Utmp(5) - 0.5D0/max(Utmp(1),subsub_dfloor)*(Utmp(2)**2 + Utmp(3)**2 + Utmp(4)**2))*(gamma - 1.0D0)
          if(pnew .lt. subsub_pfloor) okay=.false.

          if(.not.okay) dfac = dfac * 0.5D0

          if(okay) exit
          if(dfac .lt. 1.0D-12)then !! fail to update this cell
            do ivar=1,subsub_nhydro
              Utmp(ivar) = subsub_hydro(ii,ivar)
            enddo
            exit
          endif
        enddo
      endif

      !! Update
      do ivar=1,subsub_nhydro
        subsub_hydro(ii,ivar) = Utmp(ivar)
      enddo

    enddo
    enddo
    !$omp end do
  enddo

  !(right boundary)
    !$omp do collapse(2) private(ii, ivar, izu, izd, ix, iy, Utmp, okay, pnew, dfac, bc_vv, bc_cs, bc_flux, bc_cons, amaxL, amaxR, FL, FR)
    do iy=1, subsub_ngrid
    do ix=1, subsub_ngrid
  
      ii = (subsub_ngrid-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
      izu = subsub_ngrid+1
      izd = subsub_ngrid-1

      bc_vv = hbc(2,3,4,2)
      bc_cs  = csarr_bc(2,3)
      bc_flux(:) = hbc(2,3,:,5)
      bc_cons(:) = hbc(2,3,:,1)

      !! DEBUG MODE FOR SELF-GRAVITY TEST
      !! )) DEBUGG HYDRO ((         <- this is for grep
      if(subsub_dev_hydroonly .eq. 1)then
        izu = 1
        bc_vv = subsub_hdummy(ix,iy,izu,4,2)
        bc_cs = subsub_csarr(ix,iy,izu)
        bc_flux(:) = subsub_hdummy(ix,iy,izu,:,5)
        bc_cons(:) = subsub_hdummy(ix,iy,izu,:,1)
      endif


      amaxL = max(abs(subsub_hdummy(ix,iy,subsub_ngrid,4,2)) + subsub_csarr(ix,iy,subsub_ngrid), abs(subsub_hdummy(ix,iy,izd,4,2)) + subsub_csarr(ix,iy,izd))
      amaxR = max(abs(subsub_hdummy(ix,iy,subsub_ngrid,4,2)) + subsub_csarr(ix,iy,subsub_ngrid), abs(bc_vv) + bc_cs)

      do ivar=1,subsub_nhydro
        FL(ivar) = 0.5D0*(subsub_hdummy(ix,iy,subsub_ngrid,ivar,5) + subsub_hdummy(ix,iy,izd,ivar,5)) - &
          0.5D0*amaxL*(subsub_hdummy(ix,iy,subsub_ngrid,ivar,1) - subsub_hdummy(ix,iy,izd,ivar,1))
        FR(ivar) = 0.5D0*(bc_flux(ivar) + subsub_hdummy(ix,iy,subsub_ngrid,ivar,5)) - &
          0.5D0*amaxR*(bc_cons(ivar) - subsub_hdummy(ix,iy,subsub_ngrid,ivar,1))
      enddo


      !! Before update, check the negative density & pressure
      do ivar=1,subsub_nhydro
        Utmp(ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
      enddo

      okay = .true.
      if(Utmp(1) .lt. subsub_dfloor) okay=.false.
      pnew = (Utmp(5) - 0.5D0/max(Utmp(1),subsub_dfloor)*(Utmp(2)**2 + Utmp(3)**2 + Utmp(4)**2))*(gamma - 1.0D0)
      if(pnew .lt. subsub_pfloor) okay=.false.

      if(.not. okay)then
        !!----- RHEE -----
        !! Not physical, but limit negative pressure and density
        dfac = 0.5D0

        do
          do ivar=1,subsub_nhydro
            Utmp(ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))*dfac
          enddo

          okay = .true.
          if(Utmp(1) .lt. subsub_dfloor) okay=.false.
          pnew = (Utmp(5) - 0.5D0/max(Utmp(1),subsub_dfloor)*(Utmp(2)**2 + Utmp(3)**2 + Utmp(4)**2))*(gamma - 1.0D0)
          if(pnew .lt. subsub_pfloor) okay=.false.

          if(.not.okay) dfac = dfac * 0.5D0

          if(okay) exit
          if(dfac .lt. 1.0D-12)then !! fail to update this cell
            do ivar=1,subsub_nhydro
              Utmp(ivar) = subsub_hydro(ii,ivar)
            enddo
            exit
          endif
        enddo
      endif

      !! Update
      do ivar=1,subsub_nhydro
        subsub_hydro(ii,ivar) = Utmp(ivar)
      enddo

    enddo
    enddo
    !$omp end do


  return

!!----- RHEE -----
!! Solver in computationally efficient version below by minimizing flux calculations
!! Test required
!!----------------
!!------------------------------------------------------
!! X-axis
!!------------------------------------------------------
  !! (left boundary)
  !$omp do collapse(2) private(ii, xu, ivar, iy, iz, amaxL, amaxR, FL, FR, Utmp)
  do iy=1, subsub_ngrid
  do iz=1, subsub_ngrid
    ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + 1
    xu = ii + 1

    amaxL = max(abs(subsub_hdummy(1,iy,iz,2,2)) + subsub_csarr(1,iy,iz), abs(hbc(1,1,2,2)) + csarr_bc(1,1))
    amaxR = max(abs(subsub_hdummy(1,iy,iz,2,2)) + subsub_csarr(1,iy,iz), abs(subsub_hdummy(2,iy,iz,2,2)) + subsub_csarr(2,iy,iz))

    do ivar=1,subsub_nhydro
      FL(ivar) = 0.5D0*(hbc(1,1,ivar,3) + subsub_hdummy(1,iy,iz,ivar,3)) - &
        0.5D0*amaxL*(subsub_hdummy(1,iy,iz,ivar,1) - hbc(1,1,ivar,1))
      FR(ivar) = 0.5D0*(subsub_hdummy(1,iy,iz,ivar,3) + subsub_hdummy(2,iy,iz,ivar,3)) - &
        0.5D0*amaxR*(subsub_hdummy(2,iy,iz,ivar,1) - subsub_hdummy(1,iy,iz,ivar,1))
    enddo

    
    do ivar=1,subsub_nhydro
      subsub_hydro(ii,ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
      subsub_hydro(xu,ivar) = subsub_hydro(xu,ivar) + lam * FR(ivar)
    enddo
  enddo
  enddo
  !$omp end do
 
  !! (interior) 
  do ix=3, subsub_ngrid-1, 2
    !$omp do collapse(2) private(ii, xu, xd, ivar, iy, iz, amaxL, amaxR, FL, FR)
    do iy=1, subsub_ngrid
    do iz=1, subsub_ngrid
      ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
      xu = ii + 1
      xd = ii - 1

      amaxL = max(abs(subsub_hdummy(ix,iy,iz,2,2)) + subsub_csarr(ix,iy,iz), abs(subsub_hdummy(ix-1,iy,iz,2,2)) + subsub_csarr(ix-1,iy,iz))
      amaxR = max(abs(subsub_hdummy(ix,iy,iz,2,2)) + subsub_csarr(ix,iy,iz), abs(subsub_hdummy(ix+1,iy,iz,2,2)) + subsub_csarr(ix+1,iy,iz))

      do ivar=1,subsub_nhydro
        FL(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,3) + subsub_hdummy(ix-1,iy,iz,ivar,3)) - &
          0.5D0*amaxL*(subsub_hdummy(ix,iy,iz,ivar,1) - subsub_hdummy(ix-1,iy,iz,ivar,1))
        FR(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,3) + subsub_hdummy(ix+1,iy,iz,ivar,3)) - &
          0.5D0*amaxR*(subsub_hdummy(ix+1,iy,iz,ivar,1) - subsub_hdummy(ix,iy,iz,ivar,1))
      enddo

      do ivar=1,subsub_nhydro
        subsub_hydro(xd,ivar) = subsub_hydro(xd,ivar) - lam * FL(ivar)
        subsub_hydro(ii,ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
        subsub_hydro(xu,ivar) = subsub_hydro(xu,ivar) + lam * FR(ivar)
      enddo
    enddo
    enddo
    !$omp end do
  enddo

  !! (right boundary)
  !$omp do collapse(2) private(ii, xd, ivar, iy, iz, amaxL, amaxR, FL, FR)
  do iy=1, subsub_ngrid
  do iz=1, subsub_ngrid
    ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + subsub_ngrid
    xd = ii - 1

    amaxL = max(abs(subsub_hdummy(subsub_ngrid,iy,iz,2,2)) + subsub_csarr(subsub_ngrid,iy,iz), abs(subsub_hdummy(subsub_ngrid-1,iy,iz,2,2)) + subsub_csarr(subsub_ngrid-1,iy,iz))
    amaxR = max(abs(subsub_hdummy(subsub_ngrid,iy,iz,2,2)) + subsub_csarr(subsub_ngrid,iy,iz), abs(hbc(2,1,2,2)) + csarr_bc(2,1))

    do ivar=1, 5
      FL(ivar) = 0.5D0*(subsub_hdummy(subsub_ngrid,iy,iz,ivar,3) + subsub_hdummy(subsub_ngrid-1,iy,iz,ivar,3)) - &
        0.5D0*amaxL*(subsub_hdummy(subsub_ngrid,iy,iz,ivar,1) - subsub_hdummy(subsub_ngrid-1,iy,iz,ivar,1))
      FR(ivar) = 0.5D0*(subsub_hdummy(subsub_ngrid,iy,iz,ivar,3) + hbc(2,1,ivar,3)) - &
        0.5D0*amaxR*(hbc(2,1,ivar,1) - subsub_hdummy(subsub_ngrid,iy,iz,ivar,1))
    enddo

    do ivar=1, 5
      subsub_hydro(ii,ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
      if(isodd)then
        subsub_hydro(xd,ivar) = subsub_hydro(xd,ivar) - lam*FL(ivar)
      endif
    enddo
  enddo
  enddo
  !$omp end do

!!------------------------------------------------------
!! Y-axis
!!------------------------------------------------------
  !! (left boundary)
  !$omp do collapse(2) private(ii, yu, ivar, ix, iz, amaxL, amaxR, FL, FR)
  do ix=1, subsub_ngrid
  do iz=1, subsub_ngrid
    ii = (iz-1)*subsub_ngrid2 + ix
    yu = ii + subsub_ngrid

    amaxL = max(abs(subsub_hdummy(ix,1,iz,3,2)) + subsub_csarr(ix,1,iz), abs(hbc(1,2,3,2)) + csarr_bc(1,2))
    amaxR = max(abs(subsub_hdummy(ix,1,iz,3,2)) + subsub_csarr(ix,1,iz), abs(subsub_hdummy(ix,2,iz,3,2)) + subsub_csarr(ix,2,iz))

    do ivar=1,subsub_nhydro
      FL(ivar) = 0.5D0*(hbc(1,2,ivar,4) + subsub_hdummy(ix,1,iz,ivar,4)) - &
        0.5D0*amaxL*(subsub_hdummy(ix,1,iz,ivar,1) - hbc(1,2,ivar,1))
      FR(ivar) = 0.5D0*(subsub_hdummy(ix,1,iz,ivar,4) + subsub_hdummy(ix,2,iz,ivar,4)) - &
        0.5D0*amaxR*(subsub_hdummy(ix,2,iz,ivar,1) - subsub_hdummy(ix,1,iz,ivar,1))
    enddo

    
    do ivar=1,subsub_nhydro
      subsub_hydro(ii,ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
      subsub_hydro(yu,ivar) = subsub_hydro(yu,ivar) + lam * FR(ivar)
    enddo
  enddo
  enddo
  !$omp end do
 
  !! (interior) 
  do iy=3, subsub_ngrid-1, 2
    !$omp do collapse(2) private(ii, yu, yd, ix, iz, ivar, amaxL, amaxR, FL, FR)
    do iz=1, subsub_ngrid
    do ix=1, subsub_ngrid
      ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
      yu = ii + subsub_ngrid
      yd = ii - subsub_ngrid

      amaxL = max(abs(subsub_hdummy(ix,iy,iz,3,2)) + subsub_csarr(ix,iy,iz), abs(subsub_hdummy(ix,iy-1,iz,3,2)) + subsub_csarr(ix,iy-1,iz))
      amaxR = max(abs(subsub_hdummy(ix,iy,iz,3,2)) + subsub_csarr(ix,iy,iz), abs(subsub_hdummy(ix,iy+1,iz,3,2)) + subsub_csarr(ix,iy+1,iz))

      do ivar=1,subsub_nhydro
        FL(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,4) + subsub_hdummy(ix,iy-1,iz,ivar,4)) - &
          0.5D0*amaxL*(subsub_hdummy(ix,iy,iz,ivar,1) - subsub_hdummy(ix,iy-1,iz,ivar,1))
        FR(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,4) + subsub_hdummy(ix,iy+1,iz,ivar,4)) - &
          0.5D0*amaxR*(subsub_hdummy(ix,iy+1,iz,ivar,1) - subsub_hdummy(ix,iy,iz,ivar,1))
      enddo

      do ivar=1,subsub_nhydro
        subsub_hydro(yd,ivar) = subsub_hydro(yd,ivar) - lam * FL(ivar)
        subsub_hydro(ii,ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
        subsub_hydro(yu,ivar) = subsub_hydro(yu,ivar) + lam * FR(ivar)
      enddo
    enddo
    enddo
    !$omp end do
  enddo

  !! (right boundary)
  !$omp do collapse(2) private(ii, yd, ivar, ix, iz, amaxL, amaxR, FL, FR)
  do iz=1, subsub_ngrid
  do ix=1, subsub_ngrid
    ii = (iz-1)*subsub_ngrid2 + (subsub_ngrid-1)*subsub_ngrid + ix
    yd = ii - subsub_ngrid

    amaxL = max(abs(subsub_hdummy(ix,subsub_ngrid,iz,3,2)) + subsub_csarr(ix,subsub_ngrid,iz), abs(subsub_hdummy(ix,subsub_ngrid-1,iz,3,2)) + subsub_csarr(ix,subsub_ngrid-1,iz))
    amaxR = max(abs(subsub_hdummy(ix,subsub_ngrid,iz,3,2)) + subsub_csarr(ix,subsub_ngrid,iz), abs(hbc(2,2,3,2)) + csarr_bc(2,2))

    do ivar=1, 5
      FL(ivar) = 0.5D0*(subsub_hdummy(ix,subsub_ngrid,iz,ivar,4) + subsub_hdummy(ix,subsub_ngrid-1,iz,ivar,4)) - &
        0.5D0*amaxL*(subsub_hdummy(ix,subsub_ngrid,iz,ivar,1) - subsub_hdummy(ix,subsub_ngrid-1,iz,ivar,1))
      FR(ivar) = 0.5D0*(subsub_hdummy(ix,subsub_ngrid,iz,ivar,4) + hbc(2,2,ivar,4)) - &
        0.5D0*amaxR*(hbc(2,2,ivar,1) - subsub_hdummy(ix,subsub_ngrid,iz,ivar,1))
    enddo

    do ivar=1, 5
      subsub_hydro(ii,ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
      if(isodd)then
        subsub_hydro(yd,ivar) = subsub_hydro(yd,ivar) - lam*FL(ivar)
      endif
    enddo
  enddo
  enddo
  !$omp end do

!!------------------------------------------------------
!! Z-axis
!!------------------------------------------------------
  !! (left boundary)
  !$omp do collapse(2) private(ii, zu, ivar, ix, iy, amaxL, amaxR, FL, FR)
  do iy=1, subsub_ngrid
  do ix=1, subsub_ngrid
    ii = (iy-1)*subsub_ngrid + ix
    zu = ii + subsub_ngrid2

    amaxL = max(abs(subsub_hdummy(ix,iy,1,4,2)) + subsub_csarr(ix,iy,1), abs(hbc(1,3,4,2)) + csarr_bc(1,3))
    amaxR = max(abs(subsub_hdummy(ix,iy,1,4,2)) + subsub_csarr(ix,iy,1), abs(subsub_hdummy(ix,iy,2,4,2)) + subsub_csarr(ix,iy,2))

    do ivar=1,subsub_nhydro
      FL(ivar) = 0.5D0*(hbc(1,3,ivar,5) + subsub_hdummy(ix,iy,1,ivar,5)) - &
        0.5D0*amaxL*(subsub_hdummy(ix,iy,1,ivar,1) - hbc(1,3,ivar,1))
      FR(ivar) = 0.5D0*(subsub_hdummy(ix,iy,1,ivar,5) + subsub_hdummy(ix,iy,2,ivar,5)) - &
        0.5D0*amaxR*(subsub_hdummy(ix,iy,2,ivar,1) - subsub_hdummy(ix,iy,1,ivar,1))
    enddo

    
    do ivar=1,subsub_nhydro
      subsub_hydro(ii,ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
      subsub_hydro(zu,ivar) = subsub_hydro(zu,ivar) + lam * FR(ivar)
    enddo
  enddo
  enddo
  !$omp end do
 
  !! (interior) 
  do iz=3, subsub_ngrid-1, 2
    !$omp do collapse(2) private(ii, zu, zd, ivar, ix, iy, amaxL, amaxR, FL, FR)
    do iy=1, subsub_ngrid
    do ix=1, subsub_ngrid
      ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
      zu = ii + subsub_ngrid2
      zd = ii - subsub_ngrid2

      amaxL = max(abs(subsub_hdummy(ix,iy,iz,4,2)) + subsub_csarr(ix,iy,iz), abs(subsub_hdummy(ix,iy,iz-1,4,2)) + subsub_csarr(ix,iy,iz-1))
      amaxR = max(abs(subsub_hdummy(ix,iy,iz,4,2)) + subsub_csarr(ix,iy,iz), abs(subsub_hdummy(ix,iy,iz+1,4,2)) + subsub_csarr(ix,iy,iz+1))

      do ivar=1,subsub_nhydro
        FL(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,5) + subsub_hdummy(ix,iy,iz-1,ivar,5)) - &
          0.5D0*amaxL*(subsub_hdummy(ix,iy,iz,ivar,1) - subsub_hdummy(ix,iy,iz-1,ivar,1))
        FR(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,5) + subsub_hdummy(ix,iy,iz+1,ivar,5)) - &
          0.5D0*amaxR*(subsub_hdummy(ix,iy,iz+1,ivar,1) - subsub_hdummy(ix,iy,iz,ivar,1))
      enddo

      do ivar=1,subsub_nhydro
        subsub_hydro(zd,ivar) = subsub_hydro(zd,ivar) - lam * FL(ivar)
        subsub_hydro(ii,ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
        subsub_hydro(zu,ivar) = subsub_hydro(zu,ivar) + lam * FR(ivar)
      enddo
    enddo
    enddo
    !$omp end do
  enddo

  !! (right boundary)
  !$omp do collapse(2) private(ii, zd, ivar, ix, iy, amaxL, amaxR, FL, FR)
  do iy=1, subsub_ngrid
  do ix=1, subsub_ngrid
    ii = (subsub_ngrid-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
    zd = ii - subsub_ngrid2

    amaxL = max(abs(subsub_hdummy(ix,iy,subsub_ngrid,4,2)) + subsub_csarr(ix,iy,subsub_ngrid), abs(subsub_hdummy(ix,iy,subsub_ngrid-1,4,2)) + subsub_csarr(ix,iy,subsub_ngrid-1))
    amaxR = max(abs(subsub_hdummy(ix,iy,subsub_ngrid,4,2)) + subsub_csarr(ix,iy,subsub_ngrid), abs(hbc(2,3,4,2)) + csarr_bc(2,3))

    do ivar=1, 5
      FL(ivar) = 0.5D0*(subsub_hdummy(ix,iy,subsub_ngrid,ivar,5) + subsub_hdummy(ix,iy,subsub_ngrid-1,ivar,5)) - &
        0.5D0*amaxL*(subsub_hdummy(ix,iy,subsub_ngrid,ivar,1) - subsub_hdummy(ix,iy,subsub_ngrid-1,ivar,1))
      FR(ivar) = 0.5D0*(subsub_hdummy(ix,iy,subsub_ngrid,ivar,5) + hbc(2,3,ivar,5)) - &
        0.5D0*amaxR*(hbc(2,3,ivar,1) - subsub_hdummy(ix,iy,subsub_ngrid,ivar,1))
    enddo

    do ivar=1, 5
      subsub_hydro(ii,ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
      if(isodd)then
        subsub_hydro(zd,ivar) = subsub_hydro(zd,ivar) - lam*FL(ivar)
      endif
    enddo
  enddo
  enddo
  !$omp end do
end subroutine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_enforce_floors(U)
  use subsub_commons
  use subsub_parameters
  use hydro_parameters, ONLY: gamma
  implicit none
  real(dp), intent(inout), dimension(1:subsub_nhydro) :: U
  
  real(dp) :: rho, vx, vy, vz, ekin

  if(U(1) .lt. subsub_dfloor)then
    U(1) = subsub_dfloor
    U(2) = 0.0D0
    U(3) = 0.0D0
    U(4) = 0.0D0
    U(5) = subsub_pfloor/(gamma-1.0D0)
  endif
end subroutine subsub_enforce_floors
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_hydroRiemann_Rusanov_periodicBC(dt)!, hbc)
  use amr_commons
  use subsub_commons
  use subsub_parameters
  use hydro_parameters, ONLY:gamma
  use, intrinsic :: ieee_arithmetic
  implicit none

  real(dp) :: dt
  real(dp), dimension(1:2, 1:ndim, 1:subsub_nhydro, 1:5) :: hbc
  !! Local variables
  integer :: ivar, i, j
  integer :: ix, iy, iz
  integer :: ii, xu, xd, yu, yd, zu, zd
  integer :: ixd, ixu, iyd, iyu, izd, izu
  real(dp) :: amaxL, amaxR, lam, ekin
  
  real(dp), dimension(1:subsub_nhydro) :: FL, FR, Utmp
  real(dp), dimension(1:2, 1:ndim) :: csarr_bc
  logical :: isodd


  !1 conservative old
  !2 primitive old
  !3 Fx
  !4 Fy
  !5 Fz


  lam = dt / subsub_dx
  isodd = mod(subsub_ngrid,2) .eq. 1

    !!----- Compute sound speed first
  !$omp do collapse(2) private(ix, iy, iz)
  do iz=1, subsub_ngrid
  do iy=1, subsub_ngrid
  do ix=1, subsub_ngrid
    if(subsub_hdummy(ix,iy,iz,1,2) .gt. subsub_dfloor) then
      subsub_csarr(ix,iy,iz) = sqrt(gamma * subsub_hdummy(ix,iy,iz,5,2) / subsub_hdummy(ix,iy,iz,1,2))
    else
      subsub_csarr(ix,iy,iz) = sqrt(gamma * subsub_pfloor / subsub_dfloor)
    endif
  enddo
  enddo
  enddo
  !$omp end do

  !! DEBUG MODE FOR SELF-GRAVITY TEST
  !! )) DEBUGG GRAV ((         <- this is for grep
  if(subsub_dev_gravonly .eq. 1)then
    !$omp do collapse(2) private(ix, iy, iz)
    do iz=1, subsub_ngrid
    do iy=1, subsub_ngrid
    do ix=1, subsub_ngrid
      subsub_csarr(ix,iy,iz) = 0.0D0
    enddo
    enddo
    enddo
    !$omp end do
  endif


  !!----- Compute Flux
  !$omp do collapse(2) private(ix, iy, iz)
  do iz=1, subsub_ngrid
  do iy=1, subsub_ngrid
  do ix=1, subsub_ngrid
    subsub_hdummy(ix,iy,iz,1,3) = subsub_hdummy(ix,iy,iz,2,1)
    subsub_hdummy(ix,iy,iz,2,3) = subsub_hdummy(ix,iy,iz,2,1)*subsub_hdummy(ix,iy,iz,2,2) + subsub_hdummy(ix,iy,iz,5,2)
    subsub_hdummy(ix,iy,iz,3,3) = subsub_hdummy(ix,iy,iz,3,1)*subsub_hdummy(ix,iy,iz,2,2)
    subsub_hdummy(ix,iy,iz,4,3) = subsub_hdummy(ix,iy,iz,4,1)*subsub_hdummy(ix,iy,iz,2,2)
    subsub_hdummy(ix,iy,iz,5,3) = (subsub_hdummy(ix,iy,iz,5,1) + subsub_hdummy(ix,iy,iz,5,2)) * subsub_hdummy(ix,iy,iz,2,2)

    subsub_hdummy(ix,iy,iz,1,4) = subsub_hdummy(ix,iy,iz,3,1)
    subsub_hdummy(ix,iy,iz,2,4) = subsub_hdummy(ix,iy,iz,2,1)*subsub_hdummy(ix,iy,iz,3,2)
    subsub_hdummy(ix,iy,iz,3,4) = subsub_hdummy(ix,iy,iz,3,1)*subsub_hdummy(ix,iy,iz,3,2) + subsub_hdummy(ix,iy,iz,5,2)
    subsub_hdummy(ix,iy,iz,4,4) = subsub_hdummy(ix,iy,iz,4,1)*subsub_hdummy(ix,iy,iz,3,2)
    subsub_hdummy(ix,iy,iz,5,4) = (subsub_hdummy(ix,iy,iz,5,1) + subsub_hdummy(ix,iy,iz,5,2)) * subsub_hdummy(ix,iy,iz,3,2)

    subsub_hdummy(ix,iy,iz,1,5) = subsub_hdummy(ix,iy,iz,4,1)
    subsub_hdummy(ix,iy,iz,2,5) = subsub_hdummy(ix,iy,iz,2,1)*subsub_hdummy(ix,iy,iz,4,2)
    subsub_hdummy(ix,iy,iz,3,5) = subsub_hdummy(ix,iy,iz,3,1)*subsub_hdummy(ix,iy,iz,4,2)
    subsub_hdummy(ix,iy,iz,4,5) = subsub_hdummy(ix,iy,iz,4,1)*subsub_hdummy(ix,iy,iz,4,2) + subsub_hdummy(ix,iy,iz,5,2)
    subsub_hdummy(ix,iy,iz,5,5) = (subsub_hdummy(ix,iy,iz,5,1) + subsub_hdummy(ix,iy,iz,5,2)) * subsub_hdummy(ix,iy,iz,4,2)
  enddo
  enddo
  enddo
  !$omp end do

  !! Solver
  do ix=1, subsub_ngrid
    !$omp do collapse(2) private(ii, ixu, ixd, ivar, iy, iz, amaxL, amaxR, FL, FR, Utmp)
    do iz=1, subsub_ngrid
    do iy=1, subsub_ngrid
  
      ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
      ixu = ix+1
      ixd = ix-1

      if(ix.eq.1) ixd = subsub_ngrid
      if(ix.eq.subsub_ngrid) ixu = 1

      amaxL = max(abs(subsub_hdummy(ix,iy,iz,2,2)) + subsub_csarr(ix,iy,iz), abs(subsub_hdummy(ixd,iy,iz,2,2)) + subsub_csarr(ixd,iy,iz))
      amaxR = max(abs(subsub_hdummy(ix,iy,iz,2,2)) + subsub_csarr(ix,iy,iz), abs(subsub_hdummy(ixu,iy,iz,2,2)) + subsub_csarr(ixu,iy,iz))

      do ivar=1,subsub_nhydro
        FL(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,3) + subsub_hdummy(ixd,iy,iz,ivar,3)) - &
          0.5D0*amaxL*(subsub_hdummy(ix,iy,iz,ivar,1) - subsub_hdummy(ixd,iy,iz,ivar,1))
        FR(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,3) + subsub_hdummy(ixu,iy,iz,ivar,3)) - &
          0.5D0*amaxR*(subsub_hdummy(ixu,iy,iz,ivar,1) - subsub_hdummy(ix,iy,iz,ivar,1))
      enddo

      do ivar=1,subsub_nhydro
        subsub_hydro(ii,ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
      enddo
    enddo
    enddo
    !$omp end do
  enddo

  do iy=1, subsub_ngrid
    !$omp do collapse(2) private(ii, iyu, iyd, ivar, ix, iz, amaxL, amaxR, FL, FR, Utmp)
    do iz=1, subsub_ngrid
    do ix=1, subsub_ngrid
  
      ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
      iyu = iy+1
      iyd = iy-1

      if(iy.eq.1) iyd = subsub_ngrid
      if(iy.eq.subsub_ngrid) iyu = 1

      amaxL = max(abs(subsub_hdummy(ix,iy,iz,3,2)) + subsub_csarr(ix,iy,iz), abs(subsub_hdummy(ix,iyd,iz,3,2)) + subsub_csarr(ix,iyd,iz))
      amaxR = max(abs(subsub_hdummy(ix,iy,iz,3,2)) + subsub_csarr(ix,iy,iz), abs(subsub_hdummy(ix,iyu,iz,3,2)) + subsub_csarr(ix,iyu,iz))

      do ivar=1,subsub_nhydro
        FL(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,4) + subsub_hdummy(ix,iyd,iz,ivar,4)) - &
          0.5D0*amaxL*(subsub_hdummy(ix,iy,iz,ivar,1) - subsub_hdummy(ix,iyd,iz,ivar,1))
        FR(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,4) + subsub_hdummy(ix,iyu,iz,ivar,4)) - &
          0.5D0*amaxR*(subsub_hdummy(ix,iyu,iz,ivar,1) - subsub_hdummy(ix,iy,iz,ivar,1))
      enddo

      do ivar=1,subsub_nhydro
        subsub_hydro(ii,ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
      enddo
    enddo
    enddo
    !$omp end do
  enddo

  do iz=1, subsub_ngrid
    !$omp do collapse(2) private(ii, izu, izd, ivar, ix, iy, amaxL, amaxR, FL, FR, Utmp)
    do iy=1, subsub_ngrid
    do ix=1, subsub_ngrid
  
      ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
      izu = iz+1
      izd = iz-1

      if(iz.eq.1) izd = subsub_ngrid
      if(iz.eq.subsub_ngrid) izu = 1

      amaxL = max(abs(subsub_hdummy(ix,iy,iz,4,2)) + subsub_csarr(ix,iy,iz), abs(subsub_hdummy(ix,iy,izd,4,2)) + subsub_csarr(ix,iy,izd))
      amaxR = max(abs(subsub_hdummy(ix,iy,iz,4,2)) + subsub_csarr(ix,iy,iz), abs(subsub_hdummy(ix,iy,izu,4,2)) + subsub_csarr(ix,iy,izu))

      do ivar=1,subsub_nhydro
        FL(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,5) + subsub_hdummy(ix,iy,izd,ivar,5)) - &
          0.5D0*amaxL*(subsub_hdummy(ix,iy,iz,ivar,1) - subsub_hdummy(ix,iy,izd,ivar,1))
        FR(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,5) + subsub_hdummy(ix,iy,izu,ivar,5)) - &
          0.5D0*amaxR*(subsub_hdummy(ix,iy,izu,ivar,1) - subsub_hdummy(ix,iy,iz,ivar,1))
      enddo

      do ivar=1,subsub_nhydro
        subsub_hydro(ii,ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
      enddo
    enddo
    enddo
    !$omp end do
  enddo
  return


!!------------------------------------------------------
!! X-axis
!!------------------------------------------------------
  !! (left boundary)
  !$omp do collapse(2) private(ii, xu, xd, ivar, iy, iz, amaxL, amaxR, FL, FR, Utmp)
  do iy=1, subsub_ngrid
  do iz=1, subsub_ngrid
    ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + 1
    xu = ii + 1
    xd = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + subsub_ngrid

    amaxL = max(abs(subsub_hdummy(1,iy,iz,2,2)) + subsub_csarr(1,iy,iz), abs(subsub_hdummy(subsub_ngrid,iy,iz,2,2)) + subsub_csarr(subsub_ngrid,iy,iz))
    amaxR = max(abs(subsub_hdummy(1,iy,iz,2,2)) + subsub_csarr(1,iy,iz), abs(subsub_hdummy(2,iy,iz,2,2)) + subsub_csarr(2,iy,iz))

    do ivar=1,subsub_nhydro
      FL(ivar) = 0.5D0*(subsub_hdummy(subsub_ngrid,iy,iz,ivar,3) + subsub_hdummy(1,iy,iz,ivar,3)) - &
        0.5D0*amaxL*(subsub_hdummy(1,iy,iz,ivar,1) - subsub_hdummy(subsub_ngrid,iy,iz,ivar,1))
      FR(ivar) = 0.5D0*(subsub_hdummy(1,iy,iz,ivar,3) + subsub_hdummy(2,iy,iz,ivar,3)) - &
        0.5D0*amaxR*(subsub_hdummy(2,iy,iz,ivar,1) - subsub_hdummy(1,iy,iz,ivar,1))
    enddo

    
    do ivar=1,subsub_nhydro
      subsub_hydro(xd,ivar) = subsub_hydro(xd,ivar) - lam * FL(ivar)
      subsub_hydro(ii,ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
      subsub_hydro(xu,ivar) = subsub_hydro(xu,ivar) + lam * FR(ivar)
    enddo
  enddo
  enddo
  !$omp end do


  !! (interior) 
  do ix=3, subsub_ngrid-1, 2
    !$omp do collapse(2) private(ii, xu, xd, ivar, iy, iz, amaxL, amaxR, FL, FR)
    do iy=1, subsub_ngrid
    do iz=1, subsub_ngrid
      ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
      xu = ii + 1
      xd = ii - 1

      amaxL = max(abs(subsub_hdummy(ix,iy,iz,2,2)) + subsub_csarr(ix,iy,iz), abs(subsub_hdummy(ix-1,iy,iz,2,2)) + subsub_csarr(ix-1,iy,iz))
      amaxR = max(abs(subsub_hdummy(ix,iy,iz,2,2)) + subsub_csarr(ix,iy,iz), abs(subsub_hdummy(ix+1,iy,iz,2,2)) + subsub_csarr(ix+1,iy,iz))

      do ivar=1,subsub_nhydro
        FL(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,3) + subsub_hdummy(ix-1,iy,iz,ivar,3)) - &
          0.5D0*amaxL*(subsub_hdummy(ix,iy,iz,ivar,1) - subsub_hdummy(ix-1,iy,iz,ivar,1))
        FR(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,3) + subsub_hdummy(ix+1,iy,iz,ivar,3)) - &
          0.5D0*amaxR*(subsub_hdummy(ix+1,iy,iz,ivar,1) - subsub_hdummy(ix,iy,iz,ivar,1))
      enddo

      do ivar=1,subsub_nhydro
        subsub_hydro(xd,ivar) = subsub_hydro(xd,ivar) - lam * FL(ivar)
        subsub_hydro(ii,ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
        subsub_hydro(xu,ivar) = subsub_hydro(xu,ivar) + lam * FR(ivar)
if(ieee_is_nan(subsub_hydro(xd,ivar)))then
  write(*,*) ix, iy, iz, 'xd', myid
  write(*,*) 'FL = ', FL(ivar)
  write(*,*) 'amaxL = ', amaxL
  write(*,*) 'subsub_hdummy = ', subsub_hdummy(ix,iy,iz,1,1), subsub_hdummy(ix-1,iy,iz,1,1), subsub_hdummy(ix+1,iy,iz,1,1)
  write(*,*) 'subsub_hdummy = ', subsub_hdummy(ix,iy,iz,2,2), subsub_hdummy(ix-1,iy,iz,2,2), subsub_hdummy(ix+1,iy,iz,2,2)
  write(*,*) 'subsub_csarr = ', subsub_csarr(ix,iy,iz), subsub_csarr(ix-1,iy,iz), subsub_csarr(ix+1,iy,iz)
  write(*,*) 'FR = ', FR(ivar)
  stop
endif

      enddo
    enddo
    enddo
    !$omp end do
  enddo

  !! (right boundary)
  !$omp do collapse(2) private(ii, xu, xd, ivar, iy, iz, amaxL, amaxR, FL, FR)
  do iy=1, subsub_ngrid
  do iz=1, subsub_ngrid
    ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + subsub_ngrid
    xd = ii - 1
    xu = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + 1

    amaxL = max(abs(subsub_hdummy(subsub_ngrid,iy,iz,2,2)) + subsub_csarr(subsub_ngrid,iy,iz), abs(subsub_hdummy(subsub_ngrid-1,iy,iz,2,2)) + subsub_csarr(subsub_ngrid-1,iy,iz))
    amaxR = max(abs(subsub_hdummy(subsub_ngrid,iy,iz,2,2)) + subsub_csarr(subsub_ngrid,iy,iz), abs(subsub_hdummy(1,iy,iz,2,2)) + subsub_csarr(1,iy,iz))

    do ivar=1, 5
      FL(ivar) = 0.5D0*(subsub_hdummy(subsub_ngrid,iy,iz,ivar,3) + subsub_hdummy(subsub_ngrid-1,iy,iz,ivar,3)) - &
        0.5D0*amaxL*(subsub_hdummy(subsub_ngrid,iy,iz,ivar,1) - subsub_hdummy(subsub_ngrid-1,iy,iz,ivar,1))
      FR(ivar) = 0.5D0*(subsub_hdummy(subsub_ngrid,iy,iz,ivar,3) + subsub_hdummy(1,iy,iz,ivar,3)) - &
        0.5D0*amaxR*(subsub_hdummy(1,iy,iz,ivar,1) - subsub_hdummy(subsub_ngrid,iy,iz,ivar,1))
    enddo

    do ivar=1, 5
      subsub_hydro(xu,ivar) = subsub_hydro(xu,ivar) + lam * FR(ivar)
      subsub_hydro(ii,ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
      if(isodd)then
        subsub_hydro(xd,ivar) = subsub_hydro(xd,ivar) - lam*FL(ivar)
      endif
    enddo
  enddo
  enddo
  !$omp end do

!!------------------------------------------------------
!! Y-axis
!!------------------------------------------------------
  !! (left boundary)
  !$omp do collapse(2) private(ii, yu, yd, ivar, ix, iz, amaxL, amaxR, FL, FR)
  do ix=1, subsub_ngrid
  do iz=1, subsub_ngrid
    ii = (iz-1)*subsub_ngrid2 + ix
    yu = ii + subsub_ngrid
    yd = (iz-1)*subsub_ngrid2 + ix + (subsub_ngrid-1)*subsub_ngrid

    amaxL = max(abs(subsub_hdummy(ix,1,iz,3,2)) + subsub_csarr(ix,1,iz), abs(subsub_hdummy(ix,subsub_ngrid,iz,3,2)) + subsub_csarr(ix,subsub_ngrid,iz))
    amaxR = max(abs(subsub_hdummy(ix,1,iz,3,2)) + subsub_csarr(ix,1,iz), abs(subsub_hdummy(ix,2,iz,3,2)) + subsub_csarr(ix,2,iz))

    do ivar=1,subsub_nhydro
      FL(ivar) = 0.5D0*(subsub_hdummy(ix,subsub_ngrid,iz,ivar,4) + subsub_hdummy(ix,1,iz,ivar,4)) - &
        0.5D0*amaxL*(subsub_hdummy(ix,1,iz,ivar,1) - subsub_hdummy(ix,subsub_ngrid,iz,ivar,1))
      FR(ivar) = 0.5D0*(subsub_hdummy(ix,1,iz,ivar,4) + subsub_hdummy(ix,2,iz,ivar,4)) - &
        0.5D0*amaxR*(subsub_hdummy(ix,2,iz,ivar,1) - subsub_hdummy(ix,1,iz,ivar,1))
    enddo

    
    do ivar=1,subsub_nhydro
      subsub_hydro(yd,ivar) = subsub_hydro(yd,ivar) - lam * FL(ivar)
      subsub_hydro(ii,ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
      subsub_hydro(yu,ivar) = subsub_hydro(yu,ivar) + lam * FR(ivar)
    enddo
  enddo
  enddo
  !$omp end do
 
  !! (interior) 
  do iy=3, subsub_ngrid-1, 2
    !$omp do collapse(2) private(ii, yu, yd, ix, iz, ivar, amaxL, amaxR, FL, FR)
    do iz=1, subsub_ngrid
    do ix=1, subsub_ngrid
      ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
      yu = ii + subsub_ngrid
      yd = ii - subsub_ngrid

      amaxL = max(abs(subsub_hdummy(ix,iy,iz,3,2)) + subsub_csarr(ix,iy,iz), abs(subsub_hdummy(ix,iy-1,iz,3,2)) + subsub_csarr(ix,iy-1,iz))
      amaxR = max(abs(subsub_hdummy(ix,iy,iz,3,2)) + subsub_csarr(ix,iy,iz), abs(subsub_hdummy(ix,iy+1,iz,3,2)) + subsub_csarr(ix,iy+1,iz))

      do ivar=1,subsub_nhydro
        FL(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,4) + subsub_hdummy(ix,iy-1,iz,ivar,4)) - &
          0.5D0*amaxL*(subsub_hdummy(ix,iy,iz,ivar,1) - subsub_hdummy(ix,iy-1,iz,ivar,1))
        FR(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,4) + subsub_hdummy(ix,iy+1,iz,ivar,4)) - &
          0.5D0*amaxR*(subsub_hdummy(ix,iy+1,iz,ivar,1) - subsub_hdummy(ix,iy,iz,ivar,1))
      enddo

      do ivar=1,subsub_nhydro
        subsub_hydro(yd,ivar) = subsub_hydro(yd,ivar) - lam * FL(ivar)
        subsub_hydro(ii,ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
        subsub_hydro(yu,ivar) = subsub_hydro(yu,ivar) + lam * FR(ivar)
      enddo
    enddo
    enddo
    !$omp end do
  enddo

  !! (right boundary)
  !$omp do collapse(2) private(ii, yu, yd, ivar, ix, iz, amaxL, amaxR, FL, FR)
  do iz=1, subsub_ngrid
  do ix=1, subsub_ngrid
    ii = (iz-1)*subsub_ngrid2 + (subsub_ngrid-1)*subsub_ngrid + ix
    yd = ii - subsub_ngrid
    yu = (iz-1)*subsub_ngrid2 + ix

    amaxL = max(abs(subsub_hdummy(ix,subsub_ngrid,iz,3,2)) + subsub_csarr(ix,subsub_ngrid,iz), abs(subsub_hdummy(ix,subsub_ngrid-1,iz,3,2)) + subsub_csarr(ix,subsub_ngrid-1,iz))
    amaxR = max(abs(subsub_hdummy(ix,subsub_ngrid,iz,3,2)) + subsub_csarr(ix,subsub_ngrid,iz), abs(subsub_hdummy(ix,1,iz,3,2)) + subsub_csarr(ix,1,iz))

    do ivar=1, 5
      FL(ivar) = 0.5D0*(subsub_hdummy(ix,subsub_ngrid,iz,ivar,4) + subsub_hdummy(ix,subsub_ngrid-1,iz,ivar,4)) - &
        0.5D0*amaxL*(subsub_hdummy(ix,subsub_ngrid,iz,ivar,1) - subsub_hdummy(ix,subsub_ngrid-1,iz,ivar,1))
      FR(ivar) = 0.5D0*(subsub_hdummy(ix,subsub_ngrid,iz,ivar,4) + subsub_hdummy(ix,1,iz,ivar,4)) - &
        0.5D0*amaxR*(subsub_hdummy(ix,1,iz,ivar,1) - subsub_hdummy(ix,subsub_ngrid,iz,ivar,1))
    enddo

    do ivar=1, 5
      subsub_hydro(ii,ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
      if(isodd)then
        subsub_hydro(yd,ivar) = subsub_hydro(yd,ivar) - lam*FL(ivar)
      endif
      subsub_hydro(yu,ivar) = subsub_hydro(yu,ivar) + lam*FR(ivar)
    enddo
  enddo
  enddo
  !$omp end do


!!------------------------------------------------------
!! Z-axis
!!------------------------------------------------------
  !! (left boundary)
  !$omp do collapse(2) private(ii, zu, zd, ivar, ix, iy, amaxL, amaxR, FL, FR)
  do iy=1, subsub_ngrid
  do ix=1, subsub_ngrid
    ii = (iy-1)*subsub_ngrid + ix
    zu = ii + subsub_ngrid2
    zd = (subsub_ngrid-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix

    amaxL = max(abs(subsub_hdummy(ix,iy,1,4,2)) + subsub_csarr(ix,iy,1), abs(subsub_hdummy(ix,iy,subsub_ngrid,4,2)) + subsub_csarr(ix,iy,subsub_ngrid))
    amaxR = max(abs(subsub_hdummy(ix,iy,1,4,2)) + subsub_csarr(ix,iy,1), abs(subsub_hdummy(ix,iy,2,4,2)) + subsub_csarr(ix,iy,2))

    do ivar=1,subsub_nhydro
      FL(ivar) = 0.5D0*(subsub_hdummy(ix,iy,subsub_ngrid,ivar,5) + subsub_hdummy(ix,iy,1,ivar,5)) - &
        0.5D0*amaxL*(subsub_hdummy(ix,iy,1,ivar,1) - subsub_hdummy(ix,iy,subsub_ngrid,ivar,1))
      FR(ivar) = 0.5D0*(subsub_hdummy(ix,iy,1,ivar,5) + subsub_hdummy(ix,iy,2,ivar,5)) - &
        0.5D0*amaxR*(subsub_hdummy(ix,iy,2,ivar,1) - subsub_hdummy(ix,iy,1,ivar,1))
    enddo

    
    do ivar=1,subsub_nhydro
      subsub_hydro(zd,ivar) = subsub_hydro(zd,ivar) - lam * FL(ivar)
      subsub_hydro(ii,ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
      subsub_hydro(zu,ivar) = subsub_hydro(zu,ivar) + lam * FR(ivar)
    enddo
  enddo
  enddo
  !$omp end do

  !! (interior) 
  do iz=3, subsub_ngrid-1, 2
    !$omp do collapse(2) private(ii, zu, zd, ivar, ix, iy, amaxL, amaxR, FL, FR)
    do iy=1, subsub_ngrid
    do ix=1, subsub_ngrid
      ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
      zu = ii + subsub_ngrid2
      zd = ii - subsub_ngrid2

      amaxL = max(abs(subsub_hdummy(ix,iy,iz,4,2)) + subsub_csarr(ix,iy,iz), abs(subsub_hdummy(ix,iy,iz-1,4,2)) + subsub_csarr(ix,iy,iz-1))
      amaxR = max(abs(subsub_hdummy(ix,iy,iz,4,2)) + subsub_csarr(ix,iy,iz), abs(subsub_hdummy(ix,iy,iz+1,4,2)) + subsub_csarr(ix,iy,iz+1))

      do ivar=1,subsub_nhydro
        FL(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,5) + subsub_hdummy(ix,iy,iz-1,ivar,5)) - &
          0.5D0*amaxL*(subsub_hdummy(ix,iy,iz,ivar,1) - subsub_hdummy(ix,iy,iz-1,ivar,1))
        FR(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,5) + subsub_hdummy(ix,iy,iz+1,ivar,5)) - &
          0.5D0*amaxR*(subsub_hdummy(ix,iy,iz+1,ivar,1) - subsub_hdummy(ix,iy,iz,ivar,1))
      enddo

      do ivar=1,subsub_nhydro
        subsub_hydro(zd,ivar) = subsub_hydro(zd,ivar) - lam * FL(ivar)
        subsub_hydro(ii,ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
        subsub_hydro(zu,ivar) = subsub_hydro(zu,ivar) + lam * FR(ivar)
      enddo
    enddo
    enddo
    !$omp end do
  enddo

  !! (right boundary)
  !$omp do collapse(2) private(ii, zu, zd, ivar, ix, iy, amaxL, amaxR, FL, FR)
  do iy=1, subsub_ngrid
  do ix=1, subsub_ngrid
    ii = (subsub_ngrid-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
    zd = ii - subsub_ngrid2
    zu = (iy-1)*subsub_ngrid + ix

    amaxL = max(abs(subsub_hdummy(ix,iy,subsub_ngrid,4,2)) + subsub_csarr(ix,iy,subsub_ngrid), abs(subsub_hdummy(ix,iy,subsub_ngrid-1,4,2)) + subsub_csarr(ix,iy,subsub_ngrid-1))
    amaxR = max(abs(subsub_hdummy(ix,iy,subsub_ngrid,4,2)) + subsub_csarr(ix,iy,subsub_ngrid), abs(subsub_hdummy(ix,iy,1,4,2)) + subsub_csarr(ix,iy,1))

    do ivar=1, 5
      FL(ivar) = 0.5D0*(subsub_hdummy(ix,iy,subsub_ngrid,ivar,5) + subsub_hdummy(ix,iy,subsub_ngrid-1,ivar,5)) - &
        0.5D0*amaxL*(subsub_hdummy(ix,iy,subsub_ngrid,ivar,1) - subsub_hdummy(ix,iy,subsub_ngrid-1,ivar,1))
      FR(ivar) = 0.5D0*(subsub_hdummy(ix,iy,subsub_ngrid,ivar,5) + subsub_hdummy(ix,iy,1,ivar,5)) - &
        0.5D0*amaxR*(subsub_hdummy(ix,iy,1,ivar,1) - subsub_hdummy(ix,iy,subsub_ngrid,ivar,1))
    enddo

    do ivar=1, 5
      subsub_hydro(ii,ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
      if(isodd)then
        subsub_hydro(zd,ivar) = subsub_hydro(zd,ivar) - lam*FL(ivar)
      endif
      subsub_hydro(zu,ivar) = subsub_hydro(zu,ivar) + lam*FR(ivar)
    enddo
  enddo
  enddo
  !$omp end do

end subroutine
!################################################################
!################################################################
!################################################################
!################################################################
