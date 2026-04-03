allocate(subsub_hdummy(1:subsub_ngrid, 1:subsub_ngrid, 1:subsub_ngrid, 1:subsub_nhydro, 1:5))

!1 conservative old
!2 primitive old
!3 Fx
!4 Fy
!5 Fz

!subsub_hfx (subsub_ngrid+1, subsub_ngrid, subsub_ngrid+1, subsub_nhydro)


subroutine subsub_drift(dt)
  use subsub_commons
  use subsub_parameters
  use hydro_parameters, ONLY:gamma
  implicit none

  real(dp) :: dt

  !! Local variables
  integer :: i, ivar
  integer :: ii, ix, iy, iz
  integer :: subsub_ngrid2
  !integer :: ix, iy, iz, ii, xu, xd, yu, yd, zu, zd
  !real(dp) :: deltarho
  !real(dp), dimension(1:2, 1:ndim) :: vgdummy, fluxdummy
  !real(dp), dimension(1:subsub_nhydro) :: Fx_p, Fx_m, Fy_p, Fy_m, Fz_p, Fz_m
  !real(dp), dimension(1:subsub_nhydro) :: Uc, Utmp, UX_p, UX_m, UY_p, UY_m, UZ_p, UZ_m

  !real(dp), dimension(1:subsub_nhydro) :: Uc, Un, Ub, Utmp
  !real(dp), dimension(1:subsub_nhydro) :: Fx_p, Fx_m, Fy_p, Fy_m, Fz_p, Fz_m
  !real(dp), dimension(1:subsub_nhydro) :: Qc, Qb
  real(dp) :: lam, ekin
  real(dp), dimension(1:subsub_nhydro) :: Utmp
  real(dp), dimension(1:2, 1:ndim, 1:subsub_nhydro, 1:5) :: hbc

  !$omp single
  subsub_ngrid2 = subsub_ngrid * subsub_ngrid
  !$omp end single

  !$omp single
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
  !$omp end single

  !lam = dt/subsub_dx
  !! copy to old arrays & save primitive variables
  !$omp do collapse(2) private(ekin, ii)
  do iz=1, subsub_ngrid
  do iy=1, subsub_ngrid
  do ix=1, subsub_ngrid
    ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
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
    subsub_hdummy(ix,iy,iz,5,2) = (gamma - 1.0D0) * (subsub_hdummy(ix,iy,iz,5,1) - ekin)
  enddo
  enddo
  enddo
  !$omp end do

  select case(subsub_RiemannType)
  case(1)
!!-----------------------------------------------------------------
!! Rusanov Solver
!!-----------------------------------------------------------------
  call subsub_hydroRiemann_Rusanov(dt, hbc)
  end select

  !$omp do private(Utmp)
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














subroutine subsub_hydroRiemann_Rusanov(dt, hbc)
  use subsub_commons
  use subsub_parameters
  use hydro_parameters, ONLY:gamma
  implicit none

  real(dp) :: dt
  real(dp), dimension(1:2, 1:ndim, 1:subsub_nhydro, 1:5) :: hbc
  !! Local variables
  integer :: ivar
  integer :: ix, iy, iz, subsub_ngrid2
  integer :: ii, xu, xd, yu, yd, zu, zd
  real(dp) :: amaxL, amaxR, lam
  !real(dp), dimension(1:subsub_nhydro) :: PL, PR

  !real(dp) :: csL, csR
  real(dp), dimension(1:subsub_nhydro) :: FL, FR
  real(dp), dimension(1:subsub_ngrid, 1:subsub_ngrid, 1:subsub_ngrid) :: csarr
  real(dp), dimension(1:2, 1:ndim) :: csarr_bc
  logical :: isodd
  !integer :: ivar

  !1 conservative old
  !2 primitive old
  !3 Fx
  !4 Fy
  !5 Fz

  !$omp single
  subsub_ngrid2 = subsub_ngrid * subsub_ngrid
  lam = dt / subsub_dx
  isodd = mod(subsub_ngrid,2) .eq. 1
  !$omp end single

  !!----- Compute sound speed first
  !$omp do collapse(2)
  do iz=1, subsub_ngrid
  do iy=1, subsub_ngrid
  do ix=1, subsub_ngrid
    csarr(ix, iy, iz) = sqrt(gamma * subsub_hdummy(ix, iy, iz, 5, 2) / max(subsub_hdummy(ix, iy, iz, 1, 2),subsub_dfloor))
  enddo
  enddo
  enddo
  !$omp end do

  !$omp single
  csarr_bc(1,1) = sqrt(gamma*hbc(1,1,5,2) / max(hbc(1,1,1,2),subsub_dfloor))
  csarr_bc(2,1) = sqrt(gamma*hbc(2,1,5,2) / max(hbc(2,1,1,2),subsub_dfloor))

  csarr_bc(1,2) = sqrt(gamma*hbc(1,2,5,2) / max(hbc(1,2,1,2),subsub_dfloor))
  csarr_bc(2,2) = sqrt(gamma*hbc(2,2,5,2) / max(hbc(2,2,1,2),subsub_dfloor))

  csarr_bc(1,3) = sqrt(gamma*hbc(1,3,5,2) / max(hbc(1,3,1,2),subsub_dfloor))
  csarr_bc(2,3) = sqrt(gamma*hbc(2,3,5,2) / max(hbc(2,3,1,2),subsub_dfloor))
  !$omp end single

  !!----- Compute Flux
  !$omp do collapse(2)
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

  !$omp single
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
  !$omp end single

  
  !! Solver
!!------------------------------------------------------
!! X-axis
!!------------------------------------------------------
  !! (left boundary)
  !$omp do collapse(2) private(ii, xu, amaxL, amaxR, FL, FR)
  do iy=1, subsub_ngrid
  do iz=1, subsub_ngrid
    ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + 1
    xu = ii + 1

    amaxL = max(abs(subsub_hdummy(1,iy,iz,2,2)) + csarr(1,iy,iz), abs(hbc(1,1,2,2)) + csarr_bc(1,1))
    amaxR = max(abs(subsub_hdummy(1,iy,iz,2,2)) + csarr(1,iy,iz), abs(subsub_hdummy(2,iy,iz,2,2)) + csarr(2,iy,iz))

    do ivar=1,5
      FL(ivar) = 0.5D0*(hbc(1,1,ivar,3) + subsub_hdummy(1,iy,iz,ivar,3)) - &
        0.5D0*amaxL*(subsub_hdummy(1,iy,iz,ivar,1) - hbc(1,1,ivar,1))
      FR(ivar) = 0.5D0*(subsub_hdummy(1,iy,iz,ivar,3) + subsub_hdummy(2,iy,iz,ivar,3)) - &
        0.5D0*amaxR*(subsub_hdummy(2,iy,iz,ivar,1) - subsub_hdummy(1,iy,iz,ivar,1))
    enddo

    
    do ivar=1,5
      subsub_hydro(ii,ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
      subsub_hydro(xu,ivar) = subsub_hydro(xu,ivar) + lam * FR(ivar)
    enddo
  enddo
  enddo
  !$omp end do
 
  !! (interior) 
  do ix=3, subsub_ngrid-1, 2
    !$omp do collapse(2) private(ii, xu, xd, amaxL, amaxR, FL, FR)
    do iy=1, subsub_ngrid
    do iz=1, subsub_ngrid
      ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
      xu = ii + 1
      xd = ii - 1

      amaxL = max(abs(subsub_hdummy(ix,iy,iz,2,2)) + csarr(ix,iy,iz), abs(subsub_hdummy(ix-1,iy,iz,2,2)) + csarr(ix-1,iy,iz))
      amaxR = max(abs(subsub_hdummy(ix,iy,iz,2,2)) + csarr(ix,iy,iz), abs(subsub_hdummy(ix+1,iy,iz,2,2)) + csarr(ix+1,iy,iz))

      do ivar=1,5
        FL(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,3) + subsub_hdummy(ix-1,iy,iz,ivar,3)) - &
          0.5D0*amaxL*(subsub_hdummy(ix,iy,iz,ivar,1) - subsub_hdummy(ix-1,iy,iz,ivar,1))
        FR(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,3) + subsub_hdummy(ix+1,iy,iz,ivar,3)) - &
          0.5D0*amaxR*(subsub_hdummy(ix+1,iy,iz,ivar,1) - subsub_hdummy(ix,iy,iz,ivar,1))
      enddo

      do ivar=1,5
        subsub_hydro(xd,ivar) = subsub_hydro(xd,ivar) - lam * FL(ivar)
        subsub_hydro(ii,ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
        subsub_hydro(xu,ivar) = subsub_hydro(xu,ivar) + lam * FR(ivar)
      enddo
    enddo
    enddo
    !$omp end do
  enddo

  !! (right boundary)
  !$omp do collapse(2) private(ii, xd, amaxL, amaxR, FL, FR)
  do iy=1, subsub_ngrid
  do iz=1, subsub_ngrid
    ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + subsub_ngrid
    xd = ii - 1

    amaxL = max(abs(subsub_hdummy(subsub_ngrid,iy,iz,2,2)) + csarr(subsub_ngrid,iy,iz), abs(subsub_hdummy(subsub_ngrid-1,iy,iz,2,2)) + csarr(subsub_ngrid-1,iy,iz))
    amaxR = max(abs(subsub_hdummy(subsub_ngrid,iy,iz,2,2)) + csarr(subsub_ngrid,iy,iz), abs(hbc(2,1,2,2)) + csarr_bc(2,1))

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
  !$omp do collapse(2) private(ii, yu, amaxL, amaxR, FL, FR)
  do ix=1, subsub_ngrid
  do iz=1, subsub_ngrid
    ii = (iz-1)*subsub_ngrid2 + ix
    yu = ii + subsub_ngrid

    amaxL = max(abs(subsub_hdummy(ix,1,iz,2,2)) + csarr(ix,1,iz), abs(hbc(1,2,2,2)) + csarr_bc(1,2))
    amaxR = max(abs(subsub_hdummy(ix,1,iz,2,2)) + csarr(ix,1,iz), abs(subsub_hdummy(ix,2,iz,2,2)) + csarr(ix,2,iz))

    do ivar=1,5
      FL(ivar) = 0.5D0*(hbc(1,2,ivar,3) + subsub_hdummy(ix,1,iz,ivar,3)) - &
        0.5D0*amaxL*(subsub_hdummy(ix,1,iz,ivar,1) - hbc(1,2,ivar,1))
      FR(ivar) = 0.5D0*(subsub_hdummy(ix,1,iz,ivar,3) + subsub_hdummy(ix,2,iz,ivar,3)) - &
        0.5D0*amaxR*(subsub_hdummy(ix,2,iz,ivar,1) - subsub_hdummy(ix,1,iz,ivar,1))
    enddo

    
    do ivar=1,5
      subsub_hydro(ii,ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
      subsub_hydro(yu,ivar) = subsub_hydro(yu,ivar) + lam * FR(ivar)
    enddo
  enddo
  enddo
  !$omp end do
 
  !! (interior) 
  do iy=3, subsub_ngrid-1, 2
    !$omp do collapse(2) private(ii, yu, yd, amaxL, amaxR, FL, FR)
    do iz=1, subsub_ngrid
    do ix=1, subsub_ngrid
      ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
      yu = ii + subsub_ngrid
      yd = ii - subsub_ngrid

      amaxL = max(abs(subsub_hdummy(ix,iy,iz,2,2)) + csarr(ix,iy,iz), abs(subsub_hdummy(ix,iy-1,iz,2,2)) + csarr(ix,iy-1,iz))
      amaxR = max(abs(subsub_hdummy(ix,iy,iz,2,2)) + csarr(ix,iy,iz), abs(subsub_hdummy(ix,iy+1,iz,2,2)) + csarr(ix,iy+1,iz))

      do ivar=1,5
        FL(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,3) + subsub_hdummy(ix,iy-1,iz,ivar,3)) - &
          0.5D0*amaxL*(subsub_hdummy(ix,iy,iz,ivar,1) - subsub_hdummy(ix,iy-1,iz,ivar,1))
        FR(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,3) + subsub_hdummy(ix,iy+1,iz,ivar,3)) - &
          0.5D0*amaxR*(subsub_hdummy(ix,iy+1,iz,ivar,1) - subsub_hdummy(ix,iy,iz,ivar,1))
      enddo

      do ivar=1,5
        subsub_hydro(yd,ivar) = subsub_hydro(yd,ivar) - lam * FL(ivar)
        subsub_hydro(ii,ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
        subsub_hydro(yu,ivar) = subsub_hydro(yu,ivar) + lam * FR(ivar)
      enddo
    enddo
    enddo
    !$omp end do
  enddo

  !! (right boundary)
  !$omp do collapse(2) private(ii, yd, amaxL, amaxR, FL, FR)
  do iz=1, subsub_ngrid
  do ix=1, subsub_ngrid
    ii = (iz-1)*subsub_ngrid2 + (subsub_ngrid-1)*subsub_ngrid + ix
    yd = ii - subsub_ngrid

    amaxL = max(abs(subsub_hdummy(ix,subsub_ngrid,iz,2,2)) + csarr(ix,subsub_ngrid,iz), abs(subsub_hdummy(ix,subsub_ngrid-1,iz,2,2)) + csarr(ix,subsub_ngrid-1,iz))
    amaxR = max(abs(subsub_hdummy(ix,subsub_ngrid,iz,2,2)) + csarr(ix,subsub_ngrid,iz), abs(hbc(2,2,2,2)) + csarr_bc(2,2))

    do ivar=1, 5
      FL(ivar) = 0.5D0*(subsub_hdummy(ix,subsub_ngrid,iz,ivar,3) + subsub_hdummy(ix,subsub_ngrid-1,iz,ivar,3)) - &
        0.5D0*amaxL*(subsub_hdummy(ix,subsub_ngrid,iz,ivar,1) - subsub_hdummy(ix,subsub_ngrid-1,iz,ivar,1))
      FR(ivar) = 0.5D0*(subsub_hdummy(ix,subsub_ngrid,iz,ivar,3) + hbc(2,2,ivar,3)) - &
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
  !$omp do collapse(2) private(ii, zu, amaxL, amaxR, FL, FR)
  do iy=1, subsub_ngrid
  do ix=1, subsub_ngrid
    ii = (iy-1)*subsub_ngrid + ix
    zu = ii + subsub_ngrid2

    amaxL = max(abs(subsub_hdummy(ix,iy,1,2,2)) + csarr(ix,iy,1), abs(hbc(1,3,2,2)) + csarr_bc(1,3))
    amaxR = max(abs(subsub_hdummy(ix,iy,1,2,2)) + csarr(ix,iy,1), abs(subsub_hdummy(ix,iy,2,2,2)) + csarr(ix,iy,2))

    do ivar=1,5
      FL(ivar) = 0.5D0*(hbc(1,3,ivar,3) + subsub_hdummy(ix,iy,1,ivar,3)) - &
        0.5D0*amaxL*(subsub_hdummy(ix,iy,1,ivar,1) - hbc(1,3,ivar,1))
      FR(ivar) = 0.5D0*(subsub_hdummy(ix,iy,1,ivar,3) + subsub_hdummy(ix,iy,2,ivar,3)) - &
        0.5D0*amaxR*(subsub_hdummy(ix,iy,2,ivar,1) - subsub_hdummy(ix,iy,1,ivar,1))
    enddo

    
    do ivar=1,5
      subsub_hydro(ii,ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
      subsub_hydro(zu,ivar) = subsub_hydro(zu,ivar) + lam * FR(ivar)
    enddo
  enddo
  enddo
  !$omp end do
 
  !! (interior) 
  do iz=3, subsub_ngrid-1, 2
    !$omp do collapse(2) private(ii, zu, zd, amaxL, amaxR, FL, FR)
    do iy=1, subsub_ngrid
    do ix=1, subsub_ngrid
      ii = (iz-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
      zu = ii + subsub_ngrid2
      zd = ii - subsub_ngrid2

      amaxL = max(abs(subsub_hdummy(ix,iy,iz,2,2)) + csarr(ix,iy,iz), abs(subsub_hdummy(ix,iy,iz-1,2,2)) + csarr(ix,iy,iz-1))
      amaxR = max(abs(subsub_hdummy(ix,iy,iz,2,2)) + csarr(ix,iy,iz), abs(subsub_hdummy(ix,iy,iz+1,2,2)) + csarr(ix,iy,iz+1))

      do ivar=1,5
        FL(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,3) + subsub_hdummy(ix,iy,iz-1,ivar,3)) - &
          0.5D0*amaxL*(subsub_hdummy(ix,iy,iz,ivar,1) - subsub_hdummy(ix,iy,iz-1,ivar,1))
        FR(ivar) = 0.5D0*(subsub_hdummy(ix,iy,iz,ivar,3) + subsub_hdummy(ix,iy,iz+1,ivar,3)) - &
          0.5D0*amaxR*(subsub_hdummy(ix,iy,iz+1,ivar,1) - subsub_hdummy(ix,iy,iz,ivar,1))
      enddo

      do ivar=1,5
        subsub_hydro(zd,ivar) = subsub_hydro(zd,ivar) - lam * FL(ivar)
        subsub_hydro(ii,ivar) = subsub_hydro(ii,ivar) - lam * (FR(ivar) - FL(ivar))
        subsub_hydro(zu,ivar) = subsub_hydro(zu,ivar) + lam * FR(ivar)
      enddo
    enddo
    enddo
    !$omp end do
  enddo

  !! (right boundary)
  !$omp do collapse(2) private(ii, zd, amaxL, amaxR, FL, FR)
  do iy=1, subsub_ngrid
  do ix=1, subsub_ngrid
    ii = (subsub_ngrid-1)*subsub_ngrid2 + (iy-1)*subsub_ngrid + ix
    zd = ii - subsub_ngrid2

    amaxL = max(abs(subsub_hdummy(ix,iy,subsub_ngrid,2,2)) + csarr(ix,iy,subsub_ngrid), abs(subsub_hdummy(ix,iy,subsub_ngrid-1,2,2)) + csarr(ix,iy,subsub_ngrid-1))
    amaxR = max(abs(subsub_hdummy(ix,iy,subsub_ngrid,2,2)) + csarr(ix,iy,subsub_ngrid), abs(hbc(2,3,2,2)) + csarr_bc(2,3))

    do ivar=1, 5
      FL(ivar) = 0.5D0*(subsub_hdummy(ix,iy,subsub_ngrid,ivar,3) + subsub_hdummy(ix,iy,subsub_ngrid-1,ivar,3)) - &
        0.5D0*amaxL*(subsub_hdummy(ix,iy,subsub_ngrid,ivar,1) - subsub_hdummy(ix,iy,subsub_ngrid-1,ivar,1))
      FR(ivar) = 0.5D0*(subsub_hdummy(ix,iy,subsub_ngrid,ivar,3) + hbc(2,3,ivar,3)) - &
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