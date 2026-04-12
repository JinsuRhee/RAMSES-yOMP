module subsub_parameters
  use amr_parameters
  implicit none
  logical::subsub_on = .false. ! Flag for SUB2

  integer::subsub_ngrid = 4    ! # of cells in one direction
  integer::subsub_nhydro = 5   ! (rho, px, py, pz, E)
  
  !! I/O
  logical :: subsub_savecoarse = .false.
  integer :: subsub_mpidblpren = 3

  !! Mass in/out controller
  real(dp) :: subsub_massinfac = 2.0D0

  !! Density update
  real(dp) :: subsub_cfl = 0.5D0
  integer :: subsub_inflowtype = 1 ! 1 as spherical inflow
  real(dp) :: subsub_dfloor = 1.0D-10
  real(dp) :: subsub_pfloor = 1.0D-10

  !! Poisson update
  integer  :: subsub_poisson_type = 1 ! 1 as 6-Jac // 2 as CG // 3
  real(dp) :: subsub_poisson_softening = 0.0D0
  integer  :: subsub_poisson_niter = 1000
  real(dp) :: subsub_poisson_tolerance = 1.0D-4

  !! Hydro
  integer :: subsub_RiemannType = 1 ! 1 as Rusanov
  !! ETC
  real(dp) :: subsub_smallr = 1.0D-10
  
  !! DEV
  integer :: subsub_dev_icsphere = 0
  integer :: subsub_dev_gravonly = 1
  integer :: subsub_dev_hydroonly = 1

end module subsub_parameters