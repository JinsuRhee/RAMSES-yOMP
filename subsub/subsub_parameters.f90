module subsub_parameters
  use amr_parameters
  implicit none
  logical::subsub_on = .false. ! Flag for SUB2

  integer::subsub_ngrid = 4    ! # of cells in one direction
  integer::subsub_nhydro = 1

  !! I/O
  logical :: subsub_savecoarse = .false.

  !! Density update
  real(dp) :: subsub_cfl = 0.5D0
  integer :: subsub_inflowtype = 1 ! 1 as spherical inflow
  real(dp) :: subsub_densityfloor = 1.0D-10

  !! Poission update
  real(dp) :: subsub_softening = 0.0D0
  integer :: subsub_poission_niter = 1000
  real(dp) ::subsub_poissiontolerance = 1.0D-10

  !! ETC
  real(dp) :: subsub_smallr = 1.0D-10
  

end module subsub_parameters