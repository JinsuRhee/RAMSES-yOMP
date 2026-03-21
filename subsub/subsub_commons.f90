module subsub_commons
  use amr_parameters
  use subsub_parameters

  type subsub_type
     integer                               :: sink_ind
     integer                               :: sink_id
     real(dp)                              :: mass_tot
     real(dp)                              :: vxc, vyc, vzc
     integer                               :: clevel
     !real(dp), dimension(:,:), allocatable :: xg
     real(dp), dimension(:,:), allocatable :: vg !! nlevel^3 X ndim
     real(dp), dimension(:,:), allocatable :: hydro !! nlevel^3 X subsub_nhydro
  end type subsub_type

  type(subsub_type), dimension(:), allocatable :: subsub_obj    !! has the size of nisnkmax in each cpu
  integer :: subsub_end !! last index indicator

  integer :: subsub_nsink !! total number of sinks (should be equal to nsink if a sink is not created)
  integer :: subsub_mpidp

contains
  subroutine subsub_precision_mpi()
    integer :: ierr

    if (dp == kind(1.0E0)) then
      subsub_mpidp = MPI_REAL
    else if (dp == kind(1.0D0)) then
      subsub_mpidp = MPI_DOUBLE_PRECISION
    endif
  end subroutine subsub_precision_mpi
end module subsub_commons