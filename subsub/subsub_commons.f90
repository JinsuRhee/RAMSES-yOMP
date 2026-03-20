module subsub_commons
  use amr_parameters
  use subsub_parameters

  integer :: subsub_nobj

  type subsub_type
     integer                               :: nlevel
     integer                               :: sink_ind
     real(dp)                              :: mass_tot
     real(dp)                              :: vxc, vyc, vzc
     integer                               :: clevel
     !real(dp), dimension(:,:), allocatable :: xg
     real(dp), dimension(:,:), allocatable :: vg
     real(dp), dimension(:,:), allocatable :: hydro
  end type subsub_type

  type(subsub_type), dimension(:), allocatable :: subsub_obj
  integer :: subsub_end = 0


end module subsub_commons