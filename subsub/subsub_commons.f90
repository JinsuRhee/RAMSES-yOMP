module subsub_commons
  use amr_parameters
  use subsub_parameters

  real(dp)::subsub_testa=1.0D0
  real(dp)::subsub_testb=2.0D0
  real(dp)::subsub_testc=3.0D0

  

  type subsub_type
     integer                               :: nlevel
     integer                               :: sink_id
     real(dp), dimension(:,:), allocatable :: xg
     real(dp), dimension(:,:), allocatable :: vg
     real(dp), dimension(:,:), allocatable :: hydro
  end type subsub_type


end module subsub_commons