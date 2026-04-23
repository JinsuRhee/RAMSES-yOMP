module subsub_commons
  use amr_parameters
  use subsub_parameters

  type subsub_type
     integer                               :: sink_ind
     integer                               :: sink_id
     real(dp)                              :: sink_mass
     real(dp)                              :: mass_tot
     real(dp)                              :: mass_cell
     integer                               :: clevel
     integer                               :: domain
     real(dp), dimension(:,:), allocatable :: uold
        !! 0 myself
        !! 1 (x-) 2 (x+) 3(y-) 4(y+) 5(z-) 6(z+)
     !real(dp), dimension(:,:), allocatable :: xg
     !real(dp), dimension(:,:), allocatable :: vg !! ngrid^3 X ndim
     real(dp), dimension(:,:), allocatable :: hydro !! nlevel^3 X subsub_nhydro
     
     real(dp), dimension(:), allocatable :: phi !! ngird^3
     real(dp), dimension(:), allocatable :: phi_bh !! ngird^3

     real(dp), dimension(:,:), allocatable :: edgeBC
     real(dp), dimension(:,:,:,:,:), allocatable :: faceBC

     real(dp) :: v2, cs2
     logical :: new
     !! Not used at the moment
     !real(dp), dimension(:,:), allocatable :: fg !! ngrid^3 X ndim
  end type subsub_type

  type(subsub_type), dimension(:), allocatable :: subsub_obj    !! has the size of nisnkmax in each cpu
  integer :: subsub_end !! last index indicator

  integer :: subsub_nsink !! total number of sinks (should be equal to nsink if a sink is not created)
  integer :: subsub_mpidp
  integer :: subsub_nn
  integer :: subsub_nnface
  real(dp) :: subsub_boxlen
  real(dp) :: subsub_dx

  integer :: subsub_ngrid2

  real(dp), dimension(:), allocatable :: subsub_phi, subsub_phibh
  real(dp), dimension(:,:), allocatable :: subsub_fg
  real(dp), dimension(:), allocatable :: subsub_dd
  integer, dimension(:), allocatable :: subsub_faceind
  integer, dimension(:), allocatable :: subsub_faceindx
  integer, dimension(:), allocatable :: subsub_faceindy
  integer, dimension(:), allocatable :: subsub_faceindz

  real(dp), dimension(:,:), allocatable :: subsub_hydro
  real(dp), dimension(:,:,:,:,:), allocatable :: subsub_hdummy
  real(dp), dimension(:,:,:), allocatable :: subsub_csarr
  !real(dp), dimension(:,:,:,:), allocatable :: subsub_hfx, subsub_hfy, subsub_hfz
  !real(dp), dimension(:,:), allocatable :: subsub_vg_up, subsub_vg_down
  !real(dp), dimension(:,:), allocatable :: subsub_flux_up, subsub_flux_down

  real(dp), dimension(:), allocatable :: subsub_cgrhs
  real(dp), dimension(:), allocatable :: subsub_cgLphi
  real(dp), dimension(:), allocatable :: subsub_cgRes
  real(dp), dimension(:), allocatable :: subsub_cgp, subsub_cgLp
  real(dp) :: subsub_rhs2, subsub_rrold, subsub_rrnew, subsub_pLp, subsub_relres
  real(dp) :: subsub_alpha, subsub_beta
  logical :: subsub_skipcg

  real(dp), dimension(:,:,:,:,:), allocatable :: subsub_hydrobc

  integer :: subsub_debugn, subsub_debugtag
  real(dp) :: subsub_tcheck_cg(40), subsub_tcheck_cg_global(40)
  integer :: subsub_ncheck_cg(10), subsub_ncheck_cg_global(10)
  real(dp) :: subsub_tstart, subsub_tend, subsub_howlong

  real(dp), dimension(:,:,:), allocatable :: subsub_eBC, subsub_eBCall
  
  !! sink properties
  real(dp), dimension(:), allocatable :: subsub_v2sink, subsub_cs2sink

  !! BC related
  !integer :: subsub_ncloudmax
  
  !integer, dimension(:), allocatable :: subsub_clouds_ind
contains
  subroutine subsub_precision_mpi()
    use mpi_mod
    implicit none
    integer :: ierr
    

    if (dp == kind(1.0E0)) then
      subsub_mpidp = MPI_REAL
    else if (dp == kind(1.0D0)) then
      subsub_mpidp = MPI_DOUBLE_PRECISION
    endif
  end subroutine subsub_precision_mpi
end module subsub_commons