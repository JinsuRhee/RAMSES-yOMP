subroutine read_subsub_params(nml_ok)
  use subsub_commons
  use mpi_mod
  implicit none
  logical::nml_ok
  !--------------------------------------------------
  ! Local variables
  !--------------------------------------------------

#ifdef SOLVERmhd
#endif


#ifdef SOLVERmhd

#if NVAR>8+NENER

#endif
#else
#if NVAR>NDIM+2+NENER

#endif
#endif
#if NENER>0

#endif

  namelist/subsub_params/subsub_on &
       & ,subsub_ngrid &
       & ,subsub_savecoarse &
       & ,subsub_poisson_type &
       & ,subsub_poisson_niter &
       & ,subsub_poisson_tolerance &
       & , subsub_poisson_softening


  ! Read namelist file
  rewind(1)
  read(1,NML=subsub_params,END=101)
101 continue

end subroutine
       