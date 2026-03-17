subroutine init_subsub
  use amr_commons
  use subsub_commons
  use mpi_mod
  implicit none
  !!-----
  !! This routine works to initialize the subsub variable to sink particles in myid
  !!-----

  !! Local Varaibles


  !!-----
  !! Check
  !!-----
  if(subsub_on .and. .not. sink) then
    if(myid .eq. 1) then
      call subsub_log('sink should be .true.', 'init_subsub')
    endif
  endif

end subroutine