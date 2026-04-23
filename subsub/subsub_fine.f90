!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_fine(ilevel)
  !!-----
  !! This routine creates/removes a new subsub obj based on the update sink list
  !!-----
  use pm_commons
  use amr_commons
  use subsub_parameters
  use subsub_commons
  use mpi_mod
  implicit none

  integer :: ilevel

  !! Local variables
  integer :: i, j, k, subsub_sinkinmyid, info, subsub_nsink_old
  integer :: ind, i0, id0, iend, ismy, isokay
  integer, dimension(1:nsink) :: sink_hash, sink_hash_next, sink_ismatch, sink_ismatch_all
  integer, dimension(1:nsink) :: ismysink, ismysink_dump
  integer, dimension(1:nsink) :: nowmysink, nowmysink_dump
  integer, dimension(:), allocatable :: subsub_idmatch
  integer :: nx_loc
  real(dp):: scale
  character(LEN=5) :: nchar
  character(LEN=80) :: filename, dirname
  real(dp) :: ttsta, ttend, tcheck(10)

  if(nsink .eq. 0) return
  if(.not. sink) return
  if(ilevel < nlevelmax_current) return


  if(myid.eq.1) then
    subsub_tstart=MPI_WTIME()
  endif

  !!----- Update Domain
  !! This should be ahead of subsub_compute
  !! Initialize cells if it is newly formed
  call subsub_updatedomain(ilevel)

  tcheck(1) = MPI_WTIME()

  !!----- Update Sink property arrays
  call subsub_sinkprop


  !!----- Compute Boundary conditions using cloud particles
  call subsub_edgeBC

  tcheck(2) = MPI_WTIME()

  !write(*,*) 'edge done = ', myid, tcheck(2) - tcheck(1)

  !call MPI_BARRIER(MPI_COMM_WORLD, i0)
  tcheck(3) = MPI_WTIME()
  call subsub_faceBC

!! velocity to be zero
if(subsub_end .gt. 0)then
  do i=1,subsub_end
    do j=1, subsub_ngrid
      do k=1, subsub_ngrid
        subsub_obj(i)%faceBC(:,:,j,k,2) = 0.0D0
        subsub_obj(i)%faceBC(:,:,j,k,3) = 0.0D0
        subsub_obj(i)%faceBC(:,:,j,k,4) = 0.0D0
      enddo
    enddo
  enddo
  
!    if(subsub_obj(i)%sink_ind .ne. 3) cycle
!      write(*,*) subsub_obj(i)%edgeBC(1,1)
!      write(*,*) subsub_obj(i)%edgeBC(5,1)
!      write(*,*) subsub_obj(i)%edgeBC(3,1)
!      write(*,*) subsub_obj(i)%edgeBC(7,1)
!
!      write(*,*) maxval(subsub_obj(i)%faceBC(1,1,:,:,1)), minval(subsub_obj(i)%faceBC(1,1,:,:,1))
!      write(*,*) maxval(subsub_obj(i)%faceBC(1,1,:,:,2)), minval(subsub_obj(i)%faceBC(1,1,:,:,2))
!      write(*,*) maxval(subsub_obj(i)%faceBC(1,1,:,:,3)), minval(subsub_obj(i)%faceBC(1,1,:,:,3))
!      write(*,*) maxval(subsub_obj(i)%faceBC(1,1,:,:,4)), minval(subsub_obj(i)%faceBC(1,1,:,:,4))
!      write(*,*) maxval(subsub_obj(i)%faceBC(1,1,:,:,5)), minval(subsub_obj(i)%faceBC(1,1,:,:,5))
!    
!  enddo
endif
tcheck(4) = MPI_WTIME()
  !tcheck(4) = MPI_WTIME()

  !write(*,*) 'face done = ', myid, tcheck(4) - tcheck(3)
  !call MPI_BARRIER(MPI_COMM_WORLD, i0)



!  call subsub_test
  !!----- Compute BC
  !call subsub_computeBC(ilevel)
  !! Make BC
  

!do i=1, subsub_end
!if(subsub_obj(i)%sink_id .eq. 3)then
!  
!  write(*,*)'     niter        = ', subsub_debugtag
!  
!  write(*,*)'     maxrho       = ', maxval(subsub_obj(i)%hydro(:,1))
!  write(*,*)'     minrho       = ', minval(subsub_obj(i)%hydro(:,1))
!
!  write(*,*)'     maxE       = ', maxval(subsub_obj(i)%hydro(:,5))
!  write(*,*)'     minE       = ', minval(subsub_obj(i)%hydro(:,5))
!
!  write(*,*)'     maxPx       = ', maxval(subsub_obj(i)%hydro(:,2))
!  write(*,*)'     minPx       = ', minval(subsub_obj(i)%hydro(:,2))
!
!  write(*,*)'     maxPy       = ', maxval(subsub_obj(i)%hydro(:,3))
!  write(*,*)'     minPy       = ', minval(subsub_obj(i)%hydro(:,3))
!
!  write(*,*)'     maxPz       = ', maxval(subsub_obj(i)%hydro(:,4))
!  write(*,*)'     minPz       = ', minval(subsub_obj(i)%hydro(:,4))
!
!  subsub_debugtag = 0
!  subsub_tcheck_cg_global(1) = 0.0D0
!  !write(*,*)'     fact2       = ', fact2
!endif
!enddo

  !!----- Update Property
tcheck(5) = MPI_WTIME()
  call subsub_compute(ilevel)
tcheck(6) = MPI_WTIME()
  !write(*,*) 'me here : ', myid


!CALL MPI_BARRIER(MPI_COMM_WORLD, i0)
!
!do i=1, ncpu
!  if(myid .eq. i) then
!    write(*,*) ' ============ '
!    write(*,*) '    myid = ', myid
!    if(subsub_end .gt. 0) then
!      do j=1, subsub_end
!        write(*,*) '        ', subsub_obj(j)%sink_id
!      enddo
!
!      write(*,*) '    Total T = ', tcheck(6)-tcheck(1)
!      write(*,*) '      eBC ', tcheck(2)-tcheck(1)
!      write(*,*) '      fBC ', tcheck(4)-tcheck(3)
!      write(*,*) '      com ', tcheck(6)-tcheck(5), ' with ', subsub_ncheck_cg(1), ' interations'
!      write(*,*) '   '
!    endif
!    do j=1, 100000
!    enddo
!  endif
!  CALL MPI_BARRIER(MPI_COMM_WORLD, i0)
!enddo  

!do i=1, subsub_end
!if(subsub_obj(i)%sink_id .eq. 3)then
!  
!  write(*,*)'     niter        = ', subsub_debugtag
!  
!  write(*,*)'     maxrho       = ', maxval(subsub_obj(i)%hydro(:,1))
!  write(*,*)'     minrho       = ', minval(subsub_obj(i)%hydro(:,1))
!
!  write(*,*)'     maxE       = ', maxval(subsub_obj(i)%hydro(:,5))
!  write(*,*)'     minE       = ', minval(subsub_obj(i)%hydro(:,5))
!
!  write(*,*)'     maxPx       = ', maxval(subsub_obj(i)%hydro(:,2))
!  write(*,*)'     minPx       = ', minval(subsub_obj(i)%hydro(:,2))
!
!  write(*,*)'     maxPy       = ', maxval(subsub_obj(i)%hydro(:,3))
!  write(*,*)'     minPy       = ', minval(subsub_obj(i)%hydro(:,3))
!
!  write(*,*)'     maxPz       = ', maxval(subsub_obj(i)%hydro(:,4))
!  write(*,*)'     minPz       = ', minval(subsub_obj(i)%hydro(:,4))
!
!  subsub_debugtag = 0
!  subsub_tcheck_cg_global(1) = 0.0D0
!  !write(*,*)'     fact2       = ', fact2  
!endif
!enddo


  subsub_debugtag = subsub_debugtag + subsub_ncheck_cg(1)
  subsub_tcheck_cg_global(1) = subsub_tcheck_cg_global(1) + dtnew(ilevel)
  if(myid.eq.1)then
    subsub_tend=MPI_WTIME()
    subsub_howlong = subsub_howlong + subsub_tend - subsub_tstart
  endif

  call MPI_ALLREDUCE(subsub_ncheck_cg(1), subsub_ncheck_cg(2), 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, info)

!! For runtime investigation
!  call MPI_ALLREDUCE(subsub_tcheck_cg(1), subsub_tcheck_cg_global(1), 1, subsub_mpidp, MPI_MAX, MPI_COMM_WORLD, info)
!  call MPI_ALLREDUCE(subsub_tcheck_cg(2), subsub_tcheck_cg_global(2), 1, subsub_mpidp, MPI_MAX, MPI_COMM_WORLD, info)
!  call MPI_ALLREDUCE(subsub_tcheck_cg(3), subsub_tcheck_cg_global(3), 1, subsub_mpidp, MPI_MAX, MPI_COMM_WORLD, info)
!  call MPI_ALLREDUCE(subsub_tcheck_cg(4), subsub_tcheck_cg_global(4), 1, subsub_mpidp, MPI_MAX, MPI_COMM_WORLD, info)
!  call MPI_ALLREDUCE(subsub_tcheck_cg(5), subsub_tcheck_cg_global(5), 1, subsub_mpidp, MPI_MAX, MPI_COMM_WORLD, info)
!
!  call MPI_ALLREDUCE(subsub_ncheck_cg(1), subsub_ncheck_cg_global(1), 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, info)
!  call MPI_ALLREDUCE(subsub_ncheck_cg(2), subsub_ncheck_cg_global(2), 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, info)
!  call MPI_ALLREDUCE(subsub_ncheck_cg(3), subsub_ncheck_cg_global(3), 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, info)

!if(myid.eq.1) then
!  write(*,*)'     poisson      = ', subsub_tcheck_cg_global(1)
!  write(*,*)'     force        = ', subsub_tcheck_cg_global(2)
!  write(*,*)'     dt           = ', subsub_tcheck_cg_global(3)
!  write(*,*)'     kick         = ', subsub_tcheck_cg_global(4)
!  write(*,*)'     drift        = ', subsub_tcheck_cg_global(5)
!  write(*,*)'     poisson_nall = ', subsub_ncheck_cg_global(1)
!  write(*,*)'     poisson_n1   = ', subsub_ncheck_cg_global(2)
!  write(*,*)'     poisson_n2   = ', subsub_ncheck_cg_global(3)
!endif



end subroutine subsub_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_compute(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use subsub_commons
  use subsub_parameters
  use mpi_mod
  implicit none

  integer :: ilevel
  !! Local variables
  integer :: i, j
  !real(dp):: scale, mtot_old, mtot_new, mbh_old, mbh_new
  !integer :: ind_cell, ind_grid, ind_level, nx_loc, sinkind
  logical :: subflag
  !real(dp):: dx, tff, tffnew, threepi2, fourpi, vffnew
  !real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v


  if(subsub_end .eq. 0) return

  !! constants
  !nx_loc=(icoarse_max-icoarse_min+1)
  !scale=boxlen/dble(nx_loc)

  !!----- RHEE -----
  !! If the memory usage becomes problematic, we may think of reclying arrays below?
  !!----------------
  
  !1 conservative old
  !2 primitive old
  !3 Fx
  !4 Fy
  !5 Fz

  !! update sink by sink
  do i=1, subsub_end

    
!write(*,*) '%112233 beff: ', myid, i, subsub_obj(i)%mass_cell, subsub_obj(i)%mass_tot
    !! update mass_tot
    !mtot_old = subsub_obj(i)%mass_tot
    !mbh_old = subsub_obj(i)%sink_mass
    !mtot_new = (subsub_boxlen**ndim) * max(subsub_obj(i)%uold(0,1), subsub_smallr)
    !mbh_new = msink(sinkind)

   !! update density by the mass change
   !!----- RHEE ------
   !! - Delta M = mass_inflow - mass_outflow
   !! - mass_in:
   !! 		spherical inflow / rotating inflow / anisotropic inflow
   !! - mass_out:
   !! 		by jet or outflow?
   !!		by tidal mass loss?
   !!
   !! - At the moment,  if Delta M > 0: all by spherical inflow (v_in || vec(r) )
   !! 					if Delta M < 0: all by spherical outflow (v_out || vec(r) )
   !! BC as the mass flux
   !! Mass flux as the uniform mass in/outflow
   !! Infall velocity of v_ff
!   if(idsink(sinkind).eq.829) then
!     subsub_obj(i)%hydro(:,1) = mtot_old / subsub_boxlen**3
!     !write(*,*) '%456456 rho_min = ', minval(subsub_obj(i)%hydro(:,1))
!     !write(*,*) '%456456 rho_max = ', maxval(subsub_obj(i)%hydro(:,1))
!     !write(*,*) '%456456 new rho = ', uold(ind_cell,1)
!     !write(*,*) '%456456 BH mass = ', msink(subsub_obj(i)%sink_ind) 
!     !write(*,*) '%456456 BH mass = ', msink(subsub_obj(i)%sink_ind) * scale_d*scale_l**3 / 2d33 / 1e6
!
!     call subsub_iget3ind(j, 2, 2, 2)
!     write(*,*) '%456456 den(2,2,2)', subsub_obj(i)%hydro(j,1)
!     write(*,*) '%456456 den(1,1,1)', subsub_obj(i)%hydro(1,1)
!   endif

   

   call subsub_computefine(i, ilevel)!, mtot_new, mbh_new)

   !! update some properties
    subsub_obj(i)%mass_cell = (subsub_boxlen**ndim) * max(subsub_obj(i)%uold(0,1),subsub_dfloor)
    subsub_obj(i)%sink_mass = msink(subsub_obj(i)%sink_ind)

!if(mtot_new-mtot_old .gt. 0) write(*,*) 'good sink = ', idsink(sinkind)
!   if(idsink(sinkind) .eq. 829) then
!     write(*,*) '%456456 mass_old = ', mtot_old
!     write(*,*) '%456456 mass_new = ', mtot_new
!     write(*,*) '%456456 dm =', mtot_new - mtot_old
!     write(*,*) '%456456 dmbh =', mbh_new - mbh_old
!     write(*,*) '%456546 dmtot =', mtot_new - mtot_old + mbh_new - mbh_old
!     write(*,*) '%456456 by ct = ', sum(subsub_obj(i)%hydro(:,1)) * (subsub_boxlen/subsub_ngrid)**ndim
!     call subsub_iget3ind(j, 2, 2, 2)
!     write(*,*) '%456456 den(2,2,2)', subsub_obj(i)%hydro(j,1)
!     write(*,*) '%456456 den(1,1,1)', subsub_obj(i)%hydro(1,1)
!     !stop
!   endif
   !rho_flux = ((mtot_new - mtot_old) / dx**2) / 6.0D0 !! Mass/area
   !rho_flux = rho_flux / (dx / real(subsub_level,kind=dp))**3.0D0 !! from ghost cells (rho/area)

!! neglect mass conversion (starts from initial mass)
   !subsub_obj(i)%mass_tot = mtot_new
   !subsub_obj(i)%sink_mass = msink(sinkind)
!write(*,*) '%112233 done ', myid, i, ' / ', subsub_end
!write(*,*) '%112233 done: ', myid, i, subsub_obj(i)%sink_id, subsub_obj(i)%mass_cell, subsub_obj(i)%mass_tot
  enddo

end subroutine subsub_compute
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_computeBC(sink_ind, edgeBC)

  !!-----
  !! This routine finds a cell which (x, y, z) is located at
  !! (x,y,z) should be in a code unit
  !!-----
  use amr_commons
  use subsub_commons
  use pm_parameters
  use pm_commons
  implicit none
  
  integer :: sink_ind
  real(dp), dimension(1:twotondim, 1:subsub_nhydro) :: edgeBC

  ! Local variables
  integer :: i
  integer :: ncloud_max, ncloud
  integer :: sinkid
  integer :: nx_loc
  real(dp), dimension(:,:), allocatable :: clouds
  real(dp) :: dx_min, rmax, scale
  real(dp) :: x, y, z

  if(ndim.ne.3)then
    write(*,*) 'ndim = 3 is implemented only'
    stop
  endif
return

  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)

  !!----- Allocate clouds
  ncloud_max = (2*ir_cloud+1)**ndim
  allocate(clouds(1:ncloud_max, 1:subsub_nhydro+ndim))


  !! Retrieve cloud particle properties with collected cells
  x = xsink(sink_ind,1)/scale
  y = xsink(sink_ind,2)/scale
  z = xsink(sink_ind,3)/scale

  call subsub_getcloud(idsink(sink_ind), x, y, z, clouds, ncloud, ncloud_max)
    
  !! Interpolate to 8 edges
  !call subsub_edgeBC(edgeBC, clouds, ncloud, ncloud_max)


  deallocate(clouds)
end subroutine subsub_computeBC
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_edgeBC!(edgeBC, clouds, ncloud, ncloud_max)
  !!-----
  !! This routine computes contributions of cloud particles in each domain to the corresponding edge BCs
  !!-----
  use amr_commons
  use pm_commons
  use subsub_commons
  use subsub_parameters
  use hydro_parameters, ONLY: gamma
  use hydro_commons
  use cooling_module, ONLY: twopi
  use mpi_mod
  implicit none

!  real(dp), dimension(1:twotondim, 1:subsub_nhydro) :: edgeBC
!  integer :: ncloud, ncloud_max
!  real(dp), dimension(1:ncloud_max, 1:subsub_nhydro+ndim) :: clouds
!
  !! Local variables
  integer :: ilevel, igrid, i, ipart, jpart, ivar
  integer :: ncache, ngrid, npart1
  integer :: index_grid, index_cell, sink_ind
  integer :: ix, iy, iz, ind
  logical :: okay
  real(dp) :: rho, u, v, w, p, ekin
  real(dp) :: rx, ry, rz, d2, scale
  integer :: nx_loc, info


  real(dp) :: dx, xx, xx2, dx_min, pi, factG
  real(dp), dimension(1:twotondim, ndim) :: xc
  real(dp), dimension(1:3) :: skip_loc
  integer :: ncloud, ncloudall

  real(dp), dimension(1:twotondim, 1:ndim) :: edge
  real(dp) :: halfbox
  
  integer :: nn

  real(dp), dimension(1:nsink) :: subsub_r2sink, subsub_r2k

!  integer :: i, j, k, idim, ivar
!  real(dp) :: halfbox
!  real(dp) :: r, u, v, w, p, ekin
!  real(dp), dimension(1:twotondim, 1:ndim) :: edge
!
!  real(dp) :: smalle2, dd2, ww
!  real(dp), dimension(1:twotondim, 1:subsub_nhydro) :: sumhydro
!  real(dp), dimension(1:twotondim) :: sumweight

!sink_ind
!ind_grid

  !!----- Allocate nsink arrays
  allocate(subsub_eBC(1:nsink, 1:twotondim, 0:subsub_nhydro+1))
  allocate(subsub_eBCall(1:nsink, 1:twotondim, 0:subsub_nhydro+1))

  
  subsub_eBC(:,:,:) = 0.0D0
  subsub_eBCall(:,:,:) = 0.0D0

  !!----- Edge index
  halfbox = subsub_boxlen/2.0D0
  edge(1,:) = (/-halfbox, -halfbox, -halfbox/)
  edge(2,:) = (/+halfbox, -halfbox, -halfbox/)
  edge(3,:) = (/-halfbox, +halfbox, -halfbox/)
  edge(4,:) = (/+halfbox, +halfbox, -halfbox/)
  edge(5,:) = (/-halfbox, -halfbox, +halfbox/)
  edge(6,:) = (/+halfbox, -halfbox, +halfbox/)
  edge(7,:) = (/-halfbox, +halfbox, +halfbox/)
  edge(8,:) = (/+halfbox, +halfbox, +halfbox/)


  do ind=1,twotondim
    iz=(ind-1)/4
    iy=(ind-1-4*iz)/2
    ix=(ind-1-2*iy-4*iz)
    xc(ind,1)=(dble(ix)-0.5D0)
    xc(ind,2)=(dble(iy)-0.5D0)
    xc(ind,3)=(dble(iz)-0.5D0)
  end do

  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)

  
  dx_min=scale*0.5d0**(nlevelmax-nlevelsheld)/aexp

  pi=twopi/2d0
  factG=1.0D0
  if(cosmo)factG=3d0/8d0/pi*omega_m*aexp

  !!----- Compute r2sink first
  do i=1, nsink
    subsub_r2sink(i) = (factG * msink(i) / (subsub_v2sink(i) + subsub_cs2sink(i)))**2
    subsub_r2k(i) = min( max(subsub_r2sink(i), (dx_min/4.0D0)**2), (2.0*dx_min)**2 )
  enddo


  !!----- Get Cloud particles in this domain
  ncloud = 0

  do ilevel=levelmin, nlevelmax
    ncache = active(ilevel)%ngrid
    do igrid=1, ncache, nvector
      ngrid = MIN(nvector, ncache-igrid+1)
      do i=1, ngrid
        index_grid = active(ilevel)%igrid(igrid+i-1)

        npart1 = numbp(index_grid)
        if(npart1 .le. 0) cycle
    

    
        ipart = headp(index_grid)
        do jpart=1, npart1
          if(is_cloud(typep(ipart))) then

            !!----- Find cells
            ix=0; iy=0; iz=0
            if(xp(ipart,1)/scale .gt. xg(index_grid,1)) ix=1
            if(xp(ipart,2)/scale .gt. xg(index_grid,2)) iy=1
            if(xp(ipart,3)/scale .gt. xg(index_grid,3)) iz=1
            ind = 1 + ix + 2*iy + 4*iz

            index_cell = ncoarse + (ind-1)*ngridmax + index_grid

            !!----- Retrieve cell properties (primitive)

            rho = MAX(uold(index_cell,1), subsub_dfloor)
            u = uold(index_cell,2)/rho
            v = uold(index_cell,3)/rho
            w = uold(index_cell,4)/rho
            ekin = 0.5D0 * (u**2 + v**2 + w**2) * rho
            p = (uold(index_cell,5) - ekin)*(gamma-1.0D0)

            u = u - vp(ipart,1)
            v = v - vp(ipart,2)
            w = w - vp(ipart,3)

            !!----- Compute Weight
            sink_ind = -idp(ipart)

            ncloud = ncloud + 1
            do ind=1, twotondim
              rx = xsink(sink_ind,1)/scale + edge(ind,1)
              ry = xsink(sink_ind,2)/scale + edge(ind,2)
              rz = xsink(sink_ind,3)/scale + edge(ind,3)

              d2 = (xp(ipart,1)/scale - rx)**2 + (xp(ipart,2)/scale - ry)**2 + (xp(ipart,3)/scale - rz)**2! + subsub_smallr**2
              d2 = exp(-d2 / subsub_r2k(sink_ind))
              !d2 = 1.0D0 / d2

            
              !subsub_eBC(sink_ind,ind,6) = subsub_eBC(sink_ind,ind,6) + 1.0D0
              subsub_eBC(sink_ind,ind,0) = subsub_eBC(sink_ind,ind,0) + d2
              subsub_eBC(sink_ind,ind,1) = subsub_eBC(sink_ind,ind,1) + rho*d2
              subsub_eBC(sink_ind,ind,2) = subsub_eBC(sink_ind,ind,2) + u*d2
              subsub_eBC(sink_ind,ind,3) = subsub_eBC(sink_ind,ind,3) + v*d2
              subsub_eBC(sink_ind,ind,4) = subsub_eBC(sink_ind,ind,4) + w*d2
              subsub_eBC(sink_ind,ind,5) = subsub_eBC(sink_ind,ind,5) + p*d2



            enddo
          endif
          ipart = nextp(ipart)
        enddo
      enddo
    enddo
  enddo

  !! MPI Communicate
#ifndef WITHOUTMPI
  nn = nsink * twotondim * (2+subsub_nhydro)
  
  call MPI_ALLREDUCE(subsub_eBC(1,1,0),subsub_eBCall(1,1,0),nn,subsub_mpidp,MPI_SUM,MPI_COMM_WORLD,info)

  !call MPI_ALLREDUCE(ncloud,ncloudall,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  !call MPI_ALLREDUCE(xx,xx2,1,subsub_mpidp,MPI_SUM,MPI_COMM_WORLD,info)
#endif

  


  !! Input to subsub obj
  if(subsub_end .gt. 0) then
    do i=1, subsub_end
      sink_ind = subsub_obj(i)%sink_ind

      do ind=1, twotondim
        do ivar=1, subsub_nhydro
          !if(ivar.eq.1 .or. ivar.eq. subsub_nhydro)then
            subsub_obj(i)%edgeBC(ind,ivar) = subsub_eBCall(sink_ind,ind,ivar) / subsub_eBCall(sink_ind,ind,0)
          !else
          !  subsub_obj(i)%edgeBC(ind,ivar) = subsub_eBCall(sink_ind,ind,ivar) / subsub_eBCall(sink_ind,ind,6)
          !endif
        enddo
      enddo
    enddo
  endif

  !!---- Deallocate
  deallocate(subsub_eBC)
  deallocate(subsub_eBCall)
  

end subroutine subsub_edgeBC
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_faceBC
  use subsub_commons
  use amr_commons
  use hydro_parameters, ONLY:gamma
  implicit none
  
  !! Local Variables
  integer :: objind
  integer :: ud, idim, ivar, ind1, ind2

  integer,dimension(1:2,1:ndim,1:4) :: edge_index
  
  real(dp) :: frac, v1, v2, w1, w2
  integer :: i0, i1, j0, j1
  real(dp) :: rho, u, v, w, p, ekin

  if(subsub_end .eq. 0) return

  edge_index(1,1,:) = (/1,3,5,7/) 
  edge_index(2,1,:) = (/2,4,6,8/)

  edge_index(1,2,:) = (/1,2,5,6/)
  edge_index(2,2,:) = (/3,4,7,8/)

  edge_index(1,3,:) = (/1,2,3,4/)
  edge_index(2,3,:) = (/5,6,7,8/)

  frac = subsub_dx /subsub_boxlen

  
  !$omp parallel private(objind, ud, idim, i0, i1, j0, j1, ind1, ind2, v1, v2, w1, w2, ivar, rho, u, v, w, p, ekin)
  do objind=1, subsub_end
    do ud=1, 2 !! up or down
      do idim=1,ndim !! x, y, z
  
        i0 = edge_index(ud,idim,1)
        i1 = edge_index(ud,idim,2)
        j0 = edge_index(ud,idim,3)
        j1 = edge_index(ud,idim,4)
  
        !$omp do collapse(2)
        do ind1=1, subsub_ngrid
        do ind2=1, subsub_ngrid
          v1 = (dble(ind1)-0.5D0)*frac
          v2 = (dble(ind2)-0.5D0)*frac
  
          do ivar=1, subsub_nhydro
            w1 = v1 * (subsub_obj(objind)%edgeBC(i1,ivar) - subsub_obj(objind)%edgeBC(i0,ivar)) + subsub_obj(objind)%edgeBC(i0,ivar)
            w2 = v1 * (subsub_obj(objind)%edgeBC(j1,ivar) - subsub_obj(objind)%edgeBC(j0,ivar)) + subsub_obj(objind)%edgeBC(j0,ivar)
  
            subsub_obj(objind)%faceBC(ud,idim,ind1,ind2,ivar) = v2 * (w2-w1) + w1 
          enddo
          !! Convert to conservative
          rho = max(subsub_obj(objind)%faceBC(ud,idim,ind1,ind2,1),subsub_dfloor)
          u = subsub_obj(objind)%faceBC(ud,idim,ind1,ind2,2)
          v = subsub_obj(objind)%faceBC(ud,idim,ind1,ind2,3)
          w = subsub_obj(objind)%faceBC(ud,idim,ind1,ind2,4)
          p = max(subsub_obj(objind)%faceBC(ud,idim,ind1,ind2,5),subsub_pfloor)
          ekin = 0.5D0*(u**2 + v**2 + w**2)*rho

          subsub_obj(objind)%faceBC(ud,idim,ind1,ind2,1) = rho
          subsub_obj(objind)%faceBC(ud,idim,ind1,ind2,2) = u*rho
          subsub_obj(objind)%faceBC(ud,idim,ind1,ind2,3) = v*rho
          subsub_obj(objind)%faceBC(ud,idim,ind1,ind2,4) = w*rho
          subsub_obj(objind)%faceBC(ud,idim,ind1,ind2,5) = p/(gamma-1.0D0) + ekin
          
        enddo
        enddo
        !$omp end do
      enddo
    enddo
  enddo
  !$omp end parallel


end subroutine subsub_faceBC
!################################################################
!################################################################
!################################################################
!################################################################
subroutine subsub_sinkprop
  use subsub_commons
  use hydro_parameters, ONLY : gamma
  use pm_commons
  use mpi_mod
  implicit none

  !! Local variables
  integer :: myobj, i
  integer :: info
  real(dp), dimension(1:nsink) :: v2_local, cs2_local

  real(dp) :: rhotot, v2tmp, cs2tmp
  real(dp) :: rho, u, v, w, p, ekin

  !! Initialize
  allocate(subsub_v2sink(1:nsink))
  allocate(subsub_cs2sink(1:nsink))
  
  v2_local(:) = 0.0D0
  cs2_local(:) = 0.0D0
  subsub_v2sink(:) = 0.0D0
  subsub_cs2sink(:) = 0.0D0

  !! Retrieve
  if(subsub_end .gt. 0) then
    do myobj=1, subsub_end
      !! never computed; compute here
      if(subsub_obj(myobj)%v2 .lt. 0)then
        rhotot = 0.0D0
        v2tmp = 0.0D0
        cs2tmp = 0.0D0
        !$omp parallel do private(rho, u, v, w, p, ekin) reduction(+:rhotot, v2tmp, cs2tmp)
        do i=1, subsub_nn
          !! Mass-weighted v2 and cs2
          rho = max(subsub_obj(myobj)%hydro(i,1), subsub_dfloor)
          u = subsub_obj(myobj)%hydro(i,2) / rho
          v = subsub_obj(myobj)%hydro(i,3) / rho
          w = subsub_obj(myobj)%hydro(i,4) / rho
          ekin = (u**2 + v**2 + w**2)*0.5D0 * rho
          p = max((subsub_obj(myobj)%hydro(i,5) - ekin) / (gamma - 1.0D0), subsub_pfloor)

          rhotot = rhotot + rho
          v2tmp = v2tmp + ekin*2.0D0
          cs2tmp = cs2tmp + rho * (p/rho)
        enddo
        !$omp end parallel do

        subsub_obj(myobj)%v2 = v2tmp
        subsub_obj(myobj)%cs2 = cs2tmp
      endif
      v2_local(subsub_obj(myobj)%sink_ind) = subsub_obj(myobj)%v2
      cs2_local(subsub_obj(myobj)%sink_ind) = subsub_obj(myobj)%cs2


    enddo
  endif

  !! Communicate
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(v2_local(1), subsub_v2sink(1), nsink, subsub_mpidp, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(cs2_local(1), subsub_cs2sink(1), nsink, subsub_mpidp, MPI_SUM, MPI_COMM_WORLD, info)
#else
  subsub_v2sink(:) = v2_local(:)
  subsub_cs2sink(:) = subsub_cs2sink(:)
#endif
  
  !! Deallocation
  deallocate(subsub_v2sink)
  deallocate(subsub_cs2sink)
end subroutine subsub_sinkprop
!################################################################
!################################################################
!################################################################
!################################################################