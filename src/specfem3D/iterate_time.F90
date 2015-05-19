  subroutine iterate_time()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_movie
  use gravity_perturbation, only : gravity_timeseries, GRAVITY_SIMULATION

  implicit none

  time_start = wtime()



  do it = 1,1
    ! simulation status output and stability check
    if( mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5 .or. it == NSTEP ) &
      call check_stability()

    ! simulation status output and stability check
    ! updates wavefields using Newmark time scheme

    ! calculates stiffness term
    if( .not. GPU_MODE )then
      ! wavefields on CPU

       ! forward simulations

        ! 1. acoustic domain
        if( ACOUSTIC_SIMULATION ) call compute_forces_acoustic()
        ! 2. elastic domain
        if( ELASTIC_SIMULATION ) call compute_forces_viscoelastic()

      ! poroelastic solver
      if( POROELASTIC_SIMULATION ) call compute_forces_poroelastic()

    else
      ! wavefields on GPU
      ! acoustic solver
      if( ACOUSTIC_SIMULATION ) call compute_forces_acoustic_GPU()
      ! elastic solver
      ! (needs to be done first, before poroelastic one)
      if( ELASTIC_SIMULATION ) call compute_forces_viscoelastic_GPU()
    endif

   if( MOVIE_SIMULATION ) then
      call write_movie_output()
    endif

    ! first step of noise tomography, i.e., save a surface movie at every time step
!
!---- end of time iteration loop
!
  enddo   ! end of main time loop

  call it_print_elapsed_time()

  ! Transfer fields from GPU card to host for further analysis
  if( GPU_MODE ) call it_transfer_from_GPU()

!----  close energy file
  if( OUTPUT_ENERGY .and. myrank == 0 ) close(IOUT_ENERGY)

  end subroutine iterate_time


!=====================================================================


  subroutine it_read_forward_arrays()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  implicit none

  integer :: ier

! restores last time snapshot saved for backward/reconstruction of wavefields
! note: this is done here after the Newmark time scheme, otherwise the indexing for sources
!          and adjoint sources will become more complicated
!          that is, index it for adjoint sources will match index NSTEP - 1 for backward/reconstructed wavefields
  if (ADIOS_FOR_FORWARD_ARRAYS) then
    call read_forward_arrays_adios()
  else
    ! reads in wavefields
    open(unit=IIN,file=trim(prname)//'save_forward_arrays.bin',status='old',&
          action='read',form='unformatted',iostat=ier)
    if( ier /= 0 ) then
      print*,'error: opening save_forward_arrays'
      print*,'path: ',trim(prname)//'save_forward_arrays.bin'
      call exit_mpi(myrank,'error open file save_forward_arrays.bin')
    endif

    if( ACOUSTIC_SIMULATION ) then
      read(IIN) b_potential_acoustic
      read(IIN) b_potential_dot_acoustic
      read(IIN) b_potential_dot_dot_acoustic
    endif

    ! elastic wavefields
    if( ELASTIC_SIMULATION ) then
      read(IIN) b_displ
      read(IIN) b_veloc
      read(IIN) b_accel
      ! memory variables if attenuation
      if( ATTENUATION ) then
        if(FULL_ATTENUATION_SOLID) read(IIN) b_R_trace
        read(IIN) b_R_xx
        read(IIN) b_R_yy
        read(IIN) b_R_xy
        read(IIN) b_R_xz
        read(IIN) b_R_yz
        if(FULL_ATTENUATION_SOLID) read(IIN) b_epsilondev_trace
        read(IIN) b_epsilondev_xx
        read(IIN) b_epsilondev_yy
        read(IIN) b_epsilondev_xy
        read(IIN) b_epsilondev_xz
        read(IIN) b_epsilondev_yz
      endif ! ATTENUATION
    endif

    ! poroelastic wavefields
    if( POROELASTIC_SIMULATION ) then
      read(IIN) b_displs_poroelastic
      read(IIN) b_velocs_poroelastic
      read(IIN) b_accels_poroelastic
      read(IIN) b_displw_poroelastic
      read(IIN) b_velocw_poroelastic
      read(IIN) b_accelw_poroelastic
    endif

    close(IIN)
  endif

  if(GPU_MODE) then
    if( ACOUSTIC_SIMULATION ) then
    ! transfers fields onto GPU
      call transfer_b_fields_ac_to_device(NGLOB_AB,b_potential_acoustic, &
                                          b_potential_dot_acoustic,      &
                                          b_potential_dot_dot_acoustic,  &
                                          Mesh_pointer)
    endif
    ! elastic wavefields
    if( ELASTIC_SIMULATION ) then
      ! puts elastic wavefield to GPU
      call transfer_b_fields_to_device(NDIM*NGLOB_AB,b_displ,b_veloc,b_accel,Mesh_pointer)
      ! memory variables if attenuation
      if( ATTENUATION ) then
        call transfer_b_fields_att_to_device(Mesh_pointer,                    &
                           b_R_xx,b_R_yy,b_R_xy,b_R_xz,b_R_yz,                &
                           size(b_R_xx),                                      &
                           b_epsilondev_xx,b_epsilondev_yy,b_epsilondev_xy,   &
                           b_epsilondev_xz,b_epsilondev_yz,                   &
                           size(b_epsilondev_xx))
      endif
    endif
  endif

  end subroutine it_read_forward_arrays

!=====================================================================

  subroutine it_store_attenuation_arrays()

! resetting d/v/a/R/eps for the backward reconstruction with attenuation

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic

  implicit none

  integer :: ier

  if( it > 1 .and. it < NSTEP) then
    ! adjoint simulations

! note backward/reconstructed wavefields:
!       storing wavefield displ() at time step it, corresponds to time (it-1)*DT - t0 (see routine write_seismograms_to_file )
!       reconstucted wavefield b_displ() at it corresponds to time (NSTEP-it-1)*DT - t0
!       we read in the reconstructed wavefield at the end of the time iteration loop, i.e. after the Newmark scheme,
!       thus, indexing is NSTEP-it (rather than something like NSTEP-(it-1) )

    if (SIMULATION_TYPE == 3 .and. mod(NSTEP-it,NSTEP_Q_SAVE) == 0) then
      ! reads files content
      write(outputname,"('save_Q_arrays_',i6.6,'.bin')") NSTEP-it
      open(unit=IIN,file=trim(prname_Q)//trim(outputname),status='old',&
            action='read',form='unformatted',iostat=ier)
      if( ier /= 0 ) then
        print*,'error: opening save_Q_arrays'
        print*,'path: ',trim(prname_Q)//trim(outputname)
        call exit_mpi(myrank,'error open file save_Q_arrays_***.bin for reading')
      endif

      if( ELASTIC_SIMULATION ) then
        ! reads arrays from disk files
        read(IIN) b_displ
        read(IIN) b_veloc
        read(IIN) b_accel

        ! puts elastic fields onto GPU
        if(GPU_MODE) then
          ! wavefields
          call transfer_b_fields_to_device(NDIM*NGLOB_AB,b_displ,b_veloc,b_accel, Mesh_pointer)
        endif

        if(FULL_ATTENUATION_SOLID) read(IIN) b_R_trace
        read(IIN) b_R_xx
        read(IIN) b_R_yy
        read(IIN) b_R_xy
        read(IIN) b_R_xz
        read(IIN) b_R_yz
        if(FULL_ATTENUATION_SOLID) read(IIN) b_epsilondev_trace
        read(IIN) b_epsilondev_xx
        read(IIN) b_epsilondev_yy
        read(IIN) b_epsilondev_xy
        read(IIN) b_epsilondev_xz
        read(IIN) b_epsilondev_yz

        ! puts elastic fields onto GPU
        if(GPU_MODE) then
          ! attenuation arrays
          call transfer_b_fields_att_to_device(Mesh_pointer, &
                  b_R_xx,b_R_yy,b_R_xy,b_R_xz,b_R_yz,size(b_R_xx), &
!!!               b_R_trace,b_R_xx,b_R_yy,b_R_xy,b_R_xz,b_R_yz,size(b_R_xx), &
                  b_epsilondev_xx,b_epsilondev_yy,b_epsilondev_xy,b_epsilondev_xz,b_epsilondev_yz, &
!!!               b_epsilondev_trace,b_epsilondev_xx,b_epsilondev_yy,b_epsilondev_xy,b_epsilondev_xz,b_epsilondev_yz, &
                  size(b_epsilondev_xx))
        endif
      endif

      if( ACOUSTIC_SIMULATION ) then
        ! reads arrays from disk files
        read(IIN) b_potential_acoustic
        read(IIN) b_potential_dot_acoustic
        read(IIN) b_potential_dot_dot_acoustic

        ! puts acoustic fields onto GPU
        if(GPU_MODE) &
          call transfer_b_fields_ac_to_device(NGLOB_AB,b_potential_acoustic, &
                              b_potential_dot_acoustic, b_potential_dot_dot_acoustic, Mesh_pointer)

      endif
      close(IIN)

    else if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD .and. mod(it,NSTEP_Q_SAVE) == 0) then
      ! stores files content
      write(outputname,"('save_Q_arrays_',i6.6,'.bin')") it
      open(unit=IOUT,file=trim(prname_Q)//trim(outputname),status='unknown',&
           action='write',form='unformatted',iostat=ier)
      if( ier /= 0 ) then
        print*,'error: opening save_Q_arrays'
        print*,'path: ',trim(prname_Q)//trim(outputname)
        call exit_mpi(myrank,'error open file save_Q_arrays_***.bin for writing')
      endif

      if( ELASTIC_SIMULATION ) then
        ! gets elastic fields from GPU onto CPU
        if(GPU_MODE) then
          call transfer_fields_el_from_device(NDIM*NGLOB_AB,displ,veloc, accel, Mesh_pointer)
        endif

        ! writes to disk file
        write(IOUT) displ
        write(IOUT) veloc
        write(IOUT) accel

        if(GPU_MODE) then
          ! attenuation arrays
          call transfer_fields_att_from_device(Mesh_pointer, &
                     R_xx,R_yy,R_xy,R_xz,R_yz,size(R_xx), &
!!!                  R_trace,R_xx,R_yy,R_xy,R_xz,R_yz,size(R_xx), &
                     epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz, &
!!!                  epsilondev_trace,epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz, &
                     size(epsilondev_xx))
        endif

        if(FULL_ATTENUATION_SOLID) write(IOUT) R_trace
        write(IOUT) R_xx
        write(IOUT) R_yy
        write(IOUT) R_xy
        write(IOUT) R_xz
        write(IOUT) R_yz
        if(FULL_ATTENUATION_SOLID) write(IOUT) epsilondev_trace
        write(IOUT) epsilondev_xx
        write(IOUT) epsilondev_yy
        write(IOUT) epsilondev_xy
        write(IOUT) epsilondev_xz
        write(IOUT) epsilondev_yz
      endif

      if( ACOUSTIC_SIMULATION ) then
        ! gets acoustic fields from GPU onto CPU
        if(GPU_MODE) then
          call transfer_fields_ac_from_device(NGLOB_AB,potential_acoustic, &
                                              potential_dot_acoustic, potential_dot_dot_acoustic, Mesh_pointer)
        endif

        ! writes to disk file
        write(IOUT) potential_acoustic
        write(IOUT) potential_dot_acoustic
        write(IOUT) potential_dot_dot_acoustic
      endif
      close(IOUT)

    endif ! SIMULATION_TYPE
  endif ! it

  end subroutine it_store_attenuation_arrays

!=====================================================================

  subroutine it_print_elapsed_time()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic

  implicit none

  ! local parameters
  double precision :: tCPU
  integer :: ihours,iminutes,iseconds,int_tCPU

  if( myrank == 0 ) then
    ! elapsed time since beginning of the simulation
    tCPU = wtime() - time_start
    int_tCPU = int(tCPU)
    ihours = int_tCPU / 3600
    iminutes = (int_tCPU - 3600*ihours) / 60
    iseconds = int_tCPU - 3600*ihours - 60*iminutes
    write(IMAIN,*) 'Time-Loop Complete. Timing info:'
    write(IMAIN,*) 'Total elapsed time in seconds = ',tCPU
    write(IMAIN,"(' Total elapsed time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
  endif

  end subroutine it_print_elapsed_time

!=====================================================================

  subroutine it_transfer_from_GPU()

! transfers fields on GPU back onto CPU

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic

  implicit none

  ! to store forward wave fields
  if( SIMULATION_TYPE == 1 .and. SAVE_FORWARD ) then

    ! acoustic potentials
    if( ACOUSTIC_SIMULATION ) &
      call transfer_fields_ac_from_device(NGLOB_AB,potential_acoustic, &
                                          potential_dot_acoustic, potential_dot_dot_acoustic, &
                                          Mesh_pointer)

    ! elastic wavefield
    if( ELASTIC_SIMULATION ) then
      call transfer_fields_el_from_device(NDIM*NGLOB_AB,displ,veloc,accel,Mesh_pointer)

      if (ATTENUATION) &
        call transfer_fields_att_from_device(Mesh_pointer, &
                    R_xx,R_yy,R_xy,R_xz,R_yz,size(R_xx), &
!!!                 R_trace,R_xx,R_yy,R_xy,R_xz,R_yz,size(R_xx), &
                    epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz, &
!!!                 epsilondev_trace,epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz, &
                    size(epsilondev_xx))

    endif
  else if( SIMULATION_TYPE == 3 ) then

    ! to store kernels
    ! acoustic domains
    if( ACOUSTIC_SIMULATION ) then
      ! only in case needed...
      !call transfer_b_fields_ac_from_device(NGLOB_AB,b_potential_acoustic, &
      !                      b_potential_dot_acoustic, b_potential_dot_dot_acoustic, Mesh_pointer)

      ! acoustic kernels
      call transfer_kernels_ac_to_host(Mesh_pointer,rho_ac_kl,kappa_ac_kl,NSPEC_AB)
    endif

    ! elastic domains
    if( ELASTIC_SIMULATION ) then
      ! only in case needed...
      !call transfer_b_fields_from_device(NDIM*NGLOB_AB,b_displ,b_veloc,b_accel, Mesh_pointer)

      ! elastic kernels
      call transfer_kernels_el_to_host(Mesh_pointer,rho_kl,mu_kl,kappa_kl,cijkl_kl,NSPEC_AB)
    endif

    ! specific noise strength kernel
    if( NOISE_TOMOGRAPHY == 3 ) then
      call transfer_kernels_noise_to_host(Mesh_pointer,Sigma_kl,NSPEC_AB)
    endif

    ! approximative hessian for preconditioning kernels
    if( APPROXIMATE_HESS_KL ) then
      if( ELASTIC_SIMULATION ) call transfer_kernels_hess_el_tohost(Mesh_pointer,hess_kl,NSPEC_AB)
      if( ACOUSTIC_SIMULATION ) call transfer_kernels_hess_ac_tohost(Mesh_pointer,hess_ac_kl,NSPEC_AB)
    endif

  endif

  ! frees allocated memory on GPU
  call prepare_cleanup_device(Mesh_pointer, &
                              ACOUSTIC_SIMULATION,ELASTIC_SIMULATION, &
                              STACEY_ABSORBING_CONDITIONS,NOISE_TOMOGRAPHY,COMPUTE_AND_STORE_STRAIN, &
                              ATTENUATION,ANISOTROPY,APPROXIMATE_OCEAN_LOAD, &
                              APPROXIMATE_HESS_KL)

  end subroutine it_transfer_from_GPU
