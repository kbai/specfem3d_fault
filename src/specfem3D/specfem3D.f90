  subroutine specfem3D()

  use specfem_par
! ************** PROGRAM STARTS HERE **************

  ! force Flush-To-Zero if available to avoid very slow Gradual Underflow trapping
!  do while(Debugwait) 
!  enddo
  call force_ftz()
  
  ! reads in parameters
  call initialize_simulation()

  ! reads in external mesh
  if (ADIOS_FOR_MESH) then
    call read_mesh_databases_adios()
  else
    call read_mesh_databases()
  endif
    !call read_mesh_databases()
    !call read_mesh_databases_adios()


! sets up reference element GLL points/weights/derivatives
  call setup_GLL_points()

  
! detects surfaces
  call detect_mesh_surfaces()


! prepares sources and receivers
  call setup_sources_receivers()


! sets up and precomputes simulation arrays
  call prepare_timerun()


! steps through time iterations
  call iterate_time_static_solver()

! saves last time frame and finishes kernel calculations
  call finalize_simulation()

  end subroutine specfem3D

