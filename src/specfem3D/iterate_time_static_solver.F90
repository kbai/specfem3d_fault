  subroutine iterate_time_static_solver()
  
  use conjugate_gradient 
  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_movie
  use gravity_perturbation, only : gravity_timeseries, GRAVITY_SIMULATION
  use fault_solver_dynamic, only : bc_dynflt_set3d_all,SIMULATION_TYPE_DYN
  use fault_solver_kinematic, only : bc_kinflt_set_all,SIMULATION_TYPE_KIN


  implicit none

  type(CG_data) ::  CG_problem
  type(CG_Vector) :: CG_size

  CG_size%NDIM = 3
  CG_size%NELE = NGLOB_AB


  time_start = wtime()


  if (SIMULATION_TYPE_DYN) call bc_dynflt_set3d_all(load,veloc,displ)
  if (SIMULATION_TYPE_KIN) call bc_kinflt_set_all(load,veloc,displ)

  load(1,:) = load(1,:)
  load(2,:) = load(2,:)
  load(3,:) = load(3,:)
  write(*,*) 'compute the load vector'
  call CG_initialize(CG_problem,CG_size,displ,load)
  write(*,*) 'successfully initialize the CG_Problem'
!  write(*,*) size(CG_problem%Pdirection)
  do it = 1,10000
   
  write(*,*) 'successfully get into the loop!'
  call update_value_direction(CG_problem)
  write(*,*) 'successfully get the ',it,'step'
  write(*,*) maxval(abs(CG_problem%Residue))
  enddo   ! end of main time loop

   if( MOVIE_SIMULATION ) then
      call write_movie_output()
    endif

  write(*,*) maxval(displ(1,:)),minval(CG_problem%X(1,:))
  call it_print_elapsed_time()

  end subroutine iterate_time_static_solver



