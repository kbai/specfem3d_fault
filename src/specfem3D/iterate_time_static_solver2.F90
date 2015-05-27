  subroutine iterate_time_static_solver2()
  
  use conjugate_gradient 
  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_movie
  use gravity_perturbation, only : gravity_timeseries, GRAVITY_SIMULATION
  use fault_solver_dynamic, only : bc_dynflt_set3d_all,SIMULATION_TYPE_DYN,faults
  use fault_solver_kinematic, only : bc_kinflt_set_all,SIMULATION_TYPE_KIN


  implicit none

  type(CG_data) ::  CG_problem
  type(CG_Vector) :: CG_size
  logical,dimension(3,NGLOB_AB) :: MASK_default

  MASK_default(:,:) = .true.


  CG_size%NDIM = 3
  CG_size%NELE = NGLOB_AB


  time_start = wtime()


!  if (SIMULATION_TYPE_DYN) call bc_dynflt_set3d_all(load,veloc,displ)
!  if (SIMULATION_TYPE_KIN) call bc_kinflt_set_all(load,veloc,displ)
!  load(1,:) = load(1,:)
!  load(2,:) = load(2,:)
!  load(3,:) = load(3,:)
  load(:,:)=0.0_CUSTOM_REAL
  write(*,*) 'compute the load vector'
  displ(1,faults(1)%ibulk1) =  displ(1,faults(1)%ibulk1) + 1.0_CUSTOM_REAL
  displ(1,faults(1)%ibulk2) =  displ(1,faults(1)%ibulk2) - 1.0_CUSTOM_REAL
!  displ(:,:) = 0.0_CUSTOM_REAL
  
  call compute_AX(load , displ , MASK_default, MASK_default)
  load(:,faults(1)%ibulk1) = 0.0_CUSTOM_REAL
  load(:,faults(1)%ibulk2) = 0.0_CUSTOM_REAL

  call CG_initialize(CG_problem,CG_size,displ,load)
!  write(*,*) faults(1)%ibulk1 
  call CG_mask(CG_problem,faults(1)%ibulk1,faults(1)%ibulk1)
  call CG_mask(CG_problem,faults(1)%ibulk2,faults(1)%ibulk2)
  write(*,*) CG_problem%MASKX 
  write(*,*) 'successfully initialize the CG_Problem'
  write(*,*) CG_problem%CG_size%NELE
!  write(*,*) size(CG_problem%Pdirection)
  do it = 1,1000
   
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

  end subroutine iterate_time_static_solver2



