  subroutine iterate_time_static_solver2()
  
  use conjugate_gradient 
  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_movie
  use gravity_perturbation, only : gravity_timeseries, GRAVITY_SIMULATION
!  use fault_solver_dynamic, only : bc_dynflt_set3d_all,SIMULATION_TYPE_DYN,faults
!  use fault_solver_kinematic, only : bc_kinflt_set_all,SIMULATION_TYPE_KIN
  use fault_solver_qstatic, only: bc_qstaticflt_set3d_all,faults


  implicit none

  type(CG_data) ::  CG_problem
  type(CG_Vector) :: CG_size
  logical,dimension(3,NGLOB_AB) :: MASK_default
  real(kind=CUSTOM_REAL) :: Max_error
  real(kind=CUSTOM_REAL) :: Max_error_all
  MASK_default(:,:) = .true.


  CG_size%NDIM = 3
  CG_size%NELE = NGLOB_AB


  time_start = wtime()


 load(:,:)=0.0_CUSTOM_REAL
!  write(*,*) 'compute the load vector'
 
!  call compute_AX(load , displ , MASK_default, MASK_default)
  if(faults(1)%nglob>0) then
  load(:,faults(1)%ibulk1) = 0.0_CUSTOM_REAL
  load(:,faults(1)%ibulk2) = 0.0_CUSTOM_REAL
  load(1,faults(1)%ibulk1) = load(1,faults(1)%ibulk1)+10.0e6_CUSTOM_REAL*faults(1)%B(:)
  load(1,faults(1)%ibulk2) = load(1,faults(1)%ibulk2)-10.0e6_CUSTOM_REAL*faults(1)%B(:)
  endif

  call CG_initialize(CG_problem,CG_size,displ,load)
!  write(*,*) faults(1)%ibulk1 
!  call CG_mask(CG_problem,faults(1)%ibulk1,faults(1)%ibulk1)
!  call CG_mask(CG_problem,faults(1)%ibulk2,faults(1)%ibulk2)
!  write(*,*) CG_problem%MASKX 
  write(*,*) 'successfully initialize the CG_Problem',myrank
  write(*,*) CG_problem%CG_size%NELE
!  write(*,*) size(CG_problem%Pdirection)
!  call CG_initialize_preconditioner(CG_problem,rmassx,rmassy,rmassz)
  do it = 1,1000
   
!  write(*,*) 'successfully get into the loop!'
  call update_value_direction(CG_problem)
!  write(*,*) 'successfully get the ',it,'step'
  write(*,*) maxval(abs(CG_problem%Residue)),'at',it,'myrank=',myrank   
  call it_print_elapsed_time()

  Max_error = maxval(abs(CG_problem%Residue))
 call  max_all_cr(Max_error,Max_error_all)
  if(myrank ==0) & 
  write(IMAIN,*)  'max_error=',Max_error_all
  
  
  
  enddo   ! end of main time loop
  !displ(1,faults(1)%ibulk1) =  displ(1,faults(1)%ibulk1) + 1.0_CUSTOM_REAL
  !displ(1,faults(1)%ibulk2) =  displ(1,faults(1)%ibulk2) - 1.0_CUSTOM_REAL
!  
  call compute_AX(load , displ , MASK_default, MASK_default)
 
  call bc_qstaticflt_set3d_all(load,veloc,displ)
  
  displ = displ + veloc 


  call CG_mask(CG_problem,faults(1)%ibulk1,faults(1)%ibulk1)
  call CG_mask(CG_problem,faults(1)%ibulk2,faults(1)%ibulk2) 
  call compute_AX(load , veloc , MASK_default, MASK_default)
  write(*,*) 'max v',maxval(veloc(1,:))
  write(*,*) 'max load',maxval(load(1,:))
  load(:,faults(1)%ibulk1) = 0.0_CUSTOM_REAL
  load(:,faults(1)%ibulk2) = 0.0_CUSTOM_REAL
  call updateResidue(CG_problem,load)

  CG_problem%Residue(:,faults(1)%ibulk1) = 0.0_CUSTOM_REAL
  CG_problem%Residue(:,faults(1)%ibulk2) = 0.0_CUSTOM_REAL


  do it = 1001,2000
!  write(*,*) 'successfully get into the loop!'
  call update_value_direction(CG_problem)
!  write(*,*) 'successfully get the ',it,'step'
!  write(*,*) maxval(abs(CG_problem%Residue))
  enddo   ! end of main time loop
   
  call compute_AX(load , displ , MASK_default, MASK_default)
 
  call bc_qstaticflt_set3d_all(load,veloc,displ)
 





    if( MOVIE_SIMULATION ) then
      call write_movie_output()
    endif
    write(*,*) maxval(rmassx(:))
!  write(*,*) maxval(displ(1,:)),minval(CG_problem%X(1,:))
  call it_print_elapsed_time()

  end subroutine iterate_time_static_solver2



