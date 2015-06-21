#define XPOSITION 250000.0e0_CUSTOM_REAL
#define YPOSITION 250000.0e0_CUSTOM_REAL
#define ZPOSITION 0.0e0_CUSTOM_REAL
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

    type(CG_data) ::  CGC,CG
    type(CG_Vector) :: CG_size
    logical,dimension(3,NGLOB_AB) :: MASK_default
    real(kind=CUSTOM_REAL) :: Max_error
    real(kind=CUSTOM_REAL) :: Max_error_all
    real(kind=CUSTOM_REAL),dimension(3,NGLOB_AB) :: sload,displ2,displ3,precon
    integer :: max_loc,ii
    integer :: NGLOB_AB_ALL
    integer,dimension(1) :: Boundary_loc
    write(*,*) "ready to compute diagonal!"
    write(IMAIN,*) "xigll",xigll
    write(IMAIN,*) "wxgll",wxgll
    call print_computer_name()
    call make_load()
    call get_boundary()
    call prepare_MPI
    call prepare_CGC()

 !   displ(:,:)= 0.0_CUSTOM_REAL
 !   if(myrank == 10)    displ(1,ibool(1,5,1,980))=1.0
 !   if(myrank == 10)    write(*,*) "boundaryelement:", phase_ispec_inner_elastic(:,1)



 !   call compute_AX(load,displ,MASK_default,MASK_default)
 !   call compute_diagonal(veloc)
 !   if(myrank == 10) write(*,*) "checkdiagonal:", veloc(1,ibool(1,5,1,980)) ,load(1,ibool(1,5,1,980))
 !   displ(:,:)= 0.0_CUSTOM_REAL
 !   if(myrank == 10)    displ(1,ibool(1,5,1,980))=1.0
 !   if(myrank == 10)    write(*,*) "statistic:",nspec_inner_elastic,nspec_outer_elastic


  !  call compute_AX2(load,displ,MASK_default,MASK_default)
  !  call compute_diagonal2(veloc)
  !  if(myrank == 10) write(*,*) "checkdiagonal:",veloc(1,ibool(1,5,1,980)) ,load(1,ibool(1,5,1,980))




    write(IMAIN,*) maxval(sload)
  do it = 1,10000
   !     call update_value_direction(CGC)
        call Jacobi_method(CGC)
        Max_error = CGC%Norm_old
!        call compute_Anorm(max_error,CGC)
        if(myrank ==0) &
            write(IMAIN,*)  'max_error=',Max_error
    enddo
    write(*,*) maxval(displ)
    call smoother(NSPEC_AB,NGLOB_AB,displ2,displ,ibool)
    displ = displ2
     call write_movie_output()


    call prepare_CG()
!    displ = displ2
 !   call write_movie_output()
    do it = 1,1000
        call update_value_direction(CG)
        Max_error = CG%Norm_old
!        call compute_Anorm(max_error,CG)

        if(myrank ==0) &
            write(IMAIN,*)  'max_error=',Max_error

    enddo

    displ = displ2

    write(IMAIN,*) "maxdispl",max_all(displ),max_all(displ2)
    write(IMAIN,*) "mindispl",min_all(displ),min_all(displ2)



    call restriction_call(CGC%Residue,CG%Residue,MASK_default,MASK_default)
    CGC%Residue = CGC%Residue * MASK_default
    CGC%X = 0.0_CUSTOM_REAL
    write(*,*) "check max displ:",maxval(abs(displ))
    call Reinitialize(CGC)
    do it = 1,1000
        call update_value_direction(CGC)
        Max_error = CGC%Norm_old
!        call compute_Anorm(max_error,CGC)

        if(myrank ==0) &
            write(IMAIN,*)  'max_error=',Max_error
    enddo

    call smoother(NSPEC_AB,NGLOB_AB,displ3,displ,ibool)
    call move_x(CG,displ3)
    do it = 1,1000
        call update_value_direction(CG)
        Max_error = CG%Norm_old
!        call compute_Anorm(max_error,CG)

        if(myrank ==0) &
            write(IMAIN,*)  'max_error=',Max_error

    enddo





   ! do ii = 1,4
   ! do it = 1,10000
   !     call update_value_direction(CGC)
   !     call it_print_elapsed_time()
!        Max_error = sum(abs(CGC%Residue))
   !     Max_error = CGC%Norm_old
!            call sum_all_cr(Max_error,Max_error_all)
    !        if(myrank ==0) &
    !        write(IMAIN,*)  'max_error=',Max_error_all
  !  enddo   ! end of main time loop


  !  call smoother(NSPEC_AB, NGLOB_AB ,displ ,ibool)
  !  if( MOVIE_SIMULATION ) then
  !      call write_movie_output()
  !  endif
  !  displ2 = displ2 + displ
   ! displ = 0.0_CUSTOM_REAL

    !displ = load
  !     if( MOVIE_SIMULATION ) then
  !      call write_movie_output()
  !  endif

   ! write(*,*) "diagonal computed!"

   ! call sync_all()
    ! write(*,*) "maxdispl:",maxval(abs(displ))
   ! do it = 1001,2000:q


    !    call update_value_direction(CG)




    !      Max_error = maxval(abs(veloc-displ)/abs(veloc+displ))

    !   max_loc = maxloc(abs(CG%Residue(1,:)),1)
    !     call  max_all_cr(Max_error,Max_error_all)

    !    if(myrank ==0) &
        !   write(IMAIN,*)  'max_error=',Max_error_all



    !  do it  =  1,NGLOB_AB
    !   if(zstore(it) > -1.0) then
    !    write(*,1000) xstore(it) , ystore(it) , veloc(1,it), displ(1,it)
    ! endif
    ! enddo

    !    if(mod(it,1000) == 0) call Reinitialize(CG)


    ! displ = displ + veloc


    !  call CG_mask(CGC,faults(1)%ibulk1,faults(1)%ibulk1)
    !  call CG_mask(CGC,faults(1)%ibulk2,faults(1)%ibulk2)
    !  call compute_AX(load , veloc , MASK_default, MASK_default)
    !  write(*,*) 'max v',maxval(veloc(1,:))
    !  write(*,*) 'max load',maxval(load(1,:))
    !  load(:,faults(1)%ibulk1) = 0.0_CUSTOM_REAL
    !  load(:,faults(1)%ibulk2) = 0.0_CUSTOM_REAL
    !  call updateResidue(CGC,load)

    !  CGC%Residue(:,faults(1)%ibulk1) = 0.0_CUSTOM_REAL
    !  CGC%Residue(:,faults(1)%ibulk2) = 0.0_CUSTOM_REAL


    !do it = 1001,2000
    !  write(*,*) 'successfully get into the loop!'
    !call update_value_direction(CGC)
    !  write(*,*) 'successfully get the ',it,'step'
    !  write(*,*) maxval(abs(CGC%Residue))
    !enddo   ! end of main time loop

    ! call compute_AX(load , displ , MASK_default, MASK_default)

    call bc_qstaticflt_set3d_all(load,displ,veloc)






    if( MOVIE_SIMULATION ) then
        !call write_movie_output()
    endif
    write(*,*) maxval(rmassx(:))
    !  write(*,*) maxval(displ(1,:)),minval(CGC%X(1,:))
    call it_print_elapsed_time()
    1000 format ('datacome',4e15.7)

contains
subroutine make_load()
     time_start = wtime()
    load(:,:)=0.0_CUSTOM_REAL
    MASK_default(:,:) = .true.

       !  call compute_AX(load , displ , MASK_default, MASK_default)
        if(faults(1)%nglob>0) then
            load(:,faults(1)%ibulk1) = 0.0_CUSTOM_REAL
            load(:,faults(1)%ibulk2) = 0.0_CUSTOM_REAL
            load(1,faults(1)%ibulk1) = load(1,faults(1)%ibulk1)+10.0e6_CUSTOM_REAL*faults(1)%B(:)
            load(1,faults(1)%ibulk2) = load(1,faults(1)%ibulk2)-10.0e6_CUSTOM_REAL*faults(1)%B(:)
                    write(*,*) "sum of B!",myrank,sum(faults(1)%B(:))

        endif
        write(IMAIN,*) maxval(load)
        call restriction_call(sload,load,MASK_default,MASK_default)
         CG_size%NDIM = 3
         CG_size%NELE = NGLOB_AB
end subroutine make_load
!================================================================================
subroutine prepare_CGC()

    displ2 = 0.0_CUSTOM_REAL
    call CG_initialize(CGC,CG_size,displ,sload,.true.,MASK_default,MASK_default)
    call compute_Diagonal2(precon)

    do ii = 1,NGLOB_AB
        if (precon(1,ii)>0.0_CUSTOM_REAL) then
    rmassx(ii) = 1.0_CUSTOM_REAL/precon(1,ii)
    rmassy(ii) = 1.0_CUSTOM_REAL/precon(2,ii)
    rmassz(ii) = 1.0_CUSTOM_REAL/precon(3,ii)
    else
    rmassx(ii) = 0.0_CUSTOM_REAL
    rmassy(ii) = 0.0_CUSTOM_REAL
    rmassz(ii) = 0.0_CUSTOM_REAL
    endif
    enddo

    call CG_initialize_preconditioner(CGC,rmassx,rmassy,rmassz)

 !   call CG_initialize(CG,CG_size,displ2,load,.FALSE.)
 !   call compute_Diagonal(precon)
 !    rmassx = 1.0_CUSTOM_REAL/precon(1,:)
  !  rmassy = 1.0_CUSTOM_REAL/precon(2,:)
  !  rmassz = 1.0_CUSTOM_REAL/precon(3,:)
   ! call CG_initialize_preconditioner(CG,rmassx,rmassy,rmassz)

end subroutine prepare_CGC

!================================================================================


subroutine prepare_CG()
    call smoother(NSPEC_AB,NGLOB_AB,displ2,displ,ibool)
 !   displ = displ2
 !   displ2 = 0.0_CUSTOM_REAL
  !  displ2 = displ2 + displ     ! this is for testing with precomputing
!    displ2 = 0.0_CUSTOM_REAL   ! this is for testing no precomputing
    displ2(:,:) = displ2(:,:)*MASK_default(:,:)*(-1.0_CUSTOM_REAL)
 !   write(*,*) "load:",load(:,1)

   ! if (abs(minval(abs(xstore-XPOSITION)+abs(ystore-YPOSITION)+abs(zstore-ZPOSITION)))<1.0e-5_CUSTOM_REAL) then
   !     Boundary_loc(1) = minloc(abs(xstore-XPOSITION)+abs(ystore-YPOSITION)+abs(zstore-ZPOSITION),1)
   !     write(*,*) "findpoint:",myrank,Boundary_loc(1)
   !     displ2(:,Boundary_loc(1)) = 0.0_CUSTOM_REAL
        call CG_initialize(CG,CG_size,displ2,load,.FALSE.,MASK_default,MASK_default)
    !else
    !    call CG_initialize(CG,CG_size,displ2,load,.FALSE.,null(),null())

    !endif
    call compute_Diagonal(precon)
    rmassx = 1.0_CUSTOM_REAL/precon(1,:)
    rmassy = 1.0_CUSTOM_REAL/precon(2,:)
    rmassz = 1.0_CUSTOM_REAL/precon(3,:)

    call CG_initialize_preconditioner(CG,rmassx,rmassy,rmassz)
end subroutine prepare_CG

!=================================================================================

subroutine get_boundary()
    integer igll,counter,ispec,i,j,k
    do counter = 1,num_abs_boundary_faces
    ispec = abs_boundary_ispec(counter)
    do igll = 1,NGLLSQUARE
    i = abs_boundary_ijk(1,igll,counter)
    j = abs_boundary_ijk(2,igll,counter)
    k = abs_boundary_ijk(3,igll,counter)
    MASK_default(:,ibool(i,j,k,ispec)) = .false.
    enddo
    enddo
end subroutine get_boundary


end subroutine iterate_time_static_solver2





