subroutine compute_AX2(AX,X,MASKX,MASKAX)

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use pml_par
  use fault_solver_dynamic, only : bc_dynflt_set3d_all,SIMULATION_TYPE_DYN
  use fault_solver_kinematic, only : bc_kinflt_set_all,SIMULATION_TYPE_KIN

  implicit none

  integer:: iphase
  logical:: phase_is_inner
  integer:: iface,ispec,igll,i,j,k,iglob
  real(kind=CUSTOM_REAL),dimension(3,NGLOB_AB),intent(in) :: X
  real(kind=CUSTOM_REAL),dimension(3,NGLOB_AB),intent(out) :: AX
!  logical,dimension(3,NGLOB_AB),optional :: MASKXin,MASKAXin
  logical,dimension(3,NGLOB_AB),intent(in) :: MASKX,MASKAX
   AX(:,:)=0.0_CUSTOM_REAL

    do iphase=1,2

    !first for points on MPI interfaces
    if( iphase == 1 ) then
      phase_is_inner = .false.
    else
      phase_is_inner = .true.
    endif

! elastic term
    if(USE_DEVILLE_PRODUCTS) then
      ! uses Deville (2002) optimizations
      call compute_forces_viscoelastic_Dev_sim1(iphase)

    else
      ! no optimizations used
      call compute_forces_viscoelastic_noDev2(iphase,NSPEC_AB,NGLOB_AB, &
                        X,AX, MASKX, MASKAX,  &
                        xix2,xiy2,xiz2,etax2,etay2,etaz2,gammax2,gammay2,gammaz2, &
                        hprime_xx2,hprime_yy2,hprime_zz2, &
                        hprimewgll_xx2,hprimewgll_yy2,hprimewgll_zz2, &
                        wgllwgll_xy2,wgllwgll_xz2,wgllwgll_yz2, &
                        kappastore,mustore,jacobian2,ibool, &
                        num_phase_ispec_elastic,nspec_inner_elastic,nspec_outer_elastic, &
                        phase_ispec_inner_elastic)

    endif
   if (phase_is_inner .eqv. .false.) then
       ! sends accel values to corresponding MPI interface neighbors
       call assemble_MPI_vector_async_send(NPROC,NGLOB_AB,AX, &
               buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
               num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
               nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
               my_neighbours_ext_mesh, &
               request_send_vector_ext_mesh,request_recv_vector_ext_mesh)

    else
      ! waits for send/receive requests to be completed and assembles values
      call assemble_MPI_vector_async_w_ord(NPROC,NGLOB_AB,AX, &
                            buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh,&
                            max_nibool_interfaces_ext_mesh, &
                            nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                            request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
                            my_neighbours_ext_mesh,myrank)
endif



   enddo


end subroutine compute_AX2



