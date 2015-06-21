subroutine restriction_call(AX,X,MASKX,MASKAX)

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use pml_par
  use fault_solver_dynamic, only : bc_dynflt_set3d_all,SIMULATION_TYPE_DYN
  use fault_solver_kinematic, only : bc_kinflt_set_all,SIMULATION_TYPE_KIN
  use fault_solver_qstatic, only : faults

  implicit none
  integer:: iface,ispec,igll,i,j,k,iglob
  real(kind=CUSTOM_REAL),dimension(3,NGLOB_AB),intent(in) :: X
  real(kind=CUSTOM_REAL),dimension(3,NGLOB_AB),intent(out) :: AX
  real(kind=CUSTOM_REAL),dimension(3,NGLOB_AB) :: Xtmp

!  logical,dimension(3,NGLOB_AB),optional :: MASKXin,MASKAXin
  logical,dimension(3,NGLOB_AB),intent(in) :: MASKX,MASKAX
   AX(:,:)=0.0_CUSTOM_REAL
   Xtmp(:,:)=0.0_CUSTOM_REAL
   do i = 1,3
   Xtmp(i,faults(1)%ibulk1) = X(i,faults(1)%ibulk1)/faults(1)%B
   Xtmp(i,faults(1)%ibulk2) = X(i,faults(1)%ibulk2)/faults(1)%B
   enddo


      call restriction2(NSPEC_AB,NGLOB_AB, &
                        Xtmp,AX, MASKX, MASKAX,  &
                        xigll,yigll,zigll,&
                        wxgll,wygll,wzgll)
       ! sends accel values to corresponding MPI interface neighbors
       !call assemble_MPI_vector_async_send(NPROC,NGLOB_AB,AX, &
       !        buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
       !        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
       !        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
       !        my_neighbours_ext_mesh, &
       !        request_send_vector_ext_mesh,request_recv_vector_ext_mesh)
      ! waits for send/receive requests to be completed and assembles values
      !call assemble_MPI_vector_async_w_ord(NPROC,NGLOB_AB,AX,  &
      !                      buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh,&
      !                      max_nibool_interfaces_ext_mesh, &
      !                      nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
      !                      request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
      !                      my_neighbours_ext_mesh,myrank)
      call assemble_MPI_vector_blocking_ord(NPROC,NGLOB_AB,AX, &
                            num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                            nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                            my_neighbours_ext_mesh,myrank)






write(*,*) "maxval:",MAXVAL(abs(AX))
end subroutine restriction_call



