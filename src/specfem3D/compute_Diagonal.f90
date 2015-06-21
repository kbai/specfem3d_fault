subroutine compute_Diagonal(AX)

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  implicit none

  integer:: iphase
  logical:: phase_is_inner
  integer:: iface,ispec,igll,i,j,k,iglob
  real(kind=CUSTOM_REAL),dimension(3,NGLOB_AB),intent(out) :: AX
!  logical,dimension(3,NGLOB_AB),optional :: MASKXin,MASKAXin
!  write(*,*) present(MASKXin)
!  if(present(MASKXin))  then
!     MASKX = MASKXin
!  else
!     allocate(MASKX(3,NGLOB_AB))
!      write(*,*) 'set truth of the mask'
!      MASKX(:,:)=.true.
!  endif
!  write(*,*) present(MASKAXin),'maskin defined?'
!  if(present(MASKAXin))  then
!       MASKAX = MASKAXin
!  else
!     allocate(MASKAX(3,NGLOB_AB))
!     MASKAX(:,:) = .true.
!  endif
!
! ! write(*,*) MASKX
!
!  write(*,*) 'entering computeAX'
!  !write(*,*) X
    AX(:,:)=0.0_CUSTOM_REAL

    do iphase=1,2
    write(*,*) iphase
    !first for points on MPI interfaces
    if( iphase == 1 ) then
      phase_is_inner = .false.
    else
      phase_is_inner = .true.
    endif

! elastic term
     ! no optimizations used
      call compute_diagonal_K(iphase,NSPEC_AB,NGLOB_AB, &
                        AX,  &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz, &
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        kappastore,mustore,jacobian,ibool, &
                        deltat, &
                        SIMULATION_TYPE, &
                        NSPEC_BOUN, &
                       num_phase_ispec_elastic,nspec_inner_elastic,nspec_outer_elastic, &
                        phase_ispec_inner_elastic)
   write(*,*) "Diagonal finished!"
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
   AX=-AX    ! note that the force is negative


end subroutine compute_Diagonal


