subroutine compute_AX(X,AX)

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
      call compute_forces_viscoelastic_noDev(iphase,NSPEC_AB,NGLOB_AB, &
                        X,veloc,AX,load, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz, &
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        kappastore,mustore,jacobian,ibool, &
                        ATTENUATION,deltat,PML_CONDITIONS, &
                        one_minus_sum_beta,factor_common, &
                        one_minus_sum_beta_kappa,factor_common_kappa, &
                        alphaval,betaval,gammaval,&
                        NSPEC_ATTENUATION_AB,NSPEC_ATTENUATION_AB_kappa, &
                        R_trace,R_xx,R_yy,R_xy,R_xz,R_yz, &
                        epsilondev_trace,epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                        epsilondev_xz,epsilondev_yz,epsilon_trace_over_3, &
                        ANISOTROPY,NSPEC_ANISO, &
                        c11store,c12store,c13store,c14store,c15store,c16store, &
                        c22store,c23store,c24store,c25store,c26store,c33store, &
                        c34store,c35store,c36store,c44store,c45store,c46store, &
                        c55store,c56store,c66store, &
                        SIMULATION_TYPE,COMPUTE_AND_STORE_STRAIN,NSPEC_STRAIN_ONLY, &
                        NSPEC_BOUN,NSPEC2D_MOHO,NSPEC_ADJOINT, &
                        is_moho_top,is_moho_bot, &
                        dsdx_top,dsdx_bot, &
                        ispec2D_moho_top,ispec2D_moho_bot, &
                        num_phase_ispec_elastic,nspec_inner_elastic,nspec_outer_elastic, &
                        phase_ispec_inner_elastic,.false.)

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


end subroutine compute_AX


subroutine compute_XTAX(XTAX_ALL,X)

use specfem_par
use specfem_par_acoustic
use specfem_par_elastic
use specfem_par_poroelastic
use pml_par
use fault_solver_dynamic, only : bc_dynflt_set3d_all,SIMULATION_TYPE_DYN
use fault_solver_kinematic, only : bc_kinflt_set_all,SIMULATION_TYPE_KIN

implicit none

integer:: iphase,num_iteration
logical:: phase_is_inner
integer:: iface,ispec,igll,i,j,k,iglob
real(kind=CUSTOM_REAL),dimension(3,NGLOB_AB),intent(in) :: X
real(kind=CUSTOM_REAL),dimension(3,NGLOB_AB) :: AX
real(kind=CUSTOM_REAL),intent(out) :: XTAX_ALL
real(kind=CUSTOM_REAL) :: XTAX

AX(:,:)=0.0_CUSTOM_REAL
XTAX = 0.0_CUSTOM_REAL


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
call compute_forces_viscoelastic_noDev(iphase,NSPEC_AB,NGLOB_AB, &
X,veloc,AX,load, &
xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
hprime_xx,hprime_yy,hprime_zz, &
hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
kappastore,mustore,jacobian,ibool, &
ATTENUATION,deltat,PML_CONDITIONS, &
one_minus_sum_beta,factor_common, &
one_minus_sum_beta_kappa,factor_common_kappa, &
alphaval,betaval,gammaval,&
NSPEC_ATTENUATION_AB,NSPEC_ATTENUATION_AB_kappa, &
R_trace,R_xx,R_yy,R_xy,R_xz,R_yz, &
epsilondev_trace,epsilondev_xx,epsilondev_yy,epsilondev_xy, &
epsilondev_xz,epsilondev_yz,epsilon_trace_over_3, &
ANISOTROPY,NSPEC_ANISO, &
c11store,c12store,c13store,c14store,c15store,c16store, &
c22store,c23store,c24store,c25store,c26store,c33store, &
c34store,c35store,c36store,c44store,c45store,c46store, &
c55store,c56store,c66store, &
SIMULATION_TYPE,COMPUTE_AND_STORE_STRAIN,NSPEC_STRAIN_ONLY, &
NSPEC_BOUN,NSPEC2D_MOHO,NSPEC_ADJOINT, &
is_moho_top,is_moho_bot, &
dsdx_top,dsdx_bot, &
ispec2D_moho_top,ispec2D_moho_bot, &
num_phase_ispec_elastic,nspec_inner_elastic,nspec_outer_elastic, &
phase_ispec_inner_elastic,.false.)

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

do i=1,3
do j=1,NGLOB_AB
XTAX = XTAX + X(i,j)*AX(i,j)
enddo
enddo

call sum_all_all_cr(XTAX,XTAX_ALL)


end subroutine compute_XTAX


