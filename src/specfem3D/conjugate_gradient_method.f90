! Base module for kinematic and dynamic fault solvers
!
! Authors:
! Percy Galvez, Surendra Somala, Jean-Paul Ampuero

module conjugate_gradient

  use constants
  use specfem_par
  implicit none


  type CG_Vector
    integer :: NDIM,NELE
  end type

  type CG_data
    logical :: COARSE
    real(kind=CUSTOM_REAL), dimension(:,:), pointer :: X=>null(),L=>null()
    real(kind=CUSTOM_REAL), dimension(:,:), pointer :: Residue=>null() , Pdirection=>null(), M=>null(), M2=>null()
    logical , dimension(:,:), pointer :: MASKX=>null() , MASKAX=>null()
!
    type(CG_Vector) :: CG_size !the dimensions and element numbers in the vector
    real(kind=CUSTOM_REAL) :: Norm_old=0.0_CUSTOM_REAL


  end type CG_data
!!!!! DK DK  private
  logical,dimension(:),ALLOCATABLE  :: MPI_repeat



contains

!---------------------------------------------------------------------

subroutine CG_initialize(CG,CG_dim,Displ,Load,MG,X_setfalse,AX_setfalse)

!! DK DK use type(bc_dynandkinflt_type) instead of class(fault_type) for compatibility with some current compilers
  use specfem_par
  implicit none
  integer :: ier
  type(CG_data),intent(inout) :: CG
  real(kind=CUSTOM_REAL),target,dimension(:,:),intent(in) :: Displ,Load
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: W,tmp
  logical , dimension(:,:) :: X_setfalse,AX_setfalse
  type(CG_Vector),intent(in) :: CG_dim
  logical :: MG

  CG%X=>Displ
  CG%L=>Load
  CG%CG_size = CG_dim

 ! write(*,*) CG%CG_size%NELE,'internal of the program'
  allocate(W(3,NGLOB_AB))
  allocate(tmp(3,NGLOB_AB))
  if(.not. associated(CG%MASKAX)) then
  allocate(CG%MASKX(CG%CG_size%NDIM,CG%CG_size%NELE))
  allocate(CG%MASKAX(CG%CG_size%NDIM,CG%CG_size%NELE))
  CG%MASKX(:,:) = X_setfalse
  CG%MASKAX(:,:) = AX_setfalse
  endif
  allocate(CG%Pdirection(CG%CG_size%NDIM,CG%CG_size%NELE))
  allocate(CG%Residue(CG%CG_size%NDIM,CG%CG_size%NELE))
  allocate(CG%M(CG%CG_size%NDIM,CG%CG_size%NELE))
  allocate(CG%M2(CG%CG_size%NDIM,CG%CG_size%NELE))
!  allocate(CG%MPI_repeat(CG%CG_size%NELE))
!  call CG_mask(CG,X_setfalse,AX_setfalse)
  CG%COARSE = MG

  write(*,*) 'maxX:',maxval(abs(CG%X))
  call Axx(W,CG%X,CG%MASKX,CG%MASKAX,CG%COARSE)
  write(*,*) 'check right:',maxval(abs(W)),maxval(abs(Load))
  CG%Residue = Load - W
  CG%Pdirection= CG%Residue
  CG%M(:,:)=1.0_CUSTOM_REAL
!  call Prepare_MPI(CG%MPI_repeat)
  call Vector_Multiplication_Weight(CG%Norm_old,CG%Residue,CG%Residue,CG%M)
  call compute_diagonal2(CG%M2)
 ! call P(tmp,CG%Residue,CG%MASKX,CG%MASKAX,CG%M2)
  !call Vector_Multiplication(CG%Norm_old,CG%Residue,tmp)
  !CG%Pdirection= tmp
  write(*,*) "normold:",CG%Norm_old
end subroutine CG_initialize
!=================================================================
subroutine Prepare_MPI()


    use specfem_par
    implicit none

    integer :: iinterface, ipoin

     allocate(MPI_repeat(NGLOB_AB))
     MPI_repeat(:) = .true.
      ! partition border copy into the buffer
    do iinterface = 1, num_interfaces_ext_mesh
    if(my_neighbours_ext_mesh(iinterface) > myrank) then
      do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
        MPI_repeat(ibool_interfaces_ext_mesh(ipoin,iinterface)) = .false.
      enddo
      endif
    enddo
    write(*,*) MPI_repeat
end subroutine Prepare_MPI



subroutine CG_initialize_preconditioner(CG,Mx,My,Mz)

  type(CG_data),intent(inout) :: CG
  real(kind=CUSTOM_REAL),target,dimension(:),intent(in) :: Mx,My,Mz
  CG%M(1,:) = Mx
  CG%M(2,:) = My
  CG%M(3,:) = Mz
  call Vector_Multiplication_Weight(CG%Norm_old,CG%Residue,CG%Residue,CG%M)
  CG%Pdirection = CG%Pdirection * CG%M
  write(*,*) "precon:" , CG%Norm_old
end subroutine CG_initialize_preconditioner
!========================================================

!==================================================================
subroutine CG_mask(CG,X_setfalse,AX_setfalse)

  implicit none

  type(CG_data),intent(inout) :: CG
  integer,dimension(:),intent(in):: X_setfalse,AX_setfalse



  CG%MASKX(:,X_setfalse) = .false.
  CG%MASKAX(:,AX_setfalse) = .false.

end subroutine
!====================================================================
subroutine Jacobi_method(CG)
   use specfem_par
   implicit none

   type(CG_data),intent(inout) :: CG
   real(kind=CUSTOM_REAL),dimension(3,NGLOB_AB) :: q,z
   real(kind=CUSTOM_REAL),dimension(3,NGLOB_AB) :: Dx
   real(kind=CUSTOM_REAL) :: alpha,rz,zq


   z = CG%Residue*CG%M

   write(*,*)  "minM:",maxval(abs(CG%M)) ,maxloc(abs(CG%M))
   write(*,*)  "minmaxR",maxval(abs(CG%Residue)),maxval(abs(CG%Residue))
   write(*,*)  "maxDx:",maxval(abs(Dx)) , maxloc(abs(Dx))

   call Axx(q,z,CG%MASKX,CG%MASKAX,CG%COARSE)
   call Vector_Multiplication(rz,CG%Residue,z)
   call Vector_Multiplication(zq,z,q)
   alpha = rz/zq
   CG%X = CG%X + alpha * z
   CG%Residue = CG%Residue - alpha * q

   call Vector_Multiplication_Weight(CG%Norm_old,CG%Residue,CG%Residue,CG%M)

end subroutine Jacobi_method

subroutine update_value_direction(CG)

  use specfem_par
  implicit none

  real(kind=CUSTOM_REAL) :: PTW , Norm_new , alpha , beta
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: W,tmp
  type(CG_data),intent(inout) :: CG

    allocate(W(CG%CG_size%NDIM,CG%CG_size%NELE))
    allocate(tmp(CG%CG_size%NDIM,CG%CG_size%NELE))
   ! write(*,*) size(W),size(CG%Pdirection)
    call Axx(W,CG%Pdirection,CG%MASKX,CG%MASKAX,CG%COARSE)
   ! write(*,*) maxval(W),minval(W)

    call Vector_Multiplication(PTW,CG%Pdirection,W)
!    write(*,*) PTW,CG%Norm_old
    alpha = CG%Norm_old/PTW

    write(*,*) "alpha=",alpha
    CG%X(:,:) = CG%X(:,:) + alpha * CG%Pdirection

    CG%Residue(:,:) = CG%Residue(:,:) - alpha*W(:,:)

    call Vector_Multiplication_Weight(Norm_new,CG%Residue,CG%Residue,CG%M)
    write(*,*) "nmn:",Norm_new
   ! write(*,*) "mpirepeat:",CG%MPI_repeat
  !  call P(tmp,CG%Residue,CG%MASKX,CG%MASKAX,CG%M2)

   ! call Vector_Multiplication(Norm_new,CG%Residue,tmp)

    beta = Norm_new/CG%Norm_old

    CG%Norm_old = Norm_new

    CG%Pdirection = CG%Residue*CG%M + beta * CG%Pdirection
   ! CG%Pdirection = tmp + beta * CG%Pdirection

    call Axx(W,CG%X,CG%MASKX,CG%MASKAX,CG%COARSE)

    W=-W+CG%L

!call Vector_Multiplication(Norm_old,W,W)

!    write(*,*) 'true residue=',maxval(abs(W))
 !   write(*,*) 'pseudo residue=',maxval(abs(CG%Residue))
  !  deallocate(W)
end subroutine update_value_direction

!================================================================================================

subroutine Axx(AX ,X, A ,B ,COARSE)

use specfem_par

implicit none

    real(kind=CUSTOM_REAL),dimension(3,NGLOB_AB),intent(in) :: X
    real(kind=CUSTOM_REAL),dimension(3,NGLOB_AB),intent(out) :: AX
    logical,dimension(3,NGLOB_AB),intent(in) :: A
    logical,dimension(3,NGLOB_AB),intent(in) :: B
    logical :: COARSE
    if(COARSE) then
    call compute_AX2(AX , X ,A,B)    ! THIS IS AN INTERFACE AND CAN BE REPLACED BY ANY SUBROUTINE
    else
        call compute_AX(AX,X,A,B)
!        write(*,*) "fine one"
    endif
    Ax=-Ax
end subroutine Axx
!==============================================================================
subroutine Vector_Multiplication(SUM_ALL,A,B)

use specfem_par
implicit none

    real(kind=CUSTOM_REAL),dimension(:,:),intent(in) :: A , B
    real(kind=CUSTOM_REAL) :: SUM
    real(kind=CUSTOM_REAL),intent(inout) :: SUM_ALL


    integer :: i,j

    SUM = 0.0_CUSTOM_REAL
    SUM_ALL = 0.0_CUSTOM_REAL
    do i=1,3
      do j=1,NGLOB_AB
        if(MPI_repeat(j)) SUM = SUM + A(i,j)*B(i,j)
      enddo
    enddo

    call sum_all_all_cr(SUM,SUM_ALL)

end subroutine Vector_Multiplication
!==========================================================================================
subroutine Vector_Multiplication_Weight(SUM_ALL,A,B,Weight)

use specfem_par
implicit none

    real(kind=CUSTOM_REAL),dimension(:,:),intent(in) :: A , B , Weight
    real(kind=CUSTOM_REAL) :: SUM
    real(kind=CUSTOM_REAL),intent(inout) :: SUM_ALL

    integer :: i,j

    SUM = 0.0_CUSTOM_REAL
    SUM_ALL = 0.0_CUSTOM_REAL
    do i=1,3
      do j=1,NGLOB_AB
        if(MPI_repeat(j)) SUM = SUM + A(i,j)*B(i,j)*Weight(i,j)
      enddo
    enddo

    call sum_all_all_cr(SUM,SUM_ALL)

end subroutine Vector_Multiplication_Weight
!=======================================================================================
subroutine Updateresidue(CG, delta_r)

  real(kind=CUSTOM_REAL),target,dimension(:,:),intent(in) :: delta_r
  type(CG_data),intent(inout) :: CG

  CG%Residue = CG%Residue + delta_r
  CG%Pdirection = CG%Residue


  call Vector_Multiplication_Weight(CG%Norm_old,CG%Residue,CG%Residue,CG%M)


end subroutine Updateresidue
!===============================================================================
subroutine move_x(CG,dX)

    real(kind=CUSTOM_REAL),dimension(:,:),intent(in) :: dX
    real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: tmp
    type(CG_data),intent(inout) :: CG

    allocate(tmp(CG%CG_size%NDIM,CG%CG_size%NELE))
    call Axx(tmp,dX,CG%MASKX,CG%MASKAX,CG%COARSE)
    CG%Residue = CG%Residue - tmp
    call Vector_Multiplication_Weight(CG%Norm_old,CG%Residue,CG%Residue,CG%M)
    call Reinitialize(CG)
end subroutine move_x
!========================================================================
subroutine resetresidue(CG, delta_r)

  real(kind=CUSTOM_REAL),target,dimension(:,:),intent(in) :: delta_r
  type(CG_data),intent(inout) :: CG

  CG%Residue =  delta_r
  CG%Pdirection = CG%Residue


  call Vector_Multiplication_Weight(CG%Norm_old,CG%Residue,CG%Residue,CG%M)


end subroutine resetresidue
!==================================================================
subroutine Reinitialize(CG)

  type(CG_data),intent(inout) :: CG

  CG%Pdirection = CG%Residue*CG%M
  call Vector_Multiplication_Weight(CG%Norm_old,CG%Residue,CG%Residue,CG%M)
end subroutine Reinitialize
!===========================================================================================
subroutine P(PX , X , A ,B,M2 )

use specfem_par

implicit none

    real(kind=CUSTOM_REAL),dimension(3,NGLOB_AB),intent(in) :: X,M2
    real(kind=CUSTOM_REAL),dimension(3,NGLOB_AB),intent(out) :: PX
    real(kind=CUSTOM_REAL),dimension(3,NGLOB_AB) :: AX
    logical,dimension(3,NGLOB_AB),intent(in) :: A,B
    integer :: counter

    PX(:,:)=0.0_CUSTOM_REAL
    call compute_AX(PX,X,A,B)
    !call restriction_call(PX,X,A,B)
    !do counter = 1,100
    !call compute_AX2(AX,PX ,A,B)
    !Ax=-Ax
    !PX=(X-AX)/max(M2,1.0_CUSTOM_REAL)+PX
    !enddo
    ! call smoother(NSPEC_AB,NGLOB_AB,PX,ibool)
    !PX=-PX
    PX=X
end subroutine P
!===============================================================================
subroutine Destruction(CG)

  type(CG_data),intent(inout) :: CG

  deallocate(CG%MASKX)
  deallocate(CG%MASKAX)
  deallocate(CG%Pdirection)
  deallocate(CG%Residue)

end subroutine Destruction
!===============================================================================
subroutine compute_Anorm(PTW,CG)

    type(CG_data) :: CG
    real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: W
    real(kind=CUSTOM_REAL)  ::  PTW
    allocate(W(CG%CG_size%NDIM,CG%CG_size%NELE))
    call Axx(W,CG%X,CG%MASKX,CG%MASKAX,CG%COARSE)
    call Vector_Multiplication(PTW,CG%X,W)

end subroutine compute_Anorm
!================================================================================
double precision function max_all(CGvariable)
        real(kind=CUSTOM_REAL),dimension(:,:) :: CGvariable
        real(kind=CUSTOM_REAL) :: local_max,global_max
        local_max = maxval(CGvariable)
        call max_all_all_cr(local_max,global_max)
        max_all = global_max
        return
end function max_all
!================================================================================
double precision function min_all(CGvariable)
        real(kind=CUSTOM_REAL),dimension(:,:) :: CGvariable
        real(kind=CUSTOM_REAL) :: local_min,global_min
        local_min = minval(CGvariable)
        call min_all_all_cr(local_min,global_min)
        min_all = global_min
        return
end function min_all

end module conjugate_gradient
