! Base module for kinematic and dynamic fault solvers
!
! Authors:
! Percy Galvez, Surendra Somala, Jean-Paul Ampuero

module conjugate_gradient
 
  use constants

  implicit none

  type CG_Vector
    integer :: NDIM,NELE   
  end type

  type CG_data
    real(kind=CUSTOM_REAL), dimension(:,:), pointer :: X=>null(),L=>null()
    real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: Residue , Pdirection
    logical , dimension(:,:), allocatable :: MASKX , MASKAX
    type(CG_Vector) :: CG_size !the dimensions and element numbers in the vector
  end type CG_data
!!!!! DK DK  private

contains

!---------------------------------------------------------------------

subroutine CG_initialize(CG_variable,CG_dim,Displ,Load)

!! DK DK use type(bc_dynandkinflt_type) instead of class(fault_type) for compatibility with some current compilers
  use specfem_par
  implicit none
  integer :: ier
  type(CG_data),intent(out) :: CG_variable
  real(kind=CUSTOM_REAL),target,dimension(:,:),intent(in) :: Displ,Load
  type(CG_Vector),intent(in) :: CG_dim
  CG_variable%X=>Displ
  CG_variable%L=>Load
  CG_variable%CG_size = CG_dim
!  write(*,*) CG_variable%CG_size%NELE,'internal of the program'
  allocate(CG_variable%MASKX(CG_variable%CG_size%NDIM,CG_variable%CG_size%NELE))
  allocate(CG_variable%MASKAX(CG_variable%CG_size%NDIM,CG_variable%CG_size%NELE))
!  write(*,*) 'ier',ier
  allocate(CG_variable%Pdirection(CG_variable%CG_size%NDIM,CG_variable%CG_size%NELE))
!  write(*,*),'size',size(CG_variable%Pdirection)
  allocate(CG_variable%Residue(CG_variable%CG_size%NDIM,CG_variable%CG_size%NELE))
!  allocate(Norm_old(CG_variable%CG_size%NDIM,CG_variable%CG_size%NELE))
!  write(*,*) 'ier',ier,size(CG_variable%Pdirection),size(CG_variable%MASKX)
  CG_variable%Residue(:,:) = Load(:,:)
  CG_variable%Pdirection(:,:) = Load(:,:)
  CG_variable%MASKX(:,:) =.true.
  CG_variable%MASKAX(:,:) = .true.
!  call Ax(CG_variable)
!  write(*,*) CG_variable%MASKX
end subroutine CG_initialize
!==================================================================
subroutine CG_mask(CG_variable,X_setfalse,AX_setfalse)
 
  implicit none

  type(CG_data),intent(inout) :: CG_variable
  
  integer,dimension(:),intent(in):: X_setfalse,AX_setfalse

  CG_variable%MASKX(:,X_setfalse) = .false.
  CG_variable%MASKAX(:,X_setfalse) = .false.

end subroutine
!====================================================================

subroutine update_value_direction(CG_variable)

  use specfem_par
  implicit none

  real(kind=CUSTOM_REAL) :: PTW , Norm_old , Norm_new , alpha , beta
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: W
  type(CG_data),intent(inout) :: CG_variable

  !write(*,*) 'started to call the AXX'
    write(*,*) size(CG_variable%Pdirection)
    write(*,*) NGLOB_AB

    allocate(W(CG_variable%CG_size%NDIM,CG_variable%CG_size%NELE))
    call Axx(W,CG_variable%Pdirection,CG_variable%MASKX,CG_variable%MASKAX)
  ! w=Ap
    write(*,*) maxval(W),minval(W)

    call Vector_Multiplication(PTW,CG_variable%Pdirection,W)
    ! ptw = p'Ap
    write(*,*) 'calculate PTW=',PTW
    call Vector_Multiplication(Norm_old,CG_variable%Residue,CG_variable%Residue)
    !write(*,*) 'calculate Vector Multiplication'
    alpha = Norm_old/PTW

    CG_variable%X(:,:) = CG_variable%X(:,:) + alpha * CG_variable%Pdirection

    CG_variable%Residue(:,:) = CG_variable%Residue(:,:) - alpha*W(:,:)

    call Vector_Multiplication(Norm_new,CG_variable%Residue,CG_variable%Residue)

    beta = Norm_new/Norm_old

    CG_variable%Pdirection(:,:) = CG_variable%Residue(:,:) + beta * CG_variable%Pdirection(:,:)

    call Axx(W,CG_variable%X,CG_variable%MASKX,CG_variable%MASKAX)

    W=-W+CG_variable%L

!call Vector_Multiplication(Norm_old,W,W)

    write(*,*) 'true residue=',maxval(abs(W))

end subroutine update_value_direction

!================================================================================================

subroutine Axx(AX ,X, A ,B)

use specfem_par

implicit none

    real(kind=CUSTOM_REAL),dimension(3,NGLOB_AB),intent(in) :: X
    real(kind=CUSTOM_REAL),dimension(3,NGLOB_AB),intent(out) :: AX
    logical,dimension(3,NGLOB_AB),intent(in) :: A
    logical,dimension(3,NGLOB_AB),intent(in) :: B
    
    call compute_AX(AX , X ,A,B)    ! THIS IS AN INTERFACE AND CAN BE REPLACED BY ANY SUBROUTINE
    Ax=-Ax
end subroutine Axx


subroutine Vector_Multiplication(SUM_ALL,A,B)

use specfem_par
implicit none

    real(kind=CUSTOM_REAL),dimension(:,:),intent(in) :: A , B
    real(kind=CUSTOM_REAL) :: SUM
    real(kind=CUSTOM_REAL),intent(out) :: SUM_ALL
    integer :: i,j

    SUM = 0.0_CUSTOM_REAL
    SUM_ALL = 0.0_CUSTOM_REAL
    do i=1,3
      do j=1,NGLOB_AB
        SUM = SUM + A(i,j)*B(i,j)
      enddo
    enddo

    call sum_all_all_cr(SUM,SUM_ALL)

end subroutine Vector_Multiplication

end module conjugate_gradient
