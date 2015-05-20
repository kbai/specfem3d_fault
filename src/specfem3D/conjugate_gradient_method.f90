! Base module for kinematic and dynamic fault solvers
!
! Authors:
! Percy Galvez, Surendra Somala, Jean-Paul Ampuero

module conjugate_gradient

  use constants

  implicit none

  type conjugate_gradient
    real(kind=CUSTOM_REAL), dimension(:,:),   pointer :: X=>null()
    real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: Residue , Pdirection
  end type conjugate_gradient
!!!!! DK DK  private

contains

!---------------------------------------------------------------------

subroutine CG_initialize(CG_variable,Displ,Load)

!! DK DK use type(bc_dynandkinflt_type) instead of class(fault_type) for compatibility with some current compilers
  use specfem_par

  type(conjugate_gradient),intent(in) :: CG_variable
  real(kind=CUSTOM_REAL),dimension(3,NGLOB_AB),intent(in) :: Displ,Load

  CG_variable%X => Displ

  allocate(CG_variable%Pdirection(3,NGLOB_AB))
  allocate(CG_variable%Residue(3,NGLOB_AB))
!  allocate(Norm_old(3,NGLOB_AB))
 CG_variable%Residue(:,:) = Load(:,:)
 CG_variable%Pdirection(:,:) = Load(:,:)
!  call Ax(CG_variable)

end subroutine CG_initialize
!==================================================================

subroutine update_value_direction(CG_variable)

real(kind=CUSTOM_REAL) :: PTW , Norm_old , Norm_new , alpha , beta
real(kind=CUSTOM_REAL),dimension(3,NGLOB_AB) :: W
type(conjugate_gradient),intent(in) :: CG_variable

call Ax(W,CG_variable%Pdirection)
! w=Ap
call Vector_Multiplication(CG_variable%Pdirection,W,PTW)
! ptw = p'Ap
call Vector_Multiplication(Norm_old,CG_variable%Residue,CG_variable%Residue)

alpha = Norm_old/PTW

CG_variable%X(:,:) = CG_variable%X(:,:) + alpha * CG_variable%Pdirection

CG_variable%Residue(:,:) = CG_variable%Residue(:,:) - alpha*W(:,:)

call Vector_Multiplication(Norm_new,CG_variable%Residue,CG_variable%Residue)

beta = Norm_new/Norm_old

CG_variable%Pdirection(:,:) = CG_variable%Residue(:,:) + beta * CG_variable%Pdirection(:,:)


end subroutine update_value_direction



subroutine Ax(AX ,X)

real(kind=CUSTOM_REAL),dimension(3,NGLOB_AB),intent(in) :: X
real(kind=CUSTOM_REAL),dimension(3,NGLOB_AB),intent(out) :: AX

call compute_forces_viscoelastic(AX , X)    ! THIS IS AN INTERFACE AND CAN BE REPLACED BY ANY SUBROUTINE

end subroutine Ax


subroutine Vector_Multiplication(SUM_ALL,A,B)

real(kind=CUSTOM_REAL),dimension(3,NGLOB_AB),intent(in) :: A , B
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
