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
    real(kind=CUSTOM_REAL), dimension(:,:), pointer :: Residue=>null() , Pdirection=>null(), M=>null()
    logical , dimension(:,:), pointer :: MASKX=>null() , MASKAX=>null()
    type(CG_Vector) :: CG_size !the dimensions and element numbers in the vector
    real(kind=CUSTOM_REAL) :: Norm_old=0.0_CUSTOM_REAL

  end type CG_data
!!!!! DK DK  private

contains

!---------------------------------------------------------------------

subroutine CG_initialize(CG,CG_dim,Displ,Load)

!! DK DK use type(bc_dynandkinflt_type) instead of class(fault_type) for compatibility with some current compilers
  use specfem_par
  implicit none
  integer :: ier
  type(CG_data),intent(out) :: CG
  real(kind=CUSTOM_REAL),target,dimension(:,:),intent(in) :: Displ,Load
  type(CG_Vector),intent(in) :: CG_dim

  CG%X=>Displ
  CG%L=>Load
  CG%CG_size = CG_dim

  write(*,*) CG%CG_size%NELE,'internal of the program'
  allocate(CG%MASKX(CG%CG_size%NDIM,CG%CG_size%NELE))
  allocate(CG%MASKAX(CG%CG_size%NDIM,CG%CG_size%NELE))
  write(*,*) 'ier',ier
  allocate(CG%Pdirection(CG%CG_size%NDIM,CG%CG_size%NELE))
  write(*,*),'size',size(CG%Pdirection)
  allocate(CG%Residue(CG%CG_size%NDIM,CG%CG_size%NELE))
  allocate(CG%M(CG%CG_size%NDIM,CG%CG_size%NELE))
  CG%Residue(:,:) = Load(:,:)
  CG%Pdirection(:,:) = Load(:,:)
  CG%MASKX(:,:) =.true.
  CG%MASKAX(:,:) = .true.
  CG%M(:,:)=1.0_CUSTOM_REAL
  call Vector_Multiplication_Weight(CG%Norm_old,CG%Residue,CG%Residue,CG%M)
  write(*,*) "normold:",CG%Norm_old
end subroutine CG_initialize
!=================================================================
subroutine CG_initialize_preconditioner(CG,Mx,My,Mz)

  type(CG_data),intent(inout) :: CG
  real(kind=CUSTOM_REAL),target,dimension(:),intent(in) :: Mx,My,Mz
  CG%M(1,:) = Mx
  CG%M(2,:) = My
  CG%M(3,:) = Mz
  call Vector_Multiplication_Weight(CG%Norm_old,CG%Residue,CG%Residue,CG%M)
  CG%Pdirection = CG%Pdirection * CG%M

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

subroutine update_value_direction(CG)

  use specfem_par
  implicit none

  real(kind=CUSTOM_REAL) :: PTW , Norm_new , alpha , beta
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: W
  type(CG_data),intent(inout) :: CG

    allocate(W(CG%CG_size%NDIM,CG%CG_size%NELE))
    write(*,*) size(W),size(CG%Pdirection)
    call Axx(W,CG%Pdirection,CG%MASKX,CG%MASKAX)
    write(*,*) maxval(W),minval(W)

    call Vector_Multiplication(PTW,CG%Pdirection,W)
    write(*,*) PTW,CG%Norm_old
    alpha = CG%Norm_old/PTW
    write(*,*) "alpha=",alpha
    CG%X(:,:) = CG%X(:,:) + alpha * CG%Pdirection

    CG%Residue(:,:) = CG%Residue(:,:) - alpha*W(:,:)

    call Vector_Multiplication_Weight(Norm_new,CG%Residue,CG%Residue,CG%M)

    beta = Norm_new/CG%Norm_old

    CG%Norm_old = Norm_new

    CG%Pdirection = CG%Residue*CG%M + beta * CG%Pdirection

    call Axx(W,CG%X,CG%MASKX,CG%MASKAX)

    W=-W+CG%L

!call Vector_Multiplication(Norm_old,W,W)

    write(*,*) 'true residue=',maxval(abs(W))
    write(*,*) 'pseudo residue=',maxval(abs(CG%Residue))
  !  deallocate(W)
end subroutine update_value_direction

!================================================================================================

subroutine Axx(AX ,X, A ,B)

use specfem_par

implicit none

    real(kind=CUSTOM_REAL),dimension(3,NGLOB_AB),intent(in) :: X
    real(kind=CUSTOM_REAL),dimension(3,NGLOB_AB),intent(out) :: AX
    logical,dimension(3,NGLOB_AB),intent(in) :: A
    logical,dimension(3,NGLOB_AB),intent(in) :: B

    call compute_AX2(AX , X ,A,B)    ! THIS IS AN INTERFACE AND CAN BE REPLACED BY ANY SUBROUTINE
    Ax=-Ax
end subroutine Axx


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
        SUM = SUM + A(i,j)*B(i,j)
      enddo
    enddo

    call sum_all_all_cr(SUM,SUM_ALL)

end subroutine Vector_Multiplication



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
        SUM = SUM + A(i,j)*B(i,j)*Weight(i,j)
      enddo
    enddo

    call sum_all_all_cr(SUM,SUM_ALL)

end subroutine Vector_Multiplication_Weight






subroutine Updateresidue(CG, delta_r)

  real(kind=CUSTOM_REAL),target,dimension(:,:),intent(in) :: delta_r
  type(CG_data),intent(inout) :: CG

  CG%Residue = CG%Residue + delta_r
  CG%Pdirection = CG%Residue


  call Vector_Multiplication_Weight(CG%Norm_old,CG%Residue,CG%Residue,CG%M)


end subroutine Updateresidue

subroutine Reinitialize(CG)

  type(CG_data),intent(inout) :: CG

  CG%Pdirection = CG%Residue*CG%M




end subroutine Reinitialize


subroutine Destruction(CG)

  type(CG_data),intent(inout) :: CG

  deallocate(CG%MASKX)
  deallocate(CG%MASKAX)
  deallocate(CG%Pdirection)
  deallocate(CG%Residue)

end subroutine Destruction

end module conjugate_gradient
