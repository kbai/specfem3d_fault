
subroutine restriction(NSPEC_AB,NGLOB_AB,X ,AX,ibool)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,ONE_THIRD,FOUR_THIRDS,PI  !IMAIN WAS ADDED BY Kangchen

  use specfem_par, only : PI,FULL_ATTENUATION_SOLID, xigll, yigll, zigll, ystore

  implicit none


  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: AX
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(out) :: X
  integer :: NSPEC_AB, NGLOB_AB
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool





  integer :: ispec,iglob,ispec_p,num_elements
  integer :: i,j,k,l,ia
 ! real(kind=CUSTOM_REAL),dimension(8) :: shape3D ! local anisotropy parameters
 ! real(kind=CUSTOM_REAL),dimension(3,8) :: Corner_value

 X(:,:) = 0.0_CUSTOM_REAL

 do ispec_p = 1,NSPEC_AB
 ! write(*,*) ibool(1    ,1    ,1    ,ispec_p)
 X(:,ibool(1    ,1    ,1    ,ispec_p))= AX(:,ibool(1    ,1    ,1    ,ispec_p))
 X(:,ibool(NGLLX,1    ,1    ,ispec_p))= AX(:,ibool(NGLLX,1    ,1    ,ispec_p))
 X(:,ibool(NGLLX,NGLLY,1    ,ispec_p))= AX(:,ibool(NGLLX,NGLLY,1    ,ispec_p))
 X(:,ibool(1    ,NGLLY,1    ,ispec_p))= AX(:,ibool(1    ,NGLLY,1    ,ispec_p))
 X(:,ibool(1    ,1    ,NGLLZ,ispec_p))= AX(:,ibool(1    ,1    ,NGLLZ,ispec_p))
 X(:,ibool(NGLLX,1    ,NGLLZ,ispec_p))= AX(:,ibool(NGLLX,1    ,NGLLZ,ispec_p))
 X(:,ibool(NGLLX,NGLLY,NGLLZ,ispec_p))= AX(:,ibool(NGLLX,NGLLY,NGLLZ,ispec_p))
 X(:,ibool(1    ,NGLLY,NGLLZ,ispec_p))= AX(:,ibool(1    ,NGLLY,NGLLZ,ispec_p))

  enddo
end subroutine restriction

