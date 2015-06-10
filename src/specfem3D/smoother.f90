
subroutine smoother(NSPEC_AB,NGLOB_AB ,AX,ibool)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,ONE_THIRD,FOUR_THIRDS,PI  !IMAIN WAS ADDED BY Kangchen

  use specfem_par, only : PI,FULL_ATTENUATION_SOLID, xigll, yigll, zigll, ystore

  implicit none


  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: AX

  integer :: NSPEC_AB, NGLOB_AB
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool





  integer :: ispec,iglob,ispec_p,num_elements
  integer :: i,j,k,l,ia
  real(kind=CUSTOM_REAL),dimension(8) :: shape3D ! local anisotropy parameters
  real(kind=CUSTOM_REAL),dimension(3,8) :: Corner_value


 do ispec_p = 1,NSPEC_AB

          Corner_value(:,1) = AX(:,ibool(1    ,1    ,1    ,ispec_p))
          Corner_value(:,2) = AX(:,ibool(NGLLX,1    ,1    ,ispec_p))
          Corner_value(:,3) = AX(:,ibool(NGLLX,NGLLY,1    ,ispec_p))
          Corner_value(:,4) = AX(:,ibool(1    ,NGLLY,1    ,ispec_p))
          Corner_value(:,5) = AX(:,ibool(1    ,1    ,NGLLZ,ispec_p))
          Corner_value(:,6) = AX(:,ibool(NGLLX,1    ,NGLLZ,ispec_p))
          Corner_value(:,7) = AX(:,ibool(NGLLX,NGLLY,NGLLZ,ispec_p))
          Corner_value(:,8) = AX(:,ibool(1    ,NGLLY,NGLLZ,ispec_p))


          do k=2,NGLLZ-1
            do j=2,NGLLY-1
              do i=2,NGLLX-1
               iglob = ibool(i,j,k,ispec_p)
        
               call eval_shape3D_single(99999,shape3D,xigll(i),yigll(j),zigll(k),8)
               AX(:,iglob) = 0.0_CUSTOM_REAL
              do ia=1,8
                AX(:,iglob) = AX(:,iglob) + shape3D(ia)*Corner_value(:,ia)
              enddo

                                      
            enddo
          enddo
        enddo
  enddo
end subroutine smoother

