
subroutine smoother(NSPEC_AB,NGLOB_AB ,AX,X,ibool)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,ONE_THIRD,FOUR_THIRDS,PI  !IMAIN WAS ADDED BY Kangchen

  use specfem_par, only : PI,FULL_ATTENUATION_SOLID, xigll, yigll, zigll, wxgll, wygll,wzgll, jacobian

  implicit none


  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: AX,X

  integer :: NSPEC_AB, NGLOB_AB
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool





  integer :: ispec,iglob,ispec_p,num_elements
  integer :: i1,j1,k1,i2,j2,k2,iglob1,iglob2
  integer, dimension(2) :: M
  real(kind=CUSTOM_REAL),dimension(2) :: P
  real(kind=CUSTOM_REAL),dimension(NDIM) :: newtemp
  real(kind=CUSTOM_REAL) :: temp
  real(kind=CUSTOM_REAL),dimension(8) :: shape3D ! local anisotropy parameters
  real(kind=CUSTOM_REAL),dimension(2,2,2,3) :: Corner_value

 M = (/1,NGLLX/)
 P = (/1.0_CUSTOM_REAL,-1.0_CUSTOM_REAL/)
 do ispec = 1,NSPEC_AB

                 do k2 = 1,2
                    do j2 = 1,2
                        do i2 = 1,2
                            Corner_value(i2,j2,k2,:) = X(:,ibool(M(i2),M(j2),M(k2),ispec))
                        enddo
                    enddo
                enddo

          do k1=1,NGLLZ
    !         temp = temp * wzgll(k1)
            do j1=1,NGLLY
    !            temp = temp * wygll(j1)
              do i1=1,NGLLX

                  iglob1 = ibool(i1,j1,k1,ispec)
                  AX(:,iglob1) = 0.0_CUSTOM_REAL
                  do k2 = 1,2
                    do j2 = 1,2
                        do i2 = 1,2
                            iglob2 = ibool(M(i2),M(j2),M(k2),ispec)
                            newtemp = Corner_value(i2,j2,k2,:) * abs((xigll(i1)-P(i2)) * (yigll(j1)-P(j2))*(zigll(k1)-P(k2)))*0.125e0_CUSTOM_REAL
                            AX(:,iglob1) = AX(:,iglob1) + newtemp
                        enddo
                    enddo
                enddo

              enddo
            enddo
          enddo


  enddo
end subroutine smoother

