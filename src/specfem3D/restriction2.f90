
subroutine restriction2(NSPEC_AB,NGLOB_AB,X,AX, MASKX, MASKAX, &
                        xigll,yigll,zigll,&
                        wxgll,wygll,wzgll)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,N_SLS,SAVE_MOHO_MESH,ONE_THIRD,FOUR_THIRDS,PI,NGLLSQUARE  !IMAIN WAS ADDED BY Kangchen
  use fault_solver_common, only : ibool1,jacobian2Dw
  use fault_solver_qstatic, only : faults
  use specfem_par , only : myrank

  implicit none

  integer :: NSPEC_AB,NGLOB_AB

! displacement, velocity and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: X,AX,delta_accel
  logical, dimension(NDIM,NGLOB_AB),intent(in) :: MASKX, MASKAX
  real(kind=CUSTOM_REAL), dimension(NGLLX) :: xigll,yigll,zigll,wxgll,wygll,wzgll

 integer :: e ,i1, j1 ,i2, j2, iglob1, iglob2, igll
 integer , dimension(NGLLX,NGLLY) :: ij



  real(kind=CUSTOM_REAL),dimension(NDIM) :: newtemp

  integer, dimension(2) :: M,P


  M = (/1,NGLLX/)
  P = (/1,-1/)
  write(*,*) myrank,"entering loop!"
  write(*,*) myrank,faults(1)%nspec
  if(faults(1)%nspec>0)  then
  igll = 1
  do i1 = 1,NGLLX
    do j1 = 1,NGLLY
        ij(i1,j1)=igll
        igll = igll + 1
        write(*,*) "igll",igll
    enddo
enddo

do e = 1,faults(1)%nspec

  igll = 1

    do i1 = 1,NGLLX
        do j1 = 1,NGLLY
            write(*,*) "iboolsize:",size(ibool1)
            iglob1 = ibool1(ij(i1,j1),e)
            write(*,*) "iglob1",iglob1
            write(*,*) "iglob1_further:",faults(1)%ibulk1(iglob1),"size X",size(X)

            do i2 = 1,2
                do j2 = 1,2
                    iglob2 = ibool1(ij(M(i2),M(j2)),e)
                     newtemp = X(:,faults(1)%ibulk1(iglob1)) * &
                     abs((xigll(i1)-P(i2)) * (yigll(j1)-P(j2)))*0.25e0_CUSTOM_REAL*&
                     jacobian2Dw(ij(i1,j1),e)
                     AX(:,faults(1)%ibulk1(iglob2)) = AX(:,faults(1)%ibulk1(iglob2)) + newtemp
                     newtemp = X(:,faults(1)%ibulk2(iglob1)) * &
                     abs((xigll(i1)-P(i2)) * (yigll(j1)-P(j2)))*0.25e0_CUSTOM_REAL*&
                     jacobian2Dw(ij(i1,j1),e)
                     AX(:,faults(1)%ibulk2(iglob2)) = AX(:,faults(1)%ibulk2(iglob2)) + newtemp
                     write(*,*) myrank, iglob1 , iglob2
                enddo
            enddo
        enddo
    enddo
enddo
  endif
  write(*,*) myrank,"exiting the loop!"



        ! use first order Taylor expansion of displacement for local storage of stresses
        ! at this current time step, to fix attenuation in a consistent way




          !  update memory variables based upon the Runge-Kutta scheme
    ! spectral element loop
!write(*,*) 'finish the loop in single processor!'

end subroutine restriction2

