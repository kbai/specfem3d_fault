
subroutine compute_forces_viscoelastic_noDev2(iphase, &
                        NSPEC_AB,NGLOB_AB,X,AX, MASKX, MASKAX, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz, &
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz,&
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        kappastore,mustore,jacobian,ibool, &
                        num_phase_ispec_elastic,nspec_inner_elastic,nspec_outer_elastic, &
                        phase_ispec_inner_elastic)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,N_SLS,SAVE_MOHO_MESH,ONE_THIRD,FOUR_THIRDS,PI  !IMAIN WAS ADDED BY Kangchen

  use specfem_par, only : PI,FULL_ATTENUATION_SOLID, xigll, yigll, zigll, ystore
  use specfem_par_elastic, only : rmassx ,rmassy ,rmassz

  implicit none

  integer :: NSPEC_AB,NGLOB_AB

! displacement, velocity and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: X,AX,delta_accel
  logical, dimension(NDIM,NGLOB_AB),intent(in) :: MASKX, MASKAX
! time step

! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool
  real(kind=CUSTOM_REAL), dimension(2,2,2,NSPEC_AB) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
        kappastore,mustore

! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(2,2) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(2,2) :: hprime_yy,hprimewgll_yy
  real(kind=CUSTOM_REAL), dimension(2,2) :: hprime_zz,hprimewgll_zz
  real(kind=CUSTOM_REAL), dimension(2,2) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(2,2) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(2,2) :: wgllwgll_yz

! memory variables and standard linear solids for attenuation

! anisotropy


  integer :: iphase
  integer :: num_phase_ispec_elastic,nspec_inner_elastic,nspec_outer_elastic
  integer, dimension(num_phase_ispec_elastic,2) :: phase_ispec_inner_elastic



  ! moho kernel

  integer :: ispec,iglob,ispec_p,num_elements
  integer :: i,j,k,l

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl

  real(kind=CUSTOM_REAL) :: duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) :: duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  real(kind=CUSTOM_REAL) :: sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy

  real(kind=CUSTOM_REAL) :: hp1,hp2,hp3
  real(kind=CUSTOM_REAL) :: fac1,fac2,fac3

  real(kind=CUSTOM_REAL) :: lambdal,mul,lambdalplus2mul
  real(kind=CUSTOM_REAL) :: kappal




! local parameters
  real(kind=CUSTOM_REAL), dimension(2,2,2) :: dummyx_loc,dummyy_loc,dummyz_loc, &
    newtempx1,newtempx2,newtempx3,newtempy1,newtempy2,newtempy3,newtempz1,newtempz2,newtempz3
  real(kind=CUSTOM_REAL), dimension(2,2,2) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3
  integer, dimension(2) :: M

  real(kind=CUSTOM_REAL) duxdxl_att,duxdyl_att,duxdzl_att,duydxl_att
  real(kind=CUSTOM_REAL) duydyl_att,duydzl_att,duzdxl_att,duzdyl_att,duzdzl_att
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl_att,duzdxl_plus_duxdzl_att,duzdyl_plus_duydzl_att



  integer :: fault_direction
  ! local C-PML absorbing boundary conditions parameters

  ! choses inner/outer elements
  if( iphase == 1 ) then
    num_elements = nspec_outer_elastic
  else
    num_elements = nspec_inner_elastic
  endif

  M = (/1,NGLLX/)
!  WRITE(*,*) num_elements
!  write(*,*) X(1,3)
  do ispec_p = 1,num_elements
        ! returns element id from stored element list
        ispec = phase_ispec_inner_elastic(ispec_p,iphase)

        ! adjoint simulations: moho kernel

       ! Kelvin Voigt damping: artificial viscosity around dynamic faults

        ! stores displacment values in local array



          do k=1,2
            do j=1,2
              do i=1,2


                 iglob = ibool(M(i),M(j),M(k),ispec)
!                    write(*,*) MASKX(1,iglob)
                dummyx_loc(i,j,k) = X(1,iglob)
                dummyy_loc(i,j,k) = X(2,iglob)
                dummyz_loc(i,j,k) = X(3,iglob)
                dummyx_loc(i,j,k) = dummyx_loc(i,j,k)*MASKX(1,iglob) * (-1.0e0_CUSTOM_REAL)
                dummyy_loc(i,j,k) = dummyy_loc(i,j,k)*MASKX(2,iglob) * (-1.0e0_CUSTOM_REAL)
                dummyz_loc(i,j,k) = dummyz_loc(i,j,k)*MASKX(3,iglob) * (-1.0e0_CUSTOM_REAL)



              enddo
            enddo
          enddo



        ! use first order Taylor expansion of displacement for local storage of stresses
        ! at this current time step, to fix attenuation in a consistent way


     do k=1,2
        do j=1,2
          do i=1,2

          tempx1(i,j,k) = 0._CUSTOM_REAL
          tempx2(i,j,k) = 0._CUSTOM_REAL
          tempx3(i,j,k) = 0._CUSTOM_REAL

          tempy1(i,j,k) = 0._CUSTOM_REAL
          tempy2(i,j,k) = 0._CUSTOM_REAL
          tempy3(i,j,k) = 0._CUSTOM_REAL

          tempz1(i,j,k) = 0._CUSTOM_REAL
          tempz2(i,j,k) = 0._CUSTOM_REAL
          tempz3(i,j,k) = 0._CUSTOM_REAL

          do l=1,2
            hp1 = hprime_xx(i,l)
            tempx1(i,j,k) = tempx1(i,j,k) + dummyx_loc(l,j,k)*hp1
            tempy1(i,j,k) = tempy1(i,j,k) + dummyy_loc(l,j,k)*hp1
            tempz1(i,j,k) = tempz1(i,j,k) + dummyz_loc(l,j,k)*hp1

            !!! can merge these loops because NGLLX = NGLLY = NGLLZ
            hp2 = hprime_yy(j,l)
            tempx2(i,j,k) = tempx2(i,j,k) + dummyx_loc(i,l,k)*hp2
            tempy2(i,j,k) = tempy2(i,j,k) + dummyy_loc(i,l,k)*hp2
            tempz2(i,j,k) = tempz2(i,j,k) + dummyz_loc(i,l,k)*hp2

            !!! can merge these loops because NGLLX = NGLLY = NGLLZ
            hp3 = hprime_zz(k,l)
            tempx3(i,j,k) = tempx3(i,j,k) + dummyx_loc(i,j,l)*hp3
            tempy3(i,j,k) = tempy3(i,j,k) + dummyy_loc(i,j,l)*hp3
            tempz3(i,j,k) = tempz3(i,j,k) + dummyz_loc(i,j,l)*hp3
          enddo

              xixl = xix(i,j,k,ispec)
              xiyl = xiy(i,j,k,ispec)
              xizl = xiz(i,j,k,ispec)
              etaxl = etax(i,j,k,ispec)
              etayl = etay(i,j,k,ispec)
              etazl = etaz(i,j,k,ispec)
              gammaxl = gammax(i,j,k,ispec)
              gammayl = gammay(i,j,k,ispec)
              gammazl = gammaz(i,j,k,ispec)
              jacobianl = jacobian(i,j,k,ispec)

              duxdxl = xixl*tempx1(i,j,k) + etaxl*tempx2(i,j,k) + gammaxl*tempx3(i,j,k)
              duxdyl = xiyl*tempx1(i,j,k) + etayl*tempx2(i,j,k) + gammayl*tempx3(i,j,k)
              duxdzl = xizl*tempx1(i,j,k) + etazl*tempx2(i,j,k) + gammazl*tempx3(i,j,k)

              duydxl = xixl*tempy1(i,j,k) + etaxl*tempy2(i,j,k) + gammaxl*tempy3(i,j,k)
              duydyl = xiyl*tempy1(i,j,k) + etayl*tempy2(i,j,k) + gammayl*tempy3(i,j,k)
              duydzl = xizl*tempy1(i,j,k) + etazl*tempy2(i,j,k) + gammazl*tempy3(i,j,k)

              duzdxl = xixl*tempz1(i,j,k) + etaxl*tempz2(i,j,k) + gammaxl*tempz3(i,j,k)
              duzdyl = xiyl*tempz1(i,j,k) + etayl*tempz2(i,j,k) + gammayl*tempz3(i,j,k)
              duzdzl = xizl*tempz1(i,j,k) + etazl*tempz2(i,j,k) + gammazl*tempz3(i,j,k)

              ! stores derivatives of ux, uy and uz with respect to x, y and z


              ! save strain on the Moho boundary

              duxdxl_plus_duydyl = duxdxl + duydyl
              duxdxl_plus_duzdzl = duxdxl + duzdzl
              duydyl_plus_duzdzl = duydyl + duzdzl
              duxdyl_plus_duydxl = duxdyl + duydxl
              duzdxl_plus_duxdzl = duzdxl + duxdzl
              duzdyl_plus_duydzl = duzdyl + duydzl



              kappal = kappastore(M(i),M(j),M(k),ispec)
              mul = mustore(M(i),M(j),M(k),ispec)

              ! attenuation



  ! isotropic case
                lambdalplus2mul = kappal + FOUR_THIRDS * mul
                lambdal = lambdalplus2mul - 2.*mul

                ! compute stress sigma
                sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl
                sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl
                sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl

                sigma_xy = mul*duxdyl_plus_duydxl
                sigma_xz = mul*duzdxl_plus_duxdzl
                sigma_yz = mul*duzdyl_plus_duydzl
     !           write(IMAIN,*) 'sigma_xy',sigma_xy,'mul',mul   !added by Kangchen
               ! ANISOTROPY

              ! subtract memory variables if attenuation


                 ! define symmetric components of sigma
                 sigma_yx = sigma_xy
                 sigma_zx = sigma_xz
                 sigma_zy = sigma_yz

                 ! form dot product with test vector, non-symmetric form (which is useful in the case of PML)
                 tempx1(i,j,k) = jacobianl * (sigma_xx*xixl + sigma_yx*xiyl + sigma_zx*xizl) ! this goes to accel_x
                 tempy1(i,j,k) = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_zy*xizl) ! this goes to accel_y
                 tempz1(i,j,k) = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl) ! this goes to accel_z

                 tempx2(i,j,k) = jacobianl * (sigma_xx*etaxl + sigma_yx*etayl + sigma_zx*etazl) ! this goes to accel_x
                 tempy2(i,j,k) = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_zy*etazl) ! this goes to accel_y
                 tempz2(i,j,k) = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl) ! this goes to ccel_z

                 tempx3(i,j,k) = jacobianl * (sigma_xx*gammaxl + sigma_yx*gammayl + sigma_zx*gammazl) ! this goes to accel_x
                 tempy3(i,j,k) = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_zy*gammazl) ! this goes to accel_y
                 tempz3(i,j,k) = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl) ! this goes to accel_z

         enddo
      enddo
    enddo


    ! second double-loop over GLL to compute all the terms
    do k=1,2
       do j=1,2
          do i=1,2

          newtempx1(i,j,k) = 0._CUSTOM_REAL
          newtempy1(i,j,k) = 0._CUSTOM_REAL
          newtempz1(i,j,k) = 0._CUSTOM_REAL

          newtempx2(i,j,k) = 0._CUSTOM_REAL
          newtempy2(i,j,k) = 0._CUSTOM_REAL
          newtempz2(i,j,k) = 0._CUSTOM_REAL

          newtempx3(i,j,k) = 0._CUSTOM_REAL
          newtempy3(i,j,k) = 0._CUSTOM_REAL
          newtempz3(i,j,k) = 0._CUSTOM_REAL
!          write(*,*) 'BKpoint 1'
          do l=1,2
            fac1 = hprimewgll_xx(l,i)
            newtempx1(i,j,k) = newtempx1(i,j,k) + tempx1(l,j,k)*fac1
            newtempy1(i,j,k) = newtempy1(i,j,k) + tempy1(l,j,k)*fac1
            newtempz1(i,j,k) = newtempz1(i,j,k) + tempz1(l,j,k)*fac1

            !!! can merge these loops because NGLLX = NGLLY = NGLLZ
            fac2 = hprimewgll_yy(l,j)
            newtempx2(i,j,k) = newtempx2(i,j,k) + tempx2(i,l,k)*fac2
            newtempy2(i,j,k) = newtempy2(i,j,k) + tempy2(i,l,k)*fac2
            newtempz2(i,j,k) = newtempz2(i,j,k) + tempz2(i,l,k)*fac2

            !!! can merge these loops because NGLLX = NGLLY = NGLLZ
            fac3 = hprimewgll_zz(l,k)
            newtempx3(i,j,k) = newtempx3(i,j,k) + tempx3(i,j,l)*fac3
            newtempy3(i,j,k) = newtempy3(i,j,k) + tempy3(i,j,l)*fac3
            newtempz3(i,j,k) = newtempz3(i,j,k) + tempz3(i,j,l)*fac3
          enddo

          fac1 = wgllwgll_yz(j,k)
          fac2 = wgllwgll_xz(i,k)
          fac3 = wgllwgll_xy(i,j)

          ! sum contributions from each element to the global mesh using indirect addressing
          iglob = ibool(M(i),M(j),M(k),ispec)
          delta_accel(1,iglob) = - fac1*newtempx1(i,j,k) - &
                                fac2*newtempx2(i,j,k) - fac3*newtempx3(i,j,k)
          delta_accel(2,iglob) = - fac1*newtempy1(i,j,k) - &
                                fac2*newtempy2(i,j,k) - fac3*newtempy3(i,j,k)
          delta_accel(3,iglob) = - fac1*newtempz1(i,j,k) - &
                                fac2*newtempz2(i,j,k) - fac3*newtempz3(i,j,k)

          AX(1,iglob) = AX(1,iglob) + delta_accel(1,iglob) * MASKAX(1,iglob) * (-1.0e0_CUSTOM_REAL)
          AX(2,iglob) = AX(2,iglob) + delta_accel(2,iglob) * MASKAX(2,iglob) * (-1.0e0_CUSTOM_REAL)
          AX(3,iglob) = AX(3,iglob) + delta_accel(3,iglob) * MASKAX(3,iglob) * (-1.0e0_CUSTOM_REAL)


 !         X(1,iglob) = X(1,iglob) + load(1,iglob) + delta_accel(1,iglob)*deltat*deltat*rmassx(iglob)
 !         X(2,iglob) = X(2,iglob) + load(2,iglob) + delta_accel(2,iglob)*deltat*deltat*rmassy(iglob)
 !         X(3,iglob) = X(3,iglob) + load(3,iglob) + delta_accel(3,iglob)*deltat*deltat*rmassz(iglob)
        enddo
     enddo
  enddo
! write(*,*) 'BKpoint2 '

          !  update memory variables based upon the Runge-Kutta scheme
  enddo  ! spectral element loop
!write(*,*) 'finish the loop in single processor!'

end subroutine compute_forces_viscoelastic_noDev2

