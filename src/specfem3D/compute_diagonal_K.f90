
subroutine compute_diagonal_K(iphase, &
                        NSPEC_AB,NGLOB_AB,AX,  &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz, &
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz,&
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        kappastore,mustore,jacobian,ibool, &
                        deltat, &
                        SIMULATION_TYPE,&
                        NSPEC_BOUN,&
                        num_phase_ispec_elastic,nspec_inner_elastic,nspec_outer_elastic, &
                        phase_ispec_inner_elastic)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,ONE_THIRD,FOUR_THIRDS,PI  !IMAIN WAS ADDED BY Kangchen

  use fault_solver_dynamic, only : Kelvin_Voigt_eta!, KV_direction
  use specfem_par, only : PI,FULL_ATTENUATION_SOLID, xigll, yigll, zigll, ystore
  use specfem_par_elastic, only : rmassx ,rmassy ,rmassz

  implicit none

  integer :: NSPEC_AB,NGLOB_AB

! displacement, velocity and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: AX,delta_accel

! time step
  real(kind=CUSTOM_REAL) :: deltat

! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
        kappastore,mustore,jacobian

! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy,hprimewgll_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz

! memory variables and standard linear solids for attenuation


! anisotropy
  integer :: iphase
  integer :: num_phase_ispec_elastic,nspec_inner_elastic,nspec_outer_elastic
  integer, dimension(num_phase_ispec_elastic,2) :: phase_ispec_inner_elastic

! adjoint simulations
  integer :: SIMULATION_TYPE
  integer :: NSPEC_BOUN

  ! moho kernel


! local parameters
  integer :: ispec,iglob,ispec_p,num_elements
  integer :: i,j,k,l,m,t,u,v

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl

  real(kind=CUSTOM_REAL) :: duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) :: duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  real(kind=CUSTOM_REAL) :: sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy

  real(kind=CUSTOM_REAL) :: hp1,hp2,hp3
  real(kind=CUSTOM_REAL) :: fac1,fac2,fac3

  real(kind=CUSTOM_REAL) :: lambdal,mul,lambdalplus2mul
  real(kind=CUSTOM_REAL) :: kappal

  ! local anisotropy parameters


  ! local attenuation parameters


! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc, &
    newtempx1,newtempx2,newtempx3,newtempy1,newtempy2,newtempy3,newtempz1,newtempz2,newtempz3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

  ! local C-PML absorbing boundary conditions parameters

!  write(*,*) 'entering compute forces viscoelastic noDev'
  ! choses inner/outer elements
  if( iphase == 1 ) then
    num_elements = nspec_outer_elastic
  else
    num_elements = nspec_inner_elastic
  endif
!  WRITE(*,*) num_elements
!  write(*,*) X(1,3)
  do ispec_p = 1,num_elements
!        write(*,*) ispec_p
        ! returns element id from stored element list
        ispec = phase_ispec_inner_elastic(ispec_p,iphase)

        ! adjoint simulations: moho kernel
       ! Kelvin Voigt damping: artificial viscosity around dynamic faults

        ! stores displacment values in local array
           dummyx_loc(:,:,:) = 0.0_CUSTOM_REAL
           dummyy_loc(:,:,:) = 0.0_CUSTOM_REAL
           dummyz_loc(:,:,:) = 0.0_CUSTOM_REAL

        do l=1,3
          do k=1,NGLLZ
            do j=1,NGLLY
              do i=1,NGLLX



                 iglob = ibool(i,j,k,ispec)
!                    write(*,*) MASKX(1,iglob)
                dummyx_loc(i,j,k) = 0.0_CUSTOM_REAL
                dummyy_loc(i,j,k) = 0.0_CUSTOM_REAL
                dummyz_loc(i,j,k) = 0.0_CUSTOM_REAL


          if(l==1) dummyx_loc(i,j,k) = 1.0_CUSTOM_REAL
          if(l==2) dummyy_loc(i,j,k) = 1.0_CUSTOM_REAL
          if(l==3) dummyz_loc(i,j,k) = 1.0_CUSTOM_REAL


          tempx1(:,:,:) = 0._CUSTOM_REAL
          tempx2(:,:,:) = 0._CUSTOM_REAL
          tempx3(:,:,:) = 0._CUSTOM_REAL

          tempy1(:,:,:) = 0._CUSTOM_REAL
          tempy2(:,:,:) = 0._CUSTOM_REAL
          tempy3(:,:,:) = 0._CUSTOM_REAL

          tempz1(:,:,:) = 0._CUSTOM_REAL
          tempz2(:,:,:) = 0._CUSTOM_REAL
          tempz3(:,:,:) = 0._CUSTOM_REAL


        do v=1,NGLLX

            hp1 = hprime_xx(v,i)
            tempx1(v,j,k) = tempx1(v,j,k) + dummyx_loc(i,j,k)*hp1
            tempy1(v,j,k) = tempy1(v,j,k) + dummyy_loc(i,j,k)*hp1
            tempz1(v,j,k) = tempz1(v,j,k) + dummyz_loc(i,j,k)*hp1

            !!! can merge these loops because NGLLX = NGLLY = NGLLZ
            hp2 = hprime_yy(v,j)
            tempx2(i,v,k) = tempx2(i,v,k) + dummyx_loc(i,j,k)*hp2
            tempy2(i,v,k) = tempy2(i,v,k) + dummyy_loc(i,j,k)*hp2
            tempz2(i,v,k) = tempz2(i,v,k) + dummyz_loc(i,j,k)*hp2

            !!! can merge these loops because NGLLX = NGLLY = NGLLZ
            hp3 = hprime_zz(v,k)
            tempx3(i,j,v) = tempx3(i,j,v) + dummyx_loc(i,j,k)*hp3
            tempy3(i,j,v) = tempy3(i,j,v) + dummyy_loc(i,j,k)*hp3
            tempz3(i,j,v) = tempz3(i,j,v) + dummyz_loc(i,j,k)*hp3
        enddo



        do t=1,NGLLZ
            do u=1,NGLLY
              do v=1,NGLLX

              ! get derivatives of ux, uy and uz with respect to x, y and z
              xixl = xix(v,u,t,ispec)
              xiyl = xiy(v,u,t,ispec)
              xizl = xiz(v,u,t,ispec)
              etaxl = etax(v,u,t,ispec)
              etayl = etay(v,u,t,ispec)
              etazl = etaz(v,u,t,ispec)
              gammaxl = gammax(v,u,t,ispec)
              gammayl = gammay(v,u,t,ispec)
              gammazl = gammaz(v,u,t,ispec)
              jacobianl = jacobian(v,u,t,ispec)

              duxdxl = xixl*tempx1(v,u,t) + etaxl*tempx2(v,u,t) + gammaxl*tempx3(v,u,t)
              duxdyl = xiyl*tempx1(v,u,t) + etayl*tempx2(v,u,t) + gammayl*tempx3(v,u,t)
              duxdzl = xizl*tempx1(v,u,t) + etazl*tempx2(v,u,t) + gammazl*tempx3(v,u,t)

              duydxl = xixl*tempy1(v,u,t) + etaxl*tempy2(v,u,t) + gammaxl*tempy3(v,u,t)
              duydyl = xiyl*tempy1(v,u,t) + etayl*tempy2(v,u,t) + gammayl*tempy3(v,u,t)
              duydzl = xizl*tempy1(v,u,t) + etazl*tempy2(v,u,t) + gammazl*tempy3(v,u,t)
              duzdxl = xixl*tempz1(v,u,t) + etaxl*tempz2(v,u,t) + gammaxl*tempz3(v,u,t)
              duzdyl = xiyl*tempz1(v,u,t) + etayl*tempz2(v,u,t) + gammayl*tempz3(v,u,t)
              duzdzl = xizl*tempz1(v,u,t) + etazl*tempz2(v,u,t) + gammazl*tempz3(v,u,t)


              duxdxl_plus_duydyl = duxdxl + duydyl
              duxdxl_plus_duzdzl = duxdxl + duzdzl
              duydyl_plus_duzdzl = duydyl + duzdzl
              duxdyl_plus_duydxl = duxdyl + duydxl
              duzdxl_plus_duxdzl = duzdxl + duxdzl
              duzdyl_plus_duydzl = duzdyl + duydzl



              kappal = kappastore(v,u,t,ispec)
              mul = mustore(v,u,t,ispec)



                lambdalplus2mul = kappal + FOUR_THIRDS * mul
                lambdal = lambdalplus2mul - 2.*mul

                ! compute stress sigma
                sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl
                sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl
                sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl

                sigma_xy = mul*duxdyl_plus_duydxl
                sigma_xz = mul*duzdxl_plus_duxdzl
                sigma_yz = mul*duzdyl_plus_duydzl

                 sigma_yx = sigma_xy
                 sigma_zx = sigma_xz
                 sigma_zy = sigma_yz

                 ! form dot product with test vector, non-symmetric form (which is useful in the case of PML)
                 tempx1(v,u,t) = jacobianl * (sigma_xx*xixl + sigma_yx*xiyl + sigma_zx*xizl) ! this goes to accel_x
                 tempy1(v,u,t) = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_zy*xizl) ! this goes to accel_y
                 tempz1(v,u,t) = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl) ! this goes to accel_z

                 tempx2(v,u,t) = jacobianl * (sigma_xx*etaxl + sigma_yx*etayl + sigma_zx*etazl) ! this goes to accel_x
                 tempy2(v,u,t) = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_zy*etazl) ! this goes to accel_y
                 tempz2(v,u,t) = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl) ! this goes to ccel_z

                 tempx3(v,u,t) = jacobianl * (sigma_xx*gammaxl + sigma_yx*gammayl + sigma_zx*gammazl) ! this goes to accel_x
                 tempy3(v,u,t) = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_zy*gammazl) ! this goes to accel_y
                 tempz3(v,u,t) = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl) ! this goes to accel_z

              ! stores derivatives of ux, uy and uz with respect to x, y and z

            enddo
          enddo
        enddo
              ! save strain on the Moho boundary
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
         do v = 1,NGLLX
            fac1 = hprimewgll_xx(v,i)
            newtempx1(i,j,k) = newtempx1(i,j,k) + tempx1(v,j,k)*fac1
            newtempy1(i,j,k) = newtempy1(i,j,k) + tempy1(v,j,k)*fac1
            newtempz1(i,j,k) = newtempz1(i,j,k) + tempz1(v,j,k)*fac1

            !!! can merge these loops because NGLLX = NGLLY = NGLLZ
            fac2 = hprimewgll_yy(v,j)
            newtempx2(i,j,k) = newtempx2(i,j,k) + tempx2(i,v,k)*fac2
            newtempy2(i,j,k) = newtempy2(i,j,k) + tempy2(i,v,k)*fac2
            newtempz2(i,j,k) = newtempz2(i,j,k) + tempz2(i,v,k)*fac2

            !!! can merge these loops because NGLLX = NGLLY = NGLLZ
            fac3 = hprimewgll_zz(v,k)
            newtempx3(i,j,k) = newtempx3(i,j,k) + tempx3(i,j,v)*fac3
            newtempy3(i,j,k) = newtempy3(i,j,k) + tempy3(i,j,v)*fac3
            newtempz3(i,j,k) = newtempz3(i,j,k) + tempz3(i,j,v)*fac3

        enddo
          fac1 = wgllwgll_yz(j,k)
          fac2 = wgllwgll_xz(i,k)
          fac3 = wgllwgll_xy(i,j)

          ! sum contributions from each element to the global mesh using indirect addressing
          iglob = ibool(i,j,k,ispec)

          delta_accel(1,iglob) = - fac1*newtempx1(i,j,k) - &
                                fac2*newtempx2(i,j,k) - fac3*newtempx3(i,j,k)

 !       if(l==1)  write(*,*) "ss:",fac1,fac2,fac3,newtempx1(i,j,k),newtempx2(i,j,k),newtempx3(i,j,k)

          delta_accel(2,iglob) = - fac1*newtempy1(i,j,k) - &
                                fac2*newtempy2(i,j,k) - fac3*newtempy3(i,j,k)
          delta_accel(3,iglob) = - fac1*newtempz1(i,j,k) - &
                                fac2*newtempz2(i,j,k) - fac3*newtempz3(i,j,k)

        if(l==1)  then
            AX(1,iglob) = AX(1,iglob) + delta_accel(1,iglob)
        endif

        if(l==2)   AX(2,iglob) = AX(2,iglob) + delta_accel(2,iglob)
        if(l==3)   AX(3,iglob) = AX(3,iglob) + delta_accel(3,iglob)


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
enddo

end subroutine compute_diagonal_K

