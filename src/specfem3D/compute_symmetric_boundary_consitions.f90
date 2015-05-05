!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

! for elastic solver



  subroutine compute_symmetric_boundary_conditions(NSPEC_AB,NGLOB_AB, &
                                    ibool,accel, &
                                    free_surface_normal,free_surface_ijk,free_surface_ispec, &
                                    num_free_surface_faces)

! updates acceleration with ocean load term:
! approximates ocean-bottom continuity of pressure & displacement for longer period waves (> ~20s ),
! assuming incompressible fluid column above bathymetry ocean bottom

  implicit none

  include 'constants.h'

  integer :: NSPEC_AB,NGLOB_AB

  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB_AB),intent(inout) :: accel

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

  ! free surface
  integer :: num_free_surface_faces
  real(kind=CUSTOM_REAL) :: free_surface_normal(NDIM,NGLLSQUARE,num_free_surface_faces)
  integer :: free_surface_ijk(3,NGLLSQUARE,num_free_surface_faces)
  integer :: free_surface_ispec(num_free_surface_faces)

! local parameters
  real(kind=CUSTOM_REAL) :: nx,ny,nz
  real(kind=CUSTOM_REAL) :: force_normal_comp
  integer :: i,j,k,ispec,iglob
  integer :: igll,iface
!  logical,dimension(NGLOB_AB) :: updated_dof_ocean_load

  !   initialize the updates
!  updated_dof_ocean_load(:) = .false.

  ! for surface elements exactly at the top of the model (ocean bottom)
  do iface = 1,num_free_surface_faces

    ispec = free_surface_ispec(iface)
    do igll = 1, NGLLSQUARE
      i = free_surface_ijk(1,igll,iface)
      j = free_surface_ijk(2,igll,iface)
      k = free_surface_ijk(3,igll,iface)

      ! get global point number
      iglob = ibool(i,j,k,ispec)

      ! only update once
!      if(.not. updated_dof_ocean_load(iglob)) then
!
!        ! get normal
!        nx = free_surface_normal(1,igll,iface)
!        ny = free_surface_normal(2,igll,iface)
!        nz = free_surface_normal(3,igll,iface)
!
!        ! make updated component of right-hand side
!        ! we divide by rmass() which is 1 / M
!        ! we use the total force which includes the Coriolis term above
!        force_normal_comp = accel(1,iglob)*nx / rmassx(iglob) &
!                            + accel(2,iglob)*ny / rmassy(iglob) &
!                            + accel(3,iglob)*nz / rmassz(iglob)
!
!        accel(1,iglob) = accel(1,iglob) &
!          + (rmass_ocean_load(iglob) - rmassx(iglob)) * force_normal_comp * nx
!        accel(2,iglob) = accel(2,iglob) &
!          + (rmass_ocean_load(iglob) - rmassy(iglob)) * force_normal_comp * ny
!        accel(3,iglob) = accel(3,iglob) &
!          + (rmass_ocean_load(iglob) - rmassz(iglob)) * force_normal_comp * nz
!
!        ! done with this point
!        updated_dof_ocean_load(iglob) = .true.

!      endif
       accel(3,iglob)=0._CUSTOM_REAL
       
    enddo ! igll
  enddo ! iface

  end subroutine compute_symmetric_boundary_conditions
!
!-------------------------------------------------------------------------------------------------
!


