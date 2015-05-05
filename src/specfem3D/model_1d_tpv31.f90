!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                            July 2012
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

!--------------------------------------------------------------------------------------------------
!
! 1D Southern California model
!
! model is the standard model used in southern California:
!   Kanamori and Hadley (1975), Dreger and Helmberger (1990), Wald-Hutton,Given (1995)
!
!--------------------------------------------------------------------------------------------------

  subroutine model_1D_tpv31(xmesh,ymesh,zmesh,rho,vp,vs,qmu_atten)

! given a GLL point, returns super-imposed velocity model values

  !use create_regions_mesh_ext_par
  use specfem_par
  implicit none

  ! GLL point location
  real(kind=CUSTOM_REAL), intent(in) :: xmesh,ymesh,zmesh

  ! density, Vp and Vs
  real(kind=CUSTOM_REAL),intent(inout) :: vp,vs,rho,qmu_atten

  ! local parameters
  real(kind=CUSTOM_REAL) :: depth,x,y,z

  ! mesh point location
  x = xmesh
  y = ymesh
  z = zmesh

  ! depth in m
  depth = -zmesh

  ! assigns model parameters
  if( depth >= 10000.0 ) then
    ! moho
    vp=6.5d0
    vs=3.8d0
    rho=3.0d0
  else if( depth > 5000.0 ) then
    vp=5.75d0
    vs=3.45d0
    rho=2.72d0
  else if( depth > 2400.0 ) then
    vp=4.45d0+(5.2d0-4.45d0)*(depth-2400.0)/2600.0
    vs=2.55d0+(3.05d0-2.55d0)*(depth-2400.0)/2600.0
    rho=2.60d0+(2.62d0-2.60d0)*(depth-2400.0)/2600.0
  else
    vp=4.05d0
    vs=2.25d0
    rho=2.58d0
  endif
!  write(*,*) 'function:',depth,' ',vs 
  ! scale to standard units
  vp = vp * 1000.d0
  vs = vs * 1000.d0
  rho = rho * 1000.d0

  ! no attenuation information
  qmu_atten = 0.d0

  end subroutine model_1D_tpv31
