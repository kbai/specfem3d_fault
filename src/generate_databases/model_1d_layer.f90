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

  subroutine model_1D_layer(xmesh,ymesh,zmesh,rho,vp,vs,qmu_atten)

! given a GLL point, returns super-imposed velocity model values

  use create_regions_mesh_ext_par
! use specfem_par
  implicit none

  ! GLL point location
  double precision,intent(in) :: xmesh,ymesh,zmesh

  ! density, Vp and Vs
  real(kind=CUSTOM_REAL) :: vp,vs,rho,qmu_atten,alpha

  ! local parameters
  real(kind=CUSTOM_REAL):: depth,x,y,z

  integer :: counter

  real(kind=CUSTOM_REAL),dimension(11,4) :: TABLE


!  TABLE=reshape((/(/0,500,1000,1600,2400,3600,5000,9000,11000,15000,60000/),&
!(/2200,3000,3600,4400,4800,5250,5500,5750,6100,6300,6300/),&
!(/1050,1400,1950,2500,2800,3100,3250,3450,3600,3700,3700/),&
!(/2200,2450,2550,2600,2600,2620,2650,2720,2750,2900,2900/)/),(/11,4/))
  TABLE=reshape((/(/0,4000,16000,30000,42000,45000,80000,80000,80000,80000,80000/),&
(/5440,6250,6530,6800,7500,8110,8110,8110,8110,8110,8110/),&
(/3000,3450,3600,3900,4300,4490,4490,4490,4490,4490,4490/),&
(/2500,2600,2700,2900,2900,3300,3300,3300,3300,3300,3300/)/),(/11,4/))
!  TABLE=reshape((/(/0,2399,2401,4999,5001,9999,10001,15000,16000,17000,60000/),&
!(/4050,4050,4450,5200,5750,5750,6500,6500,6500,6500,6500/),&
!(/2250,2250,2550,3050,3450,3450,3800,3800,3800,3800,3800/),&
!(/2580,2580,2600,2620,2720,2720,3000,3000,3000,3000,3000/)/),(/11,4/))
 

  ! mesh point location
  x = xmesh
  y = ymesh
  z = zmesh

  ! depth in m
  depth = -zmesh
!  write(*,*) zmesh,TABLE(:,1)
  ! assigns model parameters
  do counter=1,size(TABLE,1)
    if (depth>=TABLE(counter,1) .and. depth<=TABLE(counter+1,1)) then
     alpha=(depth-TABLE(counter,1))/(TABLE(counter+1,1)-TABLE(counter,1));
     Vp=(1.0-alpha)*TABLE(counter,2)+(alpha)*TABLE(counter+1,2)
     Vs=(1.0-alpha)*TABLE(counter,3)+(alpha)*TABLE(counter+1,3)
     Rho=(1.0-alpha)*TABLE(counter,4)+(alpha)*TABLE(counter+1,4)
     exit
    endif
  enddo 
!  if( depth >= 10000.0 ) then
!    ! moho
!    vp=6.5d0
!    vs=3.8d0
!    rho=3.0d0
!  else if( depth > 5000.0 ) then
!    vp=5.75d0
!    vs=3.45d0
!    rho=2.72d0
!  else if( depth > 2400.0 ) then
!    vp=4.45d0+(5.2d0-4.45d0)*(depth-2400.0)/2600.0
!    vs=2.55d0+(3.05d0-2.55d0)*(depth-2400.0)/2600.0
!    rho=2.60d0+(2.62d0-2.60d0)*(depth-2400.0)/2600.0
!  else
!    vp=4.05d0
!    vs=2.25d0
!   rho=2.58d0
!  endif

  ! scale to standard units
!  vp = vp * 1000.d0
!  vs = vs * 1000.d0
!  rho = rho * 1000.d0

  ! no attenuation information
  qmu_atten = 0.d0

  end subroutine model_1D_layer
