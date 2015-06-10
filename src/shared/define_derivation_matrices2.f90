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

  subroutine define_derivation_matrices2(xigll,yigll,zigll,wxgll,wygll,wzgll, &
         hprime_xx,hprime_yy,hprime_zz, &
         hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
         wgllwgll_xy,wgllwgll_xz,wgllwgll_yz)

  implicit none

  include "constants.h"

! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(2) :: xigll,wxgll
  double precision, dimension(2) :: yigll,wygll
  double precision, dimension(2) :: zigll,wzgll

! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(2,2) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(2,2) :: hprime_yy,hprimewgll_yy
  real(kind=CUSTOM_REAL), dimension(2,2) :: hprime_zz,hprimewgll_zz
  real(kind=CUSTOM_REAL), dimension(2,2) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(2,2) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(2,2) :: wgllwgll_yz

! function for calculating derivatives of Lagrange polynomials
  double precision, external :: lagrange_deriv_GLL

  integer i,j,k,i1,i2,j1,j2,k1,k2

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,2,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,2,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,2,GAUSSALPHA,GAUSSBETA)

! if number of points is odd, the middle abscissa is exactly ZERO
  if(mod(2,2) /= 0) xigll((2-1)/2+1) = ZERO
  if(mod(2,2) /= 0) yigll((2-1)/2+1) = ZERO
  if(mod(2,2) /= 0) zigll((2-1)/2+1) = ZERO

! distinguish between single and double precision for reals
  if(CUSTOM_REAL == SIZE_REAL) then

! calculate derivatives of the Lagrange polynomials
! and precalculate some products in double precision
! hprime(i,j) = h'_j(xigll_i) by definition of the derivation matrix
  do i1=1,2
    do i2=1,2
      hprime_xx(i2,i1) = sngl(lagrange_deriv_GLL(i1-1,i2-1,xigll,2))
      hprimewgll_xx(i2,i1) = sngl(lagrange_deriv_GLL(i1-1,i2-1,xigll,2)*wxgll(i2))
    enddo
  enddo

  do j1=1,2
    do j2=1,2
      hprime_yy(j2,j1) = sngl(lagrange_deriv_GLL(j1-1,j2-1,yigll,2))
      hprimewgll_yy(j2,j1) = sngl(lagrange_deriv_GLL(j1-1,j2-1,yigll,2)*wygll(j2))
    enddo
  enddo

  do k1=1,2
    do k2=1,2
      hprime_zz(k2,k1) = sngl(lagrange_deriv_GLL(k1-1,k2-1,zigll,2))
      hprimewgll_zz(k2,k1) = sngl(lagrange_deriv_GLL(k1-1,k2-1,zigll,2)*wzgll(k2))
    enddo
  enddo

  do i=1,2
    do j=1,2
      wgllwgll_xy(i,j) = sngl(wxgll(i)*wygll(j))
    enddo
  enddo

  do i=1,2
    do k=1,2
      wgllwgll_xz(i,k) = sngl(wxgll(i)*wzgll(k))
    enddo
  enddo

  do j=1,2
    do k=1,2
      wgllwgll_yz(j,k) = sngl(wygll(j)*wzgll(k))
    enddo
  enddo

  else  ! double precision version

! calculate derivatives of the Lagrange polynomials
! and precalculate some products in double precision
! hprime(i,j) = h'_j(xigll_i) by definition of the derivation matrix
  do i1=1,2
    do i2=1,2
      hprime_xx(i2,i1) = lagrange_deriv_GLL(i1-1,i2-1,xigll,2)
      hprimewgll_xx(i2,i1) = lagrange_deriv_GLL(i1-1,i2-1,xigll,2)*wxgll(i2)
    enddo
  enddo

  do j1=1,2
    do j2=1,2
      hprime_yy(j2,j1) = lagrange_deriv_GLL(j1-1,j2-1,yigll,2)
      hprimewgll_yy(j2,j1) = lagrange_deriv_GLL(j1-1,j2-1,yigll,2)*wygll(j2)
    enddo
  enddo

  do k1=1,2
    do k2=1,2
      hprime_zz(k2,k1) = lagrange_deriv_GLL(k1-1,k2-1,zigll,2)
      hprimewgll_zz(k2,k1) = lagrange_deriv_GLL(k1-1,k2-1,zigll,2)*wzgll(k2)
    enddo
  enddo

  do i=1,2
    do j=1,2
      wgllwgll_xy(i,j) = wxgll(i)*wygll(j)
    enddo
  enddo

  do i=1,2
    do k=1,2
      wgllwgll_xz(i,k) = wxgll(i)*wzgll(k)
    enddo
  enddo

  do j=1,2
    do k=1,2
      wgllwgll_yz(j,k) = wygll(j)*wzgll(k)
    enddo
  enddo

  endif

  end subroutine define_derivation_matrices2

