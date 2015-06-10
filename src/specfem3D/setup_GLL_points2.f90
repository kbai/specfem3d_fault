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
!
! United States and French Government Sponsorship Acknowledged.

  subroutine setup_GLL_points2()

  use specfem_par
  implicit none
  integer :: i,j,k,ier

  if(myrank == 0) then
    write(IMAIN,*) '******************************************'
    write(IMAIN,*) 'There is a total of ',NPROC,' slices'
    write(IMAIN,*) '******************************************'
    write(IMAIN,*)
  endif

! set up GLL points, weights and derivation matrices for reference element (between -1,1)
  call define_derivation_matrices2(xigll2,yigll2,zigll2,wxgll2,wygll2,wzgll2, &
                                hprime_xx2,hprime_yy2,hprime_zz2, &
                                hprimewgll_xx2,hprimewgll_yy2,hprimewgll_zz2, &
                                wgllwgll_xy2,wgllwgll_xz2,wgllwgll_yz2)

 end subroutine setup_GLL_points2

