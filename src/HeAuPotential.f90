!!**********************************************************************************
!!    Consejo Superior de Investigaciones Científicas                              !
!!    Departamento de Física Atómica, Molecular y de Agregados                     !
!!    http://www.iff.csic.es/fama/                                                 !
!!                                                                                 !
!!    Authors:                                                                     !
!!    (2011-2012) Néstor F. Aguirre                                                !
!!                nfaguirrec@iff.csic.es                                           !
!!                nfaguirrec@gmail.com                                             !
!!    (2011-2012) María Pilar de Lara-Castells                                     !
!!                delara@iff.csic.es                                               !
!!                                                                                 !
!!    This program is free software; you can redistribute it and/or modify         !
!!    it under the terms of the GNU General Public License as published by         !
!!    the Free Software Foundation; either version 2 of the License, or            !
!!    (at your option) any later version.                                          !
!!                                                                                 !
!!    this program is distributed in the hope that it will be useful,              !
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of               !
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                !
!!    GNU General Public License for more details.                                 !
!!                                                                                 !
!!    You should have received a copy of the GNU General Public License            !
!!    along with thisPtr program. If not, write to the Free Software Foundation,   !
!!    Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.              !
!!                                                                                 !
!!**********************************************************************************

! Potential parameters for Murrell-Sorbie funtion from:
! Theoretical characteristics of the bound states of M-X complexes (M = Cu, Ag, and Au, and X = He, Ne, and Ar)
! Xiao-Fei Tong, Chuan-Lu Yang, Yi-Peng An, Mei-Shan Wang, Xiao-Guang Ma, and De-Hua Wang
! J. Chem. Phys. 131, 244304 (2009)
!
! Numerical derivatives from:
! http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/central-differences/
module HeAuPotential_
	implicit none
	private
	
	public :: &
		HeAuPotential_test
	
	real(8), private, parameter :: tocm1 = 0.6950356_8
	
	real(8), private, parameter :: De = 14.110_8 ! cm1
	real(8), private, parameter :: Re = 4.124_8  ! angs
	real(8), private, parameter :: ai(9) = [ 1.6417_8, -0.6164_8, 0.3770_8, -0.0549_8, 0.0018_8, &
						 0.0046_8, -0.0010_8, 0.0001_8,  2.700d-6 ]  ! angs^(-i), i=pos
		
	type, public :: HeAuPotential
		contains
			procedure :: init
			procedure :: V
			procedure :: dV
			procedure :: NdV
	end type HeAuPotential
	
	contains
	
	!>
	!! @brief Contructor
	!!
	subroutine init( this, model )
			class(HeAuPotential) :: this
			integer, optional, intent(in) :: model
			
	end subroutine init
	
	function V( this, R ) result( output )
		class(HeAuPotential), intent(in) :: this
		real(8), intent(in) :: R
		real(8) :: output
		
		integer :: i
		
		output = 1.0_8
		do i=1,size(ai)
			output = output + ai(i)*(R-Re)**i
		end do
		
		output = -De*output*exp(-ai(1)*(R-Re))

! 		real(8), parameter :: amu = 1822.88853_8
! 		real(8), parameter :: angs = 1.0_8/0.52917726_8
! 		real(8), parameter :: cm1 = 1.0_8/219474.63068_8
! 		
! 		real(8) :: mass = 4.0026_8*amu
! 		real(8) :: omega = 30.0_8*cm1
! 		real(8) :: effDe = De*cm1
! 		real(8) :: effRe = Re*angs
! 		
! 		output = 0.5_8*mass*omega**2*(r*angs-effRe)**2-effDe
! 		output = output/cm1
	end function V
	
	function dV( this, R ) result( output )
		class(HeAuPotential), intent(in) :: this
		real(8), intent(in) :: R
		real(8) :: output
		
		integer :: i
		
		output = 0.0_8
		do i=1,size(ai)
			output = output + i*ai(i)*(R-Re)**(i-1)
		end do
		
		output = -De*output*exp(-ai(1)*(R-Re))-this.V(R)*ai(1)

! 		real(8), parameter :: amu = 1822.88853_8
! 		real(8), parameter :: angs = 1.0_8/0.52917726_8
! 		real(8), parameter :: cm1 = 1.0_8/219474.63068_8
! 		
! 		real(8) :: mass = 4.0026_8*amu
! 		real(8) :: omega = 30.0_8*cm1
! 		real(8) :: effDe = De*cm1
! 		real(8) :: effRe = Re*angs
! 		
! 		output = mass*omega**2*(r*angs-effRe)
! 		output = output/(cm1/angs)
	end function dV
	
	function NdV( this, R, nPoints, stepSize ) result( output )
		class(HeAuPotential), intent(in) :: this
		real(8), intent(in) :: R
		integer, optional, intent(in) :: nPoints
		real(8), optional, intent(in) :: stepSize
		real(8) :: output
		
		integer :: nPointsEff
		real(8) :: h
		
		nPointsEff = 5
		if( present(nPoints) ) nPointsEff = nPoints
		
		h = 0.00001_8
		if( present(stepSize) ) h = stepSize
		
#define f(i) this.V(R+i##.0_8*h)
		select case( nPointsEff )
			case( 3 )
				output = ( f(1)-f(-1) )/(2.0_8*h)
			case( 5 )
				output = ( f(-2)-8.0_8*f(-1)+8.0_8*f(1)-f(2) )/(12.0_8*h)
			case( 7 )
				output = ( -f(-3)+9.0_8*f(-2)-45.0_8*f(-1)+45.0_8*f(1)-9.0_8*f(2)+f(3) )/(60.0_8*h)
			case( 9 )
				output = ( 3.0_8*f(-4)-32.0_8*f(-3)+168.0_8*f(-2) &
					  -672.0_8*f(-1)+672.0_8*f(1)-168.0_8*f(2)+32.0_8*f(3)-3.0_8*f(4) )/(840.0_8*h)
			case default
				write(*,*) "### ERROR ### This formula is not implemented, it's only available nPoints=3,5,7,9"
				stop
		end select
#undef f(i)
	end function NdV
	
	subroutine HeAuPotential_test()
		real(8) :: r
		type(HeAuPotential) :: potential
		
		do r = 1.0,10.0,0.01
			write(*,"(4F15.7)") r, potential.V( r ), potential.dV( r ), potential.NdV( r )
		end do
		
	end subroutine HeAuPotential_test

end module HeAuPotential_
