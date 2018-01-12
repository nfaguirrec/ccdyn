!!**********************************************************************************
!!    Consejo Superior de Investigaciones Científicas                              !
!!    Departamento de Física Atómica, Molecular y de Agregados                     !
!!    http://www.iff.csic.es/fama/                                                 !
!!                                                                                 !
!!    Authors:                                                                     !
!!    (2015-2015) Néstor F. Aguirre                                                !
!!                nfaguirrec@gmail.com                                             !
!!    (2015-2015) María Pilar de Lara-Castells                                     !
!!                pilar.delara.castells@csic.es                                    !
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

! http://physics.nist.gov/PhysRefData/PES/RefData/Mg_pot/Gr.html
!
! Numerical derivatives from:
! http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/central-differences/
module MgMgPotential_
	implicit none
	private
	
	public :: &
		MgMgPotential_test
	
	real(8), private, parameter :: De = 434.022_8     ! cm1
	real(8), private, parameter :: alpha = 1.04869_8  ! angs-1
	real(8), private, parameter :: Re = 3.91652_8     ! angs
		
	type, public :: MgMgPotential
		contains
			procedure :: init
			procedure :: V
			procedure :: dV
			procedure :: NdV
	end type MgMgPotential
	
	contains
	
	!>
	!! @brief Contructor
	!!
	subroutine init( this, model )
			class(MgMgPotential) :: this
			integer, optional, intent(in) :: model
			
	end subroutine init
	
	function V( this, r ) result( output )
		class(MgMgPotential), intent(in) :: this
		real(8), intent(in) :: r
		real(8) :: output
		
		integer :: i
		
		output = De*( exp(-2.0*alpha*(r-Re))-2.0*exp(-alpha*(r-Re)) )
	end function V
	
	function dV( this, r ) result( output )
		class(MgMgPotential), intent(in) :: this
		real(8), intent(in) :: r
		real(8) :: output
		
		integer :: i
		
		output = De*(2.0*alpha*exp(-alpha*(r-Re))-2.0*alpha*exp(-2.0*alpha*(r-Re)))
	end function dV
	
	function NdV( this, r, nPoints, stepSize ) result( output )
		class(MgMgPotential), intent(in) :: this
		real(8), intent(in) :: r
		integer, optional, intent(in) :: nPoints
		real(8), optional, intent(in) :: stepSize
		real(8) :: output
		
		integer :: nPointsEff
		real(8) :: h
		
		nPointsEff = 5
		if( present(nPoints) ) nPointsEff = nPoints
		
		h = 0.00001_8
		if( present(stepSize) ) h = stepSize
		
#define f(i) this.V(r+i##.0_8*h)
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
	
	subroutine MgMgPotential_test()
		real(8) :: r
		type(MgMgPotential) :: potential
		
		do r = 3.0,12.0,0.01
			write(*,"(4F15.7)") r, potential.V( r ), potential.dV( r ), potential.NdV( r )
		end do
		
	end subroutine MgMgPotential_test

end module MgMgPotential_
