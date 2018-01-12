!!**********************************************************************************
!!                                                                                 !
!!    Authors:                                                                     !
!!    (2016-2016) Néstor F. Aguirre                                                !
!!                nfaguirrec@gmail.com                                             !
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

! Potential parameters from:
! Davis, Jellinek and Berry JCP 86 (1987) 6456

! Numerical derivatives from:
! http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/central-differences/
module ArArPotential_
	use UnitsConverter_
	implicit none
	private
	
	public :: &
		ArArPotential_test
	
	integer, public, parameter :: JELLINEK_BERRY = 1
	
	! Parameters for JELLINEK_BERRY
	real(8), private :: De
	real(8), private :: R0
	
	type, public :: ArArPotential
		integer :: model
		
		contains
			procedure :: init
			procedure :: V
			procedure :: dV
			procedure :: NdV
	end type ArArPotential
	
	contains
	
	!>
	!! @brief Contructor
	!!
	subroutine init( this, model )
			class(ArArPotential) :: this
			integer, optional, intent(in) :: model
			
			this.model = JELLINEK_BERRY
			if( present(model) ) this.model = model
			
			select case( this.model )
				case( JELLINEK_BERRY )
! 					De = 0.000383049983726041_8  ! a.u.
! 					R0 = 6.42506882352121_8  ! a.u.
					De = 119.8*Kelvin  ! from http://dx.doi.org.sci-hub.io/10.1103/PhysRevA.11.1068
					R0 = 3.4*angs      ! from http://dx.doi.org.sci-hub.io/10.1103/PhysRevA.11.1068
			end select

	end subroutine init
	
	function V( this, R ) result( output )
		class(ArArPotential), intent(in) :: this
		real(8), intent(in) :: R
		real(8) :: output
		
		select case( this.model )
			case( JELLINEK_BERRY )
				output = 4.0_8*De*( (R0/R)**12 - (R0/R)**6 )
		end select
		
	end function V
	
	function dV( this, R ) result( output )
		class(ArArPotential), intent(in) :: this
		real(8), intent(in) :: R
		real(8) :: output
		
		select case( this.model )
			case( JELLINEK_BERRY )
				output = 4.0_8*De*( 6.0_8*R0**6/R**7 - 12.0_8*R0**12/R**13 )
		end select
		
	end function dV
	
	function NdV( this, R, nPoints, stepSize ) result( output )
		class(ArArPotential), intent(in) :: this
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
	
	subroutine ArArPotential_test()
		real(8) :: r
		type(ArArPotential) :: potential
		
		call potential.init( JELLINEK_BERRY )
		
		do r = 1.0*angs,10.0*angs,0.01*angs
			write(*,"(4E15.7)") r/angs, potential.V( r )/cm1, potential.dV( r )/(cm1/angs), potential.NdV( r )/(cm1/angs)
		end do
		
	end subroutine ArArPotential_test

end module ArArPotential_
