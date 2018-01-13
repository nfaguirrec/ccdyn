!!**********************************************************************************
!!    Consejo Superior de Investigaciones Científicas                              !
!!    Departamento de Física Atómica, Molecular y de Agregados                     !
!!    http://www.iff.csic.es/fama/                                                 !
!!                                                                                 !
!!    Authors:                                                                     !
!!    (2018-2018) Néstor F. Aguirre                                                !
!!                nfaguirrec@gmail.com                                             !
!!    (2018-2018) María Pilar de Lara-Castells                                     !
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

! Potential parameters from:
! Calculation of the Parameters of the LennardJones Potential for
!       Pairs of Identical Atoms Based on the Properties of Solid Substances
! V. P. Filippovaa, S. A. Kunavinb, and M. S. Pugachevc
! Inorganic Materials: Applied Research, 2015, Vol. 6, No. 1, pp. 1–4. 

! Numerical derivatives from:
! http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/central-differences/
module NiNiPotential_
	use UnitsConverter_
	implicit none
	private
	
	public :: &
		NiNiPotential_test
	
	integer, public, parameter :: LENNARD_JONES = 1
	
	! Parameters for LENNARD_JONES
	real(8), private :: De
	real(8), private :: R0
	
	type, public :: NiNiPotential
		integer :: model
		
		contains
			procedure :: init
			procedure :: V
			procedure :: dV
			procedure :: NdV
	end type NiNiPotential
	
	contains
	
	!>
	!! @brief Contructor
	!!
	subroutine init( this, model )
			class(NiNiPotential) :: this
			integer, optional, intent(in) :: model
			
			this.model = LENNARD_JONES
			if( present(model) ) this.model = model
			
			select case( this.model )
				case( LENNARD_JONES )
					De = 0.50722_8*eV
					R0 =  0.257366_8*nm
			end select

	end subroutine init
	
	function V( this, R ) result( output )
		class(NiNiPotential), intent(in) :: this
		real(8), intent(in) :: R
		real(8) :: output
		
		select case( this.model )
			case( LENNARD_JONES )
				output = 4.0_8*De*( (R0/R)**12 - (R0/R)**6 )
		end select
		
	end function V
	
	function dV( this, R ) result( output )
		class(NiNiPotential), intent(in) :: this
		real(8), intent(in) :: R
		real(8) :: output
		
		select case( this.model )
			case( LENNARD_JONES )
				output = 4.0_8*De*( 6.0_8*R0**6/R**7 - 12.0_8*R0**12/R**13 )
		end select
		
	end function dV
	
	function NdV( this, R, nPoints, stepSize ) result( output )
		class(NiNiPotential), intent(in) :: this
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
	
	subroutine NiNiPotential_test()
		real(8) :: r
		type(NiNiPotential) :: potential
		
		call potential.init( LENNARD_JONES )
		
		do r = 1.0*angs,10.0*angs,0.01*angs
			write(*,"(4E15.7)") r/angs, potential.V( r )/cm1, potential.dV( r )/(cm1/angs), potential.NdV( r )/(cm1/angs)
		end do
		
	end subroutine NiNiPotential_test

end module NiNiPotential_
