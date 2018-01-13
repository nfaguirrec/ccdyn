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
! Electronic States of the NiO Molecule
! Stephen P. Walch and W. A. Goddard
! J. Am. Chem. Soc. 100:5 / March I, I978
! Electronic state: X^3Sigma^-(1)


! Numerical derivatives from:
! http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/central-differences/
module NiOPotential_
	use UnitsConverter_
	implicit none
	private
	
	public :: &
		NiOPotential_test
	
	integer, public, parameter :: MORSE = 1
	
	! Parameters for MORSE
	real(8), private :: De
	real(8), private :: Re
	real(8), private :: we
	real(8), private :: rMass
	real(8), private :: alpha

	type, public :: NiOPotential
		integer :: model
		
		contains
			procedure :: init
			procedure :: V
			procedure :: dV
			procedure :: NdV
	end type NiOPotential
	
	contains
	
	!>
	!! @brief Contructor
	!!
	subroutine init( this, model )
			class(NiOPotential) :: this
			integer, optional, intent(in) :: model
			
			real(8) :: massNi, massO
			
! 			massNi = 57.9353_8
! 			massO = 15.9949_8
			
			massNi = 58.0_8
			massO = 16.0_8
			
			this.model = MORSE
			if( present(model) ) this.model = model
			
			select case( this.model )
				case( MORSE )
					De = 3.95_8*eV
					Re = 1.60*angs
					we = 841.0_8*cm1
					rMass = ( massNi*massO/(massNi+massO) )*amu
					alpha = we*sqrt(0.5*rMass/De)
			end select

	end subroutine init
	
	function V( this, R ) result( output )
		class(NiOPotential), intent(in) :: this
		real(8), intent(in) :: R
		real(8) :: output
		
		real(8) :: damp, morsev, pauli, disp
		
		select case( this.model )
			case( MORSE )
				output = De*( exp(-2.0_8*alpha*(R*angs-Re)) - 2.0_8*exp(-alpha*(R*angs-Re)) )
		end select
		
		output = output
	end function V
	
	function dV( this, R ) result( output )
		class(NiOPotential), intent(in) :: this
		real(8), intent(in) :: R
		real(8) :: output
		
		select case( this.model )
			case( MORSE )
				output = 2.0_8*De*alpha*( exp(alpha*(Re-R*angs))-exp(2*alpha*(Re-R*angs)) )
		end select
		
 		output = output
	end function dV
	
	function NdV( this, R, nPoints, stepSize ) result( output )
		class(NiOPotential), intent(in) :: this
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
	
	subroutine NiOPotential_test()
		real(8) :: r
		type(NiOPotential) :: potential
		
		call potential.init( MORSE )
		
		do r = 1.0,10.0,0.01
			write(*,"(4E15.7)") r, potential.V( r ), potential.dV( r ), potential.NdV( r )
		end do
		
	end subroutine NiOPotential_test

end module NiOPotential_
