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

! Potential parameters from:
! Resonant photoelectron spectroscopy of Au2 via a Feshbach state using high- resolution photoelectron imaging
! Iker León, Zheng Yang, and Lai-Sheng Wang
! J. Chem. Phys. 139, 194306 (2013)


! Numerical derivatives from:
! http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/central-differences/
module AuAuPotential_
	use UnitsConverter_
	implicit none
	private
	
	public :: &
		AuAuPotential_test
	
	integer, public, parameter :: MORSE = 1
	integer, public, parameter :: GAMESS_CRENBS = 2
	
	! Parameters for MORSE
	real(8), private :: Re
	real(8), private :: we
	real(8), private :: wexe
	real(8), private :: rMass
	real(8), private :: De
	real(8), private :: alpha

	! Parameters for GAMESS_CRENBS
	real(8), private :: a    
	real(8), private :: b    
	real(8), private :: R0   
	real(8), private :: w    
	real(8), private :: c1
	
	type, public :: AuAuPotential
		integer :: model
		
		contains
			procedure :: init
			procedure :: V
			procedure :: dV
			procedure :: NdV
	end type AuAuPotential
	
	contains
	
	!>
	!! @brief Contructor
	!!
	subroutine init( this, model )
			class(AuAuPotential) :: this
			integer, optional, intent(in) :: model
			
			this.model = GAMESS_CRENBS
			if( present(model) ) this.model = model
			
			select case( this.model )
				case( MORSE )
					Re = 2.4715_8*angs
					we = 190.9_8*cm1
					wexe = 0.42_8*cm1
					rMass = (196.966569_8/2.0_8)*amu
					De = 0.25_8*we**2.0_8/wexe
					alpha = we*sqrt(0.5*rMass/De)
					
				case( GAMESS_CRENBS )
					De    = 21.5407_8
					alpha = 0.288568_8/angs
					Re    = -14.5539_8*angs
					a     = 1396.93_8
					b     = 3.72037_8/angs
					R0    = 1.67588_8*angs
					w     = 15.3537_8*angs
					c1    = 0.0521176_8*angs
			end select

	end subroutine init
	
	function V( this, R ) result( output )
		class(AuAuPotential), intent(in) :: this
		real(8), intent(in) :: R
		real(8) :: output
		
		real(8) :: damp, morsev, pauli, disp
		
		select case( this.model )
			case( MORSE )
				output = De*( exp(-2.0_8*alpha*(R*angs-Re)) - 2.0_8*exp(-alpha*(R*angs-Re)) )
			case( GAMESS_CRENBS )
				damp = 0.5_8*(1.0_8+tanh(6.0_8*(R*angs-R0)/w))
				morsev = De*( exp(-2.0_8*alpha*(R*angs-Re)) - 2.0_8*exp(-alpha*(R*angs-Re)) )
				pauli = a*exp(-b*R*angs)
				disp = abs(c1)/(R*angs)
				
				output = (1.0_8-damp)*( morsev+pauli ) - damp*disp
		end select
		
		output = output/cm1
	end function V
	
	function dV( this, R ) result( output )
		class(AuAuPotential), intent(in) :: this
		real(8), intent(in) :: R
		real(8) :: output
		
		select case( this.model )
			case( MORSE )
				output = 2.0_8*De*alpha*( exp(alpha*(Re-R*angs))-exp(2*alpha*(Re-R*angs)) )
			case( GAMESS_CRENBS )
				output = this.NdV( R )*cm1
		end select
		
 		output = output/(cm1/angs)
	end function dV
	
	function NdV( this, R, nPoints, stepSize ) result( output )
		class(AuAuPotential), intent(in) :: this
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
	
	subroutine AuAuPotential_test()
		real(8) :: r
		type(AuAuPotential) :: potential
		
		call potential.init( GAMESS_CRENBS )
		
		do r = 1.0,10.0,0.01
			write(*,"(4E15.7)") r, potential.V( r ), potential.dV( r ), potential.NdV( r )
		end do
		
	end subroutine AuAuPotential_test

end module AuAuPotential_
