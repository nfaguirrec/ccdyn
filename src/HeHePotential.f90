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

module HeHePotential_
	implicit none
	private
	
        public :: &
                HeHePotential_test
	
	real(8), private, parameter :: tocm1 = 0.6950356
	real(8), private, parameter :: pi = acos(-1.0_8)
	
	real(8), private, parameter :: epsilon = 10.97_8
	real(8), private, parameter :: Rm = 2.9695_8
	real(8), private, parameter :: A = 189635.353_8
	real(8), private, parameter :: alpha = 10.70203539_8
	real(8), private, parameter :: beta = 1.90740649_8
	real(8), private, parameter :: c6 = 1.34687065_8
	real(8), private, parameter :: c8 = 0.41308398_8
	real(8), private, parameter :: c10 = 0.17060159_8
	real(8), private, parameter :: D = 1.4088_8
	real(8), private, parameter :: Aa = 0.0026_8
	real(8), private, parameter :: zeta1 = 1.003535949_8
	real(8), private, parameter :: zeta2 = 1.454790369_8
	
	type, public :: HeHePotential
		contains
			procedure :: init
			procedure :: V
			procedure :: dV
			procedure :: NdV
	end type HeHePotential
	
	contains
	
	!>
	!! @brief Contructor
	!!
	subroutine init( this, model )
			class(HeHePotential) :: this
			integer, optional, intent(in) :: model
			
	end subroutine init
	
	function V( this, R ) result( output )
		class(HeHePotential), intent(in) :: this
		real(8), intent(in) :: R
		real(8) :: output
		
		output = epsilon*( Vb( R/Rm ) + Va( R/Rm ) )*tocm1
	end function V
	
	function dV( this, R ) result( output )
		class(HeHePotential), intent(in) :: this
		real(8), intent(in) :: R
		real(8) :: output
		
		output = (epsilon/Rm)*( dVb( R/Rm ) + dVa( R/Rm ) )*tocm1
	end function dV
	
	function NdV( this, R ) result( output )
		class(HeHePotential), intent(in) :: this
		real(8), intent(in) :: R
		real(8) :: output
		
		real(8) :: h = 0.00001_8
		
! http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/central-differences/
		output = ( this.V(R-2.0_8*h)-8.0_8*this.V(R-h)+8.0_8*this.V(R+h)-this.V(R+2.0_8*h) )/(12.0_8*h)
! 		output = ( -this.V(R-3.0_8*h)+9.0_8*this.V(R-2.0_8*h)-45.0_8*this.V(R-h) &
! 			   +45.0_8*this.V(R+h)-9.0_8*this.V(R+2.0_8*h)+this.V(R+3.0_8*h) )/(60.0_8*h)
	end function NdV
	
	function Va( zeta ) result( output )
		real(8), intent(in) :: zeta
		real(8) :: output
		
		if( zeta1 <= zeta .and. zeta <= zeta2  ) then
			output = Aa*( sin( 2.0_8*pi*(zeta-zeta1)/(zeta2-zeta1) - 0.5_8*pi) + 1.0_8 )
		else
			output = 0.0_8
		end if
	end function Va
	
	function dVa( zeta ) result( output )
		real(8), intent(in) :: zeta
		real(8) :: output
		
		if( zeta1 <= zeta .and. zeta <= zeta2  ) then
			output = (2.0_8*pi*Aa/(zeta2-zeta1))*sin( 2.0_8*pi*(zeta-zeta1)/(zeta2-zeta1) )
		else
			output = 0.0_8
		end if
	end function dVa
	
	function Vb( zeta ) result( output )
		real(8), intent(in) :: zeta
		real(8) :: output
		
		output = A*exp( -alpha*zeta-beta*zeta**2.0_8 ) &
			   -( c6/zeta**6.0_8 + c8/zeta**8.0_8 + c10/zeta**10.0_8 )*F(zeta)
	end function Vb
	
	function dVb( zeta ) result( output )
		real(8), intent(in) :: zeta
		real(8) :: output
		
		output = -A*(alpha+2.0_8*beta*zeta)*exp( -alpha*zeta-beta*zeta**2.0_8 ) &
			   +( 6.0_8*c6/zeta**7.0_8 + 8.0_8*c8/zeta**9.0_8 + 10.0_8*c10/zeta**11.0_8 )*F(zeta) &
			   -( c6/zeta**6.0_8 + c8/zeta**8.0_8 + c10/zeta**10.0_8 )*dF(zeta)
	end function dVb
	
	function F( zeta ) result( output )
		real(8), intent(in) :: zeta
		real(8) :: output
		
		if( zeta <= D ) then
			output = exp( -(D/zeta-1.0_8)**2.0_8 )
		else
			output = 1.0_8
		end if
	end function F
	
	function dF( zeta ) result( output )
		real(8), intent(in) :: zeta
		real(8) :: output
		
		if( zeta <= D ) then
			output = 2.0_8*D*(D/zeta-1.0_8)*exp( -(D/zeta-1.0_8)**2.0_8 )/zeta**2.0_8
		else
			output = 0.0_8
		end if
	end function dF
	
! 	function V3( r1, r2, r3 ) result( output )
! 		real(8), intent(in) :: r1, r2, r3
! 		real(8) :: output
! 		
! 		real(8) :: Q1, Q2, Q3
! 		
! 		Q1 = (r1+r2+r3)/sqrt(3.0_8)
! 		Q2 = (r2-r3)/sqrt(2.0_8)
! 		Q3 = (2.0_8*r1-r2-r3)/sqrt(6.0_8)
! 		
! 		
! 	end function V3
	
	subroutine HeHePotential_test()
		real(8) :: r
		type(HeHePotential) :: potential
		
! 		call potential.init()
				
		do r = 1.0,20.0,0.05
			write(*,"(3F15.7)") r, potential.V( r ), potential.dV( r )
		end do
		
		write(*,*) ""
		write(*,*) ""
		
! 		call potential.destroy()
		
	end subroutine HeHePotential_test

end module HeHePotential_
