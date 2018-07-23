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
! Asynchronous multicanonical basin hopping method and its application to cobalt nanoclusters
! Lixin Zhan, Jeff Z. Y. Chen, Wing-Ki Liu, and S. K. Lai
! J. Chem. Phys. 122, 244707 (2005)
!
! Reference Clusters from:
! http://www-wales.ch.cam.ac.uk/~wales/CCD/CoCCD/cobalt.html

! Numerical derivatives from:
! http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/central-differences/
module CoCoPotential_
	use UnitsConverter_
	use IntegerList_
	implicit none
	private
	
	public :: &
		CoCoPotential_test
	
	integer, public, parameter :: MORSE = 1
	integer, public, parameter :: GUPTA = 2
	
	! Parameters for MORSE
	real(8), private :: De
	real(8), private :: Re
	real(8), private :: we
	real(8), private :: rMass
	real(8), private :: alpha
	
	! Parameters for GUPTA
	real(8), private :: zeta
	real(8), private :: q
	real(8), private :: xi
	real(8), private :: p
	real(8), private :: r0

	type, public :: CoCoPotential
		integer :: model
		
		contains
			procedure :: init
			procedure :: V
			procedure :: dV
			procedure :: NdV
			procedure :: Vnl
			procedure :: dVnl
			procedure, private, NOPASS :: beta
			procedure, private, NOPASS :: dbeta
	end type CoCoPotential
	
	contains
	
	!>
	!! @brief Contructor
	!!
	subroutine init( this, model )
			class(CoCoPotential) :: this
			integer, optional, intent(in) :: model
			
			real(8) :: massCo
			
			massCo = 58.9332_8
			
			this.model = MORSE
			if( present(model) ) this.model = model
			
			select case( this.model )
				case( MORSE )
					De = 3.95_8*eV
					Re = 1.60*angs
					we = 841.0_8*cm1
					rMass = ( massCo*massCo/(massCo+massCo) )*amu
					alpha = we*sqrt(0.5*rMass/De)
				case( GUPTA )
					zeta = 2*1.4880_8*eV
					q = 2.286_8
					xi = 2*0.0950_8*eV
					p = 11.604_8
					r0 = 2.497*angs
			end select

	end subroutine init
	
	!>
	!! @brief Returns the value of the potential at R in a.u.
	!!
	function V( this, R ) result( output )
		class(CoCoPotential), intent(in) :: this
		real(8), intent(in) :: R
		real(8) :: output
		
		real(8) :: damp, morsev, pauli, disp
		
		select case( this.model )
			case( MORSE )
				output = De*( exp(-2.0_8*alpha*(R*angs-Re)) - 2.0_8*exp(-alpha*(R*angs-Re)) )
			case( GUPTA )
				output = beta( xi, p, r0, R )
		end select
		
		output = output
	end function V
	
	!>
	!! @brief Returns the derivative of the potential at R in a.u.
	!!
	function dV( this, R ) result( output )
		class(CoCoPotential), intent(in) :: this
		real(8), intent(in) :: R
		real(8) :: output
		
		select case( this.model )
			case( MORSE )
				output = 2.0_8*De*alpha*( exp(alpha*(Re-R*angs))-exp(2*alpha*(Re-R*angs)) )
			case( GUPTA )
				output = dbeta( xi, p, r0, R )
		end select
		
 		output = output
	end function dV
	
	!>
	!! @brief Returns the numerical derivative of the potential at R in a.u.
	!!
	function NdV( this, R, nPoints, stepSize ) result( output )
		class(CoCoPotential), intent(in) :: this
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
	
	!>
	!! @brief Returns the nonlocal component of the potential for the particle i in a.u.
	!!
	function Vnl( this, i, positions, neighbourList ) result( output )
		class(CoCoPotential), intent(in) :: this
		integer :: i
		real(8), allocatable :: positions(:,:)
		type(IntegerList), optional, allocatable :: neighbourList(:)
		real(8) :: output
		
		class(IntegerListIterator), pointer :: iter
		integer :: j, N
		real(8) ::  rij(3)
		real(8) ::  d
		
		if ( this.model /= GUPTA ) return
		
		N = size(positions,dim=2)
		
		output = 0.0_8
		
		if( present(neighbourList) ) then
			if( i /= N ) then
				iter => neighbourList(i).begin
				do while( associated(iter) )
					j = iter.data
					
					if( j>i ) then
						rij(:) = positions(:,i)-positions(:,j)
						d = norm2(rij)
						
						output = output + beta( zeta**2, 2.0_8*q, r0, d )
					end if
					
					iter => iter.next
				end do
			end if
		else
			if( i /= N ) then
				do j=i+1,N
					
					if( i/=j ) then
					rij(:) = positions(:,i)-positions(:,j)
					d = norm2(rij)
					
					output = output + beta( zeta**2, 2.0_8*q, r0, d )
					end if
				end do
			end if
		end if
		
		output = -sqrt( output )
		
	end function Vnl
	
	!>
	!! @brief Returns the derivative of the nonlocal component of the potential for the particle i in a.u.
	!!
	function dVnl( this, i, positions, neighbourList, Vnl ) result( output )
		class(CoCoPotential), intent(in) :: this
		integer :: i
		real(8), allocatable :: positions(:,:)
		type(IntegerList), optional, allocatable :: neighbourList(:)
		real(8), optional :: Vnl
		real(8) :: output
		
		class(IntegerListIterator), pointer :: iter
		integer :: j, N
		real(8) ::  rij(3)
		real(8) ::  d
		
		if ( this.model /= GUPTA ) return
		
		N = size(positions,dim=2)
		
		output = 0.0_8
		
		if( present(neighbourList) ) then
			if( i /= N ) then
				iter => neighbourList(i).begin
				do while( associated(iter) )
					j = iter.data
					
					if( j>i ) then
						rij(:) = positions(:,i)-positions(:,j)
						d = norm2(rij)
						
						output = output + dbeta( zeta**2, 2.0_8*q, r0, d )
					end if
					
					iter => iter.next
				end do
			end if
		else
			if( i /= N ) then
				do j=i+1,N
				
					rij(:) = positions(:,i)-positions(:,j)
					d = norm2(rij)
				
					output = output + dbeta( zeta**2, 2.0_8*q, r0, d )
				end do
			end if
		end if
		
		if( present(Vnl) ) then
			output = output/2.0_8/Vnl
		else
			output = output/2.0_8/this.Vnl( i, positions, neighbourList )
		end if
		
	end function dVnl
	
	!>
	!! @brief Auxiliar function beta
	!!
	function beta( A, p, r0, r ) result( output )
		real(8) :: A
		real(8) :: p
		real(8) :: r0
		real(8) :: r
		real(8) :: output
		
		output = A*exp(-p*(r/r0-1.0_8))
	end function beta
	
	!>
	!! @brief Auxiliar function derivative of beta
	!!
	function dbeta( A, p, r0, r ) result( output )
		real(8) :: A
		real(8) :: p
		real(8) :: r0
		real(8) :: r
		real(8) :: output
		
		output = -p*A*exp(-p*(r/r0-1.0_8))/r0
	end function dbeta
	
	!>
	!! @brief Test
	!!
	subroutine CoCoPotential_test()
		real(8) :: r
		type(CoCoPotential) :: potential
		real(8), allocatable :: positions(:,:)
		
		call potential.init( GUPTA )
		
		allocate( positions(3,2) )
		positions(:,1) = 0.0_8
		positions(:,2) = 0.0_8
		
		do r = 1.0,15.0,0.01
			positions(3,2) = r
			write(*,"(4E15.7)") r, &
				potential.V( r )+potential.Vnl( 1, positions ), &
				potential.dV( r )+potential.dVnl( 1, positions ), &
				potential.NdV( r )
		end do
		
		deallocate( positions )
		
	end subroutine CoCoPotential_test

end module CoCoPotential_
