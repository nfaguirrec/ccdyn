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
! 					zeta = 2*1.4880_8*eV
! 					zeta = sqrt(2.0_8)*1.4880_8*eV
					zeta = 1.4880_8*eV
					q = 2.286_8
					xi = 0.0950_8*eV
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
				output = 2.0_8*beta( xi, p, r0, R )
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
				output = 2.0_8*dbeta( xi, p, r0, R )
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
			iter => neighbourList(i).begin
			do while( associated(iter) )
				j = iter.data
				
				if( i/=j ) then
					rij(:) = positions(:,i)-positions(:,j)
					d = norm2(rij)
					
					output = output + beta( 1.0_8, 2.0_8*q, r0, d )
				end if
				
				iter => iter.next
			end do
		else
			do j=1,N
				if( i/=j ) then
					
					rij(:) = positions(:,i)-positions(:,j)
					d = norm2(rij)
					
					output = output + beta( 1.0_8, 2.0_8*q, r0, d )
				end if
			end do
		end if
		
		output = -zeta*sqrt( output )
		
	end function Vnl
	
	!>
	!! @brief Returns the derivative of the nonlocal component of the potential for the particle i in a.u.
	!!
	function dVnl( this, i, j, positions, neighbourList, Vnl ) result( output )
		class(CoCoPotential), intent(in) :: this
		integer :: i, k
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
		
! 		if( present(neighbourList) ) then
! 			if( i /= N ) then
! 				iter => neighbourList(i).begin
! 				do while( associated(iter) )
! 					j = iter.data
! 					
! 					if( j>i ) then
! 						rij(:) = positions(:,i)-positions(:,j)
! 						d = norm2(rij)
! 						
! 						output = output + dbeta( zeta**2, 2.0_8*q, r0, d )
! 					end if
! 					
! 					iter => iter.next
! 				end do
! 			end if
! 		else
! 			do j=1,N
! 				if( i/=j ) then
					rij(:) = positions(:,i)-positions(:,j)
					d = norm2(rij)
				
					output = zeta**2*dbeta( 1.0_8, 2.0_8*q, r0, d )
! 				end if
! 				end if
! 			end do
! 		end if
		
		if( present(Vnl) ) then
			output = output/2.0_8/Vnl
		else
			output = output*( 1.0_8/this.Vnl( i, positions ) + 1.0_8/this.Vnl( j, positions ) )/2.0_8
! 			output = output/2.0_8/this.Vnl( i, positions, neighbourList )
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
	!! @brief
	!!
	subroutine getNeighbourList( positions, radius, neighbourList )
		real(8), allocatable, intent(in) :: positions(:,:)
		real(8), intent(in) :: radius
		type(IntegerList), allocatable, intent(inout) :: neighbourList(:)
		
		integer :: i, j, N
		
		N = size(positions,dim=2)
		
		if( allocated( neighbourList ) ) deallocate( neighbourList )
		allocate( neighbourList(N) )
		
		do i=1,N
			call neighbourList(i).clear()
			
			do j=1,N
				if( i /= j .and. norm2( positions(:,i) - positions(:,j) ) < radius ) then
					call neighbourList(i).append( j )
				end if
			end do
		end do
	end subroutine getNeighbourList
	
	!>
	!! @brief
	!!
	subroutine evaluatePotential( potential, positions, V, F )
		type(CoCoPotential), intent(in) :: potential
		real(8), allocatable, intent(in) :: positions(:,:)
		real(8), intent(out) :: V
		real(8), allocatable :: F(:,:)
		
		type(IntegerList), allocatable :: neighbourList(:)
		integer :: i, j, k, N
		real(8) :: d, rij(3)
		
		N = size(positions,dim=2)
		
! 		call getNeighbourList( positions, 10.0_8*angs, neighbourList )
		
		V = 0.0_8
		
		do i=1,N
			
			F(:,i) = 0.0_8
			do j=1,N
				if( i /= j ) then
					rij(:) = positions(:,i)-positions(:,j)
					d = norm2(rij)
					
					F(:,i) = F(:,i) - (rij(:)/d)*potential.dV( d )
					
! 					do k=1,N
						F(:,i) = F(:,i) - (rij(:)/d)*potential.dVnl( i, j, positions )
! 					end do
				end if
			end do
			
			if( i /= N ) then
				do j=i+1,N
				
					rij(:) = positions(:,i)-positions(:,j)
					d = norm2(rij)
				
					V = V + potential.V( d )
					
				end do
				V = V + potential.Vnl( i, positions )
! 				V = V + potential.Vnl( i, positions, neighbourList )
				
			end if
		end do
		V = V + potential.Vnl( N, positions )
! 		V = V + potential.Vnl( N, positions, neighbourList )
		
	end subroutine evaluatePotential
	
	!>
	!! @brief Test
	!!
	subroutine CoCoPotential_test()
		real(8) :: r, rij(3), d, r0
		type(CoCoPotential) :: potential
		real(8), allocatable :: positions(:,:)
		integer :: i, j, N
		real(8) :: V
		real(8), allocatable :: forces(:,:)
		
		r0 = 2.497_8*angs
		
		call potential.init( GUPTA )
		
! 		write(*,*) "#============="
! 		write(*,*) "# Testing N=2"
! 		write(*,*) "#============="
! 		allocate( positions(3,2) )
! 		positions(:,1) = 0.0_8
! 		positions(:,2) = 0.0_8
! 		
! 		do r = 1.0,15.0,0.01
! 			positions(3,2) = r
! 			write(*,"(4E15.7)") r, &
! 				potential.V( r )+potential.Vnl( 1, positions ), &
! 				potential.dV( r )+potential.dVnl( 1, positions ), &
! 				potential.NdV( r )
! 		end do
! 		
! 		deallocate( positions )
		
		write(*,*) ""
		write(*,*) ""
		write(*,*) "#============="
		write(*,*) "# Testing N=2"
		write(*,*) "#============="
		write(*,*) ""
		N = 2
		allocate( positions(3,N) )
		allocate( forces(3,N) )
		
		positions(:,1) = [ -0.337244_8,  0.114706_8,  0.257493_8 ]
		positions(:,2) = [  0.337244_8, -0.114706_8, -0.257493_8 ]
		
		positions = positions*r0
		
		call evaluatePotential( potential, positions, V, forces )
		
		write(*,"(A,F10.5,A)") "Energy   = ", V/eV, " eV"
		write(*,"(A,F10.5,A)") "Expected = ", -3.15065179, " eV"
		write(*,"(A,F10.5,A)") "Error    = ", V/eV-(-3.15065179_8), " eV"
		write(*,*) ""
		write(*,"(A,F10.5,A)") "Gradient = ", sqrt(sum(forces**2)/real(N,8))/(eV/r0), " eV/r0"
		
		deallocate( positions )
		deallocate( forces )
		
		write(*,*) ""
		write(*,*) "#============="
		write(*,*) "# Testing N=3"
		write(*,*) "#============="
		write(*,*) ""
		N = 3
		allocate( positions(3,N) )
		allocate( forces(3,N) )
		
		positions(:,1) = [  0.506198_8, -0.065114_8,  0.139191_8 ]
		positions(:,2) = [ -0.132223_8,  0.374877_8, -0.349048_8 ]
		positions(:,3) = [ -0.373975_8, -0.309763_8,  0.209857_8 ]
		
		positions = positions*r0

		call evaluatePotential( potential, positions, V, forces )
		
		write(*,"(A,F10.5,A)") "Energy   = ", V/eV, " eV"
		write(*,"(A,F10.5,A)") "Expected = ", -6.13875890, " eV"
		write(*,"(A,F10.5,A)") "Error    = ", V/eV-(-6.13875890_8), " eV"
		write(*,*) ""
		write(*,"(A,F10.5,A)") "Gradient = ", sqrt(sum(forces**2)/real(N,8))/(eV/r0), " eV/r0"
		
		deallocate( positions )
		deallocate( forces )
		
		write(*,*) ""
		write(*,*) "#============="
		write(*,*) "# Testing N=4"
		write(*,*) "#============="
		write(*,*) ""
		N = 4
		allocate( positions(3,N) )
		allocate( forces(3,N) )
		
		positions(:,1) = [ -0.329316,  0.227551,  0.411986 ]
		positions(:,2) = [  0.276943,  0.419839, -0.277485 ]
		positions(:,3) = [ -0.329784, -0.289447, -0.370707 ]
		positions(:,4) = [  0.382157, -0.357943,  0.236207 ]
		
		positions = positions*r0

		call evaluatePotential( potential, positions, V, forces )
		
		write(*,"(A,F10.5,A)") "Energy   = ", V/eV, " eV"
		write(*,"(A,F10.5,A)") "Expected = ", -9.53815918, " eV"
		write(*,"(A,F10.5,A)") "Error    = ", V/eV-(-9.53815918_8), " eV"
		write(*,*) ""
		write(*,"(A,F10.5,A)") "Gradient = ", sqrt(sum(forces**2)/real(N,8))/(eV/r0), " eV/r0"
		
		deallocate( positions )
		deallocate( forces )
		
		write(*,*) ""
		write(*,*) "#============="
		write(*,*) "# Testing N=5"
		write(*,*) "#============="
		write(*,*) ""
		N = 5
		allocate( positions(3,N) )
		allocate( forces(3,N) )
		
		positions(:,1) = [ -0.459800,  0.273487, -0.140265 ]
		positions(:,2) = [ -0.412354, -0.468448,  0.438350 ]
		positions(:,3) = [  0.291420,  0.157810,  0.442782 ]
		positions(:,4) = [  0.412353,  0.468449, -0.438353 ]
		positions(:,5) = [  0.168382, -0.431298, -0.302514 ]
		
		positions = positions*r0

		call evaluatePotential( potential, positions, V, forces )
		
		write(*,"(A,F10.5,A)") "Energy   = ", V/eV, " eV"
		write(*,"(A,F10.5,A)") "Expected = ", -12.80265501, " eV"
		write(*,"(A,F10.5,A)") "Error    = ", V/eV-(-12.80265501_8), " eV"
		write(*,*) ""
		write(*,"(A,F10.5,A)") "Gradient = ", sqrt(sum(forces**2)/real(N,8))/(eV/r0), " eV/r0"
		
		deallocate( positions )
		deallocate( forces )
		
		write(*,*) ""
		write(*,*) "#============="
		write(*,*) "# Testing N=8"
		write(*,*) "#============="
		write(*,*) ""
		N = 8
		allocate( positions(3,N) )
		allocate( forces(3,N) )
		
		positions(:,1) = [ -0.444567_8, -0.574217_8,  0.512405_8 ]
		positions(:,2) = [  0.445263_8, -0.417923_8,  0.207859_8 ]
		positions(:,3) = [ -0.545912_8,  0.305109_8,  0.158166_8 ]
		positions(:,4) = [  0.575406_8, -0.066947_8, -0.674060_8 ]
		positions(:,5) = [  0.063241_8,  0.146782_8,  0.874288_8 ]
		positions(:,6) = [ -0.276722_8, -0.407944_8, -0.416074_8 ]
		positions(:,7) = [ -0.194085_8,  0.494373_8, -0.712637_8 ]
		positions(:,8) = [  0.377376_8,  0.520768_8,  0.050054_8 ]
		
		positions = positions*r0
		
		call evaluatePotential( potential, positions, V, forces )
		
		write(*,"(A,F10.5,A)") "Energy   = ", V/eV, " eV"
		write(*,"(A,F10.5,A)") "Expected = ", -23.02759620, " eV"
		write(*,"(A,F10.5,A)") "Error    = ", V/eV-(-23.02759620_8), " eV"
		write(*,*) ""
		write(*,"(A,F10.5,A)") "Gradient = ", sqrt(sum(forces**2)/real(N,8))/(eV/r0), " eV/r0"
		
		deallocate( positions )
		deallocate( forces )
		
		write(*,*) ""
		write(*,*) "#=============="
		write(*,*) "# Testing N=30"
		write(*,*) "#=============="
		write(*,*) ""
		N = 30
		allocate( positions(3,N) )
		allocate( forces(3,N) )
		
		positions(:, 1) = [ -0.988126,  0.231577, -0.799667 ]
		positions(:, 2) = [ -0.152336, -0.427255,  1.060260 ]
		positions(:, 3) = [  0.571757,  1.035630,  0.796560 ]
		positions(:, 4) = [  0.310957, -0.454614,  0.235485 ]
		positions(:, 5) = [  1.184030,  0.036395, -0.710340 ]
		positions(:, 6) = [ -0.252440, -1.093760, -0.240743 ]
		positions(:, 7) = [  0.691723,  0.047620,  1.029380 ]
		positions(:, 8) = [ -1.125760, -0.746885, -0.523744 ]
		positions(:, 9) = [  0.588852, -0.754443, -0.651890 ]
		positions(:,10) = [ -0.076144,  0.494513,  1.356450 ]
		positions(:,11) = [  0.959675,  0.306719,  0.179666 ]
		positions(:,12) = [ -0.234201,  0.386507, -1.397530 ]
		positions(:,13) = [ -0.897472, -0.799054,  0.475069 ]
		positions(:,14) = [ -0.864839,  0.143141,  0.804319 ]
		positions(:,15) = [  0.543083, -0.169066, -1.413170 ]
		positions(:,16) = [ -0.658696,  1.081170,  0.921053 ]
		positions(:,17) = [  0.602339,  1.061560, -0.224682 ]
		positions(:,18) = [  0.627196,  0.765405, -1.140950 ]
		positions(:,19) = [  0.257293,  0.154844, -0.518307 ]
		positions(:,20) = [ -0.272115,  0.932492, -0.621631 ]
		positions(:,21) = [ -1.463250, -0.058290,  0.057038 ]
		positions(:,22) = [ -0.306105, -0.483913, -0.995046 ]
		positions(:,23) = [ -0.022393,  0.468934,  0.380754 ]
		positions(:,24) = [ -0.524156, -0.153576, -0.086838 ]
		positions(:,25) = [ -0.918658,  0.754363,  0.048313 ]
		positions(:,26) = [ -0.086085, -1.296690,  0.684414 ]
		positions(:,27) = [  1.240230, -0.602643,  0.079997 ]
		positions(:,28) = [  0.774346, -0.906656,  0.927025 ]
		positions(:,29) = [ -0.157927,  1.421760,  0.209492 ]
		positions(:,30) = [  0.649225, -1.375800,  0.079270 ]
		
		positions = positions*r0
		
		call evaluatePotential( potential, positions, V, forces )
		
		write(*,"(A,F10.5,A)") "Energy   = ", V/eV, " eV"
		write(*,"(A,F10.5,A)") "Expected = ", -104.1809983, " eV"
		write(*,"(A,F10.5,A)") "Error    = ", V/eV-(-104.1809983_8), " eV"
		write(*,*) ""
		write(*,"(A,F10.5,A)") "Gradient = ", sqrt(sum(forces**2)/real(N,8))/(eV/r0), " eV/r0"
		
		deallocate( positions )
		deallocate( forces )
		
		write(*,*) ""
		write(*,*) "#=============="
		write(*,*) "# Testing N=56"
		write(*,*) "#=============="
		write(*,*) ""
		N = 56
		allocate( positions(3,N) )
		allocate( forces(3,N) )
		
		positions(:, 1) = [ -0.903625, -0.101556, -0.024650 ]
		positions(:, 2) = [ -1.279710, -0.719799,  0.616673 ]
		positions(:, 3) = [ -0.434606,  0.396516,  0.693802 ]
		positions(:, 4) = [  0.385426,  0.595684, -0.664566 ]
		positions(:, 5) = [  0.217956, -1.613580, -0.152113 ]
		positions(:, 6) = [ -0.005269,  0.390654,  1.566230 ]
		positions(:, 7) = [  0.457834, -0.011575,  0.814862 ]
		positions(:, 8) = [  1.058680, -1.532780,  0.371505 ]
		positions(:, 9) = [ -0.480983,  0.770649, -0.220269 ]
		positions(:,10) = [  0.965259,  0.112565, -0.025248 ]
		positions(:,11) = [  1.322180,  0.952945,  0.301643 ]
		positions(:,12) = [  1.910340,  0.220979, -0.027352 ]
		positions(:,13) = [ -0.828198,  0.039428, -1.709050 ]
		positions(:,14) = [ -0.322058, -0.586123,  0.616602 ]
		positions(:,15) = [ -1.730920, -0.303799,  1.340380 ]
		positions(:,16) = [ -0.128357,  1.389090, -0.879219 ]
		positions(:,17) = [  0.171607, -1.229800, -1.084960 ]
		positions(:,18) = [ -0.964873,  1.190150,  0.508215 ]
		positions(:,19) = [ -0.663528, -1.440350,  0.301547 ]
		positions(:,20) = [ -0.788369, -0.197859,  1.369470 ]
		positions(:,21) = [  0.990682, -0.793004,  1.031510 ]
		positions(:,22) = [ -0.111230,  1.240740,  1.038480 ]
		positions(:,23) = [ -0.299316, -0.815608, -0.344876 ]
		positions(:,24) = [  0.542311, -0.759713,  0.172659 ]
		positions(:,25) = [ -0.630948, -1.645260, -0.671076 ]
		positions(:,26) = [ -0.926011,  0.804005, -1.081850 ]
		positions(:,27) = [  0.870075,  0.210519, -1.411030 ]
		positions(:,28) = [ -0.396294,  0.022363, -0.862294 ]
		positions(:,29) = [  0.970240, -0.777405, -1.467370 ]
		positions(:,30) = [  0.030797,  0.005413, -0.023852 ]
		positions(:,31) = [ -1.260820, -0.942664, -0.357066 ]
		positions(:,32) = [  0.805167,  0.826620,  1.158750 ]
		positions(:,33) = [ -0.695545, -1.173590,  1.267680 ]
		positions(:,34) = [  0.693501,  1.655450,  0.621988 ]
		positions(:,35) = [  0.730557,  1.445890, -0.350650 ]
		positions(:,36) = [  1.030510, -1.172880, -0.556398 ]
		positions(:,37) = [ -1.361560, -0.087050, -0.876801 ]
		positions(:,38) = [ -0.996322,  1.544620, -0.420440 ]
		positions(:,39) = [ -1.446530,  0.678626, -0.229589 ]
		positions(:,40) = [  0.497630, -0.384045, -0.741536 ]
		positions(:,41) = [  0.188800, -1.378840,  0.832689 ]
		positions(:,42) = [  1.461940, -0.282646, -0.758430 ]
		positions(:,43) = [ -0.042221,  0.625164, -1.534480 ]
		positions(:,44) = [  0.072276, -0.374439, -1.613020 ]
		positions(:,45) = [  1.507620, -0.666202,  0.174417 ]
		positions(:,46) = [  0.110290, -0.618478,  1.486950 ]
		positions(:,47) = [ -0.740567, -0.814883, -1.209060 ]
		positions(:,48) = [  0.360660,  0.826244,  0.297255 ]
		positions(:,49) = [  1.421620,  0.097628,  0.829766 ]
		positions(:,50) = [ -1.850230, -0.211508, -0.005516 ]
		positions(:,51) = [ -0.918918,  0.777533,  1.420980 ]
		positions(:,52) = [  1.347410,  0.717440, -0.679872 ]
		positions(:,53) = [ -1.397640,  0.309736,  0.697718 ]
		positions(:,54) = [ -0.153001,  1.624950,  0.102330 ]
		positions(:,55) = [  0.891634, -0.028392,  1.661160 ]
		positions(:,56) = [  0.744650,  1.192250, -1.312630 ]
		
		positions = positions*r0
		
		call evaluatePotential( potential, positions, V, forces )
		
		write(*,"(A,F10.5,A)") "Energy   = ", V/eV, " eV"
		write(*,"(A,F10.5,A)") "Expected = ", -206.0150172, " eV"
		write(*,"(A,F10.5,A)") "Error    = ", V/eV-(-206.0150172_8), " eV"
		write(*,*) ""
		write(*,"(A,F10.5,A)") "Gradient = ", sqrt(sum(forces**2)/real(N,8))/(eV/r0), " eV/r0"
		
		deallocate( positions )
		deallocate( forces )
		
	end subroutine CoCoPotential_test

end module CoCoPotential_
