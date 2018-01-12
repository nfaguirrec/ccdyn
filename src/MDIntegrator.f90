!!**********************************************************************************
!!    Consejo Superior de Investigaciones Científicas                              !
!!    Departamento de Física Atómica, Molecular y de Agregados                     !
!!    http://www.iff.csic.es/fama/                                                 !
!!                                                                                 !
!!    Authors:                                                                     !
!!    (2011-2012) Néstor F. Aguirre                                                !
!!                nfaguirrec@gmail.com                                             !
!!                nfaguirrec@iff.csic.es                                           !
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

!*
! @brief
!*
module MDIntegrator_
	implicit none
	private
	
	integer, public, parameter :: VELOCITY_VERLET = 1
	integer, public, parameter :: BEEMAN = 2
	integer, public, parameter :: BEEMAN_PREDCORR = 3
	integer, public, parameter :: BEEMAN_PREDCORR_AM = 4
	
	public :: &
		MDIntegrator_charIDtoID
	
	type, public :: MDIntegrator
		real(8) :: timeStep
		integer :: maxIter
		integer :: ttype
		
		integer, pointer :: nParticles
		real(8), pointer :: positions(:,:)
		real(8), pointer :: velocities(:,:)
		real(8), pointer :: accelerations(:,:)
		real(8), pointer :: forces(:,:)
		
		integer :: step
		
		contains
			procedure :: init
			procedure :: destroy
			
			procedure :: iterate
	end type MDIntegrator
	
	real(8), allocatable :: posBufferT0(:,:)
	real(8), allocatable :: accelBufferT0(:,:)
	real(8), allocatable :: accelBufferTm1(:,:)
	
	interface
		subroutine prototypeComputeForces()
		end subroutine prototypeComputeForces
	end interface
	
	contains
	
	!*
	! @brief Default constructor
	!*
	subroutine init( this, nParticles, positions, velocities, accelerations, forces )
		class(MDIntegrator) :: this
		integer, target :: nParticles
		real(8), target :: positions(:,:)
		real(8), target :: velocities(:,:)
		real(8), target :: accelerations(:,:)
		real(8), target :: forces(:,:)
		
		this.nParticles => nParticles
		this.positions => positions
		this.velocities => velocities
		this.accelerations => accelerations
		this.forces => forces
		
		this.timeStep = 1e-3
		this.maxIter = 1000
		this.step = 0
		this.ttype = VELOCITY_VERLET
		
		if( allocated(posBufferT0) ) deallocate(posBufferT0)
		if( allocated(accelBufferT0) ) deallocate(accelBufferT0)
		if( allocated(accelBufferTm1) ) deallocate(accelBufferTm1)
		
		allocate( posBufferT0(3,nParticles) )
		allocate( accelBufferT0(3,nParticles) )
		allocate( accelBufferTm1(3,nParticles) )
		
		posBufferT0 = 0.0_8
		accelBufferT0 = 0.0_8
		accelBufferTm1 = 0.0_8
	end subroutine init
	
	!*
	! @brief Destructor
	!*
	subroutine destroy( this )
		class(MDIntegrator) :: this
		
		nullify( this.nParticles )
		nullify( this.positions )
		nullify( this.velocities )
		nullify( this.accelerations )
		nullify( this.forces )
		
		this.timeStep = 0.0_8
		this.maxIter = 0
		this.step = 0
		
		deallocate( posBufferT0 )
		deallocate( accelBufferT0 )
		deallocate( accelBufferTm1 )
	end subroutine destroy
	
	!*
	! @brief
	!*
	subroutine control( this, maxIter, timeStep )
		class(MDIntegrator) :: this
		real(8), intent(in) :: maxIter
		real(8), intent(in) :: timeStep
		
	end subroutine control
	
	!*
	! @brief
	!*
	subroutine iterate( this, computeForces )
		class(MDIntegrator) :: this
		procedure(prototypeComputeForces) :: computeForces
		
#define x0 posBufferT0
#define x this.positions
#define v this.velocities
#define a this.accelerations
#define a0 accelBufferT0
#define am1 accelBufferTm1
#define dt this.timeStep
		
		select case( this.ttype )
			case( VELOCITY_VERLET )
				x = x + v*dt + 0.5_8*a*dt*dt
				a0 = a
				call computeForces() ! a <--
				v = v + 0.5_8*(a0+a)*dt
				
			case( BEEMAN )
				x = x + v*dt + (1.0_8/6.0_8)*( 4.0_8*a0 - am1 )*dt*dt
				a0 = a
				call computeForces() ! a <--
				v = v + (1.0_8/6.0_8)*( 2.0_8*a + 5.0_8*a0 - am1 )*dt
				am1 = a0
				
			case( BEEMAN_PREDCORR )
				x0 = x
				x = x + v*dt + (1.0_8/6.0_8)*( 4.0_8*a0 - am1 )*dt*dt
				a0 = a
				call computeForces() ! a <--
				x = x0 + v*dt + (1.0_8/6.0_8)*( a + 2.0_8*a0 )*dt*dt
				v = ( x - x0 )/dt + (1.0_8/6.0_8)*( 2.0_8*a + a0 )*dt
				am1 = a0
				
			case( BEEMAN_PREDCORR_AM )
				x0 = x
				x = x + v*dt + (1.0_8/6.0_8)*( 4.0_8*a0 - am1 )*dt*dt
				a0 = a
				call computeForces() ! a <--
				x = x0 + v*dt + (1.0_8/6.0_8)*( a + 2.0_8*a0 )*dt*dt
				v = v + (1.0_8/12.0_8)*( 5.0_8*a + 8.0_8*a0 - am1 )*dt
				am1 = a0
				
		end select

#undef x0
#undef x
#undef v
#undef a
#undef a0
#undef am1
#undef dt
		
		this.step = this.step + 1
	end subroutine iterate
	
	!*
	! @brief Converts the character ID to integer ID
	!*
	function MDIntegrator_charIDtoID( charID ) result( id )
		character(*), intent(in) :: charID
		integer :: id
		
		select case( charID )
			case( "VelocityVerlet" )
				id = VELOCITY_VERLET
			case( "VV" )
				id = VELOCITY_VERLET
			case( "Beeman" )
				id = BEEMAN
			case( "B0" )
				id = BEEMAN
			case( "BeemanPredCorr" )
				id = BEEMAN_PREDCORR
			case( "B1" )
				id = BEEMAN_PREDCORR
			case( "BeemanPredCorrAM" )
				id = BEEMAN_PREDCORR_AM
			case( "B2" )
				id = BEEMAN_PREDCORR_AM
		end select
	end function MDIntegrator_charIDtoID
	
end module MDIntegrator_
