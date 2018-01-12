!!**********************************************************************************
!!    Consejo Superior de Investigaciones Científicas                              !
!!    Departamento de Física Atómica, Molecular y de Agregados                     !
!!    http://www.iff.csic.es/fama/                                                 !
!!                                                                                 !
!!    Authors:                                                                     !
!!    (2011-2015) Néstor F. Aguirre                                                !
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

!*
! @brief
!*
module ClassicalDynamics_
	use UnitsConverter_
	use IntegerList_
	
	use HeHePotential_
	use HeAuPotential_
	use AuAuPotential_
	use MgMgPotential_
	use HeMgPotential_
	use ArArPotential_
	use RandomSampler_
	use HeDroplet_
	use MDIntegrator_
	use AtomicElementsDB_
	use EmbeddingInterface_
	implicit none
	private
	
	!**
	! @brief Public parameters
	!**
	integer, public, parameter :: XYZ = 1
	
	public :: &
		ClassicalDynamics_test
		
	type, public :: ClassicalDynamics
		integer :: nAtoms
		character(2), allocatable :: symbols(:)
		real(8), allocatable :: positions(:,:)
		real(8), allocatable :: velocities(:,:)
		real(8), allocatable :: accelerations(:,:)
		real(8), allocatable :: forces(:,:)
		
		logical :: useNeighbourList
		real(8) :: neighbourRadius
		type(IntegerList), allocatable :: neighbourList(:)
		
		real(8) :: timeStep
		real(8) :: potentialEnergy
		real(8) :: kineticEnergy
		real(8) :: totalEnergy
		real(8) :: averEnergy
		
		integer :: historyStep
		real(8) :: historyTimeStep
		logical :: externalPotential
		
		type(EmbeddingInterface), pointer :: embedding
		
		contains
			generic :: init => fromHeDroplet, fromFile
			procedure :: fromHeDroplet
			procedure :: fromFile
			procedure :: destroy
			
			procedure :: centerOfMass
			procedure :: updateNeighbourList
			procedure :: temperature
			
			procedure :: setExternalField
			
			procedure :: center
			procedure :: push
			procedure :: move
			procedure :: thermalize
			procedure :: freeEvolution
			
			procedure :: saveGeometry
			procedure :: saveVelocities
! 			procedure :: savePDF
! 			procedure :: saveLDF
! 			procedure :: saveADF
	end type ClassicalDynamics
	
	!! This is neccesary only for RandomSampler in ClassicalDynamics.init
	type(HeDroplet), pointer, private :: dropletPtr
	
	!! This is neccesary only for calculateForces
	type(ClassicalDynamics), pointer, private :: cDropletPtr
	procedure(protVext), pointer, private :: Vext
	procedure(protdVext), pointer, private :: dVext
	
	!! This is neccesary for potentials
	type(HeHePotential), private :: HeHePot
	type(HeAuPotential), private :: HeAuPot
	type(MgMgPotential), private :: MgMgPot
	type(AuAuPotential), private :: AuAuPot
	type(HeMgPotential), private :: HeMgPot
	type(ArArPotential), private :: ArArPot
	
	interface
		function protVext( symbol, x, y, z ) result( V )
			character(*), intent(in) :: symbol
			real(8), intent(in) :: x
			real(8), intent(in) :: y
			real(8), intent(in) :: z
			
			real(8) :: V
		end function protVext
		
		function protdVext( symbol, x, y, z, i ) result( dV )
			character(*), intent(in) :: symbol
			real(8), intent(in) :: x
			real(8), intent(in) :: y
			real(8), intent(in) :: z
			integer, intent(in) :: i
			
			real(8) :: dV
		end function protdVext
	end interface
	
	contains
	
	!>
	!! This is only neccesary for RandomSampler in ClassicalDynamics.init
	!!
	function dropletDist( r ) result( output )
		real(8), intent(in) :: r(:)
		real(8) :: output
		
		output = dropletPtr.density( r )
	end function dropletDist
	
	!>
	!! @brief Constructor
	!!
	subroutine fromHeDroplet( this, nAtoms, droplet, useNeighbourList, neighbourRadius, embedding )
		class(ClassicalDynamics), target :: this 
		integer, intent(in) :: nAtoms
		type(HeDroplet), target, intent(in) :: droplet
		logical, optional :: useNeighbourList
		real(8), optional :: neighbourRadius
		type(EmbeddingInterface), target, optional :: embedding
		
		type(RandomSampler) :: rs
		
		real(8) :: cm(3)
		integer :: i
		
		this.useNeighbourList = .true.
		if( present(useNeighbourList) ) this.useNeighbourList = useNeighbourList
		
		this.neighbourRadius = 10.0_8
		if( present(neighbourRadius) ) this.neighbourRadius = neighbourRadius
		
		this.nAtoms = nAtoms
		allocate( this.symbols(nAtoms) )
		allocate( this.positions(3,nAtoms) )
		allocate( this.velocities(3,nAtoms) )
		allocate( this.accelerations(3,nAtoms) )
		allocate( this.forces(3,nAtoms) )
		
		allocate( this.neighbourList(nAtoms) )
		
		this.symbols = "He"
		this.symbols(1) = "Au"
		
		write(6,'(A)', advance="no") "# Building positions ... "
		dropletPtr => droplet
		
		call rs.init( nDim=3 )
		call rs.setRange( 1, [-droplet.xmax,droplet.xmax] )
		call rs.setRange( 2, [-droplet.ymax,droplet.ymax] )
		call rs.setRange( 3, [-droplet.zmax,droplet.zmax] )
		call rs.buildSample( this.positions, dropletDist, delta=20.0_8, rMin=3.0_8 )
		
		this.positions(:,1) = 0.0_8
		this.positions = this.positions*angs
		
		dropletPtr => null()
		write(6,'(A)') "OK"
		
		call this.center()
		
		this.velocities = 0.0_8
		this.accelerations = 0.0_8
		this.forces = 0.0_8
		
		this.timeStep = 1.0e-3
		this.kineticEnergy = 0.0_8
		this.potentialEnergy = 0.0_8
		this.totalEnergy = 0.0_8
		this.averEnergy = 0.0_8
		
		this.historyStep = 0
		this.historyTimeStep = 0.0_8
		this.externalPotential = .false.
		
		if( present(embedding) ) then
			this.embedding => embedding
		end if
		
		cDropletPtr => this
		
		call HeHePot.init()
		call HeAuPot.init()
		call MgMgPot.init()
		call AuAuPot.init()
		call HeMgPot.init()
		call ArArPot.init()
	end subroutine fromHeDroplet
	
	!*
	! @brief Destructor
	!*
	subroutine fromFile( this, geomFileName, velFileName, useNeighbourList, neighbourRadius, embedding )
		class(ClassicalDynamics), target :: this
		character(*), intent(in) :: geomFileName
		character(*), optional, intent(in) :: velFileName
		logical, optional :: useNeighbourList
		real(8), optional :: neighbourRadius
		type(EmbeddingInterface), target, optional :: embedding
		
		real(8) :: cm(3)
		integer :: i
		integer :: vNAtoms
		character(5) :: vSymbol
		
		this.useNeighbourList = .true.
		if( present(useNeighbourList) ) this.useNeighbourList = useNeighbourList
		
		this.neighbourRadius = 10.0_8
		if( present(neighbourRadius) ) this.neighbourRadius = neighbourRadius
		
		!---------------------------------------------------------------
		write(6,'(A)', advance="no") "# Reading positions ... "
		
		open( 10, file=geomFileName, status="old" )
		read(10,'(I10)') this.nAtoms
		
		allocate( this.symbols(this.nAtoms) )
		allocate( this.positions(3,this.nAtoms) )
		allocate( this.velocities(3,this.nAtoms) )
		allocate( this.accelerations(3,this.nAtoms) )
		allocate( this.forces(3,this.nAtoms) )
		
		allocate( this.neighbourList(this.nAtoms) )
		
		read(10,'(X)')
		do i=1,this.nAtoms
			read(10,*) this.symbols(i), this.positions(1,i), this.positions(2,i), this.positions(3,i)
		end do
		close(10)
		
		this.positions = this.positions*angs
		
		write(6,'(A)') "OK"
		!---------------------------------------------------------------
		
		if( present(velFileName) ) then
			write(6,'(A)', advance="no") "# Reading velocities ... "
			
			open( 10, file=velFileName, status="old" )
			read(10,'(I10)') vNAtoms
			
			if( this.nAtoms /= vNAtoms ) then
				write(*,"(A)") "### ERROR ### Number of atoms in "//trim(velFileName)//"is inconsistent with the number of atoms in"//trim(geomFileName)
				stop
			end if
			
			read(10,'(X)')
			do i=1,this.nAtoms
				read(10,*) vSymbol, this.velocities(1,i), this.velocities(2,i), this.velocities(3,i)
				
				if( trim(vSymbol) /= trim(this.symbols(i)) ) then
					write(*,"(A)") "### ERROR ### Inconsistent identity of atoms"
					stop
				end if
			end do
			close(10)
			
			this.velocities = this.velocities*(angs/ps)
			
			write(6,'(A)') "OK"
		end if
		
		write(6,'(A)', advance="no") "# Fixing center of mass ... "
		
		cm = this.centerOfMass()
		do i=1,this.nAtoms
			this.positions(:,i) = this.positions(:,i) - cm(:)
		end do
		
		write(6,'(A,3F10.4,A)') "OK (", this.centerOfMass(), "  )"
		
		this.velocities = 0.0_8
		this.accelerations = 0.0_8
		this.forces = 0.0_8
		
		this.timeStep = 1.0e-3
		this.kineticEnergy = 0.0_8
		this.potentialEnergy = 0.0_8
		this.totalEnergy = 0.0_8
		this.averEnergy = 0.0_8
		
		this.historyStep = 0
		this.historyTimeStep = 0.0_8
		
		if( present(embedding) ) then
			this.embedding => embedding
		end if
		
		cDropletPtr => this
		
		call HeHePot.init()
		call HeAuPot.init()
		call MgMgPot.init()
		call AuAuPot.init()
		call HeMgPot.init()
		call ArArPot.init()
	end subroutine fromFile
	
	!*
	! @brief Destructor
	!*
	subroutine destroy( this )
		class(ClassicalDynamics) :: this
		
		integer :: i
		
		deallocate( this.symbols )
		deallocate( this.positions )
		deallocate( this.velocities )
		deallocate( this.accelerations )
		deallocate( this.forces )
		
		do i=1,this.nAtoms
			call this.neighbourList(i).clear()
		end do
		deallocate( this.neighbourList )
		
		dropletPtr => null()
		cDropletPtr => null()
	end subroutine destroy
	
	!*
	! @brief Returns the center of mass of the droplet
	!*
	function centerOfMass( this ) result( cm )
		class(ClassicalDynamics) :: this
		real(8) :: cm(3)
		
		integer :: i
		
		cm = 0.0_8
		do i=1,this.nAtoms
			cm(:) = cm(:) + this.positions(:,i)
		end do
		cm(:) = cm(:)/this.nAtoms
	end function centerOfMass

	!*
	! @brief
	!*
	subroutine updateNeighbourList( this )
		class(ClassicalDynamics) :: this
		
		integer :: i, j
		
#define pos this.positions	
		if( this.useNeighbourList ) then
			do i=1,this.nAtoms
				call this.neighbourList(i).clear()
				
				do j=1,this.nAtoms
					if( i /= j .and. norm2( pos(:,i) - pos(:,j) ) < this.neighbourRadius ) then
						call this.neighbourList(i).append( j )
					end if
				end do
			end do
		end if
#undef pos
	end subroutine updateNeighbourList
	
	!>
	!! @brief Returns the temperature in atomic units
	!! Taken from: ftp://ftp.rush.edu/users/molebio/Jay_Bardhan/Hunenberger05%20thermostat%20algorithms%20for%20molecular%20dynamics%20simulations.pdf
	!!
	function temperature( this ) result( T )
		class(ClassicalDynamics) :: this
		real(8) :: T
		
		integer :: i
		
#define sym cDropletPtr.symbols
		T = 0.0_8
		do i=1,this.nAtoms
			T = T + AtomicElementsDB_instance.atomicMass(sym(i))*sum(this.velocities(:,i)**2)
		end do
		T = T/real(3.0*this.nAtoms-3,8) ! 3 constrains, angular momentum conservation
#undef sym
	end function temperature
	
	!*
	! @brief
	!*
	subroutine setExternalField( this, V, dV )
		class(ClassicalDynamics) :: this
		procedure(protVext) :: V
		procedure(protdVext) :: dV
		
		Vext => V
		dVext => dV
		this.externalPotential = .true.
	end subroutine setExternalField
	
	!*
	! @brief
	!*
	function potBySym( symbol1, symbol2, r ) result( output )
		character(*), intent(in) :: symbol1
		character(*), intent(in) :: symbol2
		real(8), intent(in) :: r
		real(8) :: output
		
		output = 0.0_8
! 		if( trim(symbol1) == "He" .and. trim(symbol2) == "He" ) then
! 			
! 			output = HeHePot.V( r/angs )*cm1
! 		
! 		else if( trim(symbol1) == "Au" .and. trim(symbol2) == "Au" ) then
! 		
! 			output = AuAuPot.V( r/angs )*cm1
! 			
! 		else if( trim(symbol1) == "Mg" .and. trim(symbol2) == "Mg" ) then
! 		
! 			output = MgMgPot.V( r/angs )*cm1
! 			
! 		else if( ( trim(symbol1) == "He" .and. trim(symbol2) == "Au" ) .or. &
! 			 ( trim(symbol1) == "Au" .and. trim(symbol2) == "He" ) ) then
! 			
! 			output = HeAuPot.V( r/angs )*cm1
! 			
! 		else if( ( trim(symbol1) == "He" .and. trim(symbol2) == "Mg" ) .or. &
! 			 ( trim(symbol1) == "Mg" .and. trim(symbol2) == "He" ) ) then
! 			
! 			output = HeMgPot.V( r/angs )*cm1
! 			
! 		else if( trim(symbol1) == "Ar" .and. trim(symbol2) == "Ar" ) then
			
			output = ArArPot.V( r )
		
! 		else
! 			write(*,*) "### ERROR ### Potential pair (",symbol1,",",symbol2,") is not implemented"
! 			stop
! 		end if
		
	end function potBySym
	
	!*
	! @brief
	!*
	function dpotBySym( symbol1, symbol2, r ) result( output )
		character(*), intent(in) :: symbol1
		character(*), intent(in) :: symbol2
		real(8), intent(in) :: r
		real(8) :: output
		
		output = 0.0_8
! 		if( trim(symbol1) == "He" .and. trim(symbol2) == "He" ) then
! 			
! 			output = HeHePot.dV( r/angs )*(cm1/angs)
! 		
! 		else if( trim(symbol1) == "Au" .and. trim(symbol2) == "Au" ) then
! 		
! 			output = AuAuPot.dV( r/angs )*(cm1/angs)
! 			
! 		else if( trim(symbol1) == "Mg" .and. trim(symbol2) == "Mg" ) then
! 		
! 			output = MgMgPot.dV( r/angs )*(cm1/angs)
! 			
! 		else if( ( trim(symbol1) == "He" .and. trim(symbol2) == "Au" ) .or. &
! 			 ( trim(symbol1) == "Au" .and. trim(symbol2) == "He" ) ) then
! 			
! 			output = HeAuPot.dV( r/angs )*(cm1/angs)
! 			
! 		else if( ( trim(symbol1) == "He" .and. trim(symbol2) == "Mg" ) .or. &
! 			 ( trim(symbol1) == "Mg" .and. trim(symbol2) == "He" ) ) then
! 			
! 			output = HeMgPot.dV( r/angs )*(cm1/angs)
! 			
! 		else if( trim(symbol1) == "Ar" .and. trim(symbol2) == "Ar" ) then
! 			
			output = ArArPot.dV( r )
			
! 		else
! 			write(*,*) "### ERROR ### Potential pair (",symbol1,",",symbol2,") is not implemented"
! 			stop
! 		end if
		
	end function dpotBySym
	
	!>
	!! @brief
	!!
	subroutine computeForces()
		if( associated(cDropletPtr.embedding) ) then
			call computeForcesQuantum()
		else
			call computeForcesClassic()
		end if
	end subroutine computeForces
	
	!>
	!! @brief
	!!
	subroutine computeForcesClassic()
		integer i, j, k
		real(8) ::  rij(3)
		real(8) ::  d
		class(IntegerListIterator), pointer :: iter
		
#define N cDropletPtr.nAtoms
#define sym cDropletPtr.symbols
#define pos cDropletPtr.positions
#define vel cDropletPtr.velocities
#define acc cDropletPtr.accelerations
#define F cDropletPtr.forces
#define K cDropletPtr.kineticEnergy
#define U cDropletPtr.potentialEnergy
#define E cDropletPtr.totalEnergy
			
		K = 0.0_8
		U = 0.0_8
		E = 0.0_8
		
		do i=1,N
			F(:,i) = 0.0_8
			
			! Computes the forces on the i-particle induced
			! from the rest of the particles
			if( cDropletPtr.useNeighbourList ) then
				iter => cDropletPtr.neighbourList(i).begin
				do while( associated(iter) )
					j = iter.data
					
					if( i /= j ) then
						rij(:) = pos(:,i)-pos(:,j)
						d = norm2(rij)
						
						F(:,i) = F(:,i) - (rij(:)/d)*dpotBySym( sym(i), sym(j), d )
					end if
					
					iter => iter.next
				end do
			else
				do j=1,N
					if( i /= j ) then
						rij(:) = pos(:,i)-pos(:,j)
						d = norm2(rij)
						
						F(:,i) = F(:,i) - (rij(:)/d)*dpotBySym( sym(i), sym(j), d )
					end if
				end do
			end if
			
			! Computes the potential on the i-particle induced
			! from the rest of the particles
			if( cDropletPtr.useNeighbourList ) then
				if( i /= N ) then
					iter => cDropletPtr.neighbourList(i).begin
					do while( associated(iter) )
						j = iter.data
						
						if( j>i ) then    ! Si i>j y j es clasico
							rij(:) = pos(:,i)-pos(:,j)
							d = norm2(rij)
							
							U = U + potBySym( sym(i), sym(j), d )
						end if
						
						iter => iter.next
					end do
				end if
			else
				if( i /= N ) then
					do j=i+1,N
					
						rij(:) = pos(:,i)-pos(:,j)
						d = norm2(rij)
					
						U = U + potBySym( sym(i), sym(j), d )
						
! 		do d = 1.0*angs,10.0*angs,0.01*angs
! 			write(*,"(4E15.7)") d/angs, potBySym( "Ar", "Ar", d )/cm1, dpotBySym( "Ar", "Ar", d )/(cm1/angs)
! 		end do
! 		stop
					end do
				end if
			end if
				
			! Computes the forces on the i-particle induced
			! from the external potential
			if( cDropletPtr.externalPotential ) then
				do k=1,3
					F(k,i) = F(k,i) - dVext( sym(i), pos(1,i), pos(2,i), pos(3,i), k )
				end do
			end if
			
			! Computes the potential on the i-particle induced
			! from the external potential
			if( cDropletPtr.externalPotential ) then
				U = U + Vext( sym(i), pos(1,i), pos(2,i), pos(3,i) )
			end if
			
			acc(:,i) = F(:,i)/AtomicElementsDB_instance.atomicMass(sym(i))  ! las unidades de la masas ya se pusieron en F(k,i)
			K = K + 0.5_8*AtomicElementsDB_instance.atomicMass(sym(i))*sum(vel(:,i)**2)
			
			! L = A
			! V = A/ps
			! U = amu*A^2/ps^2
			! K = amu*A^2/ps^2
		end do
		
		E = K + U
		
#undef N
#undef pos
#undef vel
#undef acc
#undef F
#undef K
#undef U
#undef E
		
	end subroutine computeForcesClassic
	
	!>
	!! @brief
	!!
	subroutine computeForcesQuantum()
		integer i, j, k
		real(8) ::  rij(3)
		real(8) ::  d
		class(IntegerListIterator), pointer :: iter
		
#define N cDropletPtr.nAtoms
#define sym cDropletPtr.symbols
#define pos cDropletPtr.positions
#define vel cDropletPtr.velocities
#define acc cDropletPtr.accelerations
#define F cDropletPtr.forces
#define K cDropletPtr.kineticEnergy
#define U cDropletPtr.potentialEnergy
#define E cDropletPtr.totalEnergy
			
		K = 0.0_8
		U = 0.0_8
		E = 0.0_8
		
		call cDropletPtr.embedding.update( sym(cDropletPtr.embedding.targetAtoms(:)), pos(:,cDropletPtr.embedding.targetAtoms(:)), debug=.true. )
		
		do i=1,N
			F(:,i) = 0.0_8
			
 			if( any(i==cDropletPtr.embedding.targetAtoms) ) then  ! Si i es cuántico
					
				! Computes the forces on the i-particle induced
				! from the rest of the particles
				if( cDropletPtr.useNeighbourList ) then
					iter => cDropletPtr.neighbourList(i).begin
					do while( associated(iter) )
						j = iter.data
						if ( i /= j ) then
							if( all(j/=cDropletPtr.embedding.targetAtoms) ) then    ! Si j es clasico
								rij(:) = pos(:,i)-pos(:,j)
								d = norm2(rij)
								
								F(:,i) = F(:,i) - (rij(:)/d)*dpotBySym( sym(i), sym(j), d )
							end if
						end if
						
						iter => iter.next
					end do
				else
					do j=1,N
						if ( i /= j ) then
							if( all(j/=cDropletPtr.embedding.targetAtoms) ) then    ! Si j es clasico
								rij(:) = pos(:,i)-pos(:,j)
								d = norm2(rij)
								
								F(:,i) = F(:,i) - (rij(:)/d)*dpotBySym( sym(i), sym(j), d )
							end if
						end if
					end do
				end if
				
				F(:,i) = F(:,i) - cDropletPtr.embedding.gradient(:,i)
				
				! Computes the potential on the i-particle induced
				! from the rest of the particles
				if( cDropletPtr.useNeighbourList ) then
					if( i /= N ) then
						iter => cDropletPtr.neighbourList(i).begin
						do while( associated(iter) )
							j = iter.data
							
							if( j>i ) then
								if( all(j/=cDropletPtr.embedding.targetAtoms) ) then    ! Si j es clasico
									rij(:) = pos(:,i)-pos(:,j)
									d = norm2(rij)
									
									U = U + potBySym( sym(i), sym(j), d )
								end if
							end if
							
							iter => iter.next
						end do
					end if
				else
					if( i /= N ) then
						do j=i+1,N
							if( all(j/=cDropletPtr.embedding.targetAtoms) ) then    ! Si j es clasico
								rij(:) = pos(:,i)-pos(:,j)
								d = norm2(rij)
							
								U = U + potBySym( sym(i), sym(j), d )
							end if
						end do
					end if
				end if
				
				U = U + cDropletPtr.embedding.energy
				
 			else ! Si i no es cuántico
				
				! Computes the forces on the i-particle induced
				! from the rest of the particles
				if( cDropletPtr.useNeighbourList ) then
					iter => cDropletPtr.neighbourList(i).begin
					do while( associated(iter) )
						j = iter.data
						if( i /= j ) then
							if( all(j/=cDropletPtr.embedding.targetAtoms) ) then      ! Si i/=j y j es clasico
								rij(:) = pos(:,i)-pos(:,j)
								d = norm2(rij)
							
								F(:,i) = F(:,i) - (rij(:)/d)*dpotBySym( sym(i), sym(j), d )
							else ! Esto es redundante, pero es para decir que tambien hay interacción por pares para clasico-cuantico
								rij(:) = pos(:,i)-pos(:,j)
								d = norm2(rij)
								
								F(:,i) = F(:,i) - (rij(:)/d)*dpotBySym( sym(i), sym(j), d )
							end if
						end if
						
						iter => iter.next
					end do
				else
					do j=1,N
						if( associated(cDropletPtr.embedding) ) then
							if( i /= j ) then
								if( all(j/=cDropletPtr.embedding.targetAtoms) ) then     ! Si i/=j y j es clasico
									rij(:) = pos(:,i)-pos(:,j)
									d = norm2(rij)
									
									F(:,i) = F(:,i) - (rij(:)/d)*dpotBySym( sym(i), sym(j), d )
								else ! Esto es redundante, pero es para decir que tambien hay interacción por pares para clasico-cuantico
									rij(:) = pos(:,i)-pos(:,j)
									d = norm2(rij)
									
									F(:,i) = F(:,i) - (rij(:)/d)*dpotBySym( sym(i), sym(j), d )
								end if
							end if
						end if
					end do
				end if
				
				! Computes the potential on the i-particle induced
				! from the rest of the particles
				if( cDropletPtr.useNeighbourList ) then
					if( i /= N ) then
						iter => cDropletPtr.neighbourList(i).begin
						do while( associated(iter) )
							j = iter.data
							
							if( j>i ) then
								if( all(j/=cDropletPtr.embedding.targetAtoms) ) then    ! Si j es clasico
									rij(:) = pos(:,i)-pos(:,j)
									d = norm2(rij)
									
									U = U + potBySym( sym(i), sym(j), d )
								else ! Esto es redundante, pero es para decir que tambien hay interacción por pares para clasico-cuantico
									rij(:) = pos(:,i)-pos(:,j)
									d = norm2(rij)
									
									U = U + potBySym( sym(i), sym(j), d )
								end if
							end if
							
							iter => iter.next
						end do
					end if
				else
					if( i /= N ) then
						do j=i+1,N
						
							if( all(j/=cDropletPtr.embedding.targetAtoms) ) then    ! Si j es clasico
								rij(:) = pos(:,i)-pos(:,j)
								d = norm2(rij)
							
								U = U + potBySym( sym(i), sym(j), d )
							else ! Esto es redundante, pero es para decir que tambien hay interacción por pares para clasico-cuantico
								rij(:) = pos(:,i)-pos(:,j)
								d = norm2(rij)
							
								U = U + potBySym( sym(i), sym(j), d )
							end if
							
						end do
					end if
				end if
				
 			end if
 			
			! Computes the forces on the i-particle induced
			! from the external potential
			if( cDropletPtr.externalPotential ) then
				do k=1,3
					F(k,i) = F(k,i) - dVext( sym(i), pos(1,i), pos(2,i), pos(3,i), k )
				end do
			end if
			
			! Computes the potential on the i-particle induced
			! from the external potential
			if( cDropletPtr.externalPotential ) then
				U = U + Vext( sym(i), pos(1,i), pos(2,i), pos(3,i) )
			end if
			
			acc(:,i) = F(:,i)/AtomicElementsDB_instance.atomicMass(sym(i))  ! las unidades de la masas ya se pusieron en F(k,i)
			K = K + 0.5_8*AtomicElementsDB_instance.atomicMass(sym(i))*sum(vel(:,i)**2)
			
			! L = A
			! V = A/ps
			! U = amu*A^2/ps^2
			! K = amu*A^2/ps^2
		end do
		
		E = K + U
		
! 		write(*,*) "------------------------------"
! 		write(*,"(X,A,F10.5)") "Energy = ", E
! 		write(*,*) "Final gradient"
! 		do i=1,N
! 			if( trim(sym(i)) == "Au" ) then
! 				write(*,"(I5,A5,3F10.5)") i, trim(sym(i)), F(1,i), F(2,i), F(3,i)
! 			end if
! 		end do
! 		write(*,*) "------------------------------"
! 		stop
		
#undef N
#undef pos
#undef vel
#undef acc
#undef F
#undef K
#undef U
#undef E
		
	end subroutine computeForcesQuantum
	
	!>
	!! @brief
	!!
	subroutine center( this )
		class(ClassicalDynamics) :: this
		
		integer :: i
		real(8) :: cm(3)
		
		write(6,'(A)', advance="no") "# Fixing center of mass ... "

		cm = this.centerOfMass()
		do i=1,this.nAtoms
			this.positions(:,i) = this.positions(:,i) - cm(:)
		end do
		
		write(6,'(A,3F10.4,A)') "OK (", this.centerOfMass(), "  )"
	end subroutine center
	
	!*
	! @brief
	!*
	subroutine push( this, velocity, targetAtomId )
		class(ClassicalDynamics) :: this
		real(8), intent(in) :: velocity(3)
		integer, intent(in) :: targetAtomId
		
		integer :: i
		
		write(*,"(A,I8,A)") "# targetAtomId =", targetAtomId, " ("//sym(targetAtomId)//")"
		write(*,"(A,3F8.2,A)") "#     velocity =", velocity/(angs/ps), " angs/ps"
		
		if( targetAtomId == 0 ) then
			do i=1,this.nAtoms
				this.velocities(:, i) = this.velocities(:, i) + velocity
			end do
		else
			this.velocities(:,targetAtomId) = this.velocities(:,targetAtomId) + velocity
		end if
		
	end subroutine push
	
	!*
	! @brief
	!*
	subroutine move( this, position, targetAtomId )
		class(ClassicalDynamics) :: this
		real(8), intent(in) :: position(3)
		integer, intent(in) :: targetAtomId
		
		integer :: i
		
		write(*,"(A,I8,A)") "# targetAtomId =", targetAtomId, " ("//sym(targetAtomId)//")"
		write(*,"(A,3F8.2,A)") "#     position =", position/angs, " angs"
		
		if( targetAtomId == 0 ) then
			do i=1,this.nAtoms
				this.positions(:, i) = this.positions(:, i) + position
			end do
		else
			this.positions(:,targetAtomId) = position
		end if
		
	end subroutine move
	
	!*
	! @brief
	!*
	subroutine thermalize( this, temperature, tau, nSteps, integ, showFreq, confinementRadius, geomOFile, velOFile, status )
		class(ClassicalDynamics) :: this
		real(8), intent(in) :: temperature
		real(8), optional, intent(in) :: tau
		integer, optional, intent(in) :: nSteps
		integer, optional, intent(in) :: integ
		integer, optional, intent(in) :: showFreq
		real(8), optional, intent(in) :: confinementRadius
		character(*), optional :: geomOFile
		character(*), optional :: velOFile
		integer, optional :: status
		
		type(MDIntegrator) :: solver
		integer :: i
		real(8) :: lambda
		
		real(8) :: tauEff
		integer :: nStepsEff
		integer :: integEff
		integer :: showFreqEff
		real(8) :: confinementRadiusEff
		character(:), allocatable :: geomOFileEff
		character(:), allocatable :: velOFileEff
		real(8) :: energyConservationIndex
		
		integer :: j
		real(8) :: d0
		
		if( this.nAtoms < 2 ) then
			write(*,*) "### ERROR ### Termalization process error, nAtoms < 2"
			stop
		end if
		
		tauEff = 2.0_8
		if( present(tau) ) tauEff = tau
		
		nStepsEff = 1
		if( present(nSteps) ) nStepsEff = nSteps
		
		integEff = VELOCITY_VERLET
		if( present(integ) ) integEff = integ
		
		showFreqEff = 1
		if( present(showFreq) ) showFreqEff = showFreq
		
		confinementRadiusEff = -1.0_8
		if( present(confinementRadius) ) confinementRadiusEff = confinementRadius
		
		geomOFileEff = "geom.xyz"
		if( present(geomOFile) ) geomOFileEff = geomOFile
		
		velOFileEff = "vel.xyz"
		if( present(velOFile) ) velOFileEff = velOFile
		
		if( this.historyStep == 0 ) then
			call system( "rm -rf "//trim(geomOFile)//" 2> /dev/null" )
			call system( "rm -rf "//trim(velOFile)//" 2> /dev/null" )
		end if
		
		call solver.init( this.nAtoms, this.positions, this.velocities, this.accelerations, this.forces )
		solver.timeStep = this.timeStep
		solver.ttype = integEff
		
		write(*,*)  ""
		write(*,"(A15,A20,3A20,A10,4X,A10,3X,3A8)") "step", "t(ps)", "Ek(amuA^2/ps^2)", "Ep(amuA^2/ps^2)", &
			"Et(amuA^2/ps^2)", "ci", "T(K)", "Xcm(A)", "Ycm(A)", "Zcm(A)"
			
		write(*,"(A15,A20,3A20,A10,4X,A10,3X,3A8)") "----", "-----", "---------------", "---------------", &
			"---------------", "--", "----", "------", "------", "------"
		
		do i=0,nStepsEff
			
			call this.updateNeighbourList()
			
			if( mod(i,showFreqEff) == 0 ) then
			
				if( i/=0 ) then
					this.averEnergy = ( this.averEnergy*real((i-1)/showFreqEff,8) + this.totalEnergy/this.nAtoms )/real(i/showFreqEff,8)
					energyConservationIndex = abs(this.totalEnergy/this.nAtoms-this.averEnergy)/this.totalEnergy/this.nAtoms
				else
					call computeForces()
					this.averEnergy = this.totalEnergy/this.nAtoms
					energyConservationIndex = 0.0_8
				end if
				
				write(*,"(I15,F20.5,3F20.7,E10.1,4X,F10.3,3X,3F8.3,F8.3)") &
					solver.step, &
					this.historyTimeStep/ps, &
					this.kineticEnergy/this.nAtoms/(amu*angs**2/ps**2), &
					this.potentialEnergy/this.nAtoms/(amu*angs**2/ps**2), &
					this.totalEnergy/this.nAtoms/(amu*angs**2/ps**2), &
					energyConservationIndex, &
					this.temperature()/kelvin, &
					this.centerOfMass()/angs, lambda
					
				if( this.historyStep == 0 ) then
					call this.saveGeometry( geomOFileEff )
					call this.saveVelocities( velOFileEff )
				else
					call this.saveGeometry( geomOFileEff, .true. )
					call this.saveVelocities( velOFileEff, .true. )
				end if
				
			end if
			
			if( solver.step > 0 ) then
				if( temperature > this.temperature() ) then
					lambda = sqrt(1.0_8+(solver.timeStep/tauEff)*(temperature/this.temperature()-1.0_8))
				else
					lambda = 0.99_8
				end if
			else
				lambda = 1.0_8
			end if
			
			this.velocities = lambda*this.velocities
			
			do j=1,this.nAtoms
				d0 = sqrt( sum( this.positions(:,j)**2 ) )
				if( confinementRadiusEff > 0.0_8 .and. d0 > confinementRadiusEff ) then
 					this.positions(:,j) = 0.9*confinementRadiusEff*this.positions(:,j)/norm2(this.positions(:,j))
					this.velocities(:,j) = -0.9*this.velocities(:,j)
				end if
			end do
			
			call solver.iterate( computeForces )
			
! 			if( this.kineticEnergy/this.nAtoms/(amu*angs**2/ps**2) > 1000.0_8 ) then
! 				status = 1
! 				return
! 			end if
			
			this.historyStep = this.historyStep + 1
			this.historyTimeStep = this.historyTimeStep + this.timeStep
		end do
		
		call solver.destroy()
		status = 0
	end subroutine thermalize
	
	!*
	! @brief
	!*
	subroutine freeEvolution( this, nSteps, integ, showFreq, geomOFile, velOFile )
		class(ClassicalDynamics) :: this
		integer, optional, intent(in) :: nSteps
		integer, optional, intent(in) :: integ
		integer, optional, intent(in) :: showFreq
		character(*), optional :: geomOFile
		character(*), optional :: velOFile
		
		type(MDIntegrator) :: solver
		integer :: i
		
		integer :: nStepsEff
		integer :: integEff
		integer :: showFreqEff
		character(:), allocatable :: geomOFileEff
		character(:), allocatable :: velOFileEff
		real(8) :: energyConservationIndex
		
		if( present(nSteps) ) then
			nStepsEff = nSteps
		else
			nStepsEff = 1
		end if
		
		if( present(integ) ) then
			integEff = integ
		else
			integEff = VELOCITY_VERLET
		end if
		
		if( present(showFreq) ) then
			showFreqEff = showFreq
		else
			showFreqEff = 1
		end if
		
		if( present(geomOFile) ) then
			geomOFileEff = geomOFile
		else
			geomOFileEff = "geom.xyz"
		end if
		
		if( present(velOFile) ) then
			velOFileEff = velOFile
		else
			velOFileEff = "vel.xyz"
		end if
		
		if( this.historyStep == 0 ) then
			call system( "rm -rf "//trim(geomOFile)//" 2> /dev/null" )
			call system( "rm -rf "//trim(velOFile)//" 2> /dev/null" )
		end if
		
		call solver.init( this.nAtoms, this.positions, this.velocities, this.accelerations, this.forces )
		solver.timeStep = this.timeStep
		solver.ttype = integEff
		
		write(*,*)  ""
		write(*,"(A15,A20,3A20,A10,4X,A10,3X,3A8)") "step", "t(ps)", "Ek(amuA^2/ps^2)", "Ep(amuA^2/ps^2)", &
			"Et(amuA^2/ps^2)", "ci", "T(K)", "Xcm(A)", "Ycm(A)", "Zcm(A)"
			
		write(*,"(A15,A20,3A20,A10,4X,A10,3X,3A8)") "----", "-----", "---------------", "---------------", &
			"---------------", "--", "----", "------", "------", "------"
		
		do i=0,nStepsEff
			
			call this.updateNeighbourList()
			
			if( mod(i,showFreqEff) == 0 .and. i/=0 ) then
				
				if( i/=0 ) then
					this.averEnergy = ( this.averEnergy*real((i-1)/showFreqEff,8) + this.totalEnergy/this.nAtoms )/real(i/showFreqEff,8)
					energyConservationIndex = abs(this.totalEnergy/this.nAtoms-this.averEnergy)/this.totalEnergy/this.nAtoms
				else
					call computeForces()
					this.averEnergy = this.totalEnergy/this.nAtoms
					energyConservationIndex = 0.0_8
				end if
				
				write(*,"(I15,F20.5,3F20.7,E10.1,4X,F10.3,3X,3F8.3)") &
					solver.step, &
					this.historyTimeStep/ps, &
					this.kineticEnergy/this.nAtoms/(amu*angs**2/ps**2), &
					this.potentialEnergy/this.nAtoms/(amu*angs**2/ps**2), &
					this.totalEnergy/this.nAtoms/(amu*angs**2/ps**2), &
					energyConservationIndex, &
					this.temperature()/kelvin, &
					this.centerOfMass()/angs
					
				if( this.historyStep == 0 ) then
					call this.saveGeometry( geomOFileEff )
					call this.saveVelocities( velOFileEff )
				else
					call this.saveGeometry( geomOFileEff, .true. )
					call this.saveVelocities( velOFileEff, .true. )
				end if
				
			end if
			
			call solver.iterate( computeForces )
			
			this.historyStep = this.historyStep + 1
			this.historyTimeStep = this.historyTimeStep + this.timeStep
		end do
		
		call solver.destroy()
		
	end subroutine freeEvolution
	
	!>
	!! @brief Saves the geometry
	!!
	subroutine saveGeometry( this, fileName, append )
		class(ClassicalDynamics) :: this
		character(*), intent(in) :: fileName
		logical, optional, intent(in) :: append
		
		character(10) :: accessEff
		integer :: i
		
		accessEff = "sequential"
		
		if( present(append) .and. append ) then
			accessEff = "append"
		end if
		
		open( 10, file=fileName, access=accessEff, status="unknown" )
		write(10,'(I10)') size(this.positions,dim=2)
		write(10,'(F15.7,A10,I10,A15,F15.7)') this.totalEnergy, "step = ", this.historyStep, "time = ", this.historyTimeStep/ps
		do i=1,size(this.positions,dim=2)
			write(10,'(A5,3F20.5)') this.symbols(i), this.positions(1,i)/angs, this.positions(2,i)/angs, this.positions(3,i)/angs
		end do
		close(10)
	end subroutine saveGeometry
	
	!>
	!! @brief Saves the velocities
	!!
	subroutine saveVelocities( this, fileName, append )
		class(ClassicalDynamics) :: this
		character(*), intent(in) :: fileName
		logical, optional, intent(in) :: append
		
		character(10) :: accessEff
		integer :: i
		
		accessEff = "sequential"
		
		if( present(append) .and. append ) then
			accessEff = "append"
		end if
		
		open( 10, file=fileName, access=accessEff, status="unknown" )
		write(10,'(I10)') size(this.positions,dim=2)
		write(10,'(F15.7,A10,I10,A15,F15.7)') this.totalEnergy, "step = ", this.historyStep, "time = ", this.historyTimeStep/ps
		do i=1,size(this.velocities,dim=2)
			write(10,'(A3,3F20.5)') this.symbols(i), this.velocities(1,i)/(angs/ps), this.velocities(2,i)/(angs/ps), this.velocities(3,i)/(angs/ps)
		end do
		close(10)
	end subroutine saveVelocities
	
! 	!*
! 	! @brief Saves the pair distribution function
! 	!*
! 	subroutine savePDF( this, fileName, dr, rMax )
! 		class(ClassicalDynamics) :: this
! 		character(*), intent(in) :: fileName
! 		real(8), optional, intent(in) :: dr
! 		real(8), optional, intent(in) :: rMax
! 		
! 		real(8) :: drEff
! 		real(8) :: rMaxEff
! 		
! 		integer :: nGroups
! 		real(8), allocatable :: gk(:)
! 		real(8) :: ri, rij
! 		integer :: k_ij
! 		integer :: i, j
! 		real(8) :: sumgk
! 		
! 		real(8) :: pi = acos(-1.0_8)
! 		
! 		if( present(dr) ) then
! 			drEff = dr
! 		else
! 			drEff = 0.5_8
! 		end if
! 		
! 		if( present(rMax) ) then
! 			rMaxEff = rMax
! 		else
! 			rMaxEff = 50.0_8
! 		end if
! 		
! 		nGroups = int(rMaxEff/drEff)
! 		allocate( gk( nGroups ) )
! 		gk = 0.0_8
! 		
! 		do i=1,this.nAtoms
! 			do j=i+1,this.nAtoms
! 				rij = sqrt(sum((this.positions(:,j)-this.positions(:,i))**2))
! 				k_ij = ceiling( rij/drEff )
! 				
! 				if( k_ij <= nGroups ) then
! 					gk( k_ij ) = gk( k_ij ) + 1.0_8
! 				end if
! 			end do
! 		end do
! 		
! 		open( 10, file=fileName )
! 		do i=1,nGroups
! 			ri = (real(i,8)+0.5)*drEff
! 			write(10,'(2F20.5)') ri, gk(i)/ri**2
! 		end do
! 		close(10)
! 		
! 		deallocate( gk )
! 	end subroutine savePDF
! 	
! 	!*
! 	! @brief Saves lateral distribution function
! 	!*
! 	subroutine saveLDF( this, fileName, dr, rMin, rMax )
! 		class(ClassicalDynamics) :: this
! 		character(*), intent(in) :: fileName
! 		real(8), optional, intent(in) :: dr
! 		real(8), optional, intent(in) :: rMin
! 		real(8), optional, intent(in) :: rMax
! 		
! 		real(8) :: drEff
! 		real(8) :: rMinEff
! 		real(8) :: rMaxEff
! 		
! 		integer :: nGroups
! 		real(8), allocatable :: D(:)
! 		integer :: z_i
! 		integer :: k_i
! 		integer :: i
! 		real(8) :: sumD
! 		
! 		if( present(dr) ) then
! 			drEff = dr
! 		else
! 			drEff = 3.0_8
! 		end if
! 		
! 		if( present(rMin) ) then
! 			rMinEff = rMin
! 		else
! 			rMinEff = 0.0_8
! 		end if
! 		
! 		if( present(rMax) ) then
! 			rMaxEff = rMax
! 		else
! 			rMaxEff = 100.0_8
! 		end if
! 		
! 		nGroups = int((rMaxEff-rMinEff)/drEff)
! 		
! 		allocate( D( nGroups ) )
! 		D = 0.0_8
! 		
! 		do i=1,this.nAtoms
! 			z_i = this.positions(3,i)
! 			k_i = ceiling( ( z_i-rMinEff )/drEff )
! 			
! 			if( k_i <= nGroups ) then
! 				D( k_i ) = D( k_i ) + 1.0_8
! 			end if
! 		end do
! 		
! 		sumD = sum(D)*drEff
! 		D = D/sumD
! 		
! 		open( 10, file=fileName )
! 		do i=1,nGroups
! 			z_i = rMinEff + ( real(i-1,8)+0.5_8 )*drEff
! 			write(10,'(2F20.5)') z_i, D(i)
! 		end do
! 		close(10)
! 		
! 		deallocate( D )
! 	end subroutine saveLDF
! 	
! 	!*
! 	! @brief Saves angular distribution function
! 	!*
! 	subroutine saveADF( this, fileName, dtheta, thetaMin, thetaMax, origin )
! 		class(ClassicalDynamics) :: this
! 		character(*), intent(in) :: fileName
! 		real(8), optional, intent(in) :: dtheta
! 		real(8), optional, intent(in) :: thetaMin
! 		real(8), optional, intent(in) :: thetaMax
! 		real(8), optional, intent(in) :: origin(3)
! 		
! 		real(8) :: dthetaEff
! 		real(8) :: thetaMinEff
! 		real(8) :: thetaMaxEff
! 		real(8) :: originEff(3)
! 		
! 		real(8) :: theta_i
! 		integer :: nGroups
! 		real(8), allocatable :: D(:)
! 		real(8), allocatable :: Dsin(:)
! 		real(8) :: r_i
! 		integer :: k_i
! 		integer :: i
! 		real(8) :: sumDsin
! 		
! 		real(8) :: pi = acos(-1.0_8)
! 		
! 		if( present(dtheta) ) then
! 			dthetaEff = dtheta
! 		else
! 			dthetaEff = pi/40.0_8
! 		end if
! 		
! 		if( present(thetaMin) ) then
! 			thetaMinEff = thetaMin
! 		else
! 			thetaMinEff = 0.0_8
! 		end if
! 		
! 		if( present(thetaMax) ) then
! 			thetaMaxEff = thetaMax
! 		else
! 			thetaMaxEff = pi
! 		end if
! 		
! 		if( present(origin) ) then
! 			originEff = origin
! 		else
! 			originEff = [ 0.0_8, 0.0_8, -27.4_8 ]
! ! 			originEff = [ 0.0_8, 0.0_8, 0.0_8 ]
! 		end if
! 		
! 		nGroups = int((thetaMaxEff-thetaMinEff)/dthetaEff)
! 		
! 		allocate( D( nGroups ) )
! 		D = 0.0_8
! 		
! 		allocate( Dsin( nGroups ) )
! 		Dsin = 0.0_8
! 		
! 		do i=1,this.nAtoms
! 			r_i = sqrt( sum(( this.positions(:,i)-originEff )**2) )
! 			theta_i = acos( ( this.positions(3,i)-originEff(3) )/r_i )
! 			k_i = ceiling( theta_i/dthetaEff )
! 			
! 			if( k_i <= nGroups ) then
! 				D( k_i ) = D( k_i ) + 1.0_8/sin(theta_i)
! 				Dsin( k_i ) = D( k_i )*sin(theta_i)
! 			end if
! 		end do
! 		
! 		sumDsin = sum(Dsin)*dthetaEff
! 		D = D/sumDsin
! 		
! 		open( 10, file=fileName )
! 		write(10,'(A,3F20.5)') "#", originEff
! 		do i=1,nGroups
! 			theta_i = thetaMinEff + ( real(i-1,8)+0.5_8 )*dthetaEff
! 			write(10,'(2F20.5)') theta_i, D(i)
! 		end do
! 		close(10)
! 		
! 		deallocate( D )
! 	end subroutine saveADF
	
	!*
	! @brief Test method
	!*
	subroutine ClassicalDynamics_test()
		type(HeDroplet) :: droplet
		type(ClassicalDynamics) :: cDroplet
		type(HeHePotential) :: HeHePot
		
! 		call cDroplet.init( 300, droplet )
! 		call cDroplet.init( "init.xyz", XYZ )
! 		call cDroplet.saveGeometry( "geometry0.xyz" )
! 		call cDroplet.saveRadialDistribution( "radialDist.dat" )
! 		call cDroplet.thermalize( 0.107_8, tau=10.0_8, nSteps=100000, showFreq=500, geomOFile="salida.xyz" )
! 		call cDroplet.thermalize( 0.107_8, tau=10.0_8, nSteps=10000, showFreq=500, geomOFile="salida.xyz" )
! 		call cDroplet.freeEvolution( nSteps=500000, showFreq=1000, geomOFile="salida.xyz" )
! 		call cDroplet.saveRadialDistribution( "radialDistEnd.dat" )
! 		call cDroplet.init( "end.xyz", XYZ )
		
	end subroutine ClassicalDynamics_test
	
end module ClassicalDynamics_
