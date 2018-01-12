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

program ccdyn_
	use UnitsConverter_
	
	use HeDroplet_
	use ClassicalDynamics_
	use MDIntegrator_
	use HeTiO2Potential_
	use AuTiO2Potential_
	use EmbeddingInterface_
	implicit none
	
	integer, parameter :: SIZE_NML_VEC = 10
	integer, parameter :: SIZE_NML_STRING = 255
	integer, parameter :: INPUT_UNIT = 8
	
	type(HeDroplet) :: droplet
	type(ClassicalDynamics) :: cDroplet
	type(HeTiO2Potential) :: HeTiO2Pot
	type(AuTiO2Potential) :: AuTiO2Pot
	
	character(SIZE_NML_STRING) :: inputFileName
	character(SIZE_NML_STRING) :: outputPrefix
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! data
	integer :: idLabel
	character(SIZE_NML_STRING) :: initialDensity
	integer :: nAtoms
	character(SIZE_NML_STRING) :: initialGeom
	character(SIZE_NML_STRING) :: initialVel
	character(SIZE_NML_STRING) :: geometryHistoryFile
	character(SIZE_NML_STRING) :: velocityHistoryFile
	character(SIZE_NML_STRING) :: actions
	logical :: neighbourList
	real(8) :: neighbourRadius
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Thermalizing and FreeEvolution
	integer :: nSteps(SIZE_NML_VEC)
	real(8) :: timeStep(SIZE_NML_VEC)
	integer :: showFrequency(SIZE_NML_STRING)
	
	real(8) :: temperature(SIZE_NML_STRING)
	real(8) :: tau(SIZE_NML_STRING)
	character(SIZE_NML_STRING) :: integrator(SIZE_NML_STRING)
	real(8) :: confinementRadius
	
	integer :: currentThermalizingStep
	integer :: currentFreeEvolutionStep
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Push and move
	integer :: targetAtomId(SIZE_NML_VEC)
	real(8) :: position(3*SIZE_NML_VEC)
	real(8) :: velocity(3*SIZE_NML_VEC)
	
	integer :: currentPushStep
	integer :: currentMoveStep
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Center
	character(SIZE_NML_STRING) :: targetSymbol
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! ExternalField
	character(2) :: targetAtom(SIZE_NML_VEC)
	character(SIZE_NML_STRING) :: name(SIZE_NML_VEC)
	real(8) :: origin(3*SIZE_NML_VEC)

	integer :: nTargetAtoms
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Save
	character(SIZE_NML_STRING) :: geometry(SIZE_NML_VEC)
	character(SIZE_NML_STRING) :: velocities(SIZE_NML_VEC)
	character(SIZE_NML_STRING) :: PDF(SIZE_NML_VEC)  ! Pair distribution function
	character(SIZE_NML_STRING) :: LDF(SIZE_NML_VEC)  ! Lateral distribution function
	character(SIZE_NML_STRING) :: ADF(SIZE_NML_VEC)  ! Angular distribution function
	
	integer :: currentSaveStep
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Embedding
	integer :: targetAtoms(SIZE_NML_VEC)
	real(8) :: energyReference
	character(1000) :: command
	character(1000) :: energyFilter
	character(1000) :: gradientFilter
	character(1000) :: geometryConverter
	integer :: sizeTargetAtoms
	type(EmbeddingInterface) :: embed
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	integer :: i, j, effSize, intBuffer, status
	character(255) :: strBuffer
	
	if( command_argument_count() < 1 ) then
		write(*,*) "Usage: ccdyn inputFile.in"
		stop
	end if
	
	call get_command_argument( 1, inputFileName )
	open( INPUT_UNIT, file=trim(inputFileName) )
	
	write(*,*) ""
	write(*,*) " inputFileName = ", trim(inputFileName)
	write(*,*) ""
	
	!---------------------------------
	targetAtoms = 0
	command = ""
	gradientFilter = ""
	geometryConverter = ""
	
	namelist /embedding/ targetAtoms, command, energyFilter, gradientFilter, geometryConverter, energyReference
	read( INPUT_UNIT, embedding )
	rewind( INPUT_UNIT )
	
	if( any(targetAtoms /= 0) ) then
		write(*,*) ""
		write(*,*) "> EMBEDDING"
		write(*,*) ""
		
		sizeTargetAtoms=0
		do i=1,SIZE_NML_VEC
			if( targetAtoms(i) > 0 ) &
				sizeTargetAtoms = sizeTargetAtoms + 1
		end do
		
		if( size(targetAtoms) > 1 ) then
			write(*,"(2X,A,<sizeTargetAtoms>I)") " targetAtoms = ", targetAtoms(1:sizeTargetAtoms)
			write(*,"(2X,A,A)") " command = ", trim(command)
			write(*,"(2X,A,A)") " energyFilter = ", trim(energyFilter)
			write(*,"(2X,A,A)") " gradientFilter = ", trim(gradientFilter)
			write(*,"(2X,A,A)") " geometryConverter = ", trim(geometryConverter)
			write(*,"(2X,A,F20.5)") " energyReference = ", energyReference
		end if
		
		write(*,*) ""
		
		call embed.init( sizeTargetAtoms, targetAtoms, command, energyFilter, gradientFilter, geometryConverter, energyReference )
		
		write(*,*) ""
	end if
	!---------------------------------
	
	idLabel = 0
	initialDensity = ""
	initialGeom = ""
	initialVel = ""
	actions = "F"
	neighbourList = .true.
	neighbourRadius = 10.0_8
	namelist /data/ idLabel, initialDensity, nAtoms, initialGeom, initialVel, actions, neighbourList, neighbourRadius
	read( INPUT_UNIT, data )
	rewind( INPUT_UNIT )
	
	neighbourRadius = neighbourRadius*angs
	
	if( idLabel == 0 ) then
		outputPrefix = ""
	else
		write(*,*) " idLabel = ", idLabel
	
		write(strBuffer,'(I10.4)') idLabel
		outputPrefix = trim(strBuffer)//"/"
		call system( "if [ ! -d "//trim(outputPrefix)//" ]; then mkdir "//trim(outputPrefix)//"; fi" )
	end if
	
	write(*,*) " neighbourList = ", neighbourList
	write(*,*) " neighbourRadius = ", neighbourRadius/angs, " angs"
	write(*,*) ""
	write(*,*) "> INITIALING SYSTEM"
	write(*,*) ""
	if( len_trim(initialDensity) > 0 ) then
	
		write(*,"(2X,A,10I)") " nAtoms = ", nAtoms
		write(*,"(2X,A,30A)") " initialDensity = ", trim(initialDensity)
		
		call droplet.init( trim(initialDensity) )
		
		if( any(targetAtoms /= 0) ) then
			call cDroplet.init( nAtoms, droplet, neighbourList, neighbourRadius )
		else
			call cDroplet.init( nAtoms, droplet, neighbourList, neighbourRadius, embed )
		end if
		
	else if( len_trim(initialGeom) > 0 ) then
		write(*,*) " initialGeom = ", trim(initialGeom)
		
		if( any(targetAtoms /= 0) ) then
			if( len_trim(initialVel) > 0 ) then
				write(*,*) " initialVel = ", trim(initialVel)
				call cDroplet.init( trim(initialGeom), velFileName=trim(initialVel), useNeighbourList=neighbourList, &
										neighbourRadius=neighbourRadius, embedding=embed )
			else
				call cDroplet.init( trim(initialGeom), useNeighbourList=neighbourList, &
										neighbourRadius=neighbourRadius, embedding=embed )
			end if
		else
			if( len_trim(initialVel) > 0 ) then
				write(*,*) " initialVel = ", trim(initialVel)
				call cDroplet.init( trim(initialGeom), velFileName=trim(initialVel), useNeighbourList=neighbourList, &
										neighbourRadius=neighbourRadius )
			else
				call cDroplet.init( trim(initialGeom), useNeighbourList=neighbourList, &
										neighbourRadius=neighbourRadius )
			end if
		end if
	end if
	write(*,*) " actions = ", trim(actions)
	!---------------------------------
	
	currentThermalizingStep = 1
	currentFreeEvolutionStep = 1
	currentPushStep = 1
	currentMoveStep = 1
	currentSaveStep = 1
	
	do i=1,len(trim(actions))
		select case( actions(i:i) )
			case('T')
				nSteps = 0
				nSteps(1) = 1000
				timeStep = 0.1*ps
				temperature = 375.0_8*kelvin
				tau = 0.1_8*ps
				showFrequency = 10
				confinementRadius = -1.0_8
				geometryHistoryFile = "history.xyz"
				velocityHistoryFile = "history.vxyz"
				integrator = "VelocityVerlet"
				
				namelist /thermalize/ nSteps, timeStep, temperature, tau, showFrequency, confinementRadius, geometryHistoryFile, velocityHistoryFile, integrator
				read( INPUT_UNIT, thermalize )
				rewind( INPUT_UNIT )
				
				timeStep = timeStep*ps
				temperature = temperature*kelvin
! 				tau = tau
				confinementRadius = confinementRadius*angs
				
				do while( .true. )

					write(*,*) ""
					write(*,*) "> THERMALIZE"
					write(*,*) ""
					write(*,"(2X,A30,I10)") "nSteps = ", nSteps( currentThermalizingStep )
					write(*,"(2X,A30,F10.5,A5)") "timeStep = ", timeStep( currentThermalizingStep )/ps, "ps"
					write(*,"(2X,A30,F10.5,A5)") "temperature = ", temperature( currentThermalizingStep )/kelvin, "K"
					write(*,"(2X,A30,E10.5)") "tau = ", tau( currentThermalizingStep )
					write(*,"(2X,A30,I10)") "showFrequency = ", showFrequency( currentThermalizingStep )
					write(*,"(2X,A30,F10.5,A5)") "confinementRadius = ", confinementRadius/angs, "A"
					write(*,"(2X,A30,A30)") "geometryHistoryFile = ", trim(geometryHistoryFile)
					write(*,"(2X,A30,A30)") "velocityHistoryFile = ", trim(velocityHistoryFile)
					write(*,"(2X,A30,A30,A,I1,A)") "integrator = ", trim(integrator( currentThermalizingStep )), &
						"(", MDIntegrator_charIDtoID( trim(integrator( currentThermalizingStep )) ), ")"
					
					cDroplet.timeStep = timeStep( currentThermalizingStep )
				
					call cDroplet.thermalize( temperature( currentThermalizingStep ), &
						tau( currentThermalizingStep ), nSteps( currentThermalizingStep ), &
						MDIntegrator_charIDtoID( trim(integrator( currentThermalizingStep )) ), &
						showFrequency( currentThermalizingStep ), &
						confinementRadius, &
						trim(outputPrefix)//trim(geometryHistoryFile), &
						trim(outputPrefix)//trim(velocityHistoryFile), status )
						
					if( status == 0 ) then
						exit
					else
						write(*,*) ""
						write(*,*) "> Unstable geometry. New trial"
						write(*,*) ""
					end if
					
				end do
				
				currentThermalizingStep = currentThermalizingStep + 1
				
			case('F')
				nSteps = 0
				nSteps(1) = 1000
				showFrequency = 100
				geometryHistoryFile = "history.xyz"
				velocityHistoryFile = "history.vxyz"
				integrator = "VelocityVerlet"
				
				namelist /freeEvolution/ nSteps, timeStep, showFrequency, geometryHistoryFile, velocityHistoryFile, integrator
				read( INPUT_UNIT, freeEvolution )
				rewind( INPUT_UNIT )
				
				timeStep = timeStep*ps
				
				write(*,*) ""
				write(*,*) "> FREE EVOLUTION"
				write(*,*) ""
				write(*,"(2X,A30,I10)") "nSteps = ", nSteps( currentFreeEvolutionStep )
				write(*,"(2X,A30,F10.5,A5)") "timeStep = ", timeStep( currentFreeEvolutionStep )/ps, "ps"
				write(*,"(2X,A30,I10)") "showFrequency = ", showFrequency( currentFreeEvolutionStep )
				write(*,"(2X,A30,A30)") "geometryHistoryFile = ", trim(geometryHistoryFile)
				write(*,"(2X,A30,A30)") "velocityHistoryFile = ", trim(velocityHistoryFile)
				write(*,"(2X,A30,A30,A,I1,A)") "integrator = ", trim(integrator( currentFreeEvolutionStep )), &
					"(", MDIntegrator_charIDtoID( trim(integrator( currentFreeEvolutionStep )) ), ")"
					
				cDroplet.timeStep = timeStep( currentFreeEvolutionStep )
				call cDroplet.freeEvolution( nSteps( currentFreeEvolutionStep ), &
					MDIntegrator_charIDtoID( trim(integrator( currentFreeEvolutionStep )) ), &
					showFrequency( currentFreeEvolutionStep ), trim(outputPrefix)//trim(geometryHistoryFile), &
					trim(outputPrefix)//trim(velocityHistoryFile) )
				
				currentFreeEvolutionStep = currentFreeEvolutionStep + 1
				
			case('E')
				targetAtom = "XX"
				name = "XX"
				origin = 0.0_8
				namelist /externalField/ targetAtom, name, origin
				read( INPUT_UNIT, externalField )
				rewind( INPUT_UNIT )
				
				origin = origin*angs
				
				do j=1,SIZE_NML_STRING
					if( targetAtom(j) == "XX" ) then
						nTargetAtoms = j-1
						exit
					end if
				end do
				
				write(*,*) ""
				write(*,*) "> EXTERNAL FIELD APPLIED ( ", nTargetAtoms, " )"
				write(*,*) ""
				do j=1,nTargetAtoms
					write(*,"(2X,A,A10)") "targetAtom = ", targetAtom(j)
					write(*,"(2X,A,A)") "name = ", trim(name(j))
					write(*,"(2X,A,3F10.5,A5)") "origin = ", origin(3*(j-1)+1:3*(j-1)+3)/angs, "A"
					write(*,*) ""
					
					call initExternalPotential( targetAtom(j), name(j) )
				end do
				
				call cDroplet.setExternalField( V, dV )
				
			case('C')
				write(*,*) ""
				write(*,*) "> CENTER"
				write(*,*) ""
				
				call cDroplet.center()
				
			case('P')
				targetAtomId = 0
				velocity = 0.0_8
				namelist /push/ targetAtomId, velocity
				read( INPUT_UNIT, push )
				rewind( INPUT_UNIT )
				
				velocity = velocity*(angs/ps)
				
				write(*,*) ""
				write(*,*) "> PUSH"
				write(*,*) ""
				write(*,"(2X,A,I10)") "targetAtomId = ", targetAtomId(currentPushStep)
				write(*,"(2X,A,3F10.5,A5)") "    velocity = ", velocity(3*(currentPushStep-1)+1:3*(currentPushStep-1)+3)/(angs/ps), "A/ps"
				write(*,*) ""
				
				call cDroplet.push( velocity(3*(currentPushStep-1)+1:3*(currentPushStep-1)+3), targetAtomId(currentPushStep) )
				
				currentPushStep = currentPushStep + 1
				
			case('M')
				targetAtomId = 0
				position = 0.0_8
				namelist /move/ targetAtomId, position
				read( INPUT_UNIT, move )
				rewind( INPUT_UNIT )
				
				position = position*angs
				
				write(*,*) ""
				write(*,*) "> MOVE"
				write(*,*) ""
				write(*,"(2X,A,3F10.5,A5)") "position = ", position(3*(currentMoveStep-1)+1:3*(currentMoveStep-1)+3)/angs, "A"
				
				call cDroplet.move( position(3*(currentMoveStep-1)+1:3*(currentMoveStep-1)+3), targetAtomId(currentMoveStep) )
				
				currentMoveStep = currentMoveStep + 1
				
			case('S')
				geometry = ""
				geometry(1) = "geometry00.xyz"
				velocities = ""
				PDF = ""
				ADF = ""
				LDF = ""
				
				namelist /save/ geometry, velocities, PDF, LDF, ADF
				read( INPUT_UNIT, save )
				rewind( INPUT_UNIT )
				
				write(*,*) ""
				write(*,*) "> SAVE GEOMETRY"
				write(*,*) ""
				
				write(*,"(2X,20A,A20)") "geometry = ", trim(geometry(currentSaveStep))
				call cDroplet.saveGeometry( trim(outputPrefix)//trim(geometry(currentSaveStep)) )
				
				if( len_trim(velocities(currentSaveStep)) /= 0 ) then
					write(*,"(2X,20A,A20)") "velocities = ", trim(velocities(currentSaveStep))
					call cDroplet.saveVelocities( trim(outputPrefix)//trim(velocities(currentSaveStep)) )
				end if
					
! 				if( len_trim(PDF(currentSaveStep)) /= 0 ) then
! 					write(*,"(2X,20A,A20)") "PDF = ", trim(PDF(currentSaveStep))
! 					call cDroplet.savePDF( trim(outputPrefix)//trim(PDF(currentSaveStep)) )
! 				end if
! 					
! 				if( len_trim(LDF(currentSaveStep)) /= 0 ) then
! 					write(*,"(2X,20A,A20)") "LDF = ", trim(LDF(currentSaveStep))
! 					call cDroplet.saveLDF( trim(outputPrefix)//trim(LDF(currentSaveStep)) )
! 				end if
! 				
! 				if( len_trim(ADF(currentSaveStep)) /= 0 ) then
! 					write(*,"(2X,20A,A20)") "ADF = ", trim(ADF(currentSaveStep))
! 					call cDroplet.saveADF( trim(outputPrefix)//trim(ADF(currentSaveStep)) )
! 				end if
				
				currentSaveStep = currentSaveStep + 1
		end select
	end do
	
	contains
	
	subroutine initExternalPotential( symbol, name )
		character(*), intent(in) :: symbol
		character(*), intent(in) :: name
		
		if( trim(symbol) == "He" ) then
			select case( trim(name) )
				case( "HeTiO2_LAP" )
					call HeTiO2Pot.init( LAP )
				case( "HeTiO2_CMP" )
					call HeTiO2Pot.init( CMP )
				case( "HeTiO2_M3DP" )
					call HeTiO2Pot.init( M3DP )
				case( "HeTiO2_LAPC3LJ" )
					call HeTiO2Pot.init( LAPC3LJ )
				case( "HeTiO2_LAPC3LJv2" )
					call HeTiO2Pot.init( LAPC3LJv2 )
				case default
					write(*,*) "### ERROR ### External potential ", trim(name) , " for ", trim(symbol)," is not implemented"
					stop
			end select
		else if( trim(symbol) == "Au" ) then
			select case( trim(name) )
				case( "AuTiO2_SM" )
					call AuTiO2Pot.init( SM )
				case( "AuTiO2_SM1Dz" )
					call AuTiO2Pot.init( SM1Dz )
				case( "AuTiO2_SMSP" )
					call AuTiO2Pot.init( SMSP )
				case( "AuTiO2_SMSP1Dz" )
					call AuTiO2Pot.init( SMSP1Dz )
				case( "AuTiO2_SMSPDisp1Dz" )
					call AuTiO2Pot.init( SMSPDisp1Dz )
				case default
					write(*,*) "### ERROR ### External potential ", trim(name) , " for ", trim(symbol)," is not implemented"
					stop
			end select
		end if
	end subroutine initExternalPotential
	
	function V( symbol, x, y, z ) result( output )
		character(*), intent(in) :: symbol
		real(8), intent(in) :: x
		real(8), intent(in) :: y
		real(8), intent(in) :: z
		real(8) :: output
		
		if( trim(symbol) == "He" ) then
			output = HeTiO2Pot.V( (x-origin(1))/angs, (y-origin(2))/angs, (z-origin(3))/angs )*cm1
		else if( trim(symbol) == "Au" ) then
			output = AuTiO2Pot.V( (x-origin(1))/angs, (y-origin(2))/angs, (z-origin(3))/angs )*cm1
		else
			write(*,*) "### ERROR ### External potential for ", trim(symbol)," is not implemented"
			stop
		end if
	end function V

	function dV( symbol, x, y, z, coord ) result( output )
		character(*), intent(in) :: symbol
		real(8), intent(in) :: x
		real(8), intent(in) :: y
		real(8), intent(in) :: z
		integer, intent(in) :: coord
		real(8) :: output
		
		if( trim(symbol) == "He" ) then
			output = HeTiO2Pot.dV( (x-origin(1))/angs, (y-origin(2))/angs, (z-origin(3))/angs, coord )*(cm1/angs)
		else if( trim(symbol) == "Au" ) then
			output = AuTiO2Pot.dV( (x-origin(1))/angs, (y-origin(2))/angs, (z-origin(3))/angs, coord )*(cm1/angs)
		else
			write(*,*) "### ERROR ### External potential for ", trim(symbol)," is not implemented"
			stop
		end if
	end function dV

end program ccdyn_
