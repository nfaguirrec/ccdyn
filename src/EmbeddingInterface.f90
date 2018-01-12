!!**********************************************************************************
!!                                                                                 !
!!    Authors:                                                                     !
!!    (2015-2015) Néstor F. Aguirre                                                !
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

!>
!! @brief
!!
module EmbeddingInterface_
	use UnitsConverter_
	
	implicit none
	private
	
	type, public :: EmbeddingInterface
		integer, allocatable :: targetAtoms(:)
		real(8) :: energy
		real(8), allocatable :: gradient(:,:)
		character(1000) :: command
		character(1000) :: energyFilter
		character(1000) :: gradientFilter
		character(1000) :: geometryConverter
		real(8) :: energyReference
		
		contains
			generic :: init => initEmbedding
			generic :: assignment(=) => copyEmbedding
			
			procedure :: initEmbedding
			procedure :: copyEmbedding
			final :: destroyEmbedding
			
			procedure :: update
	end type EmbeddingInterface
	
	contains
	
	!>
	!! @brief Constructor
	!!
	subroutine initEmbedding( this, sizeTargetAtoms, targetAtoms, command, energyFilter, gradientFilter, geometryConverter, energyReference )
		class(EmbeddingInterface) :: this 
		integer :: sizeTargetAtoms
		integer :: targetAtoms(*)
		character(1000) :: command
		character(1000) :: energyFilter
		character(1000) :: gradientFilter
		character(1000) :: geometryConverter
		real(8) :: energyReference
		
		integer :: i, pos
		
		allocate( this.targetAtoms(sizeTargetAtoms) )
		allocate( this.gradient(3,sizeTargetAtoms) )
		this.targetAtoms = targetAtoms(1:sizeTargetAtoms)
		this.command = command
		this.energyFilter = energyFilter
		this.gradientFilter = gradientFilter
		this.geometryConverter = geometryConverter
		this.energyReference = energyReference
		
! 		call this.update( debug=.true. )
	end subroutine initEmbedding
	
	!>
	!! @brief Copy constructor
	!!
	subroutine copyEmbedding( this, other )
		class(EmbeddingInterface), intent(inout) :: this
		class(EmbeddingInterface), intent(in) :: other

		stop "### ERROR ### EmbeddingInterface.copy is not implemented yet"
! 		this.targetAtoms(:)
! 		this.command
! 		this.gradientFilter
	end subroutine copyEmbedding
	
	!>
	!! @brief Destructor
	!!
	subroutine destroyEmbedding( this )
		type(EmbeddingInterface) :: this
		
		if( allocated(this.targetAtoms) ) deallocate( this.targetAtoms )
	end subroutine destroyEmbedding
	
	!>
	!! @brief
	!!
	subroutine update( this, symbols, geometry, debug )
		class(EmbeddingInterface) :: this
		character(2) :: symbols(:)
		real(8) :: geometry(:,:)
		logical, optional :: debug
		
		integer :: i, pos
		character(5) :: sBuffer
		character(1000) :: oFileName
		
		oFileName = "interface.emb"
		
		open( 31, file="geom.emb", status="unknown" )
		do i=1,size(geometry(1,:))
			write(31,"(A5,3F15.5)") symbols(i), geometry(:,i)/angs
		end do
		close(31)
		
		call system( trim(this.geometryConverter) )
		
		call system( trim(this.command) )
		
		call system( trim(this.energyFilter)//" > "//trim(oFileName) )
		
		! Loading energy
		!open( 31, file=trim(oFileName), status="old" )
		open( 31, file=trim(oFileName), status="unknown" )
		read( 31, * ) this.energy
		close(31)
		
		this.energy = this.energy - this.energyReference
		
		call system( trim(this.gradientFilter)//" > "//trim(oFileName) )
		
		! Loading gradient
		open( 31, file=trim(oFileName), status="unknown" )
		!open( 31, file=trim(oFileName), status="old" )
		do i=1,size(this.targetAtoms)
			!read( 31, * ) pos, this.gradient(:,pos)
			read( 31, * ) sBuffer, this.gradient(:,i)   !<<< OJO TIENEN QUE SE CONSECUTIVOS. HAY QUE ARREGLAR ESTO
		end do
		close(31)
		
		call system( "rm "//trim(oFileName) )
		
		if( present(debug) .and. debug ) then
			write(*,"(A,F10.5)") "   ENERGY >>", this.energy
		
			do i=1,size(this.targetAtoms)
				write(*,"(A,I5,A5,3F10.5)") "   GRADIENT >>", i, symbols(i), this.gradient(:,i)
			end do
		end if
		
	end subroutine update
	
end module EmbeddingInterface_
