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

program main
	use UnitsConverter_
	use Atom_
	use Molecule_
	implicit none
	
	integer :: nAtoms
	type(Atom) :: atomBuffer
	type(Molecule) :: mol0, mol
	real(8) :: rmax
	integer :: i
! 	type(RandomSampler) :: rs
	
	nAtoms = 20
	rmax = 10.0_8
	
	call mol0.init( "benceno.xyz" )
! 	call mol0.show( formatted=.true. )
	
	call mol.init( nAtoms, "Geometry generated with HeDropletInitialGuess" )
	
	do i=1,nAtoms
		if( i <= mol0.nAtoms() ) then
			atomBuffer = mol0.atoms(i)
		else
			call atomBuffer.init( "He", 0.0_8, 0.0_8, 0.0_8 )
		end if
		
		mol.atoms(i) = atomBuffer
	end do
	
! 	call mol.randomGeometry( radius=1000.0_8*angs, maxIter=100000, fixedAtoms=[1,2,3,4,5,6,7,8,9,10] )
! 	call mol.randomGeometry( maxIter=100000, overlappingRadius=0.2_8*angs, fixedAtoms=[1,2,3,4,5,6,7,8,9,10] )
	call mol.show( formatted=.true. )
	
! 	allocate( this.symbols(nAtoms) )
! 	allocate( this.positions(3,nAtoms) )
! 	
! 	this.symbols = "He"
! 	
! 	call rs.init( nDim=3 )
! 	call rs.setRange( 1, [-rmax,rmax] )
! 	call rs.setRange( 2, [-rmax,rmax] )
! 	call rs.setRange( 3, [-rmax,rmax] )
! 	call rs.buildSample( this.positions, dropletDist, delta=20.0_8, rMin=3.0_8 )
! 	
! 	this.positions(:,1) = 0.0_8
! 	this.positions = this.positions*angs
	
	contains
	
	function dropletDist( r ) result( output )
		real(8), intent(in) :: r(:)
		real(8) :: output
		
! 		output = dropletPtr.density( r )
		output = 0.0_8
	end function dropletDist
	
end program main
