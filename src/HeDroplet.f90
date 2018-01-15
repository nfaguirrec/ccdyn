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

!*
! @brief
!*
module HeDroplet_
	implicit none
	private
	
	public :: &
		HeDroplet_test
	
	type, public :: HeDroplet
		real(8) :: xmax, ymax, zmax
		real(8) :: hx, hy, hz
		complex(8) , allocatable :: psi(:,:,:)
		real(8) , allocatable :: den(:,:,:)
		real(8) , allocatable :: x(:)
		real(8) , allocatable :: y(:)
		real(8) , allocatable :: z(:)
		
		integer(4) :: nx,ny,nz
		
		contains
			procedure :: init
			procedure :: destroy
			generic :: density => densityR, densityXYZ
			
			procedure, private :: load
			procedure, private :: densityR
			procedure, private :: densityXYZ
	end type HeDroplet
	
	contains
	
	!*
	! @brief Constructor
	!*
	subroutine init( this, fileName )
		implicit none
		class(HeDroplet) :: this 
		character(*), intent(in) :: fileName
		
		call this.load( fileName )
	end subroutine init
	
	!*
	! @brief Copy constructor
	!*
	subroutine copy( this, other )
		implicit none
		class(HeDroplet) :: this
		type(HeDroplet), intent(in) :: other

	end subroutine copy
	
	!*
	! @brief Destructor
	!*
	subroutine destroy( this )
		implicit none
		class(HeDroplet) :: this
		
		deallocate(this.psi)
		deallocate(this.den)
		deallocate(this.x)
		deallocate(this.y)
		deallocate(this.z)
		
	end subroutine destroy
	
	!*
	! @brief
	!*
	subroutine load( this, fileName )
		class(HeDroplet) :: this
		character(*) :: fileName
		
		real(8) :: xmax,ymax,zmax,hx,hy,hz,ximp,yimp,zimp
		integer(4) :: ix,iy,iz
		integer(4) :: nx,ny,nz
		logical :: limp
		
		open( unit=1, file=fileName, status='old' )
		read(1,*) xmax,ymax,zmax,hx,hy,hz,nx,ny,nz,limp,ximp,yimp,zimp
		
		this.xmax = xmax
		this.ymax = ymax
		this.zmax = zmax
		
		this.hx = hx
		this.hy = hy
		this.hz = hz
		
		allocate(this.psi(nx,ny,nz))
		allocate(this.den(nx,ny,nz))
		allocate(this.x(nx))
		allocate(this.y(ny))
		allocate(this.z(nz))
		
		! Reads the wave function
		read(1,*) this.psi
		close(1)
		
		! Build grid in real space
		do ix=1,nx
			this.x(ix) = -xmax+hx*(ix-1)
		end do
		
		do iy=1,ny
			this.y(iy) = -ymax+hy*(iy-1)
		end do
		
		do iz=1,nz
			this.z(iz) = -zmax+hz*(iz-1)
		end do
		
		this.den = this.psi*conjg(this.psi)
	end subroutine load
	
	function densityR( this, r ) result( output )
		class(HeDroplet) :: this
		real(8), intent(in) :: r(:)
		real(8) :: output
		
		integer(4) :: ix,iy,iz
		
		ix = int((r(1)+this.xmax)/this.hx)+1
		iy = int((r(2)+this.ymax)/this.hy)+1
		iz = int((r(3)+this.zmax)/this.hz)+1
		
		output = this.den(ix,iy,iz)
	end function densityR
	
	function densityXYZ( this, x, y, z ) result( output )
		class(HeDroplet) :: this
		real(8), intent(in) :: x
		real(8), intent(in) :: y
		real(8), intent(in) :: z
		real(8) :: output
		
		integer(4) :: ix,iy,iz
		
		ix = int((x+this.xmax)/this.hx)+1
		iy = int((y+this.ymax)/this.hy)+1
		iz = int((z+this.zmax)/this.hz)+1
		
		output = this.den(ix,iy,iz)
	end function densityXYZ
	
	!*
	! @brief Test method
	!*
	subroutine HeDroplet_test()
		implicit none
		
		type(HeDroplet) :: droplet
		real(8) :: maxDen
		real(8) :: z
		
! 		call droplet.init( "density.001.dat" )
		call droplet.init( "/home/nestor/Develop/drop-WorkSpace/N300/density.N300.dat" )
! 		call droplet.init( "/home/nestor/Develop/drop-WorkSpace/N600/density.N600.dat" )
! 		call droplet.init( "/home/nestor/Develop/drop-WorkSpace/N1000/density.N1000.dat" )
		
		maxDen = maxval(droplet.den)
		
		z=-25.0_8
		do while( z<=25.0_8 )
! 			write(*,"(2F20.7)") z, droplet.density(0.0_8, 0.0_8, z)/0.021836  ! bulk density from Phys. Rev. B. 52 (1995) 1193
			write(*,"(2F20.7)") z, droplet.density(0.0_8, -z, 0.0_8)/maxDen  ! bulk density from Phys. Rev. B. 52 (1995) 1193
! 			write(*,"(2F20.7)") z, droplet.density(0.0_8, 0.0_8, z)/maxDen  ! bulk density from Phys. Rev. B. 52 (1995) 1193
			z = z + 0.75_8
		end do
		
	end subroutine HeDroplet_test
	
end module HeDroplet_
