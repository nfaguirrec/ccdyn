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

module HeTiO2Potential_
	implicit none
	private
	
	public :: &
		HeTiO2Potential_test
	
	real(8), private, parameter :: pi = acos(-1.0_8)
	real(8), private, parameter :: aCell = 2.9850_8
	real(8), private, parameter :: bCell = 6.49548_8
	
	integer, public, parameter :: M3DP = 1
	integer, public, parameter :: LAP = 2
	integer, public, parameter :: CMP = 3
	integer, public, parameter :: LAPC3LJ = 4
	integer, public, parameter :: LAPC3LJv2 = 5
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Morse 3D potential (M3DP) parameters
	real(8), private :: D_M3DP(3,4)
	real(8), private :: alpha_M3DP(3,4)
	real(8), private :: Re_M3DP(3,4)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Laterally Averaged Potential (LAP) parameters
	real(8), private :: De_LAP
	real(8), private :: alpha_LAP
	real(8), private :: Re_LAP
	
	! Corrección C3LJ a LAP
	real(8), private :: c3LJ_LAP
	real(8), private :: z0LJ_LAP
	real(8), private :: wLJ_LAP
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! C3 Lennard-Jones correction to LAP (LAPC3LJ)
	real(8), private :: De_LAPC3LJ
	real(8), private :: alpha_LAPC3LJ
	real(8), private :: Re_LAPC3LJ
	real(8), private :: c3_LAPC3LJ
	real(8), private :: z0_LAPC3LJ
	real(8), private :: w_LAPC3LJ
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Corrugated Morse Potential (CMP) parameters
	real(8), private :: De_CMP
	real(8), private :: alpha_CMP
	real(8), private :: Re_CMP
	real(8), private :: nu0_CMP
	real(8), private :: zeta_CMP(3,4)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! MorseC3 LAP (LAPC3LJv2)
	real(8), private :: De_LAPC3LJv2
	real(8), private :: alpha_LAPC3LJv2
	real(8), private :: Re_LAPC3LJv2
	real(8), private :: c3_LAPC3LJv2
	real(8), private :: z0_LAPC3LJv2
	real(8), private :: w_LAPC3LJv2
	
	!*
	! @brief Clase que representa el potencial de interaccion HeTiO2(110)
	!*
	type, public :: HeTiO2Potential
		integer, private :: model
		logical, private :: grimme
		
		contains
			procedure :: init
			procedure :: destroy
			procedure :: V
			procedure :: NdV
			procedure :: dV
	end type HeTiO2Potential
	
	contains
	
	!*
	! @brief Contructor
	!*
	subroutine init( this, model, grimme )
		class(HeTiO2Potential) :: this
		integer, optional, intent(in) :: model
		logical, optional, intent(in) :: grimme
		
		integer :: modelEff
		integer :: grimmeEff
		
		if( present(model) ) then
			modelEff = model
		else
			modelEff = M3DP
		end if
		
		if( present(grimme) ) then
			grimmeEff = grimme
		else
			grimmeEff = .false.
		end if
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Morse 3D potential (M3DP) parameters
		D_M3DP = reshape( &
			[ 54.671500_8, -2.228480_8, 0.585194_8, &
			  13.211000_8,  4.842040_8, 0.363525_8, &
			   4.614530_8,  3.747220_8, 0.333350_8, &
			  -1.182030_8,  3.222690_8, 0.587320_8  ], &
			[3,4] )
		
		alpha_M3DP = reshape( &
			[  1.492420_8, -0.010921_8, 0.003189_8, &
			  -0.168927_8, -0.097618_8,-0.024780_8, &
			   0.062441_8, -0.034535_8, 0.012846_8, &
			   0.016801_8, 0.030015_8, -0.011699_8  ], &
			[3,4] )
		
		Re_M3DP = reshape( &
			[ 4.073990_8,  0.057537_8, -0.001500_8, &
			 -0.568990_8, -0.115578_8, -0.018127_8, &
			 -0.090046_8, -0.075153_8, -0.023349_8, &
			  0.008590_8, -0.051210_8, -0.017931_8  ], &
			[3,4] )
			
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Laterally Averaged Potential (LAP) parameters
		De_LAP = 39.95_8
		alpha_LAP = 1.721_8
		Re_LAP = 4.3619_8
		
		! Corrección LJ a LAP
! 		De_LAP = 39.95_8*0.95_8
		c3LJ_LAP = 12.562529_8
		z0LJ_LAP = Re_LAP+2.0_8*log(2.0_8)/alpha_LAP
		wLJ_LAP = 4.0_8
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! C3 Lennard-Jones correction to LAP (LAPC3LJ)
		De_LAPC3LJ = 39.95_8*0.95_8
		alpha_LAPC3LJ = alpha_LAP
		Re_LAPC3LJ = Re_LAP
		c3_LAPC3LJ = 12.562529_8
		z0_LAPC3LJ = Re_LAP+2.0_8*log(2.0_8)/alpha_LAP
		w_LAPC3LJ = 4.0_8
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Corrugated Morse Potential (CMP) parameters
		De_CMP = De_LAP
		alpha_CMP = alpha_LAP
		Re_CMP = Re_LAP
		zeta_CMP = reshape( &
			[ 3.272747_8,  0.043029_8, -0.018600_8, &
			 -0.619662_8, -0.146484_8, -0.012232_8, &
			 -0.046390_8, -0.081815_8, -0.009533_8, &
			  0.025380_8, -0.011625_8, -0.009877_8  ], &
			[3,4] )
		nu0_CMP = 2.5036300_8*exp(2.0_8*alpha_CMP*zeta_CMP(1,1))
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! MorseC3 LAP (LAPC3LJv2)
		De_LAPC3LJv2 = 74.0325066877981_8
		alpha_LAPC3LJv2 = 1.49864795196324_8
		Re_LAPC3LJv2 = 4.32863209306751_8
		c3_LAPC3LJv2 = -15.4243338712553_8
		z0_LAPC3LJv2 = 4.92457636012702_8
		w_LAPC3LJv2 = 5.09690288557883_8
		
		this.model = modelEff
		this.grimme = grimmeEff
	end subroutine init
	
	!*
	! @brief Destructor
	!*
	subroutine destroy( this )
		class(HeTiO2Potential), intent(in) :: this
		
	end subroutine destroy
	
	!*
	! @brief Funcion generica que se encarga de llamar a la funcion
	!        de potencial correcto, segun el atributo model
	!
	! @param x
	!        Coordenada x en angstroms
	! @param y
	!        Coordenada y en angstroms
	! @param z
	!        Coordenada z en angstroms
	!*
	function V( this, x, y, z ) result( output )
		class(HeTiO2Potential), intent(in) :: this
		real(8), intent(in) :: x
		real(8), intent(in) :: y
		real(8), intent(in) :: z
		real(8) :: output
		
		select case( this.model )
			case( M3DP )
				output = V_M3DP( x, y, z )
				
				if( this.grimme ) then
					output = output + EdispGrimme( x, y, z )
				end if
			case( LAP )
				output = V_LAP( x, y, z )
				
				if( this.grimme ) then
					output = output + EdispGrimmeLA( z )
				end if
			case( CMP )
				output = V_CMP( x, y, z )
				
				if( this.grimme ) then
					output = output + EdispGrimme( x, y, z )
				end if
			case( LAPC3LJ )
				output = V_LAPC3LJ( x, y, z )
				
				if( this.grimme ) then
					write(*,*) "## WARNING ## Grimme correction is not available for LAPC3LJ model"
				end if
			case( LAPC3LJv2 )
				output = V_LAPC3LJv2( x, y, z )
				
				if( this.grimme ) then
					write(*,*) "## WARNING ## Grimme correction is not available for LAPC3LJ model"
				end if
		end select
	end function V
	
	
	function NdV( this, x, y, z, coord, stepSize, nPoints ) result( output )
		class(HeTiO2Potential), intent(in) :: this
		real(8), intent(in) :: x
		real(8), intent(in) :: y
		real(8), intent(in) :: z
		integer, intent(in) :: coord
		real(8), optional, intent(in) :: stepSize
		integer, optional, intent(in) :: nPoints
		real(8) :: output
		
		integer :: nPointsEff
		real(8) :: h
		
		nPointsEff = 5
		if( present(nPoints) ) nPointsEff = nPoints
		
		h = 0.0001_8
		if( present(stepSize) ) h = stepSize

#define fx(i) this.V(x+i##.0_8*h,y,z)
#define fy(i) this.V(x,y+i##.0_8*h,z)
#define fz(i) this.V(x,y,z+i##.0_8*h)
		output = 0.0_8
		select case( nPointsEff )
			case( 3 )
				select case( coord )
					case( 1 )
						output = ( fx(1)-fx(-1) )/(2.0_8*h)
					case( 2 )
						output = ( fy(1)-fy(-1) )/(2.0_8*h)
					case( 3 )
						output = ( fz(1)-fz(-1) )/(2.0_8*h)
				end select
			case( 5 )
				select case( coord )
					case( 1 )
						output = ( fx(-2)-8.0_8*fx(-1)+8.0_8*fx(1)-fx(2) )/(12.0_8*h)
					case( 2 )
						output = ( fy(-2)-8.0_8*fy(-1)+8.0_8*fy(1)-fy(2) )/(12.0_8*h)
					case( 3 )
						output = ( fz(-2)-8.0_8*fz(-1)+8.0_8*fz(1)-fz(2) )/(12.0_8*h)
				end select
! 			case( 7 )
! 				output = ( -fz(-3)+9.0_8*fz(-2)-45.0_8*fz(-1)+45.0_8*fz(1)-9.0_8*fz(2)+fz(3) )/(60.0_8*h)
! 			case( 9 )
! 				output = ( 3.0_8*fz(-4)-32.0_8*fz(-3)+168.0_8*fz(-2) &
! 					  -672.0_8*fz(-1)+672.0_8*fz(1)-168.0_8*fz(2)+32.0_8*fz(3)-3.0_8*fz(4) )/(840.0_8*h)
			case default
				write(*,*) "### ERROR ### This formula is not implemented, it's only available nPoints=3,5,7,9"
				stop
		end select
#undef fx
#undef fy
#undef fz
	end function NdV
	
	!*
	! @brief Funcion generica que se encarga de llamar a la funcion
	!        de derivada de potencial correcto, segun el atributo model
	!
	! @param x
	!        Coordenada x en angstroms
	! @param y
	!        Coordenada y en angstroms
	! @param z
	!        Coordenada z en angstroms
	! @param coord
	!        Coordenada respecto a la cual se deriva [ ( 1, 2, 3 ) => ( x, y, z ) ]
	!*
	function dV( this, x, y, z, coord ) result( output )
		class(HeTiO2Potential), intent(in) :: this
		real(8), intent(in) :: x
		real(8), intent(in) :: y
		real(8), intent(in) :: z
		integer, intent(in) :: coord
		real(8) :: output
		
		select case( this.model )
			case( LAP )
				output = dV_LAP( x, y, z, coord )
				
				if( this.grimme ) then
! 					output = output + EdispGrimmeLA( z )
					write(*,*) "## WARNING ## Derivative for Grimme correction is not implemented"
				end if
			case( CMP )
				output = dV_CMP( x, y, z, coord )
				
				if( this.grimme ) then
! 					output = output + EdispGrimme( x, y, z )
					write(*,*) "## WARNING ## Derivative for Grimme correction is not implemented"
				end if
			case( M3DP )
				output = dV_M3DP( x, y, z, coord )
				
				if( this.grimme ) then
! 					output = output + EdispGrimme( x, y, z )
					write(*,*) "## WARNING ## Derivative for Grimme correction is not implemented"
				end if
			case( LAPC3LJ )
				output = dV_LAPC3LJ( x, y, z, coord )
				
				if( this.grimme ) then
! 					output = output + EdispGrimmeLA( z )
					write(*,*) "## WARNING ## Derivative for Grimme correction is not available"
				end if
			case( LAPC3LJv2 )
				output = dV_LAPC3LJv2( x, y, z, coord )
				
				if( this.grimme ) then
					write(*,*) "## WARNING ## Derivative for Grimme correction is not available"
				end if
		end select
	end function dV
	
	!*
	! @brief Potencial segun el modelo Morse3D
	!*
	function V_M3DP( x, y, z ) result( output )
		real(8), intent(in) :: x
		real(8), intent(in) :: y
		real(8), intent(in) :: z
		real(8) :: output
		
#define Dxy F(x,y,1)
#define Axy F(x,y,2)
#define Zxy F(x,y,3)
		
		output = Dxy*( exp(-2.0_8*Axy*( z-Zxy ))-2.0_8*exp(-Axy*(z-Zxy)) )
		
#undef Vxyz
#undef Dxy
#undef Axy
	end function V_M3DP
	
	!*
	! @brief Derivada del potencial segun el modelo Morse3D
	!*
	function dV_M3DP( x, y, z, coord ) result( output )
		real(8), intent(in) :: x
		real(8), intent(in) :: y
		real(8), intent(in) :: z
		integer, intent(in) :: coord
		real(8) :: output
		
#define Vxyz V_M3DP(x,y,z)
#define Dxy F(x,y,1)
#define Axy F(x,y,2)
#define Zxy F(x,y,3)
#define dDxydx dF(x,y,1,1)
#define dDxydy dF(x,y,1,2)
#define dAxydx dF(x,y,2,1)
#define dAxydy dF(x,y,2,2)
#define dZxydx dF(x,y,3,1)
#define dZxydy dF(x,y,3,2)

		select case ( coord )
			case( 1 )
				output = dDxydx*Vxyz/Dxy &
					-2.0_8*Dxy*(exp(-2.0_8*Axy)-exp(-Axy))*( dAxydx*(z-Zxy)-Axy*dZxydx )
			case( 2 )
				output = dDxydy*Vxyz/Dxy &
					-2.0_8*Dxy*(exp(-2.0_8*Axy)-exp(-Axy))*( dAxydy*(z-Zxy)-Axy*dZxydy )
			case( 3 )
				output = 2.0_8*Dxy*Axy*( exp(-Axy*(z-Zxy)) - exp(-2.0_8*Axy*(z-Zxy)) )
		end select
	
#undef Vxyz
#undef Dxy
#undef Axy
#undef Zxy
#undef dDxydx
#undef dDxydy
#undef dAxydx
#undef dAxydy
#undef dZxydx
#undef dZxydy
	end function dV_M3DP
	
	!*
	! @brief Potencial segun el modelo de potencial
	!        lateralmente promediado
	!*
	function V_LAP( x, y, z ) result( output )
		real(8), intent(in) :: x
		real(8), intent(in) :: y
		real(8), intent(in) :: z
		real(8) :: output
		
		output = De_LAP*( exp(-2.0_8*alpha_LAP*(z-Re_LAP)) - 2.0_8*exp(-alpha_LAP*(z-Re_LAP)) ) !+ EdispC3LJ_LAP( z )
	end function V_LAP
	
	!*
	! @brief Derivada del potencial segun el modelo de potencial
	!        lateralmente promediado
	!*
	function dV_LAP( x, y, z, coord ) result( output )
		real(8), intent(in) :: x
		real(8), intent(in) :: y
		real(8), intent(in) :: z
		integer, intent(in) :: coord
		real(8) :: output
		
		select case ( coord )
			case( 1 )
				output = 0.0_8
			case( 2 )
				output = 0.0_8
			case( 3 )
				output = 2.0_8*De_LAP*alpha_LAP*( exp(-alpha_LAP*(z-Re_LAP)) - exp(-2.0_8*alpha_LAP*(z-Re_LAP)) ) !+ dEdispC3LJ_LAP( z )
		end select
	end function dV_LAP
	
	!*
	! @brief Potencial segun el modelo de potencial
	!        de Morse corrugado
	!*
	function V_CMP( x, y, z ) result( output )
		real(8), intent(in) :: x
		real(8), intent(in) :: y
		real(8), intent(in) :: z
		real(8) :: output
		
		real(8) :: preexp
		
#define zetaxy F(x,y,4)
		
		preexp = exp(2.0_8*alpha_CMP*zetaxy)/nu0_CMP
		output = De_CMP*( preexp*exp(-2.0_8*alpha_CMP*(z-Re_CMP)) - 2.0_8*exp(-alpha_CMP*(z-Re_CMP)) )
		
#undef Dxy
	end function V_CMP
	
	!*
	! @brief Derivada del potencial segun el modelo de potencial
	!        de Morse corrugado
	!*
	function dV_CMP( x, y, z, coord ) result( output )
		real(8), intent(in) :: x
		real(8), intent(in) :: y
		real(8), intent(in) :: z
		integer, intent(in) :: coord
		real(8) :: output
		
#define zetaxy F(x,y,4)
#define dzetaxydx dF(x,y,4,1)
#define dzetaxydy dF(x,y,4,2)
			
		select case ( coord )
			case( 1 )
				output = 2.0_8*De_CMP*exp(-2.0_8*alpha_CMP*(z-Re_CMP-zetaxy))*alpha_CMP*dzetaxydx
			case( 2 )
				output = 2.0_8*De_CMP*exp(-2.0_8*alpha_CMP*(z-Re_CMP-zetaxy))*alpha_CMP*dzetaxydy
			case( 3 )
				output = 2.0_8*De_CMP*alpha_CMP*( exp(-alpha_CMP*(z-Re_CMP)) - exp(-2.0_8*alpha_CMP*(z-Re_CMP-zetaxy)) )
		end select
			
#undef zetaxy
#undef dzetaxydx
#undef dzetaxydy
	end function dV_CMP
	
	!*
	! @brief Funcion auxiliar generica
	! param = [ ( 1, 2, 3, 4 ) => ( D_M3DP, alpha_M3DP, Re_M3DP, zeta_CMP ) ]
	!*
	pure function F( x, y, param ) result( output )
		real(8), intent(in) :: x
		real(8), intent(in) :: y
		integer, intent(in) :: param
		real(8) :: output
		
		real(8) :: Lx(3,1)
		real(8) :: Ly(4,1)
		real(8) :: rbuffer(1,1)
		
		Lx = transpose(reshape( &
			cos( 2.0_8*pi*[ 0.0_8, 1.0_8, 2.0_8 ]*x/aCell ), &
		     [1,3] ))
		
		Ly = transpose(reshape( &
			cos( 2.0_8*pi*[ 0.0_8, 1.0_8, 2.0_8, 3.0_8 ]*y/bCell ), &
		     [1,4] ))
		
		select case ( param )
			case( 1 )
				rbuffer = matmul( matmul( transpose(Lx), D_M3DP ), Ly )
			case( 2 )
				rbuffer = matmul( matmul( transpose(Lx), alpha_M3DP ), Ly )
			case( 3 )
				rbuffer = matmul( matmul( transpose(Lx), Re_M3DP ), Ly )
			case( 4 )
				rbuffer = matmul( matmul( transpose(Lx), zeta_CMP ), Ly )
		end select
		
		output = rbuffer(1,1)
	end function F
	
	!*
	! @brief Derivada de la funcion auxiliar generica
	! param = [ ( 1, 2, 3, 4 ) => ( D_M3DP, alpha_M3DP, Re_M3DP, zeta_CMP ) ]
	! coord = [ ( 1, 2 ) => ( x, y ) ]
	!*
	pure function dF( x, y, param, coord ) result( output )
		real(8), intent(in) :: x
		real(8), intent(in) :: y
		integer, intent(in) :: param
		integer, intent(in) :: coord
		real(8) :: output
		
		real(8) :: Lx(3,1)
		real(8) :: Ly(4,1)
		real(8) :: Sx(3,1)
		real(8) :: Sy(4,1)
		
		real(8) :: Work(3,4)
		real(8) :: rbuffer(1,1)
		
		select case ( param )
			case( 1 )
				Work = D_M3DP
			case( 2 )
				Work = alpha_M3DP
			case( 3 )
				Work = Re_M3DP
			case( 4 )
				Work = zeta_CMP
		end select
		
		select case ( coord )
			case( 1 )
				Sx = transpose(reshape( &
					[0.0_8, 1.0_8, 2.0_8]*sin( 2.0_8*pi*[ 0.0_8, 1.0_8, 2.0_8 ]*x/aCell ), &
					[1,3] ))
					
				Ly = transpose(reshape( &
					cos( 2.0_8*pi*[ 0.0_8, 1.0_8, 2.0_8, 3.0_8 ]*y/bCell ), &
					[1,4] ))
				
				rbuffer = matmul( matmul( transpose(Sx), Work ), Ly )
				output = -(2.0_8*pi/aCell)*rbuffer(1,1)
			case( 2 )
				Lx = transpose(reshape( &
					cos( 2.0_8*pi*[ 0.0_8, 1.0_8, 2.0_8 ]*x/aCell ), &
					[1,3] ))
				
				Sy = transpose(reshape( &
					[ 0.0_8, 1.0_8, 2.0_8, 3.0_8 ]*sin( 2.0_8*pi*[ 0.0_8, 1.0_8, 2.0_8, 3.0_8 ]*y/bCell ), &
					[1,4] ))
				
				rbuffer = matmul( matmul( transpose(Lx), Work ), Sy )
				output = -(2.0_8*pi/bCell)*rbuffer(1,1)
		end select
	end function dF
	
	!*
	! @brief Potencial segun el modelo de potencial lateralmente
	!        promediado corregido por c3 de Lennard-Jones
	!*
	function V_LAPC3LJ( x, y, z ) result( output )
		real(8), intent(in) :: x
		real(8), intent(in) :: y
		real(8), intent(in) :: z
		real(8) :: output
		
		output = De_LAPC3LJ*( exp(-2.0_8*alpha_LAPC3LJ*(z-Re_LAPC3LJ)) &
			- 2.0_8*exp(-alpha_LAPC3LJ*(z-Re_LAPC3LJ)) ) + EdispC3LJ( z, z0_LAPC3LJ, w_LAPC3LJ, c3_LAPC3LJ )
			
	end function V_LAPC3LJ
	
	!*
	! @brief Derivada del potencial segun el modelo de potencial
	!        lateralmente promediado corregido por c3 de Lennard-Jones
	!*
	function dV_LAPC3LJ( x, y, z, coord ) result( output )
		real(8), intent(in) :: x
		real(8), intent(in) :: y
		real(8), intent(in) :: z
		integer, intent(in) :: coord
		real(8) :: output
		
		select case ( coord )
			case( 1 )
				output = 0.0_8
			case( 2 )
				output = 0.0_8
			case( 3 )
				output = 2.0_8*De_LAPC3LJ*alpha_LAPC3LJ*( exp(-alpha_LAPC3LJ*(z-Re_LAPC3LJ)) &
							- exp(-2.0_8*alpha_LAPC3LJ*(z-Re_LAPC3LJ)) ) + dEdispC3LJ( z, z0_LAPC3LJ, w_LAPC3LJ, c3_LAPC3LJ )
		end select
	end function dV_LAPC3LJ
	
	!*
	! @brief Funcion de damping ( función tipo sigmoide )
	!
	! @param z
	!        Coordenada independiente
	! @param z0
	!        Parámetro que indica donde la forma sigmoide de la función
	!        vale justo 0.5
	! @param w
	!        Parámetro que indica el ancho de la sigmoide
	! @param der
	!        Opcional: si es .true. retorna el valor de la primera
	!        derivada respecto a z
	!*
	real(8) function fDamping( z, z0, w, der )
		real(8), intent(in) :: z
		real(8), intent(in) :: z0
		real(8), intent(in) :: w
		logical, optional, intent(in) :: der
		
		integer :: nout
		
! 		if( present(der) .and. der == .true. ) then
! 			fDamping = 0.5_8*( 1.0_8 + tanh(6.0_8*(z-z0)/w) )
! 		else
! 			! to avoid floating overflow error in cosh, cosh(100.0) = 1.344058570908068E+043
! 			if( abs(6.0_8*(z-z0)/w) < 100.0_8 ) then
! 				fDamping = (3.0_8/w)*( 1.0_8/cosh(6.0_8*(z-z0)/w) )**2
! 			else
! 				fDamping = 0.0_8
! 			end if
! 		end if

		if( present(der) .and. der == .true. ) then
			! to avoid floating overflow error in cosh, cosh(100.0) = 1.344058570908068E+043
			if( abs(6.0_8*(z-z0)/w) < 100.0_8 ) then
				fDamping = (3.0_8/w)*( 1.0_8/cosh(6.0_8*(z-z0)/w) )**2
			else
				fDamping = 0.0_8
			end if
		else
			fDamping = 0.5_8*( 1.0_8 + tanh(6.0_8*(z-z0)/w) )
		end if

	end function fDamping
	
	!*
	! @brief Energía de dispersión promediada lateralmente proveniente de suponer que
	!        cada átomo de la superficie interacciona con el átomo de He a través de
	!        un potencial de pares tipo Lennard-Jones (12,6)
	!
	! @param zHe
	!        Posición de átomo de He relativo al Ti(5f) de la superficie
	!*
	real(8) function EdispC3LJ( zHe, z0, w, c3 )
		real(8), intent(in) :: zHe
		real(8), intent(in) :: z0
		real(8), intent(in) :: w
		real(8), intent(in) :: c3
		
		EdispC3LJ = -fDamping( zHe, z0, w )*( c3/zHe )**3
	end function EdispC3LJ
	
	!*
	! @brief Derivada respecto a z de la energía de dispersión tipo Lennard-Jones
	!        promediada lateralmente
	!
	! @param zHe
	!        Posición de átomo de He relativo al Ti(5f) de la superficie
	!*
	real(8) function dEdispC3LJ( zHe, z0, w, c3 )
		real(8), intent(in) :: zHe
		real(8), intent(in) :: z0
		real(8), intent(in) :: w
		real(8), intent(in) :: c3
		
		dEdispC3LJ = ( 3.0_8*fDamping( zHe, z0, w )/zHe - fDamping(zHe, z0, w, der=.true.) )*( c3/zHe )**3
	end function dEdispC3LJ
	
	!*
	! @brief Promedio lateral de la energía de dispersión de Grimme
	!
	! @param zHe
	!        Posición de átomo de He relativo al Ti(5f) de la superficie
	!*
	real(8) function EdispGrimmeLA( zHe )
		real(8), intent(in) :: zHe
		
		integer, parameter :: nPointsX = 10
		integer, parameter :: nPointsY = 5
		
		integer :: i, j
		real(8) :: xi, yj
		real(8) :: dx, dy
		real(8) :: sum
		
		dx = aCell/nPointsX
		dy = bCell/nPointsY
		
		sum = 0.0_8
		do i=1,nPointsX
			xi = real(i,8)*dx
			do j=1,nPointsY
				yj = real(j,8)*dy
				sum = sum + EdispGrimme( xi, yj, zHe )
			end do
		end do
		
! 		write(*,*) xi, yj, nPointsX*dx, nPointsY*dy
		
		EdispGrimmeLA = sum/nPointsX/nPointsY
	end function EdispGrimmeLA
	
	!*
	! @brief Energía de dispersión de Grimme
	!
	! @param xHe, yHe, zHe
	!        Posición de átomo de He relativo al Ti(5f) de la superficie
	!*
	real(8) function EdispGrimme( xHe, yHe, zHe )
		real(8), intent(in) :: xHe
		real(8), intent(in) :: yHe
		real(8), intent(in) :: zHe
		
		real(8), parameter :: Htocm1 = 219474.63068_8
		real(8), parameter :: unit_J_to_H = 2.2937128d17
		real(8), parameter :: unit_nm_to_a0 = 18.89726134d0
		real(8), parameter :: unit_pm_to_a0 = unit_nm_to_a0/1d3
		real(8), parameter :: unit_mol_to_au = 6.0221417930d23
		real(8), parameter :: unit_Ang_to_a0 = 1d0/0.529177208
		
		integer, parameter :: nAtomsO = 3811
		integer, parameter :: nAtomsTi = 1984
		
		real(8) :: c6O, rvdWO
		real(8) :: c6Ti, rvdWTi
		real(8) :: c6He, rvdWHe
		real(8) :: d, s6_PBE
		
		real(8) :: c6OHe
		real(8) :: c6TiHe
		real(8) :: rvdWOHe
		real(8) :: rvdWTiHe
		real(8) :: dampF
		integer :: i, j
		real(8) :: r
		
		real(8) :: xO(nAtomsO), yO(nAtomsO), zO(nAtomsO)
		real(8) :: xTi(nAtomsTi), yTi(nAtomsTi), zTi(nAtomsTi)
		
		xTi = [ &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		2.98500,   2.98500,   2.98500,   2.98500,   2.98500, &
		2.98500,   2.98500,   2.98500,   2.98500,   2.98500, &
		5.97000,   5.97000,   5.97000,   5.97000,   5.97000, &
		5.97000,   5.97000,   5.97000,   5.97000,   5.97000, &
		8.95500,   8.95500,   8.95500,   8.95500,   8.95500, &
		8.95500,   8.95500,   8.95500,   8.95500,   8.95500, &
		11.94000,  11.94000,  11.94000,  11.94000,  11.94000, &
		11.94000,  11.94000,  11.94000,  11.94000,  11.94000, &
		14.92500,  14.92500,  14.92500,  14.92500,  14.92500, &
		14.92500,  14.92500,  14.92500,  14.92500,  14.92500, &
		17.91000,  17.91000,  17.91000,  17.91000,  17.91000, &
		17.91000,  17.91000,  17.91000,  17.91000,  17.91000, &
		20.89500,  20.89500,  20.89500,  20.89500,  20.89500, &
		20.89500,  20.89500,  20.89500,  20.89500,  20.89500, &
		23.88000,  23.88000,  23.88000,  23.88000,  23.88000, &
		23.88000,  23.88000,  23.88000,  23.88000,  23.88000, &
		26.86500,  26.86500,  26.86500,  26.86500,  26.86500, &
		26.86500,  26.86500,  26.86500,  26.86500,  26.86500, &
		-29.85000, -29.85000, -29.85000, -29.85000, -29.85000, &
		-29.85000, -29.85000, -29.85000, -29.85000, -29.85000, &
		-26.86500, -26.86500, -26.86500, -26.86500, -26.86500, &
		-26.86500, -26.86500, -26.86500, -26.86500, -26.86500, &
		-23.88000, -23.88000, -23.88000, -23.88000, -23.88000, &
		-23.88000, -23.88000, -23.88000, -23.88000, -23.88000, &
		-20.89500, -20.89500, -20.89500, -20.89500, -20.89500, &
		-20.89500, -20.89500, -20.89500, -20.89500, -20.89500, &
		-17.91000, -17.91000, -17.91000, -17.91000, -17.91000, &
		-17.91000, -17.91000, -17.91000, -17.91000, -17.91000, &
		-14.92500, -14.92500, -14.92500, -14.92500, -14.92500, &
		-14.92500, -14.92500, -14.92500, -14.92500, -14.92500, &
		-11.94000, -11.94000, -11.94000, -11.94000, -11.94000, &
		-11.94000, -11.94000, -11.94000, -11.94000, -11.94000, &
		-8.95500,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-8.95500,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-5.97000,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-5.97000,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-2.98500,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		-2.98500,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		-1.49250,  -1.49250,  -1.49250,  -1.49250,  -1.49250, &
		-1.49250,  -1.49250,  -1.49250,  -1.49250,  -1.49250, &
		1.49250,   1.49250,   1.49250,   1.49250,   1.49250, &
		1.49250,   1.49250,   1.49250,   1.49250,   1.49250, &
		4.47750,   4.47750,   4.47750,   4.47750,   4.47750, &
		4.47750,   4.47750,   4.47750,   4.47750,   4.47750, &
		7.46250,   7.46250,   7.46250,   7.46250,   7.46250, &
		7.46250,   7.46250,   7.46250,   7.46250,   7.46250, &
		10.44750,  10.44750,  10.44750,  10.44750,  10.44750, &
		10.44750,  10.44750,  10.44750,  10.44750,  10.44750, &
		13.43250,  13.43250,  13.43250,  13.43250,  13.43250, &
		13.43250,  13.43250,  13.43250,  13.43250,  13.43250, &
		16.41750,  16.41750,  16.41750,  16.41750,  16.41750, &
		16.41750,  16.41750,  16.41750,  16.41750,  16.41750, &
		19.40250,  19.40250,  19.40250,  19.40250,  19.40250, &
		19.40250,  19.40250,  19.40250,  19.40250,  19.40250, &
		22.38750,  22.38750,  22.38750,  22.38750,  22.38750, &
		22.38750,  22.38750,  22.38750,  22.38750,  22.38750, &
		25.37250,  25.37250,  25.37250,  25.37250,  25.37250, &
		25.37250,  25.37250,  25.37250,  25.37250,  25.37250, &
		28.35750,  28.35750,  28.35750,  28.35750,  28.35750, &
		28.35750,  28.35750,  28.35750,  28.35750,  28.35750, &
		-28.35750, -28.35750, -28.35750, -28.35750, -28.35750, &
		-28.35750, -28.35750, -28.35750, -28.35750, -28.35750, &
		-25.37250, -25.37250, -25.37250, -25.37250, -25.37250, &
		-25.37250, -25.37250, -25.37250, -25.37250, -25.37250, &
		-22.38750, -22.38750, -22.38750, -22.38750, -22.38750, &
		-22.38750, -22.38750, -22.38750, -22.38750, -22.38750, &
		-19.40250, -19.40250, -19.40250, -19.40250, -19.40250, &
		-19.40250, -19.40250, -19.40250, -19.40250, -19.40250, &
		-16.41750, -16.41750, -16.41750, -16.41750, -16.41750, &
		-16.41750, -16.41750, -16.41750, -16.41750, -16.41750, &
		-13.43250, -13.43250, -13.43250, -13.43250, -13.43250, &
		-13.43250, -13.43250, -13.43250, -13.43250, -13.43250, &
		-10.44750, -10.44750, -10.44750, -10.44750, -10.44750, &
		-10.44750, -10.44750, -10.44750, -10.44750, -10.44750, &
		-7.46250,  -7.46250,  -7.46250,  -7.46250,  -7.46250, &
		-7.46250,  -7.46250,  -7.46250,  -7.46250,  -7.46250, &
		-4.47750,  -4.47750,  -4.47750,  -4.47750,  -4.47750, &
		-4.47750,  -4.47750,  -4.47750,  -4.47750,  -4.47750, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		2.98500,   2.98500,   2.98500,   2.98500,   2.98500, &
		2.98500,   2.98500,   2.98500,   2.98500,   2.98500, &
		5.97000,   5.97000,   5.97000,   5.97000,   5.97000, &
		5.97000,   5.97000,   5.97000,   5.97000,   5.97000, &
		8.95500,   8.95500,   8.95500,   8.95500,   8.95500, &
		8.95500,   8.95500,   8.95500,   8.95500,   8.95500, &
		11.94000,  11.94000,  11.94000,  11.94000,  11.94000, &
		11.94000,  11.94000,  11.94000,  11.94000,  11.94000, &
		14.92500,  14.92500,  14.92500,  14.92500,  14.92500, &
		14.92500,  14.92500,  14.92500,  14.92500,  14.92500, &
		17.91000,  17.91000,  17.91000,  17.91000,  17.91000, &
		17.91000,  17.91000,  17.91000,  17.91000,  17.91000, &
		20.89500,  20.89500,  20.89500,  20.89500,  20.89500, &
		20.89500,  20.89500,  20.89500,  20.89500,  20.89500, &
		23.88000,  23.88000,  23.88000,  23.88000,  23.88000, &
		23.88000,  23.88000,  23.88000,  23.88000,  23.88000, &
		26.86500,  26.86500,  26.86500,  26.86500,  26.86500, &
		26.86500,  26.86500,  26.86500,  26.86500,  26.86500, &
		-29.85000, -29.85000, -29.85000, -29.85000, -29.85000, &
		-29.85000, -29.85000, -29.85000, -29.85000, -29.85000, &
		-26.86500, -26.86500, -26.86500, -26.86500, -26.86500, &
		-26.86500, -26.86500, -26.86500, -26.86500, -26.86500, &
		-23.88000, -23.88000, -23.88000, -23.88000, -23.88000, &
		-23.88000, -23.88000, -23.88000, -23.88000, -23.88000, &
		-20.89500, -20.89500, -20.89500, -20.89500, -20.89500, &
		-20.89500, -20.89500, -20.89500, -20.89500, -20.89500, &
		-17.91000, -17.91000, -17.91000, -17.91000, -17.91000, &
		-17.91000, -17.91000, -17.91000, -17.91000, -17.91000, &
		-14.92500, -14.92500, -14.92500, -14.92500, -14.92500, &
		-14.92500, -14.92500, -14.92500, -14.92500, -14.92500, &
		-11.94000, -11.94000, -11.94000, -11.94000, -11.94000, &
		-11.94000, -11.94000, -11.94000, -11.94000, -11.94000, &
		-8.95500,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-8.95500,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-5.97000,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-5.97000,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-2.98500,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		-2.98500,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		-1.49250,  -1.49250,  -1.49250,  -1.49250,  -1.49250, &
		-1.49250,  -1.49250,  -1.49250,  -1.49250,  -1.49250, &
		1.49250,   1.49250,   1.49250,   1.49250,   1.49250, &
		1.49250,   1.49250,   1.49250,   1.49250,   1.49250, &
		4.47750,   4.47750,   4.47750,   4.47750,   4.47750, &
		4.47750,   4.47750,   4.47750,   4.47750,   4.47750, &
		7.46250,   7.46250,   7.46250,   7.46250,   7.46250, &
		7.46250,   7.46250,   7.46250,   7.46250,   7.46250, &
		10.44750,  10.44750,  10.44750,  10.44750,  10.44750, &
		10.44750,  10.44750,  10.44750,  10.44750,  10.44750, &
		13.43250,  13.43250,  13.43250,  13.43250,  13.43250, &
		13.43250,  13.43250,  13.43250,  13.43250,  13.43250, &
		16.41750,  16.41750,  16.41750,  16.41750,  16.41750, &
		16.41750,  16.41750,  16.41750,  16.41750,  16.41750, &
		19.40250,  19.40250,  19.40250,  19.40250,  19.40250, &
		19.40250,  19.40250,  19.40250,  19.40250,  19.40250, &
		22.38750,  22.38750,  22.38750,  22.38750,  22.38750, &
		22.38750,  22.38750,  22.38750,  22.38750,  22.38750, &
		25.37250,  25.37250,  25.37250,  25.37250,  25.37250, &
		25.37250,  25.37250,  25.37250,  25.37250,  25.37250, &
		28.35750,  28.35750,  28.35750,  28.35750,  28.35750, &
		28.35750,  28.35750,  28.35750,  28.35750,  28.35750, &
		-28.35750, -28.35750, -28.35750, -28.35750, -28.35750, &
		-28.35750, -28.35750, -28.35750, -28.35750, -28.35750, &
		-25.37250, -25.37250, -25.37250, -25.37250, -25.37250, &
		-25.37250, -25.37250, -25.37250, -25.37250, -25.37250, &
		-22.38750, -22.38750, -22.38750, -22.38750, -22.38750, &
		-22.38750, -22.38750, -22.38750, -22.38750, -22.38750, &
		-19.40250, -19.40250, -19.40250, -19.40250, -19.40250, &
		-19.40250, -19.40250, -19.40250, -19.40250, -19.40250, &
		-16.41750, -16.41750, -16.41750, -16.41750, -16.41750, &
		-16.41750, -16.41750, -16.41750, -16.41750, -16.41750, &
		-13.43250, -13.43250, -13.43250, -13.43250, -13.43250, &
		-13.43250, -13.43250, -13.43250, -13.43250, -13.43250, &
		-10.44750, -10.44750, -10.44750, -10.44750, -10.44750, &
		-10.44750, -10.44750, -10.44750, -10.44750, -10.44750, &
		-7.46250,  -7.46250,  -7.46250,  -7.46250,  -7.46250, &
		-7.46250,  -7.46250,  -7.46250,  -7.46250,  -7.46250, &
		-4.47750,  -4.47750,  -4.47750,  -4.47750,  -4.47750, &
		-4.47750,  -4.47750,  -4.47750,  -4.47750,  -4.47750, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		2.98500,   2.98500,   2.98500,   2.98500,   2.98500, &
		2.98500,   2.98500,   2.98500,   2.98500,   2.98500, &
		5.97000,   5.97000,   5.97000,   5.97000,   5.97000, &
		5.97000,   5.97000,   5.97000,   5.97000,   5.97000, &
		8.95500,   8.95500,   8.95500,   8.95500,   8.95500, &
		8.95500,   8.95500,   8.95500,   8.95500,   8.95500, &
		11.94000,  11.94000,  11.94000,  11.94000,  11.94000, &
		11.94000,  11.94000,  11.94000,  11.94000,  11.94000, &
		14.92500,  14.92500,  14.92500,  14.92500,  14.92500, &
		14.92500,  14.92500,  14.92500,  14.92500,  14.92500, &
		17.91000,  17.91000,  17.91000,  17.91000,  17.91000, &
		17.91000,  17.91000,  17.91000,  17.91000,  17.91000, &
		20.89500,  20.89500,  20.89500,  20.89500,  20.89500, &
		20.89500,  20.89500,  20.89500,  20.89500,  20.89500, &
		23.88000,  23.88000,  23.88000,  23.88000,  23.88000, &
		23.88000,  23.88000,  23.88000,  23.88000,  23.88000, &
		26.86500,  26.86500,  26.86500,  26.86500,  26.86500, &
		26.86500,  26.86500,  26.86500,  26.86500,  26.86500, &
		-29.85000, -29.85000, -29.85000, -29.85000, -29.85000, &
		-29.85000, -29.85000, -29.85000, -29.85000, -29.85000, &
		-26.86500, -26.86500, -26.86500, -26.86500, -26.86500, &
		-26.86500, -26.86500, -26.86500, -26.86500, -26.86500, &
		-23.88000, -23.88000, -23.88000, -23.88000, -23.88000, &
		-23.88000, -23.88000, -23.88000, -23.88000, -23.88000, &
		-20.89500, -20.89500, -20.89500, -20.89500, -20.89500, &
		-20.89500, -20.89500, -20.89500, -20.89500, -20.89500, &
		-17.91000, -17.91000, -17.91000, -17.91000, -17.91000, &
		-17.91000, -17.91000, -17.91000, -17.91000, -17.91000, &
		-14.92500, -14.92500, -14.92500, -14.92500, -14.92500, &
		-14.92500, -14.92500, -14.92500, -14.92500, -14.92500, &
		-11.94000, -11.94000, -11.94000, -11.94000, -11.94000, &
		-11.94000, -11.94000, -11.94000, -11.94000, -11.94000, &
		-8.95500,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-8.95500,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-5.97000,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-5.97000,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-2.98500,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		-2.98500,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		-1.49250,  -1.49250,  -1.49250,  -1.49250,  -1.49250, &
		-1.49250,  -1.49250,  -1.49250,  -1.49250,  -1.49250, &
		1.49250,   1.49250,   1.49250,   1.49250,   1.49250, &
		1.49250,   1.49250,   1.49250,   1.49250,   1.49250, &
		4.47750,   4.47750,   4.47750,   4.47750,   4.47750, &
		4.47750,   4.47750,   4.47750,   4.47750,   4.47750, &
		7.46250,   7.46250,   7.46250,   7.46250,   7.46250, &
		7.46250,   7.46250,   7.46250,   7.46250,   7.46250, &
		10.44750,  10.44750,  10.44750,  10.44750,  10.44750, &
		10.44750,  10.44750,  10.44750,  10.44750,  10.44750, &
		13.43250,  13.43250,  13.43250,  13.43250,  13.43250, &
		13.43250,  13.43250,  13.43250,  13.43250,  13.43250, &
		16.41750,  16.41750,  16.41750,  16.41750,  16.41750, &
		16.41750,  16.41750,  16.41750,  16.41750,  16.41750, &
		19.40250,  19.40250,  19.40250,  19.40250,  19.40250, &
		19.40250,  19.40250,  19.40250,  19.40250,  19.40250, &
		22.38750,  22.38750,  22.38750,  22.38750,  22.38750, &
		22.38750,  22.38750,  22.38750,  22.38750,  22.38750, &
		25.37250,  25.37250,  25.37250,  25.37250,  25.37250, &
		25.37250,  25.37250,  25.37250,  25.37250,  25.37250, &
		28.35750,  28.35750,  28.35750,  28.35750,  28.35750, &
		28.35750,  28.35750,  28.35750,  28.35750,  28.35750, &
		-28.35750, -28.35750, -28.35750, -28.35750, -28.35750, &
		-28.35750, -28.35750, -28.35750, -28.35750, -28.35750, &
		-25.37250, -25.37250, -25.37250, -25.37250, -25.37250, &
		-25.37250, -22.38750, -22.38750, -22.38750, -22.38750, &
		-22.38750, -22.38750, -19.40250, -19.40250, -19.40250, &
		-19.40250, -19.40250, -19.40250, -19.40250, -19.40250, &
		-19.40250, -16.41750, -16.41750, -16.41750, -13.43250, &
		-13.43250, -13.43250, -13.43250, -13.43250, -13.43250, &
		-13.43250, -13.43250, -13.43250, -13.43250, -10.44750, &
		-10.44750, -10.44750, -10.44750, -10.44750, -10.44750, &
		-10.44750, -10.44750, -10.44750, -10.44750,  -7.46250, &
		-7.46250,  -7.46250,  -7.46250,  -7.46250,  -7.46250, &
		-7.46250,  -7.46250,  -7.46250,  -7.46250,  -4.47750, &
		-4.47750,  -4.47750,  -4.47750,  -4.47750,  -4.47750, &
		-4.47750,  -4.47750,  -4.47750,  -4.47750,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   2.98500, &
		2.98500,   2.98500,   2.98500,   2.98500,   2.98500, &
		2.98500,   2.98500,   2.98500,   2.98500,   5.97000, &
		5.97000,   5.97000,   5.97000,   5.97000,   5.97000, &
		5.97000,   5.97000,   5.97000,   5.97000,   8.95500, &
		8.95500,   8.95500,   8.95500,   8.95500,   8.95500, &
		8.95500,   8.95500,   8.95500,   8.95500,  11.94000, &
		11.94000,  11.94000,  11.94000,  11.94000,  11.94000, &
		11.94000,  11.94000,  11.94000,  11.94000,  14.92500, &
		14.92500,  14.92500,  14.92500,  14.92500,  14.92500, &
		14.92500,  14.92500,  14.92500,  14.92500,  17.91000, &
		17.91000,  17.91000,  17.91000,  17.91000,  17.91000, &
		17.91000,  17.91000,  17.91000,  17.91000,  20.89500, &
		20.89500,  20.89500,  20.89500,  20.89500,  20.89500, &
		20.89500,  20.89500,  20.89500,  20.89500,  23.88000, &
		23.88000,  23.88000,  23.88000,  23.88000,  23.88000, &
		23.88000,  23.88000,  23.88000,  23.88000,  26.86500, &
		26.86500,  26.86500,  26.86500,  26.86500,  26.86500, &
		26.86500,  26.86500,  26.86500,  26.86500, -29.85000, &
		-29.85000, -29.85000, -29.85000, -29.85000, -29.85000, &
		-29.85000, -29.85000, -29.85000, -29.85000, -26.86500, &
		-26.86500, -26.86500, -26.86500, -26.86500, -26.86500, &
		-26.86500, -26.86500, -26.86500, -26.86500, -23.88000, &
		-23.88000, -23.88000, -23.88000, -23.88000, -23.88000, &
		-23.88000, -23.88000, -23.88000, -23.88000, -20.89500, &
		-20.89500, -20.89500, -20.89500, -20.89500, -20.89500, &
		-20.89500, -20.89500, -20.89500, -20.89500, -17.91000, &
		-17.91000, -17.91000, -17.91000, -17.91000, -17.91000, &
		-17.91000, -17.91000, -17.91000, -17.91000, -14.92500, &
		-14.92500, -14.92500, -14.92500, -14.92500, -14.92500, &
		-14.92500, -14.92500, -14.92500, -14.92500, -11.94000, &
		-11.94000, -11.94000, -11.94000, -11.94000, -11.94000, &
		-11.94000, -11.94000, -11.94000, -11.94000,  -8.95500, &
		-8.95500,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-8.95500,  -8.95500,  -8.95500,  -8.95500,  -5.97000, &
		-5.97000,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-5.97000,  -5.97000,  -5.97000,  -5.97000,  -2.98500, &
		-2.98500,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		-2.98500,  -2.98500,  -2.98500,  -2.98500,  -1.49250, &
		-1.49250,  -1.49250,  -1.49250,  -1.49250,  -1.49250, &
		-1.49250,  -1.49250,  -1.49250,  -1.49250,   1.49250, &
		1.49250,   1.49250,   1.49250,   1.49250,   1.49250, &
		1.49250,   1.49250,   1.49250,   1.49250,   4.47750, &
		4.47750,   4.47750,   4.47750,   4.47750,   4.47750, &
		4.47750,   4.47750,   4.47750,   4.47750,   7.46250, &
		7.46250,   7.46250,   7.46250,   7.46250,   7.46250, &
		7.46250,   7.46250,   7.46250,   7.46250,  10.44750, &
		10.44750,  10.44750,  10.44750,  10.44750,  10.44750, &
		10.44750,  10.44750,  10.44750,  10.44750,  13.43250, &
		13.43250,  13.43250,  13.43250,  13.43250,  13.43250, &
		13.43250,  13.43250,  13.43250,  13.43250,  16.41750, &
		16.41750,  16.41750,  16.41750,  16.41750,  16.41750, &
		16.41750,  16.41750,  16.41750,  16.41750,  19.40250, &
		19.40250,  19.40250,  19.40250,  19.40250,  19.40250, &
		19.40250,  19.40250,  19.40250,  19.40250,  22.38750, &
		22.38750,  22.38750,  22.38750,  22.38750,  22.38750, &
		22.38750,  22.38750,  22.38750,  22.38750,  25.37250, &
		25.37250,  25.37250,  25.37250,  25.37250,  25.37250, &
		25.37250,  25.37250,  25.37250,  25.37250,  28.35750, &
		28.35750,  28.35750,  28.35750,  28.35750,  28.35750, &
		28.35750,  28.35750,  28.35750,  28.35750, -28.35750, &
		-28.35750, -28.35750, -28.35750, -28.35750, -28.35750, &
		-28.35750, -28.35750, -28.35750, -28.35750, -25.37250, &
		-25.37250, -25.37250, -25.37250, -25.37250, -25.37250, &
		-25.37250, -25.37250, -25.37250, -25.37250, -22.38750, &
		-22.38750, -22.38750, -22.38750, -22.38750, -22.38750, &
		-22.38750, -22.38750, -22.38750, -22.38750, -19.40250, &
		-19.40250, -19.40250, -19.40250, -19.40250, -19.40250, &
		-19.40250, -19.40250, -19.40250, -19.40250, -16.41750, &
		-16.41750, -16.41750, -16.41750, -16.41750, -16.41750, &
		-16.41750, -16.41750, -16.41750, -16.41750, -13.43250, &
		-13.43250, -13.43250, -13.43250, -13.43250, -13.43250, &
		-13.43250, -13.43250, -13.43250, -13.43250, -10.44750, &
		-10.44750, -10.44750, -10.44750, -10.44750, -10.44750, &
		-10.44750, -10.44750, -10.44750, -10.44750,  -7.46250, &
		-7.46250,  -7.46250,  -7.46250,  -7.46250,  -7.46250, &
		-7.46250,  -7.46250,  -7.46250,  -7.46250,  -4.47750, &
		-4.47750,  -4.47750,  -4.47750,  -4.47750,  -4.47750, &
		-4.47750,  -4.47750,  -4.47750,  -4.47750,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   2.98500, &
		2.98500,   2.98500,   2.98500,   2.98500,   2.98500, &
		2.98500,   2.98500,   2.98500,   2.98500,   5.97000, &
		5.97000,   5.97000,   5.97000,   5.97000,   5.97000, &
		5.97000,   5.97000,   5.97000,   5.97000,   8.95500, &
		8.95500,   8.95500,   8.95500,   8.95500,   8.95500, &
		8.95500,   8.95500,   8.95500,   8.95500,  11.94000, &
		11.94000,  11.94000,  11.94000,  11.94000,  11.94000, &
		11.94000,  11.94000,  11.94000,  11.94000,  14.92500, &
		14.92500,  14.92500,  14.92500,  14.92500,  14.92500, &
		14.92500,  14.92500,  14.92500,  14.92500,  17.91000, &
		17.91000,  17.91000,  17.91000,  17.91000,  17.91000, &
		17.91000,  17.91000,  17.91000,  17.91000,  20.89500, &
		20.89500,  20.89500,  20.89500,  20.89500,  20.89500, &
		20.89500,  20.89500,  20.89500,  20.89500,  23.88000, &
		23.88000,  23.88000,  23.88000,  23.88000,  23.88000, &
		23.88000,  23.88000,  23.88000,  23.88000,  26.86500, &
		26.86500,  26.86500,  26.86500,  26.86500,  26.86500, &
		26.86500,  26.86500,  26.86500,  26.86500, -29.85000, &
		-29.85000, -29.85000, -29.85000, -29.85000, -29.85000, &
		-29.85000, -29.85000, -29.85000, -29.85000, -26.86500, &
		-26.86500, -26.86500, -26.86500, -26.86500, -26.86500, &
		-26.86500, -26.86500, -26.86500, -26.86500, -23.88000, &
		-23.88000, -23.88000, -23.88000, -23.88000, -23.88000, &
		-23.88000, -23.88000, -23.88000, -23.88000, -20.89500, &
		-20.89500, -20.89500, -20.89500, -20.89500, -20.89500, &
		-20.89500, -20.89500, -20.89500, -20.89500, -17.91000, &
		-17.91000, -17.91000, -17.91000, -17.91000, -17.91000, &
		-17.91000, -17.91000, -17.91000, -17.91000, -14.92500, &
		-14.92500, -14.92500, -14.92500, -14.92500, -14.92500, &
		-14.92500, -14.92500, -14.92500, -14.92500, -11.94000, &
		-11.94000, -11.94000, -11.94000, -11.94000, -11.94000, &
		-11.94000, -11.94000, -11.94000, -11.94000,  -8.95500, &
		-8.95500,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-8.95500,  -8.95500,  -8.95500,  -8.95500,  -5.97000, &
		-5.97000,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-5.97000,  -5.97000,  -5.97000,  -5.97000,  -2.98500, &
		-2.98500,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		-2.98500,  -2.98500,  -2.98500,  -2.98500,  -1.49250, &
		-1.49250,  -1.49250,  -1.49250,  -1.49250,  -1.49250, &
		-1.49250,  -1.49250,  -1.49250,  -1.49250,   1.49250, &
		1.49250,   1.49250,   1.49250,   1.49250,   1.49250, &
		1.49250,   1.49250,   1.49250,   1.49250,   4.47750, &
		4.47750,   4.47750,   4.47750,   4.47750,   4.47750, &
		4.47750,   4.47750,   4.47750,   4.47750,   7.46250, &
		7.46250,   7.46250,   7.46250,   7.46250,   7.46250, &
		7.46250,   7.46250,   7.46250,   7.46250,  10.44750, &
		10.44750,  10.44750,  10.44750,  10.44750,  10.44750, &
		10.44750,  10.44750,  10.44750,  10.44750,  13.43250, &
		13.43250,  13.43250,  13.43250,  13.43250,  13.43250, &
		13.43250,  13.43250,  13.43250,  13.43250,  16.41750, &
		16.41750,  16.41750,  16.41750,  16.41750,  16.41750, &
		16.41750,  16.41750,  16.41750,  16.41750,  19.40250, &
		19.40250,  19.40250,  19.40250,  19.40250,  19.40250, &
		19.40250,  19.40250,  19.40250,  19.40250,  22.38750, &
		22.38750,  22.38750,  22.38750,  22.38750,  22.38750, &
		22.38750,  22.38750,  22.38750,  22.38750,  25.37250, &
		25.37250,  25.37250,  25.37250,  25.37250,  25.37250, &
		25.37250,  25.37250,  25.37250,  25.37250,  28.35750, &
		28.35750,  28.35750,  28.35750,  28.35750,  28.35750, &
		28.35750,  28.35750,  28.35750,  28.35750, -28.35750, &
		-28.35750, -28.35750, -28.35750, -28.35750, -28.35750, &
		-28.35750, -28.35750, -28.35750, -28.35750, -25.37250, &
		-25.37250, -25.37250, -25.37250, -25.37250, -25.37250, &
		-25.37250, -25.37250, -25.37250, -25.37250, -22.38750, &
		-22.38750, -22.38750, -22.38750, -22.38750, -22.38750, &
		-22.38750, -22.38750, -22.38750, -22.38750, -19.40250, &
		-19.40250, -19.40250, -19.40250, -19.40250, -19.40250, &
		-19.40250, -19.40250, -19.40250, -19.40250, -16.41750, &
		-16.41750, -16.41750, -16.41750, -16.41750, -16.41750, &
		-16.41750, -16.41750, -16.41750, -16.41750, -13.43250, &
		-13.43250, -13.43250, -13.43250, -13.43250, -13.43250, &
		-13.43250, -13.43250, -13.43250, -13.43250, -10.44750, &
		-10.44750, -10.44750, -10.44750, -10.44750, -10.44750, &
		-10.44750, -10.44750, -10.44750, -10.44750,  -7.46250, &
		-7.46250,  -7.46250,  -7.46250,  -7.46250,  -7.46250, &
		-7.46250,  -7.46250,  -7.46250,  -7.46250,  -4.47750, &
		-4.47750,  -4.47750,  -4.47750,  -4.47750,  -4.47750, &
		-4.47750,  -4.47750,  -4.47750,  -4.47750 &
		]

		yTi = [ &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		-16.23871,  -3.24774,  22.73419, -29.22967, -22.73419, &
		-16.23871,  -9.74322,  -3.24774,   3.24774,   9.74322, &
		16.23871,  22.73419,  29.22967, -29.22967, -22.73419, &
		-16.23871,   9.74322,  22.73419, -16.23871,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-25.98193, -19.48645, -12.99097,  -6.49548,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322 &
		]

		zTi = [ &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.38548,   6.38548,   6.38548,   6.38548,   6.38548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		6.74548,   6.74548,   6.74548,   6.74548,   6.74548, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.43774,   3.43774,   3.43774,   3.43774,   3.43774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		3.16774,   3.16774,   3.16774,   3.16774,   3.16774, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.43774, &
		-3.43774,  -3.43774,  -3.43774,  -3.43774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -3.16774, &
		-3.16774,  -3.16774,  -3.16774,  -3.16774,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.38548, &
		-6.38548,  -6.38548,  -6.38548,  -6.38548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548,  -6.74548, &
		-6.74548,  -6.74548,  -6.74548,  -6.74548 &
		]
		
		xO = [ &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		2.98500,   2.98500,   2.98500,   2.98500,   2.98500, &
		2.98500,   2.98500,   2.98500,   2.98500,   2.98500, &
		5.97000,   5.97000,   5.97000,   5.97000,   5.97000, &
		5.97000,   5.97000,   5.97000,   5.97000,   5.97000, &
		8.95500,   8.95500,   8.95500,   8.95500,   8.95500, &
		8.95500,   8.95500,   8.95500,   8.95500,   8.95500, &
		11.94000,  11.94000,  11.94000,  11.94000,  11.94000, &
		11.94000,  11.94000,  11.94000,  11.94000,  11.94000, &
		14.92500,  14.92500,  14.92500,  14.92500,  14.92500, &
		14.92500,  14.92500,  14.92500,  14.92500,  14.92500, &
		17.91000,  17.91000,  17.91000,  17.91000,  17.91000, &
		17.91000,  17.91000,  17.91000,  17.91000,  17.91000, &
		20.89500,  20.89500,  20.89500,  20.89500,  20.89500, &
		20.89500,  20.89500,  20.89500,  20.89500,  20.89500, &
		23.88000,  23.88000,  23.88000,  23.88000,  23.88000, &
		23.88000,  23.88000,  23.88000,  23.88000,  23.88000, &
		26.86500,  26.86500,  26.86500,  26.86500,  26.86500, &
		26.86500,  26.86500,  26.86500,  26.86500,  26.86500, &
		-29.85000, -29.85000, -29.85000, -29.85000, -29.85000, &
		-29.85000, -29.85000, -29.85000, -29.85000, -29.85000, &
		-26.86500, -26.86500, -26.86500, -26.86500, -26.86500, &
		-26.86500, -26.86500, -26.86500, -26.86500, -26.86500, &
		-23.88000, -23.88000, -23.88000, -23.88000, -23.88000, &
		-23.88000, -23.88000, -23.88000, -23.88000, -23.88000, &
		-20.89500, -20.89500, -20.89500, -20.89500, -20.89500, &
		-20.89500, -20.89500, -20.89500, -20.89500, -20.89500, &
		-17.91000, -17.91000, -17.91000, -17.91000, -17.91000, &
		-17.91000, -17.91000, -17.91000, -17.91000, -17.91000, &
		-14.92500, -14.92500, -14.92500, -14.92500, -14.92500, &
		-14.92500, -14.92500, -14.92500, -14.92500, -14.92500, &
		-11.94000, -11.94000, -11.94000, -11.94000, -11.94000, &
		-11.94000, -11.94000, -11.94000, -11.94000, -11.94000, &
		-8.95500,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-8.95500,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-5.97000,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-5.97000,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-2.98500,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		-2.98500,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		1.49250,   1.49250,   1.49250,   1.49250,   1.49250, &
		1.49250,   1.49250,   1.49250,   1.49250,   1.49250, &
		4.47750,   4.47750,   4.47750,   4.47750,   4.47750, &
		4.47750,   4.47750,   4.47750,   4.47750,   4.47750, &
		7.46250,   7.46250,   7.46250,   7.46250,   7.46250, &
		7.46250,   7.46250,   7.46250,   7.46250,   7.46250, &
		10.44750,  10.44750,  10.44750,  10.44750,  10.44750, &
		10.44750,  10.44750,  10.44750,  10.44750,  10.44750, &
		13.43250,  13.43250,  13.43250,  13.43250,  13.43250, &
		13.43250,  13.43250,  13.43250,  13.43250,  13.43250, &
		16.41750,  16.41750,  16.41750,  16.41750,  16.41750, &
		16.41750,  16.41750,  16.41750,  16.41750,  16.41750, &
		19.40250,  19.40250,  19.40250,  19.40250,  19.40250, &
		19.40250,  19.40250,  19.40250,  19.40250,  19.40250, &
		22.38750,  22.38750,  22.38750,  22.38750,  22.38750, &
		22.38750,  22.38750,  22.38750,  22.38750,  22.38750, &
		25.37250,  25.37250,  25.37250,  25.37250,  25.37250, &
		25.37250,  25.37250,  25.37250,  25.37250,  25.37250, &
		28.35750,  28.35750,  28.35750,  28.35750,  28.35750, &
		28.35750,  28.35750,  28.35750,  28.35750,  28.35750, &
		-28.35750, -28.35750, -28.35750, -28.35750, -28.35750, &
		-28.35750, -28.35750, -28.35750, -28.35750, -28.35750, &
		-25.37250, -25.37250, -25.37250, -25.37250, -25.37250, &
		-25.37250, -25.37250, -25.37250, -25.37250, -25.37250, &
		-22.38750, -22.38750, -22.38750, -22.38750, -22.38750, &
		-22.38750, -22.38750, -22.38750, -22.38750, -22.38750, &
		-19.40250, -19.40250, -19.40250, -19.40250, -19.40250, &
		-19.40250, -19.40250, -19.40250, -19.40250, -19.40250, &
		-16.41750, -16.41750, -16.41750, -16.41750, -16.41750, &
		-16.41750, -16.41750, -16.41750, -16.41750, -16.41750, &
		-13.43250, -13.43250, -13.43250, -13.43250, -13.43250, &
		-13.43250, -13.43250, -13.43250, -13.43250, -13.43250, &
		-10.44750, -10.44750, -10.44750, -10.44750, -10.44750, &
		-10.44750, -10.44750, -10.44750, -10.44750, -10.44750, &
		-7.46250,  -7.46250,  -7.46250,  -7.46250,  -7.46250, &
		-7.46250,  -7.46250,  -7.46250,  -7.46250,  -7.46250, &
		-4.47750,  -4.47750,  -4.47750,  -4.47750,  -4.47750, &
		-4.47750,  -4.47750,  -4.47750,  -4.47750,  -4.47750, &
		-1.49250,  -1.49250,  -1.49250,  -1.49250,  -1.49250, &
		-1.49250,  -1.49250,  -1.49250,  -1.49250,  -1.49250, &
		-1.49250,  -1.49250,  -1.49250,  -1.49250,  -1.49250, &
		-1.49250,  -1.49250,  -1.49250,  -1.49250,  -1.49250, &
		1.49250,   1.49250,   1.49250,   1.49250,   1.49250, &
		1.49250,   1.49250,   1.49250,   1.49250,   1.49250, &
		4.47750,   4.47750,   4.47750,   4.47750,   4.47750, &
		4.47750,   4.47750,   4.47750,   4.47750,   4.47750, &
		7.46250,   7.46250,   7.46250,   7.46250,   7.46250, &
		7.46250,   7.46250,   7.46250,   7.46250,   7.46250, &
		10.44750,  10.44750,  10.44750,  10.44750,  10.44750, &
		10.44750,  10.44750,  10.44750,  10.44750,  10.44750, &
		13.43250,  13.43250,  13.43250,  13.43250,  13.43250, &
		13.43250,  13.43250,  13.43250,  13.43250,  13.43250, &
		16.41750,  16.41750,  16.41750,  16.41750,  16.41750, &
		16.41750,  16.41750,  16.41750,  16.41750,  16.41750, &
		19.40250,  19.40250,  19.40250,  19.40250,  19.40250, &
		19.40250,  19.40250,  19.40250,  19.40250,  19.40250, &
		22.38750,  22.38750,  22.38750,  22.38750,  22.38750, &
		22.38750,  22.38750,  22.38750,  22.38750,  22.38750, &
		25.37250,  25.37250,  25.37250,  25.37250,  25.37250, &
		25.37250,  25.37250,  25.37250,  25.37250,  25.37250, &
		28.35750,  28.35750,  28.35750,  28.35750,  28.35750, &
		28.35750,  28.35750,  28.35750,  28.35750,  28.35750, &
		-28.35750, -28.35750, -28.35750, -28.35750, -28.35750, &
		-28.35750, -28.35750, -28.35750, -28.35750, -28.35750, &
		-25.37250, -25.37250, -25.37250, -25.37250, -25.37250, &
		-25.37250, -25.37250, -25.37250, -25.37250, -25.37250, &
		-22.38750, -22.38750, -22.38750, -22.38750, -22.38750, &
		-22.38750, -22.38750, -22.38750, -22.38750, -22.38750, &
		-19.40250, -19.40250, -19.40250, -19.40250, -19.40250, &
		-19.40250, -19.40250, -19.40250, -19.40250, -19.40250, &
		-16.41750, -16.41750, -16.41750, -16.41750, -16.41750, &
		-16.41750, -16.41750, -16.41750, -16.41750, -16.41750, &
		-13.43250, -13.43250, -13.43250, -13.43250, -13.43250, &
		-13.43250, -13.43250, -13.43250, -13.43250, -13.43250, &
		-10.44750, -10.44750, -10.44750, -10.44750, -10.44750, &
		-10.44750, -10.44750, -10.44750, -10.44750, -10.44750, &
		-7.46250,  -7.46250,  -7.46250,  -7.46250,  -7.46250, &
		-7.46250,  -7.46250,  -7.46250,  -7.46250,  -7.46250, &
		-4.47750,  -4.47750,  -4.47750,  -4.47750,  -4.47750, &
		-4.47750,  -4.47750,  -4.47750,  -4.47750,  -4.47750, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		2.98500,   2.98500,   2.98500,   2.98500,   2.98500, &
		2.98500,   2.98500,   2.98500,   2.98500,   2.98500, &
		5.97000,   5.97000,   5.97000,   5.97000,   5.97000, &
		5.97000,   5.97000,   5.97000,   5.97000,   5.97000, &
		8.95500,   8.95500,   8.95500,   8.95500,   8.95500, &
		8.95500,   8.95500,   8.95500,   8.95500,   8.95500, &
		11.94000,  11.94000,  11.94000,  11.94000,  11.94000, &
		11.94000,  11.94000,  11.94000,  11.94000,  11.94000, &
		14.92500,  14.92500,  14.92500,  14.92500,  14.92500, &
		14.92500,  14.92500,  14.92500,  14.92500,  14.92500, &
		17.91000,  17.91000,  17.91000,  17.91000,  17.91000, &
		17.91000,  17.91000,  17.91000,  17.91000,  17.91000, &
		20.89500,  20.89500,  20.89500,  20.89500,  20.89500, &
		20.89500,  20.89500,  20.89500,  20.89500,  20.89500, &
		23.88000,  23.88000,  23.88000,  23.88000,  23.88000, &
		23.88000,  23.88000,  23.88000,  23.88000,  23.88000, &
		26.86500,  26.86500,  26.86500,  26.86500,  26.86500, &
		26.86500,  26.86500,  26.86500,  26.86500,  26.86500, &
		-29.85000, -29.85000, -29.85000, -29.85000, -29.85000, &
		-29.85000, -29.85000, -29.85000, -29.85000, -29.85000, &
		-26.86500, -26.86500, -26.86500, -26.86500, -26.86500, &
		-26.86500, -26.86500, -26.86500, -26.86500, -26.86500, &
		-23.88000, -23.88000, -23.88000, -23.88000, -23.88000, &
		-23.88000, -23.88000, -23.88000, -23.88000, -23.88000, &
		-20.89500, -20.89500, -20.89500, -20.89500, -20.89500, &
		-20.89500, -20.89500, -20.89500, -20.89500, -20.89500, &
		-17.91000, -17.91000, -17.91000, -17.91000, -17.91000, &
		-17.91000, -17.91000, -17.91000, -17.91000, -17.91000, &
		-14.92500, -14.92500, -14.92500, -14.92500, -14.92500, &
		-14.92500, -14.92500, -14.92500, -14.92500, -14.92500, &
		-11.94000, -11.94000, -11.94000, -11.94000, -11.94000, &
		-11.94000, -11.94000, -11.94000, -11.94000, -11.94000, &
		-8.95500,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-8.95500,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-5.97000,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-5.97000,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-2.98500,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		-2.98500,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		2.98500,   2.98500,   2.98500,   2.98500,   2.98500, &
		2.98500,   2.98500,   2.98500,   2.98500,   2.98500, &
		5.97000,   5.97000,   5.97000,   5.97000,   5.97000, &
		5.97000,   5.97000,   5.97000,   5.97000,   5.97000, &
		8.95500,   8.95500,   8.95500,   8.95500,   8.95500, &
		8.95500,   8.95500,   8.95500,   8.95500,   8.95500, &
		11.94000,  11.94000,  11.94000,  11.94000,  11.94000, &
		11.94000,  11.94000,  11.94000,  11.94000,  11.94000, &
		14.92500,  14.92500,  14.92500,  14.92500,  14.92500, &
		14.92500,  14.92500,  14.92500,  14.92500,  14.92500, &
		17.91000,  17.91000,  17.91000,  17.91000,  17.91000, &
		17.91000,  17.91000,  17.91000,  17.91000,  17.91000, &
		20.89500,  20.89500,  20.89500,  20.89500,  20.89500, &
		20.89500,  20.89500,  20.89500,  20.89500,  20.89500, &
		23.88000,  23.88000,  23.88000,  23.88000,  23.88000, &
		23.88000,  23.88000,  23.88000,  23.88000,  23.88000, &
		26.86500,  26.86500,  26.86500,  26.86500,  26.86500, &
		26.86500,  26.86500,  26.86500,  26.86500,  26.86500, &
		-29.85000, -29.85000, -29.85000, -29.85000, -29.85000, &
		-29.85000, -29.85000, -29.85000, -29.85000, -29.85000, &
		-26.86500, -26.86500, -26.86500, -26.86500, -26.86500, &
		-26.86500, -26.86500, -26.86500, -26.86500, -26.86500, &
		-23.88000, -23.88000, -23.88000, -23.88000, -23.88000, &
		-23.88000, -23.88000, -23.88000, -23.88000, -23.88000, &
		-20.89500, -20.89500, -20.89500, -20.89500, -20.89500, &
		-20.89500, -20.89500, -20.89500, -20.89500, -20.89500, &
		-17.91000, -17.91000, -17.91000, -17.91000, -17.91000, &
		-17.91000, -17.91000, -17.91000, -17.91000, -17.91000, &
		-14.92500, -14.92500, -14.92500, -14.92500, -14.92500, &
		-14.92500, -14.92500, -14.92500, -14.92500, -14.92500, &
		-11.94000, -11.94000, -11.94000, -11.94000, -11.94000, &
		-11.94000, -11.94000, -11.94000, -11.94000, -11.94000, &
		-8.95500,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-8.95500,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-5.97000,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-5.97000,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-2.98500,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		-2.98500,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		-1.49250,  -1.49250,  -1.49250,  -1.49250,  -1.49250, &
		-1.49250,  -1.49250,  -1.49250,  -1.49250,  -1.49250, &
		1.49250,   1.49250,   1.49250,   1.49250,   1.49250, &
		1.49250,   1.49250,   1.49250,   1.49250,   1.49250, &
		4.47750,   4.47750,   4.47750,   4.47750,   4.47750, &
		4.47750,   4.47750,   4.47750,   4.47750,   4.47750, &
		7.46250,   7.46250,   7.46250,   7.46250,   7.46250, &
		7.46250,   7.46250,   7.46250,   7.46250,   7.46250, &
		10.44750,  10.44750,  10.44750,  10.44750,  10.44750, &
		10.44750,  10.44750,  10.44750,  10.44750,  10.44750, &
		13.43250,  13.43250,  13.43250,  13.43250,  13.43250, &
		13.43250,  13.43250,  13.43250,  13.43250,  13.43250, &
		16.41750,  16.41750,  16.41750,  16.41750,  16.41750, &
		16.41750,  16.41750,  16.41750,  16.41750,  16.41750, &
		19.40250,  19.40250,  19.40250,  19.40250,  19.40250, &
		19.40250,  19.40250,  19.40250,  19.40250,  19.40250, &
		22.38750,  22.38750,  22.38750,  22.38750,  22.38750, &
		22.38750,  22.38750,  22.38750,  22.38750,  22.38750, &
		25.37250,  25.37250,  25.37250,  25.37250,  25.37250, &
		25.37250,  25.37250,  25.37250,  25.37250,  25.37250, &
		28.35750,  28.35750,  28.35750,  28.35750,  28.35750, &
		28.35750,  28.35750,  28.35750,  28.35750,  28.35750, &
		-28.35750, -28.35750, -28.35750, -28.35750, -28.35750, &
		-28.35750, -28.35750, -28.35750, -28.35750, -28.35750, &
		-25.37250, -25.37250, -25.37250, -25.37250, -25.37250, &
		-25.37250, -25.37250, -25.37250, -25.37250, -25.37250, &
		-22.38750, -22.38750, -22.38750, -22.38750, -22.38750, &
		-22.38750, -22.38750, -22.38750, -22.38750, -22.38750, &
		-19.40250, -19.40250, -19.40250, -19.40250, -19.40250, &
		-19.40250, -19.40250, -19.40250, -19.40250, -19.40250, &
		-16.41750, -16.41750, -16.41750, -16.41750, -16.41750, &
		-16.41750, -16.41750, -16.41750, -16.41750, -16.41750, &
		-13.43250, -13.43250, -13.43250, -13.43250, -13.43250, &
		-13.43250, -13.43250, -13.43250, -13.43250, -13.43250, &
		-10.44750, -10.44750, -10.44750, -10.44750, -10.44750, &
		-10.44750, -10.44750, -10.44750, -10.44750, -10.44750, &
		-7.46250,  -7.46250,  -7.46250,  -7.46250,  -7.46250, &
		-7.46250,  -7.46250,  -7.46250,  -7.46250,  -7.46250, &
		-4.47750,  -4.47750,  -4.47750,  -4.47750,  -4.47750, &
		-4.47750,  -4.47750,  -4.47750,  -4.47750,  -4.47750, &
		1.49250,   1.49250,   1.49250,   1.49250,   1.49250, &
		1.49250,   1.49250,   4.47750,   4.47750,   4.47750, &
		4.47750,   4.47750,   4.47750,   7.46250,   7.46250, &
		7.46250,   7.46250,   7.46250,   7.46250,  10.44750, &
		10.44750,  10.44750,  10.44750,  10.44750,  10.44750, &
		13.43250,  13.43250,  13.43250,  13.43250,  13.43250, &
		13.43250,  16.41750,  16.41750,  16.41750,  16.41750, &
		16.41750,  16.41750,  19.40250,  19.40250,  19.40250, &
		19.40250,  19.40250,  19.40250,  22.38750,  22.38750, &
		22.38750,  22.38750,  22.38750,  22.38750,  25.37250, &
		25.37250,  25.37250,  25.37250,  28.35750,  28.35750, &
		28.35750,  28.35750,  28.35750,  28.35750, -28.35750, &
		-28.35750, -28.35750, -28.35750, -25.37250, -25.37250, &
		-25.37250, -25.37250, -25.37250, -25.37250, -22.38750, &
		-22.38750, -22.38750, -22.38750, -22.38750, -22.38750, &
		-19.40250, -19.40250, -19.40250, -19.40250, -19.40250, &
		-19.40250, -16.41750, -16.41750, -16.41750, -16.41750, &
		-16.41750, -16.41750, -13.43250, -13.43250, -13.43250, &
		-13.43250, -13.43250, -13.43250, -13.43250, -10.44750, &
		-10.44750, -10.44750, -10.44750, -10.44750, -10.44750, &
		-10.44750,  -7.46250,  -7.46250,  -7.46250,  -7.46250, &
		-7.46250,  -7.46250,  -4.47750,  -4.47750,  -4.47750, &
		-4.47750,  -4.47750,  -4.47750,  -1.49250,  -1.49250, &
		-1.49250,  -1.49250,  -1.49250,  -1.49250,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   2.98500,   2.98500,   2.98500,   2.98500, &
		2.98500,   2.98500,   2.98500,   5.97000,   5.97000, &
		5.97000,   5.97000,   5.97000,   5.97000,   5.97000, &
		5.97000,   8.95500,   8.95500,   8.95500,   8.95500, &
		8.95500,   8.95500,  11.94000,  11.94000,  11.94000, &
		11.94000,  11.94000,  11.94000,  11.94000,  14.92500, &
		14.92500,  14.92500,  14.92500,  14.92500,  14.92500, &
		14.92500,  17.91000,  17.91000,  17.91000,  17.91000, &
		17.91000,  17.91000,  20.89500,  20.89500,  20.89500, &
		20.89500,  20.89500,  20.89500,  20.89500,  23.88000, &
		23.88000,  23.88000,  23.88000,  23.88000,  23.88000, &
		23.88000,  26.86500,  26.86500,  26.86500,  26.86500, &
		26.86500,  26.86500,  26.86500,  26.86500,  26.86500, &
		26.86500, -29.85000, -29.85000, -29.85000, -29.85000, &
		-29.85000, -29.85000, -29.85000, -29.85000, -29.85000, &
		-29.85000, -26.86500, -26.86500, -26.86500, -26.86500, &
		-26.86500, -26.86500, -26.86500, -26.86500, -26.86500, &
		-26.86500, -23.88000, -23.88000, -23.88000, -23.88000, &
		-23.88000, -23.88000, -23.88000, -23.88000, -23.88000, &
		-23.88000, -20.89500, -20.89500, -20.89500, -20.89500, &
		-20.89500, -20.89500, -20.89500, -20.89500, -20.89500, &
		-20.89500, -17.91000, -17.91000, -17.91000, -17.91000, &
		-17.91000, -17.91000, -17.91000, -17.91000, -17.91000, &
		-17.91000, -14.92500, -14.92500, -14.92500, -14.92500, &
		-14.92500, -14.92500, -14.92500, -14.92500, -14.92500, &
		-14.92500, -11.94000, -11.94000, -11.94000, -11.94000, &
		-11.94000, -11.94000, -11.94000, -11.94000, -11.94000, &
		-11.94000,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-8.95500,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-8.95500,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-5.97000,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-5.97000,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		-2.98500,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		-2.98500,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   2.98500,   2.98500,   2.98500,   2.98500, &
		2.98500,   2.98500,   2.98500,   2.98500,   2.98500, &
		2.98500,   5.97000,   5.97000,   5.97000,   5.97000, &
		5.97000,   5.97000,   5.97000,   5.97000,   5.97000, &
		5.97000,   8.95500,   8.95500,   8.95500,   8.95500, &
		8.95500,   8.95500,   8.95500,   8.95500,   8.95500, &
		8.95500,  11.94000,  11.94000,  11.94000,  11.94000, &
		11.94000,  11.94000,  11.94000,  11.94000,  11.94000, &
		11.94000,  14.92500,  14.92500,  14.92500,  14.92500, &
		14.92500,  14.92500,  14.92500,  14.92500,  14.92500, &
		14.92500,  17.91000,  17.91000,  17.91000,  17.91000, &
		17.91000,  17.91000,  17.91000,  17.91000,  17.91000, &
		17.91000,  20.89500,  20.89500,  20.89500,  20.89500, &
		20.89500,  20.89500,  20.89500,  20.89500,  20.89500, &
		20.89500,  23.88000,  23.88000,  23.88000,  23.88000, &
		23.88000,  23.88000,  23.88000,  23.88000,  23.88000, &
		23.88000,  26.86500,  26.86500,  26.86500,  26.86500, &
		26.86500,  26.86500,  26.86500,  26.86500,  26.86500, &
		26.86500, -29.85000, -29.85000, -29.85000, -29.85000, &
		-29.85000, -29.85000, -29.85000, -29.85000, -29.85000, &
		-29.85000, -26.86500, -26.86500, -26.86500, -26.86500, &
		-26.86500, -26.86500, -26.86500, -26.86500, -26.86500, &
		-26.86500, -23.88000, -23.88000, -23.88000, -23.88000, &
		-23.88000, -23.88000, -23.88000, -23.88000, -23.88000, &
		-23.88000, -20.89500, -20.89500, -20.89500, -20.89500, &
		-20.89500, -20.89500, -20.89500, -20.89500, -20.89500, &
		-20.89500, -17.91000, -17.91000, -17.91000, -17.91000, &
		-17.91000, -17.91000, -17.91000, -17.91000, -17.91000, &
		-17.91000, -14.92500, -14.92500, -14.92500, -14.92500, &
		-14.92500, -14.92500, -14.92500, -14.92500, -14.92500, &
		-14.92500, -11.94000, -11.94000, -11.94000, -11.94000, &
		-11.94000, -11.94000, -11.94000, -11.94000, -11.94000, &
		-11.94000,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-8.95500,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-8.95500,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-5.97000,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-5.97000,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		-2.98500,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		-2.98500,  -1.49250,  -1.49250,  -1.49250,  -1.49250, &
		-1.49250,  -1.49250,  -1.49250,  -1.49250,  -1.49250, &
		-1.49250,   1.49250,   1.49250,   1.49250,   1.49250, &
		1.49250,   1.49250,   1.49250,   1.49250,   1.49250, &
		1.49250,   4.47750,   4.47750,   4.47750,   4.47750, &
		4.47750,   4.47750,   4.47750,   4.47750,   4.47750, &
		4.47750,   7.46250,   7.46250,   7.46250,   7.46250, &
		7.46250,   7.46250,   7.46250,   7.46250,   7.46250, &
		7.46250,  10.44750,  10.44750,  10.44750,  10.44750, &
		10.44750,  10.44750,  10.44750,  10.44750,  10.44750, &
		10.44750,  13.43250,  13.43250,  13.43250,  13.43250, &
		13.43250,  13.43250,  13.43250,  13.43250,  13.43250, &
		13.43250,  16.41750,  16.41750,  16.41750,  16.41750, &
		16.41750,  16.41750,  16.41750,  16.41750,  16.41750, &
		16.41750,  19.40250,  19.40250,  19.40250,  19.40250, &
		19.40250,  19.40250,  19.40250,  19.40250,  19.40250, &
		19.40250,  22.38750,  22.38750,  22.38750,  22.38750, &
		22.38750,  22.38750,  22.38750,  22.38750,  22.38750, &
		22.38750,  25.37250,  25.37250,  25.37250,  25.37250, &
		25.37250,  25.37250,  25.37250,  25.37250,  25.37250, &
		25.37250,  28.35750,  28.35750,  28.35750,  28.35750, &
		28.35750,  28.35750,  28.35750,  28.35750,  28.35750, &
		28.35750, -28.35750, -28.35750, -28.35750, -28.35750, &
		-28.35750, -28.35750, -28.35750, -28.35750, -28.35750, &
		-28.35750, -25.37250, -25.37250, -25.37250, -25.37250, &
		-25.37250, -25.37250, -25.37250, -25.37250, -25.37250, &
		-25.37250, -22.38750, -22.38750, -22.38750, -22.38750, &
		-22.38750, -22.38750, -22.38750, -22.38750, -22.38750, &
		-22.38750, -19.40250, -19.40250, -19.40250, -19.40250, &
		-19.40250, -19.40250, -19.40250, -19.40250, -19.40250, &
		-19.40250, -16.41750, -16.41750, -13.43250, -13.43250, &
		-10.44750, -10.44750, -10.44750, -10.44750, -10.44750, &
		-10.44750, -10.44750, -10.44750, -10.44750, -10.44750, &
		-7.46250,  -7.46250,  -7.46250,  -7.46250,  -7.46250, &
		-7.46250,  -7.46250,  -7.46250,  -7.46250,  -7.46250, &
		-4.47750,  -4.47750,  -4.47750,  -4.47750,  -4.47750, &
		-4.47750,  -4.47750,  -4.47750,  -4.47750,  -4.47750, &
		-1.49250,  -1.49250,   1.49250,   1.49250,   1.49250, &
		1.49250,   1.49250,   1.49250,   1.49250,   1.49250, &
		1.49250,   1.49250,   4.47750,   4.47750,   4.47750, &
		4.47750,   4.47750,   4.47750,   4.47750,   4.47750, &
		4.47750,   4.47750,   7.46250,   7.46250,   7.46250, &
		7.46250,   7.46250,   7.46250,   7.46250,   7.46250, &
		7.46250,   7.46250,  10.44750,  10.44750,  10.44750, &
		10.44750,  10.44750,  10.44750,  10.44750,  10.44750, &
		13.43250,  13.43250,  13.43250,  13.43250,  16.41750, &
		16.41750,  16.41750,  16.41750,  16.41750,  16.41750, &
		16.41750,  16.41750,  16.41750,  16.41750,  19.40250, &
		19.40250,  19.40250,  19.40250,  19.40250,  19.40250, &
		19.40250,  19.40250,  19.40250,  19.40250,  22.38750, &
		22.38750,  22.38750,  22.38750,  22.38750,  22.38750, &
		22.38750,  22.38750,  22.38750,  22.38750,  25.37250, &
		25.37250,  25.37250,  25.37250,  25.37250,  25.37250, &
		28.35750,  28.35750,  28.35750,  28.35750,  28.35750, &
		28.35750, -28.35750, -28.35750, -28.35750, -28.35750, &
		-28.35750, -28.35750, -28.35750, -28.35750, -28.35750, &
		-28.35750, -25.37250, -25.37250, -25.37250, -25.37250, &
		-25.37250, -25.37250, -25.37250, -25.37250, -25.37250, &
		-25.37250, -22.38750, -22.38750, -22.38750, -22.38750, &
		-22.38750, -22.38750, -22.38750, -22.38750, -22.38750, &
		-22.38750, -19.40250, -19.40250, -19.40250, -19.40250, &
		-16.41750, -16.41750, -16.41750, -16.41750, -16.41750, &
		-16.41750, -16.41750, -16.41750, -13.43250, -13.43250, &
		-13.43250, -13.43250, -13.43250, -13.43250, -13.43250, &
		-13.43250, -13.43250, -13.43250, -10.44750, -10.44750, &
		-10.44750, -10.44750, -10.44750, -10.44750, -10.44750, &
		-10.44750, -10.44750, -10.44750,  -7.46250,  -7.46250, &
		-7.46250,  -7.46250,  -7.46250,  -7.46250,  -7.46250, &
		-7.46250,  -7.46250,  -7.46250,  -4.47750,  -4.47750, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		2.98500,   2.98500,   2.98500,   2.98500,   2.98500, &
		2.98500,   2.98500,   2.98500,   2.98500,   2.98500, &
		5.97000,   5.97000,   5.97000,   5.97000,   5.97000, &
		5.97000,   5.97000,   5.97000,   5.97000,   5.97000, &
		8.95500,   8.95500,   8.95500,   8.95500,   8.95500, &
		8.95500,   8.95500,   8.95500,   8.95500,   8.95500, &
		11.94000,  11.94000,  14.92500,  14.92500,  14.92500, &
		14.92500,  14.92500,  14.92500,  14.92500,  14.92500, &
		14.92500,  14.92500,  17.91000,  17.91000,  17.91000, &
		17.91000,  17.91000,  17.91000,  17.91000,  17.91000, &
		17.91000,  17.91000,  20.89500,  20.89500,  20.89500, &
		20.89500,  20.89500,  20.89500,  20.89500,  20.89500, &
		20.89500,  20.89500,  23.88000,  23.88000,  23.88000, &
		23.88000,  23.88000,  23.88000,  23.88000,  23.88000, &
		26.86500,  26.86500,  26.86500,  26.86500, -29.85000, &
		-29.85000, -29.85000, -29.85000, -29.85000, -29.85000, &
		-29.85000, -29.85000, -29.85000, -29.85000, -26.86500, &
		-26.86500, -26.86500, -26.86500, -26.86500, -26.86500, &
		-26.86500, -26.86500, -26.86500, -26.86500, -23.88000, &
		-23.88000, -23.88000, -23.88000, -23.88000, -23.88000, &
		-23.88000, -23.88000, -23.88000, -23.88000, -20.89500, &
		-20.89500, -20.89500, -20.89500, -20.89500, -20.89500, &
		-17.91000, -17.91000, -17.91000, -17.91000, -17.91000, &
		-17.91000, -14.92500, -14.92500, -14.92500, -14.92500, &
		-14.92500, -14.92500, -14.92500, -14.92500, -14.92500, &
		-14.92500, -11.94000, -11.94000, -11.94000, -11.94000, &
		-11.94000, -11.94000, -11.94000, -11.94000, -11.94000, &
		-11.94000,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-8.95500,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-8.95500,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-5.97000,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-5.97000,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		-2.98500,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		-2.98500,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   2.98500,   2.98500,   2.98500,   2.98500, &
		2.98500,   2.98500,   2.98500,   2.98500,   2.98500, &
		2.98500,   5.97000,   5.97000,   5.97000,   5.97000, &
		5.97000,   5.97000,   5.97000,   5.97000,   5.97000, &
		5.97000,   8.95500,   8.95500,   8.95500,   8.95500, &
		8.95500,   8.95500,   8.95500,   8.95500,   8.95500, &
		8.95500,  11.94000,  11.94000,  11.94000,  11.94000, &
		11.94000,  11.94000,  11.94000,  11.94000,  11.94000, &
		11.94000,  14.92500,  14.92500,  14.92500,  14.92500, &
		14.92500,  14.92500,  14.92500,  14.92500,  14.92500, &
		14.92500,  17.91000,  17.91000,  17.91000,  17.91000, &
		17.91000,  17.91000,  17.91000,  17.91000,  17.91000, &
		17.91000,  20.89500,  20.89500,  20.89500,  20.89500, &
		20.89500,  20.89500,  20.89500,  20.89500,  20.89500, &
		20.89500,  23.88000,  23.88000,  23.88000,  23.88000, &
		23.88000,  23.88000,  23.88000,  23.88000,  23.88000, &
		23.88000,  26.86500,  26.86500,  26.86500,  26.86500, &
		26.86500,  26.86500,  26.86500,  26.86500,  26.86500, &
		26.86500, -29.85000, -29.85000, -29.85000, -29.85000, &
		-29.85000, -29.85000, -29.85000, -29.85000, -29.85000, &
		-29.85000, -26.86500, -26.86500, -26.86500, -26.86500, &
		-26.86500, -26.86500, -26.86500, -26.86500, -26.86500, &
		-26.86500, -23.88000, -23.88000, -23.88000, -23.88000, &
		-23.88000, -23.88000, -23.88000, -23.88000, -23.88000, &
		-23.88000, -20.89500, -20.89500, -20.89500, -20.89500, &
		-20.89500, -20.89500, -20.89500, -20.89500, -20.89500, &
		-20.89500, -17.91000, -17.91000, -17.91000, -17.91000, &
		-17.91000, -17.91000, -17.91000, -17.91000, -17.91000, &
		-17.91000, -14.92500, -14.92500, -14.92500, -14.92500, &
		-14.92500, -14.92500, -14.92500, -14.92500, -14.92500, &
		-14.92500, -11.94000, -11.94000, -11.94000, -11.94000, &
		-11.94000, -11.94000, -11.94000, -11.94000, -11.94000, &
		-11.94000,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-8.95500,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-8.95500,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-5.97000,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-5.97000,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		-2.98500,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		-2.98500,   1.49250,   1.49250,   1.49250,   1.49250, &
		1.49250,   1.49250,   1.49250,   1.49250,   1.49250, &
		1.49250,   4.47750,   4.47750,   4.47750,   4.47750, &
		4.47750,   4.47750,   4.47750,   4.47750,   4.47750, &
		4.47750,   7.46250,   7.46250,   7.46250,   7.46250, &
		7.46250,   7.46250,   7.46250,   7.46250,   7.46250, &
		7.46250,  10.44750,  10.44750,  10.44750,  10.44750, &
		10.44750,  10.44750,  10.44750,  10.44750,  10.44750, &
		10.44750,  13.43250,  13.43250,  13.43250,  13.43250, &
		13.43250,  13.43250,  13.43250,  13.43250,  13.43250, &
		13.43250,  16.41750,  16.41750,  16.41750,  16.41750, &
		16.41750,  16.41750,  16.41750,  16.41750,  16.41750, &
		16.41750,  19.40250,  19.40250,  19.40250,  19.40250, &
		19.40250,  19.40250,  19.40250,  19.40250,  19.40250, &
		19.40250,  22.38750,  22.38750,  22.38750,  22.38750, &
		22.38750,  22.38750,  22.38750,  22.38750,  22.38750, &
		22.38750,  25.37250,  25.37250,  25.37250,  25.37250, &
		25.37250,  25.37250,  25.37250,  25.37250,  25.37250, &
		25.37250,  28.35750,  28.35750,  28.35750,  28.35750, &
		28.35750,  28.35750,  28.35750,  28.35750,  28.35750, &
		28.35750, -28.35750, -28.35750, -28.35750, -28.35750, &
		-28.35750, -28.35750, -28.35750, -28.35750, -28.35750, &
		-28.35750, -25.37250, -25.37250, -25.37250, -25.37250, &
		-25.37250, -25.37250, -25.37250, -25.37250, -25.37250, &
		-25.37250, -22.38750, -22.38750, -22.38750, -22.38750, &
		-22.38750, -22.38750, -22.38750, -22.38750, -22.38750, &
		-22.38750, -19.40250, -19.40250, -19.40250, -19.40250, &
		-19.40250, -19.40250, -19.40250, -19.40250, -19.40250, &
		-19.40250, -16.41750, -16.41750, -16.41750, -16.41750, &
		-16.41750, -16.41750, -16.41750, -16.41750, -16.41750, &
		-16.41750, -13.43250, -13.43250, -13.43250, -13.43250, &
		-13.43250, -13.43250, -13.43250, -13.43250, -13.43250, &
		-13.43250, -10.44750, -10.44750, -10.44750, -10.44750, &
		-10.44750, -10.44750, -10.44750, -10.44750, -10.44750, &
		-10.44750,  -7.46250,  -7.46250,  -7.46250,  -7.46250, &
		-7.46250,  -7.46250,  -7.46250,  -7.46250,  -7.46250, &
		-7.46250,  -4.47750,  -4.47750,  -4.47750,  -4.47750, &
		-4.47750,  -4.47750,  -4.47750,  -4.47750,  -4.47750, &
		-4.47750,  -1.49250,  -1.49250,  -1.49250,  -1.49250, &
		-1.49250,  -1.49250,  -1.49250,  -1.49250,  -1.49250, &
		-1.49250,  -1.49250,  -1.49250,  -1.49250,  -1.49250, &
		-1.49250,  -1.49250,  -1.49250,  -1.49250,  -1.49250, &
		-1.49250,   1.49250,   1.49250,   1.49250,   1.49250, &
		1.49250,   1.49250,   1.49250,   1.49250,   1.49250, &
		1.49250,   4.47750,   4.47750,   4.47750,   4.47750, &
		4.47750,   4.47750,   4.47750,   4.47750,   4.47750, &
		4.47750,   7.46250,   7.46250,   7.46250,   7.46250, &
		7.46250,   7.46250,   7.46250,   7.46250,   7.46250, &
		7.46250,  10.44750,  10.44750,  10.44750,  10.44750, &
		10.44750,  10.44750,  10.44750,  10.44750,  10.44750, &
		10.44750,  13.43250,  13.43250,  13.43250,  13.43250, &
		13.43250,  13.43250,  13.43250,  13.43250,  13.43250, &
		13.43250,  16.41750,  16.41750,  16.41750,  16.41750, &
		16.41750,  16.41750,  16.41750,  16.41750,  16.41750, &
		16.41750,  19.40250,  19.40250,  19.40250,  19.40250, &
		19.40250,  19.40250,  19.40250,  19.40250,  19.40250, &
		19.40250,  22.38750,  22.38750,  22.38750,  22.38750, &
		22.38750,  22.38750,  22.38750,  22.38750,  22.38750, &
		22.38750,  25.37250,  25.37250,  25.37250,  25.37250, &
		25.37250,  25.37250,  25.37250,  25.37250,  25.37250, &
		25.37250,  28.35750,  28.35750,  28.35750,  28.35750, &
		28.35750,  28.35750,  28.35750,  28.35750,  28.35750, &
		28.35750, -28.35750, -28.35750, -28.35750, -28.35750, &
		-28.35750, -28.35750, -28.35750, -28.35750, -28.35750, &
		-28.35750, -25.37250, -25.37250, -25.37250, -25.37250, &
		-25.37250, -25.37250, -25.37250, -25.37250, -25.37250, &
		-25.37250, -22.38750, -22.38750, -22.38750, -22.38750, &
		-22.38750, -22.38750, -22.38750, -22.38750, -22.38750, &
		-22.38750, -19.40250, -19.40250, -19.40250, -19.40250, &
		-19.40250, -19.40250, -19.40250, -19.40250, -19.40250, &
		-19.40250, -16.41750, -16.41750, -16.41750, -16.41750, &
		-16.41750, -16.41750, -16.41750, -16.41750, -16.41750, &
		-16.41750, -13.43250, -13.43250, -13.43250, -13.43250, &
		-13.43250, -13.43250, -13.43250, -13.43250, -13.43250, &
		-13.43250, -10.44750, -10.44750, -10.44750, -10.44750, &
		-10.44750, -10.44750, -10.44750, -10.44750, -10.44750, &
		-10.44750,  -7.46250,  -7.46250,  -7.46250,  -7.46250, &
		-7.46250,  -7.46250,  -7.46250,  -7.46250,  -7.46250, &
		-7.46250,  -4.47750,  -4.47750,  -4.47750,  -4.47750, &
		-4.47750,  -4.47750,  -4.47750,  -4.47750,  -4.47750, &
		-4.47750,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   2.98500,   2.98500,   2.98500,   2.98500, &
		2.98500,   2.98500,   2.98500,   2.98500,   2.98500, &
		2.98500,   5.97000,   5.97000,   5.97000,   5.97000, &
		5.97000,   5.97000,   5.97000,   5.97000,   5.97000, &
		5.97000,   8.95500,   8.95500,   8.95500,   8.95500, &
		8.95500,   8.95500,   8.95500,   8.95500,   8.95500, &
		8.95500,  11.94000,  11.94000,  11.94000,  11.94000, &
		11.94000,  11.94000,  11.94000,  11.94000,  11.94000, &
		11.94000,  14.92500,  14.92500,  14.92500,  14.92500, &
		14.92500,  14.92500,  14.92500,  14.92500,  14.92500, &
		14.92500,  17.91000,  17.91000,  17.91000,  17.91000, &
		17.91000,  17.91000,  17.91000,  17.91000,  17.91000, &
		17.91000,  20.89500,  20.89500,  20.89500,  20.89500, &
		20.89500,  20.89500,  20.89500,  20.89500,  20.89500, &
		20.89500,  23.88000,  23.88000,  23.88000,  23.88000, &
		23.88000,  23.88000,  23.88000,  23.88000,  23.88000, &
		23.88000,  26.86500,  26.86500,  26.86500,  26.86500, &
		26.86500,  26.86500,  26.86500,  26.86500,  26.86500, &
		26.86500, -29.85000, -29.85000, -29.85000, -29.85000, &
		-29.85000, -29.85000, -29.85000, -29.85000, -29.85000, &
		-29.85000, -26.86500, -26.86500, -26.86500, -26.86500, &
		-26.86500, -26.86500, -26.86500, -26.86500, -26.86500, &
		-26.86500, -23.88000, -23.88000, -23.88000, -23.88000, &
		-23.88000, -23.88000, -23.88000, -23.88000, -23.88000, &
		-23.88000, -20.89500, -20.89500, -20.89500, -20.89500, &
		-20.89500, -20.89500, -20.89500, -20.89500, -20.89500, &
		-20.89500, -17.91000, -17.91000, -17.91000, -17.91000, &
		-17.91000, -17.91000, -17.91000, -17.91000, -17.91000, &
		-17.91000, -14.92500, -14.92500, -14.92500, -14.92500, &
		-14.92500, -14.92500, -14.92500, -14.92500, -14.92500, &
		-14.92500, -11.94000, -11.94000, -11.94000, -11.94000, &
		-11.94000, -11.94000, -11.94000, -11.94000, -11.94000, &
		-11.94000,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-8.95500,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-8.95500,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-5.97000,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-5.97000,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		-2.98500,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		-2.98500,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   2.98500,   2.98500,   2.98500,   2.98500, &
		2.98500,   2.98500,   2.98500,   2.98500,   2.98500, &
		2.98500,   5.97000,   5.97000,   5.97000,   5.97000, &
		5.97000,   5.97000,   5.97000,   5.97000,   5.97000, &
		5.97000,   8.95500,   8.95500,   8.95500,   8.95500, &
		8.95500,   8.95500,   8.95500,   8.95500,   8.95500, &
		8.95500,  11.94000,  11.94000,  11.94000,  11.94000, &
		11.94000,  11.94000,  11.94000,  11.94000,  11.94000, &
		11.94000,  14.92500,  14.92500,  14.92500,  14.92500, &
		14.92500,  14.92500,  14.92500,  14.92500,  14.92500, &
		14.92500,  17.91000,  17.91000,  17.91000,  17.91000, &
		17.91000,  17.91000,  17.91000,  17.91000,  17.91000, &
		17.91000,  20.89500,  20.89500,  20.89500,  20.89500, &
		20.89500,  20.89500,  20.89500,  20.89500,  20.89500, &
		20.89500,  23.88000,  23.88000,  23.88000,  23.88000, &
		23.88000,  23.88000,  23.88000,  23.88000,  23.88000, &
		23.88000,  26.86500,  26.86500,  26.86500,  26.86500, &
		26.86500,  26.86500,  26.86500,  26.86500,  26.86500, &
		26.86500, -29.85000, -29.85000, -29.85000, -29.85000, &
		-29.85000, -29.85000, -29.85000, -29.85000, -29.85000, &
		-29.85000, -26.86500, -26.86500, -26.86500, -26.86500, &
		-26.86500, -26.86500, -26.86500, -26.86500, -26.86500, &
		-26.86500, -23.88000, -23.88000, -23.88000, -23.88000, &
		-23.88000, -23.88000, -23.88000, -23.88000, -23.88000, &
		-23.88000, -20.89500, -20.89500, -20.89500, -20.89500, &
		-20.89500, -20.89500, -20.89500, -20.89500, -20.89500, &
		-20.89500, -17.91000, -17.91000, -17.91000, -17.91000, &
		-17.91000, -17.91000, -17.91000, -17.91000, -17.91000, &
		-17.91000, -14.92500, -14.92500, -14.92500, -14.92500, &
		-14.92500, -14.92500, -14.92500, -14.92500, -14.92500, &
		-14.92500, -11.94000, -11.94000, -11.94000, -11.94000, &
		-11.94000, -11.94000, -11.94000, -11.94000, -11.94000, &
		-11.94000,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-8.95500,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-8.95500,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-5.97000,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-5.97000,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		-2.98500,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		-2.98500,  -1.49250,  -1.49250,  -1.49250,  -1.49250, &
		-1.49250,  -1.49250,  -1.49250,  -1.49250,  -1.49250, &
		-1.49250,   1.49250,   1.49250,   1.49250,   1.49250, &
		1.49250,   1.49250,   1.49250,   1.49250,   1.49250, &
		1.49250,   4.47750,   4.47750,   4.47750,   4.47750, &
		4.47750,   4.47750,   4.47750,   4.47750,   4.47750, &
		4.47750,   7.46250,   7.46250,   7.46250,   7.46250, &
		7.46250,   7.46250,   7.46250,   7.46250,   7.46250, &
		7.46250,  10.44750,  10.44750,  10.44750,  10.44750, &
		10.44750,  10.44750,  10.44750,  10.44750,  10.44750, &
		10.44750,  13.43250,  13.43250,  13.43250,  13.43250, &
		13.43250,  13.43250,  13.43250,  13.43250,  13.43250, &
		13.43250,  16.41750,  16.41750,  16.41750,  16.41750, &
		16.41750,  16.41750,  16.41750,  16.41750,  16.41750, &
		16.41750,  19.40250,  19.40250,  19.40250,  19.40250, &
		19.40250,  19.40250,  19.40250,  19.40250,  19.40250, &
		19.40250,  22.38750,  22.38750,  22.38750,  22.38750, &
		22.38750,  22.38750,  22.38750,  22.38750,  22.38750, &
		22.38750,  25.37250,  25.37250,  25.37250,  25.37250, &
		25.37250,  25.37250,  25.37250,  25.37250,  25.37250, &
		25.37250,  28.35750,  28.35750,  28.35750,  28.35750, &
		28.35750,  28.35750,  28.35750,  28.35750,  28.35750, &
		28.35750, -28.35750, -28.35750, -28.35750, -28.35750, &
		-28.35750, -28.35750, -28.35750, -28.35750, -28.35750, &
		-28.35750, -25.37250, -25.37250, -25.37250, -25.37250, &
		-25.37250, -25.37250, -25.37250, -25.37250, -25.37250, &
		-25.37250, -22.38750, -22.38750, -22.38750, -22.38750, &
		-22.38750, -22.38750, -22.38750, -22.38750, -22.38750, &
		-22.38750, -19.40250, -19.40250, -19.40250, -19.40250, &
		-19.40250, -19.40250, -19.40250, -19.40250, -19.40250, &
		-19.40250, -16.41750, -16.41750, -16.41750, -16.41750, &
		-16.41750, -16.41750, -16.41750, -16.41750, -16.41750, &
		-16.41750, -13.43250, -13.43250, -13.43250, -13.43250, &
		-13.43250, -13.43250, -13.43250, -13.43250, -13.43250, &
		-13.43250, -10.44750, -10.44750, -10.44750, -10.44750, &
		-10.44750, -10.44750, -10.44750, -10.44750, -10.44750, &
		-10.44750,  -7.46250,  -7.46250,  -7.46250,  -7.46250, &
		-7.46250,  -7.46250,  -7.46250,  -7.46250,  -7.46250, &
		-7.46250,  -4.47750,  -4.47750,  -4.47750,  -4.47750, &
		-4.47750,  -4.47750,  -4.47750,  -4.47750,  -4.47750, &
		-4.47750,   1.49250,   1.49250,   1.49250,   1.49250, &
		1.49250,   1.49250,   1.49250,   1.49250,   1.49250, &
		1.49250,   4.47750,   4.47750,   4.47750,   4.47750, &
		4.47750,   4.47750,   4.47750,   4.47750,   4.47750, &
		4.47750,   7.46250,   7.46250,   7.46250,   7.46250, &
		7.46250,   7.46250,   7.46250,   7.46250,   7.46250, &
		7.46250,  10.44750,  10.44750,  10.44750,  10.44750, &
		10.44750,  10.44750,  10.44750,  10.44750,  10.44750, &
		10.44750,  13.43250,  13.43250,  13.43250,  13.43250, &
		13.43250,  13.43250,  13.43250,  13.43250,  13.43250, &
		13.43250,  16.41750,  16.41750,  16.41750,  16.41750, &
		16.41750,  16.41750,  16.41750,  16.41750,  16.41750, &
		16.41750,  19.40250,  19.40250,  19.40250,  19.40250, &
		19.40250,  19.40250,  19.40250,  19.40250,  19.40250, &
		19.40250,  22.38750,  22.38750,  22.38750,  22.38750, &
		22.38750,  22.38750,  22.38750,  22.38750,  22.38750, &
		22.38750,  25.37250,  25.37250,  25.37250,  25.37250, &
		25.37250,  25.37250,  25.37250,  25.37250,  25.37250, &
		25.37250,  28.35750,  28.35750,  28.35750,  28.35750, &
		28.35750,  28.35750,  28.35750,  28.35750,  28.35750, &
		28.35750, -28.35750, -28.35750, -28.35750, -28.35750, &
		-28.35750, -28.35750, -28.35750, -28.35750, -28.35750, &
		-28.35750, -25.37250, -25.37250, -25.37250, -25.37250, &
		-25.37250, -25.37250, -25.37250, -25.37250, -25.37250, &
		-25.37250, -22.38750, -22.38750, -22.38750, -22.38750, &
		-22.38750, -22.38750, -22.38750, -22.38750, -22.38750, &
		-22.38750, -19.40250, -19.40250, -19.40250, -19.40250, &
		-19.40250, -19.40250, -19.40250, -19.40250, -19.40250, &
		-19.40250, -16.41750, -16.41750, -16.41750, -16.41750, &
		-16.41750, -16.41750, -16.41750, -16.41750, -16.41750, &
		-16.41750, -13.43250, -13.43250, -13.43250, -13.43250, &
		-13.43250, -13.43250, -13.43250, -13.43250, -13.43250, &
		-13.43250, -10.44750, -10.44750, -10.44750, -10.44750, &
		-10.44750, -10.44750, -10.44750, -10.44750, -10.44750, &
		-10.44750,  -7.46250,  -7.46250,  -7.46250,  -7.46250, &
		-7.46250,  -7.46250,  -7.46250,  -7.46250,  -7.46250, &
		-7.46250,  -4.47750,  -4.47750,  -4.47750,  -4.47750, &
		-4.47750,  -4.47750,  -4.47750,  -4.47750,  -4.47750, &
		-4.47750,  -1.49250,  -1.49250,  -1.49250,  -1.49250, &
		-1.49250,  -1.49250,  -1.49250,  -1.49250,  -1.49250, &
		-1.49250,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   2.98500,   2.98500,   2.98500,   2.98500, &
		2.98500,   2.98500,   2.98500,   2.98500,   2.98500, &
		2.98500,   5.97000,   5.97000,   5.97000,   5.97000, &
		5.97000,   5.97000,   5.97000,   5.97000,   5.97000, &
		5.97000,   8.95500,   8.95500,   8.95500,   8.95500, &
		8.95500,   8.95500,   8.95500,   8.95500,   8.95500, &
		8.95500,  11.94000,  11.94000,  11.94000,  11.94000, &
		11.94000,  11.94000,  11.94000,  11.94000,  11.94000, &
		11.94000,  14.92500,  14.92500,  14.92500,  14.92500, &
		14.92500,  14.92500,  14.92500,  14.92500,  14.92500, &
		14.92500,  17.91000,  17.91000,  17.91000,  17.91000, &
		17.91000,  17.91000,  17.91000,  17.91000,  17.91000, &
		17.91000,  20.89500,  20.89500,  20.89500,  20.89500, &
		20.89500,  20.89500,  20.89500,  20.89500,  20.89500, &
		20.89500,  23.88000,  23.88000,  23.88000,  23.88000, &
		23.88000,  23.88000,  23.88000,  23.88000,  23.88000, &
		23.88000,  26.86500,  26.86500,  26.86500,  26.86500, &
		26.86500,  26.86500,  26.86500,  26.86500,  26.86500, &
		26.86500, -29.85000, -29.85000, -29.85000, -29.85000, &
		-29.85000, -29.85000, -29.85000, -29.85000, -29.85000, &
		-29.85000, -26.86500, -26.86500, -26.86500, -26.86500, &
		-26.86500, -26.86500, -26.86500, -26.86500, -26.86500, &
		-26.86500, -23.88000, -23.88000, -23.88000, -23.88000, &
		-23.88000, -23.88000, -23.88000, -23.88000, -23.88000, &
		-23.88000, -20.89500, -20.89500, -20.89500, -20.89500, &
		-20.89500, -20.89500, -20.89500, -20.89500, -20.89500, &
		-20.89500, -17.91000, -17.91000, -17.91000, -17.91000, &
		-17.91000, -17.91000, -17.91000, -17.91000, -17.91000, &
		-17.91000, -14.92500, -14.92500, -14.92500, -14.92500, &
		-14.92500, -14.92500, -14.92500, -14.92500, -14.92500, &
		-14.92500, -11.94000, -11.94000, -11.94000, -11.94000, &
		-11.94000, -11.94000, -11.94000, -11.94000, -11.94000, &
		-11.94000,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-8.95500,  -8.95500,  -8.95500,  -8.95500,  -8.95500, &
		-8.95500,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-5.97000,  -5.97000,  -5.97000,  -5.97000,  -5.97000, &
		-5.97000,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		-2.98500,  -2.98500,  -2.98500,  -2.98500,  -2.98500, &
		-2.98500 &
		]

		yO = [ &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		1.25772,   7.75321,  14.24869,  20.74417,  27.23965, &
		-31.21969, -24.72421, -18.22873, -11.73324,  -5.23776, &
		1.25772,   7.75321,  14.24869,  20.74417,  27.23965, &
		-31.21969, -24.72421, -18.22873, -11.73324,  -5.23776, &
		1.25772,   7.75321,  14.24869,  20.74417,  27.23965, &
		-31.21969, -24.72421, -18.22873, -11.73324,  -5.23776, &
		1.25772,   7.75321,  14.24869,  20.74417,  27.23965, &
		-31.21969, -24.72421, -18.22873, -11.73324,  -5.23776, &
		1.25772,   7.75321,  14.24869,  20.74417,  27.23965, &
		-31.21969, -24.72421, -18.22873, -11.73324,  -5.23776, &
		1.25772,   7.75321,  14.24869,  20.74417,  27.23965, &
		-31.21969, -24.72421, -18.22873, -11.73324,  -5.23776, &
		1.25772,   7.75321,  14.24869,  20.74417,  27.23965, &
		-31.21969, -24.72421, -18.22873, -11.73324,  -5.23776, &
		1.25772,   7.75321,  14.24869,  20.74417,  27.23965, &
		-31.21969, -24.72421, -18.22873, -11.73324,  -5.23776, &
		1.25772,   7.75321,  14.24869,  20.74417,  27.23965, &
		-31.21969, -24.72421, -18.22873, -11.73324,  -5.23776, &
		1.25772,   7.75321,  14.24869,  20.74417,  27.23965, &
		-31.21969, -24.72421, -18.22873, -11.73324,  -5.23776, &
		1.25772,   7.75321,  14.24869,  20.74417,  27.23965, &
		-31.21969, -24.72421, -18.22873, -11.73324,  -5.23776, &
		1.25772,   7.75321,  14.24869,  20.74417,  27.23965, &
		-31.21969, -24.72421, -18.22873, -11.73324,  -5.23776, &
		1.25772,   7.75321,  14.24869,  20.74417,  27.23965, &
		-31.21969, -24.72421, -18.22873, -11.73324,  -5.23776, &
		1.25772,   7.75321,  14.24869,  20.74417,  27.23965, &
		-31.21969, -24.72421, -18.22873, -11.73324,  -5.23776, &
		1.25772,   7.75321,  14.24869,  20.74417,  27.23965, &
		-31.21969, -24.72421, -18.22873, -11.73324,  -5.23776, &
		1.25772,   7.75321,  14.24869,  20.74417,  27.23965, &
		-31.21969, -24.72421, -18.22873, -11.73324,  -5.23776, &
		1.25772,   7.75321,  14.24869,  20.74417,  27.23965, &
		-31.21969, -24.72421, -18.22873, -11.73324,  -5.23776, &
		1.25772,   7.75321,  14.24869,  20.74417,  27.23965, &
		-31.21969, -24.72421, -18.22873, -11.73324,  -5.23776, &
		1.25772,   7.75321,  14.24869,  20.74417,  27.23965, &
		-31.21969, -24.72421, -18.22873, -11.73324,  -5.23776, &
		1.25772,   7.75321,  14.24869,  20.74417,  27.23965, &
		-31.21969, -24.72421, -18.22873, -11.73324,  -5.23776, &
		-1.25772,   5.23776,  11.73324,  18.22873,  24.72421, &
		31.21969, -27.23965, -20.74417, -14.24869,  -7.75321, &
		-1.25772,   5.23776,  11.73324,  18.22873,  24.72421, &
		31.21969, -27.23965, -20.74417, -14.24869,  -7.75321, &
		-1.25772,   5.23776,  11.73324,  18.22873,  24.72421, &
		31.21969, -27.23965, -20.74417, -14.24869,  -7.75321, &
		-1.25772,   5.23776,  11.73324,  18.22873,  24.72421, &
		31.21969, -27.23965, -20.74417, -14.24869,  -7.75321, &
		-1.25772,   5.23776,  11.73324,  18.22873,  24.72421, &
		31.21969, -27.23965, -20.74417, -14.24869,  -7.75321, &
		-1.25772,   5.23776,  11.73324,  18.22873,  24.72421, &
		31.21969, -27.23965, -20.74417, -14.24869,  -7.75321, &
		-1.25772,   5.23776,  11.73324,  18.22873,  24.72421, &
		31.21969, -27.23965, -20.74417, -14.24869,  -7.75321, &
		-1.25772,   5.23776,  11.73324,  18.22873,  24.72421, &
		31.21969, -27.23965, -20.74417, -14.24869,  -7.75321, &
		-1.25772,   5.23776,  11.73324,  18.22873,  24.72421, &
		31.21969, -27.23965, -20.74417, -14.24869,  -7.75321, &
		-1.25772,   5.23776,  11.73324,  18.22873,  24.72421, &
		31.21969, -27.23965, -20.74417, -14.24869,  -7.75321, &
		-1.25772,   5.23776,  11.73324,  18.22873,  24.72421, &
		31.21969, -27.23965, -20.74417, -14.24869,  -7.75321, &
		-1.25772,   5.23776,  11.73324,  18.22873,  24.72421, &
		31.21969, -27.23965, -20.74417, -14.24869,  -7.75321, &
		-1.25772,   5.23776,  11.73324,  18.22873,  24.72421, &
		31.21969, -27.23965, -20.74417, -14.24869,  -7.75321, &
		-1.25772,   5.23776,  11.73324,  18.22873,  24.72421, &
		31.21969, -27.23965, -20.74417, -14.24869,  -7.75321, &
		-1.25772,   5.23776,  11.73324,  18.22873,  24.72421, &
		31.21969, -27.23965, -20.74417, -14.24869,  -7.75321, &
		-1.25772,   5.23776,  11.73324,  18.22873,  24.72421, &
		31.21969, -27.23965, -20.74417, -14.24869,  -7.75321, &
		-1.25772,   5.23776,  11.73324,  18.22873,  24.72421, &
		31.21969, -27.23965, -20.74417, -14.24869,  -7.75321, &
		-1.25772,   5.23776,  11.73324,  18.22873,  24.72421, &
		31.21969, -27.23965, -20.74417, -14.24869,  -7.75321, &
		-1.25772,   5.23776,  11.73324,  18.22873,  24.72421, &
		31.21969, -27.23965, -20.74417, -14.24869,  -7.75321, &
		-1.25772,   5.23776,  11.73324,  18.22873,  24.72421, &
		31.21969, -27.23965, -20.74417, -14.24869,  -7.75321, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-0.00000,   6.49548,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,  -6.49548, &
		-2.03002,   4.46547,  10.96095,  17.45643,  23.95191, &
		30.44740, -28.01195, -21.51646, -15.02098,  -8.52550, &
		-2.03002,   4.46547,  10.96095,  17.45643,  23.95191, &
		30.44740, -28.01195, -21.51646, -15.02098,  -8.52550, &
		-2.03002,   4.46547,  10.96095,  17.45643,  23.95191, &
		30.44740, -28.01195, -21.51646, -15.02098,  -8.52550, &
		-2.03002,   4.46547,  10.96095,  17.45643,  23.95191, &
		30.44740, -28.01195, -21.51646, -15.02098,  -8.52550, &
		-2.03002,   4.46547,  10.96095,  17.45643,  23.95191, &
		30.44740, -28.01195, -21.51646, -15.02098,  -8.52550, &
		-2.03002,   4.46547,  10.96095,  17.45643,  23.95191, &
		30.44740, -28.01195, -21.51646, -15.02098,  -8.52550, &
		-2.03002,   4.46547,  10.96095,  17.45643,  23.95191, &
		30.44740, -28.01195, -21.51646, -15.02098,  -8.52550, &
		-2.03002,   4.46547,  10.96095,  17.45643,  23.95191, &
		30.44740, -28.01195, -21.51646, -15.02098,  -8.52550, &
		-2.03002,   4.46547,  10.96095,  17.45643,  23.95191, &
		30.44740, -28.01195, -21.51646, -15.02098,  -8.52550, &
		-2.03002,   4.46547,  10.96095,  17.45643,  23.95191, &
		30.44740, -28.01195, -21.51646, -15.02098,  -8.52550, &
		-2.03002,   4.46547,  10.96095,  17.45643,  23.95191, &
		30.44740, -28.01195, -21.51646, -15.02098,  -8.52550, &
		-2.03002,   4.46547,  10.96095,  17.45643,  23.95191, &
		30.44740, -28.01195, -21.51646, -15.02098,  -8.52550, &
		-2.03002,   4.46547,  10.96095,  17.45643,  23.95191, &
		30.44740, -28.01195, -21.51646, -15.02098,  -8.52550, &
		-2.03002,   4.46547,  10.96095,  17.45643,  23.95191, &
		30.44740, -28.01195, -21.51646, -15.02098,  -8.52550, &
		-2.03002,   4.46547,  10.96095,  17.45643,  23.95191, &
		30.44740, -28.01195, -21.51646, -15.02098,  -8.52550, &
		-2.03002,   4.46547,  10.96095,  17.45643,  23.95191, &
		30.44740, -28.01195, -21.51646, -15.02098,  -8.52550, &
		-2.03002,   4.46547,  10.96095,  17.45643,  23.95191, &
		30.44740, -28.01195, -21.51646, -15.02098,  -8.52550, &
		-2.03002,   4.46547,  10.96095,  17.45643,  23.95191, &
		30.44740, -28.01195, -21.51646, -15.02098,  -8.52550, &
		-2.03002,   4.46547,  10.96095,  17.45643,  23.95191, &
		30.44740, -28.01195, -21.51646, -15.02098,  -8.52550, &
		-2.03002,   4.46547,  10.96095,  17.45643,  23.95191, &
		30.44740, -28.01195, -21.51646, -15.02098,  -8.52550, &
		2.03002,  15.02098,  21.51646,  28.01195, -30.44740, &
		-23.95191, -17.45643,   2.03002,   8.52550,  15.02098, &
		21.51646,  28.01195, -30.44740,   2.03002,   8.52550, &
		15.02098,  21.51646,  28.01195, -30.44740,   2.03002, &
		8.52550,  15.02098,  21.51646, -10.96095,  -4.46547, &
		2.03002,   8.52550,  15.02098,  21.51646, -10.96095, &
		-4.46547,   2.03002,   8.52550, -23.95191, -17.45643, &
		-10.96095,  -4.46547,   2.03002,   8.52550, -23.95191, &
		-17.45643, -10.96095,  -4.46547,  28.01195, -30.44740, &
		-23.95191, -17.45643, -10.96095,  -4.46547,  28.01195, &
		-30.44740, -23.95191, -17.45643,  15.02098,  21.51646, &
		28.01195, -30.44740, -23.95191, -17.45643,  15.02098, &
		21.51646,  28.01195, -30.44740,   2.03002,   8.52550, &
		15.02098,  21.51646,  28.01195, -30.44740,   2.03002, &
		8.52550,  15.02098,  21.51646, -10.96095,  -4.46547, &
		2.03002,   8.52550,  15.02098,  21.51646, -10.96095, &
		-4.46547,   2.03002,   8.52550, -23.95191, -17.45643, &
		-10.96095,  -4.46547,   2.03002,   8.52550,  28.01195, &
		-23.95191, -17.45643, -10.96095,  -4.46547,   2.03002, &
		28.01195, -30.44740, -23.95191, -17.45643, -10.96095, &
		-4.46547,  15.02098,  28.01195, -30.44740, -23.95191, &
		-17.45643, -10.96095,  15.02098,  21.51646,  28.01195, &
		-30.44740, -23.95191, -17.45643,   2.03002,  15.02098, &
		21.51646,  28.01195, -30.44740, -23.95191,   0.00000, &
		6.49548,  12.99097,  19.48645,  25.98193, -32.47741, &
		-12.99097,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -12.99097,  -6.49548,   0.00000,   6.49548, &
		12.99097,  19.48645,  25.98193, -25.98193, -12.99097, &
		-6.49548,   0.00000,  12.99097, -25.98193, -19.48645, &
		-12.99097,  -6.49548,   0.00000,   6.49548,  12.99097, &
		25.98193, -25.98193, -19.48645, -12.99097,   0.00000, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,  12.99097,  25.98193, -32.47741, &
		-25.98193, -12.99097,  12.99097,  19.48645,  25.98193, &
		-32.47741, -25.98193, -19.48645, -12.99097,   0.00000, &
		12.99097,  19.48645,  25.98193, -25.98193, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   1.26772,   7.76321,  14.25869,  20.75417, &
		27.24965, -31.20969, -24.71421, -18.21873, -11.72324, &
		-5.22776,   1.26772,   7.76321,  14.25869,  20.75417, &
		27.24965, -31.20969, -24.71421, -18.21873, -11.72324, &
		-5.22776,   1.26772,   7.76321,  14.25869,  20.75417, &
		27.24965, -31.20969, -24.71421, -18.21873, -11.72324, &
		-5.22776,   1.26772,   7.76321,  14.25869,  20.75417, &
		27.24965, -31.20969, -24.71421, -18.21873, -11.72324, &
		-5.22776,   1.26772,   7.76321,  14.25869,  20.75417, &
		27.24965, -31.20969, -24.71421, -18.21873, -11.72324, &
		-5.22776,   1.26772,   7.76321,  14.25869,  20.75417, &
		27.24965, -31.20969, -24.71421, -18.21873, -11.72324, &
		-5.22776,   1.26772,   7.76321,  14.25869,  20.75417, &
		27.24965, -31.20969, -24.71421, -18.21873, -11.72324, &
		-5.22776,   1.26772,   7.76321,  14.25869,  20.75417, &
		27.24965, -31.20969, -24.71421, -18.21873, -11.72324, &
		-5.22776,   1.26772,   7.76321,  14.25869,  20.75417, &
		27.24965, -31.20969, -24.71421, -18.21873, -11.72324, &
		-5.22776,   1.26772,   7.76321,  14.25869,  20.75417, &
		27.24965, -31.20969, -24.71421, -18.21873, -11.72324, &
		-5.22776,   1.26772,   7.76321,  14.25869,  20.75417, &
		27.24965, -31.20969, -24.71421, -18.21873, -11.72324, &
		-5.22776,   1.26772,   7.76321,  14.25869,  20.75417, &
		27.24965, -31.20969, -24.71421, -18.21873, -11.72324, &
		-5.22776,   1.26772,   7.76321,  14.25869,  20.75417, &
		27.24965, -31.20969, -24.71421, -18.21873, -11.72324, &
		-5.22776,   1.26772,   7.76321,  14.25869,  20.75417, &
		27.24965, -31.20969, -24.71421, -18.21873, -11.72324, &
		-5.22776,   1.26772,   7.76321,  14.25869,  20.75417, &
		27.24965, -31.20969, -24.71421, -18.21873, -11.72324, &
		-5.22776,   1.26772,   7.76321, -11.72324,  -5.22776, &
		1.26772,   7.76321,  14.25869,  20.75417,  27.24965, &
		-31.20969, -24.71421, -18.21873, -11.72324,  -5.22776, &
		1.26772,   7.76321,  14.25869,  20.75417,  27.24965, &
		-31.20969, -24.71421, -18.21873, -11.72324,  -5.22776, &
		1.26772,   7.76321,  14.25869,  20.75417,  27.24965, &
		-31.20969, -24.71421, -18.21873, -11.72324,  -5.22776, &
		-14.25869,  -7.76321,  -1.26772,   5.22776,  11.72324, &
		18.21873,  24.71421,  31.20969, -27.24965, -20.75417, &
		-14.25869,  -7.76321,  -1.26772,   5.22776,  11.72324, &
		18.21873,  24.71421,  31.20969, -27.24965, -20.75417, &
		-14.25869,  -7.76321,  -1.26772,   5.22776,  11.72324, &
		18.21873,  24.71421,  31.20969, -27.24965, -20.75417, &
		-14.25869,  -7.76321,  -1.26772,   5.22776,  11.72324, &
		18.21873,  24.71421,  31.20969, -27.24965, -20.75417, &
		-27.24965, -20.75417, -14.25869,  -7.76321,  -1.26772, &
		5.22776,  11.72324,  18.21873,  24.71421,  31.20969, &
		-27.24965, -20.75417, -14.25869,  -7.76321,  -1.26772, &
		5.22776,  11.72324,  18.21873,  24.71421,  31.20969, &
		-27.24965, -20.75417, -14.25869,  -7.76321,  -1.26772, &
		5.22776,  11.72324,  18.21873,  24.71421,  31.20969, &
		-27.24965, -20.75417, -14.25869,  -7.76321,  -1.26772, &
		5.22776,  11.72324,  18.21873,  24.71421,  31.20969, &
		24.71421,  31.20969, -27.24965, -20.75417, -14.25869, &
		-7.76321,  -1.26772,   5.22776,  11.72324,  18.21873, &
		24.71421,  31.20969, -27.24965, -20.75417, -14.25869, &
		-7.76321,  -1.26772,   5.22776,  11.72324,  18.21873, &
		24.71421,  31.20969, -27.24965, -20.75417, -14.25869, &
		-7.76321,  -1.26772,   5.22776,  11.72324,  18.21873, &
		24.71421,  31.20969, -27.24965, -20.75417, -14.25869, &
		-7.76321,  -1.26772,   5.22776,  11.72324,  18.21873, &
		11.72324,  18.21873,  24.71421,  31.20969, -27.24965, &
		-20.75417, -14.25869,  -7.76321,  -1.26772,   5.22776, &
		11.72324,  18.21873,  24.71421,  31.20969, -27.24965, &
		-20.75417, -14.25869,  -7.76321,  -1.26772,   5.22776, &
		11.72324,  18.21873,  24.71421,  31.20969, -27.24965, &
		-20.75417, -14.25869,  -7.76321,  -1.26772,   5.22776, &
		11.72324,  18.21873,  24.71421,  31.20969, -27.24965, &
		-20.75417, -14.25869,  -7.76321,  -1.26772,   5.22776, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-3.24774,   3.24774,   9.74322,  16.23871,  22.73419, &
		29.22967, -29.22967, -22.73419, -16.23871,  -9.74322, &
		-16.23871,  -9.74322,  -3.24774,   3.24774,   9.74322, &
		16.23871,  22.73419,  29.22967, -29.22967, -22.73419, &
		-16.23871,  -9.74322,  -3.24774,   3.24774,   9.74322, &
		16.23871,  22.73419,  29.22967, -29.22967, -22.73419, &
		-16.23871,  -9.74322,  -3.24774,   3.24774,   9.74322, &
		16.23871,  22.73419,  29.22967, -29.22967, -22.73419, &
		-16.23871,  -9.74322,  -3.24774,   3.24774,   9.74322, &
		16.23871,  22.73419,  29.22967, -29.22967, -22.73419, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		-29.22967, -22.73419, -16.23871,  -9.74322,  -3.24774, &
		3.24774,   9.74322,  16.23871,  22.73419,  29.22967, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,  -2.03002,   4.46547,  10.96095,  17.45643, &
		23.95191,  30.44740, -28.01195, -21.51646, -15.02098, &
		-8.52550,  -2.03002,   4.46547,  10.96095,  17.45643, &
		23.95191,  30.44740, -28.01195, -21.51646, -15.02098, &
		-8.52550,  -2.03002,   4.46547,  10.96095,  17.45643, &
		23.95191,  30.44740, -28.01195, -21.51646, -15.02098, &
		-8.52550,  -2.03002,   4.46547,  10.96095,  17.45643, &
		23.95191,  30.44740, -28.01195, -21.51646, -15.02098, &
		-8.52550,  -2.03002,   4.46547,  10.96095,  17.45643, &
		23.95191,  30.44740, -28.01195, -21.51646, -15.02098, &
		-8.52550,  -2.03002,   4.46547,  10.96095,  17.45643, &
		23.95191,  30.44740, -28.01195, -21.51646, -15.02098, &
		-8.52550,  -2.03002,   4.46547,  10.96095,  17.45643, &
		23.95191,  30.44740, -28.01195, -21.51646, -15.02098, &
		-8.52550,  -2.03002,   4.46547,  10.96095,  17.45643, &
		23.95191,  30.44740, -28.01195, -21.51646, -15.02098, &
		-8.52550,  -2.03002,   4.46547,  10.96095,  17.45643, &
		23.95191,  30.44740, -28.01195, -21.51646, -15.02098, &
		-8.52550,  -2.03002,   4.46547,  10.96095,  17.45643, &
		23.95191,  30.44740, -28.01195, -21.51646, -15.02098, &
		-8.52550,  -2.03002,   4.46547,  10.96095,  17.45643, &
		23.95191,  30.44740, -28.01195, -21.51646, -15.02098, &
		-8.52550,  -2.03002,   4.46547,  10.96095,  17.45643, &
		23.95191,  30.44740, -28.01195, -21.51646, -15.02098, &
		-8.52550,  -2.03002,   4.46547,  10.96095,  17.45643, &
		23.95191,  30.44740, -28.01195, -21.51646, -15.02098, &
		-8.52550,  -2.03002,   4.46547,  10.96095,  17.45643, &
		23.95191,  30.44740, -28.01195, -21.51646, -15.02098, &
		-8.52550,  -2.03002,   4.46547,  10.96095,  17.45643, &
		23.95191,  30.44740, -28.01195, -21.51646, -15.02098, &
		-8.52550,  -2.03002,   4.46547,  10.96095,  17.45643, &
		23.95191,  30.44740, -28.01195, -21.51646, -15.02098, &
		-8.52550,  -2.03002,   4.46547,  10.96095,  17.45643, &
		23.95191,  30.44740, -28.01195, -21.51646, -15.02098, &
		-8.52550,  -2.03002,   4.46547,  10.96095,  17.45643, &
		23.95191,  30.44740, -28.01195, -21.51646, -15.02098, &
		-8.52550,  -2.03002,   4.46547,  10.96095,  17.45643, &
		23.95191,  30.44740, -28.01195, -21.51646, -15.02098, &
		-8.52550,  -2.03002,   4.46547,  10.96095,  17.45643, &
		23.95191,  30.44740, -28.01195, -21.51646, -15.02098, &
		-8.52550,   2.03002,   8.52550,  15.02098,  21.51646, &
		28.01195, -30.44740, -23.95191, -17.45643, -10.96095, &
		-4.46547,   2.03002,   8.52550,  15.02098,  21.51646, &
		28.01195, -30.44740, -23.95191, -17.45643, -10.96095, &
		-4.46547,   2.03002,   8.52550,  15.02098,  21.51646, &
		28.01195, -30.44740, -23.95191, -17.45643, -10.96095, &
		-4.46547,   2.03002,   8.52550,  15.02098,  21.51646, &
		28.01195, -30.44740, -23.95191, -17.45643, -10.96095, &
		-4.46547,   2.03002,   8.52550,  15.02098,  21.51646, &
		28.01195, -30.44740, -23.95191, -17.45643, -10.96095, &
		-4.46547,   2.03002,   8.52550,  15.02098,  21.51646, &
		28.01195, -30.44740, -23.95191, -17.45643, -10.96095, &
		-4.46547,   2.03002,   8.52550,  15.02098,  21.51646, &
		28.01195, -30.44740, -23.95191, -17.45643, -10.96095, &
		-4.46547,   2.03002,   8.52550,  15.02098,  21.51646, &
		28.01195, -30.44740, -23.95191, -17.45643, -10.96095, &
		-4.46547,   2.03002,   8.52550,  15.02098,  21.51646, &
		28.01195, -30.44740, -23.95191, -17.45643, -10.96095, &
		-4.46547,   2.03002,   8.52550,  15.02098,  21.51646, &
		28.01195, -30.44740, -23.95191, -17.45643, -10.96095, &
		-4.46547,   2.03002,   8.52550,  15.02098,  21.51646, &
		28.01195, -30.44740, -23.95191, -17.45643, -10.96095, &
		-4.46547,   2.03002,   8.52550,  15.02098,  21.51646, &
		28.01195, -30.44740, -23.95191, -17.45643, -10.96095, &
		-4.46547,   2.03002,   8.52550,  15.02098,  21.51646, &
		28.01195, -30.44740, -23.95191, -17.45643, -10.96095, &
		-4.46547,   2.03002,   8.52550,  15.02098,  21.51646, &
		28.01195, -30.44740, -23.95191, -17.45643, -10.96095, &
		-4.46547,   2.03002,   8.52550,  15.02098,  21.51646, &
		28.01195, -30.44740, -23.95191, -17.45643, -10.96095, &
		-4.46547,   2.03002,   8.52550,  15.02098,  21.51646, &
		28.01195, -30.44740, -23.95191, -17.45643, -10.96095, &
		-4.46547,   2.03002,   8.52550,  15.02098,  21.51646, &
		28.01195, -30.44740, -23.95191, -17.45643, -10.96095, &
		-4.46547,   2.03002,   8.52550,  15.02098,  21.51646, &
		28.01195, -30.44740, -23.95191, -17.45643, -10.96095, &
		-4.46547,   2.03002,   8.52550,  15.02098,  21.51646, &
		28.01195, -30.44740, -23.95191, -17.45643, -10.96095, &
		-4.46547,   2.03002,   8.52550,  15.02098,  21.51646, &
		28.01195, -30.44740, -23.95191, -17.45643, -10.96095, &
		-4.46547,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,   0.00000,   6.49548,  12.99097,  19.48645, &
		25.98193, -32.47741, -25.98193, -19.48645, -12.99097, &
		-6.49548,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,   1.25772,   7.75321,  14.24869,  20.74417, &
		27.23965, -31.21969, -24.72421, -18.22873, -11.73324, &
		-5.23776,   1.25772,   7.75321,  14.24869,  20.74417, &
		27.23965, -31.21969, -24.72421, -18.22873, -11.73324, &
		-5.23776,   1.25772,   7.75321,  14.24869,  20.74417, &
		27.23965, -31.21969, -24.72421, -18.22873, -11.73324, &
		-5.23776,   1.25772,   7.75321,  14.24869,  20.74417, &
		27.23965, -31.21969, -24.72421, -18.22873, -11.73324, &
		-5.23776,   1.25772,   7.75321,  14.24869,  20.74417, &
		27.23965, -31.21969, -24.72421, -18.22873, -11.73324, &
		-5.23776,   1.25772,   7.75321,  14.24869,  20.74417, &
		27.23965, -31.21969, -24.72421, -18.22873, -11.73324, &
		-5.23776,   1.25772,   7.75321,  14.24869,  20.74417, &
		27.23965, -31.21969, -24.72421, -18.22873, -11.73324, &
		-5.23776,   1.25772,   7.75321,  14.24869,  20.74417, &
		27.23965, -31.21969, -24.72421, -18.22873, -11.73324, &
		-5.23776,   1.25772,   7.75321,  14.24869,  20.74417, &
		27.23965, -31.21969, -24.72421, -18.22873, -11.73324, &
		-5.23776,   1.25772,   7.75321,  14.24869,  20.74417, &
		27.23965, -31.21969, -24.72421, -18.22873, -11.73324, &
		-5.23776,   1.25772,   7.75321,  14.24869,  20.74417, &
		27.23965, -31.21969, -24.72421, -18.22873, -11.73324, &
		-5.23776,   1.25772,   7.75321,  14.24869,  20.74417, &
		27.23965, -31.21969, -24.72421, -18.22873, -11.73324, &
		-5.23776,   1.25772,   7.75321,  14.24869,  20.74417, &
		27.23965, -31.21969, -24.72421, -18.22873, -11.73324, &
		-5.23776,   1.25772,   7.75321,  14.24869,  20.74417, &
		27.23965, -31.21969, -24.72421, -18.22873, -11.73324, &
		-5.23776,   1.25772,   7.75321,  14.24869,  20.74417, &
		27.23965, -31.21969, -24.72421, -18.22873, -11.73324, &
		-5.23776,   1.25772,   7.75321,  14.24869,  20.74417, &
		27.23965, -31.21969, -24.72421, -18.22873, -11.73324, &
		-5.23776,   1.25772,   7.75321,  14.24869,  20.74417, &
		27.23965, -31.21969, -24.72421, -18.22873, -11.73324, &
		-5.23776,   1.25772,   7.75321,  14.24869,  20.74417, &
		27.23965, -31.21969, -24.72421, -18.22873, -11.73324, &
		-5.23776,   1.25772,   7.75321,  14.24869,  20.74417, &
		27.23965, -31.21969, -24.72421, -18.22873, -11.73324, &
		-5.23776,   1.25772,   7.75321,  14.24869,  20.74417, &
		27.23965, -31.21969, -24.72421, -18.22873, -11.73324, &
		-5.23776,  -1.25772,   5.23776,  11.73324,  18.22873, &
		24.72421,  31.21969, -27.23965, -20.74417, -14.24869, &
		-7.75321,  -1.25772,   5.23776,  11.73324,  18.22873, &
		24.72421,  31.21969, -27.23965, -20.74417, -14.24869, &
		-7.75321,  -1.25772,   5.23776,  11.73324,  18.22873, &
		24.72421,  31.21969, -27.23965, -20.74417, -14.24869, &
		-7.75321,  -1.25772,   5.23776,  11.73324,  18.22873, &
		24.72421,  31.21969, -27.23965, -20.74417, -14.24869, &
		-7.75321,  -1.25772,   5.23776,  11.73324,  18.22873, &
		24.72421,  31.21969, -27.23965, -20.74417, -14.24869, &
		-7.75321,  -1.25772,   5.23776,  11.73324,  18.22873, &
		24.72421,  31.21969, -27.23965, -20.74417, -14.24869, &
		-7.75321,  -1.25772,   5.23776,  11.73324,  18.22873, &
		24.72421,  31.21969, -27.23965, -20.74417, -14.24869, &
		-7.75321,  -1.25772,   5.23776,  11.73324,  18.22873, &
		24.72421,  31.21969, -27.23965, -20.74417, -14.24869, &
		-7.75321,  -1.25772,   5.23776,  11.73324,  18.22873, &
		24.72421,  31.21969, -27.23965, -20.74417, -14.24869, &
		-7.75321,  -1.25772,   5.23776,  11.73324,  18.22873, &
		24.72421,  31.21969, -27.23965, -20.74417, -14.24869, &
		-7.75321,  -1.25772,   5.23776,  11.73324,  18.22873, &
		24.72421,  31.21969, -27.23965, -20.74417, -14.24869, &
		-7.75321,  -1.25772,   5.23776,  11.73324,  18.22873, &
		24.72421,  31.21969, -27.23965, -20.74417, -14.24869, &
		-7.75321,  -1.25772,   5.23776,  11.73324,  18.22873, &
		24.72421,  31.21969, -27.23965, -20.74417, -14.24869, &
		-7.75321,  -1.25772,   5.23776,  11.73324,  18.22873, &
		24.72421,  31.21969, -27.23965, -20.74417, -14.24869, &
		-7.75321,  -1.25772,   5.23776,  11.73324,  18.22873, &
		24.72421,  31.21969, -27.23965, -20.74417, -14.24869, &
		-7.75321,  -1.25772,   5.23776,  11.73324,  18.22873, &
		24.72421,  31.21969, -27.23965, -20.74417, -14.24869, &
		-7.75321,  -1.25772,   5.23776,  11.73324,  18.22873, &
		24.72421,  31.21969, -27.23965, -20.74417, -14.24869, &
		-7.75321,  -1.25772,   5.23776,  11.73324,  18.22873, &
		24.72421,  31.21969, -27.23965, -20.74417, -14.24869, &
		-7.75321,  -1.25772,   5.23776,  11.73324,  18.22873, &
		24.72421,  31.21969, -27.23965, -20.74417, -14.24869, &
		-7.75321,  -1.25772,   5.23776,  11.73324,  18.22873, &
		24.72421,  31.21969, -27.23965, -20.74417, -14.24869, &
		-7.75321,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322,  -3.24774,   3.24774,   9.74322,  16.23871, &
		22.73419,  29.22967, -29.22967, -22.73419, -16.23871, &
		-9.74322 &
		]

		zO = [ &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		7.86321,   7.86321,   7.86321,   7.86321,   7.86321, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		6.66548,   6.66548,   6.66548,   6.66548,   6.66548, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		5.29776,   5.29776,   5.29776,   5.29776,   5.29776, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		4.51546,   4.51546,   4.51546,   4.51546,   4.51546, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   3.28774, &
		3.28774,   3.28774,   3.28774,   3.28774,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.99002,   1.99002,   1.99002,   1.99002, &
		1.99002,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   1.26772,   1.26772,   1.26772,   1.26772, &
		1.26772,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		0.00000,   0.00000,   0.00000,   0.00000,   0.00000, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.26772,  -1.26772,  -1.26772,  -1.26772, &
		-1.26772,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -1.99002,  -1.99002,  -1.99002,  -1.99002, &
		-1.99002,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -3.28774,  -3.28774,  -3.28774,  -3.28774, &
		-3.28774,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -4.51546,  -4.51546,  -4.51546,  -4.51546, &
		-4.51546,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -5.29776,  -5.29776,  -5.29776,  -5.29776, &
		-5.29776,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -6.66548,  -6.66548,  -6.66548,  -6.66548, &
		-6.66548,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321,  -7.86321,  -7.86321,  -7.86321,  -7.86321, &
		-7.86321 &
		]
		
		! Sets the origin on Ti5f
		zO = zO - 6.38548
		zTi = zTi - 6.38548
		
		! C6 values and vdW radii from
		! S. Grimme, J Comput Chem 27 (2006) 1787-1799
		! Units c6=[J nm^6 mol^{-1}], r=[Angstrom]
		c6He = 0.08_8*unit_J_to_H*unit_nm_to_a0**6/unit_mol_to_au
		c6O = 0.70_8*unit_J_to_H*unit_nm_to_a0**6/unit_mol_to_au
		c6Ti = 10.80_8*unit_J_to_H*unit_nm_to_a0**6/unit_mol_to_au
		
		rvdWHe = 1.012_8*unit_Ang_to_a0
		rvdWO = 1.342_8*unit_Ang_to_a0
		rvdWTi = 1.562_8*unit_Ang_to_a0
		
		d = 20.0_8
		s6_PBE = -0.75_8
		
		c6OHe = sqrt( c6O*c6He )
		c6TiHe = sqrt( c6Ti*c6He )
		rvdWOHe = rvdWO+rvdWHe
		rvdWTiHe = rvdWTi+rvdWHe
		
		EdispGrimme = 0.0_8
		
		! Ti-He contribution
		do i=1, nAtomsTi
			if( zTi(i) > 1.5-6.38548 ) then
			r = sqrt( (xTi(i)-xHe)**2 + (yTi(i)-yHe)**2 + (zTi(i)-zHe)**2 )
			r = r*unit_Ang_to_a0
			
			dampF = 1.0_8/( 1.0_8+exp(-d*(r/rvdWTiHe-1.0)) )
			EdispGrimme = EdispGrimme + s6_PBE*dampF*c6TiHe/r**6
			end if
		end do
		
		! O-He contribution
		do i=1, nAtomsO
			if( zTi(i) > 1.5-6.38548 ) then
			r = sqrt( (xO(i)-xHe)**2 + (yO(i)-yHe)**2 + (zO(i)-zHe)**2 )
			r = r*unit_Ang_to_a0
			
			dampF = 1.0_8/( 1.0_8+exp(-d*(r/rvdWOHe-1.0)) )
			EdispGrimme = EdispGrimme + s6_PBE*dampF*c6OHe/r**6
			end if
		end do
		
		EdispGrimme = EdispGrimme*Htocm1
	end function EdispGrimme
	
	!*
	! @brief Potencial segun el modelo de potencial lateralmente
	!        promediado corregido por c3
	!*
	function V_LAPC3LJv2( x, y, z ) result( output )
		real(8), intent(in) :: x
		real(8), intent(in) :: y
		real(8), intent(in) :: z
		real(8) :: output
		
		output = De_LAPC3LJv2*( exp(-2.0_8*alpha_LAPC3LJv2*(z-Re_LAPC3LJv2)) &
				  - 2.0_8*exp(-alpha_LAPC3LJv2*(z-Re_LAPC3LJv2)) ) &
				  + EdispC3LJ( z, z0_LAPC3LJv2, w_LAPC3LJv2, c3_LAPC3LJv2 )
				  
	end function V_LAPC3LJv2
	
	!*
	! @brief Derivada del potencial segun el modelo de potencial
	!        lateralmente promediado corregido por c3
	!*
	function dV_LAPC3LJv2( x, y, z, coord ) result( output )
		real(8), intent(in) :: x
		real(8), intent(in) :: y
		real(8), intent(in) :: z
		integer, intent(in) :: coord
		real(8) :: output
		
		select case ( coord )
			case( 1 )
				output = 0.0_8
			case( 2 )
				output = 0.0_8
			case( 3 )
				output = 2.0_8*De_LAPC3LJv2*alpha_LAPC3LJv2*( exp(-alpha_LAPC3LJv2*(z-Re_LAPC3LJv2)) &
				          - exp(-2.0_8*alpha_LAPC3LJv2*(z-Re_LAPC3LJv2)) ) &
				          + dEdispC3LJ( z, z0_LAPC3LJv2, w_LAPC3LJv2, c3_LAPC3LJv2 )
		end select
	end function dV_LAPC3LJv2
	
	!*
	! @brief Subrutina encargada de varificar el coorecto
	!        funcionamiento de la clase
	!*
	subroutine HeTiO2Potential_test()
		real(8) :: x, y, z
		real(8) :: zVec(17)
		real(8) :: zPot(5)
		integer :: i, j
		integer :: g1, g2, maxg
		integer :: nx, ny
		real(8) :: sum
		real(8) :: aprox, exact
		real(8), allocatable :: vg(:,:)
		real(8) :: estimated1, estimated2, delta, delta00
		type(HeTiO2Potential) :: potential
		
! 		zVec = [ 9.58548, 10.33548, 11.08548, 11.83548, 12.58548, 13.33548, 14.08548, 14.83548, 15.58548 ]
! 		zVec = zVec - 6.38548
! 		
! 		call potential.init( M3DP, grimme=.true. )
		
! 		zVec = [ 2.48838806, 2.77838806, 2.92338806, 3.06838806, 3.21338806, 3.35838806, 3.50338806, 3.64838806, 3.93838806, 4.22838806, 4.51838806, 4.80838806, 5.09838806, 5.38838806, 5.67838806, 6.54838806, 6.83838806 ]
! 		do i=1,size(zVec)
! 			write(*,"(F5.3,F30.13)") zVec(i), EdispGrimme( 0.0_8, 0.0_8, zVec(i) )
! 		end do
		
! 		call potential.init( LAP )
! 		
! 		z = 2.0_8
! 		do i=0,100
! 			write(*,*) z, potential.V( 0.0_8, 0.0_8, z ), potential.dV( 0.0_8, 0.0_8, z, 3 )
! 			
! 			z = z + 0.1_8
! 		end do
! 		
! 		call potential.destroy()
		
! 		maxg = 3
!  		
! 		allocate( vg(maxg+1,maxg+1) )
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Ejemplo de calculo de lo factores nu_G para el
		! modelo Morse corrugado
! #define zetaxy F(x,y,4)
! 		call potential.init( CMP )
! 		
! 		do g1=0,maxg
! 			do g2=0,maxg
! 				nx = 1000
! 				ny = 1000
! 				
! 				sum = 0.0_8
! 				do x = -aCell/2.0_8,aCell/2.0_8,aCell/real(nx,8)
! 					do y = -bCell/2.0_8,bCell/2.0_8,bCell/real(ny,8)
! 						sum = sum + exp(2.0_8*alpha_CMP*zetaxy)*cos(2.0_8*pi*g1*x/aCell+2.0_8*pi*g2*y/bCell)
! 					end do
! 				end do
! 				
! 				estimated1 =  sum/real(nx*ny,8)
! 				
! 				nx = nx+10
! 				ny = ny+10
! 				
! 				sum = 0.0_8
! 				do x = -aCell/2.0_8,aCell/2.0_8,aCell/real(nx,8)
! 					do y = -bCell/2.0_8,bCell/2.0_8,bCell/real(ny,8)
! 						sum = sum + exp(2.0_8*alpha_CMP*zetaxy)*cos(2.0_8*pi*g1*x/aCell+2.0_8*pi*g2*y/bCell)
! 					end do
! 				end do
! 				
! 				estimated2 =  sum/real(nx*ny,8)
! 				
! 				delta = abs( estimated2-estimated1 )
! 				
! 				vg(g1+1,g2+1) = estimated2
! 				
! 				if( g1==0 .and. g2==0 ) then
! 					delta00 = delta
! 					write(*,"(2I6,F20.8,A4,F10.8,F20.8,A4,F10.8)") &
! 						g1, g2, &
! 						estimated2, "+/-", abs(delta), &
! 						1.00_8, "+/-", 0.00_8
! 				else
! 					write(*,"(2I6,F20.8,A4,F10.8,F20.8,A4,F10.8)") &
! 						g1, g2, &
! 						estimated2, "+/-", abs(delta), &
! 						estimated2/vg(1,1), "+/-", abs(estimated2*sqrt( (delta00/vg(1,1))**2+(delta/estimated2)**2 ))
! 				end if
! 			end do
! 		end do
! 			
! 		write(*,*) ""
! 		write(*,*) ""
! 		
! 		call potential.destroy()
! #undef zetaxy
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Comparacion de los factores V_G calculados
		! numericamente comparados con su resultado
		! analitico
! 		call potential.init( CMP )
! 		
! 		do g1=0,maxg
! 			do g2=0,maxg
! 				i=1
! 				write(*,"(A,3I5)") "#", i, g1, g2
! 				
! 				do z = 2.5,10.0,0.1
! 					sum = 0.0_8
! 					do x = -aCell/2.0_8,aCell/2.0_8,aCell/real(nx,8)
! 						do y = -bCell/2.0_8,bCell/2.0_8,bCell/real(ny,8)
! 							sum = sum + potential.V( x, y, z )*cos(2.0_8*pi*g1*x/aCell+2.0_8*pi*g2*y/bCell)
! 						end do
! 					end do
! 						
! 					aprox = sum/real(nx*ny,8)
! 						
! 					if( g1==0 .and. g2==0 ) then
! 						exact = De_CMP*( exp(-2.0_8*alpha_CMP*(z-Re_CMP))-2.0_8*exp(-alpha_CMP*(z-Re_CMP)) )
! 					else
! 						exact = De_CMP*( vg(g1+1,g2+1)/vg(1,1) )*exp(-2.0_8*alpha_CMP*(z-Re_CMP))
! 					end if
! 					
! 					write(*,"(4F20.7)") z, aprox, exact, abs(aprox-exact)
! 				end do
! 				
! 				write(*,*) ""
! 				write(*,*) ""
! 				i=i+1
! 			end do
! 		end do
! 		
! 		call potential.destroy()
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Plano de potencial para el modelo de Morse
		! corrugado
! 		zPot = [ 3.2_8, 4.0_8, 4.2_8, 4.6_8, 5.0_8 ]
! 		call potential.init( CMP )
! 		
! 		do i=1,5
! 			z = zPot(i)
! 			do x = -aCell,aCell,2.0_8*aCell/100.0_8
! 				do y = -2.0_8*bCell,2.0_8*bCell,4.0_8*bCell/200.0_8
! 					write(*,"(4F20.7)") x, y, potential.V( x, y, z )
! 				end do
! 				write(*,*) ""
! 			end do
! 			
! 			write(*,*) ""
! 			write(*,*) ""
! 		end do
! 		
! 		call potential.destroy()
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Calculo de los factores V_G para el potencial
		! de Morse 3D
! 		call potential.init( CMP )
! 		
! 		nx = 100
! 		ny = 100
! 		
! 		do g1=0,maxg
! 			do g2=0,maxg
! 				i=1
! 				write(*,"(A,3I5)") "#", i, g1, g2
! 				
! 				do z = 2.5,10.0,0.1
! 					sum = 0.0_8
! 					do x = -aCell/2.0_8,aCell/2.0_8,aCell/real(nx,8)
! 						do y = -bCell/2.0_8,bCell/2.0_8,bCell/real(ny,8)
! 							sum = sum + potential.V( x, y, z )*cos(2.0_8*pi*g1*x/aCell+2.0_8*pi*g2*y/bCell)
! 						end do
! 					end do
! 						
! 					aprox = sum/real(nx*ny,8)
! 					
! 					write(*,"(4F20.7)") z, aprox
! 				end do
! 				
! 				write(*,*) ""
! 				write(*,*) ""
! 				i=i+1
! 			end do
! 		end do
! 		
! 		call potential.destroy()

		call potential.init( LAPC3LJv2 )
! 		call potential.init( LAPC3LJ )
! 		call potential.init( LAP )
		
		write(*,"(1A,4A20)") "#", "z", "V(x,y,z)", "dV(x,y,z)/dz", "N[dV(x,y,z)/dz]"
!  		y = bCell/3.0_8
		y = 0.0_8
		z = 2.0_8
		do while( z <= 10.0_8 )
! 		z = -25.0_8
! 		do while( z <= -15.0_8 )
! 			write(*,"(6F20.7)") z, potential.V( 0.0_8, y, z+25.0_8 ), &
! 				potential.dV( 0.0_8, y, z+25.0_8, 3 ), potential.NdV( 0.0_8, y, z+25.0_8, 3 )
			write(*,"(6F20.7)") z, potential.V( 0.0_8, y, z ), &
				potential.dV( 0.0_8, y, z, 3 ), potential.NdV( 0.0_8, y, z, 3 )
			
			z = z + 0.01_8
		end do
		
		call potential.destroy()
		
	end subroutine HeTiO2Potential_test
	
end module HeTiO2Potential_
