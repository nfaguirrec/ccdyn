!!**********************************************************************************
!!    Consejo Superior de Investigaciones Científicas                              !
!!    Departamento de Física Atómica, Molecular y de Agregados                     !
!!    http://www.iff.csic.es/fama/                                                 !
!!                                                                                 !
!!    Authors:                                                                     !
!!    (2011-2015) Néstor F. Aguirre                                                !
!!                nfaguirrec@iff.csic.es                                           !
!!                nfaguirrec@gmail.com                                             !
!!    (2011-2015) María Pilar de Lara-Castells                                     !
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

module AuTiO2Potential_
	implicit none
	private
	
	public :: &
		AuTiO2Potential_test
		
	real(8), private, parameter :: angs = 1.0_8/0.52917726_8
	real(8), private, parameter :: eV = 1.0_8/27.211396132_8
	real(8), private, parameter :: cm1 = 1.0_8/219474.63068_8
		
	integer, public, parameter :: SM = 1
	integer, public, parameter :: SM1Dz = 2
	integer, public, parameter :: SMSP = 3
	integer, public, parameter :: SMSP1Dz = 4
	integer, public, parameter :: SMSPDisp1Dz = 5
	
	real(8), private, parameter :: pi = acos(-1.0_8)
	real(8), private, parameter :: aCell = 2.953_8
	real(8), private, parameter :: bCell = 6.486_8
	
	real(8), private :: aDe
	real(8), private :: bDe
	real(8), private :: cDe
	real(8), private :: gammaDe
	real(8), private :: tauDe
	
	real(8), private :: aAlpha
	real(8), private :: bAlpha
	real(8), private :: cAlpha
	real(8), private :: gammaAlpha
	real(8), private :: tauAlpha
	
	real(8), private :: aBeta
	real(8), private :: bBeta
	real(8), private :: cBeta
	real(8), private :: gammaBeta
	real(8), private :: tauBeta
	
	real(8), private :: aZe
	real(8), private :: bZe
	real(8), private :: cZe
	real(8), private :: gammaZe
	real(8), private :: tauZe
	
	! Parameters for SMSPDisp1Dz model
	real(8), private :: rcut
	real(8), private :: maxV
	real(8), private :: mindV
	real(8), private :: C3
	real(8), private :: C4
	real(8), private :: C6
	real(8), private :: C8
	real(8), private :: C10
	real(8), private :: C12
	real(8), private :: C14
	
	!*
	! @brief Clase que representa el potencial de interaccion AuTiO2(110)
	!*
	type, public :: AuTiO2Potential
		integer, private :: model
		logical, private :: grimme
		
		contains
			procedure :: init
			procedure :: destroy
			procedure :: V
			procedure :: dV
			procedure, private :: V1D
			procedure, private :: V3D
			procedure, private :: NdV
			procedure, private :: dV1D
			procedure, private :: dV3D
	end type AuTiO2Potential
	
	contains
	
	!*
	! @brief Contructor
	!*
	subroutine init( this, model, grimme )
		class(AuTiO2Potential) :: this
		integer, optional, intent(in) :: model
		logical, optional, intent(in) :: grimme
		
		integer :: modelEff
		integer :: grimmeEff
		
		if( present(model) ) then
			modelEff = model
		else
			modelEff = SMSP
		end if
		
		if( present(grimme) ) then
			grimmeEff = grimme
		else
			grimmeEff = .false.
		end if
		
		if( modelEff == SM .or. modelEff == SM1Dz ) then
		
			aDe = 0.405721_8
			bDe = 0.012854_8
			cDe = -0.034879_8
			gammaDe = 0.620000_8
			tauDe = -0.160000_8
			
			aAlpha = 2.86538_8
			bAlpha = 1.46072_8
			cAlpha = 0.640971_8
			gammaAlpha = 1.000000_8
			tauAlpha = -0.127564_8
			
			aBeta = 15.9827_8
			bBeta = -7.93087_8
			cBeta = -0.242313_8
			gammaBeta = 1.000000_8
			tauBeta = -0.112711_8
			
			aZe = 0.272486_8
			bZe = -0.014505_8
			cZe = 0.007032_8
			gammaZe = 0.733418_8
			tauZe = -0.0415664_8
			
		else if( modelEff == SMSP .or. modelEff == SMSP1Dz ) then
		
			aDe = 0.359541_8
			bDe = 0.0289728_8
			cDe = -0.045617_8
			gammaDe = 0.648537_8
			tauDe = -0.18907_8
			
			aAlpha = 1.96388_8
			bAlpha = 0.338402_8
			cAlpha = 0.464906_8
			gammaAlpha = 1.56058_8
			tauAlpha = -0.0324907_8
			
			aBeta = 28.6901_8
			bBeta = -5.24277_8
			cBeta = -7.374_8
			gammaBeta = 1.73418_8
			tauBeta = -0.176578_8
			
			aZe = 0.273244_8
			bZe = -0.0144939_8
			cZe = 0.00705536_8
			gammaZe = 0.758325_8
			tauZe = -0.047944_8
			
		else if( modelEff == SMSPDisp1Dz ) then
			
			rcut = 2.0_8
			maxV = 1.81634739890781_8
			mindV = -4.91534910424389_8
			C3 = -6.7455_8
			C4 = 82.5543_8
			C6 = -6044.46_8
			C8 = 88471.3_8
			C10 = -518462.0_8
			C12 = 1.39401D6
			C14 = -1.43623D6
			
		end if
		
		this.model = modelEff
		this.grimme = grimmeEff
	end subroutine init
	
	!*
	! @brief Destructor
	!*
	subroutine destroy( this )
		class(AuTiO2Potential), intent(in) :: this
		
	end subroutine destroy
	
	!*
	! @brief
	!*
	function V( this, x, y, z ) result( output )
		class(AuTiO2Potential), intent(in) :: this
		real(8), intent(in) :: x
		real(8), intent(in) :: y
		real(8), intent(in) :: z
		real(8) :: output
		
		if( this.model == SM1Dz .or. this.model == SMSP1Dz ) then
			output = this.V3D( 0.0_8, bCell/3.0_8, z )
		else if( this.model == SMSPDisp1Dz ) then
			output = V_SMSPDisp1Dz( z )
		else
			output = this.V3D( x, y, z )
		end if
	end function V
	
	!*
	! @brief
	!*
	function dV( this, x, y, z, coord ) result( output )
		class(AuTiO2Potential), intent(in) :: this
		real(8), intent(in) :: x
		real(8), intent(in) :: y
		real(8), intent(in) :: z
		integer, intent(in) :: coord
		real(8) :: output
		
		if( this.model == SM1Dz .or. this.model == SMSP1Dz ) then
			select case( coord )
				case( 1 )
					output = 0.0_8
				case( 2 )
					output = 0.0_8
				case( 3 )
					output = this.dV3D( 0.0_8, bCell/3.0_8, z, 3 )
			end select
		else if( this.model == SMSPDisp1Dz ) then
			select case( coord )
				case( 1 )
					output = 0.0_8
				case( 2 )
					output = 0.0_8
				case( 3 )
					output = dV_SMSPDisp1Dz( z )
			end select
		else
			output = this.dV3D( x, y, z, coord )
		end if
	end function dV
	
	!*
	! @brief Funcion generica que se encarga de llamar a la funcion
	!        de potencial correcto, segun el atributo model, pero
	!        a lo largo de la coordenada z, siendo x,y el mínimo
	!
	! @param z
	!        Coordenada z en angstroms
	!*
	function V1D( this, z ) result( output )
		class(AuTiO2Potential), intent(in) :: this
		real(8), intent(in) :: z
		real(8) :: output
		
! 		output = this.V3D( 0.0_8, 1.0377600_8, z )
		output = this.V3D( 0.0_8, bCell/3.0_8, z )
		
	end function V1D
	
	!*
	! @brief Funcion generica que se encarga de llamar a la funcion
	!        derivada del potencial correcto, segun el atributo model,
	!        pero a lo largo de la coordenada z, siendo x,y el mínimo
	!
	! @param z
	!        Coordenada z en angstroms
	!*
	function dV1D( this, z ) result( output )
		class(AuTiO2Potential), intent(in) :: this
		real(8), intent(in) :: z
		real(8) :: output
		
! 		output = this.dV3D( 0.0_8, 1.0377600_8, z, 3 )
		output = this.dV3D( 0.0_8, bCell/3.0_8, z, 3 )
		
	end function dV1D
	
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
	function V3D( this, x, y, z ) result( output )
		class(AuTiO2Potential), intent(in) :: this
		real(8), intent(in) :: x
		real(8), intent(in) :: y
		real(8), intent(in) :: z
		real(8) :: output
		
		real(8) :: yEff
		real(8) :: zEff
		
		yEff = abs(y/bCell)
		zEff = z/29.0_8+0.166318_8
		
		yEff = yEff - real( int(yEff), 8 )
		
		if( yEff >= 0.0_8 .and. yEff <= 0.5_8 ) then
			output = V_SM2D( x, 2.0_8*yEff, zEff )
		else if( yEff > 0.5_8 .and. yEff <= 1.0_8 ) then
			output = V_SM2D( x, 2.0_8*(1.0_8-yEff), zEff )
! 		else if( abs(yEff-1.0_8) < 1e-6 ) then
! 			output = V_SM2D( x, 0.0_8, zEff )
		end if
		
! 		output = output*9678.64917839637_8 ! eV to cm-1
		
#define Ze        F(yEff, aZe, bZe, cZe, gammaZe, tauZe)
		output = (1000000.0_8*exp(-(z-Ze)/0.1_8)+output)*9678.64917839637_8 ! eV to cm-1
#undef Ze
		
	end function V3D
	
	function NdV( this, x, y, z, coord, stepSize, nPoints ) result( output )
		class(AuTiO2Potential), intent(in) :: this
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
		
		h = 0.000001_8
		if( present(stepSize) ) h = stepSize

#define fx(i) this.V(x+i##.0_8*h,y,z)
#define fy(i) this.V(x,y+i##.0_8*h,z)
#define fz(i) this.V(x,y,z+i##.0_8*h)
		output = 0.0_8
		select case( nPointsEff )
			case( 3 )
				select case( coord )
					case( 1 )
! 						output = ( fx(1)-fx(-1) )/(2.0_8*h)
						output = 0.0_8
					case( 2 )
						output = ( fy(1)-fy(-1) )/(2.0_8*h)
					case( 3 )
						output = ( fz(1)-fz(-1) )/(2.0_8*h)
				end select
			case( 5 )
				select case( coord )
					case( 1 )
! 						output = ( fx(-2)-8.0_8*fx(-1)+8.0_8*fx(1)-fx(2) )/(12.0_8*h)
						output = 0.0_8
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
	function dV3D( this, x, y, z, coord ) result( output )
		class(AuTiO2Potential), intent(in) :: this
		real(8), intent(in) :: x
		real(8), intent(in) :: y
		real(8), intent(in) :: z
		integer, intent(in) :: coord
		real(8) :: output
		
		real(8) :: yEff
		real(8) :: zEff
		
		yEff = abs(y/bCell)
		zEff = z/29.0_8+0.166318_8
		
		yEff = yEff - real( int(yEff), 8 )
		
		if( coord == 1 .or. zEff > 0.4_8 ) then
			output = 0.0_8
			return
		end if
		
		output = 0.0_8
		if( y<0.0_8 ) then
			if( yEff > 0.0_8 .and. yEff <= 0.5_8 ) then
				output = -dV_SM2D( x, 2.0_8*yEff, zEff, coord )
			else if( yEff > 0.5_8 .and. yEff <= 1.0_8 ) then
				output = dV_SM2D( x, 2.0_8*(1.0_8-yEff), zEff, coord )
			else if( abs(yEff) < 1e-6 ) then
				output = dV_SM2D( x, 0.0_8, zEff, coord )
			end if
		else
			if( yEff > 0.0_8 .and. yEff <= 0.5_8 ) then
				output = dV_SM2D( x, 2.0_8*yEff, zEff, coord )
			else if( yEff > 0.5_8 .and. yEff <= 1.0_8 ) then
				output = -dV_SM2D( x, 2.0_8*(1.0_8-yEff), zEff, coord )
			else if( abs(yEff) < 1e-6 ) then
				output = dV_SM2D( x, 0.0_8, zEff, coord )
			end if
		end if
		
! 		output = output + 1000.0_8*x
! 		output = output*9678.64917839637_8 ! eV to cm-1
		
#define Ze        F(yEff, aZe, bZe, cZe, gammaZe, tauZe)
		output = (-10000000.0_8*exp(-(z-Ze)/0.1_8)+output)*9678.64917839637_8 ! eV to cm-1
#undef Ze
	end function dV3D
	
	!*
	! @brief Potencial segun el modelo superMorse2D
	!*
	function V_SM2D( x, y, z ) result( output )
		real(8), intent(in) :: x
		real(8), intent(in) :: y
		real(8), intent(in) :: z
		real(8) :: output
		
#define De    F(y, aDe, bDe, cDe, gammaDe, tauDe)
#define alpha F(y, aAlpha, bAlpha, cAlpha, gammaAlpha, tauAlpha)
		
		output = De*( exp(-2.0_8*alpha*P(y,z))-2.0_8*exp(-alpha*P(y,z)) )
#undef De
#undef alpha
	end function V_SM2D
	
	!*
	! @brief Derivada del potencial segun el modelo superMorse2D
	!*
	function dV_SM2D( x, y, z, coord ) result( output )
		real(8), intent(in) :: x
		real(8), intent(in) :: y
		real(8), intent(in) :: z
		integer, intent(in) :: coord
		real(8) :: output
		
#define Vxyz      V_SM2D(x,y,z)
#define De        F(y, aDe, bDe, cDe, gammaDe, tauDe)
#define alpha     F(y, aAlpha, bAlpha, cAlpha, gammaAlpha, tauAlpha)
#define dDedy    dF(y, aDe, bDe, cDe, gammaDe, tauDe)
#define dalphady dF(y, aAlpha, bAlpha, cAlpha, gammaAlpha, tauAlpha)
	
		output = 0.0_8
		select case ( coord )
			case( 1 )
! 				output = 1.0d6
				output = 0.0_8
			case( 2 )
				output = ( (Vxyz/De)*dDedy + 2.0_8*De*exp(-2.0_8*P(y,z)*alpha) &
				  *( exp(alpha*P(y,z))-1.0_8 )*( P(y,z)*dalphady + alpha*dP(y,z,2) ) )/(0.5*bCell)
			case( 3 )
				output = 2.0_8*De*alpha*exp(-2.0_8*P(y,z)*alpha)*( exp(alpha*P(y,z))-1.0_8 )*dP(y,z,3)/29.0_8
		end select
	
#undef Vxyz
#undef De
#undef alpha
#undef dDedy
#undef dalphady
	end function dV_SM2D
	
	!*
	! @brief Funcion auxiliar P
	!*
	function P( y, z ) result( output )
		real(8), intent(in) :: y
		real(8), intent(in) :: z
		real(8) :: output
		
#define beta      F(y, aBeta, bBeta, cBeta, gammaBeta, tauBeta)
#define Ze        F(y, aZe, bZe, cZe, gammaZe, tauZe)
		
		output = exp(beta*(z-Ze))-1.0_8
		
#undef beta
#undef Ze
	end function P
	
	!*
	! @brief Derivada de la funcion auxiliar P
        ! @param coord
        !        Coordenada respecto a la cual se deriva [ ( 1, 2, 3 ) => ( x, y, z ) ]
	!*
	function dP( y, z, coord ) result( output )
		real(8), intent(in) :: y
		real(8), intent(in) :: z
		integer, intent(in) :: coord
		real(8) :: output
		
#define beta      F(y, aBeta, bBeta, cBeta, gammaBeta, tauBeta)
#define Ze        F(y, aZe, bZe, cZe, gammaZe, tauZe)
#define dbetady  dF(y, aBeta, bBeta, cBeta, gammaBeta, tauBeta)
#define dZedy    dF(y, aZe, bZe, cZe, gammaZe, tauZe)
		
		select case( coord )
			case( 1 )
				output = 0.0_8
			case( 2 )
				output = exp(beta*(z-Ze))*( (z-Ze)*dbetady-beta*dZedy )
			case( 3 )
				output = beta*exp(beta*(z-Ze))
		end select
		
#undef beta
#undef Ze
#undef dbetady
#undef dZedy
	end function dP
	
	!*
	! @brief Funcion auxiliar F
	!*
	pure function F( y, a, b, c, gamma, tau ) result( output )
		real(8), intent(in) :: y
		real(8), intent(in) :: a
		real(8), intent(in) :: b
		real(8), intent(in) :: c
		real(8), intent(in) :: gamma
		real(8), intent(in) :: tau
		real(8) :: output
		
		real(8) :: t
		
		t = y**gamma + tau*sin(2.0_8*pi*y**gamma)
		output = a + b*cos(pi*t) + c*cos(2.0_8*pi*t)
	end function F
	
	!*
	! @brief Derivada de la funcion auxiliar F
	!*
	pure function dF( y, a, b, c, gamma, tau ) result( output )
		real(8), intent(in) :: y
		real(8), intent(in) :: a
		real(8), intent(in) :: b
		real(8), intent(in) :: c
		real(8), intent(in) :: gamma
		real(8), intent(in) :: tau
		real(8) :: output
		
		real(8) :: t
		real(8) :: dtdy
		
		t = y**gamma + tau*sin(2.0_8*pi*y**gamma)
		dtdy = ( gamma*y**(gamma-1.0_8)*(1.0_8+2.0_8*pi*tau*cos(2.0_8*pi*y**gamma)) )
		
		output = -pi*sin(pi*t)*( b+4.0_8*c*cos(pi*t) )*dtdy
! 		output = -pi*( b*sin(pi*t)+2.0_8*c*sin(2.0_8*pi*t) )*dtdy/bCell
	end function dF
	
	!*
	! @brief Potencial segun el modelo SMSPDisp1Dz
	!*
	function V_SMSPDisp1Dz( z ) result( output )
		real(8), intent(in) :: z
		real(8) :: output
		
		if( z > 2.0_8 ) then
			output = C3/z**3 + C4/z**4 + C6/z**6 + C8/z**8 + C10/z**10 + C12/z**12 + C14/z**14
		else
			output = 1.82618769531726_8
		end if

! 		output = 2.0_8*0.5_8*( 1.0_8 - tanh( 20.0_8*(z-2.0_8) ) )
		
		output = output*eV/cm1
	end function V_SMSPDisp1Dz
	
	!*
	! @brief Derivada del potencial segun el modelo SMSPDisp1Dz
	!*
	function dV_SMSPDisp1Dz( z ) result( output )
		real(8), intent(in) :: z
		real(8) :: output
		
		if( z > 2.0_8 ) then
			output = -3.0_8*C3/z**4 - 6.0_8*C6/z**7 - 8.0_8*C8/z**9 - 10.0_8*C10/z**11 - 12.0_8*C12/z**13 - 14.0_8*C14/z**15 - 4.0_8*C4/z**5
		else
			output = 0.0_8
		end if

! 		if( z>-5.0 .and. z<15.0_8 ) then
! 			output = -20.0_8/cosh( 20.0_8*(z-2.0_8) )**2
! 		else
! 			output =0.0_8
! 		end if
		
		output = output*eV/cm1
	end function dV_SMSPDisp1Dz
	
	!*
	! @brief Subrutina encargada de varificar el coorecto
	!        funcionamiento de la clase
	!*
	subroutine AuTiO2Potential_test()
		real(8) :: x, y, z
		real(8) :: sum
		type(AuTiO2Potential) :: potential
		
! 		call potential.init( model=SMSP, grimme=.false. )
! 		call potential.init( model=SM, grimme=.false. )
		call potential.init( model=SMSPDisp1Dz )
		
! 		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 		! set contours
! 		! set cntrparam levels 100
! 		! splot [] [] [-0.5:0.5] "salida" w l
! 		write(*,"(1A,A29,4A30)") "#", "x", "y", "z", "V(x,y,z)", "dV(x,y,z)/dx"
! 		
! 		z = 2.0_8
! 		do while( z <= 5.5_8 )
! 		
! ! 			x = -aCell/2.0_8
! ! 			do while( x <= aCell/2.0_8 )
! 
! ! 				y = -bCell/2.0_8
! ! 				do while( y <= bCell/2.0_8 )
! 				y = -bCell
! 				do while( y <= bCell )
! 				
! 					write(*,"(5F30.7)") x, y, z, potential.V( x, y, z ), potential.NdV( x, y, z, 2 )
! 					
! ! 					y = y + bCell/50.0_8
! 					y = y + bCell/200.0_8
! 				end do
! 				
! ! 				x = x + 0.1_8
! ! 				write(*,*) ""
! ! 			end do
! 			
! ! 			z = z + 0.1_8
! 			z = z + 0.01_8
! 			
! 			write(*,*) ""
! 		end do

		
! 		write(*,*) ""
! 		write(*,*) ""
! 		
! 		write(*,"(1A,4A20)") "#", "x", "V(x,y,z)", "dV(x,y,z)/dx", "N[dV(x,y,z)/dx]"
! 		
! 		z = 2.8_8
! 		y = bCell/2.0_8
! 		x = -2.0*aCell
! 		do while( x <= 2.001*aCell )
! 			write(*,"(4F30.7)") x, potential.V( x, y, z ), &
! 				potential.dV( x, y, z, 1 ), potential.NdV( x, y, z, 1 )
! 			
! 			x = x + aCell/101.0_8
! 		end do
! 		
! 		write(*,*) ""
! 		write(*,*) ""
! 		
! 		write(*,"(1A,4A20)") "#", "y", "V(x,y,z)", "dV(x,y,z)/dy", "N[dV(x,y,z)/dy]"
! 		
! 		z = -22.0_8
! ! 		z = 2.8_8
! 		y = -bCell
! 		do while( y <= bCell )
! 			write(*,"(4F30.7)") y, potential.V( 0.0_8, y, z+25.0_8 ), &
! 				potential.dV( 0.0_8, y, z+25.0_8, 2 ), potential.NdV( 0.0_8, y, z+25.0_8, 2 )
! ! 			write(*,"(4F20.7)") y, potential.V( 0.0_8, y, z ), &
! ! 				potential.dV( 0.0_8, y, z, 2 ), potential.NdV( 0.0_8, y, z, 2 )
! 			
! 			y = y + 2.0_8*bCell/101.0_8
! 		end do
! 		
! 		write(*,*) ""
! 		write(*,*) ""
! 		
		write(*,"(1A,4A20)") "#", "z", "V(x,y,z)", "dV(x,y,z)/dz", "N[dV(x,y,z)/dz]"
!  		y = bCell/2.0_8
 		y = bCell/3.0_8
		z = 1.0_8
		do while( z <= 7.5_8 )
! 		z = -30.0_8
! 		do while( z <= -15.0_8 )
! 			write(*,"(6F30.7)") z, potential.V( 0.0_8, y, z+25.0_8 ), &
! 				potential.dV( 0.0_8, y, z+25.0_8, 3 ), potential.NdV( 0.0_8, y, z+25.0_8, 3 )
			write(*,"(6F20.7)") z, potential.V( 0.0_8, y, z ), &
				potential.dV( 0.0_8, y, z, 3 ), potential.NdV( 0.0_8, y, z, 3 )
			
			z = z + 0.01_8
		end do
		
! 		write(*,"(1A,4A20)") "#", "z", "V(x,y,z)", "dV(x,y,z)/dz", "N[dV(x,y,z)/dz]"
! 		z = 1.8_8
! 		do while( z <= 7.0_8 )
! 			write(*,"(2F30.7)") z, potential.V( z )
! 			
! 			z = z + 0.01_8
! 		end do
		
		call potential.destroy()
		
	end subroutine AuTiO2Potential_test
	
end module AuTiO2Potential_
