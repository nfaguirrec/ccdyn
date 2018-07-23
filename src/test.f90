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

program test_
	use HeDroplet_
	use ClassicalDynamics_
	use MDIntegrator_
	use HeHePotential_
	use HeAuPotential_
	use AuAuPotential_
	use ArArPotential_
	use AuTiO2Potential_
	use HeTiO2Potential_
	use AtomicElementsDB_
	use NiOPotential_
	use NiNiPotential_
	use CoCoPotential_
	implicit none
	
! 	call HeDroplet_test()
!	call HeTiO2Potential_test()
! 	call HeHePotential_test()
! 	call HeAuPotential_test()
!  	call AuAuPotential_test()
!  	call ArArPotential_test()
!  	call AuTiO2Potential_test()
!  	call ClassicalDynamics_test()
! 	call AtomicElementsDB_test()
! 	call NiOPotential_test()
! 	call NiNiPotential_test()
	call CoCoPotential_test()
	
end program test_
