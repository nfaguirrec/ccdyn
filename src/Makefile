MAKEFLAGS = -s
FC = ifort
FCFLAGS = -traceback -mkl -fpic -w -I. -g -fpp -fpe0 -I${HOME}/Develop/scift/src
LDFLAGS = -L${HOME}/Develop/scift/src -lscift
TLIB = 

all: ArArPotential.o AuAuPotential.o AuTiO2Potential.o ClassicalDynamics.o CoCoPotential.o EmbeddingInterface.o HeAuPotential.o HeDroplet.o HeHePotential.o HeMgPotential.o HeTiO2Potential.o MDIntegrator.o MgMgPotential.o NiNiPotential.o NiOPotential.o  $(TLIB) ccdyn HeDropletInitialGuess test 

ArArPotential.o ararpotential_.mod: ArArPotential.f90 
	fmake
	echo -n 'Building ArArPotential.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) -c ArArPotential.f90

AuAuPotential.o auaupotential_.mod: AuAuPotential.f90 
	fmake
	echo -n 'Building AuAuPotential.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) -c AuAuPotential.f90

AuTiO2Potential.o autio2potential_.mod: AuTiO2Potential.f90 
	fmake
	echo -n 'Building AuTiO2Potential.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) -c AuTiO2Potential.f90

ccdyn: ccdyn.f90 HeHePotential.o AuTiO2Potential.o ArArPotential.o CoCoPotential.o MgMgPotential.o HeAuPotential.o AuAuPotential.o HeDroplet.o ClassicalDynamics.o EmbeddingInterface.o NiNiPotential.o NiOPotential.o MDIntegrator.o HeMgPotential.o HeTiO2Potential.o hehepotential_.mod autio2potential_.mod ararpotential_.mod cocopotential_.mod mgmgpotential_.mod heaupotential_.mod auaupotential_.mod hedroplet_.mod classicaldynamics_.mod embeddinginterface_.mod ninipotential_.mod niopotential_.mod mdintegrator_.mod hemgpotential_.mod hetio2potential_.mod
	fmake
	echo -n 'Building ccdyn.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) ccdyn.f90 HeHePotential.o AuTiO2Potential.o ArArPotential.o CoCoPotential.o MgMgPotential.o HeAuPotential.o AuAuPotential.o HeDroplet.o ClassicalDynamics.o EmbeddingInterface.o NiNiPotential.o NiOPotential.o MDIntegrator.o HeMgPotential.o HeTiO2Potential.o -o ccdyn $(LDFLAGS)

ClassicalDynamics.o classicaldynamics_.mod: ClassicalDynamics.f90 hehepotential_.mod ararpotential_.mod cocopotential_.mod mgmgpotential_.mod heaupotential_.mod auaupotential_.mod hedroplet_.mod embeddinginterface_.mod ninipotential_.mod niopotential_.mod mdintegrator_.mod hemgpotential_.mod
	fmake
	echo -n 'Building ClassicalDynamics.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) -c ClassicalDynamics.f90

CoCoPotential.o cocopotential_.mod: CoCoPotential.f90 
	fmake
	echo -n 'Building CoCoPotential.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) -c CoCoPotential.f90

EmbeddingInterface.o embeddinginterface_.mod: EmbeddingInterface.f90 
	fmake
	echo -n 'Building EmbeddingInterface.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) -c EmbeddingInterface.f90

HeAuPotential.o heaupotential_.mod: HeAuPotential.f90 
	fmake
	echo -n 'Building HeAuPotential.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) -c HeAuPotential.f90

HeDroplet.o hedroplet_.mod: HeDroplet.f90 
	fmake
	echo -n 'Building HeDroplet.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) -c HeDroplet.f90

HeDropletInitialGuess: HeDropletInitialGuess.f90  
	fmake
	echo -n 'Building HeDropletInitialGuess.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) HeDropletInitialGuess.f90  -o HeDropletInitialGuess $(LDFLAGS)

HeHePotential.o hehepotential_.mod: HeHePotential.f90 
	fmake
	echo -n 'Building HeHePotential.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) -c HeHePotential.f90

HeMgPotential.o hemgpotential_.mod: HeMgPotential.f90 
	fmake
	echo -n 'Building HeMgPotential.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) -c HeMgPotential.f90

HeTiO2Potential.o hetio2potential_.mod: HeTiO2Potential.f90 
	fmake
	echo -n 'Building HeTiO2Potential.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) -c HeTiO2Potential.f90

MDIntegrator.o mdintegrator_.mod: MDIntegrator.f90 
	fmake
	echo -n 'Building MDIntegrator.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) -c MDIntegrator.f90

MgMgPotential.o mgmgpotential_.mod: MgMgPotential.f90 
	fmake
	echo -n 'Building MgMgPotential.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) -c MgMgPotential.f90

NiNiPotential.o ninipotential_.mod: NiNiPotential.f90 
	fmake
	echo -n 'Building NiNiPotential.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) -c NiNiPotential.f90

NiOPotential.o niopotential_.mod: NiOPotential.f90 
	fmake
	echo -n 'Building NiOPotential.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) -c NiOPotential.f90

test: test.f90 AuTiO2Potential.o HeHePotential.o ArArPotential.o MgMgPotential.o CoCoPotential.o HeAuPotential.o AuAuPotential.o EmbeddingInterface.o HeDroplet.o ClassicalDynamics.o NiNiPotential.o NiOPotential.o MDIntegrator.o HeMgPotential.o HeTiO2Potential.o autio2potential_.mod hehepotential_.mod ararpotential_.mod mgmgpotential_.mod cocopotential_.mod heaupotential_.mod auaupotential_.mod embeddinginterface_.mod hedroplet_.mod classicaldynamics_.mod ninipotential_.mod niopotential_.mod mdintegrator_.mod hemgpotential_.mod hetio2potential_.mod
	fmake
	echo -n 'Building test.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) test.f90 AuTiO2Potential.o HeHePotential.o ArArPotential.o MgMgPotential.o CoCoPotential.o HeAuPotential.o AuAuPotential.o EmbeddingInterface.o HeDroplet.o ClassicalDynamics.o NiNiPotential.o NiOPotential.o MDIntegrator.o HeMgPotential.o HeTiO2Potential.o -o test $(LDFLAGS)

clean:
	rm -f *.o *.mod *~ .deps ccdyn HeDropletInitialGuess test  2> /dev/null

distclean:
	rm -f *.o *.mod *~ .deps Makefile ccdyn HeDropletInitialGuess test  2> /dev/null

Makefile: .deps
	fmake

.deps:
	fmake
