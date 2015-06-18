FC = ifort
OMP = -openmp 
FFLAGS = -fpp  -mkl  $(OMP) -fpe0 -fast
LFLAGS = $(OMP) -mkl
OBJECTS = commonModules.o dataLoader.o pointNeighbor.o implicit.o  \
		  mLaplace.o biconjGrad.o smoothing.o ncalcRHS.o \
		  subrutinas.o meshMove.o ns2DComp.ALE.o 

.PHONY: clean

all: ns

xhost: FFLAGS += -xHost 
xhost: ns

mkl: FFLAGS += -mkl
mkl: LFLAGS =+ -mkl
mkl: ns

debug: FFLAGS = -O0 -fpp $(OMP) -check -traceback -warn nounused -mkl
debug: ns

idb: FFLAGS = -debug -O0 -fpp $(OMP)
idb: LFLAGS = -debug $(OMP)
idb: ns

profgen: FFLAGS = -fpp -prof-gen $(OMP)
profgen: ns

profile: FFLAGS += -p
profile: ns

ns: $(OBJECTS)
	$(FC) $(LFLAGS) $(OBJECTS) -o ns 

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) ns *.mod
