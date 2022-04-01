FFLAGS := -O3
LFLAGS := -L$(ALGENCAN)/sources/algencan/lib -L$(ALGENCAN)/sources/hsl/lib -L$(ALGENCAN)/sources/blas/lib
LDFLAGS := -lalgencan -lhsl -lblas 

OBJ := drand.o geometry.o geompack2.o modamerica.o modcesaro.o vorintpols.o

all: covering.o $(OBJ)
	gfortran $^ $(LFLAGS) -o covering $(LDFLAGS)

covering.o: covering.f90 $(OBJ)
	gfortran $(FFLAGS) -c -I$(ALGENCAN)/sources/algencan/inc $<

%.o: %.f90
	gfortran $(FFLAGS) -c $<

clean:
	rm -f $(OBJ) covering.o
	rm -f modamerica.mod modcesaro.mod vorcells_polygons.mod
	rm -f covering

cleanimages:
	rm -f picture-*.mps
	rm -f picture-*.log
	rm -f picture-*.mp
	rm -f picture-*.pdf

cleanoutputs:
	rm -f output*.txt
	rm -f output.dat
	rm -f *tabline.txt
	rm -f bestrad.txt
