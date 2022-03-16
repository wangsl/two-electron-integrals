
# $Id: Makefile,v 1.1 2005/07/18 20:32:28 wangsl Exp wangsl $

.DEFAULT: .f .F .for .c .C .cpp .cc .f90 .F90
.SUFFIXES: .f .F .for .c .C .cpp .cc .f90 .F90

O = .

F77 = gfortran
CC = gcc
CCC = g++

CFLAGS = -O3 -fPIC -fopenmp -I/home/wang/PyQuante-1.6.5/libint-release-1-1-4/include

FFLAGS = -O3 -fPIC -fopenmp

Link = $(CCC) $(CFLAGS) 

LIBS = -L/home/wang/PyQuante-1.6.5/libint-release-1-1-4/lib -lint -lgfortran

EXENAME = test.x

OBJS = $(O)/test-2.o $(O)/gauss_hermite.o  $(O)/rys_chebyshev_coeffs.o \
	$(O)/indent.o  $(O)/die.o \
	$(O)/primitive_gaussian.o \
	$(O)/rys_2e_integral.o  $(O)/rys_2e_integralio.o \
	$(O)/crys.o  \
	$(O)/fmt.o  $(O)/libint-test.o #$(O)/fmtgen.o $(O)/gfloat.o $(O)/gamgen.o $(O)/fmtset.o

#$(O)/test-2.o $(O)/rys-test.o

$(EXENAME) : $(OBJS) 
	$(Link) -o $(EXENAME) $(OBJS) $(LIBS)

$(O)/%.o: %.c
	cd $(O) ; $(CC) $(CFLAGS) -c $<
$(O)/%.o: %.cc
	cd $(O) ; $(CCC) $(CFLAGS) -c $<
$(O)/%.o: %.cpp
	cd $(O) ; $(CCC) $(CFLAGS) -c $<
$(O)/%.o: %.C
	cd $(O) ; $(CCC) $(CFLAGS) -c $<
$(O)/%.o: %.F
	cd $(O) ; $(F77) $(FFLAGS) -c $<
$(O)/%.o: %.for
	cd $(O) ; $(F77) $(FFLAGS) -c $<
$(O)/%.o: %.f90
	cd $(O) ; $(F90) $(FFLAGS) -c $<

%io.C: %.h
	perl io.pl $<

clean:
	rm -f *.o *~ $(EXENAME) *.mod *io.C depend

.PRECIOUS: %io.C %.C %.c %.f %.F %.h $(O)/%.o

depend :
	$(CXX) $(CXXFLAGS) -MM *.[cC] | perl dep.pl > $@
	#bash dep.sh *.F >> $@

ifneq ($(MAKECMDGOALS), clean)
include depend
endif
