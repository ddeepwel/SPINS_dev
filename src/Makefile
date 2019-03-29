# Makefile to build SPINS, generalized for multiple systemd

# The "system.mk" file in the current directory will contain
# system-specific make variables, most notably the
# C/C++ compiler/linker, 
include system.mk
include VERSION


# Compile with debugging flags
DEBUG?=false
# Compile with optimizations
OPTIM?=true
# Compile with extra optimizations that may require significantly
# longer for compilation
SLOW_OPTIM?=false

# Compile with instrumentation support for timing
TIMINGS?=false

# If MPICXX isn't separately defined, then set it to be the same
# as CXX
ifeq ($(strip $(MPICXX)),)
   MPICXX:=$(CXX)
endif

# Flag to pass version info to code
VERSION:= -DMAJOR_VERSION=${MAJOR_VERSION} -DMINOR_VERSION=${MINOR_VERSION} -DPATCH_VERSION=${PATCH_VERSION}

# Assemble the CFLAGS
CFLAGS:=$(SYSTEM_CFLAGS) $(MPI_CFLAGS) $(SYSTEM_CXXFLAGS)
LDFLAGS:=$(SYSTEM_LDFLAGS)
ifeq ($(DEBUG),true)
   CFLAGS:=$(CFLAGS) $(DEBUG_CFLAGS)
   LDFLAGS:=$(LDFLAGS) $(DEBUG_LDFLAGS)
endif
ifeq ($(OPTIM),true)
   CFLAGS:=$(CFLAGS) $(OPTIM_CFLAGS)
   LDFLAGS:=$(LDFLAGS) $(OPTIM_LDFLAGS)
   ifeq ($(SLOW_OPTIM),true)
      CFLAGS:=$(CFLAGS) $(EXTRA_OPTIM_CFLAGS)
      LDFLAGS:=$(LDFLAGS) $(EXTRA_OPTIM_LDFLAGS)
   endif
endif

ifeq ($(TIMINGS),true)
   CFLAGS:=$(CFLAGS) -DTIMING_ENABLE
endif


INCLUDE_DIRS := -I../include $(MPI_INCDIR) $(LAPACK_INCDIR) $(BLITZ_INCDIR) $(FFTW_INCDIR) $(UMF_INCDIR)
CFLAGS := $(CFLAGS) $(INCLUDE_DIRS)

LIB_DIRS := -L../lib $(MPI_LIBDIR) $(LAPACK_LIBDIR) $(BLITZ_LIBDIR) $(FFTW_LIBDIR) $(UMF_LIBDIR)
LDLIBS := -lfftw3 $(MPI_LIB) -lumfpack -lamd -lboost_program_options $(LAPACK_LIB) -lblitz -lm -lstdc++ 
LDFLAGS := $(LDFLAGS) $(LIB_DIRS) 

.PHONY: all 
all: tests/test_deriv_x tests/test_write_x tests/test_esolve_x tests/test_heat_x tests/test_ode_x tests/test_ns_x

.PHONY: clean
clean:
	rm -f *.o tests/*.o cases/*.o cases/*.src.? tests/*.src.? Science/*.o BaseCase/*.o cases/*/*.src.? cases/*/*.o

##
## Short-names for source files
##

objfiles: $(shell ls *.cpp | sed -e 's/cpp/o/g')

# SPINS files
SPINS_BASE := TArray.o T_util.o Parformer.o ESolver.o Timestep.o NSIntegrator.o Splits.o Par_util.o Split_reader.o gmres.o gmres_1d_solver.o gmres_2d_solver.o grad.o multigrid.o Options.o BaseCase.o Sorter.o timing.o

BASECASE_CPPS := $(wildcard BaseCase/*.cpp)
BASECASE_OBJS := $(addprefix BaseCase/,$(notdir $(BASECASE_CPPS:.cpp=.o)))

SCIENCE_CPPS := $(wildcard Science/*.cpp)
SCIENCE_OBJS := $(addprefix Science/,$(notdir $(SCIENCE_CPPS:.cpp=.o)))

##
## Recipes
##

NSIntegrator.o: NSIntegrator.cpp NSIntegrator_impl.cc

tests/test%.o: tests/test%.cpp
	$(MPICXX) ${VERSION} $(CFLAGS) -o $@ -c $<

tests/test%.x: tests/test%.o tests/test%.src.o ${SCIENCE_OBJS} ${BASECASE_OBJS} ${SPINS_BASE}
	$(LD) ${VERSION} $(LDFLAGS) -o $@ $^ $(LDLIBS)

cases/%.o: cases/%.cpp NSIntegrator_impl.cc NSIntegrator.hpp
	$(MPICXX) ${VERSION} $(CFLAGS) -o $@ -c  $<

cases/%.x: cases/%.o cases/%.src.o ${SCIENCE_OBJS} ${BASECASE_OBJS} ${SPINS_BASE}
	$(LD) ${VERSION} $(LDFLAGS) -o $@ $^ $(LDLIBS)

cases/%_x: cases/%.o cases/%.src.o ${SCIENCE_OBJS} ${BASECASE_OBJS} ${SPINS_BASE}
	$(LD) ${VERSION} $(LDFLAGS) -o $@ $^ $(LDLIBS)

tests/test%.src.c: tests/test%.cpp CaseFileSource.c
	echo "const char casefilesource[] = {`xxd -i < $<`, 0x00};" > $@
	echo "const char casefilename[] = \"$<\";" >> $@
	cat CaseFileSource.c >> $@

cases/%.src.c: cases/%.cpp CaseFileSource.c
	echo "const char casefilesource[] = {`xxd -i < $<`, 0x00};" > $@
	echo "const char casefilename[] = \"$<\";" >> $@
	cat CaseFileSource.c >> $@

%.o: %.cpp	*.hpp
	$(MPICXX) ${VERSION} $(CFLAGS) -o $@ -c $< 

print-%  : ; @echo $* = $($*)

.SECONDARY:
