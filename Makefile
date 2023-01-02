# Disable the default rules
MAKEFLAGS += --no-builtin-rules --no-builtin-variables

# Project name
NAME := csv

# Configuration settings
FC := gfortran
AR := ar rcs
LD := $(FC)
RM := rm -f

# List of all source files
SRCS := adke.for \
        art_heat.f \
        art_visc.f \
        av_vel.f \
        beeman.for \
        break_dam.for \
        combined_visc.for \
        density.f \
        direct_find.f \
        eos.f \
	ex.f90 \
        external_force.f \
        ghost_part.for \
        gmod_part.for \
        grid_geom.f \
        hsml.f \
        init_grid.f \
        input.f \
        internal_force \
        kernel.f \
        link_list.f \
        lr_flux.for \
        mod_interf.for \
         mv_interf \
        outflow.for \
        output.f \
        prcorr.for \
        pressure_pd.for \
        puzir.for \
        raley_taylor.for \
        refine.for \
        reigenvalues.for \
        reigenvectors.for \
        riemann.for \
        roe_flux.for \
        shepard.for \
        shifting.for \
        single_step.for \
        sound_speed.for \
        sym_euler.f \
        sym_verlet.for \
        sym_verlet2.for \
        symplectic.for \
        tank.for \
        tank15.for \
        tensile.for \
        time_elapsed.f90 \
        time_integration.f \
        time_integration.f \
        time_print.f90 \
        treedefs.f90 \
        treeload.for \
        treesearch.for \
        verlet.for \
        virt_part.f \
        viscosity.f \
        vtk.for \
        ws_coeff.for

TEST_SRC:=sph.f
# Create lists of the build artefacts in this project
OBJS := $(addsuffix .o, $(SRCS))
TEST_OBJS:=$(addsuffix .o, $(TEST_SRCS))
LIB := $(patsubst %, lib%.a, $(NAME))
TEST_EXE := $(patsubst %.f, %.exe, $(TEST_SRCS))

# Declare all public targets
.PHONY: all clean
all: $(LIB) $(TEST_EXE)

# Create the static library from the object files
$(LIB): $(OBJS)
	$(AR) $@ $^

# Link the test executables
$(TEST_EXE): %.exe: %.f90.o $(LIB)
	$(LD) -o $@ $^

# Create object files from Fortran source
$(OBJS) $(TEST_OBJS): %.o: %
	$(FC) -c -o $@ $<

# Create lists of the build artefacts in this project
OBJS := $(addsuffix .o, $(SRCS))
TEST_OBJS := $(addsuffix .o, $(TEST_SRCS))
LIB := $(patsubst %, lib%.a, $(NAME))
TEST_EXE := $(patsubst %.f90, %.exe, $(TEST_SRCS))

# Declare all public targets
.PHONY: all clean
all: $(LIB) $(TEST_EXE)

# Create the static library from the object files
$(LIB): $(OBJS)
	$(AR) $@ $^

# Link the test executables
$(TEST_EXE): %.exe: %.f90.o $(LIB)
	$(LD) -o $@ $^

# Create object files from Fortran source
$(OBJS) $(TEST_OBJS): %.o: %
	$(FC) -c -o $@ $<

# Define all module interdependencies
ex.mod := ex.f90.o
module.f90.o: $(ex.mod)


# Cleanup, filter to avoid removing source code by accident
clean:
	$(RM) $(filter %.o, $(OBJS) $(TEST_OBJS)) $(filter %.exe, $(TEST_EXE)) $(LIB) $(wildcard *.mod)
        

