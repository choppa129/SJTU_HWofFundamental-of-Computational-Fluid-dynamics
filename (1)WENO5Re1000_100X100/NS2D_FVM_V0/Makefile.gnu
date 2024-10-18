
export FC = gfortran
export FFLAGS =  \
		-O3 \
		-fdefault-real-8 -fdefault-double-8

# List all application source files
sources := 	global.F90 \
		common_grid.F90 \
		Cart_utility.F90 \
		Cart_vector_trans.F90 \
		Cart_Convection_Diffusion.F90 \
		Cart_bc.F90 \
		Cart_ini.F90 \
		Cart_projection.F90 \
		Cart_adv.F90 \
		Cart_save.F90
%.o:%.F90
	$(FC) -c $(FFLAGS) $<
objects := $(sources:.F90=.o)


# Identify the main program
main := MAIN.F90 
mainobject := $(main:.F90=.o)


# compiles the program.
Cart: $(mainobject) $(objects)
	$(FC) $(FFLAGS)  -o $@ $^ ${PETSC_KSP_LIB}




