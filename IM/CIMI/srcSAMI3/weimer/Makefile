
# Name of the application. 
NAME = ./w05sc.x

# Source files for this application. 
FFILES = read_data.f90\
         w05sc.f90\
         test_w05sc_mod2.f90 

OBJFILES = $(FFILES:.f90=.o)

# Flags to pass on to the compiler and linker.
FFLAGS = -mcmodel=large -vec_report0 -save -O2 -r8

# Compiler
F90 = ifort 

# Creating OBJFILES
%.o: %.f90
	$(F90) $(FFLAGS) -c $*.f90

# Rules

$(NAME): 	$(OBJFILES)
	$(F90) $(FFLAGS) -o $@ $(OBJFILES)

clean:
	rm *.out *.o *.x *.mod wei05sc_*_f90.dat



