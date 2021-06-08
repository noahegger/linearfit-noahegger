COMPILER= gfortran

FLAGS = 

EXEC = nuclear_energies

SRC = $(wildcard *.f90) 

OBJ = $(SRC:.f90=.o)

$(EXEC): $(OBJ)
	$(COMPILER) $(FLAGS) -o $@ $^

types.o: types.f90
	$(COMPILER) $(FLAGS) -c $<

read_write.o: read_write.f90 types.o nuclear_model.o
	$(COMPILER) $(FLAGS) -c $<

linear_algebra.o: linear_algebra.f90 types.o
	$(COMPILER) $(FLAGS) -c $<

nuclear_model.o: nuclear_model.f90 types.o linear_algebra.o
	$(COMPILER) $(FLAGS) -c $<

main.o: main.f90 types.o read_write.o nuclear_model.o
	$(COMPILER) $(FLAGS) -c $<

clean:
	rm -rf *.o *.mod

mrproper: clean
	rm -rf $(EXEC)
