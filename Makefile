F90 = gfortran
F90FLAGS = -c -fbounds-check

PROG = masca
OBJS =   mod_precision.o mod_variables.o mod_maillage.o mod_flux.o mod_calcul.o mod_io.o  masca.o

SUFFIXES = .f90.o
.SUFFIXES : .f90 .o
$(SUFFIXES) :
	$(F90) $(F90FLAGS) $<

$(PROG) : $(OBJS)
	$(F90) $(OBJS) -o $(PROG)

clean:
	@rm -f *.o *.mod *~ a.out core $(PROG) SORTIES/visit/*
	@echo "Repertoire nettoyé"

run :
	make clean
	make
	@ ./$(PROG)
affiche :
	make clean
	@gedit *.f90 Makefile datas &
