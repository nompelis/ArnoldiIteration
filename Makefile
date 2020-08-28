 FOPTS = -Wall
 FOPTS += -fbounds-check
 FORT = gfortran

all:
	$(FORT) $(FOPTS) code.f

clean:
	rm -f *.o *.mod a.out

