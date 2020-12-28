THIS=.

COMPILER=gfortran-7

OBJECTS=gg_ttb.o \
 uub_ttb.o \
 ubu_ttb.o \
 gg_ttbg.o \
 uub_ttbg.o \
 ubu_ttbg.o \
 gu_ttbu.o \
 gub_ttbub.o \
 ug_ttbu.o \
 ubg_ttbub.o

ALL: $(OBJECTS);
	ar rcs libMatrix.a ../helas204/libdhelas.a $(OBJECTS)
	rm *.o

%.o: %.f;
	$(COMPILER) -c -o  $@ $<

clean: ;
	rm *.o *.a
