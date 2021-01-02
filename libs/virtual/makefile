THIS=.

COMPILER=gfortran-7
QCDLOOP=~/src/QCDLoop
BSYPPTT=.

OBJECTS=bsyspinors.o bsyspinit.o ampggint.o ampgg.o ampqqint.o ampqq.o bsypptt.o

VPATH=$(BSYPPTT)
ALL: $(OBJECTS);
	ar rcs libBSYPPTT.a $(OBJECTS)
	rm *.o

%.o: %.f;
	$(COMPILER) -std=legacy -c -o  $@ $<

clean: ;
	rm *.a
