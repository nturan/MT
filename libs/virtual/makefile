THIS=.

COMPILER=gfortran-7
QCDLOOP=../QCDLoop
BSYPPTT=../bsypptt/fortran

OBJECTS=bsyspinors.o bsyspinit.o ampggint.o ampgg.o ampqqint.o ampqq.o bsypptt.o

VPATH=$(BSYPPTT)
ALL: $(OBJECTS);
	ar rcs libBSYPPTT.a $(OBJECTS)
	rm *.o

%.o: %.f;
	$(COMPILER) -c -o  $@ $<

clean: ;
	rm *.a
