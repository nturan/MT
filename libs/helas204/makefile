#$Modified: Thu Oct  7 23:04:55 2010 by uwer $

include settings 
include makefile.arch 

VPATH = ../../..


AR = ar
ARFLAGS = r
RANLIB = ranlib

# object files

OBJECTS = \
	   jggxxx.o \
           gggxxx.o \
           jgggxx.o \
           ggggxx.o \
           boostx.o \
           coup1x.o \
           coup2x.o \
           coup3x.o \
           coup4x.o \
           fsixxx.o \
           fsoxxx.o \
           fvixxx.o \
           fvoxxx.o \
           hioxxx.o \
           hvsxxx.o \
           hvvxxx.o \
	   iosxxx.o \
           iovxxx.o \
           sxxxxx.o \
           ixxxxx.o \
           j3xxxx.o \
           jioxxx.o \
           jtioxx.o \
           jvsxxx.o \
           jvvxxx.o \
           jw3wxx.o \
           jwwwxx.o \
	   mom2cx.o \
           oxxxxx.o \
           rotxxx.o \
           vssxxx.o \
           vvsxxx.o \
           vvvxxx.o \
           vxxxxx.o \
           w3w3xx.o \
           wwwwxx.o \
           eaixxx.o  hsssxx.o  jeexxx.o  jvssxx.o  ssssxx.o  vvssxx.o \
           eaoxxx.o  hssxxx.o  jssxxx.o  momntx.o  sssxxx.o


# define target

TARGET = libdhelas.a

# rules


all: 
	make $(TARGET)
	make clean

$(TARGET): $(OBJECTS)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJECTS)
	$(RANLIB) $(TARGET)

clean:
	rm -f $(OBJECTS)

distclean: clean
	rm -f $(TARGET)

%.o: %.F
	$(FC) -c $(FFLAGS) $<

# dependencies

# end
