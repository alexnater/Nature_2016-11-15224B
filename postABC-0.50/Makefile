CC=		g++
CFLAGS=		-O3
PROG=		postABC
INCLUDES=	-I.
RANDLIBPATH=	./randlib1.3
MSLIBPATH=	./ms_lib
LIBPATH=	-L $(RANDLIBPATH) -L $(MSLIBPATH)
LIBS=		-lrand -lms


postABC: randlib mslib
	@echo "Compiling $@..."
	$(CC) $(CFLAGS) *.cpp $(LIBPATH) -o postABC $(LIBS)
#	rm $(MSLIBPATH)/libms.a $(RANDLIBPATH)/librand.a
	@echo "...done!"


all: $(PROG)



randlib: 
	@echo "Compiling $@..."
	gcc $(CFLAGS) -c $(RANDLIBPATH)/randlib.c $(RANDLIBPATH)/linpack.c $(RANDLIBPATH)/com.c
	ar -rcs $(RANDLIBPATH)/librand.a randlib.o linpack.o com.o
	rm randlib.o linpack.o com.o
	@echo "...done!"



mslib:
	@echo "Compiling $@..."
	gcc $(CFLAGS) -c $(MSLIBPATH)/ms.c $(MSLIBPATH)/streec.c $(MSLIBPATH)/rand1.c
	ar -rcs $(MSLIBPATH)/libms.a ms.o streec.o rand1.o
	rm ms.o streec.o rand1.o
	@echo "...done!"

