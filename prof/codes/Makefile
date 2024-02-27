CC = mpic++

RM = rm -f 
LAPACK_PATH = /opt/lapack/
BLAS_PATH = /opt/blas
GSL_PATH = /opt/gsl
MPI_PATH = /opt/mpi/gcc-4.8.5/openmpi-4.1.0
ARPACK_PATH = /home/arago/local/arpackpp

LDFLAGS =  -fopenmp -mcmodel=medium -L$(LAPACK_PATH)/lib  -llapack -L$(BLAS_PATH)/lib/  -lrefblas -L$(GSL_PATH)/lib -L$(GSL_PATH)/lib -L$(ARPACK_PATH)/lib  -lsuperlu  -lopenblas -lgfortran -lgsl -lgslcblas  -larpack 
LDFLAGS_POST =  -fopenmp -mcmodel=medium -L$(LAPACK_PATH)/lib  -llapack -L$(BLAS_PATH)/lib/  -lrefblas -L$(GSL_PATH)/lib -L$(GSL_PATH)/lib -L$(ARPACK_PATH)/lib   -lgfortran   -lgsl -lgslcblas
#LDFLAGS =  -fopenmp -mcmodel=medium -L$(LAPACK_PATH)/lib  -llapack -L$(BLAS_PATH)/lib/  -lblas -L$(GSL_PATH)/lib -L$(GSL_PATH)/lib -L$(ARPACK_PATH)/lib  -lgfortran -lgsl -lgslcblas  -larpack 

#CFLAGS = -fgnu89-inline -c -g -O2 -Wall -mcmodel=medium -fopenmp -I${HOME}/local/include/arpackpp
#CFLAGS = -c -g -O2 -std=c99 -Wall -mcmodel=medium -fopenmp -I${HOME}/local/include/arpackpp
CFLAGS = -c -g -O2 -Wall -Wno-restrict -Wno-maybe-uninitialized -Wno-misleading-indentation -Wno-format-overflow -Wno-write-strings -mcmodel=medium -fopenmp  -I${ARPACK_PATH}/include -I$(GSL_PATH)/include/ 
BUILD = t2 vert local static
OBJS = t2.o matrix.o basic.o hamil.o ls_basis.o arp.o

.PHONY: all clean  

default: t2
all: $(BUILD) 
clean: ; $(RM) *.o 

t2: $(OBJS)
	$(CC) $^    $(LDLIBS) $(LDFLAGS) -o $@
vert: vert_new.o matrix.o basic.o
	$(CC) $^    $(LDLIBS) $(LDFLAGS) -o vert_new
genG: genG.o matrix.o basic.o ls_basis.o hamil.o
	$(CC) $^    $(LDLIBS) $(LDFLAGS) -o $@
local: local.o matrix.o basic.o ls_basis.o hamil.o
	$(CC) $^    $(LDLIBS) $(LDFLAGS) -o $@

static: static.o matrix.o basic.o operators.o
	$(CC) $^    $(LDFLAGS) $(LDLIBS) -o $@


t2.o:	$($@:.o=.c)	lattice.h matrix.h matrix.c inverse.c ls_basis.h ls_basis.c basic.c ls.c hop.c memcifbathm.c ran2.c hamil.h hamil.c arp.h arp.c Makefile 
ls_basis.o:	$($@:.o=.c)	lattice.h matrix.h matrix.c ls_basis.h ls_basis.c Makefile 
arp.o:		$($@:.o=.c)	Makefile
hamil.o:	$($@:.o=.c)	Makefile
basic.o:	$($@:.o=.c)	Makefile
operators.o:	$($@:.o=.c)	Makefile
matrix.o:	$($@:.o=.c)	matrix.h	Makefile
static.o:	$($@:.o=.c)	ls_basis.c basic.h basic.c lattice.h operators.c	Makefile
vert_new.o:	$($@:.o=.c)	matrix.h matrix.c lattice.h inverse.c basic.c hop.c ls.c	Makefile
local.o:	$($@:.o=.c)	matrix.h matrix.c lattice.h inverse.c basic.c hop.c ls.c	Makefile
genG.o:		$($@:.o=.c)	matrix.h matrix.c lattice.h inverse.c basic.c hop.c ls.c	Makefile
local_wsqbroad.o:	$($@:.o=.c)	matrix.h matrix.c lattice.h inverse.c basic.c hop.c ls.c	Makefile

.c.o:
	$(CC) $(CFLAGS) -o $@ $< 
