# Set TARGET to "cuda" or "blas"
# Set FLOAT to "single" or "double"
TARGET  = cuda
FLOAT   = double

FORT    = gfortran
FOPTS   = -O3 -march=native -ffast-math -funroll-loops -fstrict-aliasing -cpp -D$(TARGET) -D$(FLOAT) -Darch$(NBITS) -Wunused
ifeq ($(TARGET),cuda)
  CC      = nvcc
  CCOPTS  = -O3 -arch sm_20
else
  CC      = g++
  CCOPTS  = -O3
endif

CUDADIR = /usr/lib/x86_64-linux-gnu
BLASDIR = /usr/local/atlas/lib

SRCDIR  = ./src
LIBDIR  = ./lib
INCDIR  = ./include
TSTDIR  = ./testing
OBJDIR  = $(SRCDIR)/obj_$(TARGET)

TSTSRC  = $(wildcard $(TSTDIR)/*.f90)
TSTEXE  = $(TSTSRC:.f90=)

ifeq ($(TARGET),cuda)
  LIBS = -L$(CUDADIR) -lcudart -lcublas
else
  LIBS = -L$(BLASDIR) -lsatlas
endif

all: lib

test: $(TSTEXE)

lib: $(LIBDIR)/libfortrix.a

build: $(OBJDIR)/fortrix.o

clean:
	rm -rf $(OBJDIR)
	rm -rf $(LIBDIR)
	rm -rf $(INCDIR)
	rm -f $(TSTEXE)

clean-all:
	rm -rf $(SRCDIR)/obj_*
	rm -rf $(LIBDIR)
	rm -rf $(INCDIR)
	rm -f $(TSTEXE)

$(TSTDIR)/%: $(TSTDIR)/%.f90 $(TSTDIR)/binary_entrywise_operation.f $(LIBDIR)/libfortrix.a
	$(FORT) $(FOPTS) -o $@ -I$(INCDIR) -I$(OBJDIR) $< -L$(LIBDIR) -lfortrix $(LIBS)

$(LIBDIR)/libfortrix.a: $(OBJDIR)/fortrix.o
	mkdir -p $(LIBDIR)
	ar cr $(LIBDIR)/libfortrix.a $(OBJDIR)/*.o
	mkdir -p $(INCDIR)
	cp -f $(OBJDIR)/fortrix.mod $(INCDIR)/

$(OBJDIR)/fortrix.o: $(SRCDIR)/fortrix.f90 $(OBJDIR)/matrix_module.o $(OBJDIR)/scalar_module.o $(OBJDIR)/dmatrix_module.o $(OBJDIR)/tmatrix_module.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/matrix_module.o: $(SRCDIR)/matrix_module.f90 $(OBJDIR)/base_matrix_module.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/dmatrix_module.o: $(SRCDIR)/dmatrix_module.f90 $(OBJDIR)/base_matrix_module.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/tmatrix_module.o: $(SRCDIR)/tmatrix_module.f90 $(OBJDIR)/base_matrix_module.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/scalar_module.o: $(SRCDIR)/scalar_module.f90 $(OBJDIR)/base_matrix_module.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/base_matrix_module.o: $(SRCDIR)/base_matrix_module.f90 $(OBJDIR)/matrix_handler_$(TARGET).o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/matrix_handler_cuda.o: $(SRCDIR)/matrix_handler_cuda.f90 $(OBJDIR)/matrix_handler.o $(OBJDIR)/cuda_interface.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/matrix_handler_blas.o: $(SRCDIR)/matrix_handler_blas.f90 $(OBJDIR)/matrix_handler.o $(OBJDIR)/arrays.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/cuda_interface.o: $(SRCDIR)/cuda_interface.f90 $(OBJDIR)/cublas_v2_fortran.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/cublas_v2_fortran.o: $(SRCDIR)/cublas_v2_fortran.c $(SRCDIR)/cublas_v2_fortran.h $(OBJDIR)/cuda_operations.o
	$(CC) $(CCOPTS) -c -I$(CUDADIR)/include -o $@ $<

$(OBJDIR)/arrays.o: $(SRCDIR)/arrays.c $(SRCDIR)/arrays.h
	mkdir -p $(OBJDIR)
	$(CC) $(CCOPTS) -c -o $@ $<

$(OBJDIR)/cuda_operations.o: $(SRCDIR)/cuda_operations.cu $(SRCDIR)/cuda_operations.h
	mkdir -p $(OBJDIR)
	$(CC) $(CCOPTS) -c -I$(CUDADIR)/include -o $@ $<

$(OBJDIR)/matrix_handler.o: $(SRCDIR)/matrix_handler.f90
	mkdir -p $(OBJDIR)
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

