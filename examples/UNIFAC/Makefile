FORT    = gfortran
FOPTS   = -O3 -march=native -ffast-math -funroll-loops -fstrict-aliasing -cpp -Wunused

FORTRIX = ../..
LIBS = -L$(FORTRIX)/lib/ -lfortrix -L/usr/local/cuda/lib64 -lcudart -lcublas

SRCDIR = ./src
OBJDIR = ./obj

all: unifac preunifac

clean:
	rm -rf $(OBJDIR)
	rm -f unifac
	rm -f preunifac

unifac: $(OBJDIR)/unifac.o  $(OBJDIR)/scilab_style.o  $(OBJDIR)/mString.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -o $@ $^ $(LIBS)

$(OBJDIR)/unifac.o: $(SRCDIR)/unifac.f90 $(OBJDIR)/scilab_style.o $(OBJDIR)/mTypes.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -I$(FORTRIX)/include -o $@ $<

$(OBJDIR)/mTypes.o: $(SRCDIR)/mTypes.f90
	mkdir -p $(OBJDIR)
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

preunifac: $(OBJDIR)/preunifac.o $(OBJDIR)/scilab_style.o $(OBJDIR)/mString.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -o $@ $^

$(OBJDIR)/preunifac.o: $(SRCDIR)/preunifac.f90 $(OBJDIR)/scilab_style.o
	$(FORT) $(FOPTS) -I$(OBJDIR) -c -o $@ $<

$(OBJDIR)/scilab_style.o: $(SRCDIR)/scilab_style.f90 $(OBJDIR)/mString.o
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mString.o: $(SRCDIR)/mString.f90
	mkdir -p $(OBJDIR)
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<


