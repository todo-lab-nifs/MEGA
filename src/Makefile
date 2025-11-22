FC = mpiifx

# FFLAGS = -O3 -qmkl -xCORE-AVX512 -Bdynamic -fpic -mcmodel=medium -qopenmp -fpp -DMKL -DPRM2048MPI -DPTCL1MPI -DSMP4 -DOPEN_MP
# LDFLAGS = -O3 -qmkl -qopenmp
FFLAGS = -O3 -qmkl -xCORE-AVX512 -Bdynamic -fpic -mcmodel=medium -fpp -DMKL -DPRM1024MPI -DPTCL1MPI -DXEON6900P -DKTI

LDFLAGS = -O3 -qmkl

LIBS = 
OBJS = mega2025_open4.o
BIN  = mega2025_open4_kti-1024mpi.go

.SUFFIXES: .f90 .o

#  default : a.out
all:$(BIN)

$(BIN) : $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

.f90.o:
	$(FC) $(FFLAGS) -c $*.f90

clean :
	rm -f *.o $(BIN)
#	rm -f $(OBJS) $(OBJS:.o=.mod) $(BIN)

