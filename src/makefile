# Simple makefile to compile MPI version of AutoHetDet_0D
# If you have the intel fortran compiler installed, runnin make in src directory should 
# produced the executable AutoHetDet_0D_Opt
# This script doesn't clean up any files, such as *.mod
#
# To run with MPI, use something like > mpiexec -machinefile ../mpi.hosts ../AutoHetDet_0D
#
# Adding this environmental variables (can be done either on command or setting the environment)
# can allow processes at MPI barriers to go idle
#
# mpiexec -genv I_MPI_THREAD_YIELD 3 -genv I_MPI_THREAD_SLEEP 10 -machinefile ../88core5Nodes.hosts ../AutoHetDet_0D_Opt
#
# Values of 100 for I_MPI_THREAD_SLEEP also work, but setting it to 1000 can slow things down.
#
# Intel Fortran MPI compilier
FC = mpiifort

# Note, -heap-arrays was needed for Windows compilation, not sure if needed here.
# FFLAGSCPU = -O3 -mkl -xHost -qopenmp -fpp -D__MPI__ -ipo -o Darwin-MEP-1D
# -qopt-matmul (this flag forces matrix multiplication, but is default for O3 (not sure about O2)
FFLAGSCPU = -O3 -mkl -xHost -qopenmp -ipo -qopt-matmul -o AutoHetDet_0D_Opt

# using BiM:
BIMDIR = ./BiM
# Module source files must be compiled in a specific order because of inter dependencies
BIMSRC  = $(BIMDIR)/bim__NoWrite_SOLOUTmod.for
BIMSRC += $(BIMDIR)/subbim_NoLinAlg.for

# using hyperBOB
BOBDIR = ./hyperBOB
BOBSRC  = $(BOBDIR)/kind_module.f90
BOBSRC += $(BOBDIR)/bobyqa.f90
BOBSRC += $(BOBDIR)/hyperBOB_V1.4.f90

# using latin hypercube
LATSRC = ./latin_random/latin_random.f90

# LINPACK to LAPACK
LIN2LAPSRC = ./linpack2lapack/Linpack2Lapack.f90

# Main program source files
DRVSRC = ./Module_Thermo_V5.1.f90 
DRVSRC += ./AutoHetDet_opt_V4.7.f90

# note, add back  $(LIN2LAPSRC) if returning to MKL version
CPU: $(BIMSRC) $(LIN2LAPSRC) $(BOBSRC) $(LATSRC) $(DRVSRC)
	$(FC) $(FFLAGSCPU) $(BIMSRC) $(LIN2LAPSRC) $(BOBSRC) $(LATSRC) $(DRVSRC)



