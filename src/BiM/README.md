BiM code can be found here:
http://web.math.unifi.it/users/brugnano/BiM/index.html

A few modifications were made to the orginal code so that it could be run on a cluster
using MPI.  Most of the modifications regard output. Also, the LINPACK calls have been 
redirected to LAPACK so that Intel's MKL can be used.  Consequently, calls to BLAS and
other routines in MKL have been removed from subbim.f and renamed here as subbim_NoLinAlg.f
