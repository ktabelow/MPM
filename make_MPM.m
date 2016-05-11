% to compile the fortran files

mex gethani.F cgethani.F
mex FFLAGS='$FFLAGS -fopenmp' -largeArrayDims vaws2.F cvaws2.F;


