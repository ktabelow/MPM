% to compile the fortran files

mex gethani.F cgethani.F
mex FFLAGS='$FFLAGS' -largeArrayDims vaws2.F cvaws2.F;
mex FFLAGS='$FFLAGS' -largeArrayDims vpaws2.F cpvaws2.F;


