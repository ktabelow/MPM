
% to compile the fortran files

mex FFLAGS='$FFLAGS' -largeArrayDims gethani.F cgethani.F
mex FFLAGS='$FFLAGS' -largeArrayDims  -v vaws2.F cvaws2.F;
mex FFLAGS='$FFLAGS' -largeArrayDims  -v vpaws2.F cpvaws2.F;
mex FFLAGS='$FFLAGS' -largeArrayDims  -v geticov.F esticov.F;

