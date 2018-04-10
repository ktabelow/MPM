
% to compile the fortran files

mex FFLAGS='$FFLAGS' -largeArrayDims gethani.F cgethani.F
mex FFLAGS='$FFLAGS' -largeArrayDims vaws2.F cvaws2.F;
mex FFLAGS='$FFLAGS' -largeArrayDims vpaws2.F cpvaws2.F;


