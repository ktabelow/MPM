
% to compile the fortran files

mex FFLAGS='$FFLAGS' -compatibleArrayDims gethani.F cgethani.F
mex FFLAGS='$FFLAGS' -compatibleArrayDims vaws2.F cvaws2.F;
mex FFLAGS='$FFLAGS' -compatibleArrayDims vpaws2.F cpvaws2.F;
mex FFLAGS='$FFLAGS' -compatibleArrayDims geticov.F esticov.F;

