% get affine transformation matrix
% 
% Input:
%   VF     - structure containing image volume information of moved image
%   VG     - structure containing image volume information of reference image 
%

function affineTransMatrix = MPM_get_coreg_matrix(VF,VG)
    
    coregflags.sep = [4 2];
    coregflags.graphics = 0;
    x = spm_coreg(VG,VF, coregflags);
    M = spm_matrix(x);
    MM= spm_get_space(VF.fname);
    affineTransMatrix = M\MM;
end