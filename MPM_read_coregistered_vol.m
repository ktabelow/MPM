% read coregistered image volume
%
% written by C. D'Alonzo
% 
% Input:
%   VF     - structure containing image volume information of moved image
%   VG     - structure containing image volume information of reference image 
%
%   optional:  
%
% affineTransMatrix - matrix with the affine transformation, if not present
%                     the function MPM_get_coreg_matrix(VF,VG) is called
% slices            - vector with slices to be read, if not present read the
%                     whole volume
% ========================================================================
function Vol = MPM_read_coregistered_vol(VF,VG,varargin)
     dm      = VG.dim;
    
    slices = [];
    affineTransMatrix = [];
    
    

    for k=1:2:length(varargin),     % overwrites default parameter
         eval([varargin{k},'=varargin{',int2str(k+1),'};']);
    end;

    if isempty(affineTransMatrix), 
        %affineTransMatrix = MPM_get_coreg_matrix(VF,VG); 
        affineTransMatrix = spm_get_space(VF.fname);
    end
    if isempty(slices), slices = 1:dm(3); end;
    
   
    Vol    = zeros([dm(1),dm(2),numel(slices)]);
    
    for p=1:numel(slices),
        matrixP = VG.mat*spm_matrix([0 0 slices(p)]);
        Vol(:,:,p) = spm_slice_vol(VF,affineTransMatrix\matrixP,dm(1:2),1);
    end
       
    
end
