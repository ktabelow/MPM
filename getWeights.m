function [wghts] = getWeights(imgFile)
V = spm_vol(imgFile);
wghts     = sqrt(sum(V(1).mat(1:3,1:3).^2)); % [yvoxel/xvoxel zvoxel/xvoxel]
% wghts   = [vxg(2)/vxg(1) vxg(3)/vxg(1)];
end