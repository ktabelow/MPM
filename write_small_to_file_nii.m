% write a nii file of a big volume, with data only from a small contained  volume
% 
function []= write_small_to_file_nii(outputdir,fileprefix, big_volume,small_volume_data,zStart, zEnd, sdim)
[pth, nam, ext] = spm_fileparts(big_volume.fname);
big_volume.fname = fullfile(outputdir,nam);
data = zeros(sdim);
data(:,:,zStart:zEnd) = small_volume_data;

my_write_vol_nii(data, big_volume(1), fileprefix);

end