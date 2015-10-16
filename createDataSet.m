% create the data structure to apply the ESTATICSmodel
%
% =========================================================================
% 2015/10/15
%
% 
%
% [dataset] = createDataSet(sdim,t1Files,pdFiles,mtFiles,maskFile,t1TR,pdTR,mtTR,t1TE,pdTE,mtTE,t1FA,pdFA, mtFA);
% 
% validate and use the inputs to assembly the necessary struct 
% to apply the ESTATICSmodel
%
% Input:
%
%  sdim     - a vector containing the 3 dimensions of the cubus
%  dir      - a string containing the folder where the data files are
%  t1Files  - a string cell array with the names of the T1 images files
%  pdFiles  - a string cell array with the names of the PD images files
%  mtFiles  - a string cell array with the names of the MT images files 
%             (can be empty)
%  maskFile - a string cell array with the name of the mask file
%             (can be empty)
%  t1TR       - a vector containing the repetition times for T1
%  pdTR       - a vector containing the repetition times for PD
%  mtTR       - a vector containing the repetition times for MT
%             (can be empty)
%  t1TE       - a vector containing the echo times for T1
%  pdTE       - a vector containing the echo times for PD
%  mtTE       - a vector containing the echo times for MT
%             (can be empty)
%  t1FA       - a vector containing the flip angles for T1
%  pdFA       - a vector containing the flip angles for PD
%  mtFA       - a vector containing the flip angles for MT
%             (can be empty)
%             
% for the function to work, all imput arguments must be passed, even if
% empty
% care must be taken that the t1** vectors have the same length as t1Files,
% same for pd and mt
%
% Output:
%
%  dataset  - a struct that contains:
%
%  sdim     - a vector containing the 3 dimensions of the cubus
%  dir      - a string containing the folder where the data files are
%  t1Files  - a string cell array with the names of the T1 images files
%  pdFiles  - a string cell array with the names of the PD images files
%  mtFiles  - a string cell array with the names of the MT images files 
%             (can be empty)
%  maskFile - a string cell array with the name of the mask file
%             (can be empty)
%  mask     - a 3D matrix with 1 for the voxels to be considered 
%  nv       - number of parameters to estimate 
%             (can be 4 -> complete model
%                  or 3 -> missing MT)
%  nFiles   - number of imput data files 
%  TR       - a vector containing all the repetition times
%  TE       - a vector containing all the echo times
%  FA       - a vector containing all the flip angles
%
% =========================================================================

function [dataset] = createDataSet(sdim,dir,t1Files,pdFiles,mtFiles,maskFile,t1TR,pdTR,mtTR,t1TE,pdTE,mtTE,t1FA,pdFA, mtFA)
  dataset.dir = dir;
  
  if isempty(sdim)
      error('need spatial dimensionality of the data'); 
  end
  if ~isnumeric(sdim)|| length(sdim)~=3
       error('need exactly 3 numbers for spatial dimension'); 
  end
  dataset.sdim = sdim;
  if isempty(t1Files)
      error('cell array of T1 files required'); 
  end
  dataset.t1Files = t1Files;
  if isempty(pdFiles)
      error('cell array of PD files required'); 
  end
  dataset.pdFiles = pdFiles;
  
  if isempty(mtFiles)
      dataset.nv = 3;
  else
      dataset.nv = 4;
      dataset.mtFiles = mtFiles;
  end
  
  dataset.maskFile = maskFile;
  if isempty(maskFile)
      dataset.mask=ones(sdim);
  else 
      slices=1:sdim(3);
      [mask(:,:,:),~] = loadImageSPM(fullfile(dir,[maskFile{1} '.nii']),'slices',slices);
      mask = round(mask(:));
      mask = reshape (mask, sdim);
      dataset.mask=mask;
      clear mask;
  end
  
  dataset.nFiles = length(t1Files)+ length(pdFiles) + length(mtFiles);
  
  if length(t1FA)~=length(t1Files) || length(t1TR)~=length(t1Files) || length(t1TE)~=length(t1Files)
     error('one of the T1 parameter is wrongly dimensioned');
  end
  
  if length(pdFA)~=length(pdFiles) || length(pdTR)~=length(pdFiles) || length(pdTE)~=length(pdFiles)
     error('one of the PD parameter is wrongly dimensioned');
  end
  
  if length(mtFA)~=length(mtFiles) || length(mtTR)~=length(mtFiles) || length(mtTE)~=length(mtFiles)
     error('one of the MT parameter is wrongly dimensioned');
  end
  
  if dataset.nv==4
      dataset.FA = [t1FA(:); mtFA(:); pdFA(:)];
      dataset.TE = [t1TE(:); mtTE(:); pdTE(:)];
      dataset.TR = [t1TR(:); mtTR(:); pdTR(:)];
  else
      dataset.FA = [t1FA(:); pdFA(:)];
      dataset.TE = [t1TE(:); pdTE(:)];
      dataset.TR = [t1TR(:); pdTR(:)];
  end
end