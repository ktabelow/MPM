% read data and estimate the 4 (or 3) 3D parameters of the ESTATICSmodel
%
% =========================================================================
% 2015/10/6
%
% matlab Implementation of the R function estimateESTATICS.r in package qMRI
% written by K. Tabelow and J. Polzehl (?)
%
% This implementation is based on the model of Lars Ruthotto to read the 
% data files and estimate the paramaters level by level and not on each
% voxel. This requires less memory and less execution time.
%
% [model] = estimateESTATICS(dataset, varargin)
% 
% reads the data from the dataset and estimates the 4 3D parameters 
% of the ESTATICSmodel using estimate.m
%
% Input:
%
%  dataset     - a struct containing the informations of the model
% 
% These informations must be contained: 
%  sdim     - a vector containing the 3 dimensions of the cubus
%  zStart   - start value on the z axis of the interest volume
%  zEnd     - end value on the z axis of the interest volume
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
% Optional Arguments:
%
%  TEScale   - factor to scale the echo times, default: 100
%  dataScale - factor to scale the data, default: 1000
%  verbose   - additional information, default: true
%  tolerance - to stop the GaussNewton algorithm
%
% Output:
%
%  model  - a struct that contains all information of the dataset and:
%
%  TEScale   - factor to scale the echo times, default: 100
%  DataScale - factor to scale the data, default: 1000
%  modelCoeff- 4 (or 3) 3D parameters of the ESTATICS model
%  invCov    - a 5D matrix containing for each voxel the inverse covariance
%              matrix (4x4 or 3x3)
%
% =========================================================================

function [model] = estimateESTATICS(dataset, varargin)
%% checks if the dataset in input is present
% it is important to pass the correct struct to this function
 

if nargin==0,
    error('There has to be at least one dataset input');    
end

%% checks if input dataset is a struct

if ~isstruct(dataset) 
    error('Wrong input data type, struct expected'); 
end
% to check if the model struct has at least the requested fields 
if ~isfield(dataset,{'sdim','zStart','zEnd','dir','nFiles','nv','t1Files','pdFiles','mtFiles', 'TE'}),
    error('Wrong input data type, struct with fields sdim, zStart, zEnd, dir, nFiles,t1Files,pdFiles,mtFiles, mask, nv, TE expected');    
end

%% sets the default parameters

TEScale = 100;
DataScale = 1000;
verbose = true;
tolerance = 1e-5;

%% overwrites default parameters, if present

for k=1:2:length(varargin),     
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;
%% a little information

if verbose,
    fprintf('   verbose = %d \n', verbose );    
    fprintf('   tolerance = %f \n', tolerance );
    fprintf('   zStart = %d \n', dataset.zStart );
    fprintf('   zEnd = %d \n', dataset.zEnd );
end

%% getting parameters from dataset

nv = dataset.nv;
t1Files = dataset.t1Files;
pdFiles = dataset.pdFiles;
if nv==4, 
    mtFiles = dataset.mtFiles;
else
    mtFiles=[];
end
    
sdim = dataset.sdim;
TE = dataset.TE;
zStart = dataset.zStart;
zEnd = dataset.zEnd;

%% preparing for reading and estimate
% until now implemented only for levels containing one layer
% numberLayers=1;
coeff=zeros(nv*sdim(1)*sdim(2),(zEnd-zStart+1));
sigma2 = zeros(nv*nv*sdim(1)*sdim(2),(zEnd-zStart+1));%
invCov=cell((zEnd-zStart+1)/1,1); %sdim(3)

if verbose,
    fprintf('Starting the ESTATICS model at %s \n',datestr(now));
    msg = sprintf('Percent done: 0.0');
    fprintf(msg);    
end

% iterating on all the layers between zStart and zEnd
% when zStart = 1 and zEnd = sdim(3) the all cubus is done
for k=1:(zEnd-zStart+1) %1:sdim(3)
    
    i= zStart+k-1;
    maski = dataset.mask(:,:,k);
  res=estimate(i,tolerance,maski,sdim,t1Files,mtFiles,pdFiles,TE,TEScale,DataScale,dataset);
%  res=estimate(i,sdim,t1Files,mtFiles,pdFiles,TE,TEScale,DataScale);
  coeff(:,k)=res.coeff;  
  invCov{k}=res.invCov;
  sigma2(:,k)=res.sig2;
  if verbose
  % Display the progress
       percentDone = 100 * k / (zEnd-zStart+1); %sdim(3) ;        
       reverseStr = repmat(sprintf('\b'), 1, length(msg));
       msg = sprintf('Percent done: %3.1f', percentDone);
       fprintf([reverseStr, msg]);
  end
end


coeff = reshape (coeff, [nv sdim(1) sdim(2) (zEnd-zStart+1)]);
sigma2 = reshape (sigma2, [nv nv sdim(1) sdim(2) (zEnd-zStart+1)]);
if verbose
    fprintf('\n');
end

if verbose
fprintf('getting the invCov matrix: \n');
end

invC = zeros([nv nv sdim(1) sdim(2) (zEnd-zStart+1)]);
vec = zeros(1,nv*nv*sdim(1)*sdim(2));

sizeV = nv*sdim(1)*sdim(2);
for i = 1: nv*sdim(1)*sdim(2)
    j = i-1;
   vec(j*nv+1) = j*sizeV+1 +nv*fix(j/nv); 
   vec(j*nv+2) = j*sizeV+2 +nv*fix(j/nv); 
   vec(j*nv+3) = j*sizeV+3 +nv*fix(j/nv); 
   if nv==4
       vec(j*nv+4) = j*sizeV+4 +nv*fix(j/nv); 
   end
end
clear sizeV

for k=1:length(invCov)  
        invC(:,:,:,:,k) = reshape(full(invCov{k}(vec)),[nv nv sdim(1) sdim(2)]);        
end


if verbose,
fprintf('finished at: %s \n',datestr(now));
end

model = dataset;
model.modelCoeff = coeff;
model.invCov = invC./sigma2;
model.TEScale = TEScale;
model.DataScale = DataScale;



end