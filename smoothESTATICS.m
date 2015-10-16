% performs adaptive smoothing on the 4 (3) 3D parameters of the ESTATICSmodel
%
% =========================================================================
% 2015/10/6
%
% matlab Implementation of the R function smoothESTATICS.r in package qMRI
% written by K. Tabelow and J. Polzehl (?)
%
% calls the Fortran subroutine vaws.f to perform the smoothing algorithm;
% vaws.f works with OpenMP, in order not to have problems 
% the environment variable OMP_NUM_THREADS has to be set
% TO DO reimplementation of vaws.F in C
%
% [modelS] = smoothESTATICS(model,varargin)
% 
% performs adaptive smoothing on the 4 3D parameters of the ESTATICSmodel
%
% Input:
%
%  model     - ESTATIC model to be smoothed (as struct)
%
% Optional Arguments:
%
%  kstar   - numbers of iteration, default: 16
%  alpha   - , default: 0.05
%  wghts   - weights, default: []
%  verbose - additional information, default: true
%
% Output:
%
%  modelS  - a struct that contains:
%
%  modelCoeff       - smoothed 4 parameters
%  bi               - 
%  smoothPar        - vector containing lambda, hakt and alpha
%
% =========================================================================

function [modelS] = smoothESTATICS(model,varargin)

%% checks if the model in input is present
% it is important to pass the correct struct to this function
 

if nargin==0,
    error('There has to be at least one input, the mpmESTATICSmodel');    
end

%% checks if input model is a struct

if ~isstruct(model) 
    error('Wrong input data type, struct expected'); 
end
% to check if the model struct has at least the requested fields 
if ~isfield(model,{'modelCoeff','invCov','sdim','nFiles','mask','nv'}),
    error('Wrong input data type, struct with fields modelCoeff, invCov, sdim, nFiles, mask, nv expected');    
end

%% sets the default parameters

kstar = 16;
alpha = 0.05;
wghts = [];
verbose = true;

%% overwrites default parameters, if present

for k=1:2:length(varargin),     
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;
%% a little information

if verbose,
    fprintf('verbose = %d \n ', verbose );
    fprintf('kstar = %d \n', kstar );
    fprintf('alpha = %f \n', alpha );
    
end
%% setting parameters

% nv describes the model, with or without mt
% length of the vector to smooth (# parameters of model)
nv = model.nv; 
% determine a suitable adaptation bandwidth
lambda = nv*finv(1-alpha, nv, model.nFiles - nv);

% adjust for non-isotropic voxel if necessary
if isempty(wghts)
    wghts = [1 1 1];   
end
% make first spatial dimension unit, and use second and third for reference
wghts = wghts(1)./wghts(2:3);

% spatial dimension and number of voxel

n1 =  model.sdim(1);
n2 =  model.sdim(2);
n3 =  model.sdim(3);
n = n1*n2*n3;

% initialisation for first step
zobj.bi = ones([1 n]);
zobj.theta = model.modelCoeff;
bi = zobj.bi;

% set the number of usable threads (cores in R)
nthreadsinuse=2*feature('numcores');
setenv('OMP_NUM_THREADS',num2str(nthreadsinuse));
mccores= cast(str2double(getenv('OMP_NUM_THREADS')),'int32');

% define the maximum bandwidth  from the number of iteration
hmax = 1.25^(kstar/3);

if verbose,
    fprintf('Starting the smoothing process ...\n\n');
    msg = sprintf('Percent done: 0.0');
    fprintf(msg);
    protocol=cell(1,kstar);
end

%% perform the iteration
k = 1;

while k<=kstar
   % determine the actual bandwidth for this step
   hakt = gethani (1, 1.25*hmax, 2, 1.25^k, wghts, 1e-4);
    
   % we need the (approx.) size of the weighting scheme array
   tmp = 2*fix(hakt./[1 wghts])+1;
   dlw = tmp(1:3);
   
   % perform the actual adaptive smoothing
   % [zobj.bi, zobj.theta, zobj.hakt] = vaws(reshape(model.modelCoeff,nv,[]), reshape(model.mask,n,1), nv, n1, n2, n3, hakt, lambda, reshape(zobj.theta,nv,[]), reshape(model.invCov,[nv nv n]), zobj.bi, zeros(nv,n), mccores, zeros(prod(dlw),1), wghts, zeros(nv,mccores) );
   [zobj.bi, zobj.theta, zobj.hakt] = vaws(reshape(model.modelCoeff,nv,[]), reshape(model.mask,1,n), nv, n1, n2, n3, hakt, lambda, reshape(zobj.theta,nv,[]), reshape(model.invCov,[nv nv n]), zobj.bi, zeros(nv,n), mccores, zeros(1,prod(dlw)), wghts, zeros(nv,mccores) );
   zobj.bi = max(bi,zobj.bi);
   bi = zobj.bi;
   
   % some verbose stuff
   if verbose
       diff=reshape(zobj.theta,1,[]) - reshape(model.modelCoeff,1,[]);
       diff2=diff.^2;
       diffabs=abs(diff);
       protocol{k} = sprintf('k= %2.d bandwith: %3.3f MSE %3.3d MAE %3.3d mean(bi) %3.3f \n',k,hakt, mean(diff2(:)),mean(diffabs(:)), mean(zobj.bi(:)) );
       % Display the progress
       percentDone = 100 * k / kstar;        
       reverseStr = repmat(sprintf('\b'), 1, length(msg));
       msg = sprintf('Percent done: %3.1f', percentDone);
       fprintf([reverseStr, msg]);
       
   end
   
   %go for the next iteration
    k = k+1;
end


%% some verbose stuff

if verbose
fprintf('\n');
for i=1:numel(protocol)
    fprintf(protocol{i});
end
end

%% preparing the output

% set correct dimension
 zobj.theta = reshape(zobj.theta, [nv n1 n2 n3]);
 zobj.bi = reshape(zobj.bi, [n1 n2 n3]);

% assign values
 modelS = model;
 modelS.modelCoeff = zobj.theta;
 modelS.bi = zobj.bi;
 modelS.smoothPar = [lambda, hakt, alpha]; 
  


end
% end function smoothESTATICS()
