% calculate R1, delta, PD and R2star of a model
%
% =========================================================================
% 2015/10/7
%
% matlab Implementation of the R function calculateQI.r in package qMRI
% written by K. Tabelow and J. Polzehl (?)
%
% 
%
% [qi] = calculateQI(model,varargin)
% 
% calculate longitudinal relaxation rate R1, effective proton density PD and magnetisation
% transfer saturation delta in a given model
%
% Input:
%
%  model     - model to be evaluate (as struct)
%  
% Optional Arguments:
%
%  b1File   - filename of the correction field
%  TR2      - repetition time R2, default: 0
%  verbose  - additional information, default: true
%
% Output:
%
%  qi  - a struct that contains:
%
%  b1Map       - correction map
%  R1          - longitudinal relaxation rate
%  R2star      - apparent transverse relaxation rate
%  PD          - effective proton density
%  delta       - magnetisation transfer saturation
%
% =========================================================================


function [qi]= calculateQI(model,varargin)
%% checks if the model in input is present
% it is important to pass the correct struct to this function
 

if nargin<3,
    error('There have to be two inputs, a model and a correction file');    
end

%% checks if input model is a struct

if ~isstruct(model) 
    error('Wrong input data type, struct expected'); 
end
% to check if the model struct has the requested fields (to adjust)
if ~isfield(model,{'modelCoeff','sdim', 'zStart','zEnd'}),
    error('Wrong input data type, struct with at least fields modelCoeff, sdim , zStart, zEnd expected');    
end

%% sets the default parameters
b1File={};
TR2 = 0;
verbose = true;

%% overwrites default parameters, if present

for k=1:2:length(varargin),     
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

%% read B1 correction field

slices = model.zStart : model.zEnd; %1:model.sdim(3);
if ~isempty(b1File) || (length(b1File)==1 && strcmp(b1File{1},''))
   if verbose,fprintf('reading correction file from %s \n',b1File{1}); end   
  [b1Map(:,:,:),~] = loadImageSPM(fullfile(b1File{1}),'slices',slices);
  b1Map = b1Map./100; % still to check if is ok
  b1Map(b1Map < 0) = 0;
else 
    if verbose, fprintf('no B1 correction\n'); end
    b1Map = ones([model.sdim(1) model.sdim(2) (model.zEnd-model.zStart+1)]);
end

%% get correct flip angles an TR times

t1FA = model.FA(1);
if model.nv==4, 
    pdFA = model.FA(length(model.t1Files)+length(model.mtFiles)+1); 
else
    pdFA = model.FA(length(model.t1Files)+1);
end
    
t1TR = model.TR(1);



%% calculate E1

if verbose, fprintf('calculating R1...'); end
alphat1 = b1Map.*t1FA/180*pi; 
alphapd = b1Map.*pdFA/180*pi;
SINalphat1 = sin(alphat1);
COSalphat1 = cos(alphat1);
SINalphapd = sin(alphapd);
COSalphapd = cos(alphapd);

clear alphat1 alphapd;

if model.nv == 4 
 enum = reshape(model.modelCoeff(1,:,:,:),[model.sdim(1) model.sdim(2) (model.zEnd-model.zStart+1)]) - SINalphat1./SINalphapd.* reshape(model.modelCoeff(3,:,:,:),[model.sdim(1) model.sdim(2) (model.zEnd-model.zStart+1)]);
 denom = reshape(model.modelCoeff(1,:,:,:),[model.sdim(1) model.sdim(2) (model.zEnd-model.zStart+1)]).*COSalphat1 - SINalphat1./SINalphapd.* reshape(model.modelCoeff(3,:,:,:),[model.sdim(1) model.sdim(2) (model.zEnd-model.zStart+1)]).* COSalphapd;
else 
    enum = reshape(model.modelCoeff(1,:,:,:),[model.sdim(1) model.sdim(2) (model.zEnd-model.zStart+1)]) - SINalphat1./SINalphapd.* reshape(model.modelCoeff(2,:,:,:),[model.sdim(1) model.sdim(2) (model.zEnd-model.zStart+1)]);
 denom = reshape(model.modelCoeff(1,:,:,:),[model.sdim(1) model.sdim(2) (model.zEnd-model.zStart+1)]).*COSalphat1 - SINalphat1./SINalphapd.* reshape(model.modelCoeff(2,:,:,:),[model.sdim(1) model.sdim(2) (model.zEnd-model.zStart+1)]).* COSalphapd;
end

E1 = enum./denom;

clear enum denom COSalphapd SINalphapd;

R1 = -log(E1)/t1TR;

if verbose, fprintf('done\n'); end


%% calculate PD

if verbose, fprintf('calculating PD...'); end
enum = (1 - COSalphat1 .* E1).*reshape(model.modelCoeff(1,:,:,:),[model.sdim(1) model.sdim(2) (model.zEnd-model.zStart+1)]).*model.DataScale;
denom = SINalphat1.*(1-E1);

PD =  enum./denom;

clear enum denom SINalphat1;

if verbose, fprintf('done\n'); end

%% calculate delta
if model.nv == 4 
if verbose, fprintf('calculating MT...'); end
mtFA = model.FA(length(model.t1Files)+1);
mtTR = model.TR(length(model.t1Files)+1);
alphamt = b1Map .* mtFA / 180 * pi;
E1mt = E1.^(mtTR/t1TR);
E2mt = E1.^(TR2/t1TR);
enum = reshape(model.modelCoeff(2,:,:,:),[model.sdim(1) model.sdim(2) (model.zEnd-model.zStart+1)]) .* model.DataScale - (1-E2mt).*sin(alphamt) .*PD; 
denom = reshape(model.modelCoeff(2,:,:,:),[model.sdim(1) model.sdim(2) (model.zEnd-model.zStart+1)]) .* model.DataScale .* cos(alphamt) .*E1mt + PD.* (E2mt-E1mt).*sin(alphamt);

delta = 1 - enum./denom;
clear alphamt enum denom
else
    delta = [];
end



if verbose, fprintf('done\n'); end
%% prepare output
qi.model = model;
qi.b1Map = b1Map;
qi.R1 = R1;
if model.nv==4
    qi.R2star = reshape(model.modelCoeff(4,:,:,:),[model.sdim(1) model.sdim(2) (model.zEnd-model.zStart+1)])/model.TEScale; 
else
    qi.R2star = reshape(model.modelCoeff(3,:,:,:),[model.sdim(1) model.sdim(2) (model.zEnd-model.zStart+1)])/model.TEScale; %teScale
end
qi.PD = PD;
qi.delta = delta;

end