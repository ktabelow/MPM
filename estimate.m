% estimate the 4 3D parameters of the ESTATICSmodel of the slice i
% has to run for all i in sdim(3) to get the complete model
%
% =========================================================================
% 2015/09/30
%
% matlab Implementation of the R function estimateESTATICS.r in package qMRI
% written by K. Tabelow and J. Polzehl (?)
%
% from Lars Ruthotto 
% functions loadImageSPM, GaussNewton, FLASHobjFctn, FLASHobjFctn2, 
% Armijo  all written by Lars Ruthotto
% http://www.mathcs.emory.edu/~lruthot/
%
% 
%
% function [res] = estimate(i, tolerance, maski, sdim, dir,t1Files,pdFiles,mtFiles,TE, TEScale, DataScale)
% 
% estimate the 4 3D parameters of the ESTATICSmodel
%
% Input:
%
%  i         - the slice on which the estimation is done
%  tolerance - to stop the GaussNewton algorithm
%  maski     - submask of the level i
%  sdim      - dimension
%  dir       - string with the name of the folder containing the files
%  t1Files   -
%  pdFiles   - input files
%  mtFile    -
%  TE        - echo times vector
%  TEScale   - scaling factor for TE
%  DataScale - scaling factor for data
%
%
%  
%
% Output:
%
%  res - struct containing 
% 
%  modelCoeff -  the coefficients of the ESTATICS model for the level,
%  invCov     -  inverse covariance matrix (unscaled) = J'*J
%                attention!: this matrix is sparse and the nv*nv 
%                single matrices are blocks on the diagonal  
%  sig2       -  vector containing the sigma2 term
%                attention!: this vector already contains nv*nv*nvoxel values
%                 
%  
%
% =========================================================================

function [res] = estimate(i, tolerance, maski,sdim, t1Files,mtFiles,pdFiles,TE, TEScale, DataScale)

if ~isempty(mtFiles)
    nv = 4;
else
    nv = 3;
end
% if the mask on one level is empty doesn't do anything (to speed things up)
if sum(maski(:))>0,
slices = i;
m = [sdim(1), sdim(2), numel(slices)];

T1 = zeros([length(t1Files),m]);
for k=1:length(t1Files),
    [T1(k,:,:,:),omega] = loadImageSPM(fullfile(t1Files{k}),'slices',slices);
end
 
if ~isempty(mtFiles)
MT = zeros([length(mtFiles),m]);
for k=1:length(mtFiles),
    [MT(k,:,:,:),omega] = loadImageSPM(fullfile(mtFiles{k}),'slices',slices);
end
end

PD = zeros([length(pdFiles),m]);
for k=1:length(pdFiles),
    [PD(k,:,:,:),omega] = loadImageSPM(fullfile(pdFiles{k}),'slices',slices);
end

TE = TE./TEScale; %[t1TE(:)./100; mtTE(:)./100; pdTE(:)./100];
%m0 = [1000/DataScale; 1000/DataScale;1000/DataScale; 0.05*TEScale];
if ~isempty(mtFiles)
    m0 = [1000; 1000;1000; 0.05*TEScale];
else
    m0 = [1000; 1000; 0.05*TEScale];
end

%% perform optimization


    % the commented part does the optimisation only on the elements in mask
    % but it ist too slow...
    %T1t = reshape(T1(:,logical(maski)), size(T1,1),[]);
    %MTt = reshape(MT(:,logical(maski)), size(MT,1),[]);
    %PDt = reshape(PD(:,logical(maski)), size(PD,1),[]);
    T1t = reshape(T1, size(T1,1),[]); %
    
    PDt = reshape(PD, size(PD,1),[]); %
    if ~isempty(mtFiles)
        MTt = reshape(MT, size(MT,1),[]); %
        data = [T1t./DataScale; MTt./DataScale; PDt./DataScale];


        indicator = [ones(length(t1Files),1);2*ones(length(mtFiles),1); 3*ones(length(pdFiles),1)];
        fctn = @(model) FLASHObjectiveFunction(model,data,TE,indicator);
    
        % m0big has to be changed if we want a voxel-dipendent starting point
        m0big = kron(ones(size(T1t,2),1),m0); 
        lower = zeros(size(m0big));
        % mi1 = GaussNewton(fctn,m0big); % old gauss without constraints
        mi1 = ProjGaussNewton(fctn,m0big,'tolJ',tolerance,'tolG',tolerance,'tolU',tolerance,'lower',lower,'maxIter',10);
        

        %maski = reshape(maski,numel(maski),[]);
        %maskbig = repelem(maski,4);
        %m1big = zeros(4*prod(m),1);

        % fill in information and set correct dimension
        %m1big(logical(maskbig),1) = mi1(:,1);
        %res.coeff=m1big; %

        res.coeff=mi1;

        %% analyze result

        [~,para,dD,~,Hstar] = fctn(mi1); % ~ is for H, Hstar contains only J'*J
        sig2 = para.Dcs/(size(TE,1)-numel(m0));


        %iind = 1: length(maski);
        %iind = iind(logical(maski));
        %dia = sparse(iind,iind,ones(length(iind),1),sdim(1)*sdim(2),sdim(1)*sdim(2)); %,sdim(1)*sdim(2),sdim(1)*sdim(2)
        %big_indexing = kron(dia, [1 1 1 1 ; 1 1 1 1; 1 1 1 1 ; 1 1 1 1]);
        %small_indexing = kron(speye(numel(iind)),[1 1 1 1 ; 1 1 1 1; 1 1 1 1 ; 1 1 1 1] );
        %small_invCov = H./sig2;
        %res.invCov = spalloc(4*sdim(1)*sdim(2),4*sdim(1)*sdim(2),numel(small_invCov));
        %res.invCov(logical(big_indexing))=small_invCov(logical(small_indexing));
        res.invCov=Hstar;
        res.sig2 = repelem(sig2,nv*nv);
        %res.invCov=J'*J./sig2;
        %res.invCov=H/sig2; %res.invCov=H./sig2;
        res.para=para;
        res.dD=dD;
        res.omega=omega;
    else 
        data = [T1t./DataScale; PDt./DataScale];


        indicator = [ones(length(t1Files),1); 2*ones(length(pdFiles),1)];
        fctn = @(model) FLASHObjectiveFunction_WithoutMT(model,data,TE,indicator);
    
        % m0big has to be changed if we want a voxel-dipendent starting point
        m0big = kron(ones(size(T1t,2),1),m0); 
        lower = zeros(size(m0big));
        % mi1 = GaussNewton(fctn,m0big);
        mi1 = ProjGaussNewton(fctn,m0big,'tolJ',tolerance,'tolG',tolerance,'tolU',tolerance,'lower',lower,'maxIter',10);

        %maski = reshape(maski,numel(maski),[]);
        %maskbig = repelem(maski,4);
        %m1big = zeros(4*prod(m),1);

        % fill in information and set correct dimension
        %m1big(logical(maskbig),1) = mi1(:,1);
        %res.coeff=m1big; %

        res.coeff=mi1;

        %% analyze result

        [~,para,dD,~,Hstar] = fctn(mi1); % ~ is for H, Hstar contains only J'*J
        sig2 = para.Dcs/(size(TE,1)-numel(m0));


        %iind = 1: length(maski);
        %iind = iind(logical(maski));
        %dia = sparse(iind,iind,ones(length(iind),1),sdim(1)*sdim(2),sdim(1)*sdim(2)); %,sdim(1)*sdim(2),sdim(1)*sdim(2)
        %big_indexing = kron(dia, [1 1 1 1 ; 1 1 1 1; 1 1 1 1 ; 1 1 1 1]);
        %small_indexing = kron(speye(numel(iind)),[1 1 1 1 ; 1 1 1 1; 1 1 1 1 ; 1 1 1 1] );
        %small_invCov = H./sig2;
        %res.invCov = spalloc(4*sdim(1)*sdim(2),4*sdim(1)*sdim(2),numel(small_invCov));
        %res.invCov(logical(big_indexing))=small_invCov(logical(small_indexing));

        res.invCov=Hstar; %
        res.sig2 = repelem(sig2,nv*nv);
        res.para=para;
        res.dD=dD;
        res.omega=omega;
    end
else
   
    res.coeff=zeros(nv*sdim(1)*sdim(2),1);
    res.invCov = sparse(nv*sdim(1)*sdim(2),nv*sdim(1)*sdim(2));
    res.sig2=zeros(nv*nv*sdim(1)*sdim(2),1);
end

 end