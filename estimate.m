% estimate the 4 3D parameters of the ESTATICSmodel of the slice i
% has to run for all i in sdim(3) to get the complete model
%
% =========================================================================
% 
% last modified 2016/03/17
%
% written by C. D'Alonzo 
%
% 
% functions loadImageSPM, GaussNewton, FLASHObjectiveFunction, 
% written by Lars Ruthotto http://www.mathcs.emory.edu/~lruthot/
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

function [res] = estimate(i, tolerance, maski,sdim, t1Files,mtFiles,pdFiles,TE, TEScale, DataScale, dataset)

if ~isempty(mtFiles)
    nv = 4;
else
    nv = 3;
end
% if the mask on one level is empty doesn't do anything (to speed things up)
if sum(maski(:))>0,
slices = i;
m = [sdim(1), sdim(2), numel(slices)];

if isfield(dataset,'t1_aTMat'),
    T1 = zeros([length(t1Files),m]);
    for k=1:length(t1Files),
        %[T1(k,:,:,:),omega] = loadImageSPM(fullfile(t1Files{k}),'slices',slices);
        %T1(k,:,:,:) = MPM_read_coregistered_vol(spm_vol(fullfile(t1Files{k})),spm_vol(dataset.Pdmean{1}),'slices',slices,'affineTransMatrix',dataset.t1_aTMat(:,:,k));
        T1(k,:,:,:) = MPM_read_coregistered_vol(spm_vol(fullfile(t1Files{k})),spm_vol(dataset.Pdmean{1}),'slices',slices,'affineTransMatrix',dataset.t1_aTMat);

    end
 
    if ~isempty(mtFiles)
        MT = zeros([length(mtFiles),m]);
        for k=1:length(mtFiles),
            %[MT(k,:,:,:),omega] = loadImageSPM(fullfile(mtFiles{k}),'slices',slices);
            MT(k,:,:,:) = MPM_read_coregistered_vol(spm_vol(fullfile(mtFiles{k})),spm_vol(dataset.Pdmean{1}),'slices',slices,'affineTransMatrix',dataset.mt_aTMat);
        end
    end

    PD = zeros([length(pdFiles),m]);
    for k=1:length(pdFiles),
        %[PD(k,:,:,:),omega] = loadImageSPM(fullfile(pdFiles{k}),'slices',slices);
        PD(k,:,:,:) = MPM_read_coregistered_vol(spm_vol(fullfile(pdFiles{k})),spm_vol(dataset.Pdmean{1}),'slices',slices,'affineTransMatrix',dataset.pd_aTMat);
    end
else
    T1 = zeros([length(t1Files),m]);
    for k=1:length(t1Files),
        %[T1(k,:,:,:),omega] = loadImageSPM(fullfile(t1Files{k}),'slices',slices);
        T1(k,:,:,:) = MPM_read_coregistered_vol(spm_vol(fullfile(t1Files{k})),spm_vol(dataset.Pdmean{1}),'slices',slices);
    end
 
    if ~isempty(mtFiles)
        MT = zeros([length(mtFiles),m]);
        for k=1:length(mtFiles),
            %[MT(k,:,:,:),omega] = loadImageSPM(fullfile(mtFiles{k}),'slices',slices);
            MT(k,:,:,:) = MPM_read_coregistered_vol(spm_vol(fullfile(mtFiles{k})),spm_vol(dataset.Pdmean{1}),'slices',slices);
        end
    end

    PD = zeros([length(pdFiles),m]);
    for k=1:length(pdFiles),
        %[PD(k,:,:,:),omega] = loadImageSPM(fullfile(pdFiles{k}),'slices',slices);
        PD(k,:,:,:) = MPM_read_coregistered_vol(spm_vol(fullfile(pdFiles{k})),spm_vol(dataset.Pdmean{1}),'slices',slices);
    end
end

TE = TE./TEScale; 
%m0 = [1000/DataScale; 1000/DataScale;1000/DataScale; 0.05*TEScale];
if ~isempty(mtFiles)
    %m0 = [1000; 1000;1000; 0.05*TEScale];
    m0 = [1000/DataScale; 1000/DataScale;1000/DataScale; 0.05*TEScale];
else
    %m0 = [1000; 1000; 0.05*TEScale];
    m0 = [1000/DataScale; 1000/DataScale; 0.05*TEScale];
end

%% perform optimization


    % the optimisation only on the elements in mask was too slow and was
    % eliminated
    
    T1t = reshape(T1, size(T1,1),[]); 
    
    PDt = reshape(PD, size(PD,1),[]); 
    if ~isempty(mtFiles)
        MTt = reshape(MT, size(MT,1),[]); 
        
        data = [T1t./DataScale; MTt./DataScale; PDt./DataScale];

        dataMasked = data(:,maski==1);
        indicator = [ones(length(t1Files),1);2*ones(length(mtFiles),1); 3*ones(length(pdFiles),1)];
        
        fctn = @(model,credit) FLASHObjectiveFunction(model,dataMasked,TE,indicator,credit);
    
        % m0big has to be changed if we want a voxel-dipendent starting point
        m0big = kron(ones(numel(find(maski)),1),m0); 
        lower = zeros(size(m0big));
        
        mi1 = zeros(numel(m0),numel(maski));
        mt = ProjGaussNewton(fctn,m0big,'verbose',0,'tolJ',tolerance,'tolG',tolerance,'tolU',tolerance,'lower',lower,'maxIter',50);
        mi1(:,maski==1) = mt;
        mi1 = mi1(:);
        
        

        res.coeff=mi1;

        %% analyze result

        
        [~,para,dD,~,Hstar] = FLASHObjectiveFunction(mi1,data,TE,indicator,[]);  % ~ is for H, Hstar contains only J'*J
        sig2 = para.Dcs/(size(TE,1)-numel(m0));

        res.invCov=Hstar;
        %res.sig2 = repelem(sig2,nv*nv);% unluckily is repelem not
        %available for matlab version older then 2015a
       
        %res.sig2 = kron(sig2,ones(1,nv*nv));
        res.sig2 = sig2;
        
        res.para=para;
        res.dD=dD;
        %res.omega=omega;
    else 
        data = [T1t./DataScale; PDt./DataScale];


        indicator = [ones(length(t1Files),1); 2*ones(length(pdFiles),1)];
        % fctn = @(model) FLASHObjectiveFunction_WithoutMT(model,data,TE,indicator);
    
        % m0big has to be changed if we want a voxel-dipendent starting point
        m0big = kron(ones(size(T1t,2),1),m0); 
        lower = zeros(size(m0big));
       

        fctn2 = @(model,credit) FLASHObjectiveFunction_WithoutMT(model,data,TE,indicator,credit);
        mi1 = ProjGaussNewton(fctn2,m0big,'tolJ',tolerance,'tolG',tolerance,'tolU',tolerance,'lower',lower,'maxIter',50,'verbose',0);
        mi1 = mi1(:);

        res.coeff=mi1;

        %% analyze result

        
        [~,para,dD,~,Hstar] = fctn2(mi1,[]);  % ~ is for H, Hstar contains only J'*J
        sig2 = para.Dcs/(size(TE,1)-numel(m0));



        res.invCov=Hstar; %
        %res.sig2 = repelem(sig2,nv*nv); % unluckily is repelem not
        %available for matlab version older then 2015a
        
        %res.sig2 = kron(sig2,ones(1,nv*nv));
        res.sig2 = sig2;
        res.para=para;
        res.dD=dD;
        %res.omega=omega;
    end
else
   
    res.coeff=zeros(nv*sdim(1)*sdim(2),1);
    res.invCov = sparse(nv*sdim(1)*sdim(2),nv*sdim(1)*sdim(2));
    %res.sig2=zeros(nv*nv*sdim(1)*sdim(2),1);
    res.sig2=zeros(sdim(1)*sdim(2),1);
end

 end