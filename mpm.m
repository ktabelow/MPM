% script to run the functions of the ESTATICS model

%% setting variables
sdim = [322 368 256];
dir = '/Home/stat/tabelow/projects/MPM/';

t1Files = {strcat(dir,'data3/T1r1/sMQ03089-0005-00001-000256-01.img'), ...
           strcat(dir,'data3/T1r1/sMQ03089-0005-00001-000512-02.img'), ...
           strcat(dir,'data3/T1r1/sMQ03089-0005-00001-000768-03.img'), ...
           strcat(dir,'data3/T1r1/sMQ03089-0005-00001-001024-04.img'), ... 
           strcat(dir,'data3/T1r1/sMQ03089-0005-00001-001280-05.img'), ... 
           strcat(dir,'data3/T1r1/sMQ03089-0005-00001-001536-06.img')};
       
mtFiles = {strcat(dir,'data3/MTr1/sMQ03089-0007-00001-000256-01.img'), ... 
           strcat(dir,'data3/MTr1/sMQ03089-0007-00001-000512-02.img'), ... 
           strcat(dir,'data3/MTr1/sMQ03089-0007-00001-000768-03.img'), ...                                                                           
           strcat(dir,'data3/MTr1/sMQ03089-0007-00001-001024-04.img'), ...                                                                           
          strcat(dir,'data3/MTr1/sMQ03089-0007-00001-001280-05.img'), ...                                                                           
           strcat(dir,'data3/MTr1/sMQ03089-0007-00001-001536-06.img')};
pdFiles = {strcat(dir,'data3/PDr1/sMQ03089-0006-00001-000256-01.img'), ...                                                                           
           strcat(dir,'data3/PDr1/sMQ03089-0006-00001-000512-02.img'), ...                                                                           
           strcat(dir,'data3/PDr1/sMQ03089-0006-00001-000768-03.img'), ...                                                                           
           strcat(dir,'data3/PDr1/sMQ03089-0006-00001-001024-04.img'), ...                                                                           
           strcat(dir,'data3/PDr1/sMQ03089-0006-00001-001280-05.img'), ...                                                                           
           strcat(dir,'data3/PDr1/sMQ03089-0006-00001-001536-06.img')};
maskFile =  strcat(dir,'data3/MSK_c1sMQ03089-0007-00001-000256-01_olsq_MT.nii');  

    
mtFA = [5, 5, 5, 5, 5, 5];
pdFA = [5, 5, 5, 5, 5, 5];
t1FA = [27, 27, 27, 27, 27, 27];


mtTE = [2.71, 5.17,  7.63, 10.09 , 12.55, 15.01];
pdTE = [2.71, 5.17,  7.63, 10.09 , 12.55, 15.01];
t1TE = [2.71, 5.17,  7.63, 10.09 , 12.55, 15.01];


mtTR = [26.1, 26.1,26.1,26.1,26.1,26.1];
pdTR = [22.5, 22.5,22.5,22.5,22.5,22.5];
t1TR = [22.5, 22.5,22.5,22.5,22.5,22.5];

% this is the B1 field map coregistered and resliced to t1Files[1]
b1File = { strcat(dir,'data3/B1r1/rsmuB1map_sMQ03089-0002-00001-000001-01.img')};

% looking for weights
try 
    wghts=getWeights(t1Files{1});
catch ME
wghts = [];
end

%% works on the all cubus, levels by levels
% defines how many planes in each level
height = 30;
% calculates how big (high) the overlapping has to be
% to assure a good smoothing
kstar = 16;
hmax = 1.25^(kstar/3);
hakt = gethani (1, 1.25*hmax, 2, 1.25^kstar, [1 1], 1e-4);
hdelta = ceil(hakt); % height of half of the overlapping
% calculates the interval between the zStart values of the different levels
interval = int64(height - 2*hdelta);

% preparing result variables
R1 = zeros(sdim);
R2star = zeros(sdim);
PD = zeros(sdim);
delta = zeros(sdim);

% prepares variable to save the whole mask
totalmask = zeros(sdim);


% start iteration on all the levels
for startLayerVoxel = 1:interval:sdim(3),
    
    zStart = double(startLayerVoxel);
    %fprintf('Starting at %d \n',zStart);
    
    if sdim(3)-(startLayerVoxel + height)> 2*hdelta,
        % in case the next starting point has enough planes after it
        zEnd = double(startLayerVoxel + height); 
    else
        % in case the next starting point has NOT enough planes after it
        zEnd = sdim(3);
        % startLayerVoxel = sdim(3)+1;
    end
    %fprintf('Ending at %d \n',zEnd);
    
% function [dataset] = createDataSet(sdim,zStart, zEnd, dir,t1Files,pdFiles,mtFiles,maskFile,t1TR,pdTR,mtTR,t1TE,pdTE,mtTE,t1FA,pdFA, mtFA)
dataset = createDataSet(sdim,zStart, zEnd,t1Files,pdFiles,mtFiles,maskFile,t1TR,pdTR,mtTR,t1TE,pdTE,mtTE,t1FA,pdFA, mtFA);

% function [model] = estimateESTATICS(dataset, varargin)
modelMPM3 = estimateESTATICS(dataset);

% function [modelS] = smoothESTATICS(model, varargin)
% modelMPM3s = smoothESTATICS(modelMPM3);

% new function smoothing only on the elements in mask
modelMPM3snew = smoothESTATICSmask(modelMPM3, 'wghts', wghts);

% function [qi] = calculateQI(model, varargin)
%qi = calculateQI(modelMPM3, 'TR2',3.6,'b1File',b1File);
%qiS = calculateQI(modelMPM3s, 'TR2',3.6,'b1File',b1File);

qiSnew = calculateQI(modelMPM3snew, 'TR2',3.6,'b1File',b1File,'verbose', false);

if zStart==1
R1(:,:,zStart:zEnd) = qiSnew.R1;
R2star(:,:,zStart:zEnd) = qiSnew.R2star;
PD(:,:,zStart:zEnd) = qiSnew.PD;
delta(:,:,zStart:zEnd) = qiSnew.delta;
totalmask(:,:,zStart:zEnd) = qiSnew.model.mask;
else 
    R1(:,:,zStart+hdelta:zEnd) = qiSnew.R1(:,:, 1+hdelta: (zEnd-zStart+1) );
R2star(:,:,zStart+hdelta:zEnd) = qiSnew.R2star(:,:, 1+hdelta: (zEnd-zStart+1) );
PD(:,:,zStart+hdelta:zEnd) = qiSnew.PD(:,:, 1+hdelta: (zEnd-zStart+1) );
delta(:,:,zStart+hdelta:zEnd) = qiSnew.delta(:,:, 1+hdelta: (zEnd-zStart+1) );
totalmask(:,:,zStart+hdelta:zEnd) = qiSnew.model.mask(:,:, 1+hdelta: (zEnd-zStart+1) );
end
if zEnd==sdim(3),
    break;
end

end

R2star(isnan(R2star))=0;
PD(isnan(PD))=0;
delta(isnan(delta))=0;
R1(isnan(R1))=0;
 % set all values of the voxels outside the mask to 0
    R1(totalmask<1)=0;
    R2star(totalmask<1)=0;
    PD(totalmask<1)=0;
    delta(totalmask<1)=0; 
