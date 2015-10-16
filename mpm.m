% script to run the functions of the ESTATICS model
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


mtTE = [2.71,  2.71,  2.71, 5.17,  7.63,10.09];
pdTE = [2.71,  2.71,  2.71, 5.17,  7.63,10.09];
t1TE = [2.71,  2.71,  2.71, 5.17,  7.63,10.09];%, 12.55, 15.01];


mtTR = [26.1, 26.1,26.1,26.1,26.1,26.1];
pdTR = [22.5, 22.5,22.5,22.5,22.5,22.5];
t1TR = [22.5, 22.5,22.5,22.5,22.5,22.5];

% this is the B1 field map coregistered and resliced to t1Files[1]
b1File = { strcat(dir,'data3/B1r1/rsmuB1map_sMQ03089-0002-00001-000001-01.img')};

zStart = 201;
zEnd = 204;

% function [dataset] = createDataSet(sdim,zStart, zEnd, dir,t1Files,pdFiles,mtFiles,maskFile,t1TR,pdTR,mtTR,t1TE,pdTE,mtTE,t1FA,pdFA, mtFA)
dataset = createDataSet(sdim,zStart, zEnd,t1Files,pdFiles,mtFiles,maskFile,t1TR,pdTR,mtTR,t1TE,pdTE,mtTE,t1FA,pdFA, mtFA);

% function [model] = estimateESTATICS(dataset, varargin)
modelMPM3 = estimateESTATICS(dataset);

% function [modelS] = smoothESTATICS(model, varargin)
modelMPM3s = smoothESTATICS(modelMPM3);

% function [qi] = calculateQI(model, varargin)
qi = calculateQI(modelMPM3, 'TR2',3.6,'b1File',b1File);
qiS = calculateQI(modelMPM3s, 'TR2',3.6,'b1File',b1File);