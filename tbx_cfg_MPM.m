function MPM = tbx_cfg_MPM
% Configuration file for toolbox MPM
%_______________________________________________________________________
% Copyright (C) 2015 WIAS Berlin

% C.D'Alonzo from K.Tabelow configuration file for aws4SMP
%
% ---------------------------------------------------------------------
% spm_file t1Files
% ---------------------------------------------------------------------
    t1Files         = cfg_files;
    t1Files.tag     = 't1Files';
    t1Files.name    = 't1Files';
    t1Files.help    = {'Select the T1 files.'};
    t1Files.filter  = 'image';
    t1Files.ufilter = '.*';
    t1Files.num     = [1 Inf];
%    --------------------------------------------------------------------
% spm_file mtFiles
% ---------------------------------------------------------------------
    mtFiles         = cfg_files;
    mtFiles.tag     = 'mtFiles';
    mtFiles.name    = 'mtFiles';
    mtFiles.help    = {'Select the MT files.'};
    mtFiles.filter  = 'image';
    mtFiles.ufilter = '.*';
    mtFiles.num     = [0 Inf];
% --------------------------------------------------------------------
% spm_file pdFiles
% ---------------------------------------------------------------------
    pdFiles         = cfg_files;
    pdFiles.tag     = 'pdFiles';
    pdFiles.name    = 'pdFiles';
    pdFiles.help    = {'Select the PD files.'};
    pdFiles.filter  = 'image';
    pdFiles.ufilter = '.*';
    pdFiles.num     = [1 Inf];
    

% ---------------------------------------------------------------------
% sdim spatial dimensions
% ---------------------------------------------------------------------
    sdim        = cfg_entry;
    sdim.tag     = 'sdim';
    sdim.name    = 'sdim';
    sdim.help    = {'Dimension of the cubus. '};
    sdim.strtype = 'e';
    sdim.num     = [1 3];
    sdim.val    = {[322 368 256]};
    
% ---------------------------------------------------------------------
% tr2 for the mt aquisition
% ---------------------------------------------------------------------
    tr2        = cfg_entry;
    tr2.tag     = 'tr2';
    tr2.name    = 'tr2';
    tr2.help    = {'TR2>0 for the MT acquisition'};
    tr2.strtype = 'e';
    tr2.num     = [1 1];
    tr2.val    = {3.6};
    
% ---------------------------------------------------------------------
% height of the interest level (especially for workstation)
% ---------------------------------------------------------------------
    height        = cfg_entry;
    height.tag     = 'height';
    height.name    = 'height';
    height.help    = {'Height of the level considered in each step. Please consider the resources of your computer when changing this number!'};
    height.strtype = 'e';
    height.num     = [1 1];
    height.val    = {30};
    
    
% ---------------------------------------------------------------------
% number of iteration of the smoothing algorithm
% ---------------------------------------------------------------------
    kstar        = cfg_entry;
    kstar.tag     = 'kstar';
    kstar.name    = 'kstar';
    kstar.help    = {'Number of iteration of the smoothing algorithm'};
    kstar.strtype = 'e';
    kstar.num     = [1 1];
    kstar.val    = {16};
    
% ---------------------------------------------------------------------
% zStart especially for workstation
% ---------------------------------------------------------------------
%    zStart        = cfg_entry;
%    zStart.tag     = 'zStart';
%    zStart.name    = 'zStart';
%    zStart.help    = {'Start level of the volume of interest '};
%    zStart.strtype = 'e';
%    zStart.num     = [1 1];
%    zStart.val    = {201};
    
% ---------------------------------------------------------------------
% zEnd especially for workstation
% ---------------------------------------------------------------------
%    zEnd        = cfg_entry;
%    zEnd.tag     = 'zEnd';
%    zEnd.name    = 'zEnd';
%    zEnd.help    = {'End level of the volume of interest '};
%    zEnd.strtype = 'e';
%    zEnd.num     = [1 1];
%    zEnd.val    = {225};

% --------------------------------------------------------------------
% spm_file maskFile
% ---------------------------------------------------------------------
    maskFile         = cfg_files;
    maskFile.tag     = 'maskFile';
    maskFile.name    = 'maskFile';
    maskFile.help    = {'Select the mask.'};
    maskFile.filter  = 'image';
    maskFile.ufilter = '.*';
    maskFile.num     = [0 1];
% --------------------------------------------------------------------
% spm_file b1File
% ---------------------------------------------------------------------
    b1File         = cfg_files;
    b1File.tag     = 'b1File';
    b1File.name    = 'b1File';
    b1File.help    = {'Select the correction field file.'};
    b1File.filter  = 'image';
    b1File.ufilter = '.*';
    b1File.num     = [0 1];
    
% ---------------------------------------------------------------------
% t1TR 
% ---------------------------------------------------------------------
    t1TR        = cfg_entry;
    t1TR.tag     = 't1TR';
    t1TR.name    = 't1TR';
    t1TR.help    = {'relaxation times '};
    t1TR.strtype = 'e';
    t1TR.num     = [1 Inf];
    t1TR.val    = {[22.5 22.5 22.5 22.5 22.5 22.5]};
    
% ---------------------------------------------------------------------
% mtTR 
% ---------------------------------------------------------------------
    mtTR        = cfg_entry;
    mtTR.tag     = 'mtTR';
    mtTR.name    = 'mtTR';
    mtTR.help    = {'magnetisation transfer '};
    mtTR.strtype = 'e';
    mtTR.num     = [1 Inf];
    mtTR.val    = {[26.1 26.1 26.1 26.1 26.1 26.1]};
% ---------------------------------------------------------------------
% pdTR 
% ---------------------------------------------------------------------
    pdTR        = cfg_entry;
    pdTR.tag     = 'pdTR';
    pdTR.name    = 'pdTR';
    pdTR.help    = {'proton density '};
    pdTR.strtype = 'e';
    pdTR.num     = [1 Inf];
    pdTR.val    = {[22.5 22.5 22.5 22.5 22.5 22.5]};
        
% ---------------------------------------------------------------------
% t1TE 
% ---------------------------------------------------------------------
    t1TE        = cfg_entry;
    t1TE.tag     = 't1TE';
    t1TE.name    = 't1TE';
    t1TE.help    = {'relaxation times echo times '};
    t1TE.strtype = 'e';
    t1TE.num     = [1 Inf];
    t1TE.val    = {[2.71, 5.17,  7.63, 10.09 , 12.55, 15.01]};
    
% ---------------------------------------------------------------------
% mtTE 
% ---------------------------------------------------------------------
    mtTE        = cfg_entry;
    mtTE.tag     = 'mtTE';
    mtTE.name    = 'mtTE';
    mtTE.help    = {'mt echo times '};
    mtTE.strtype = 'e';
    mtTE.num     = [1 Inf];
    mtTE.val    = {[2.71, 5.17,  7.63, 10.09 , 12.55, 15.01]};
% ---------------------------------------------------------------------
% pdTE 
% ---------------------------------------------------------------------
    pdTE        = cfg_entry;
    pdTE.tag     = 'pdTE';
    pdTE.name    = 'pdTE';
    pdTE.help    = {'proton density echo times '};
    pdTE.strtype = 'e';
    pdTE.num     = [1 Inf];
    pdTE.val    = {[2.71, 5.17,  7.63, 10.09 , 12.55, 15.01]};
% ---------------------------------------------------------------------
% t1FA 
% ---------------------------------------------------------------------
    t1FA        = cfg_entry;
    t1FA.tag     = 't1FA';
    t1FA.name    = 't1FA';
    t1FA.help    = {'relaxation times flip angle '};
    t1FA.strtype = 'e';
    t1FA.num     = [1 Inf];
    t1FA.val    = {[27, 27, 27, 27, 27, 27]};
    
% ---------------------------------------------------------------------
% mtFA 
% ---------------------------------------------------------------------
    mtFA        = cfg_entry;
    mtFA.tag     = 'mtFA';
    mtFA.name    = 'mtFA';
    mtFA.help    = {'mt flip angle '};
    mtFA.strtype = 'e';
    mtFA.num     = [1 Inf];
    mtFA.val    = {[5,5,5,5,5,5]};
% ---------------------------------------------------------------------
% pdFA 
% ---------------------------------------------------------------------
    pdFA        = cfg_entry;
    pdFA.tag     = 'pdFA';
    pdFA.name    = 'pdFA';
    pdFA.help    = {'proton density  flip angle '};
    pdFA.strtype = 'e';
    pdFA.num     = [1 Inf];
    pdFA.val    = {[5,5,5,5,5,5]};
        
        

% ---------------------------------------------------------------------
% fwhm Bias FWHM
% ---------------------------------------------------------------------
%     kernel       = cfg_menu;
%     kernel.tag     = 'kernel';
%     kernel.name    = 'Location kernel';
%     kernel.help    = {'Choice of location kernel'};
%     kernel.val     = {2};
%     kernel.labels = {'Uniform'
%                      'Epanechnikov'
%                      'Biweight'
%                      'Triweight'               
%                      'Gauss'}';
%     kernel.values = {0 1 2 3 4};

% ---------------------------------------------------------------------
% 
% ---------------------------------------------------------------------
    MPM         = cfg_exbranch;
    MPM.tag     = 'MPM';
    MPM.name    = 'MPM';
    MPM.val     = {t1Files mtFiles pdFiles sdim tr2 height kstar ...
                   maskFile b1File t1TR mtTR pdTR ...
                   t1TE mtTE pdTE t1FA mtFA pdFA};
    MPM.help    = {'This toolbox implements multi parameter mapping for SPM.'}';
    MPM.prog    = @spm_local_mpm;

%======================================================================
function [] = spm_local_mpm(job)

    
    if ~isdeployed, addpath(fullfile(spm('Dir'),'toolbox','MPM')); end
   
    mpmESTATICS(job);
    
end





end




