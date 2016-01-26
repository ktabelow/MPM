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
    mtFiles.num     = [1 Inf];
%    mtFiles.val     = {[]};
    
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
    sdim.help    = {'Dimension of the cubus. ' ...
        'If  the vector [0 0 0] is not changed, the spatial dimension will be read from the first t1 file. ' ...
        ' '};
    sdim.strtype = 'e';
    sdim.num     = [1 3];
    % sdim.val    = {[322 368 256]};
    sdim.val    = {[0 0 0]};
    
% ---------------------------------------------------------------------
% tr2 for the mt aquisition
% ---------------------------------------------------------------------
    tr2        = cfg_entry;
    tr2.tag     = 'tr2';
    tr2.name    = 'TR2';
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
    height.help    = {'Height of the level considered in each step. ' ...
                      'Please consider the resources of your computer when changing this number!' ...
                      ' '};
    height.strtype = 'e';
    height.num     = [1 1];
    height.val    = {30};
    
    
% ---------------------------------------------------------------------
% number of iteration of the smoothing algorithm
% ---------------------------------------------------------------------
    kstar        = cfg_entry;
    kstar.tag     = 'kstar';
    kstar.name    = 'kstar';
    kstar.help    = {'Number of iteration of the smoothing algorithm. ' ...
                     'If 0 no smoothing will be performed. ' ...
                     ' '};
    kstar.strtype = 'e';
    kstar.num     = [1 1];
    kstar.val    = {16};
    
% ---------------------------------------------------------------------
% adaptation bandwidth of the smoothing algorithm
% ---------------------------------------------------------------------
    lambda        = cfg_entry;
    lambda.tag     = 'lambda';
    lambda.name    = 'lambda';
    lambda.help    = {'Adaptation bandwidth of the smoothing algorithm' ...
                      ' '};
    lambda.strtype = 'e';
    lambda.num     = [0 Inf];
    lambda.val    = {[]};
    
% ---------------------------------------------------------------------
% precision for the estatics algorithm
% ---------------------------------------------------------------------
    tol        = cfg_entry;
    tol.tag     = 'tol';
    tol.name    = 'tol';
    tol.help    = {'Tolerance level to stop the convergence in the ESTATIC algorithm optimisation.' ...
                   'A smaller value will deliver better results, but will take longer. ' ...
                   'On the contrary a value too big will produce artifacts. ' ...
                   'We advise to leave the default value and repeat with a smaller value, in case of stripes in the output. ' ... 
                   ' '};
    tol.strtype = 'e';
    tol.num     = [1 1];
    tol.val    = {1e-5};
   
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
    maskFile.val     = {[]};
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
    b1File.val     = {[]};
    
% ---------------------------------------------------------------------
% save ESTATICS modell - if it is on, the ESTATICS model parameter and
% relevant data are written
% ---------------------------------------------------------------------
saveESTA   = cfg_menu;
saveESTA.tag     = 'saveESTA';
saveESTA.name    = 'Save the ESTATICS model?';
saveESTA.help    = {'This option enables to save the ESTATICS model. ' ...
                    'If activated, it will save a ESTATICS_model_filename.mat file ' ...
                    'with all the necessary information about the ESTATICS model.' ...
                    'Activating this option will slow down the process and ' ...
                    'requires the user to have enough free memory to save the file (up to a few GB). ' ...
                    'Once the model is saved, it can be used in the toolbox branch' ...
                    ' "Use an existing ESTATICS model" to repeat the smoothing algorithm ' ...
                    'with other kstar and lambda values.' ...
                    ' '};
saveESTA.labels = {
               'No'
               'Yes'
}';
saveESTA.values = {0 1};
saveESTA.val    = {0};
% ---------------------------------------------------------------------
% t1TR 
% ---------------------------------------------------------------------
    t1TR        = cfg_entry;
    t1TR.tag     = 't1TR';
    t1TR.name    = 't1TR';
    t1TR.help    = {'Relaxation times. Please insert a value for each t1File. ' ...
                    'If not specified, the program will try to read it from the .nii file. '};
    t1TR.strtype = 'e';
    t1TR.num     = [0 Inf];
   % t1TR.val    = {[22.5 22.5 22.5 22.5 22.5 22.5]};
    t1TR.val    = {[]};
    
% ---------------------------------------------------------------------
% mtTR 
% ---------------------------------------------------------------------
    mtTR        = cfg_entry;
    mtTR.tag     = 'mtTR';
    mtTR.name    = 'mtTR';
    mtTR.help    = {'Magnetisation transfer. Please insert a value for each mtFile. ' ...
                    'If not specified, the program will try to read it from the .nii file. '};
    mtTR.strtype = 'e';
    mtTR.num     = [0 Inf];
   % mtTR.val    = {[26.1 26.1 26.1 26.1 26.1 26.1]};
    mtTR.val    = {[]};
% ---------------------------------------------------------------------
% pdTR 
% ---------------------------------------------------------------------
    pdTR        = cfg_entry;
    pdTR.tag     = 'pdTR';
    pdTR.name    = 'pdTR';
    pdTR.help    = {'Proton density. Please insert a value for each pdFile. ' ...
                    'If not specified, the program will try to read it from the .nii file. '};
    pdTR.strtype = 'e';
    pdTR.num     = [0 Inf];
   % pdTR.val    = {[22.5 22.5 22.5 22.5 22.5 22.5]};
    pdTR.val    = {[]};
        
% ---------------------------------------------------------------------
% t1TE 
% ---------------------------------------------------------------------
    t1TE        = cfg_entry;
    t1TE.tag     = 't1TE';
    t1TE.name    = 't1TE';
    t1TE.help    = {'Relaxation times echo times. Please insert a value for each t1File. ' ...
        'If not specified, the program will try to read it from the .nii file. '};
    t1TE.strtype = 'e';
    t1TE.num     = [0 Inf];
   % t1TE.val    = {[2.71, 5.17,  7.63, 10.09 , 12.55, 15.01]};
    t1TE.val    = {[]};
    
% ---------------------------------------------------------------------
% mtTE 
% ---------------------------------------------------------------------
    mtTE        = cfg_entry;
    mtTE.tag     = 'mtTE';
    mtTE.name    = 'mtTE';
    mtTE.help    = {'Mt echo times. Please insert a value for each mtFile. ' ...
        'If not specified, the program will try to read it from the .nii file. '};
    mtTE.strtype = 'e';
    mtTE.num     = [0 Inf];
   % mtTE.val    = {[2.71, 5.17,  7.63, 10.09 , 12.55, 15.01]};
    mtTE.val    = {[]};
% ---------------------------------------------------------------------
% pdTE 
% ---------------------------------------------------------------------
    pdTE        = cfg_entry;
    pdTE.tag     = 'pdTE';
    pdTE.name    = 'pdTE';
    pdTE.help    = {'Proton density echo times.Please insert a value for each pdFile. ' ...
        'If not specified, the program will try to read it from the .nii file. '};
    pdTE.strtype = 'e';
    pdTE.num     = [0 Inf];
    % pdTE.val    = {[2.71, 5.17,  7.63, 10.09 , 12.55, 15.01]};
    pdTE.val    = {[]};
% ---------------------------------------------------------------------
% t1FA 
% ---------------------------------------------------------------------
    t1FA        = cfg_entry;
    t1FA.tag     = 't1FA';
    t1FA.name    = 't1FA';
    t1FA.help    = {'Relaxation times flip angle. Please insert a value for each t1File. ' ...
        'If not specified, the program will try to read it from the .nii file. '};
    t1FA.strtype = 'e';
    t1FA.num     = [0 Inf];
   % t1FA.val    = {[27, 27, 27, 27, 27, 27]};
    t1FA.val    = {[]};
    
% ---------------------------------------------------------------------
% mtFA 
% ---------------------------------------------------------------------
    mtFA        = cfg_entry;
    mtFA.tag     = 'mtFA';
    mtFA.name    = 'mtFA';
    mtFA.help    = {'Mt flip angle. Please insert a value for each mtFile. ' ...
        'If not specified, the program will try to read it from the .nii file. '};
    mtFA.strtype = 'e';
    mtFA.num     = [0 Inf];
   % mtFA.val    = {[5,5,5,5,5,5]};
    mtFA.val    = {[]};
% ---------------------------------------------------------------------
% pdFA 
% ---------------------------------------------------------------------
    pdFA        = cfg_entry;
    pdFA.tag     = 'pdFA';
    pdFA.name    = 'pdFA';
    pdFA.help    = {'Proton density flip angle. Please insert a value for each pdFile. ' ...
        'If not specified, the program will try to read it from the .nii file. '};
    pdFA.strtype = 'e';
    pdFA.num     = [0 Inf];
   % pdFA.val    = {[5,5,5,5,5,5]};
    pdFA.val    = {[]};
        
    
% --------------------------------------------------------------------
% .mat file containing an ESTATIC model
% ---------------------------------------------------------------------
    ESTAmodel         = cfg_files;
    ESTAmodel.tag     = 'ESTAmodel';
    ESTAmodel.name    = 'ESTATICS model';
    ESTAmodel.help    = {'Select the .mat file containing the ESTATICS model produced by the mpm toolbox.'};
    ESTAmodel.filter  = '.mat';
    ESTAmodel.ufilter = '.*';
    ESTAmodel.num     = [1 1];
    %ESTAmodel.val     = {[]};
    
% ---------------------------------------------------------------------        
% branch with MT files
% ---------------------------------------------------------------------
     MT_branch         = cfg_exbranch;
     MT_branch.tag     = 'MT_branch';
     MT_branch.name    = 'Model with t1, pd and mt files';
     MT_branch.val     = {t1Files mtFiles pdFiles maskFile tr2 height kstar lambda ...
                          tol b1File saveESTA t1TR mtTR pdTR ...
                          t1TE mtTE pdTE t1FA mtFA pdFA};
     MT_branch.help    = {'This branch implements multi parameter mapping for a dataset with T1,PD and MT files.'}';
     MT_branch.prog    = @spm_local_mpm;
     
% ---------------------------------------------------------------------        
% branch without MT files
% ---------------------------------------------------------------------
     withoutMT_branch         = cfg_exbranch;
     withoutMT_branch.tag     = 'withoutMT_branch';
     withoutMT_branch.name    = 'Model without mt files';
     withoutMT_branch.val     = {t1Files pdFiles maskFile tr2 height kstar lambda ...
                                tol b1File saveESTA t1TR pdTR ... 
                                t1TE pdTE t1FA pdFA}; 
     withoutMT_branch.help    = {'This branch implements multi parameter mapping for a dataset without MT files.'}';
     withoutMT_branch.prog    = @spm_local_mpm_noMT;

% ---------------------------------------------------------------------        
% branch using an existing ESTATICS model
% ---------------------------------------------------------------------
     ESTAmodel_branch         = cfg_exbranch;
     ESTAmodel_branch.tag     = 'ESTAmodel_branch';
     ESTAmodel_branch.name    = 'Use an existing ESTATICS model';
     ESTAmodel_branch.val     = {ESTAmodel height kstar lambda b1File tr2};
     ESTAmodel_branch.help    = {'This branch implements the smoothing and final calculation step of the mpm method given an existing ESTATICS model.'}';
     ESTAmodel_branch.prog    = @spm_local_mpm_givenESTATICS;
     
% ---------------------------------------------------------------------        
% branch producing .nii file from an existing ESTATICS model
% ---------------------------------------------------------------------
     extractNII_branch         = cfg_exbranch;
     extractNII_branch.tag     = 'extractNII_branch';
     extractNII_branch.name    = 'Get nii files from ESTATICS model';
     extractNII_branch.val     = {ESTAmodel height};
     extractNII_branch.help    = {'This branch produces the nii files of the estimated parameter ' ...
         'given an existing ESTATICS model. No smoothing applied or final parameters calculated.'}';
     extractNII_branch.prog    = @spm_local_mpm_extractNII;

% ---------------------------------------------------------------------
%   MPM toolbox
% ---------------------------------------------------------------------
    MPM         = cfg_choice;
    MPM.tag     = 'MPM';
    MPM.name    = 'MPM Multi-parameter Mapping';
    MPM.values     = {MT_branch withoutMT_branch ESTAmodel_branch extractNII_branch};
    MPM.help    = {'This toolbox implements multi parameter mapping for SPM.'}';
%    MPM.prog    = @spm_local_mpm;

%     MPM         = cfg_exbranch;
%     MPM.tag     = 'MPM';
%     MPM.name    = 'MPM';
%     MPM.val     = {t1Files mtFiles pdFiles sdim tr2 height kstar lambda ...
%                    tol maskFile b1File t1TR mtTR pdTR ...
%                    t1TE mtTE pdTE t1FA mtFA pdFA};
%     MPM.help    = {'This toolbox implements multi parameter mapping for SPM.'}';
% %     MPM.prog    = @spm_local_mpm;

%======================================================================
function [] = spm_local_mpm(job)

    
    if ~isdeployed, addpath(fullfile(spm('Dir'),'toolbox','MPM')); end
    job.sdim = [0 0 0];
    mpmESTATICS(job);
    
end

function [] = spm_local_mpm_noMT(job)

    
    if ~isdeployed, addpath(fullfile(spm('Dir'),'toolbox','MPM')); end
    job.mtFiles = [];
    job.mtTR = [];
    job.mtTE = [];
    job.mtFA = [];
    job.sdim = [0 0 0];
    mpmESTATICS(job);
    
end

function [] = spm_local_mpm_givenESTATICS(job)

    
    if ~isdeployed, addpath(fullfile(spm('Dir'),'toolbox','MPM')); end
    
    mpm_givenESTATICS(job);
    
end

function [] = spm_local_mpm_extractNII(job)

    
    if ~isdeployed, addpath(fullfile(spm('Dir'),'toolbox','MPM')); end
    
    mpm_extractNII(job);
    
end
end




