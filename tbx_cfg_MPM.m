function MPM = tbx_cfg_MPM
% Configuration file for toolbox MPM
%_______________________________________________________________________
% Copyright (C) 2015 WIAS Berlin

% Author: Chiara D'Alonzo 
% Maintainer: Karsten Tabelow 
% '
% ---------------------------------------------------------------------
% spm_file t1Files
% ---------------------------------------------------------------------
    t1Files         = cfg_files;
    t1Files.tag     = 't1Files';
    t1Files.name    = 'T1w volumes';
    t1Files.help    = {'Select the T1w volumes of the MPM sequence.'};
    t1Files.filter  = 'image';
    t1Files.ufilter = '.*';
    t1Files.num     = [1 Inf];
%    --------------------------------------------------------------------
% spm_file mtFiles
% ---------------------------------------------------------------------
    mtFiles         = cfg_files;
    mtFiles.tag     = 'mtFiles';
    mtFiles.name    = 'MTw volumes';
    mtFiles.help    = {'Select the MTw volumes of the MPM sequence.' ...
                       'If no MTw volumes were acquired in the sequence, choose the option - Model without MTw volumes - from the toolbox menu instead.'};
    mtFiles.filter  = 'image';
    mtFiles.ufilter = '.*';
    mtFiles.num     = [1 Inf];
%    mtFiles.val     = {[]};
    
% --------------------------------------------------------------------
% spm_file pdFiles
% ---------------------------------------------------------------------
    pdFiles         = cfg_files;
    pdFiles.tag     = 'pdFiles';
    pdFiles.name    = 'PDw volumes';
    pdFiles.help    = {'Select the PDw volumes of the MPM sequence.'};
    pdFiles.filter  = 'image';
    pdFiles.ufilter = '.*';
    pdFiles.num     = [1 Inf];
    
% ---------------------------------------------------------------------
% odir Output Directory
% ---------------------------------------------------------------------
    odir         = cfg_files;
    odir.tag     = 'odir';
    odir.name    = 'Output Directory';
    odir.help    = {'Select the directory where the output files should be written.'};
    odir.filter = 'dir';
    odir.ufilter = '.*';
    odir.num     = [1 1];
% ---------------------------------------------------------------------
% sdim spatial dimensions
% ---------------------------------------------------------------------
%     sdim         = cfg_entry;
%     sdim.tag     = 'sdim';
%     sdim.name    = 'dim';
%     sdim.help    = {'Spatial dimension of the data volumes.' ...
%                     'If  the default vector [0 0 0] is not changed, the spatial dimension will be defined by the first T1w volumes.'};
%     sdim.strtype = 'e';
%     sdim.num     = [1 3];
%     sdim.val     = {[0 0 0]};
    
% ---------------------------------------------------------------------
% tr2 for the mt aquisition
% ---------------------------------------------------------------------
    tr2         = cfg_entry;
    tr2.tag     = 'tr2';
    tr2.name    = 'TR2';
    tr2.help    = {'Type the second repetition time TR2 > 0 for the MTw acquisitions.'};
    tr2.strtype = 'e';
    tr2.num     = [1 1];
    tr2.val     = {3.6};
    
% ---------------------------------------------------------------------
% height of the interest level (especially for workstation)
% ---------------------------------------------------------------------
    height         = cfg_entry;
    height.tag     = 'height';
    height.name    = 'Number of slices processed';
    height.help    = {'Type a number of slices that are to be processed at once.' ...
                      'This option saves memory and is necessary when your computer has only very limited hardware resources. This comes at the cost of computation time, as smoothing is then performed on overlapping slices more then once to guarantee correct 3D results. If memory is no problem on your computer (e.g. on a large compute custer), you might be able to process the full volumes at once. However, typically 16GB memory are not much for volumes of spatial dimension 256x256x256 ...' ...
                      ' '};
    height.strtype = 'e';
    height.num     = [1 1];
    height.val     = {30};
    
    
% ---------------------------------------------------------------------
% number of iteration of the smoothing algorithm
% ---------------------------------------------------------------------
    kstar         = cfg_entry;
    kstar.tag     = 'kstar';
    kstar.name    = 'Number of iterations (kstar)';
    kstar.help    = {'Number of iteration of the smoothing algorithm.' ...
                     'If 0 is given, no smoothing will be performed.' ...
                     'A larger number of iterations leads to more smoothing in regions of homogeneous intensity. If smooth trends are in the data (which always are), the result for large values of kstar is a step-function and gives a cartoon-like impression of the images.' ...
                     'Smaller values of kstar lead to lower amount of noise reduction and reduce the cartoon-like effects.' ...
                     'For a non-adaptive version of the smoothing method (see parameter lamdba) kstar directly corresponds to the bandwidth of the non-adaptive kernel.'};
    kstar.strtype = 'e';
    kstar.num     = [1 1];
    kstar.val     = {16};
    
% ---------------------------------------------------------------------
% adaptation bandwidth of the smoothing algorithm
% ---------------------------------------------------------------------
    lambda         = cfg_entry;
    lambda.tag     = 'lambda';
    lambda.name    = 'Adaptation bandwidth of the smoothing (lambda)';
    lambda.help    = {'Adaptation bandwidth of the smoothing algorithm' ...
                      'The parameter steers the amount of adaptivity of the method. If lambda is 0 (zero) the method will not smooth at all and the resulting images are unchanged. If lambda is very large, say 1e20, then the procedure is non-adaptive and is equivalent to a non-adaptive kernel smoother with a bandwidth that corresponds to the number of iterations (kstar). A proper choice of lambda should (in theory) only depend on the type of noise (but not the SNR). In practise, values of about 10-20 give good results.'};
    lambda.strtype = 'e';
    lambda.num     = [0 Inf];
    lambda.val     = {20};
    
% ---------------------------------------------------------------------
% precision for the estatics algorithm
% ---------------------------------------------------------------------
    tol         = cfg_entry;
    tol.tag     = 'tol';
    tol.name    = 'tol';
    tol.help    = {'Tolerance level to stop the convergence in the ESTATIC algorithm optimization.' ...
                   'A smaller value will deliver better results, but will take longer.' ...
                   'Contrary, large values might result in artifacts like intensity stripes along the slices in the final R2*, R1, PD images.' ...
                   'We advise to leave the default value and repeat the computation with a smaller value, in case you see such stripes in the output files.'};
    tol.strtype = 'e';
    tol.num     = [1 1];
    tol.val     = {1e-5};
   
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
    maskFile.name    = 'Brain mask';
    maskFile.help    = {'Select a mask file.'
                        'If no file is selected, every voxel is considered in the calculations.'};
    maskFile.filter  = 'image';
    maskFile.ufilter = '.*';
    maskFile.num     = [0 1];
    maskFile.val     = {[]};

% --------------------------------------------------------------------
% spm_file b1File
% ---------------------------------------------------------------------
%     b1File         = cfg_files;
%     b1File.tag     = 'b1File';
%     b1File.name    = 'B1 correction file';
%     b1File.help    = {'Select the correction field file if available' ...
%                       'This volume has to come with the same spatial dimension as the MPM data. If no volume is selected, no correction is performed.'};
%     b1File.filter  = 'image';
%     b1File.ufilter = '.*';
%     b1File.num     = [0 1];
%     b1File.val     = {[]};
    
% --------------------------------------------------------------------
% spm_file b1File
% ---------------------------------------------------------------------
    b1FileA         = cfg_files;
    b1FileA.tag     = 'b1FileA';
    b1FileA.name    = 'B1 correction file - Amplitude Image';
    b1FileA.help    = {'Select the amplitude image file.' ...
                      'If no volume is selected, no correction is performed.'};
    b1FileA.filter  = 'image';
    b1FileA.ufilter = '.*';
    b1FileA.num     = [0 1];
    b1FileA.val     = {[]};
    
% --------------------------------------------------------------------
% spm_file b1File
% ---------------------------------------------------------------------
    b1FileP         = cfg_files;
    b1FileP.tag     = 'b1FileP';
    b1FileP.name    = 'B1 correction file - Phase Image';
    b1FileP.help    = {'Select the phase image file.' ...
                      'If no volume is selected, no correction is performed.'};
    b1FileP.filter  = 'image';
    b1FileP.ufilter = '.*';
    b1FileP.num     = [0 1];
    b1FileP.val     = {[]};
    
% ---------------------------------------------------------------------
% save ESTATICS modell - if it is on, the ESTATICS model parameter and
% relevant data are written
% ---------------------------------------------------------------------
    saveESTA         = cfg_menu;
    saveESTA.tag     = 'saveESTA';
    saveESTA.name    = 'Save the ESTATICS model?';
    saveESTA.help    = {'This option enables to save the ESTATICS model. ' ...
                        'If activated, the ESTATICS parameters are saved to a ESTATICS_model_filename.mat file with all the necessary information about the ESTATICS model.' ...
                        'Activating this option will slow down the computation and requires enough disk space to save the file (up to some GB, as the covariance matrix has to be saved as well.).' ...
                        'Once the model is saved, it can be used in the toolbox branch "Use an existing ESTATICS model" to repeat the smoothing algorithm with other adaptation bandwidth values and number of iterations.' ...
                      ' '};
    saveESTA.labels  = {'No'
                       'Yes'};
    saveESTA.values  = {0 1};
    saveESTA.val     = {0};
    
% ---------------------------------------------------------------------
% coregister images - if it is on, all the calculation are done on the
% coregistered images; if off, the matrix returned from spm_get_space is
% used
% ---------------------------------------------------------------------
    coregIM         = cfg_menu;
    coregIM.tag     = 'coregIM';
    coregIM.name    = 'Coregister the images?';
    coregIM.help    = {'This option enables to perform the coregistration of the images. Use it if your images have not been coregistered, but consider the additional time required.' ...
                      ' '};
    coregIM.labels  = {'No'
                       'Yes'};
    coregIM.values  = {0 1};
    coregIM.val     = {0};

% ---------------------------------------------------------------------
% t1TR 
% ---------------------------------------------------------------------
    t1TR         = cfg_entry;
    t1TR.tag     = 't1TR';
    t1TR.name    = 'TR values for T1 weighted volumes';
    t1TR.help    = {'Repetition times of the T1w volumes. Please insert a value for each T1 weighted volume. ' ...
                    'If not specified, the toolbox will try to read it from the description field in the header of the image files. ' ...
                    'This description has to contain a part like TR=22.5ms/TE=2.71ms/FA=27deg'};
    t1TR.strtype = 'e';
    t1TR.num     = [0 Inf];
    t1TR.val     = {[]};
    
% ---------------------------------------------------------------------
% mtTR 
% ---------------------------------------------------------------------
    mtTR         = cfg_entry;
    mtTR.tag     = 'mtTR';
    mtTR.name    = 'TR values for MT weighted volumes';
    mtTR.help    = {'Repetition times of the MTw volumes. Please insert a value for each MT weighted volume. ' ...
                    'If not specified, the toolbox will try to read it from the description field in the header of the image files. ' ...
                    'This description has to contain a part like TR=22.5ms/TE=2.71ms/FA=27deg'};
    mtTR.strtype = 'e';
    mtTR.num     = [0 Inf];
    mtTR.val     = {[]};
% ---------------------------------------------------------------------
% pdTR 
% ---------------------------------------------------------------------
    pdTR         = cfg_entry;
    pdTR.tag     = 'pdTR';
    pdTR.name    = 'TR values for PD weighted volumes';
    pdTR.help    = {'Repetition times of the PDw volumes. Please insert a value for each PD weighted volume.  ' ...
                    'If not specified, the toolbox will try to read it from the description field in the header of the image files. ' ...
                    'This description has to contain a part like TR=22.5ms/TE=2.71ms/FA=27deg'};
    pdTR.strtype = 'e';
    pdTR.num     = [0 Inf];
    pdTR.val     = {[]};
        
% ---------------------------------------------------------------------
% t1TE 
% ---------------------------------------------------------------------
    t1TE         = cfg_entry;
    t1TE.tag     = 't1TE';
    t1TE.name    = 'TE values for the T1w volumes';
    t1TE.help    = {'Echo times for the sequence of T1w images. Please insert a value for each T1 weighted volume. ' ...
                    'If not specified, the toolbox will try to read it from the description field in the header of the image files.' ...
                    'This description has to contain a part like TR=22.5ms/TE=2.71ms/FA=27deg'};
    t1TE.strtype = 'e';
    t1TE.num     = [0 Inf];
    t1TE.val     = {[]};
    
% ---------------------------------------------------------------------
% mtTE 
% ---------------------------------------------------------------------
    mtTE         = cfg_entry;
    mtTE.tag     = 'mtTE';
    mtTE.name    = 'TE values for the MTw volumes';
    mtTE.help    = {'Echo times for the sequence of MTw images. Please insert a value for each MT weighted volume. ' ...
                    'If not specified, the toolbox will try to read it from the description field in the header of the image files. ' ...
                    'This description has to contain a part like TR=22.5ms/TE=2.71ms/FA=27deg'};
    mtTE.strtype = 'e';
    mtTE.num     = [0 Inf];
    mtTE.val     = {[]};

% ---------------------------------------------------------------------
% pdTE 
% ---------------------------------------------------------------------
    pdTE         = cfg_entry;
    pdTE.tag     = 'pdTE';
    pdTE.name    = 'TE values for the PDw volumes';
    pdTE.help    = {'Echo times for the sequence of PDw images. Please insert a value for each PD weighted volume. ' ...
                    'If not specified, the toolbox will try to read it from the description field in the header of the image files. ' ...
                    'This description has to contain a part like TR=22.5ms/TE=2.71ms/FA=27deg'};
    pdTE.strtype = 'e';
    pdTE.num     = [0 Inf];
    pdTE.val     = {[]};

% ---------------------------------------------------------------------
% t1FA 
% ---------------------------------------------------------------------
    t1FA         = cfg_entry;
    t1FA.tag     = 't1FA';
    t1FA.name    = 'FA values for the T1w volumes';
    t1FA.help    = {'Flip angle for the T1w volumes. Please insert a value for each T1w volume. ' ...
                    'If not specified, the toolbox will try to read it from the description field in the header of the image files. ' ...
                    'This description has to contain a part like TR=22.5ms/TE=2.71ms/FA=27deg'};
    t1FA.strtype = 'e';
    t1FA.num     = [0 Inf];
    t1FA.val     = {[]};
    
% ---------------------------------------------------------------------
% mtFA 
% ---------------------------------------------------------------------
    mtFA         = cfg_entry;
    mtFA.tag     = 'mtFA';
    mtFA.name    = 'FA values for the MTw volumes';
    mtFA.help    = {'Flip angle for the MTw volumes. Please insert a value for each MTw volume. ' ...
                    'If not specified, the toolbox will try to read it from the description field in the header of the image files. ' ...
                    'This description has to contain a part like TR=22.5ms/TE=2.71ms/FA=27deg'};
    mtFA.strtype = 'e';
    mtFA.num     = [0 Inf];
    mtFA.val     = {[]};
% ---------------------------------------------------------------------
% pdFA 
% ---------------------------------------------------------------------
    pdFA         = cfg_entry;
    pdFA.tag     = 'pdFA';
    pdFA.name    = 'FA values for the PDw volumes';
    pdFA.help    = {'Flip angle for the PDw volumes. Please insert a value for each PDw volume. ' ...
                    'If not specified, the toolbox will try to read it from the description field in the header of the image files. ' ...
                    'This description has to contain a part like TR=22.5ms/TE=2.71ms/FA=27deg'};
    pdFA.strtype = 'e';
    pdFA.num     = [0 Inf];
    pdFA.val     = {[]};
        
    
% --------------------------------------------------------------------
% .mat file containing an ESTATIC model
% ---------------------------------------------------------------------
    ESTAmodel         = cfg_files;
    ESTAmodel.tag     = 'ESTAmodel';
    ESTAmodel.name    = 'ESTATICS model';
    ESTAmodel.help    = {'Select the .mat file containing the ESTATICS model created by this toolbox.'};
    ESTAmodel.filter  = '.mat';
    ESTAmodel.ufilter = '.*';
    ESTAmodel.num     = [1 1];
    
% ---------------------------------------------------------------------        
% branch with MT files
% ---------------------------------------------------------------------
     MT_branch         = cfg_exbranch;
     MT_branch.tag     = 'MT_branch';
     MT_branch.name    = 'Model with T1w, PDw and MTw volumes';
     MT_branch.val     = {t1Files mtFiles pdFiles maskFile odir tr2 coregIM saveESTA  ...
                          kstar lambda b1FileA b1FileP height tol  ...
                          t1TR mtTR pdTR t1TE mtTE pdTE t1FA mtFA pdFA};
     MT_branch.help    = {'This branch implements multi parameter mapping (MPM) for a dataset with T1w, PDw and MTw volumes.'};
     MT_branch.prog    = @spm_local_mpm;
     
% ---------------------------------------------------------------------        
% branch without MT files
% ---------------------------------------------------------------------
     withoutMT_branch         = cfg_exbranch;
     withoutMT_branch.tag     = 'withoutMT_branch';
     withoutMT_branch.name    = 'Model without MTw volumes';
     withoutMT_branch.val     = {t1Files pdFiles maskFile odir coregIM saveESTA  ...
                                kstar lambda b1FileA b1FileP height tol ... 
                                t1TR pdTR t1TE pdTE t1FA pdFA}; 
     withoutMT_branch.help    = {'This branch implements multi parameter mapping (MPM) for a dataset without MTw volumes.'};
     withoutMT_branch.prog    = @spm_local_mpm_noMT;

% ---------------------------------------------------------------------        
% branch using an existing ESTATICS model
% ---------------------------------------------------------------------
     ESTAmodel_branch         = cfg_exbranch;
     ESTAmodel_branch.tag     = 'ESTAmodel_branch';
     ESTAmodel_branch.name    = 'Use an existing ESTATICS model';
     ESTAmodel_branch.val     = {ESTAmodel odir kstar lambda b1FileA b1FileP tr2 height};
     ESTAmodel_branch.help    = {'This branch implements the smoothing and final calculation of the quantitative maps R1, R2*, PD (and MT) given an existing ESTATICS model.'};
     ESTAmodel_branch.prog    = @spm_local_mpm_givenESTATICS;
     
% ---------------------------------------------------------------------        
% branch producing .nii file from an existing ESTATICS model
% ---------------------------------------------------------------------
     extractNII_branch         = cfg_exbranch;
     extractNII_branch.tag     = 'extractNII_branch';
     extractNII_branch.name    = 'Get image files from ESTATICS model';
     extractNII_branch.val     = {ESTAmodel odir height};
     extractNII_branch.help    = {'This branch creates the NIfTI files of the estimated ESTATICS parameter ' ...
         'from an existing ESTATICS model. No smoothing is applied. No calculation of R1, R2*, PD (and MT) maps is performed.'};
     extractNII_branch.prog    = @spm_local_mpm_extractNII;

% ---------------------------------------------------------------------
%   MPM toolbox
% ---------------------------------------------------------------------
    MPM         = cfg_choice;
    MPM.tag     = 'MPM';
    MPM.name    = 'MPM Multi-Parameter Mapping';
    MPM.values     = {MT_branch withoutMT_branch ESTAmodel_branch extractNII_branch};
    MPM.help    = {'This toolbox implements multi parameter mapping based on the estimation of the ESTATICS model along with adaptive smoothing of the volumes.'};

%======================================================================
function [] = spm_local_mpm(job)

    
    if ~isdeployed, addpath(fullfile(spm('Dir'),'toolbox','MPM')); end
    job.sdim = [0 0 0];
    
    mpm_assert_output_directory(job.odir{1});
    
    mpmESTATICS(job);
    
end

function [] = spm_local_mpm_noMT(job)

    
    if ~isdeployed, addpath(fullfile(spm('Dir'),'toolbox','MPM')); end
    mpm_assert_output_directory(job.odir{1});
    
    job.mtFiles = [];
    job.mtTR = [];
    job.mtTE = [];
    job.mtFA = [];
    job.tr2 = 1;
    job.sdim = [0 0 0];
    mpmESTATICS(job);
    
end

function [] = spm_local_mpm_givenESTATICS(job)

    
    if ~isdeployed, addpath(fullfile(spm('Dir'),'toolbox','MPM')); end
    
    mpm_assert_output_directory(job.odir{1});
    
    mpm_givenESTATICS(job);
    
end

function [] = spm_local_mpm_extractNII(job)

    
    if ~isdeployed, addpath(fullfile(spm('Dir'),'toolbox','MPM')); end
    
    mpm_assert_output_directory(job.odir{1});
    
    mpm_extractNII(job);
    
end

function [] = mpm_assert_output_directory(odir)
        % checks if it is possible to write to the chosen output directory
    timestamp = datestr(now);
    timestamp = strrep(timestamp,' ','_');
    timestamp = strrep(timestamp,'-','_');
    timestamp = strrep(timestamp,':','_');
    try mkdir(odir,timestamp);
    catch ME
        fprintf(ME.message);
        error('You cannot write to the chosen output directory. Please choose another!');
    end
    rmdir(fullfile(odir,timestamp));
end

end




