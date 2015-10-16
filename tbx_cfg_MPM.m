function MPM = tbx_cfg_MPM
% Configuration file for toolbox MPM
%_______________________________________________________________________
% Copyright (C) 2015 WIAS Berlin

% Chiara from K.Tabelow configuration file for aws4SMP
%
% ---------------------------------------------------------------------
% spm_file t1Files.img
% ---------------------------------------------------------------------
    spm_t1Files         = cfg_files;
    spm_t1Files.tag     = 't1_file';
    spm_t1Files.name    = 'Select t1Files';
    spm_t1Files.help    = {'Select the t1 files.'};
    spm_t1Files.filter  = 'image';
    spm_t1Files.ufilter = '.*';
    spm_t1Files.num     = [1 8];
%    --------------------------------------------------------------------
% spm_file SPM.mat
% ---------------------------------------------------------------------
    spm_mtFiles         = cfg_files;
    spm_mtFiles.tag     = 'mt_file';
    spm_mtFiles.name    = 'Select mtFiles';
    spm_mtFiles.help    = {'Select the mt files.'};
    spm_mtFiles.filter  = 'image';
    spm_mtFiles.ufilter = '.*';
    spm_mtFiles.num     = [1 8];
% ---------------------------------------------------------------------
% hmax Maximal Bandwidth for test
% ---------------------------------------------------------------------
    sdim        = cfg_entry;
    sdim.tag     = 'sdim';
    sdim.name    = 'sdim';
    sdim.help    = {'Dimension of the cubus. '};
    sdim.strtype = 'e';
    sdim.num     = [1 3];
    sdim.val    = {[322 368 256]};

% ---------------------------------------------------------------------
% ladjust for test
% ---------------------------------------------------------------------
    ladjust        = cfg_entry;
    ladjust.tag     = 'ladjust';
    ladjust.name    = 'ladjust';
    ladjust.help    = {'Adjust adaptation bandwidth. '};
    ladjust.strtype = 'e';
    ladjust.num     = [1 1];
    ladjust.val    = {1};

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
    MPM.val     = {spm_t1Files spm_mtFiles sdim};
    MPM.help    = {'This toolbox implements multi parameter mapping for SPM.'}';
    MPM.prog    = @spm_local_mpm;

%======================================================================
function spm_local_mpm(job)

    % access by job.spm_t1file, job.hmax, job.ladjust
    if ~isdeployed, addpath(fullfile(spm('Dir'),'toolbox','MPM')); end

    %load(job.spm_file{:});
    fileattrib(job.spm_t1Files{:})
    
   % [I, xCon] = spm_conman(SPM, 'T&F', Inf, ...
   %         '    Select contrasts...', ' for conjunction', 1);
    % Now I contains the selected contrast in the struct xCon

    %if isempty(fieldnames(SPM.xCon))
    %  SPM.xCon = xCon;
    %end
       %createDataSet(sdim,job.spm_t1Files,job.spm_mtFiles)
    %spm_smooth_gamma(SPM, I, job.ladjust, job.hmax);
end





end




