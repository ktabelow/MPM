function PDcalculation(PMT, PA)

% Performs "PD flattening"
%
% FORMAT PDcalculation(PMT, PA, Pmask)
%_______________________________________________________________________
% INPUTS
%   PMT     Filename of the estimated MT map.
%   PA      Filename of the estimated A map.
%   Pmask   (optional) Filename of the ???
%==========================================================================

    spm_jobman('initcfg');
    disp('----- Calculating Proton Density map -----');

    if ~exist('hmri_get_defaults', 'var') %TODO: this is the wrong question as hmri_get_defaults is a MATLAB function 
        PDproc.PDmap    = 1;      % Calculation of PD maps requires a B1 map. Set to 0 if a B1 map is not available
        PDproc.WBMaskTh = 0.1;    % Threshold for calculation of whole-brain mask from TPMs
        PDproc.WMMaskTh = 0.95;   % Threshold for calculation of white-matter mask from TPMs
        PDproc.biasreg  = 10^(-5);
        PDproc.biasfwhm = 50;
        threshA         = 10^5;
        threshMT        = 5;
    else
        % get PD processing default settings etc. from hmri toolbox
        PDproc   = hmri_get_defaults('PDproc');
        threshA  = hmri_get_defaults('qMRI_maps_thresh.A');
        threshMT = hmri_get_defaults('qMRI_maps_thresh.MT');
    end
    interpcons  = -4; % interpolation parameter for SPM registration
    

    [pMT, ~, ~] = spm_fileparts(PMT);
    % Creation of whole-brain and white-matter masks
    VG = spm_vol(PMT);
    dm = VG.dim;         % get volume spatial dimension
%    MTforA_vol = ACID_read_vols(VG, VG, interpcons); % TODO: use spm function here
    MTforA_vol    = zeros(dm);
    for p=1:dm(3)
        M = VG.mat*spm_matrix([0 0 p]);
        MTforA_vol(:,:,p) = spm_slice_vol(VG, VG.mat\M, dm(1:2),interpcons);
    end
    Atmp = zeros(dm);
    % loop slicewise through the MT file
    % and mask all voxel within 5 voxel from the border
    % and cut MT value at threshMT
    for p = 1:dm(3)
        MTforA = MTforA_vol( : , : , p);
        if (p < 5) || (p > dm(3) - 5)
            MTforA = zeros(size(MTforA, 1), size(MTforA, 2));
        else
            MTforA(1:5, : ) = 0; 
            MTforA(end-5:end, : ) = 0;
            MTforA( : , 1:5) = 0; 
            MTforA( : , end-5:end) = 0;
        end
        Atmp( : , : , p) = max(min(MTforA, threshMT), -threshMT);
    end
    % write to disk as <filename of the MT file> plus suffix '_MTforA'
    VG.fname = [VG.fname(1:end-4) '_MTforA.nii'];
    VMTforA = spm_write_vol(VG, Atmp);
    
    clear matlabbatch
    matlabbatch{1}.spm.spatial.preproc.channel.vols = {VMTforA.fname};
    matlabbatch{1}.spm.spatial.preproc.channel.write = [1 0]; % saving bias field
    spm_jobman('run', matlabbatch);

    [pMTforA, fMTforA, eMTforA] = spm_fileparts(VMTforA.fname);
    TPMs = spm_read_vols(spm_vol(spm_select('FPList', spm_fileparts(VMTforA.fname), ['^c[1-6]' fMTforA eMTforA])));
    % Whole brain mask = WM + GM + bone but not the CSF! (CP: why?)
    WBmask = zeros(dm);
    WBmask(sum(cat(4, TPMs( : , : , : , 1:2), TPMs( : , : , : , end)), 4) >= PDproc.WBMaskTh) = 1;
    % White matter mask as well
    WMmask = zeros(dm);
    WMmask(squeeze(TPMs( : , : , : , 2)) >= PDproc.WMMaskTh) = 1;

    % remove all temporary files
    temp = spm_select('FPList', pMT, '^.*_MTforA');
    for counter = 1:size(temp, 1)
        delete(deblank(temp(counter, : )));
    end
    % End of creation of whole-brain and white-matter masks

    % Saves masked A map for bias-field correction later
    [pA, fA, eA] = spm_fileparts(PA);
    Vsave = spm_vol(PA);
    Amap = spm_read_vols(Vsave) .* WBmask;
    Amap(Amap == Inf) = 0;
    Amap(isnan(Amap)) = 0;
    Amap(Amap >= threshA) = 0; 
    Vsave.fname = fullfile(spm_str_manip(Vsave.fname, 'h'), ['masked_' spm_str_manip(Vsave.fname, 't')]);
    spm_write_vol(Vsave, Amap);
    Pmask = spm_select('FPList', spm_fileparts(Vsave.fname), ['^' spm_str_manip(Vsave.fname, 't')]);

    % Bias-field correction of masked A map
    clear matlabbatch
    matlabbatch{1}.spm.spatial.preproc.channel.vols = {Pmask};
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg = PDproc.biasreg;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = PDproc.biasfwhm;
    matlabbatch{1}.spm.spatial.preproc.channel.write = [1 0]; % saving bias field
    % spm_jobman('initcfg');
    spm_jobman('run', matlabbatch);

    temp = spm_select('FPList', pA, ['^(c|masked).*' fA]);
    for counter = 1:size(temp, 1)
        delete(deblank(temp(counter, : )));
    end

    % Bias field correction of A map. The bias field is calculated on
    % the masked A map but we want to apply it on the unmasked A map. We
    % therefore need to explicitly load the bias field and apply it on the original A map instead of just
    % loading the bias-field corrected A map from the previous step
    % PA=spm_select('FPList',ptmp ,'^s.*_A.(img|nii)$');
    bf = fullfile(pA, spm_select('List', pA, ['^BiasField_masked_' fA eA]));
    BF = double(spm_read_vols(spm_vol(bf)));
    Y = BF.*spm_read_vols(spm_vol(PA));

    % Calibration of flattened A map to % water content using typical white
    % matter value from the litterature (69%)

    A_WM = WMmask.*Y;
    Y = Y/mean(A_WM(A_WM~=0))*69;
    sprintf('mean White Matter intensity: %04d', mean(A_WM(A_WM~=0)))
    sprintf('SD White Matter intensity %04d', std(A_WM(A_WM~=0), [], 1))
    Y(Y>200 ) = 0;
    % MFC: Estimating Error for data set to catch bias field issues:
    errorEstimate = std(A_WM(A_WM > 0))./mean(A_WM(A_WM > 0));
    Vsave = spm_vol(PA);
    if (strfind(fA, 'A_'))
        Vsave.fname = fullfile(pA, [regexprep(fA,'^A_','PD_') eA]);
    elseif (strfind(fA, 'Alinear_'))
        Vsave.fname = fullfile(pA, [regexprep(fA,'^Alinear_','PDlinear_') eA]);
    else 
        Vsave.fname = fullfile(pA, ['PD_' fA eA]);
    end
    Vsave.descrip = ['PD Map.  Error Estimate: ', num2str(errorEstimate)];
    if errorEstimate > 0.06
        % MFC: Testing on 15 subjects showed 6% is a good cut-off:
        warning(['Error estimate is high: ', Vsave.fname]);
    end

    % V.fname = P;
    spm_write_vol(Vsave, Y);

    temp = spm_select('FPList', pA, ['^Bias.*' fA]);
    for counter = 1:size(temp, 1)
        delete(deblank(temp(counter, : )));
    end

end

