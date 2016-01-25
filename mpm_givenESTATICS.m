% wrapper function for MPM functions given an ESTATICS model
%
%
% =========================================================================
% 2016/01/08
%
% 
% function [] = mpm_givenESTATICS(job)
% 
% call all the functions of the MPM iterating on each region containing "heigth"
% planes x,y. 
%
% Input:
%  job             - a struct containing:
%  
%  ESTAmodel       - a file containing the metadata
%  tr2             -
%  height          - height of the region of interest 
%  kstar           - number of steps for the smoothing algorithm (if 0 no smoothing)          
%  lambda          - adaptation bandwidth
%  b1File          - correction field
%
%  The ESTAmodel must contain:
%
%  maskFile        - the complete name of the file containing the mask
%  mask            - the mask
%  invCov          - inverse covariance matrix
%  t1Files 
%  mtFiles         - original input files
%  pdFiles 
%  nFiles          - number of original imput files
%  sdim            - spatial dimensionality of the data 
%  nv              - number of parameters (3 or 4)
%  modelCoeff      - a 3 or 4 dimensional array containing the vauels of
%                    the parameters
%  TR
%  TE              - arrays containing the relaxation rate, flip angle and echo time  
%  FA
%
% no Output
% the 3 or 4 final variables R1, PD, R2star and in case of a model with MT also delta
% are written each in a .nii file
%                 
%  
%
% =========================================================================

function [] = mpm_givenESTATICS(job)

    fprintf('Starting the MPM at %s \n',datestr(now));


 %% works on the all volume, level by level
    
    fprintf('Loading data...\n');
    meta=matfile(char(job.ESTAmodel));
%     coeffFiles = meta.modelCoeff; %
    modelMPM.nv = meta.nv;
    modelMPM.nFiles = meta.nFiles;
    modelMPM.t1Files = meta.t1Files;
    modelMPM.pdFiles = meta.pdFiles;
    if modelMPM.nv==4,
         modelMPM.mtFiles = meta.mtFiles;
    end
    modelMPM.FA = meta.FA;
    modelMPM.TE = meta.TE;
    modelMPM.TR = meta.TR;
    modelMPM.TEScale = meta.TEScale;
    modelMPM.DataScale = meta.DataScale;
    
    Vt1 = spm_vol(modelMPM.t1Files{1}); 
    Vpd = spm_vol(modelMPM.pdFiles{1});  
    if modelMPM.nv==4, Vmt = spm_vol(modelMPM.mtFiles{1}); end 
   
    sdim = meta.sdim;
    modelMPM.sdim = sdim; 
   
 
    % calculates how big the overlapping has to be
    % to assure a good smoothing    
    hmax = 1.25^(job.kstar/3);
    hakt = gethani (1, 1.25*hmax, 2, 1.25^job.kstar, [1 1], 1e-4);
    hdelta = ceil(hakt); % height of half of the overlapping
    
    % check for the input height and adjust it if too big or too small
    if abs(job.height)>sdim(3),
        job.height = sdim(3);
    elseif job.height<0,
        job.height = abs(job.height);
    end
    
    % calculates the interval between the zStart of the different levels
    if job.height>2*hdelta,
        interval = int64(job.height - 2*hdelta);
    else
        job.height = (job.height + 2*hdelta +2);
        interval = int64(job.height - 2*hdelta);
    end
    
    
    
    % preparing result variables
    R1 = zeros(sdim);
    R2star = zeros(sdim);
    PD = zeros(sdim);
    delta = zeros(sdim);

    % prepares variable to save the whole mask (it could be avoided) 
    totalmask = zeros(sdim);
    
    % looking for weights
    try 
        wghts=getWeights(meta.t1Files{1});
    catch 
    wghts = [];
    end
    
    % prepares the data        
    %    dataset = createDataSet(job.sdim,job.t1Files,job.pdFiles,job.mtFiles,char(job.maskFile),job.t1TR,job.pdTR,job.mtTR,job.t1TE,job.pdTE,job.mtTE,job.t1FA,job.pdFA, job.mtFA);
        
    
   
    
    spm_progress_bar('Init',sdim(3),'planes completed');
    
    %% start iteration on all the levels
    for startLayerVoxel = 1:interval:sdim(3),
        
        spm_progress_bar('Set',startLayerVoxel);
        
        zStart = double(startLayerVoxel);        
%         fprintf('on level %d \n', zStart);
        if sdim(3)-(startLayerVoxel + job.height)> 2*hdelta,
            % in case the next starting point has enough planes after it
            zEnd = double(startLayerVoxel + job.height); 
        else
            % in case the next starting point has NOT enough planes after it
            zEnd = sdim(3);
        end
 
        modelMPM.zStart = zStart;
        modelMPM.zEnd = zEnd;
        
        %% read the mask, model coefficients and invCov between zStart and zEnd
%         if isempty(meta.maskFile) || (length(meta.maskFile)==1 && strcmp(meta.maskFile{1},''))
%             meta.mask=ones([sdim(1) sdim(2) zEnd-zStart+1]);
%         else 
            slices=zStart:zEnd; %1:sdim(3);
            %[mask(:,:,:),~] = loadImageSPM(fullfile(dir,[maskFile{1} '.nii']),'slices',slices);
            if isfield(modelMPM, 'mask'),
                modelMPM= rmfield(modelMPM,'mask');
            end
            mask=meta.mask(:,:,slices);
            modelMPM.mask = mask;
            clear mask;
%             [mask(:,:,:),~] = loadImageSPM(fullfile(meta.maskFile) ,'slices',slices);
%             mask = round(mask(:));
%             mask = reshape (mask, [sdim(1) sdim(2) zEnd-zStart+1]);
%             modelMPM.mask=mask;
%             clear mask;
            if isfield(modelMPM, 'modelCoeff'),
                modelMPM= rmfield(modelMPM,'modelCoeff');
            end                  
            modelCoeff(1,:,:,:) = meta.modelCoeff(1,:,:,slices);
            modelCoeff(2,:,:,:) = meta.modelCoeff(2,:,:,slices);
            modelCoeff(3,:,:,:) = meta.modelCoeff(3,:,:,slices);
%             [modelMPM.modelCoeff(1,:,:,:),~]= loadImageSPM(fullfile(coeffFiles{1}) ,'slices',slices);
%             [modelMPM.modelCoeff(2,:,:,:),~]= loadImageSPM(fullfile(coeffFiles{2}) ,'slices',slices);
%             [modelMPM.modelCoeff(3,:,:,:),~]= loadImageSPM(fullfile(coeffFiles{3}) ,'slices',slices);
            if meta.nv==4,  
                modelCoeff(4,:,:,:) = meta.modelCoeff(4,:,:,slices);
%                 [modelMPM.modelCoeff(4,:,:,:),~]= loadImageSPM(fullfile(coeffFiles{4}) ,'slices',slices);            
            end
            modelMPM.modelCoeff = modelCoeff;
            clear modelCoeff;
            
            if isfield(modelMPM, 'invCov'),
                modelMPM= rmfield(modelMPM,'invCov');
            end    
            invCov=meta.invCov(:,:,:,:,slices);            
            modelMPM.invCov=invCov;
            clear invCov;
            
%         end
        
        
        
        
        
        %% smooth and calculate the final four parameters with the functions
        % function [modelS] = smoothESTATICSmask(model, varargin)
        % function [qi] = calculateQI(model, varargin)
        % (if kstar is 0, the smoothing step isn't done)
        if job.kstar~=0, 
            modelMPMs_mask = smoothESTATICSmask(modelMPM, 'verbose', false, 'wghts', wghts, 'lambda',job.lambda); 
            qiSnew = calculateQI(modelMPMs_mask, 'TR2',job.tr2,'b1File',job.b1File , 'verbose', false);
        else
            qiSnew = calculateQI(modelMPM, 'TR2',job.tr2,'b1File',job.b1File , 'verbose', false);
        end
        
        
        %% save in the result variables
        if zStart==1
            R1(:,:,zStart:zEnd) = qiSnew.R1;
            R2star(:,:,zStart:zEnd) = qiSnew.R2star;
            PD(:,:,zStart:zEnd) = qiSnew.PD;
            if modelMPM.nv==4, delta(:,:,zStart:zEnd) = qiSnew.delta; end%
            totalmask(:,:,zStart:zEnd) = qiSnew.model.mask;
        else 
            R1(:,:,zStart+hdelta:zEnd) = qiSnew.R1(:,:, 1+hdelta: (zEnd-zStart+1) );
            R2star(:,:,zStart+hdelta:zEnd) = qiSnew.R2star(:,:, 1+hdelta: (zEnd-zStart+1) );
            PD(:,:,zStart+hdelta:zEnd) = qiSnew.PD(:,:, 1+hdelta: (zEnd-zStart+1) );
            if modelMPM.nv==4, delta(:,:,zStart+hdelta:zEnd) = qiSnew.delta(:,:, 1+hdelta: (zEnd-zStart+1) ); end %
            totalmask(:,:,zStart+hdelta:zEnd) = qiSnew.model.mask(:,:, 1+hdelta: (zEnd-zStart+1) );
            
        end
        if zEnd==sdim(3),
                break;
        end
    
    end
    spm_progress_bar('Set',sdim(3));
    spm_progress_bar('Clear');
    
    %% steps to save nii files
    % set all the NaN values to 0 (that makes easier to display correctly)
    R1(isnan(R1))=0;
    R2star(isnan(R2star))=0;
    PD(isnan(PD))=0;
    if modelMPM.nv==4, delta(isnan(delta))=0; end
    
    % set all values of the voxels outside the mask to 0
    R1(totalmask<1)=0;
    R2star(totalmask<1)=0;
    PD(totalmask<1)=0;
    if modelMPM.nv==4, delta(totalmask<1)=0; end
   
    try
    % write 3 or 4 files for R1, PD, R2star and in case delta
    % function []= write_small_to_file_nii(outputdir,filenamepr, big_volume,small_volume_data,zStart, zEnd, sdim)
    write_small_to_file_nii(pwd,'R1_', Vt1, R1, 1, sdim(3), sdim);
    write_small_to_file_nii(pwd,'R2star_', Vt1, R2star, 1, sdim(3), sdim);
    write_small_to_file_nii(pwd,'PD_', Vpd, PD, 1, sdim(3), sdim);
    if modelMPM.nv==4,  write_small_to_file_nii(pwd,'delta_', Vmt, delta, 1, sdim(3), sdim); end
    catch
        error('it was not possible to save the resulting .nii files. Check to have writing rights in the current directory')
    end
    fprintf('Ending the MPM at %s \n',datestr(now));

end