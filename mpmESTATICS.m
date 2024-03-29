% wrapper function for all the steps of the MPM toolbox
%
% =========================================================================
% 2015/12/08
% 
% by C. D'Alonzo
%
% FORTRAN code written by J. Polzehl
% 
% function [] = mpmESTATICS(job)
% 
% call all the functions of the MPM iterating on each region containing "heigth"
% planes x,y. 
%
% Input:
%  job             - a struct containing:
%
%  t1Files 
%  mtFiles         - input Files
%  pdFiles 
%  sdim            - spatial dimensionality of the data 
%                    (in default case [0 0 0] will be read from the file)
%  saveESTA        - 1 to save a.mat file containing the ESTATICS model or 0 otherwise 
%  coregIM         - 1 to coregister the image to the first T1w image, 0
%                    otherwise
%  confInt         - 1 to estimate the confidence intervals for R1 and PD,
%                    0 otherwise
%  tr2 
%  height          - height of the region of interest 
%  kstar           - number of steps for the smoothing algorithm (if 0 no smoothing)          
%  lambda          - adaptation bandwidth
%  maskFile        - mask containing the voxel to consider
%  b1FileA         - amplitude image (correction field)
%  b1FileB         - phase image (correction field)
%  odir            - output directory
%  t1TR mtTR pdTR  
%  t1TE mtTE pdTE  - parameter (read from file if empty)
%  t1FA mtFA pdFA
%
% no Output
%
% the 3 or 4 final variables R1, PD, R2star and in case of a model with MT also delta
% are written each in a .nii file
%                 
% In case the users chose to save the ESTATICS model, a .mat file with
% references and data regarding the model will be also saved.
%
% In case the user chose to estimate the confidence interval, also
% 4 .nii files each containing upper and lower boundaries for R1 and PD are written
% and a reference to them is saved in the ESTATICS model .mat file, if that option is active.
%
% =========================================================================

function [] = mpmESTATICS(job)

    fprintf('Starting the MPM at %s \n',datestr(now));


 %% works on the all volume, level by level
    
    V = spm_vol(job.t1Files{1});
    
    % read the dimension array if not given
    if sum(job.sdim)==0      
      job.sdim = V.dim;
    end  
 
    % calculates how big the overlapping has to be
    % to assure a good smoothing    
    hmax = 1.25^(job.kstar/3);
    hakt = gethani (1, 1.25*hmax, 2, 1.25^job.kstar, [1 1], 1e-4);
    hsigma = 2;
    hdelta = max(hsigma, ceil(hakt)); % height of half of the overlapping
    
    % check for the input height and adjust it if too big or too small
    if abs(job.height)>job.sdim(3)
        job.height = job.sdim(3);
    elseif job.height<0
        job.height = abs(job.height);
    end
    
    % calculates the interval between the zStart of the different levels
    if job.height>2*hdelta
        interval = int64(job.height - 2*hdelta);
    else
        job.height = (job.height + 2*hdelta +2);
        interval = int64(job.height - 2*hdelta);
    end
    
    
    
    % preparing result variables
    R1 = zeros(job.sdim);
    R2star = zeros(job.sdim);
    PD = zeros(job.sdim);
    delta = zeros(job.sdim);

    
    if job.confInt==1
     R1_low = zeros(job.sdim);
     R1_up  = zeros(job.sdim);
     PD_low = zeros(job.sdim);
     PD_up  = zeros(job.sdim);
    end
    
    % prepares variable to save the whole mask (it could be avoided) 
    totalmask = zeros(job.sdim);
    
    % looking for weights
    try 
        wghts=getWeights(job.t1Files{1});
    catch 
    wghts = [];
    end
    
    % prepares the data        
        dataset = createDataSet(job.sdim,job.t1Files,job.pdFiles,job.mtFiles,char(job.maskFile),job.t1TR,job.pdTR,job.mtTR,job.t1TE,job.pdTE,job.mtTE,job.t1FA,job.pdFA, job.mtFA);
    
    % calculates the PD average to which all images are coregistered
        pdVol = spm_vol(job.pdFiles{1});
        Pdmean = zeros(job.sdim);
        for i=1:length(job.pdFiles) 
            [I,~,~,~] = loadImageSPM(fullfile(job.pdFiles{i}));
            Pdmean = Pdmean + I;
        end
        Pdmean = Pdmean./length(job.pdFiles);
        write_small_to_file_nii(job.odir{1},'PD_average_', pdVol, Pdmean,1, job.sdim(3), job.sdim);
        [~, nam, ~] = spm_fileparts(pdVol.fname);
        job.Pdmean = {fullfile(job.odir{1},strcat('PD_average_',nam,'.nii'))};
        dataset.Pdmean = job.Pdmean;
        clear Pdmean;
        
    % if amplitude and phase image are both present, 
    % produces the b1 correction field with the same dimensionality of the data 
    
    if (isempty(job.b1FileA) || (length(job.b1FileA)==1 && strcmp(job.b1FileA{1},'')) ) || ...
            (isempty(job.b1FileP) || (length(job.b1FileP)==1 && strcmp(job.b1FileP{1},'')) )
        job.b1File = [];
    else
        try
            fprintf('\nProducing coregistered correction field... \n');
            
        b1_aTMat = MPM_get_coreg_matrix(spm_vol(job.b1FileA{1}),spm_vol(job.Pdmean{1}));
        b1Vol = MPM_read_coregistered_vol(spm_vol(job.b1FileP{1}),spm_vol(job.Pdmean{1}),'affineTransMatrix',b1_aTMat);
        fprintf('\nSaving the coregistered correction field file... ');
        write_small_to_file_nii(job.odir{1},'b1File_registeredTo_', pdVol,b1Vol,1, job.sdim(3), job.sdim);
        [~, nam, ~] = spm_fileparts(pdVol.fname);
        job.b1File = {fullfile(job.odir{1},strcat('b1File_registeredTo_',nam,'.nii'))};
        fprintf('Done. \n');
        catch ME
           fprintf(ME.message);
            fprintf('\nSomething went wrong when registering and saving the correction field file. No correction field will be applied. \n');
             job.b1File = [];
        end
    end
    % coregister the images
    % NB the files are not modified (not even the headers), but the saved
    % matrices are used later to read the datas
    if job.coregIM==1
        % not coregistering the mask. The user is responsible for creating a correct mask
        % fprintf('\nCoregistering mask... \n'); 
        % dataset.mask_aTMat = MPM_get_coreg_matrix(spm_vol(fullfile(dataset.maskFile)),spm_vol(job.pdFiles{1}));
        
        %dataset.t1_aTMat = zeros(4,4, length(job.t1Files)); % old
        %registration system
        %dataset.pd_aTMat = zeros(4,4, length(job.pdFiles));
        
        dataset.t1_aTMat = zeros(4,4);
        dataset.pd_aTMat = zeros(4,4);
        
        % calculates the T1 average to calculate the transformation matrix        
        t1mean = zeros(job.sdim);
        for i=1:length(job.t1Files) 
            [I,~,~,~] = loadImageSPM(fullfile(job.t1Files{i}));
            t1mean = t1mean + I;
        end
        t1mean = t1mean./length(job.t1Files);
        write_small_to_file_nii(job.odir{1},'T1_average_', pdVol, t1mean,1, job.sdim(3), job.sdim);
        [~, nam, ~] = spm_fileparts(pdVol.fname);
        job.t1mean = {fullfile(job.odir{1},strcat('T1_average_',nam,'.nii'))};
        dataset.t1mean = job.t1mean;
        clear t1mean;
        
        fprintf('\nCoregistering T1w files... \n');
%         for k=1:length(job.t1Files),
%             dataset.t1_aTMat(:,:,k) = MPM_get_coreg_matrix(spm_vol(fullfile(job.t1Files{k})),spm_vol(job.Pdmean{1}));
%         end
        dataset.t1_aTMat(:,:) = MPM_get_coreg_matrix(spm_vol(fullfile(job.t1mean{1})),spm_vol(job.Pdmean{1}));
         

        fprintf('\nCoregistering PDw files... \n');
%         for k=1:length(job.pdFiles),
%             dataset.pd_aTMat(:,:,k) = MPM_get_coreg_matrix(spm_vol(fullfile(job.pdFiles{k})),spm_vol(job.Pdmean{1}));
%         end  
        dataset.pd_aTMat(:,:) = MPM_get_coreg_matrix(spm_vol(fullfile(job.Pdmean{1})),spm_vol(job.Pdmean{1})); %identity?!?


        if ~isempty(job.mtFiles)
            %dataset.mt_aTMat = zeros(4,4, length(job.mtFiles));
            dataset.mt_aTMat = zeros(4,4);
            fprintf('\nCoregistering MTw files... \n');
%             for k=1:length(job.mtFiles),
%                 dataset.mt_aTMat(:,:,k) = MPM_get_coreg_matrix(spm_vol(fullfile(job.mtFiles{k})),spm_vol(job.Pdmean{1}));
%             end
            
            % calculates the MT average to calculate the transformation matrix        
            MTmean = zeros(job.sdim);
            for i=1:length(job.mtFiles) 
                [I,~,~,~] = loadImageSPM(fullfile(job.mtFiles{i}));
                MTmean = MTmean + I;
            end
            MTmean = MTmean./length(job.mtFiles);
            write_small_to_file_nii(job.odir{1},'MT_average_', pdVol, MTmean,1, job.sdim(3), job.sdim);
            [~, nam, ~] = spm_fileparts(pdVol.fname);
            job.MTmean = {fullfile(job.odir{1},strcat('MT_average_',nam,'.nii'))};
            dataset.MTmean = job.MTmean;
            clear MTmean;
            
            % calculate transfo matrix
            dataset.mt_aTMat(:,:) = MPM_get_coreg_matrix(spm_vol(fullfile(job.MTmean{1})),spm_vol(job.Pdmean{1}));
           
 
        end
    end

    
    
    
    
    % if saving ESTATICS model, prepare the .mat with the metadata
    % 
    if job.saveESTA==1
        % prepare to save the whole invCov directly in the .mat file
        sigma2i=zeros([job.sdim]);
        try
            [~,fname,~] = spm_fileparts(V.fname);
            filename = strcat(job.odir{1},filesep,'ESTATICS_model_',fname,'.mat');
            save(filename,'sigma2i','-v7.3');
            meta = matfile(filename,'Writable',true);
        
        catch ME
            fprintf(ME.message);
            fprintf('\nIt is not possible to save the ESTATICS model. Check to have writing rights in the current directory, enough memory and to be using at least version 7.3.');
            job.saveESTA = 0;
      
        end
        clear sigma2i;
    end
    
    
    % save metadata for the ESTATICS model
    if job.saveESTA==1
        
           try
           meta.sdim = dataset.sdim;
           meta.t1Files = dataset.t1Files;
           meta.pdFiles = dataset.pdFiles;
           if dataset.nv == 4
               meta.mtFiles = dataset.mtFiles;
           else
               meta.mtFiles = [];
           end
           meta.maskFile = dataset.maskFile;
           meta.TR = dataset.TR;
           meta.TE = dataset.TE;
           meta.FA = dataset.FA;
           meta.nv = dataset.nv;
           meta.nFiles = dataset.nFiles; 
           meta.P2_a = dataset.P2_a; 
           meta.P2_b = dataset.P2_b; 
           meta.modelCoeff = zeros([dataset.nv dataset.sdim]);
           catch ME
            fprintf(ME.message);   
            fprintf('\nThere was a problem saving the ESTATICS model. Check to have writing rights in the current directory and to be using at least version 7.3.\n');
            job.saveESTA = 0;      
           end
        
    end
    
    
    spm_progress_bar('Init',job.sdim(3),'planes completed');
    
    %% start iteration on all the levels
    for startLayerVoxel = 1:interval:job.sdim(3)
        
        spm_progress_bar('Set',startLayerVoxel);
        
        zStart = double(startLayerVoxel);        

        if job.sdim(3)-(startLayerVoxel + job.height)> 2*hdelta
            % in case the next starting point has enough planes after it
            zEnd = double(startLayerVoxel + job.height); 
        else
            % in case the next starting point has NOT enough planes after it
            zEnd = job.sdim(3);
        end
 
        dataset.zStart = zStart;
        dataset.zEnd = zEnd;
        
        %% read the mask between zStart and zEnd
        if isempty(dataset.maskFile) || (length(dataset.maskFile)==1 && strcmp(dataset.maskFile{1},''))
            dataset.mask=ones([dataset.sdim(1) dataset.sdim(2) zEnd-zStart+1]);
            if zStart==1, fprintf('no mask file - working on all voxels\n'); end
        else 
            if zStart==1
                V1 = spm_vol(dataset.maskFile);
                maskdim=V1.dim;
            end
            if dataset.sdim~=maskdim
                dataset.mask=ones([dataset.sdim(1) dataset.sdim(2) zEnd-zStart+1]);
                if zStart==1, fprintf('no correct dimension in mask file - working on all voxels\n'); end
            else
                if zStart==1 
                    fprintf('correct dimension in mask file\n'); 
                end
                slices=zStart:zEnd; 
                %[mask(:,:,:),~] = loadImageSPM(fullfile(dataset.maskFile) ,'slices',slices);
                % if we decided to coregister the mask in line 125
                %if isfield(dataset,'mask_aTMat'),
                %    mask(:,:,:) = MPM_read_coregistered_vol(spm_vol(fullfile(dataset.maskFile)),spm_vol(job.pdFiles{1}),'slices',slices,'affineTransMatrix',dataset.mask_aTMat);
                %else
                    mask(:,:,:) = MPM_read_coregistered_vol(spm_vol(fullfile(dataset.maskFile)),spm_vol(job.pdFiles{1}),'slices',slices);
                %end
                mask = round(mask(:));
                mask = reshape (mask, [dataset.sdim(1) dataset.sdim(2) zEnd-zStart+1]);
                dataset.mask=mask;
                clear mask;
            end
        end
        
        
        %% estimate from all the data the 4or3 parameter 
        % function [model] = estimateESTATICS(dataset, varargin)
        warning('off','all');
        modelMPM = estimateESTATICS(dataset,'verbose', false, 'tolerance', job.tol);
        warning('on','all');
        
        
        %% calculate the confidence interval
        if job.confInt==1
            [R1_1,R1_2,PD_1,PD_2] = getConfidenceIntervall(modelMPM,'b1File',job.b1File);
        end
        
        %% save the 4or3 parameters of the ESTATICS model and DataScale and TEScale in the .mat
        if job.saveESTA==1        
         try
           meta.TEScale = modelMPM.TEScale;
           meta.DataScale = modelMPM.DataScale;           
          if zStart==1
              meta.modelCoeff(1,:,:,1:zEnd) = modelMPM.modelCoeff(1,:,:,:);
              meta.modelCoeff(2,:,:,1:zEnd) = modelMPM.modelCoeff(2,:,:,:);
              if dataset.nv==4
                  meta.modelCoeff(3,:,:,1:zEnd) = modelMPM.modelCoeff(3,:,:,:);
                  meta.modelCoeff(4,:,:,1:zEnd) = modelMPM.modelCoeff(4,:,:,:);
              else
                   meta.modelCoeff(3,:,:,1:zEnd) = modelMPM.modelCoeff(3,:,:,:);
              end
            
          else 
              meta.modelCoeff(1,:,:,zStart+hdelta:zEnd) = modelMPM.modelCoeff(1,:,:, 1+hdelta: (zEnd-zStart+1)); 
              meta.modelCoeff(2,:,:,zStart+hdelta:zEnd) = modelMPM.modelCoeff(2,:,:, 1+hdelta: (zEnd-zStart+1)); 
              if dataset.nv==4
                  meta.modelCoeff(3,:,:,zStart+hdelta:zEnd) = modelMPM.modelCoeff(3,:,:, 1+hdelta: (zEnd-zStart+1)); 
                  meta.modelCoeff(4,:,:,zStart+hdelta:zEnd) = modelMPM.modelCoeff(4,:,:, 1+hdelta: (zEnd-zStart+1)) ; 
              else
                  meta.modelCoeff(3,:,:,zStart+hdelta:zEnd) = modelMPM.modelCoeff(3,:,:, 1+hdelta: (zEnd-zStart+1)); 
              end
          end
         catch ME
            fprintf(ME.message);
            fprintf('\nThere was a problem saving the ESTATICS model. Check to have enough free space and to be using at least version 7.3.\n');
            job.saveESTA = 0;      
         end        
        end
        
        
        
        %% smooth and calculate the final four parameters with the functions
        % function [modelS] = smoothESTATICSmask(model, varargin)
        % function [qi] = calculateQI(model, varargin)
        % (if kstar or lambda are 0, the smoothing step isn't done)
        if job.kstar~=0 && job.lambda~=0 
            modelMPMs_mask = smoothESTATICSmask(modelMPM, 'verbose', false, 'wghts', wghts, 'kstar',job.kstar, 'patchsize', job.patchsize,'lambda',job.lambda); 
            qiSnew = calculateQI(modelMPMs_mask, 'TR2',job.tr2,'b1File',job.b1File , 'verbose', false);
            qiSnewLinear = calculateQIlinear(modelMPMs_mask, 'TR2',job.tr2,'b1File',job.b1File , 'verbose', false);
        else
            qiSnew = calculateQI(modelMPM, 'TR2',job.tr2,'b1File',job.b1File , 'verbose', false);            
            qiSnewLinear = calculateQIlinear(modelMPM, 'TR2',job.tr2,'b1File',job.b1File , 'verbose', false);            
        end

        
        %% save in the result variables
        if zStart==1
            R1(:,:,zStart:zEnd) = qiSnew.R1;
            R1linear(:,:,zStart:zEnd) = qiSnewLinear.R1;
            R2star(:,:,zStart:zEnd) = qiSnew.R2star;
            PD(:,:,zStart:zEnd) = qiSnew.PD;
            PDlinear(:,:,zStart:zEnd) = qiSnewLinear.PD;
            if ~isempty(job.mtFiles) || (length(job.mtFiles)==1 &&~strcmp(job.mtFiles{1},'')), delta(:,:,zStart:zEnd) = qiSnew.delta; end%
            if ~isempty(job.mtFiles) || (length(job.mtFiles)==1 &&~strcmp(job.mtFiles{1},'')), deltaLinear(:,:,zStart:zEnd) = qiSnewLinear.delta; end%
            totalmask(:,:,zStart:zEnd) = qiSnew.model.mask;
            
            if job.confInt==1
                R1_low(:,:,zStart:zEnd) = R1_1;
                R1_up(:,:,zStart:zEnd) = R1_2;
                PD_low(:,:,zStart:zEnd) = PD_2;
                PD_up(:,:,zStart:zEnd) = PD_1;
            end
            
            if job.saveESTA==1        
               try
%               meta.invCov(:,:,:,:,zStart:zEnd)=modelMPM.invCov;
               meta.sigma2i(:,:,zStart:zEnd)=modelMPM.sigma2i;
               catch ME
                   fprintf(ME.message);                   
                   fprintf('\nThere was a problem saving the ESTATICS model. Check to have enough free space and to be using at least version 7.3.\n');
                   job.saveESTA = 0;      
               end        
            end          
            
        else 
            R1(:,:,zStart+hdelta:zEnd) = qiSnew.R1(:,:, 1+hdelta: (zEnd-zStart+1) );
            R1linear(:,:,zStart+hdelta:zEnd) = qiSnewLinear.R1(:,:, 1+hdelta: (zEnd-zStart+1) );
            R2star(:,:,zStart+hdelta:zEnd) = qiSnew.R2star(:,:, 1+hdelta: (zEnd-zStart+1) );
            PD(:,:,zStart+hdelta:zEnd) = qiSnew.PD(:,:, 1+hdelta: (zEnd-zStart+1) );
            PDlinear(:,:,zStart+hdelta:zEnd) = qiSnewLinear.PD(:,:, 1+hdelta: (zEnd-zStart+1) );
            if ~isempty(job.mtFiles) || (length(job.mtFiles)==1 &&~strcmp(job.mtFiles{1},'')), delta(:,:,zStart+hdelta:zEnd) = qiSnew.delta(:,:, 1+hdelta: (zEnd-zStart+1) ); end %
            if ~isempty(job.mtFiles) || (length(job.mtFiles)==1 &&~strcmp(job.mtFiles{1},'')), deltaLinear(:,:,zStart+hdelta:zEnd) = qiSnewLinear.delta(:,:, 1+hdelta: (zEnd-zStart+1) ); end %
            totalmask(:,:,zStart+hdelta:zEnd) = qiSnew.model.mask(:,:, 1+hdelta: (zEnd-zStart+1) );
            
            if job.confInt==1
                R1_low(:,:,zStart+hdelta:zEnd) = R1_1(:,:, 1+hdelta: (zEnd-zStart+1) );
                R1_up(:,:,zStart+hdelta:zEnd) = R1_2(:,:, 1+hdelta: (zEnd-zStart+1) );
                PD_low(:,:,zStart+hdelta:zEnd) = PD_2(:,:, 1+hdelta: (zEnd-zStart+1) );
                PD_up(:,:,zStart+hdelta:zEnd) = PD_1(:,:, 1+hdelta: (zEnd-zStart+1) );
            end
            
            if job.saveESTA==1        
               try
%               meta.invCov(:,:,:,:,zStart:zEnd)=modelMPM.invCov;
               meta.sigma2i(:,:,zStart:zEnd)=modelMPM.sigma2i;
               catch
               fprintf('There was a problem while saving sigma2i. Check to have enough free space and to be using at least version 7.3.');
               job.saveESTA = 0;      
               end        
            end
            
        end
        if zEnd==job.sdim(3)
                break;
        end
    %fprintf('Done the level starting at %d \n',startLayerVoxel);
    end
    spm_progress_bar('Set',job.sdim(3));
    spm_progress_bar('Clear');
    
    %% steps to save nii files
    % set all the NaN values to 0 (that makes easier to display correctly)
    R1(isnan(R1))=0;
    R1linear(isnan(R1linear))=0;
    R2star(isnan(R2star))=0;
    PD(isnan(PD))=0;
    PDlinear(isnan(PDlinear))=0;
    if ~isempty(job.mtFiles) || (length(job.mtFiles)==1 &&~strcmp(job.mtFiles{1},'')), delta(isnan(delta))=0; end
    if ~isempty(job.mtFiles) || (length(job.mtFiles)==1 &&~strcmp(job.mtFiles{1},'')), deltaLinear(isnan(deltaLinear))=0; end
    
    if job.confInt==1
        R1_low(isnan(R1_low))=0;
        R1_up(isnan(R1_up))=0;
        PD_low(isnan(PD_low))=0;
        PD_up(isnan(PD_up))=0;
    end
    
    % set all values of the voxels outside the mask to 0
    R1(totalmask<1)=0;
    R1linear(totalmask<1)=0;
    R2star(totalmask<1)=0;
    PD(totalmask<1)=0;
    PDlinear(totalmask<1)=0;
    if ~isempty(job.mtFiles) || (length(job.mtFiles)==1 &&~strcmp(job.mtFiles{1},'')), delta(totalmask<1)=0; end
    if ~isempty(job.mtFiles) || (length(job.mtFiles)==1 &&~strcmp(job.mtFiles{1},'')), deltaLinear(totalmask<1)=0; end
   
    if job.confInt==1
        R1_low(totalmask<1)=0;
        R1_up(totalmask<1)=0;
        PD_low(totalmask<1)=0;
        PD_up(totalmask<1)=0;
    end
    
    
    try
    % write 3 or 4 files for R1, PD, R2star and in case delta
    % function []= write_small_to_file_nii(outputdir,filenamepr, big_volume,small_volume_data,zStart, zEnd, sdim)
    write_small_to_file_nii(job.odir{1},'R1_', spm_vol(job.t1Files{1}), R1, 1, job.sdim(3), job.sdim);
    write_small_to_file_nii(job.odir{1},'R1linear_', spm_vol(job.t1Files{1}), R1linear, 1, job.sdim(3), job.sdim);
    write_small_to_file_nii(job.odir{1},'R2star_', spm_vol(job.t1Files{1}), R2star, 1, job.sdim(3), job.sdim);
    write_small_to_file_nii(job.odir{1},'A_', spm_vol(job.t1Files{1}), PD, 1, job.sdim(3), job.sdim);
    write_small_to_file_nii(job.odir{1},'Alinear_', spm_vol(job.t1Files{1}), PDlinear, 1, job.sdim(3), job.sdim);
    [pA, fA, eA] = spm_fileparts(job.t1Files{1});
    
    if job.confInt==1
        write_small_to_file_nii(job.odir{1},'R1_low_', spm_vol(job.t1Files{1}), R1_low, 1, job.sdim(3), job.sdim);
        write_small_to_file_nii(job.odir{1},'R1_up_', spm_vol(job.t1Files{1}), R1_up, 1, job.sdim(3), job.sdim);
        write_small_to_file_nii(job.odir{1},'A_low_', spm_vol(job.t1Files{1}), PD_low, 1, job.sdim(3), job.sdim);
        write_small_to_file_nii(job.odir{1},'A_up_', spm_vol(job.t1Files{1}), PD_up, 1, job.sdim(3), job.sdim);
        
    end
    eA = '.nii';
    if ~isempty(job.mtFiles) || (length(job.mtFiles)==1 &&~strcmp(job.mtFiles{1},''))  
        write_small_to_file_nii(job.odir{1},'MT_', spm_vol(job.t1Files{1}), delta, 1, job.sdim(3),  job.sdim); 
        write_small_to_file_nii(job.odir{1},'MTlinear_', spm_vol(job.t1Files{1}), deltaLinear, 1, job.sdim(3),  job.sdim); 
        PDcalculation(fullfile(job.odir{1}, ['MTlinear_' fA eA]), fullfile(job.odir{1}, ['Alinear_' fA eA]));
        PDcalculation(fullfile(job.odir{1}, ['MT_' fA eA]), fullfile(job.odir{1}, ['A_' fA eA]));
    else
        PDcalculation(fullfile(job.odir{1}, ['R1_' fA eA]), fullfile(job.odir{1}, ['A_' fA eA]));
        PDcalculation(fullfile(job.odir{1}, ['R1linear_' fA eA]), fullfile(job.odir{1}, ['Alinear_' fA eA]));
    end
    catch
        error('it was not possible to save the resulting .nii files. Check to have writing rights in the current directory')
    end
   % save mask in the .mat (in case we are saving the model)
    if job.saveESTA==1,   
%               clear invCov;
               try
               meta.mask=totalmask;
                
                if job.confInt==1
                    meta.R1_low = {fullfile(job.odir{1},strcat('R1_low_',job.t1Files{1}))};
                    meta.R1_up = {fullfile(job.odir{1},strcat('R1_up_',job.t1Files{1}))};
                    meta.PD_low = {fullfile(job.odir{1},strcat('A_low_',job.t1Files{1}))};
                    meta.PD_up = {fullfile(job.odir{1},strcat('A_up_',job.t1Files{1}))};
                end
               
               catch
               fprintf('There was a problem saving the ESTATICS model. Check to have enough free space and to be using at least version 7.3.');
               job.saveESTA = 0;      
               end        
    end
    fprintf('Ending the MPM at %s \n',datestr(now));

end