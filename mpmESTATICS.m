% wrapper function for all the steps of the MPM toolbox
%
% =========================================================================
% 2015/12/08
%
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
%  tr2 
%  height          - height of the region of interest 
%  kstar           - number of steps for the smoothing algorithm (if 0 no smoothing)          
%  lambda          - adaptation bandwidth
%  maskFile        - mask containing the voxel to consider
%  b1File          - correction field
%  t1TR mtTR pdTR  
%  t1TE mtTE pdTE  - parameter (read from file if empty)
%  t1FA mtFA pdFA
%
% no Output
% the 3 or 4 final variables R1, PD, R2star and in case of a model with MT also delta
% are written each in a .nii file
%                 
%  
%
% =========================================================================

function [] = mpmESTATICS(job)

    fprintf('Starting the MPM at %s \n',datestr(now));


 %% works on the all volume, level by level
    
    V = spm_vol(job.t1Files{1});
    
    % read the dimension array if not given
    if sum(job.sdim)==0,      
      job.sdim = V.dim;
    end  
 
    % calculates how big the overlapping has to be
    % to assure a good smoothing    
    hmax = 1.25^(job.kstar/3);
    hakt = gethani (1, 1.25*hmax, 2, 1.25^job.kstar, [1 1], 1e-4);
    hdelta = ceil(hakt); % height of half of the overlapping
    
    % check for the input height and adjust it if too big or too small
    if abs(job.height)>job.sdim(3),
        job.height = job.sdim(3);
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
    R1 = zeros(job.sdim);
    R2star = zeros(job.sdim);
    PD = zeros(job.sdim);
    delta = zeros(job.sdim);

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
        
    
    % if saving ESTATICS model, prepare the .mat with the metadata
    % 
    if job.saveESTA==1,
        % prepare to save the whole invCov directly in the .mat file
        invCov=zeros([dataset.nv dataset.nv job.sdim]);
        try
            [~,fname,~] = spm_fileparts(V.fname);
            filename = strcat('ESTATICS_model_',fname,'.mat');
            save(filename,'invCov','-v7.3');
            meta = matfile(filename,'Writable',true);
         % S_t1         
%          wV1        = V; 
%          dt        = [spm_type('float32'),spm_platform('bigend')];
%          dm        = V(1).dim;
%          Ni1        = nifti;
%          Ni1.mat    = V(1).mat;
%          Ni1.mat0    = V(1).mat;
%          wV1.fname    = fullfile(pwd,['ESTATICS_S_t1_' fname '.nii']);
%          Ni1.dat      = file_array(wV1.fname,dm,dt, 0,1,0);
%          create(Ni1);
%          zeroValues = zeros(dataset.sdim);
%          for i=1:dataset.sdim(3),
%          Ni1.dat(:,:,i) = zeroValues(:,:,i);
%          end
         
         % S_pd         
%          wV2        = V;          
%          Ni2        = nifti;
%          Ni2.mat    = V(1).mat;
%          Ni2.mat0    = V(1).mat;
%          wV2.fname    = fullfile(pwd,['ESTATICS_S_pd_' fname '.nii']);
%          Ni2.dat      = file_array(wV2.fname,dm,dt, 0,1,0);
%          create(Ni2);
%          for i=1:dataset.sdim(3),
%          Ni2.dat(:,:,i) = zeroValues(:,:,i);
%          end
%          
%          if dataset.nv==4,
                 % S_mt                
%                 wV3        = V;          
%                 Ni3        = nifti;
%                 Ni3.mat    = V(1).mat;
%                 Ni3.mat0    = V(1).mat;
%                 wV3.fname    = fullfile(pwd,['ESTATICS_S_mt_' fname '.nii']);
%                 Ni3.dat      = file_array(wV3.fname,dm,dt, 0,1,0);
%                 create(Ni3);
%                 for i=1:dataset.sdim(3),
%                    Ni3.dat(:,:,i) = zeroValues(:,:,i);
%                 end
%          end
         
         % R2star  
%                 wV4        = V;          
%                 Ni4        = nifti;
%                 Ni4.mat    = V(1).mat;
%                 Ni4.mat0    = V(1).mat;
%                 wV4.fname    = fullfile(pwd,['ESTATICS_R2star_' fname '.nii']);
%                 Ni4.dat      = file_array(wV4.fname,dm,dt, 0,1,0);
%                 create(Ni4);
%                 for i=1:dataset.sdim(3),
%                 Ni4.dat(:,:,i) = zeroValues(:,:,i);
%                 end
        % add the name of the file to the .mat file
%         if dataset.nv==4,
%             meta.modelCoeff = {wV1.fname, wV2.fname, wV3.fname, wV4.fname};
%         else
%             meta.modelCoeff = {wV1.fname, wV2.fname, wV4.fname};
%         end
        
        catch ME
            fprintf(ME.message);
            fprintf('\nIt is not possible to save the ESTATICS model. Check to have writing rights in the current directory, enough memory and to be using at least version 7.3.');
            job.saveESTA = 0;
      
        end
        clear invCov;
    end
    
    
    % save metadata for the ESTATICS model
    if job.saveESTA==1,
        
           try
           meta.sdim = dataset.sdim;
           meta.t1Files = dataset.t1Files;
           meta.pdFiles = dataset.pdFiles;
           if dataset.nv == 4,
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
           meta.modelCoeff = zeros([dataset.nv dataset.sdim]);
           catch ME
            fprintf(ME.message);   
            fprintf('\nThere was a problem saving the ESTATICS model. Check to have writing rights in the current directory and to be using at least version 7.3.\n');
            job.saveESTA = 0;      
           end
        
    end
    
    
    spm_progress_bar('Init',job.sdim(3),'planes completed');
    
    %% start iteration on all the levels
    for startLayerVoxel = 1:interval:job.sdim(3),
        
        spm_progress_bar('Set',startLayerVoxel);
        
        zStart = double(startLayerVoxel);        

        if job.sdim(3)-(startLayerVoxel + job.height)> 2*hdelta,
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
            if zStart==1,
                V1 = spm_vol(dataset.maskFile);
                maskdim=V1.dim;
            end
            if dataset.sdim~=maskdim,
                dataset.mask=ones([dataset.sdim(1) dataset.sdim(2) zEnd-zStart+1]);
                if zStart==1, fprintf('no correct dimension in mask file - working on all voxels\n'); end
            else
                if zStart==1, fprintf('correct dimension in mask file\n'); end
                slices=zStart:zEnd; %1:sdim(3);
                %[mask(:,:,:),~] = loadImageSPM(fullfile(dir,[maskFile{1} '.nii']),'slices',slices);
                [mask(:,:,:),~] = loadImageSPM(fullfile(dataset.maskFile) ,'slices',slices);
                mask = round(mask(:));
                mask = reshape (mask, [dataset.sdim(1) dataset.sdim(2) zEnd-zStart+1]);
                dataset.mask=mask;
                clear mask;
            end
        end
        
        
        %% estimate from all the data the 4or3 parameter 
        % function [model] = estimateESTATICS(dataset, varargin)
        modelMPM = estimateESTATICS(dataset,'verbose', false, 'tolerance', job.tol);
        
        %% save the 4or3 parameters of the ESTATICS model and DataScale and TEScale in the .mat
        if job.saveESTA==1,        
         try
           meta.TEScale = modelMPM.TEScale;
           meta.DataScale = modelMPM.DataScale;
          if zStart==1,
              fprintf('saving coefficient on level zstart\n');
%               size(meta.modelCoeff)
%               size(modelMPM.modelCoeff)
%                size(meta.modelCoeff(1,:,:,1:zEnd))
%                size(modelMPM.modelCoeff(1,:,:,:))
              meta.modelCoeff(1,:,:,1:zEnd) = modelMPM.modelCoeff(1,:,:,:);
              meta.modelCoeff(2,:,:,1:zEnd) = modelMPM.modelCoeff(2,:,:,:);
%               Ni1.dat(:,:,1:zEnd) = reshape(modelMPM.modelCoeff(1,:,:,:),[job.sdim(1), job.sdim(2), zEnd]);
%               Ni2.dat(:,:,1:zEnd) = reshape(modelMPM.modelCoeff(2,:,:,:),[job.sdim(1), job.sdim(2), zEnd]);
              if dataset.nv==4,
                  meta.modelCoeff(3,:,:,1:zEnd) = modelMPM.modelCoeff(3,:,:,:);
                  meta.modelCoeff(4,:,:,1:zEnd) = modelMPM.modelCoeff(4,:,:,:);
%                   Ni3.dat(:,:,1:zEnd) = reshape(modelMPM.modelCoeff(3,:,:,:),[job.sdim(1), job.sdim(2), zEnd]);
%                   Ni4.dat(:,:,1:zEnd) = reshape(modelMPM.modelCoeff(4,:,:,:),[job.sdim(1), job.sdim(2), zEnd]);
              else
                   meta.modelCoeff(3,:,:,1:zEnd) = modelMPM.modelCoeff(3,:,:,:);
%                   Ni4.dat(:,:,1:zEnd) = reshape(modelMPM.modelCoeff(3,:,:,:),[job.sdim(1), job.sdim(2), zEnd]);
              end
            
          else 
              fprintf('saving coefficient on level %d \n',zStart);
              meta.modelCoeff(1,:,:,zStart+hdelta:zEnd) = modelMPM.modelCoeff(1,:,:, 1+hdelta: (zEnd-zStart+1)); 
              meta.modelCoeff(2,:,:,zStart+hdelta:zEnd) = modelMPM.modelCoeff(2,:,:, 1+hdelta: (zEnd-zStart+1)); 
%               Ni1.dat(:,:,zStart+hdelta:zEnd) = reshape(modelMPM.modelCoeff(1,:,:, 1+hdelta: (zEnd-zStart+1)),[job.sdim(1), job.sdim(2), zEnd-zStart+1-hdelta] ); 
%               Ni2.dat(:,:,zStart+hdelta:zEnd) = reshape(modelMPM.modelCoeff(2,:,:, 1+hdelta: (zEnd-zStart+1)),[job.sdim(1), job.sdim(2), zEnd-zStart+1-hdelta] ); 
              if dataset.nv==4,
                  meta.modelCoeff(3,:,:,zStart+hdelta:zEnd) = modelMPM.modelCoeff(3,:,:, 1+hdelta: (zEnd-zStart+1)); 
                  meta.modelCoeff(4,:,:,zStart+hdelta:zEnd) = modelMPM.modelCoeff(4,:,:, 1+hdelta: (zEnd-zStart+1)) ; 
%                   Ni3.dat(:,:,zStart+hdelta:zEnd) = reshape(modelMPM.modelCoeff(3,:,:,1+hdelta: (zEnd-zStart+1)),[job.sdim(1), job.sdim(2), zEnd-zStart+1-hdelta] ); 
%                   Ni4.dat(:,:,zStart+hdelta:zEnd) = reshape(modelMPM.modelCoeff(4,:,:,1+hdelta: (zEnd-zStart+1)),[job.sdim(1), job.sdim(2), zEnd-zStart+1-hdelta] ); 
              else
                  meta.modelCoeff(3,:,:,zStart+hdelta:zEnd) = modelMPM.modelCoeff(3,:,:, 1+hdelta: (zEnd-zStart+1)); 
%                   Ni4.dat(:,:,zStart+hdelta:zEnd) = reshape(modelMPM.modelCoeff(3,:,:,1+hdelta: (zEnd-zStart+1)),[job.sdim(1), job.sdim(2), zEnd-zStart+1-hdelta] ); 
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
            if ~isempty(job.mtFiles) || (length(job.mtFiles)==1 &&~strcmp(job.mtFiles{1},'')), delta(:,:,zStart:zEnd) = qiSnew.delta; end%
            totalmask(:,:,zStart:zEnd) = qiSnew.model.mask;
            
            if job.saveESTA==1,        
               try
               meta.invCov(:,:,:,:,zStart:zEnd)=modelMPM.invCov;
               catch ME
                   fprintf(ME.message);                   
                   fprintf('\nThere was a problem saving the ESTATICS model. Check to have enough free space and to be using at least version 7.3.\n');
                   job.saveESTA = 0;      
               end        
            end          
            
        else 
            R1(:,:,zStart+hdelta:zEnd) = qiSnew.R1(:,:, 1+hdelta: (zEnd-zStart+1) );
            R2star(:,:,zStart+hdelta:zEnd) = qiSnew.R2star(:,:, 1+hdelta: (zEnd-zStart+1) );
            PD(:,:,zStart+hdelta:zEnd) = qiSnew.PD(:,:, 1+hdelta: (zEnd-zStart+1) );
            if ~isempty(job.mtFiles) || (length(job.mtFiles)==1 &&~strcmp(job.mtFiles{1},'')), delta(:,:,zStart+hdelta:zEnd) = qiSnew.delta(:,:, 1+hdelta: (zEnd-zStart+1) ); end %
            totalmask(:,:,zStart+hdelta:zEnd) = qiSnew.model.mask(:,:, 1+hdelta: (zEnd-zStart+1) );
            
            if job.saveESTA==1,        
               try
               meta.invCov(:,:,:,:,zStart+hdelta:zEnd)=modelMPM.invCov(:,:,:,:,1+hdelta: (zEnd-zStart+1) );
               catch
               fprintf('There was a problem saving the ESTATICS model. Check to have enough free space and to be using at least version 7.3.');
               job.saveESTA = 0;      
               end        
            end
            
        end
        if zEnd==job.sdim(3),
                break;
        end
    
    end
    spm_progress_bar('Set',job.sdim(3));
    spm_progress_bar('Clear');
    
    %% steps to save nii files
    % set all the NaN values to 0 (that makes easier to display correctly)
    R1(isnan(R1))=0;
    R2star(isnan(R2star))=0;
    PD(isnan(PD))=0;
    if ~isempty(job.mtFiles) || (length(job.mtFiles)==1 &&~strcmp(job.mtFiles{1},'')), delta(isnan(delta))=0; end
    
    % set all values of the voxels outside the mask to 0
    R1(totalmask<1)=0;
    R2star(totalmask<1)=0;
    PD(totalmask<1)=0;
    if ~isempty(job.mtFiles) || (length(job.mtFiles)==1 &&~strcmp(job.mtFiles{1},'')), delta(totalmask<1)=0; end
   
    try
    % write 3 or 4 files for R1, PD, R2star and in case delta
    % function []= write_small_to_file_nii(outputdir,filenamepr, big_volume,small_volume_data,zStart, zEnd, sdim)
    write_small_to_file_nii(pwd,'R1_', spm_vol(job.t1Files{1}), R1, 1, job.sdim(3), job.sdim);
    write_small_to_file_nii(pwd,'R2star_', spm_vol(job.t1Files{1}), R2star, 1, job.sdim(3), job.sdim);
    write_small_to_file_nii(pwd,'PD_', spm_vol(job.pdFiles{1}), PD, 1, job.sdim(3), job.sdim);
    if ~isempty(job.mtFiles) || (length(job.mtFiles)==1 &&~strcmp(job.mtFiles{1},'')),  write_small_to_file_nii(pwd,'delta_', spm_vol(job.mtFiles{1}), delta, 1, job.sdim(3),  job.sdim); end
    catch
        error('it was not possible to save the resulting .nii files. Check to have writing rights in the current directory')
    end
   % save mask in the .mat (in case we are saving the model)
    if job.saveESTA==1,   
               clear invCov;
               try
               meta.mask=totalmask;
               catch
               fprintf('There was a problem saving the ESTATICS model. Check to have enough free space and to be using at least version 7.3.');
               job.saveESTA = 0;      
               end        
    end
    fprintf('Ending the MPM at %s \n',datestr(now));

end