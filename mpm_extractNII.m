% function to get the nii files of the estimated parameter given an ESTATICS model
%
%
% =========================================================================
% 2016/01/25
%
% written by C. D'Alonzo 
%
% function [] = mpm_extractNII(job)
% 
% produces the nii files of the parameters, iterating on each region containing "heigth"
% planes x,y. 
%
% Input:
%  job             - a struct containing:
%  
%  ESTAmodel       - a file containing the metadata
%  height          - height of the region of interest 
%  odir            - output directory   
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
%  modelCoeff      - a 3 or 4 dimensional array containing the values of
%                    the parameters
%  TR
%  TE              - arrays containing the relaxation rate, flip angle and echo time  
%  FA
%
% no Output
% the 3 or 4 parameter are written each in a .nii file
%                 
%  
%
% =========================================================================

function [] = mpm_extractNII(job)

    fprintf('Starting the function at %s \n',datestr(now));


 %% works on the all volume, level by level
    
    fprintf('Loading data...\n');
    meta=matfile(char(job.ESTAmodel));
    t1Files = meta.t1Files;
    V = spm_vol(t1Files{1});
    [~,fname,~] = spm_fileparts(V.fname);
    
%     coeffFiles = meta.modelCoeff; %
    nv = meta.nv;
%     modelMPM.nFiles = meta.nFiles;
%     modelMPM.t1Files = meta.t1Files;
%     modelMPM.pdFiles = meta.pdFiles;
%     if modelMPM.nv==4,
%          modelMPM.mtFiles = meta.mtFiles;
%     end
%     modelMPM.FA = meta.FA;
%     modelMPM.TE = meta.TE;
%     modelMPM.TR = meta.TR;
    TEScale = meta.TEScale;
    DataScale = meta.DataScale;
    
%     Vt1 = spm_vol(modelMPM.t1Files{1}); 
%     Vpd = spm_vol(modelMPM.pdFiles{1});  
%     if modelMPM.nv==4, Vmt = spm_vol(modelMPM.mtFiles{1}); end 
   
    sdim = meta.sdim;
   
 % S_t1         
         wV1        = V; 
         dt        = [spm_type('float32'),spm_platform('bigend')];
         dm        = V(1).dim;
         Ni1        = nifti;
         Ni1.mat    = V(1).mat;
         Ni1.mat0    = V(1).mat;
         wV1.fname    = fullfile(job.odir{1},['ESTATICS_S_t1_' fname '.nii']);
         Ni1.dat      = file_array(wV1.fname,dm,dt, 0,1,0);
         create(Ni1);
         zeroValues = zeros(sdim);
         for i=1:sdim(3),
         Ni1.dat(:,:,i) = zeroValues(:,:,i);
         end
         
         % S_pd         
         wV2        = V;          
         Ni2        = nifti;
         Ni2.mat    = V(1).mat;
         Ni2.mat0    = V(1).mat;
         wV2.fname    = fullfile(job.odir{1},['ESTATICS_S_pd_' fname '.nii']);
         Ni2.dat      = file_array(wV2.fname,dm,dt, 0,1,0);
         create(Ni2);
         for i=1:sdim(3),
         Ni2.dat(:,:,i) = zeroValues(:,:,i);
         end
         
         if nv==4,
                 % S_mt                
                wV3        = V;          
                Ni3        = nifti;
                Ni3.mat    = V(1).mat;
                Ni3.mat0    = V(1).mat;
                wV3.fname    = fullfile(job.odir{1},['ESTATICS_S_mt_' fname '.nii']);
                Ni3.dat      = file_array(wV3.fname,dm,dt, 0,1,0);
                create(Ni3);
                for i=1:sdim(3),
                   Ni3.dat(:,:,i) = zeroValues(:,:,i);
                end
         end
         
         % R2star  
                wV4        = V;          
                Ni4        = nifti;
                Ni4.mat    = V(1).mat;
                Ni4.mat0    = V(1).mat;
                wV4.fname    = fullfile(job.odir{1},['ESTATICS_R2star_' fname '.nii']);
                Ni4.dat      = file_array(wV4.fname,dm,dt, 0,1,0);
                create(Ni4);
                for i=1:sdim(3),
                Ni4.dat(:,:,i) = zeroValues(:,:,i);
                end
    
    % check for the input height and adjust it if too big or too small
    if abs(job.height)>sdim(3),
        job.height = sdim(3);
    elseif job.height<0,
        job.height = abs(job.height);
    end
    interval = fix(job.height);
   
    
    
    
   
    
    spm_progress_bar('Init',sdim(3),'planes completed');
    
    %% start iteration on all the levels
    for startLayerVoxel = 1:interval:sdim(3),
        
        spm_progress_bar('Set',startLayerVoxel);
        
        zStart = double(startLayerVoxel);        
        
        if sdim(3)-(startLayerVoxel + interval)> job.height,
            % in case the next starting point has enough planes after it
            zEnd = double(startLayerVoxel + interval); 
        else
            % in case the next starting point has NOT enough planes after it
            zEnd = sdim(3);
        end
%         fprintf('on level %d - %d\n', zStart, zEnd);
        
        %% read the mask and model coefficients between zStart and zEnd
            slices=zStart:zEnd; %1:sdim(3);
            %[mask(:,:,:),~] = loadImageSPM(fullfile(dir,[maskFile{1} '.nii']),'slices',slices);
            try
            mask=meta.mask(:,:,slices);
            
                           
            modelCoeff(1,:,:,:) = meta.modelCoeff(1,:,:,slices);
            
            modelCoeff(2,:,:,:) = meta.modelCoeff(2,:,:,slices);
            
               
            modelCoeff(3,:,:,:) = meta.modelCoeff(3,:,:,slices);
            if nv==4,  
                modelCoeff(4,:,:,:) = meta.modelCoeff(4,:,:,slices);            
            end
            
            modelCoeff(:,mask<1) = 0;
            Ni1.dat(:,:,slices) = reshape(modelCoeff(1,:,:,:),[sdim(1), sdim(2),length(slices)]).*DataScale; % zEnd-zStart+1
            Ni2.dat(:,:,slices) = reshape(modelCoeff(2,:,:,:),[sdim(1), sdim(2), length(slices)]).*DataScale;
              if nv==4,
                  Ni3.dat(:,:,slices) = reshape(modelCoeff(3,:,:,:),[sdim(1), sdim(2),length(slices)]).*DataScale;
                  Ni4.dat(:,:,slices) = reshape(modelCoeff(4,:,:,:),[sdim(1), sdim(2),length(slices)])./TEScale;
              else
                   
                  Ni4.dat(:,:,slices) = reshape(modelCoeff(3,:,:,:),[sdim(1), sdim(2), length(slices)])./TEScale;
              end
            clear modelCoeff mask;
            
            catch ME
            fprintf(ME.message);
            fprintf('\nThere was a problem reading the ESTATICS model and/or saving the nii files. Check to have enough free space and to be using at least version 7.3.\n');
            
            end        
        
        
        
        
        if zEnd==sdim(3),
                break;
        end
    
    end
    spm_progress_bar('Set',sdim(3));
    spm_progress_bar('Clear');
    
    
    fprintf('Ending the function at %s \n',datestr(now));

end