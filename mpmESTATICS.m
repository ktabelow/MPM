function [] = mpmESTATICS(job)

 %% works on the all cubus, levels by levels
    
    % calculates how big the overlapping has to be
    % to assure a good smoothing    
    hmax = 1.25^(job.kstar/3);
    hakt = gethani (1, 1.25*hmax, 2, 1.25^job.kstar, [1 1], 1e-4);
    hdelta = ceil(hakt); % height of half of the overlapping
    % calculates the interval between the zStart of the different levels
    interval = int64(job.height - 2*hdelta);

    % preparing result variables
    R1 = zeros(job.sdim);
    R2star = zeros(job.sdim);
    PD = zeros(job.sdim);
    delta = zeros(job.sdim);

    % prepares variable to save the whole mask
    totalmask = zeros(job.sdim);
    
    % looking for weights
    try 
        wghts=getWeights(job.t1Files{1});
    catch 
    wghts = [];
    end
    
    fprintf('Starting the MPM at %s \n',datestr(now));
    spm_progress_bar('Init',job.sdim(3),'planes completed');
    
    % start iteration on all the levels
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
 
        % prepares the data        
        dataset = createDataSet(job.sdim,zStart, zEnd, job.t1Files,job.pdFiles,job.mtFiles,char(job.maskFile),job.t1TR,job.pdTR,job.mtTR,job.t1TE,job.pdTE,job.mtTE,job.t1FA,job.pdFA, job.mtFA);
        
        % estimate from all the data 4 parameter 
        % function [model] = estimateESTATICS(dataset, varargin)
        modelMPM3 = estimateESTATICS(dataset,'verbose', true);
        
        % smooth and calculate the final four parameters with the functions
        % function [modelS] = smoothESTATICSmask(model, varargin)
        % function [qi] = calculateQI(model, varargin)
        % (if kstar is 0, the smoothing step isn't done)
        if job.kstar~=0, 
            modelMPM3snew = smoothESTATICSmask(modelMPM3, 'verbose', false, 'wghts', wghts); 
            qiSnew = calculateQI(modelMPM3snew, 'TR2',job.tr2,'b1File',job.b1File , 'verbose', false);
        else
            qiSnew = calculateQI(modelMPM3, 'TR2',job.tr2,'b1File',job.b1File , 'verbose', false);
        end

        
        
        if zStart==1
            R1(:,:,zStart:zEnd) = qiSnew.R1;
            R2star(:,:,zStart:zEnd) = qiSnew.R2star;
            PD(:,:,zStart:zEnd) = qiSnew.PD;
            if ~strcmp(job.mtFiles{1},''), delta(:,:,zStart:zEnd) = qiSnew.delta; end%
            totalmask(:,:,zStart:zEnd) = qiSnew.model.mask;
        else 
                R1(:,:,zStart+hdelta:zEnd) = qiSnew.R1(:,:, 1+hdelta: (zEnd-zStart+1) );
            R2star(:,:,zStart+hdelta:zEnd) = qiSnew.R2star(:,:, 1+hdelta: (zEnd-zStart+1) );
            PD(:,:,zStart+hdelta:zEnd) = qiSnew.PD(:,:, 1+hdelta: (zEnd-zStart+1) );
            if ~strcmp(job.mtFiles{1},''), delta(:,:,zStart+hdelta:zEnd) = qiSnew.delta(:,:, 1+hdelta: (zEnd-zStart+1) ); end %
            totalmask(:,:,zStart+hdelta:zEnd) = qiSnew.model.mask(:,:, 1+hdelta: (zEnd-zStart+1) );
        end
        if zEnd==job.sdim(3),
                break;
        end
    
    end
    spm_progress_bar('Set',job.sdim(3));
    spm_progress_bar('Clear');
    
    
    % set all the NaN values to 0 (that makes easier to display correctly)
    R1(isnan(R1))=0;
    R2star(isnan(R2star))=0;
    PD(isnan(PD))=0;
    if ~strcmp(job.mtFiles{1},''), delta(isnan(delta))=0; end
    
    % set all values of the voxels outside the mask to 0
    R1(totalmask<1)=0;
    R2star(totalmask<1)=0;
    PD(totalmask<1)=0;
    if ~strcmp(job.mtFiles{1},''), delta(totalmask<1)=0; end
   
    % write 3 or 4 files for R1, PD, R2star and in case delta
    % function []= write_small_to_file_nii(outputdir,filenamepr, big_volume,small_volume_data,zStart, zEnd, sdim)
    write_small_to_file_nii(pwd,'R1_', spm_vol(job.t1Files{1}), R1, 1, job.sdim(3), job.sdim);
    write_small_to_file_nii(pwd,'R2star_', spm_vol(job.t1Files{1}), R2star, 1, job.sdim(3), job.sdim);
    write_small_to_file_nii(pwd,'PD_', spm_vol(job.pdFiles{1}), PD, 1, job.sdim(3), job.sdim);
    if ~strcmp(job.mtFiles{1},''),  write_small_to_file_nii(pwd,'delta_', spm_vol(job.mtFiles{1}), delta, 1, job.sdim(3),  job.sdim); end
    
    fprintf('Ending the MPM at %s \n',datestr(now));

end