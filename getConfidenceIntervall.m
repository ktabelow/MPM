% calculate the confidence interval for R1 and PD of the ESTATICSmodel
%
% =========================================================================
% 2016/07/27
%
% written by C. D'Alonzo 
%
%
% getConfidenceIntervall(model,varargin)
% 
% see preprint .... ISSN ... for details
%
% Input:
%
%  model     - ESTATIC model (as struct)
%
% Optional Arguments:
%
%  b1File   - filename of the correction field
%  alpha    - confidence level,  default 0.05
%
% Output:
%   
%  model     - with fields 
% TO DO
%
% =========================================================================

function [R1_1,R1_2,PD_1,PD_2] = getConfidenceIntervall(model,varargin)

%% sets the default parameters
b1File  = {};

alpha   = 0.05;

%% overwrites default parameters, if present

for k=1:2:length(varargin),     
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;


c = 4*finv(1-alpha,model.nv,model.nFiles-model.nv);

s11 = reshape(model.invCov(1,1,:,:,:),[model.sdim(1) model.sdim(2) (model.zEnd-model.zStart+1)]);  
s12 = reshape(model.invCov(1,2,:,:,:),[model.sdim(1) model.sdim(2) (model.zEnd-model.zStart+1)]);
s22 = reshape(model.invCov(2,2,:,:,:),[model.sdim(1) model.sdim(2) (model.zEnd-model.zStart+1)]);
coeff1 = reshape(model.modelCoeff(1,:,:,:),[model.sdim(1) model.sdim(2) (model.zEnd-model.zStart+1)]);
coeff2 = reshape(model.modelCoeff(2,:,:,:),[model.sdim(1) model.sdim(2) (model.zEnd-model.zStart+1)]);


[m,M] = mpm_get_theta2(s11,s12,s22,c,coeff1,coeff2);

[theta1a_m,theta1a_M] = mpm_get_theta1(s11,s12,s22,c,coeff1,m, coeff2);
[theta1b_m,theta1b_M] = mpm_get_theta1(s11,s12,s22,c,coeff1,M, coeff2);

minTheta1 = min(theta1a_m,theta1b_m);
maxTheta1 = max(theta1a_M,theta1b_M);


[R1_1,PD_1] = calculate_R1_PD(model,minTheta1,'b1File',b1File); % change to model
[R1_2,PD_2] = calculate_R1_PD(model,maxTheta1,'b1File',b1File);




function [minTheta2, maxTheta2] = mpm_get_theta2(s11,s12,s22,c,theta1hat,theta2hat)

Q_a = (s12./s11).^2 - s22./s11;
  
p_a = theta1hat + s12./s11 .* theta2hat;

den = Q_a.^2 .* theta2hat.^2 - p_a.^2 .* Q_a;

q = ((c./s11).^2 - p_a.^2*c./s11) ./ den;

p = -2 * Q_a .*theta2hat*c ./ (s11 .* den);

theta2pos = theta2hat - p/2 + sqrt(p.^2/4 - q);

theta2neg = theta2hat - p/2 - sqrt(p.^2/4 - q);

minTheta2 = min(theta2pos, theta2neg);

maxTheta2 = max(theta2pos, theta2neg);



function [minTheta1, maxTheta1] = mpm_get_theta1(s11,s12,s22,c,theta1hat,theta2,theta2hat)

delta2 = theta2 - theta2hat;

Q_a = (s12./s11).^2 - s22./s11;

theta1pos = theta1hat - s12./s11.*delta2 + sqrt(Q_a.*delta2.^2 + c/s11);

theta1neg = theta1hat - s12./s11.*delta2 - sqrt(Q_a.*delta2.^2 + c/s11);

minTheta1 = min(theta1pos, theta1neg);
maxTheta1 = max(theta1pos, theta1neg);



function [R1,PD]= calculate_R1_PD(model,theta1,varargin)
%% read B1 correction field

%% sets the default parameters
b1File  = {};



%% overwrites default parameters, if present

for k=1:2:length(varargin),     
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

slices = model.zStart : model.zEnd; %1:model.sdim(3);
if (length(b1File)==1 && strcmp(b1File{1},''))
    %if model.zStart==1, fprintf('no B1 correction file - no correction will be performed\n'); end
    %if verbose, fprintf('no B1 correction\n'); end
    b1Map = ones([model.sdim(1) model.sdim(2) (model.zEnd-model.zStart+1)]);
elseif ~isempty(b1File)   % || model.sdim==b1dim
   V1 = spm_vol(b1File{1});
   b1dim=V1.dim;
    if model.sdim==b1dim,
        %if model.zStart==1, fprintf('reading correction file from %s \n',b1File{1}); end
        %if verbose, fprintf('reading correction file from %s \n',b1File{1}); end   
        %[b1Map(:,:,:),~] = loadImageSPM(fullfile(b1File{1}),'slices',slices);
        b1Map(:,:,:) = MPM_read_coregistered_vol(spm_vol(b1File{1}),spm_vol(model.pdFiles{1}),'slices',slices);
        b1Map = b1Map./100; % still to check if is ok
        b1Map(b1Map < 0) = 0;
    else
        %if model.zStart==1, fprintf('no correct dimension in B1 correction file - no correction will be performed\n'); end
        %if verbose, fprintf('no correct B1 correction file\n'); end
        b1Map = ones([model.sdim(1) model.sdim(2) (model.zEnd-model.zStart+1)]);
    end
else 
    %if model.zStart==1, fprintf('no B1 correction file selected - no correction will be performed\n'); end
    %if verbose, fprintf('no B1 correction\n'); end
    b1Map = ones([model.sdim(1) model.sdim(2) (model.zEnd-model.zStart+1)]);
end

%% get correct flip angles an TR times

t1FA = model.FA(1);
if model.nv==4, 
    pdFA = model.FA(length(model.t1Files)+length(model.mtFiles)+1); 
else
    pdFA = model.FA(length(model.t1Files)+1);
end
    
t1TR = model.TR(1);



%% calculate E1


alphat1 = b1Map.*t1FA/180*pi; 
alphapd = b1Map.*pdFA/180*pi;
SINalphat1 = sin(alphat1);
COSalphat1 = cos(alphat1);
SINalphapd = sin(alphapd);
COSalphapd = cos(alphapd);

clear alphat1 alphapd;

if model.nv == 4 
 enum = theta1 - SINalphat1./SINalphapd.* reshape(model.modelCoeff(3,:,:,:),[model.sdim(1) model.sdim(2) (model.zEnd-model.zStart+1)]);
 denom = theta1.*COSalphat1 - SINalphat1./SINalphapd.* reshape(model.modelCoeff(3,:,:,:),[model.sdim(1) model.sdim(2) (model.zEnd-model.zStart+1)]).* COSalphapd;
else 
    enum = theta1 - SINalphat1./SINalphapd.* reshape(model.modelCoeff(2,:,:,:),[model.sdim(1) model.sdim(2) (model.zEnd-model.zStart+1)]);
    denom = theta1.*COSalphat1 - SINalphat1./SINalphapd.* reshape(model.modelCoeff(2,:,:,:),[model.sdim(1) model.sdim(2) (model.zEnd-model.zStart+1)]).* COSalphapd;
end

E1 = enum./denom;

clear enum denom COSalphapd SINalphapd;

R1 = -log(E1)/t1TR;

R1(imag(R1)>0)=nan; %
R1 = real(R1);

% RF spoiling correction Preibisch and Deichmann MRM 61 (2009) 125-135
% These coefficients depend on the sequence!! See getPolynomsP2_ab and
% MTprot in VBQ
P2_a = model.P2_a;
P2_b = model.P2_b;
R1 = R1 ./ ((P2_a(1)*b1Map.^2 + P2_a(2)*b1Map + P2_a(3)) .* R1 + (P2_b(1)*b1Map.^2 + P2_b(2)*b1Map + P2_b(3)));
E1 = exp(- R1 * t1TR);
% END spoiling correction

R1 = R1 * 1000; % TR is in ms, so R1 is now in s^(-1)



%% calculate PD


enum = (1 - COSalphat1 .* E1).*reshape(model.modelCoeff(1,:,:,:),[model.sdim(1) model.sdim(2) (model.zEnd-model.zStart+1)]).*model.DataScale;
denom = SINalphat1.*(1-E1);

PD =  enum./denom;

clear enum denom SINalphat1;






%
% calculate_R1_PD <- function(theta1,
%                                mt,
%                                verbose = TRUE) {
%   
%   b1Map <- 1
%   ## get correct flip angles and TR times
%   t1FA <- 27
%   pdFA <- 13
%   t1TR <- 22.5
%   
%   ## calculate E1
%   if (verbose) cat("calculating R1 ... ")
%   alphat1 <- b1Map * t1FA / 180 * pi
%   alphapd <- b1Map * pdFA / 180 * pi
%   SINalphat1 <- sin(alphat1)
%   COSalphat1 <- cos(alphat1)
%   SINalphapd <- sin(alphapd)
%   COSalphapd <- cos(alphapd)
%   rm(alphat1, alphapd)
%   
%   enum <- theta1 - SINalphat1/SINalphapd * mt
%   denom <- theta1 * COSalphat1 - SINalphat1/SINalphapd * mt * COSalphapd
%   
%   E1 <- enum/denom
%   rm(enum, denom, COSalphapd, SINalphapd)
%   R1 <- -log(E1)/t1TR
%   if (verbose) cat("done\n")
%   
%   ## calculate PD
%   if (verbose) cat("calculating PD ... ")
%   enum <- (1 - COSalphat1 * E1) * theta1 * 1000
%   denom <- SINalphat1 * (1 - E1)
%   PD <- enum/denom
%   rm(enum, denom, SINalphat1)
%   if (verbose) cat("done\n")
%   
%   
%   
%   
%   invisible(list(R1 = R1,
%                  PD = PD
%                  ))
% }