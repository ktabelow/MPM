function [Dc,para,dD,H] = FLASHobjFctnWithoutMT(model,data,TE,indicator)
%
% Input
% model - [T1(0); PD(0); R2*] column vector
% data  - [T1; PD], column vector
% TE - echo times, size(TE,1) = size(data,1)
% TR - readout times
% FA - flip angle
% indicator - values are {1,2} for T1,PD images

ndata = numel(TE);
model = reshape(model,3,[]);
nvoxel = size(model,2);
data  = reshape(data,ndata,[]);


doDerivative = (nargout>2);

% evaluate forward model
pred = model(indicator,:) .* exp(-TE*model(3,:));
res  = pred(:) - data(:);
Dcs  = sum(res.^2,1); % misfit for each voxel
Dc   = 0.5*sum(Dcs);

para = struct('omega',[],'m',[],'Dc',Dc,'Dcs',Dcs,'Rc',0.0,'res',res);
if not(doDerivative)
    dD = [];
    H  = [];
else
    I = 1:size(data,1);
    J = indicator;
    A1 = sparse(I,J,ones(size(data,1),1), ndata,3);
   
    I  = 1:(nvoxel*ndata);
    J  = kron(1:nvoxel, 3*ones(ndata,1));
    V = -sdiag(TE)*model(indicator,:);
    A2 = sparse(I(:),J(:),V(:),nvoxel*ndata,3*nvoxel); 
    
    A = kron(speye(nvoxel),A1) + A2;
    
    J = sdiag(exp(-TE*model(3,:)))*A;
    
    dD = J'*res;
    para.dDs = sqrt(sum(reshape(dD,3,[]).^2,1));
    H  = J'*J + 1e-4*speye(numel(model));
end
    


    
 
    
    
    
end
