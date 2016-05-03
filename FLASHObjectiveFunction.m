  
function [Dc,para,dD,H,Hstar] = FLASHObjectiveFunction(model,data,TE,indicator,credit)
%
% Input
% model - [T1(0); MT(0); PD(0); R2*] column vector
% data  - [T1; MT; PD], column vector
% TE - echo times, size(TE,1) = size(data,1)
% TR - readout times
% FA - flip angle
% indicator - values are {1,2,3} for T1,MT,PD images
% credit - specifies indices that need to be worked on

if nargin<=4 || isempty(credit)
    credit = 1:size(data,2);
end

ndata = numel(TE);
model = reshape(model,4,[]);
model = model(:,credit);
nvoxel = size(model,2);
data  = reshape(data,ndata,[]);
data  = data(:,credit);

doDerivative = (nargout>2);

% evaluate forward model
pred = model(indicator,:) .* exp(-TE*model(4,:));
res  = pred - data;
Dcs  = sum(res.^2,1); % misfit for each voxel
Dc   = 0.5*sum(Dcs);

para = struct('omega',[],'m',[],'nvoxel',nvoxel,'Dc',Dc,'Dcs',Dcs,'Rc',0.0,'res',res);
if not(doDerivative)
    dD = [];
    H  = [];
else
    I = 1:size(data,1);
    J = indicator;
    A1 = sparse(I,J,ones(size(data,1),1), ndata,4);
   
    I  = 1:(nvoxel*ndata);
    J  = kron(1:nvoxel, 4*ones(ndata,1));
    V = -sdiag(TE)*model(indicator,:);
    A2 = sparse(I(:),J(:),V(:),nvoxel*ndata,4*nvoxel); 
    
    A = kron(speye(size(model,2)),A1) + A2;
    
    J = sdiag(exp(-TE*model(4,:)))*A;
    
    dD = J'*res(:);
    para.dDs = reshape(dD,4,[]);
    Hstar = J'*J;
    H  = Hstar + 1e-4*speye(numel(model));
end
    
end
