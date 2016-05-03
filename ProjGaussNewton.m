function [uc,his] = ProjGaussNewton(fctn,uc,varargin)

if nargin==0
    help(mfilename);
    runMinimalExample;
    return;
end

maxIter      = 30;
uStop        = [];
Jstop        = [];
paraStop     = [];
tolJ         = 1e-5;            % for stopping, objective function
tolU         = 1e-4;            %   - " -     , current value
tolG         = 1e-4;            %   - " -     , norm of gradient
LSMaxIter    = 15;              % maximum number of line search iterations
LSreduction  = 0;               % minimal reduction in line search
Plots        = @(task,para) [];
verbose      = false;
upper        = Inf*ones(numel(uc),1);
lower        = -Inf*ones(numel(uc),1);
maxStep      = 1000;
for k=1:2:length(varargin),     % overwrites default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;
if verbose
    fprintf('%s %s %s\n',ones(1,30)*char('='),mfilename,ones(1,40)*char('='));
end

% evaluate objective function for stopping criteria

if isempty(uStop), uStop = uc; end
if isempty(Jstop) || isempty(paraStop),
    [Jstop,paraStop] = fctn(uStop,[]); Jstop = abs(Jstop) + (Jstop == 0);
    Jstop = paraStop.Dcs;
    Plots('stop',paraStop);
end;
nvoxel = paraStop.nvoxel;
uc    = reshape(uc,[],nvoxel);
lower = reshape(lower,[],nvoxel);
upper = reshape(upper,[],nvoxel);

if verbose
fprintf('[ maxIter=%s / tolJ=%s / tolU=%s / tolG=%s / length(yc)=%d ]\n',...
  num2str(maxIter),num2str(tolJ),num2str(tolU),num2str(tolG),length(uc));
fprintf('%4s %9s %9s %9s %9s %9s %2s %9s %9s\n','iter','Jc','Jold-Jc','norm(dJ)','Dc','Rc','LS','Active','perc.conv');
fprintf('%4d %9.2e %9.2e %9.2e %9.2e %9.2e\n',-1,sum(Jstop),0,0,paraStop.Dc,paraStop.Rc)
end
% evaluate objective function for starting values and plots
[Jc,para,dJ,H] = fctn(uc,[]); 
Plots('start',para);
iter = 0; uOld = 0*uc; Jold = para.Dcs; u0 = uc; LSiter = 0;

% initialize
nvoxel    = para.nvoxel;
STOP      = zeros(5,nvoxel);
credit    = true(1,nvoxel);
Active    = (uc<=lower)|(uc>=upper);
normU     = @(u) sqrt(sum(u).^2);

Dcs    = para.Dcs;
while 1
    Plots(iter,para);
    
    
    % check stopping rules
    STOP(1,credit) = (iter>0) & (abs(Jold(credit)-para.Dcs)   <= tolJ*(1+abs(Jstop(credit))));
    STOP(2,:)      = (iter>0) & (normU(uc-uOld)               <= tolU*(1+normU(u0)));
    STOP(3,credit) = (normU(para.dDs)                         <= tolG*(1+abs(Jstop(credit))));
    STOP(4,credit) = (norm(para.dDs)                          <= 1e6*eps);
    STOP(5,:) = (iter >= maxIter);
    
    credit = ( ~STOP(1,:) | ~STOP(2,:) | ~STOP(3,:)) & ~STOP(4,:);
    
    % update
    Active =  (uc(:,credit)<=lower(:,credit)) | (uc(:,credit)>=upper(:,credit));
    [Jc,para,dJ,H] = fctn(uc,credit); % evalute objective function
    Dcs(credit) = para.Dcs;
    
    % some output
    if verbose
        nActive = nnz((uc<=lower) | (uc>=upper));
        fprintf('%4d %9.2e %9.2e %9.2e %9.2e %9.2e %2d\t%2d\t%1.1f%%\n',...
        iter,.5*sum(Dcs),sum(.5*(Jold-Dcs)),norm(dJ),para.Dc,para.Rc,LSiter,nActive,100*nnz(not(credit))/nvoxel);
    end
    
    if isempty(credit) || iter >=maxIter, break;  end;
    
    iter = iter + 1;
    % solve the Gauss-Newton System on inactive set
    Hr  = H(not(Active),not(Active));
    Jr  = dJ(not(Active));
    dur = -(Hr\Jr);
    if max(abs(dur))>maxStep
        dur = maxStep*dur/max(abs(dur));
    end
    du = 0*uc(:,credit);
    du(not(Active)) = dur;
        
    % take gradient step on inactive set
    gl = max(0,-dJ(uc(:,credit)==lower(:,credit)));
    if not(isempty(gl)),
        if max(abs(gl))>max(abs(dur)),
            gl = max(abs(dur))*gl./max(abs(gl));
        end
        gl = max(gl,0);
        du(uc(:,credit)==lower(:,credit))=gl;    
    end
    
    gu = min(0,-dJ(uc(:,credit)==upper(:,credit)));
    if not(isempty(gu)),
        if max(abs(gu))>max(abs(dur)),
            gu = max(abs(dur))*gu./max(abs(gu));
        end
        gu = min(gu,0);
        du(uc(:,credit)==upper(:,credit)) = gu;
    end
    
    if all(normU(du)<=1e-10),
        %fprintf('\n %s - found stationary point.\n\n',mfilename)
        break;
    end
        
    
    
    % perform Armijo linesearch
    t = 1;
    LSiter = 1; creditLS = credit; idx = 1:nnz(creditLS);
    while 1
        ut(:,creditLS) = uc(:,creditLS) + t*du(:,idx);
        ut(ut<lower) = lower(ut<lower);
        ut(ut>upper) = upper(ut>upper);
        
        [~,parat] = fctn(ut,creditLS);
        % Armijo condition: fctn(uc+t*du) <= Jc + t*LSreduction*(dJ*du)
        LS = parat.Dcs < (para.Dcs(idx)); % compare
        if all(LS) || LSiter>=LSMaxIter, break; end
        creditLS(creditLS) = ~LS;
        idx(LS) = [];
        t = t/2;
        LSiter = LSiter+1;
    end  
    if not(LS),  warning('Linesearch failed!'); break; end
    
    uOld = uc; Jold(credit) = para.Dcs; uc = ut;
end
% report history
if verbose
fprintf('%s\nSTOPPING:\n',ones(1,87)*char('-'));
fprintf('%1.1f%%\t[ %-10s=%16.8e <= %-25s=%16.8e]\n',100*nnz(STOP(1,:))/nvoxel,...
  '(Jold-Jc)',max(Jold-Jc),'tolJ*(1+|Jstop|)'   ,tolJ*max(1+abs(Jstop)));
fprintf('%1.1f%%\t[ %-10s=%16.8e <= %-25s=%16.8e]\n',100*nnz(STOP(2,:))/nvoxel,...
  '|yc-yOld|',norm(uc-uOld),'tolU*(1+norm(yc)) ',tolU*(1+norm(uc)));
fprintf('%1.1f%%\t[ %-10s=%16.8e <= %-25s=%16.8e]\n',100*nnz(STOP(3,:))/nvoxel,...
  '|dJ|',norm(dJ),'tolG*(1+abs(Jstop))',min(tolG*(1+abs(Jstop))));
fprintf('%1.1f%%\t[ %-10s=%16.8e <= %-25s=%16.8e]\n',100*nnz(STOP(4,:))/nvoxel,...
  'norm(dJ)',norm(dJ),'eps',1e3*eps);
fprintf('%d\t[ %-10s=  %-14d >= %-25s=  %-14d]\n',any(STOP(5,:)),...
  'iter',iter,'maxIter',maxIter);

fprintf('%s %s : done ! %s\n',ones(1,30)*char('='),mfilename,ones(1,30)*char('='));
end


function [Jc,para,dJ,d2J] = testFun(A,b,x)
    res = A*x-b;
    Jc  = 0.5*res'*res;
    dJ  = A'*res;
    d2J = A'*A;
    para = struct('m',0,'omega',0,'Dc',Jc,'Dcs',Jc,'dDs',dJ,'Rc',0,'nvoxel',1);


function runMinimalExample
    A = speye(4);
    b = - rand(4,1);
    u = ones(4,1);
    l = zeros(4,1);
    
    fctn = @(x,credit) testFun(A,b,x);
    
    x0 = [0;0; 1;1];
    xopt = ProjGaussNewton(fctn,x0,'upper',u,'lower',l,'verbose',true);
    
    xopt
    
    


