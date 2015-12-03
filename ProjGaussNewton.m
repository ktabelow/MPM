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
solver       = [];
verbose      = false;
upper        = Inf*ones(numel(uc),1);
lower        = -Inf*ones(numel(uc),1);
maxStep      = 1000;
for k=1:2:length(varargin),     % overwrites default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;
if verbose
    fprintf('%s %s %s\n',ones(1,20)*char('='),mfilename,ones(1,20)*char('='));
end

% evaluate objective function for stopping criteria
if isempty(uStop), uStop = uc; end
if isempty(Jstop) || isempty(paraStop),
    [Jstop,paraStop] = fctn(uStop); Jstop = abs(Jstop) + (Jstop == 0);
    Plots('stop',paraStop);
end;
if verbose
fprintf('[ maxIter=%s / tolJ=%s / tolU=%s / tolG=%s / length(yc)=%d ]\n',...
  num2str(maxIter),num2str(tolJ),num2str(tolU),num2str(tolG),length(uc));
fprintf('%4s %9s %9s %9s %9s %9s %2s %9s\n','iter','Jc','Jold-Jc','norm(dJ)','Dc','Rc','LS','Active');
fprintf('%4d %9.2e %9.2e %9.2e %9.2e %9.2e\n',-1,Jstop,0,0,paraStop.Dc,paraStop.Rc)
end
% evaluate objective function for starting values and plots
[Jc,para,dJ,H] = fctn(uc); 
Plots('start',para);
iter = 0; uOld = 0*uc; Jold = Jc; u0 = uc; LSiter = 0;

% initialize
STOP = zeros(5,1);

Active = (uc<=lower)|(uc>=upper);
while 1
    Plots(iter,para);
    % some output
    if verbose
        fprintf('%4d %9.2e %9.2e %9.2e %9.2e %9.2e %2d\t%2d\n',...
        iter,Jc,Jold-Jc,norm(dJ),para.Dc,para.Rc,LSiter,nnz(Active));
    end
    % check stopping rules
    STOP(1) = (iter>0) && abs(Jold-Jc)   <= tolJ*(1+abs(Jstop));
    STOP(2) = (iter>0) && (norm(uc-uOld) <= tolU*(1+norm(u0)));
    STOP(3) = norm(dJ)                   <= tolG*(1+abs(Jstop));
    STOP(4) = norm(dJ)                   <= 1e6*eps;
    STOP(5) = (iter >= maxIter);
    
    if all(STOP(1:3)) || any(STOP(4:5)), break;  end;
    
    iter = iter + 1;
    % solve the Gauss-Newton System on inactive set
    Hr = H(not(Active),not(Active));
    Jr = dJ(not(Active));
    dur = -(Hr\Jr);
    if max(abs(dur))>maxStep
        dur = maxStep*dur/max(abs(dur));
    end
    du = 0*uc;
    du(not(Active)) = dur;
        
    % take gradient step on inactive set
    gl = -dJ(uc==lower);
    if not(isempty(gl)),
        if max(abs(gl))>max(abs(dur)),
            gl = max(abs(dur))*gl./max(abs(gl));
        end
        gl = max(gl,0);
        du(uc==lower)=gl;    
    end
    
    gu = -dJ(uc==upper);
    if not(isempty(gu)),
        if max(abs(gu))>max(abs(dur)),
            gu = max(abs(dur))*gu./max(abs(gu));
        end
        gu = min(gu,0);
        du(uc==upper) = gu;
    end
    
    if norm(du)<=1e-10,
        fprintf('\n %s - found stationary point.\n\n',mfilename)
        break;
    end
        
    
    
    % perform Armijo linesearch
    t = 1;
    LSiter = 1;
    while 1
        ut = uc + t*du;
        ut(ut<lower) = lower(ut<lower);
        ut(ut>upper) = lower(ut>upper);
        
        Jt = fctn(ut);
        % Armijo condition: fctn(uc+t*du) <= Jc + t*LSreduction*(dJ*du)
        LS = (Jt<Jc + t*LSreduction*(dJ'*du)); % compare
        if LS || LSiter>=LSMaxIter, break; end
        t = t/2;
        LSiter = LSiter+1;
    end  
    if not(LS),  warning('Linesearch failed!'); break; end
    
    % update 
    uOld = uc; Jold = Jc; uc = ut;
    Active = (uc<=lower) | (uc>=upper);
    [Jc,para,dJ,H] = fctn(uc); % evalute objective function
end
% report history
if verbose
fprintf('%s\nSTOPPING:\n',ones(1,70)*char('-'));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(1),...
  '(Jold-Jc)',(Jold-Jc),'tolJ*(1+|Jstop|)'   ,tolJ*(1+abs(Jstop)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(2),...
  '|yc-yOld|',norm(uc-uOld),'tolU*(1+norm(yc)) ',tolU*(1+norm(uc)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(3),...
  '|dJ|',norm(dJ),'tolG*(1+abs(Jstop))',tolG*(1+abs(Jstop)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(4),...
  'norm(dJ)',norm(dJ),'eps',1e3*eps);
fprintf('%d[ %-10s=  %-14d >= %-25s=  %-14d]\n',STOP(5),...
  'iter',iter,'maxIter',maxIter);

fprintf('%s %s : done ! %s\n',ones(1,20)*char('='),mfilename,ones(1,20)*char('='));
end


function [Jc,para,dJ,d2J] = testFun(A,b,x)
    res = A*x-b;
    Jc  = 0.5*res'*res;
    dJ  = A'*res;
    d2J = A'*A;
    para = struct('m',0,'omega',0,'Dc',Jc,'Rc',0);


function runMinimalExample
    A = speye(4);
    b = - rand(4,1);
    u = ones(4,1);
    l = zeros(4,1);
    
    fctn = @(x) testFun(A,b,x);
    
    x0 = [0;0; 1;1];
    xopt = ProjGaussNewton(fctn,x0,'upper',u,'lower',l,'verbose',true);
    
    xopt,
    
    
    


