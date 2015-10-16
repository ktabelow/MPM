function [uc,his] = GaussNewton(fctn,uc,varargin)

if nargin==0
    help(mfilename);
    runMinimalExample;
    return;
end

maxIter      = 10;
uStop        = [];
Jstop        = [];
paraStop     = [];
tolJ         = 1e-3;            % for stopping, objective function
tolU         = 1e-2;            %   - " -     , current value
tolG         = 1e-2;            %   - " -     , norm of gradient
LSMaxIter    = 10;              % maximum number of line search iterations
LSreduction  = 0;               % minimal reduction in line search
Plots        = @(task,para) [];
solver       = [];
verbose      = false;
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
fprintf('%4s %9s %9s %9s %9s %9s %2s\n','iter','Jc','Jold-Jc','norm(dJ)','Dc','Rc','LS');
fprintf('%4d %9.2e %9.2e %9.2e %9.2e %9.2e\n',-1,Jstop,0,0,paraStop.Dc,paraStop.Rc)
end
% evaluate objective function for starting values and plots
[Jc,para,dJ,H] = fctn(uc); 
Plots('start',para);
iter = 0; uOld = 0*uc; Jold = Jc; u0 = uc; LSiter = 0;

% initialize
STOP = zeros(5,1);

while 1
    Plots(iter,para);
    % some output
    if verbose
        fprintf('%4d %9.2e %9.2e %9.2e %9.2e %9.2e %2d\n',...
        iter,Jc,Jold-Jc,norm(dJ),para.Dc,para.Rc,LSiter);
    end
    % check stopping rules
    STOP(1) = (iter>0) && abs(Jold-Jc)   <= tolJ*(1+abs(Jstop));
    STOP(2) = (iter>0) && (norm(uc-uOld) <= tolU*(1+norm(u0)));
    STOP(3) = norm(dJ)                   <= tolG*(1+abs(Jstop));
    STOP(4) = norm(dJ)                   <= 1e6*eps;
    STOP(5) = (iter >= maxIter);
    
    if all(STOP(1:3)) || any(STOP(4:5)), break;  end;
    
    iter = iter + 1;
    % solve the Gauss-Newton System
    du =  solveGN(H,-dJ,solver);
    
    % perform Armijo linesearch
    t = 1;
    LSiter = 1;
    while 1
        ut = uc + t*du;
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
function du = solveGN(H,rhs,solver)
tolCG     = .1;
maxIterCG = 100;
if isnumeric(H),
    if isempty(solver),
        du = H\rhs;
    else
        switch solver
            case 'pcg'
                L   = tril(H); % Symmetric Gauss Seidel Preconditioning,
                D   = diag(H); % L is lower, D is diagonal, U = L'
                SGS = @(x) L\(D.*(L'\x));
                [du,flag,relres,iter] = pcg(H,rhs,tolCG,maxIterCG,SGS);
            case 'jacobi-pcg'
                D   = diag(H); % D is diagonal
                PC = @(x) D.\x;
                [du,flag,relres,iter] = pcg(H,rhs,tolCG,maxIterCG,PC);
            case 'cg'
                [du,flag,relres,iter] = pcg(H,rhs,tolCG,maxIterCG);
                
        end
    end
elseif isa(H,'function_handle'),
    du = pcg(H,rhs,tolCG,maxIterCG);
elseif isa(H,'struct'),
    if isempty(solver), solver = H.solver; end;
    switch solver
        case 'pcg'
            C = @(x) H.d2R\x;
            A = @(x) H.d2D(x) + H.d2R * x;
            
            du = pcg(A,rhs,tolCG,maxIterCG,C);
            
        case 'cg'
            A = @(x) H.d2D(x) + H.d2R * x;
            du = pcg(A,rhs,tolCG,maxIterCG);
        case 'jacobi-cg'
            C = @(x) full(diag(H.d2R)).\x;
            
            A = @(x) H.d2D(x) + H.d2R * x;
            du = pcg(A,rhs,tolCG,maxIterCG,C);
            
        otherwise error('%s- solver %s nyi!',mfilename,solver);
    end
end


function runMinimalExample
domain = [0 3 0 5];
m      = [33 33];
h      = (domain(2:2:end)-domain(1:2:end))./m;

