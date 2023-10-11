function [flux,z] = forward(model,core,tol)
% model must have fields 'S' and 'rev'
% core is a binary vector indicating whether core or non-core reaction
[m,n] = size(model.S);
rev_core = model.rev & core; % reversible core reactions
n_rev_core = sum(rev_core);% number of reversible core reactions

% objective
f = -1*[zeros(n,1);unifrnd(1,1.1,n_rev_core,1)];

% equalities
Aeq = [model.S, sparse(m,n_rev_core)];
beq = zeros(m,1);
csenseeq = repmat('E',m,1); % equality

% inequalities
temp1 = speye(n);
temp2 = speye(n_rev_core);
Aineq = [temp1(rev_core,:),-1*temp2 ];
bineq = zeros(n_rev_core,1);
csenseineq = repmat('G',n_rev_core,1); % greater than

% bounds
irr_core=model.rev==0 & core==1;
lb = model.lb;
lb(irr_core)=max([tol*ones(sum(irr_core),1),lb(irr_core)],[],2);% irreversible core reactions
ub = model.ub;
lb = [lb;-Inf(n_rev_core,1)];
ub = [ub; tol*ones(n_rev_core,1)];

% Set up LP problem
LPproblem.A=[Aeq;Aineq];
LPproblem.b=[beq;bineq];
LPproblem.lb=lb;
LPproblem.ub=ub;
LPproblem.c=f;
LPproblem.osense=1;%minimise
LPproblem.csense = [csenseeq; csenseineq];

solution = solveCobraLP(LPproblem);
if solution.stat~=1
    fprintf('%s%s\n',num2str(solution.stat),' = solution.stat')
    fprintf('%s%s\n',num2str(solution.origStat),' = solution.origStat')
    warning('LP solution may not be optimal')
end

x=solution.full;

if ~isempty(x)
    flux = x(1:n);
    z = x(n+1:end);
else
    flux=nan(n,1);
    z=nan(n_rev_core,1);
end

