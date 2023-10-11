function reacInd = findConsistentReacID(model,direction,tol)
SF=1;
[m,n] = size(model.S);
dir0 = direction==0;
n_ = sum(dir0);

% objective
f = [zeros(n,1);1*ones(n_,1)];

% equalities
Aeq = [model.S, sparse(m,n_)];
beq = zeros(m,1);
csenseeq = repmat('E',m,1); % equality

% inequalities
temp1 = speye(n);
temp2 = speye(n_);
Aineq1 = [temp1(dir0,:),temp2];
bineq1 = zeros(n_,1);
csenseineq1 = repmat('G',n_,1); % greater than

Aineq2 = [temp1(dir0,:),-1*temp2];
bineq2 = zeros(n_,1);
csenseineq2 = repmat('L',n_,1); % lesser than

% bounds
lb = model.lb;
lb(direction==1)=max([tol*ones(sum(direction==1),1),lb(direction==1)],[],2);
lb = [lb;zeros(n_,1)]*SF;
ub = model.ub;
ub(direction==-1)=-tol*ones(sum(direction==-1),1);
ub = [ub;Inf(n_,1)]*SF;

% Set up LP problem
LPproblem.A=[Aeq;Aineq1;Aineq2];
LPproblem.b=[beq;bineq1;bineq2];
LPproblem.lb=lb;
LPproblem.ub=ub;
LPproblem.c=f;
LPproblem.osense=1;%minimise
LPproblem.csense = [csenseeq; csenseineq1; csenseineq2];
solution = solveCobraLP(LPproblem);

if solution.stat~=1
    fprintf('%s%s\n',num2str(solution.stat),' = solution.stat')
    fprintf('%s%s\n',num2str(solution.origStat),' = solution.origStat')
    warning('LP solution may not be optimal')
end

x=solution.full;
reacInd = abs(x(1:n))>0;