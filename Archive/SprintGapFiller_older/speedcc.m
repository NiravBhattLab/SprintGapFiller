function [ConsReacIDS,feas,LPS] = speedcc(model,tol,feas)

if nargin < 2 || isempty(tol)
    tol=1e-4; 
end
if nargin < 3 || isempty(feas)
    feas = 0.99;
end

[m,n] = size(model.S);
IrR = model.ub<=0; % reactions irreversible in reverse direction
temp = model.ub(IrR);
SF=1;
model.S(:,IrR) = -model.S(:,IrR);
model.ub(IrR) = -1*model.lb(IrR);
model.lb(IrR) = -1*temp;
rev = ones(n,1);
rev(model.lb>=0) = 0;
model.rev=rev;
LPS=0;
model.lb=model.lb*SF; model.ub=model.ub*SF;
prev_rxns = false(numel(model.rxns),1);
temp_core = true(numel(model.rxns),1);
while sum(temp_core)~=sum(prev_rxns)
    prev_rxns = temp_core;
    LPS = LPS+2;
    [flux1,~] = forwardcc(model,temp_core,tol*SF);
    temp_core(abs(flux1)>=tol*SF*feas)=false;
    [flux2,~] = reverse(model,temp_core,tol*SF);
    temp_core(abs(flux2)>=tol*SF*feas)=false;
end

ConsReacIDS=setdiff([1:numel(model.rxns)],find(temp_core))';