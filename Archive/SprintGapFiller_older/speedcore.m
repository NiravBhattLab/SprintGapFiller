function [ConsModel,stat,LPS] = speedcore(model,core,tol)
%%INPUT
%       model: COBRA model structure. The model has to be consistent model
%
%       core: Indices of reactions that have to be present in the final
%             model
%
%       tol: tolerance level (minimum absolute flux that has to be carried
%            by all the reactions in the model)

%%OUTPUT
%       ConsModel: The consistent model with no blocked reactions and has
%                  all the core reactions in them
%       
%       LPS: Number of LPs used to get the ConsModel

%%AUTHOR
%       Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

if nargin < 3 || isempty(tol)
    tol=1e-4; 
end

SF=1e3;
[m,n] = size(model.S);
core = ismember([1:n],core)';
IrR = model.ub<=0; % reactions irreversible in reverse direction
temp = model.ub(IrR);
model.S(:,IrR) = -model.S(:,IrR);
model.ub(IrR) = -1*model.lb(IrR);
model.lb(IrR) = -1*temp;
rev = ones(n,1);
rev(model.lb>=0) = 0;
model.rev=rev;

temp_core = core;
flux = zeros(n,1);
LPS=0;
while any(temp_core)
    LPS = LPS+1;
    [flux1,~,f_stat] = forward(model,temp_core,tol);
    if sum(abs(flux))==0
        flux = flux1;
    else
        c1=unifrnd(0.45,0.55,1);
        flux = (c1*flux)+((1-c1)*flux1);
    end
    flux=flux*SF;
    temp_core(core==1 & abs(flux)>=tol) = 0; 
    if ~any(temp_core)
        break
    end
    LPS = LPS+1;
    [flux2,~,r_stat] = reverse(model,temp_core,tol);
    c1=unifrnd(0.45,0.55,1);
    flux = (c1*flux)+((1-c1)*flux2);
    flux=flux*SF;
    temp_core(core==1 & abs(flux)>=tol) = 0; 
    % If no optimal solution can be found
    if f_stat~=1 || r_stat~=1
        stat=0;
        break
    else
        stat=1;
    end
end

direction = zeros(n,1);
direction(core==1&flux>0) = 1;
direction(core==1&flux<0) = -1;
LPS = LPS+1;
if stat
    reacInd = findConsistentReacID(model,direction,tol); %LPminimal
    ConsModel = removeRxns(model, setdiff(model.rxns,model.rxns(reacInd)));
else
    ConsModel = struct();
end
% ConsModel = removeUnusedGenes(ConsModel);