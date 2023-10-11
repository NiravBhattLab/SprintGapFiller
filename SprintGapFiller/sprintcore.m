function [ConsModel,LPS] = sprintcore(model,core,tol)
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


[~,n] = size(model.S);
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
    [flux1,~] = forward(model,temp_core,tol);
    if sum(abs(flux))==0
        flux = flux1;
    else
        c1=unifrnd(0.45,0.55,1);
        flux = (c1*flux)+((1-c1)*flux1);
    end
    temp_core(core==1 & abs(flux)>=0.99*tol) = 0; 
    if ~any(temp_core)
        break
    end
    LPS = LPS+1;
    [flux2,~] = reverse(model,temp_core,tol);
    c1=unifrnd(0.45,0.55,1);
    flux = (c1*flux)+((1-c1)*flux2);
    temp_core(core==1 & abs(flux)>=0.99*tol) = 0; 
end

direction = zeros(n,1);
direction(core==1&flux>0) = 1;
direction(core==1&flux<0) = -1;
LPS = LPS+1;
reacInd = findConsistentReacID(model,direction,tol); %LPminimal

ConsModel = removeRxns(model, setdiff(model.rxns,model.rxns(reacInd)));
% ConsModel = removeUnusedGenes(ConsModel);