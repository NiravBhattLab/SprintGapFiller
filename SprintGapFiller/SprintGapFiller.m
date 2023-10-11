function [MicComModel,BlockedCoreRxns] = SprintGapFiller(model,core,tol)
%%INPUT
%       model: Universal inconsistent model with blocked reaction
%
%       core: core reactions which has to be present in the final model
%
%       tol: Minimum absolute flux required for a reaction to be unblocked

%%OUTPUT
%       MicComModel: The consistent community model
%       
%       BlockedCoreRxns: Core reactions that cannot have absolute flux 
%                        above the tol value

%%AUTHOR
%       Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras


if nargin < 2 || isempty(tol)
    tol=1e-4; 
end

[~,n] = size(model.S); % number of metabolites and reactions
IrR = model.ub<=0; % reactions irreversible in reverse direction
temp = model.ub(IrR);
model.S(:,IrR) = -model.S(:,IrR);
model.ub(IrR) = -1*model.lb(IrR);
model.lb(IrR) = -1*temp;
rev = ones(n,1);
rev(model.lb>=0) = 0;
model.rev=rev;
prev_rxns = false(n,1);
temp_core = true(n,1);
flux = zeros(n,1); % initiating the flux vector that will carry directionality info

while sum(temp_core)~=sum(prev_rxns)
    prev_rxns = temp_core;
    % maximizing number of reactions with forward flux
    [flux1,~] = forwardcc(model,temp_core,tol);
    temp_core(abs(flux1)>=tol*0.99)=false;
    if sum(abs(flux))==0
        flux = flux1;
    else
        c1=round(unifrnd(0.45,0.55,1),4);
        flux = (c1*flux)+((1-c1)*flux1);
    end
    % maximizing number of reactions with reverse flux
    [flux2,~] = reverse(model,temp_core,tol);
    temp_core(abs(flux2)>=tol*0.99)=false;
    
    c1=round(unifrnd(0.45,0.55,1),4);
    flux = (c1*flux)+((1-c1)*flux2);
    
end
BlckdRxns = find(temp_core);
BlockedCoreRxns = model.rxns(intersect(find(core),BlckdRxns));
core(temp_core) = 0; % not forcing any flux through the blocked reactions
direction = zeros(n,1);
direction(core==1&flux>0) = 1;
direction(core==1&flux<0) = -1;
if any(core==1&flux==0)
    warning('Any of the core reactions carry zero flux')
end
reacInd= findConsistentReacID(model,direction,tol); %LPminimal
MicComModel = removeRxns(model, setdiff(model.rxns,model.rxns(reacInd)));

end