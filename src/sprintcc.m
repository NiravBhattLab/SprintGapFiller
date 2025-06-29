function [ConsReacIDS,LPS] = sprintcc(model,tol,consistencyType)
% USAGE:
% [ConsReacIDS,LPS] = sprintcc(model,tol,consistencyType)
%
% INPUTS:
%    model: COBRA model structure.
% 
% OPTIONAL INPUTS:
%    tol:             Tolerance level (minimum absolute flux that has to be carried
%                     by a reaction for it to be defined as consistent).
%                     (Default: 1e-4)
%    consistencyType: Type of consistency check to make. 'stoichiometry'
%                     assumes steadystate condition and 'topology' assumes accumulation
%                     condition. (Default: 'stoichiometry')
%
% OUTPUTS:
%    ConsReacIDS: Reaction IDs corresponding to reactions in the input
%                 model that are consistent/unblocked
%    LPS:         Number of LPs solved to get the consistent reaction IDs
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

if ~exist('tol', 'var') || isempty(tol)
    tol =1e-4;    
end
if ~exist('consistencyType', 'var') || isempty(consistencyType)
    consistencyType = 'stoichiometry';    
end

[~,n] = size(model.S);
IrR = model.ub<=0; % reactions irreversible in reverse direction
temp = model.ub(IrR);
model.S(:,IrR) = -model.S(:,IrR);
model.ub(IrR) = -1*model.lb(IrR);
model.lb(IrR) = -1*temp;
rev = ones(n,1);
rev(model.lb>=0) = 0;
model.rev=rev;
LPS=0;
prev_rxns = false(numel(model.rxns),1);
temp_core = true(numel(model.rxns),1);

if strcmp(consistencyType,'stoichiometry')
    steadystate = 1;
elseif strcmp(consistencyType,'topology')
    steadystate = 0;
end

while sum(temp_core)~=sum(prev_rxns)
    prev_rxns = temp_core;
    LPS = LPS+2;
    [flux1,~] = forwardcc(model,temp_core,tol,steadystate);
    temp_core(abs(flux1)>=tol*0.99)=false;
    [flux2,~] = reverse(model,temp_core,tol,steadystate);
    temp_core(abs(flux2)>=tol*0.99)=false;
end

ConsReacIDS=setdiff(1:numel(model.rxns),find(temp_core))';