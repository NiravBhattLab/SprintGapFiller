function [ConsModel,LPS] = sprintcore(model,coreRxns,tol,gapFilltype,weights,nSol,altSolMethod,probType,solveTime,remGene,prevSols)
% USAGE:
%   [ConsModel,LPS] = sprintcore(model,coreRxns,tol,gapFilltype,weights,nSol,altSolMethod,probType,solveTime,remGene,prevSols)
%
% INPUTS:
%   model:    COBRA model structure. The model has to be a consistent model
%   coreRxns: Indices of reactions that have to be present in the final
%             model
%
% OPTIONAL INPUTS:
%   tol:          Tolerance level (minimum absolute flux that has to be carried
%                 by all the reactions in the model) (Default: 1e-4)
%   gapFilltype:  Type of gapfilling to apply. Either 'topology' or
%                 'stoichiometry' (stoichiometric matrix) based. (Default:'stoichiometry')
%   weights:      Weights for non-core reactions. More the weights, lesser
%                 the chance to get included in the final model (Default: ones)
%   nSol:         Number of alternative solutions required (Default: 1)
%   altSolMethod: Method to find the alternate solutions.
%                 accepted values: 'coreDirection', 'pathwayExclusion'.
%                 Note: 'pathwayExclusion' works only for MILP probType
%                 (Default: 'coreDirection')
%   probType:     Which method to use to find the minimal reaction set.
%                 accepted values: 'LP','MILP','DC'. (Default: 'LP')
%   solveTime:    maximum runtime for solving MILP problem, if probType=='MILP'
%                 (Default: 7200s)
%   remGene:      Bool value indicating whether to remove the unused genes
%                 or not (Default: 0 (doesn't remove the unused genes))
%   prevSols:     A cell of previosuly obtained solutions that should not
%                 be a part of any new solutions
%
% OUTPUTS:
%   ConsModel: The consistent model with no blocked reactions and has
%              all the core reactions in them (If nSol==1). A cell 
%              consisting of consistent models in them (If nSol > 1)
%   LPS:       Number of LPs used to get ConsModel, if probType is 'LP'
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

if ~exist('tol', 'var') || isempty(tol)
    tol=1e-4;     
end
if ~exist('weights', 'var') || isempty(weights)
    weights = ones(numel(model.rxns),1);  
end
if ~exist('nSol', 'var') || isempty(nSol)
    nSol=1;  
end
if ~exist('probType', 'var') || isempty(probType)
    probType='LP';  
end
if ~exist('remGene', 'var') || isempty(remGene)
    remGene=0;  
end
if ~exist('prevSols', 'var') || isempty(remGene)
    prevSols={};  
end
if ~exist('gapFilltype', 'var') || isempty(gapFilltype)
    gapFilltype='stoichiometry';  
end

[~,n] = size(model.S);
coreRxns = ismember(1:n,coreRxns)';
IrR = model.ub<=0; % reactions irreversible in reverse direction
temp = model.ub(IrR);
model.S(:,IrR) = -model.S(:,IrR);
model.ub(IrR) = -1*model.lb(IrR);
model.lb(IrR) = -1*temp;
rev = ones(n,1);
rev(model.lb>=0) = 0;
model.rev=rev;

temp_core = coreRxns;
flux = zeros(n,1);
LPS=0;

if strcmp(gapFilltype,'stoichiometry')
    steadystate = 1;
elseif strcmp(gapFilltype,'topology')
    steadystate = 0;
end

while any(temp_core)
    LPS = LPS+1;
    [flux1,~] = forward(model,temp_core,tol,steadystate);
    if sum(abs(flux))==0
        flux = flux1;
    else
        c1=unifrnd(0.45,0.55,1);
        flux = (c1*flux)+((1-c1)*flux1);
    end
    temp_core(coreRxns==1 & abs(flux)>=tol*1e-1) = 0; 
    if ~any(temp_core)
        break
    end
    LPS = LPS+1;
    [flux2,~] = reverse(model,temp_core,tol,steadystate);
    c1=unifrnd(0.45,0.55,1);
    flux = (c1*flux)+((1-c1)*flux2);
    temp_core(coreRxns==1 & abs(flux)>=tol*1e-1) = 0; 
end

direction = zeros(n,1);
direction(coreRxns==1&flux>0) = 1;
direction(coreRxns==1&flux<0) = -1;

if strcmp(probType,'MILP')
    if ~exist('solveTime', 'var') || isempty(solveTime)
        solveTime=7200;     
    end
    [reacInd,x] = findConsistentReacID(model,direction,weights,tol,steadystate,probType,solveTime,flux,prevSols);
elseif strcmp(probType,'LP')
    LPS = LPS+1;
    [reacInd,x] = findConsistentReacID(model,direction,weights,tol,steadystate,probType);
elseif strcmp(probType,'DC')
    [reacInd,x] = findConsistentReacID(model,direction,weights,tol,steadystate,probType);
end

flux = x(1:numel(model.rxns));
ConsModel = removeRxns(model, setdiff(model.rxns,model.rxns(reacInd)));
if remGene
    ConsModel = removeUnusedGenes(ConsModel);
end

if nSol>1
    ConsModel = {ConsModel};
    if strcmp(altSolMethod,'coreDirection')
        LPS = {LPS};    
        for j=2:nSol
            LPS_ = 0;
            temp_core = coreRxns;
            flux = zeros(n,1);
            while any(temp_core)
                LPS_ = LPS_+1;
                r = randi(2);
                if r==1
                    [flux1,~] = forward(model,temp_core,tol,steadystate);
                else
                    [flux1,~] = reverse(model,temp_core,tol,steadystate);
                end
                if sum(abs(flux))==0
                    flux = flux1;
                else
                    c1=unifrnd(0.45,0.55,1);
                    flux = (c1*flux)+((1-c1)*flux1);
                end
                temp_core(coreRxns==1 & abs(flux)>=tol*1e-1) = 0; 
                if ~any(temp_core)
                    break
                end
                LPS_ = LPS_+1;
                if r==1
                    [flux2,~] = reverse(model,temp_core,tol,steadystate);
                else
                    [flux2,~] = forward(model,temp_core,tol,steadystate);
                end
                c1=unifrnd(0.45,0.55,1);
                flux = (c1*flux)+((1-c1)*flux2);
                temp_core(coreRxns==1 & abs(flux)>=tol*1e-1) = 0; 
            end
        
            direction = zeros(n,1);
            direction(coreRxns==1&flux>0) = 1;
            direction(coreRxns==1&flux<0) = -1;
            if strcmp(probType,'MILP')
                if ~exist('solveTime', 'var') || isempty(solveTime)
                    solveTime=7200;
                end
                [reacInd,x] = findConsistentReacID(model,direction,weights,tol,steadystate,probType,solveTime,flux);
            elseif strcmp(probType,'LP')
                LPS_ = LPS_+1;
                [reacInd,x] = findConsistentReacID(model,direction,weights,tol,steadystate,probType);
            elseif strcmp(probType,'DC')
                [reacInd,x] = findConsistentReacID(model,direction,weights,tol,steadystate,probType);
            end
    
            flux = x(1:numel(model.rxns));
            Mod = removeRxns(model, setdiff(model.rxns,model.rxns(reacInd)));
            if remGene
                Mod = removeUnusedGenes(Mod);
            end 
            ConsModel = [ConsModel;Mod];
            LPS = [LPS;LPS_];
        end
    elseif strcmp(altSolMethod,'pathwayExclusion')
        if ~strcmp(probType,'MILP')
            error('The probType has to be MILP for using pathwayExclusion method of finding alternate models')
        end
        for j=2:nSol
            prevSols = [prevSols;find(reacInd)];
            [reacInd,~,stat] = findConsistentReacID(model,direction,weights,tol,steadystate,probType,solveTime,[],prevSols);
            if stat~=1
                break
            end
            Mod = removeRxns(model, setdiff(model.rxns,model.rxns(reacInd)));
            if remGene
                Mod = removeUnusedGenes(Mod);
            end
            ConsModel = [ConsModel;Mod];
        end
    end
end