function [ConsModel,LPS] = sprintcore(model,coreRxns,tol,gapFilltype,coreMets,weights,nSol,altSolMethod,probType,solveTime,remGene,prevSols)
% USAGE:
%   [ConsModel,LPS] = sprintcore(model,coreRxns,tol,gapFilltype,coreMets,weights,nSol,altSolMethod,probType,solveTime,remGene,prevSols)
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
%                 'consistency' (stoichiometric matrix) based. (Default:'consistency')
%   coreMets:     Indices of metabolites that the model needs to produce. (Default: {}) 
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
if ~exist('coreMets', 'var') || isempty(coreMets)
    coreMets={};  
end

if ~isempty(coreMets) % if core metabolites are provided
    metFlag=1;
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
if sum(temp_core)==0 % if no core reactions are provided atleast core mets must be provided
    sol = optimizeCbModel(model);
    flux = sol.x;
else
    flux = zeros(n,1);
end

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


if metFlag % if there are any core metabolites
    coreRxnsMets = coreRxns;
    % checking if the variable 'flux' captured the core metabolites
    no_flux_met = [];
    for CmID = 1:numel(coreMets)
        curr_met = coreMets(CmID);
        corr_rxns = find(model.S(curr_met,:));
        if any(flux(corr_rxns))
            cons_corr_rxns = intersect(find(flux),corr_rxns);
            % updating the coreRxns based on the obtained new rxns from coreMets
            coreRxnsMets(cons_corr_rxns) = 1;
        else
            no_flux_met = [no_flux_met;coreMets(CmID)]; % if any core metabolite is not captured in the 'flux'
        end
    end
    if any(no_flux_met) % for the core metabolites which are yet to be added
        rxnsTemp = find(model.S(no_flux_met(1),:));
        for tempId =2:numel(no_flux_met)
            rxnsTemp = union(rxnsTemp,find(model.S(no_flux_met(tempId),:)));
        end
        temp_core = ismember(1:numel(model.rxns),rxnsTemp)';
        while any(temp_core)
            LPS = LPS+1;
            [flux1,~] = forward(model,temp_core,tol,steadystate);
            c1=unifrnd(0.45,0.55,1);
            flux = (c1*flux)+((1-c1)*flux1);
            [temp_core,new_core,no_flux_met] = get_unblocked_metRxns(flux,model,no_flux_met,temp_core);
            coreRxnsMets(new_core) = 1;
            if ~any(temp_core)
                break
            end
            LPS = LPS+1;
            [flux2,~] = reverse(model,temp_core,tol,steadystate);
            c1=unifrnd(0.45,0.55,1);
            flux = (c1*flux)+((1-c1)*flux2);
            [temp_core,new_core,no_flux_met] = get_unblocked_metRxns(flux,model,no_flux_met,temp_core);
            coreRxnsMets(new_core) = 1;
        end
    end
end

direction = zeros(n,1);
if metFlag
    direction(coreRxnsMets==1&flux>0) = 1;
    direction(coreRxnsMets==1&flux<0) = -1;
else
    direction(coreRxns==1&flux>0) = 1;
    direction(coreRxns==1&flux<0) = -1;
end

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
            if sum(temp_core)==0 % if no core reactions are provided atleast core mets must be provided
                sol = optimizeCbModel(model);
                flux = sol.x;
            else
                flux = zeros(n,1);
            end

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
            

            if metFlag % if there are any core metabolites
                coreRxnsMets = coreRxns;
                % checking if the variable 'flux' captured the core metabolites
                no_flux_met = [];
                for CmID = 1:numel(coreMets)
                    curr_met = coreMets(CmID);
                    corr_rxns = find(model.S(curr_met,:));
                    if any(flux(corr_rxns))
                        cons_corr_rxns = intersect(find(flux),corr_rxns);
                        % updating the coreRxns based on the obtained new rxns from coreMets
                        coreRxnsMets(cons_corr_rxns) = 1;
                    else
                        no_flux_met = [no_flux_met;coreMets(CmID)]; % if any core metabolite is not captured in the 'flux'
                    end
                end
                if any(no_flux_met) % for the core metabolites which are yet to be added
                    rxnsTemp = find(model.S(no_flux_met(1),:));
                    for tempId =2:numel(no_flux_met)
                        rxnsTemp = union(rxnsTemp,find(model.S(no_flux_met(tempId),:)));
                    end
                    temp_core = ismember(1:numel(model.rxns),rxnsTemp)';
                    while any(temp_core)
                        LPS_ = LPS_+1;
                        [flux1,~] = forward(model,temp_core,tol,steadystate);
                        c1=unifrnd(0.45,0.55,1);
                        flux = (c1*flux)+((1-c1)*flux1);
                        [temp_core,new_core,no_flux_met] = get_unblocked_metRxns(flux,model,no_flux_met,temp_core);
                        coreRxnsMets(new_core) = 1;
                        if ~any(temp_core)
                            break
                        end
                        LPS_ = LPS_+1;
                        [flux2,~] = reverse(model,temp_core,tol,steadystate);
                        c1=unifrnd(0.45,0.55,1);
                        flux = (c1*flux)+((1-c1)*flux2);
                        [temp_core,new_core,no_flux_met] = get_unblocked_metRxns(flux,model,no_flux_met,temp_core);
                        coreRxnsMets(new_core) = 1;
                    end
                end
            end
        end
            direction = zeros(n,1);
            if metFlag
                direction(coreRxnsMets==1&flux>0) = 1;
                direction(coreRxnsMets==1&flux<0) = -1;
            else
                direction(coreRxns==1&flux>0) = 1;
                direction(coreRxns==1&flux<0) = -1;
            end
            
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


function [temp_core,new_core,blockMets] = get_unblocked_metRxns(flux,model,no_flux_met,temp_core)
    new_core = intersect(find(temp_core),find(flux));
    cons_rxn_ids = find(flux);
    blockMets = setdiff(no_flux_met,find(sum(model.S(:,cons_rxn_ids)~=0,2)));
    if isempty(blockMets)
        temp_core = zeros(numel(model.rxns),1);
    else
        rxnsTemp = find(model.S(blockMets(1),:));
        for tempId =2:numel(blockMets)
            rxnsTemp = union(rxnsTemp,find(model.S(blockMets(tempId),:)));
        end
        temp_core = ismember(1:numel(model.rxns),rxnsTemp);
    end
    
end