function [MicComModel,removedCoreRxns] = BuildCommunityModels(Folder_path, abbr, Umodel, minBio, ConsiderOtherTranRxn, TransferCore, media, tol)
%%INPUT
%       Folder_path: Matlab cell listing the paths to all the microbial
%                    model structures in .mat format. Lower (*lb) and upper
%                    (*ub) bounds can be provided for any reaction in the 
%                    microbial models. Reaction ID (*rxns) and metabolite 
%                    ID (*mets) should be in same format as in Umodel. 
%                    Metabolite IDs(*mets) should include the compartment 
%                    info(Eg: glc_D(e), pyr(c))
%
%       abbr: Matlab cell listing model abbrevations. All rxns and mets 
%             will have this prefix. Must be same order as in Folder_path
%
%       Umodel: COBRA model structure of Universal model. This model should
%               be a superset of all the reactions and this model should 
%               not have any blocked reactions. All the exchange reactions
%               ID should begin with 'EX_'
%               The following fields are required:
%                   * S - `m x n` Stoichiometric matrix
%                   * b  - `m x 1` change in concentration with time
%                   * c  - `n x 1` Linear objective coefficients
%                   * lb - `n x 1` Lower bounds on net flux
%                   * ub - `n x 1` Upper bounds on net flux
%                   * mets - metabolite IDs
%                   * rxns - reaction IDs
%
%       minBio: minimum flux that a biomass reaction should carry
%
%       ConsiderOtherTranRxn: Boolean value
%               1: Transport reactions that are not present in the given 
%                  microbial model but present in the universal model will
%                  also be considered for inclusion
%               0: Transport reactions that are present in the microbial 
%                  model will only be considered for inclusion
%
%       TransferCore: Boolen value
%               1: Transport reactions that are present in the microbial
%                  model will also be considered as core reactions
%               0: Transport reactions that are present in the microbial
%                  model will not be considered as core reactions
%
%       media: matlab structure with fields
%              *exc_rxns: list of exchange reactions (must be in same 
%                         format as Umodel.rxns)
%              *lb: Lower bounds of corresponding reactions in exc_rxns
%              *ub: Upper bounds of corresponding reactions in exc_rxns
%
%       tol: minimum absolute flux value that every reaction in the 
%            community model should carry (default: 1e-4)

%%OUTPUT
%       ComModel: The consistent community model
%       
%       removedRxns: Reactions that were removed because of inconsistency 
%                    (or) inability to carry the minimum flux (tol) in the given
%                    community and media conditions

%%AUTHOR
%       Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

if ~exist('tol', 'var') 
    tol = 1e-4;
end
if numel(Folder_path) ~= numel(abbr)
    error('Model names has to be of the same size as the number of models in Foler_path')
end

n_models = numel(Folder_path); % number of models
Umodel_new=Umodel;
exc_rxns=Umodel_new.rxns(startsWith(Umodel_new.rxns,'EX_'));
exc_rxnFormulas = printRxnFormula(Umodel_new,'rxnAbbrList',exc_rxns,'printFlag',false);
Umodel_new=removeRxns(Umodel_new,exc_rxns); % remove all exchange rxns
S=[];lb=[];ub=[];c=[];b=[];rxns=[];mets=[];core=[];
for i=1:n_models
    load(Folder_path{i})
    if ConsiderOtherTranRxn
        UmodelTemp = Umodel_new;
    else
        % removing all the transport reactions that are not there in the
        % given microbial model
        UmodelTemp = removeRxns(Umodel_new,setdiff(Umodel_new.rxns(getTransRxns(Umodel_new)),...
            model.rxns(getTransRxns(model))));
    end
    if TransferCore
        coreTemp = ismember(UmodelTemp.rxns,model.rxns);
    else
        coreTemp = ismember(UmodelTemp.rxns,model.rxns(setdiff([1:numel(model.rxns)],getTransRxns(model))));
    end
    
    % Adding the biomass reaction
    BioForm = printRxnFormula(model,model.rxns(find(model.c)),0);
    UmodelTemp=addReaction(UmodelTemp,'bio','reactionFormula',BioForm{1},...
        'lowerBound',minBio,'upperBound',1000);
    
    core = [core;coreTemp;1]; % one refers to the biomass reaction
    new_rxns = cellfun(@(x)rename_rxns(x,abbr{i}),UmodelTemp.rxns,'uni',false);
    rxns = [rxns;new_rxns];
    new_mets=cellfun(@(x)rename_mets(x,abbr{i}),UmodelTemp.mets,'uni',false);
    mets=[mets;new_mets];
    S = blkdiag(S,UmodelTemp.S);c=[c;UmodelTemp.c];b=[b;UmodelTemp.b];
    new_lb = UmodelTemp.lb;new_ub = UmodelTemp.ub;
    [loca,locb] = ismember(model.rxns,UmodelTemp.rxns);
    locb = locb(locb~=0);
    new_lb(locb)=model.lb(loca);new_ub(locb)=model.ub(loca);
    lb=[lb;new_lb];ub=[ub;new_ub];
end

% merging the extracellular metabolite rows and removing the extra
% metabolites
[Umets,~,ix] = unique(mets);
counts = accumarray(ix,1).';
counts = counts';
for j=1:numel(counts)
    if counts(j)>1
        ids = find(ismember(mets,Umets{j}));
        S(ids(1),:)=sum(S(ids,:),1);
        S(ids(2:end),:)=[];
        b(ids(2:end))=[];
        mets(ids(2:end))=[];
    end
end
    
ComModel=struct();
ComModel.S=S;ComModel.lb=lb;
ComModel.ub=ub;ComModel.c=c;
ComModel.b=b;ComModel.mets=mets;
ComModel.rxns=rxns;

for i=1:numel(media.exc_rxns)
    id  = find(ismember(exc_rxns,media.exc_rxns{i}));
    ComModel=addReaction(ComModel,exc_rxns{id},'reactionFormula',exc_rxnFormulas{id},...
        'lowerBound',media.lb(i),'upperBound',media.ub(i));
end
core = [core;zeros(numel(media.exc_rxns),1)];

[a,feas,LPScc]=speedcc(ComModel,tol);

removedRxns = ComModel.rxns(setdiff([1:numel(ComModel.rxns)],a));
ConsUComModel = removeRxns(ComModel,removedRxns);
removedRxns = intersect(removedRxns,ComModel.rxns(find(core)));
[MicComModel,stat,LPSco] = speedcore(ConsUComModel,find(core(a)),tol);
% If no optimal solution is found then try including few reactions with flux value less than tol
while stat~=1
    [a1,feas,LPScc]=speedcc(ComModel,tol,feas/10);
    removedRxns = ComModel.rxns(setdiff([1:numel(ComModel.rxns)],a1));
    ConsUComModel = removeRxns(ComModel,removedRxns);
    removedCoreRxns = intersect(removedRxns,ComModel.rxns(find(core)));
    % Core reaction will always be those that can carry abs flux \geq tol
    [MicComModel,stat,LPSco] = speedcore(ConsUComModel,find(core(a)),tol);
end
end

function rxns = rename_rxns(a,ABR)
rxns = [ABR,'_',a];
end

function mets = rename_mets(a,ABR)
if ~strcmp(a(end-1),'e')&~contains(a,'biomass')
    mets = [ABR,'_',a];
else
    mets = a;
end
end

function ids = getTransRxns(model)
% this part of code is adapted from findTransRxns.m in COBRA toolbox
[~,compSymbols]=arrayfun(@parseMetNames, model.mets);
ids=[];
for i = 1:numel(model.rxns)
    % get compartment symbols for each rxn
    compSymbolsTmp=compSymbols(model.S(:,i)~=0);
    % if there's more than 1 compartment involved, it's a transport rxn
    if length(unique(compSymbolsTmp))>1
        ids=[ids;i];
    end
end
end