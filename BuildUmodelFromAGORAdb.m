% To build universal model from the given AGORA models
% This code returns 2 variables: 
%       ConsUmodel: The required consistent universal model
%       TrblRxns: reactions that need attention due to conflict in bounds,
%                 stoichiometry and metids
clear
folder = './AGORA2/'; % path to the folder that has all the AGORA models

% Retrieving path to all the models
items=dir(folder); 
Path2AllModels = {};
modelNames = {};
for i=3:numel(items)
    p =[folder,items(i).name];
    modelNames=[modelNames;items(i).name];
    Path2AllModels=[Path2AllModels;p];
end

% Initializing the variables 'rxnsAttr' and 'metsAttr' that is the union of
% all the reactions and metabolites attributes respectively
load(Path2AllModels{1})
rxnFormulas=printRxnFormula(model,model.rxns,0);
rxnsAttr=struct('lb',model.lb,'ub',model.ub,'rxns',{model.rxns},'subSystems',{model.subSystems},...
    'rxnNames',{model.rxnNames},'rxnECNumbers',{model.rxnECNumbers},'rxnKEGGID',{model.rxnKEGGID});
metsAttr =struct('mets',{model.mets},'metNames',{model.metNames},'metCharges',{model.metCharges},...
    'metFormulas',{model.metFormulas},'metChEBIID',{model.metChEBIID},'metKEGGID',{model.metKEGGID},...
    'metPubChemID',{model.metPubChemID},'metHMDBID',{model.metHMDBID},'metSmiles',{model.metSmiles});
TrblRxns={}; % To store the reactions that has conflict among the models
bio1_formula={};
for i=1:numel(Path2AllModels)
    load(Path2AllModels{i})
    av_ids = find(ismember(model.rxns,rxnsAttr.rxns));
    nav_ids = setdiff([1:numel(model.rxns)],av_ids)';
    % Appending the new reactions
    AvrxnFormula = printRxnFormula(model,model.rxns(av_ids),0);
    NavrxnFormula= printRxnFormula(model,model.rxns(nav_ids),0);
    rxnFormulas =[rxnFormulas;NavrxnFormula];
    rxnsAttr = appendRxnAttr(model,nav_ids,rxnsAttr);
    % Appending the new metabolites
    nav_met_ids = find(~ismember(model.mets,metsAttr.mets));
    metsAttr = appendMetAttr(model,nav_met_ids,metsAttr);
    % Checking for conflicts in the reaction formula
    if sum(ismember(AvrxnFormula,rxnFormulas))~=numel(av_ids)
        Irxnid = av_ids(~ismember(AvrxnFormula,rxnFormulas));
        for k=1:numel(Irxnid)
            TrblRxns(end+1,1) = model.rxns(Irxnid(k));
            TrblRxns(end,2) = modelNames(i);
            TrblRxns(end,3) = {'Stoich or metid'};
        end
    end
    [C,iA,iB] = intersect(model.rxns(av_ids),rxnsAttr.rxns);
    % Checking for conflicts in the lower bound
    if sum(model.lb(av_ids(iA))==rxnsAttr.lb(iB))~=numel(av_ids)
        Irxn = C(model.lb(av_ids(iA))~=rxnsAttr.lb(iB));
        for k=1:numel(Irxn)
            TrblRxns(end+1,1) = Irxn(k);
            TrblRxns(end,2) = modelNames(i);
            TrblRxns(end,3) = {'lb'};
        end
    end
    % Checking for conflicts in the upper bound
    if sum(model.ub(iA)==rxnsAttr.ub(iB))~=numel(av_ids)
        Irxn = C(model.ub(iA)~=rxnsAttr.ub(iB));
        for k=1:numel(Irxn)
            TrblRxns{end+1,1} = Irxn(k);
            TrblRxns{end,2} = modelNames{i};
            TrblRxns{end,3} = 'ub';
        end
    end
    if strcmp(model.rxns(find(model.c)),'bio1')
        bio1_formula=[bio1_formula;printRxnFormula(model,model.rxns(find(model.c)),0)];
    end
end

% Removing the biomass reactions that are expected to vary in stoichiometry
temp_var = TrblRxns(:,1);
temp_var = TrblRxns(~contains(temp_var,'bio'),:);

% Creating a universal model
Umodel = createModel(rxnsAttr.rxns,rxnsAttr.rxnNames,rxnFormulas,'subSystemList',...
    rxnsAttr.subSystems,'lowerBoundList',rxnsAttr.lb,'upperBoundList',rxnsAttr.ub);

% Adding metabolite attributes to the universal model
[~,ord]=ismember(Umodel.mets,metsAttr.mets);
Umodel.metNames=metsAttr.metNames(ord);
Umodel.metCharges=metsAttr.metCharges(ord);
Umodel.metFormulas=metsAttr.metFormulas(ord);
Umodel.metChEBIID=metsAttr.metChEBIID(ord);
Umodel.metKEGGID=metsAttr.metKEGGID(ord);
Umodel.metPubChemID=metsAttr.metPubChemID(ord);
Umodel.metHMDBID=metsAttr.metHMDBID(ord);
Umodel.metSmiles=metsAttr.metSmiles(ord);

% Adding reaction attributes to the universal model
[~,ord]=ismember(Umodel.rxns,rxnsAttr.rxns);
Umodel.rxnKEGGID = rxnsAttr.rxnKEGGID(ord);
Umodel.rxnECNumbers = rxnsAttr.rxnECNumbers(ord);

% Adding all the bio1 reactions
for j=1:numel(bio1_formula)
    rxn_id = ['bio1',num2str(j)];
    Umodel = addReaction(Umodel,rxn_id,'reactionFormula',bio1_formula{j},'lowerBound',0,...
        'upperBound',1000);
end

% Getting the consistent universal model
a = speedcc(Umodel,1e-4);
ConsUmodel = removeRxns(Umodel,Umodel.rxns(setdiff([1:numel(Umodel.rxns)],a)));

% Removing biomass reactions
bio_rxns = ConsUmodel.rxns(contains(ConsUmodel.rxns,'bio'));
bio_rxns = setdiff(bio_rxns,{'EX_biomass(e)','pbiosynthesis'});
ConsUmodel = removeRxns(ConsUmodel,bio_rxns);


% Saving the consistent universal model
save('ConsUmodel','ConsUmodel')
% Saving the results
save('Results_BuildUmodelFromAGORAdb')

function rxnsAttr = appendRxnAttr(model,rxn_ids,rxnsAttr)
% This function appends the attributes of the reactions
rxnsAttr.lb=[rxnsAttr.lb;model.lb(rxn_ids)];
rxnsAttr.ub=[rxnsAttr.ub;model.ub(rxn_ids)];
rxnsAttr.rxns = [rxnsAttr.rxns;model.rxns(rxn_ids)];
rxnsAttr.subSystems = [rxnsAttr.subSystems;model.subSystems(rxn_ids)];
rxnsAttr.rxnNames = [rxnsAttr.rxnNames;model.rxnNames(rxn_ids)];
rxnsAttr.rxnECNumbers = [rxnsAttr.rxnECNumbers;model.rxnECNumbers(rxn_ids)];
rxnsAttr.rxnKEGGID = [rxnsAttr.rxnKEGGID;model.rxnKEGGID(rxn_ids)];
end
function metsAttr = appendMetAttr(model,met_ids,metsAttr)
% This function appends the attributes of the metabolites
metsAttr.mets=[metsAttr.mets;model.mets(met_ids)];
metsAttr.metNames=[metsAttr.metNames;model.metNames(met_ids)];
metsAttr.metCharges=[metsAttr.metCharges;model.metCharges(met_ids)];
metsAttr.metFormulas=[metsAttr.metFormulas;model.metFormulas(met_ids)];
metsAttr.metChEBIID=[metsAttr.metChEBIID;model.metChEBIID(met_ids)];
metsAttr.metKEGGID=[metsAttr.metKEGGID;model.metKEGGID(met_ids)];
metsAttr.metPubChemID=[metsAttr.metPubChemID;model.metPubChemID(met_ids)];
metsAttr.metHMDBID=[metsAttr.metHMDBID;model.metHMDBID(met_ids)];
metsAttr.metSmiles=[metsAttr.metSmiles;model.metSmiles(met_ids)];
end