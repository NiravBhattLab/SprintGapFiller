%% This code will generate a community model. The microbes in the community are randomly chosen from AGORA2 database
%% The parameter 'n' refer to number of microbes to be considered for model building
%%
clear
% initCobraToolbox(false)
% changeCobraSolver('ibm_cplex','all')

% Number of microbial models to build in the community
n=2; 
% Folder where all the microbial models are stored
folder='./AGORA2/'; 

% Listing path to all the models in the given folder
items=dir(folder); 
Path2AllModels = {};
for i=3:numel(items)
    p =[folder,items(i).name];
    Path2AllModels=[Path2AllModels;p];
end

% Loading Model names files (The reactions and metabolites will be named
% after this
load('ModelNames.mat')

% Choosing random model ids (community model will be built for only these
% models)
ids = sort(randsample(numel(ModelNames),n)); 
% ids=1:n;
Path2nModels = Path2AllModels(ids);
Models = ModelNames(ids);

% tolerance level (minimum flux required to be carried in all the models)
tol=1e-5;
% Loading the consistent universal model
load('ConsUmodel.mat')

% Media constraints given as bounds on exchange reactions
ids =startsWith(ConsUmodel.rxns,'EX_');
media=struct();
media.exc_rxns = ConsUmodel.rxns(ids);
media.lb= ConsUmodel.lb(ids);media.ub= ConsUmodel.ub(ids);

% Minimum flux to be carried by biomass reaction in all the microbes
minBio = 0.001;
% Whether or not to consider adding the trasport reactions that are not present in the microbial models
ConsiderOtherTranRxn=0;
% Whether or not to consider transport reactions as core reactions
TransferCore=0;
% Building a microbial community model
tic
Cmodel=BuildCommunityModels(Path2nModels, Models, ConsUmodel, minBio, ConsiderOtherTranRxn, TransferCore, media, tol);
toc

