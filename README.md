# SprintGapFiller
Codes for "SprintGapFiller: For fast gapfilling in microbial community models"   

### Requirements
1. MATLAB
2. [COBRA Toolbox](http://opencobra.github.io/cobratoolbox/)

# Follow the instructions below to get a consistent gap-filled microbial community model (microbes are randomly chosen from AGORA2)
### Initiating the cobra toolbox
```
initCobraToolbox(false)
```
### Choosing the solver type
```
changeCobraSolver('ibm_cplex','all')
```
### Changing the feasible tolerance value
```
changeCobraSolverParams('LP', 'feasTol', 1e-8);
```
### Parameters required to build a microbial community model
```
n = 100; % Number of microbes in the community
tol = 1e-3; % Minimum flux required for a reaction to be unblocked
folder = './AGORA2/'; % Path to AGORA model files
minBio = 0.001; % Minimum flux to be carried by biomass reaction in all the microbes
ConsiderOtherTranRxn=0; % Whether or not to consider adding the trasport reactions that are not present in the microbial models
TransferCore=0; % Whether or not to consider transport reactions as core reactions
```
### Getting path all the AGORA2 model files
```
% Listing path to all the models in the given folder
items=dir(folder); 
Path2AllModels = {};
for i=3:numel(items)
    p =[folder,items(i).name];
    Path2AllModels=[Path2AllModels;p];
end
```
### Loading the names file and the universal model
```
load('ModelNames.mat') % Loading Model names files (The reactions and metabolites will be named after this)
load('ConsUmodel.mat') % Loading the consistent universal model
```
### Setting the media constraints
```
% Media constraints given as bounds on exchange reactions
ids =startsWith(ConsUmodel.rxns,'EX_');
media=struct();
media.exc_rxns = ConsUmodel.rxns(ids);
media.lb= ConsUmodel.lb(ids);media.ub= ConsUmodel.ub(ids);
```
### Choosing random microbial models and building a community 
```
ids = sort(randsample(numel(ModelNames),n)); 
Path2nModels = Path2AllModels(ids);
Models = ModelNames(ids);
Cmodel=BuildCommunityModels(Path2nModels, Models, ConsUmodel, minBio, ConsiderOtherTranRxn, TransferCore, media, tol);
```

