function [fhandle] = mcmc_plotCustomSpecies(em, m, mi, mai, di, spl)
%% DONT NEED THIS FUNCTION> CAN JUST USE MCMC_TRAJECTORIES WITH A CUSTOM SPECIES ARGUMENT. 
%mcmc_plotCustomSpecies Take an exported model, and a parameter array (2D), and plot the
%trajectories for up to the last 5 parameter points in the parameter array
%   -- em: exported model object. 
%   -- m: parameter array, ORDERED (ie, parameters ordered by the following
%           type of code: 
%               mvarray = masterVecArray(m, mai); 
%               marrayOrd = mvarray(mi.paramMaps(mi.orderingIx, 1),:,:);
%               and then made 2D if needed. 
%           This array has dimensions: nPoints x nParam
%           last 5 points get plotted. If there are less than 5 points, then
%           all of them get plotted. 
%   -- mi: model_info struct. Must be a scalar struct. 
%   -- mai: master_info struct. must be a scalar struct. 
%   -- spl: cell array vector list of cell array vector lists of species to plot. 
%           The outer cell array is the number of subplots. The inner cell
%           array is the list of species whose trajectories get summed in
%           the plotting. This is the same specification as the
%           measuredSpecies field in the model_info struct. 
%   

% determine the subplot dimensions from the 
nParam = size(m, 2); 
[n1 n2] = twofactors(nParam);
if isprime(nParam)
    n1 = ceil(nParam/5); % 5 columns
    n2 = 5;
end
nPoints = size(m, 1);
if nPoints >4
    sIx = nPoints-4;
    nSimCurves = 5;
    
else
    sIx = 1;
    nSimCurves = nPoints;
end

tv = di(mi.dataToMapTo).timeVector;
dose = mi.dosedVals;

% simulate the model at the required points. 
[da, idxnotused] = simulatecurves(em,m, nSimCurves, dose, tv, spl);

for i = sIx:nPoints
    
    
    
    
    
% plot the results




end

