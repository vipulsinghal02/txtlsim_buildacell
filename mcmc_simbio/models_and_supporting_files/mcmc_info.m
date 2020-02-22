function mcmc_info
% mcmc_info.m - This help file describes the contents and function of the 
% mcmc_info struct. You can type help mcmc_info in the command line to view
% its contents. It does not contain any code, and is for reference only. 
% 
%
% OVERVIEW: 
% 
% The mcmc_info struct contains information on how to set up the parameter
% inference problem using the simbiology model(s), dataset(s),
% parameter(s), estimation hyperparameters, the parameter 'concurrence' or 
% 'sharing' pattern, and other information. 
% 
% 
% PARAMETER CONCURRENCE, TOPOLOGIES AND GEOMETRIES:
% 
% In general, each estimation problem can have multiple experiment-model
% pairs that inform some common set of parameters. For example, we might
% have an estimtion problem where the tetR transcription factor shows up in
% two different circuits, say, the incoherent feedforward loop (IFFL) and the
% genetic toggle switch. Then, many parameters associated with the tetR
% transcription factor will be expected to appear in in models
% corresponding to both circuits. In general, we do not expect the
% parameter value associated primarily with the of action of tetR (we will 
% call such parameter values "part specific" or "circuit specific") to be
% different between different circuits. Therefore, if we have data on both
% circuits, and want to fit a model of each circuit to its corresponding
% data, we will want to estimate the values of the tetR-associated
% parameters jointly. However, there will also be many parameters that will
% be specific to each model, and those will be estimated independently in
% the two models. The parameters that are unique to models (ie not shared)
% need not come from models of different circuits; they can be from models 
% of the same circuits as well. For example, that there are two copies of 
% the toggle circuit, one in a fast growing strain, and the second in a 
% slow growing stain, so that the dilution rate parameters (for instance) 
% in the two models of toggle switch are expected to be different. In
% general such 'environment specific' effects can lead to circuits having
% models that differ in parameters whose values would otherwise be equal. 
% 
% In general, we will have an arbitrary number of models,
% and an arbitrary pattern of paramter sharing/concurrence between the
% models. The mcmc_info file is used to specify the list of models, the
% datasets these models correspond to, the pattern of parameter concurrence
% in the parameter inference problem, the experimental setup (dosing
% information and which species get measured, for example), the MCMC
% hyperparameters (like number of steps, number of walkers, step size,
% etc). 
% 
% We also introduce some additional terminology that will be ueful at
% various points. 
% 
% TOPOLOGY: Different circuits lead to models that have different network
% topologies. We will say that these models have different topologies. 
% In the example above, there are two model
% topologies, the IFFL topology and the genetic toggle switch topology.
%
% GEOMETRY: Models with the same topology can still differ if any parameters
% are expected to have different values in different contexts. We say that
% such models are geometrically different. An example of such a difference
% was described in the example above when the genetic toggle was implemented
% in the two strains with different growth characteristics. The dilution
% parameter in the models of the genetic toggle would be expected to be
% different in the two strains. In this case we say that the genetic toggle
% topology has two geometries associated with it
% 
% In general, we can identify a model with a pair of indices, (i, j), for
% the ith topology and the jth geometry. In the above, let the IFFL be
% topology 1, and the genetic toggle be topology 2, and so the three models
% described above may be identified with (1, 1); (2, 1); (2, 1), and each 
% of these in general will have a dataset associated with it. We can
% then describe the parameter concurrence pattern in terms of these models.  
% For example, the tetR specific parameters get shared between all three 
% models, the dilution rate parameters are given different identities
% (i.e., are independently estimated) between model (2, 1) ad (2, 2). We
% will discuss how the concurrence pattern is specified when we discuss the
% paramMaps array.
% 
% The terminology "topology" was chosen to reflect the fact
% that the models corresponding to the two different circuits
% have different reaction network topologies. Similarly, the term
% "geometry" is used in contrast to topologies because the circuits that
% are geometrically different have models that are different in at most the
% parameters, and never in the set of chemcial reaction equations (i.e.,
% never in the reaction network topology). 
% 
% Each topology will have at least one geometry, and in general may have
% multiple. All the geometries within a topology can *only* differ in the
% parameter values, and not in the initial conditions specification (dosing)
% or the identities of the species being measured (the measuresSpecies property). 
% Note that there are two ways that initial conditions (dosing) can vary:
% We can vary the identities of the species that are dosed, or, even if the
% species being dosed are the same between two models, the set of initial 
% values (dose values) that these species take could be different. In both
% these cases, the two models are considered to have different topologies.
% 
% 
% Examples of files that generate mcmc_info structs for a few different
% cases are: 
% mcmc_info_constgfp3i.m -- constitutive gene expression, single topology
% and single geometry
% 
% mcmc_info_constgfp3ii.m -- constitutive gene expression, single topology,
% two geometries. Certain parameters shared between the two geometries. 
% i.e., at each iteration of the parameter inference algorithm, only a single 
% point was used for the shared parameters. Eventually, the shared parameters
% had a common distribution for both the geometries. Parameters that were
% not shared had separate values estimated. 
% 
% mcmc_info_tetR1i.m -- a single topology-geometry setup for the tetR
% repression circuit. 
% 
% mcmc_info_constgfptetR1.m -- This estimation problem involves two
% topologies, with each topology having one geometry. The topologies are
% the constitutive gene expression model, and the tetR repression model.
% This example involves a case where both these circuits are implemented in
% the same environment, so that any parameters associated with the
% environment (such as the concentration of enzymatic machinery or transcript 
% elongation rate) are expected to be the same in both models, while any
% parameters associated with circuit parts (such as tetR dimerization or 
% promoter binding rate) that are only present in one of the models are
% (trivially) unique to each model, and estimated independently of the
% parameters in the other model. 
% 
% THE MCMC_INFO STRUCT
% 
% In this section, we describe the fields within the mcmc_info struct. The
% mcmc_info struct has three fields, each of which is a matlab struct in
% itself: 
% 
% - runsim_info: contains information on the MCMC run hyperparameters like
% the stepsize, number of steps, number of walkers, whether parallelization
% is to be used, etc. 
% - model_info: contains information on the models (topologies and
% geometries) to be used in the estiamtion problem. 
% - master_info: contains information on the full set of parameters to be
% estimated (collected from all the models), including their names, the
% ranges of parameter values to use, etc.
% 
% We now describe these three structs in greater detail. 
%
% MCMC_INFO.RUNSIM_INFO
% The runsim_info struct has the fields:
% 
% - stdev: The standard deviation used in the computation of the log
% likelihood. We use a gaussian likelihood function (i.e., the L2 norm of 
% the residuals, weighted by the standard deviation, when the log of the
% likelihood is computed). 
% 
% - tightening: This is simply an additional parameter that the stdev value
% gets divided by. A large value means that the MCMC has a "tighter"
% rejection criterion, i.e., it rejects likelihood decreasing proposals 
% more easily. Recall that in the Metropolis algorithm, if a
% proposed parameter point increases the likelihood, it always gets
% accepted, but if it decreases the likelihood, it gets rejected
% probabilistically. Thus, with a higher value of the tightening parameter,
% the rejection happens with a higher probability. In the limit of large
% values of this parameter, MCMC reduces to a deterministic method which
% always accepts proposals that increase the likelihood, and always rejects
% proposals that decrease the likelihood. We generally just set this value
% to a default of 1. 
% 
% - nW: The number of MCMC walkers for the ensemble MCMC. A good value is
% 200 - 600 walkers for most 16Gb RAM machines. A larger number of walkers,
% if your machine can handle it, is always better. See the original paper
% by Goodman and Weare https://projecteuclid.org/euclid.camcos/1513731992
% and the python implementation https://arxiv.org/abs/1202.3665 and
% especially the dicumentation pages of the python implementation,
% http://dfm.io/emcee/current/ , for an in depth discussion. 
% 
% - stepsize: The step size. Do not use a value of 1, it will not work. All
% other positive values should be fine. The simulation depends on this
% hyperparameter quite critically, so test a few values and see what gives
% you a rejection percentage of less than 90%. See the links above for more
% discussion. 
%
% - nPoints: The number of times the forward model gets simulated in a
% single iteration (see nIter below). The actual number of steps reported
% by the algorithm for each iteration = nPoints / (nW*thinning). It is best
% not to set this value to be larger than 300000 (three hundred thousand)
% points, since at numbers significantly larger than this, the machines we
% tested the toolbox on tended to either freeze or not run the simulation.
% We have defined a separate hyperparameter, nIter, that allows for the
% MCMC runs to scale to an arbitrary number of points without these
% problems. 
% 
% - nIter: Number of iterations to run the MCMC for. It turns out that if
% we want to run the MCMC for 3 million points, it is better to set nPoints
% to 300k points, and nIter to 10, so that the MCMC simulations are run 10
% times, with each iteration a continuation of the previous (as it would
% have been if a 3 million point run was possible). A typical value is 10, 
% though 5 or 50 are not uncommon, depending on your needs. 
% 
% - thinnning: number of steps to skip before recording the walker
% positions. A typical value is 10. 
% 
% - parallel: A boolean valued flag, specifying if the parallel computing
% capabilities of MATLAB should be used. 
% 
% MCMC_INFO.MODEL_INFO
% 
% The model_info struct is a struct array of length nTopologies, 
% and specifies the properties of models, and the pattern of parameter 
% sharing across the topologies and geometries for the purposes of setting
% up the concurrent parameter inference problem. Each element of this array
% corresponds to one model topology. The fields of the struct define
% various properties associated with that topology. The ith topology's
% properties may be set / accessed as mcmc_info.model_info(i), with fields:
% 
% - circuitInfo: A human readable description of the model / circuit /
% topology. Use the text formatting you would use for the function fprintf.
% 
% - modelObj: A Simbiology model class object. In the terminology 
% of the concurrent parameter inference problem introduced above, this is a 
% network topology that this element of the struct array corresponds to. 
% 
% - namesUnord: A list of strings naming the parameters in the model object 
% that are set from values in the master_vector (see the 
% mcmc_info.master_info struct). Use this list to specify which parameters 
% in the model need to be set by the values from the master_vector. The 
% paramMaps array defined below is used to map the values from the 
% master_vector to the parameters listed in namesUnord. 
% 
% - paramMaps: This is one of the most important fields as far as the
% defining the concurrence feature is concerned. This is a matrix of size
% length(mcmc_info.master_info.namesUnord) x nGeometries, where nGeometries
% is the number of geometries that are to be defined for this topology.
% Each column of this matrix contains a list of indices of the elements of
% master_vector that are used to specify the values of the parameters named
% in namesUnord (in the order specified in namesUnord) for one geometry. 
% 
% 
% For example, suppose the current model has three parameters, with name
% strings 'k1', 'k2', and 'k3',
% and we have three geometries: g1, g2 and g3. For simplicity, assume that
% there is only one topology (with the straightforward generalization to an
% arbitrary number of topologies, see for example mcmc_info_tetR1i.m). 
% Rename the master_vector to be mv, and let it have 7 entries: 
% 
% mv = [p1, p2, p3, ..., p7]';
% 
% where p1, ..., p7 are numerical values. Then if a paramMaps matrix is 
% defined as paramMaps = [1  1  1 
%                 2  4  6
%                 3  5  7],
% 
% leads to an estimation problem with three geometries, with these geometries 
% having the relationship with the elements of mv depicted
% in the table below.   
% 
%           \  
%            \ Geometry
%   Parameter \
%                   g1      g2      g3
%                _____________________________
%      k1       | mv(1)    mv(1)   mv(1)
%      k2       | mv(2)    mv(4)   mv(6)
%      k3       | mv(3)    mv(5)   mv(7)
% 
% The master vector mv has 7 entries, and these may be partitioned into
% entries that are either fixed or to be estimated. For example, here we
% might want to set mv(1) to be fixed, and mv(2:7) to be estimated. Then,
% the value of mv(1) must be specified in the file that sets up the
% mcmc_info struct, while the values of elements mv(2:7) are the values that 
% are to be estimated, and are iteratively set each time the MCMC algorithm 
% tries a new point in the 6-D parameter space being searched. Thus, these 
% values do not matter at the time the mcmc_info stuct is specified, and
% can be set to anything initially. See hte MCMC_INFO.MASTER_INFO section 
% below, where this and all other specification related to the master 
% vector are made. 
% 
% Overall, the MCMC algorithm will search a 6-dimensional space, and at
% every point that is a proposal, distribute the 7 values (1 fixed, and the 
% remaining 6 specified by the 'current' proposed point) 
% into the three geometries as described in the table above, and simulate 
% these models to compute, using the dataset corresponding to each geometry, 
% the residuals and the likelihood. 
% 
% - dosedNames: A cell array of strings representing the names of the species
% that are to be 'dosed'. Dosing here means that the initial values of these
% species are varied as part of the experiment, and a set of trajectories
% are generated for a given parameter. These set of trajectories are then
% compared to corresponding data. The length of this array of strings is
% given by nDosedSpecies. See also the dosedVals field below. 
% 
% - dosedVals: A matrix of (numerical) dose values, of size nDosedSpecies by 
% nDoseCombinations. Here nDosedSpecies is the length of the dosedNames
% (cell) array of strings. nDoseCombinations is the number of different
% doses (initial conditions) the model is to be simulated for. Each dose
% combination column vector is used to set the initial values of the species
% specified in the dosedNames array, and a model trajectory is generated.
% Thus, there is one trajectory associated with each dose combination. We
% note that in the data (defined in the data_info struct array), there must
% be a corresponding number of dose combinations. 
% 
% - measuredNames: This is a cell array of (sub) cell arrays of strings 
% containing the names of the species in the model. For example, if we have
% {{'species a'}, {'species b', 'species c'}} then the first output of the
% model is the trajectory of species a, and the second output is the sum of
% the trajectores of species b and species c. These two outputs will 
% correspond to the first column and the second column of the data array
% (resp.). Note that the strings 'species x' must correspond to
% species in the model object.
% 
% - measuredSpeciesIndex: This is an array of indices mappints the cells of
% the measuredNames field to the columns of the dataArray field of the
% data_info struct's relevant entry (the relevant entry specified by the
% dataToMapTo array defined below). 
% 
% - dataToMapTo: A numerical vector array that maps the data_info struct 
% array's elements (type help data_info into the command line to read more 
% about how the data is stored using data_info structs) to the geometries
% defined by the various columns of the paramMaps matrix. It therefore has
% length nGeometries, where recall that nGeometries is the number of 
% geometries associated with the current topology, and defined by the columns 
% of the paramMaps struct. 
% 
% 
% MCMC_INFO.MASTER_INFO
% 
% The master_info struct is used to specify information about the
% master_vector, the elements of the master vector that are to be
% estimated, and other global information about the parameter inference
% problem. It is a scalar struct, and has fields: 
% 
% estNames: cell array of strings specifying the names of the parameters or
% species within the master vector that are to be estimated. You cannot
% estimte species that are also being dosed (ie, are in the dosedNames
% array in one of the elements in the mcmc_info.model_info struct array. 
% 
% masterVector: A numerical vector of parameter and species initial 
% concentration values that are to be distributed to all the models (over
% all the topologies and geometries) as specified by the various paramMaps
% matrices in the elements of the mcmc_info.model_info struct array. These
% values are log transformed, and range from -inf to +inf, even when the
% parameter values and species concentration values are restricted to be
% non-negative. Not all the values in master vector need to be estimated. 
% The entries that are fixed are specified by the fixedParams field of 
% this struct. 
% 
% paramRanges: This is a length(estParams) x 2 matrix array setting the
% upper and lower bound on the (log transformed) values of the parameters
% and species being estimated. 
% 
% - fixedParams: This is a vector of length 
% length(masterVector) - length(estParams) that contains indices of the
% elements of the master vector that are to be fixed. These elements of the
% master vector must be specified in the file that sets up the mcmc_info 
% struct. See mcmc_info_constgfp3i.m for an example. There, kf is set to a 
% fixed value, and the remaining parameters are estimated. The fixedParams
% vector is therefore just a scalar containing the index of the kf entry of
% the master vector. 
% 
% - semanticGroups: This is an experimental feature (in the sense that "we 
% are trying it out", rather than it "pertains to experimental data") that
% could be useful. We will supply more complete documentation on how to 
% use it in future versions. For now, the idea is that when the 
% integrability of the initial spead of points is checked using 
% integrableLHS_v2.m, we might want that 
% 
% See the files mcmc_info_constgfp3i.m,  mcmc_info_constgfp3ii.m and 
% mcmc_info_tetR1i.m for concrete functional examples of setting up an
% mcmc_info struct. 

% Copyright (c) 2018, Vipul Singhal, Caltech
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.




%% OLD STUFF: IGNORE. 
% 
% 
% that are 
% % defined for this topology.
% 
% 
% 
% nCols 
% % containing the indices of the elements of the data_info
% % struct that the model geometries correspond to. These are used when
% % the model predictions are compared to the data in the computation 
% % of the log likelihood during MCMC. 
% % 
% % - 
% nCols
% 
% % The measuredNames and measuredSpeciesIndex arrays are used to pick out the 
% % species coordinates of the simulation trajectories that are to be used as
% % outputs of the model, which in turn are compared to experimental data
% % trajectories. 
% % cell array of strings (defined next) then 
% 
% sets the values 
% % simulation 
% % 
% % species. It has length nDosedSpecies. 
% % 
% % 
% 
% 
% % have 
% 
% 
% 
% % 
% 
% 
% I.e., if we use the same circuit, but
% % change the species that we are varying 
% 
% 
% 
% % Furthermore, suppose the IFFL was initialized ("doesd") at a few
% % different initial condition vectors, and the corresponding output data
% % was collected. 
% 
% 
% Each circuit will have a different model and data 
% % associated with it, but the tetR dimerization rate parameters (for example) 
% % will show up in both models. In general, we do not want to estimate two
% % separate values (or *distributions* of values in the Bayesian parameter 
% % estimation scenario, especially when the parameters are non-identifiable)
% % from each model-data pair. Instead, we want to only estimate a single
% % value / distribution from both sets of information. We call models that
% % differ in the network topology (i.e., the set of chemical reactions
% % defining the model
% 
% 
% 
% 
% multiple different model
% % topologies, and each model topology might have multiple 'variants' of
% % that topology associated with it. We call these variants 'geometries', to
% % distinguish from the variation in models due to the topology of the
% % network itself. Here, there is only one network topology and geometry.
% % Type 'help mcmc_info' in the command line to learn more about the general
% % use case, and see the files mcmc_info_constgfp3tetR1 for an example of
% % the use of two distincti topologies in a single estimation problem. 
% 
% 
% Variations can For example, we might have different
% % fixed parameters associated with it. We will call these parameters
% % geometries associated with that model, since the network topology of the
% % model is the same, but the fact that certain parameters can be different
% % is analogous to the networks being geometrically different. 
% 


% In general, each estimation problem can have multiple experiment-model
% pairs that inform some common set of parameters. For example, we might
% have an estimtion problem where the tetR transcription factor shows up in
% two different circuits, say, the incoherent feedforward loop and the
% genetic toggle switch. Each circuit will have a different model and data 
% associated with it, but the tetR dimerization rate parameters (for example) 
% will show up in both models. In general, we do not want to estimate two
% separate values (or *distributions* of values in the Bayesian parameter 
% estimation scenario, especially when the parameters are non-identifiable)
% from each model-data pair. Instead, we want to only estimate a single
% value / distribution from both sets of information. We call models that
% differ in the network topology (i.e., the set of chemical reactions
% defining the model




multiple different model
% topologies, and each model topology might have multiple 'variants' of
% that topology associated with it. We call these variants 'geometries', to
% distinguish from the variation in models due to the topology of the
% network itself. Here, there is only one network topology and geometry.
% Type 'help mcmc_info' in the command line to learn more about the general
% use case, and see the files mcmc_info_constgfp3tetR1 for an example of
% the use of two distincti topologies in a single estimation problem. 


Variations can For example, we might have different
% fixed parameters associated with it. We will call these parameters
% geometries associated with that model, since the network topology of the
% model is the same, but the fact that certain parameters can be different
% is analogous to the networks being geometrically different. 




end