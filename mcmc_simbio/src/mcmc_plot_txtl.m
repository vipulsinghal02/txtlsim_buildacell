function mcmc_plot_txtl(m, p, d, varargin)
% Generate txtl_plot plot using user specified parameters, for models used
% in the mcmc_simbio toolbox (ie, those with parameters already globalized
% and bound using the Kd rules). If the parameters are not globalized, they
% *must* be specified in the p struct that is input into this function,
% along with the reaction string specifying which reaction that parameter
% is to be scoped under.
%
% This function is useful when you have results from mcmc and want to just
% plot an unexported Simbiology model object with those parameters.
%
% INPUTS:
% m: the Simbiology model object you want to simulate. This needs to have
% gone through the setup reactions phase and globalize parameters phases
% already.
%
% p: A struct specifying the parameters you want to use for the simulation.
% This has fields:
%       paramNames: name of parameter
%       paramVals: value of the parameter; this is not log transformed.
%       reactionString: This is either 'global' for globally (model object
%                  level) parameters, or a string corresponding to the
%                  reaction that parameter belongs to.
%
% d: A struct specifying the dose to use. There are two fields: 'dosedVals'
% and 'dosedNames'. dosedVals is a column vector of values to use, and
% dosedNames is a cell array of the names of the species to be dosed.
%
% varargin: optional arguments.
%
%
% Name value pair:
% 1. 'model_info', model_info(j) (One element of model_info struct array)
% inputs the jth elements of the model_info struct array to this function.
% The model_info struct contains the
% parameters to be used with the model, These get overwritten by any
% parameters that are also specified int eh parameter struct p provided to
% this function. The model info struct must be for a single topology, i.e.,
% must be a struct of length 1. For multiple topology model info struct,
% slice them as mi(j) before passing to this function.
%
% 2. 'which_geometry', g (numeric)
% If there are multiple geometries (for example, columns of paramMaps) in the model
% info struct, this one specifies which one to use to set the model
% parameters. If unspecified, then the first geometry is used in the
% model_info.
%
% 3. 'master_info', mai (struct)
% This is the master_info struct that contains the parameters used by the
% model_info struct (namesUnord is the array of names, and paramMaps is the
% map that takes the master vector to the model's parameters specified by
% the names in model_info.namesUnord
%
% Vipul Singhal, Caltech, 2018


% parse optional inputs
inp = inputParser;
inp.addParameter('model_info', [], @isstruct); %
inp.addParameter('which_geometry', 1, @isnumeric); %
inp.addParameter('master_info', [], @isstruct); %
inp.parse(varargin{:});
inp = inp.Results;


% set the parameters in the model object - first from the model_info, then
% from the parameters struct p.
if ~isempty(inp.model_info) && ~isempty(inp.master_info)
    [m, pnotset] = set_mobj_from_model_info(m, ...
        inp.model_info, ...
        inp.master_info, ...
        inp.which_geometry);
end

% now set the parameters from the supplied parameter array.
for i = 1:length(p)
    obj = sbioselect(m, 'name', p(i).paramNames);
    if ~isempty(obj)
        switch class(obj)
            case 'SimBiology.Parameter'
                obj.Value = exp(p(i).paramVals);
            case 'SimBiology.Species'
                obj.InitialAmount = exp(p(i).paramVals);
        end
    else
        % maybe the parameter is reaction level scoped.
%         rObj = [];
%         if ~strcmp(p(i).reactionString, 'global')
%             rObj = sbioselect(m,'Type','reaction','reaction',p(i).reactionString);
%             
%             %             if it is a reaction, then it could be reveresible or
%             %             irreversible. Then we need to set the parameter within the
%             %             kinetic law object. I dont have time for this. 
%             
%             
%             %             rObj.Value = log(p(i).paramVals);
%             %             !todo...
%             error('todo')
%         end
        error('paramnotset_see_what_happened.')
        % otherwise just report that it did not get set.
%         paramsNotSet = [paramsNotSet; pnames(i)];
%         warning(['Object with name ' pnames{i} ' not found. Skipping setting its value.'])
    end
%     can also look into: setparam(m, p(i).reactionString, p(i).paramNames, p(i).paramVals);
end

% now run the model under the specified dose and record the output.
for i = 1:length(d.dosedNames)
    sp = sbioselect(m, 'Name', d.dosedNames{i});
    sp.InitialAmount = d.dosedVals(i);
end
[simData] = txtl_runsim(m,14*60*60);


% finally plot stuff.
txtl_plot(simData,m);




% Save then in the current directory as a low res jpg,
% displaying the directory they were saved in in the command window.
% If a file with the current name exists, then append a numeric
% to the filename.
print(['txtlplot_' datestr(now, 'yyyymmddTHHMMSS')],'-djpeg', '-r100')

end







