function m = model_txtl_pLacdeGFP(varargin)
% model_txtl_ptetdeGFP_pLactetR_aTc: txtlsim model containing: 
% ptet-UTR1-deGFP
% pLac-UTR1-tetR
% aTc
% 
% 
% ~~~ INPUTS ~~~
% 
% Optional name-value pair input arguments: 
% 	 'timeVector': A vector of timepoints to set up the model for. Default is 1:100. 
%    'paramInfo': A param_info struct that is used to set parameter values in the model object. 
%         This struct has fields: 
%        'paramNames':         A string corresponding to the parameter name, or 
%                                  a 1 by 2 cell array of strings of forward and backward 
%                                  reaction rate parameters for reversible reactions. 
%        'paramVals':          A numerical value of the parameter, or a 1 by 2 vector 
%                                  for reversible reactions. 
%        'reactionString':     This is either the reaction string or the  
%                                  string 'global'. If the paramNames is a 1 x 2 cell array
%                                  and the paramVals is a 
%                                  then the reaction is reversible, and if the paramNames is
%                                  a string, then it is a irreversible reaction. When this 
%                                  argument is 'global', we always have a string and a scalar 
%                                  for paramNames and paramVals respectively. 
%        'reactionString':     this is either the reaction string or the  
%                                  string 'global'. 
%        'paramRanges':        Non-negative orthant values for the parameter upper and 
%                                  lower bounds. For each parameter this is either a 
%                                  2 by 1 vector (for irreversible reactions) or a 2 by 2
%                              matrix, where the first and second rows are the upper 
%                                  and lower bounds respectively, and the first and second 
%                                  columns are the bounds for the forward and reverse rate
%                                  parameters respectively. If nothing is specified, then 
%                                  the bounds for a parameter with value VAL is [VAL/10; VAL*10]
%        'paramNotes'          Human readable notes. 
%                  
% 
% ~~~ OUTPUTS ~~~
%
% m: a Simbiology Model Object
%
% emo: An exported Simbiology Model Object. 
%
% mi: An updated mcmc_info struct. The following fields are added to the mcmc_info struct 
%   names_ord: An ordered list of species and parameters to be estimated.  
%   emo: Exported model Object. 

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

p = inputParser;
p.addParameter('paramInfo', [], @isstruct);
p.addParameter('timeVector', 1:100, @isnumeric);
p.parse(varargin{:});
p = p.Results;
tv = p.timeVector;


%% ################## EDIT THIS SECTION AS NEEDED ####################

% Specify the parameter config files, dna string and initial concentration,
% inducers etc. See the TXTL Modeling toolbox documentation for more information. 

% setup model object (this model uses the txtl toolbox at 
tube1 = txtl_extract('Emcmc2018');
tube2 = txtl_buffer('Emcmc2018');
tube3 = txtl_newtube('pLacdeGFP');
txtl_add_dna(tube3, ...
  'plac(50)', 'utr1(20)', 'deGFP(1000)', 0, 'plasmid');	
m = txtl_combine([tube1, tube2, tube3]);
m.UserData.energymode = 'regeneration';
%% MOST LIKELY YOU WILL NOT NEED TO TOUCH THIS SECTION
% See: 
% https://www.mathworks.com/help/simbio/ug/selecting-absolute-
% tolerance-and-relative-tolerance-for-simulation.html
cs1 = getconfigset(m);
set(cs1.RuntimeOptions, 'StatesToLog', 'all');
set(cs1.SolverOptions, 'AbsoluteToleranceScaling', 1);
set(cs1.SolverOptions, 'AbsoluteTolerance', 1.0e-6);
set(cs1.SolverOptions, 'AbsoluteToleranceStepSize', tv(end)*1.0e-6*0.1);
set(cs1.SolverOptions, 'RelativeTolerance', 1.0e-6);
% try: AbsoluteToleranceStepSize = AbsoluteTolerance * StopTime * 0.1

tic
    [~] = txtl_runsim(m,tv(end));
toc



%
%% DONT TOUCH THIS SECTION
% modify the base model to get it ready for parameter estimation. 
% Change the scope of the reversible reaction parameters from the
% reaction to the model level scope. This is needed for adding rules
% between the parameters (like Kd rules).
% The rules are needed because in parameter esitmation when you change some
% base parameter (like Kd) the parameters that are associated should also
% change. Just like for the reversible reactions, the irreversible 
% reaction params need to be model scoped because
% for the elongation rate the parameter is tied to NTP consumption rates

globalize_params(m)

%% Optionally set parameters for in the model object
if ~isempty(p.paramInfo)
    pinf = p.paramInfo;
    % set the parameters
    for i = 1:length(pinf, 1) 
        setparam(m, pinf.reactionString(i), pinf.paramNames(i), pinf.paramVals(i));
    end
end


end

