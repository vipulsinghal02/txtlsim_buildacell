% induction.m - gene expression with an inducer
% R. M. Murray, 11 Sep 2012
%
% This file contains a simple example of setting up a TXTL simulation
% for a negatively autoregulated gene.  The constants for this example
% come from the Simbiology toolbox example page:
%
%    http://www.mathworks.com/help/toolbox/simbio/gs/fp58748.html
%
close all
clear all


% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
tube1 = txtl_extract('E30VNPRL');
tube2 = txtl_buffer('E30VNPRL');

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('circuit');

% Define the DNA strands (defines TX-TL species + reactions)
dna_tetR = txtl_add_dna(tube3, ...
  'p70(50)', 'utr1(20)', 'tetR(647)', 15, 'linear');
dna_deGFP = txtl_add_dna(tube3, ...
  'ptet(50)', 'utr1(20)', 'deGFP(1000)', 5, 'linear');


%
% Next we have to set up the reactions that describe how the circuit
% works.  Transcription and translation are already included above, so
% we just need to include protein-protein and protein-DNA interactions.
%
% Note that the commands in this section are standard Simbiology commands,
% so you can put anything you want here.
%

% No additional reactions required for this circuit
% tetR-DNA interactions are automatically included in tetR setup

%
% Describe the actual experiment that we want to run.  This includes 
% combining the various tubes and also adding any additional inducers
% or purified proteins that you want to include in the run.
%

% Put in gamS to protect linear DNA
gamS = txtl_addspecies(tube3, 'protein gamS', 100);

% Set up the plot
figure(2); clf();
count = 1;

% Do runs at different inducer levels, linearly spaced
levels = [0 2 5 10 20 40 60 80 100];
maxGFP = zeros(1, length(levels));
colors = {'r', 'b', 'g', 'c', 'm', 'y', 'k', 'r--', 'b--'};
% Mix the contents of the individual tubes
  Mobj = txtl_combine([tube1, tube2, tube3]);

for atc = levels 
  
  % Run a simulation
  configsetObj = getconfigset(Mobj, 'active');
  set(configsetObj, 'StopTime', 14*60*60);
  %set(configsetObj, 'SolverType', 'ode23'); 
  [t_ode{count}, x_ode{count}, mObj, simData] = txtl_runsim(Mobj, configsetObj);
  

  
  % Add additional inducer for the next run
  if count < size(levels,2)
  inducer = txtl_addspecies(Mobj, 'aTc', levels(count+1)-levels(count));
  count = count + 1;
  end
end

for count = 1:9
  % Plot the time trace  
  figure(2); hold on;
  iGFP = findspecies(Mobj, 'protein tetR');
  plot(t_ode{count}, x_ode{count}(:, iGFP), colors{count});
  labels{count} = [int2str(levels(count)) ' nM aTc'];
  title('protein tetR')
  % Keep track of the max expression for later plotting
  maxGFP(count) = max(x_ode{count}(:, iGFP));
end
%title('Time Responses');
lgh = legend(labels, 'Location', 'Northwest');
legend(gca, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');
for count = 1:9
  % Plot the time trace  
  figure(3); hold on;
  iGFP = findspecies(Mobj, 'protein deGFP*');
  plot(t_ode{count}, x_ode{count}(:, iGFP), colors{count});
  labels{count} = [int2str(levels(count)) ' nM aTc'];
  title('protein deGFP*')
  % Keep track of the max expression for later plotting
  maxGFP(count) = max(x_ode{count}(:, iGFP));
end
%title('Time Responses');
lgh = legend(labels, 'Location', 'Northwest');
legend(gca, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');
for count = 1:9
  % Plot the time trace  
  figure(4); hold on;
  iGFP = findspecies(Mobj, 'aTc');
  plot(t_ode{count}, x_ode{count}(:, iGFP), colors{count});
  labels{count} = [int2str(levels(count)) ' nM aTc'];
  title('aTc')
  % Keep track of the max expression for later plotting
  maxGFP(count) = max(x_ode{count}(:, iGFP));
end
%title('Time Responses');
lgh = legend(labels, 'Location', 'Northwest');
legend(gca, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');

for count = 1:9
  % Plot the time trace  
  figure(5); hold on;
  iGFP = findspecies(Mobj, 'DNA ptet--utr1--deGFP');
  plot(t_ode{count}, x_ode{count}(:, iGFP), colors{count});
  labels{count} = [int2str(levels(count)) ' nM aTc'];
  title('DNA ptet--utr1--deGFP')
  % Keep track of the max expression for later plotting
  maxGFP(count) = max(x_ode{count}(:, iGFP));
end
%title('Time Responses');
lgh = legend(labels, 'Location', 'Northwest');
legend(gca, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');

for count = 1:9
  % Plot the time trace  
  figure(6); hold on;
  iGFP = findspecies(Mobj, '2 aTc:protein tetRdimer');
  plot(t_ode{count}, x_ode{count}(:, iGFP), colors{count});
  labels{count} = [int2str(levels(count)) ' nM aTc'];
  title('2 aTc:protein tetRdimer')
  % Keep track of the max expression for later plotting
  maxGFP(count) = max(x_ode{count}(:, iGFP));
end
%title('Time Responses');
lgh = legend(labels, 'Location', 'Northwest');
legend(gca, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');

for count = 1:9
  % Plot the time trace  
  figure(7); hold on;
  iGFP = findspecies(Mobj, 'DNA ptet--utr1--deGFP:protein tetRdimer');
  plot(t_ode{count}, x_ode{count}(:, iGFP), colors{count});
  labels{count} = [int2str(levels(count)) ' nM aTc'];
  title('DNA ptet--utr1--deGFP:protein tetRdimer')
  % Keep track of the max expression for later plotting
  maxGFP(count) = max(x_ode{count}(:, iGFP));
end
%title('Time Responses');
lgh = legend(labels, 'Location', 'Northwest');
legend(gca, 'boxoff');
ylabel('Species amounts [nM]');
xlabel('Time [min]');



% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
