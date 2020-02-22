function [f, meanconc, stdconc , idxnotused] = ...
    plotEstimTraces_ver201711(m,em,ts, datmat, ds, ms, varargin)
%plotEstimTraces_ver201711 based on the function plotEstimTraces, except
%measured species is as used in the mcmc toolbox. Ie, ms is a cell array of
%cell arrays. 
% 
% m is the MCMC data, either in 3D or 2D
% ts = tspan
%   datmat = data matrix
% em = exported model object
% ds = dosing strat
% ms = measured species



p = inputParser;
p.addParameter('title',[] ,@ischar)
p.addParameter('nplot',150,@isnumeric);
p.addParameter('colorseq', [7, 4, 6, 1, 8, 9, 10],...
    @isnumeric);...
    % colorscheme: meangirls (RNA), rainforest (GFP), swamplands (CFP), 
% bog (YFP), ...
p.addParameter('Visible', 'on', @ischar)
p.parse(varargin{:});
p=p.Results;


if ndims(m) == 3
    m = m(:,:)';
end
   
idx = randperm(size(m, 1), p.nplot);


% extract dose matrix and dose names from dosing strat
nICs = length(ds(1).concentrations);
nDSP = length(ds);
dn = cell(length(ds));
for i = 1:length(ds)
    dn{i} = ds(i).species;
end
dose = zeros(nICs , nDSP);
for j = 1:nICs
    for i = 1:nDSP
        dose(j,i) = ds(i).concentrations(j);
    end
end

nts = length(ts);
nMS = length(ms);
% Compute sample trajectories, max values for axis limits, means, 
% standard deviations 
meanconc = zeros(nts, nMS, nICs);
stdconc = zeros(nts, nMS, nICs);
samp = zeros(nts, nMS, p.nplot, nICs);

idxnotused = zeros(p.nplot, nICs);
for i = 1:nICs
    for kk=1:p.nplot
        try
        sd = simulate(em, [exp(m(idx(kk),:)'); dose(i,:)']);
        sd = resample(sd, ts);
        
        % since measured species is a cell array of cell arrays of strings,
        % we need to loop over the outer cell, summing the values of the
        % species in the inner cell. 
        for ss = 1:nMS
                [~, XX] = selectbyname(sd, ms{ss}); 
                X = sum(XX, 2);
                % output the relevant data to the samples 
                samp(:,ss,kk, i) = X;
        
        end
        
        
        
        
        
        catch ME
            ME.message
            idx(kk);
            idxnotused(kk, i) = idx(kk);
        end
    end
    % compute mean and std over the p.nplot dimension, 
    % for each timepoint, IC and
    % each species
    kknotused = find(idxnotused(:,i));
    for j = 1:nMS
        meanconc(:,j, i) = ...
            mean(samp(:,j,setdiff(1:(p.nplot),kknotused), i), 3);
        stdconc(:, j, i) = ...
            std(samp(:,j,setdiff(1:(p.nplot), kknotused), i),0, 3);
    end
end

%%

%

cc = colorschemes;
f = figure('Visible', p.Visible);
% ss = get(0, 'screensize');
% set(gcf, 'Position', [50 100 ss(3)/1.1 ss(4)/1.3]);
% Compute species maxes for plotting
mxtemp = max(max(meanconc + stdconc, [], 1), [],3);
mxtemp = max(mxtemp, max(datmat)); % dm is the data matrix, 
% DIM1:time, DIM2:measured species
maxsp = mxtemp';

h = zeros(nICs, nMS);
ptch = zeros(nICs, nMS);
d = zeros(nICs, nMS);

csq = p.colorseq;

for i = 1:nICs
    for j = 1:nMS
    subplot(nICs, nMS,nMS*(i-1)+j);
    [h(i, j), ptch(i, j)] = boundedline(ts/3600, meanconc(:,j, i),...
        stdconc(:, j, i));
    set(ptch(i, j), 'FaceColor', cc{2,csq(j)}(1,:), 'FaceAlpha', 0.5);
    set(h(i, j), 'Color', cc{2,csq(j)}(2,:), 'LineStyle', '--');
    hold on 
    set(h(i, j), 'LineWidth', 2)
    d(i, j)=plot(ts/3600,datmat((i-1)*nts + 1 : (i)*nts,j),...
        'color',cc{2,csq(j)}(3,:) ,'linewidth',2);
    hold on 
    set(gca, 'Ylim', [0, round(maxsp(j)*1.1, -(order(maxsp(j)*1.1)-2))])
    xlabel('time/h')
    ylabel('conc/nM')
    
%     title(sprintf('%s conc, dosed %s = %0.2g nM', mn{j}, dn{1}, dose(i,1))) 
    % right now I can only support a single species dose. Need to come up
    % with an elegant way of putting all the dosing information in the
    % title or in floating text. 
    if i ==1 && j ==1
        ax = gca;
    end
    end
end
% if ~isempty(p.title)
% suptitle(p.title)
% end

% handles = [h(1, :), d(1,:)];
% lg1 = cell(1, nMS);
% lg2 = cell(1, nMS);
% for i = 1:nMS
%     lg1{i} = [mn{i} ' sample mean'];
%     lg2{i} = [mn{i} ' exp data'];
% end

% legstr = [lg1, lg2];
%  legend(ax, handles,legstr, 'Location', 'NorthWest');

end


