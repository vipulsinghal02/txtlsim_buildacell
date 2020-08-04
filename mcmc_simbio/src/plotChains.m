function plotChains(m, nW, legends, varargin)
% Plot markov chain trajectories. 
% m is a nParam x nWalkers x numSamples
% nW is the number of walkers to plot. 
% legs is legends
p = inputParser;
p.addParameter('Visible', 'on', @ischar)
p.parse(varargin{:});
p=p.Results;

[nParam, nWalkers, nSamples] = size(m);
[n1 n2] = twofactors(nParam);

if isprime(nParam)
    if nParam<12
        n1 = ceil(nParam/3); % 5 columns
        n2 = 3;
    elseif nParam < 20
        n1 = ceil(nParam/4); % 5 columns
        n2 = 4;
    else
        n1 = ceil(nParam/5); % 5 columns
        n2 = 5;
    end
    
end

wix = unique(ceil(rand(nW, 1)*nWalkers));
m = m(:,wix, :);
figure('Visible', p.Visible)
set(0,'Units','normalized')

set(gcf,'Units', 'normalized')
set(gcf, 'Position', [0.05, 0.1, 0.9, 0.85])
for i = 1:nParam
    subplot(n1, n2, i)

    for j = 1:length(wix)
        plot1 = plot(1:nSamples, squeeze(m(i, j, :)),...
            'LineWidth', 0.1,...
            'color', [0.2 0.7 0.1].^2);
        plot1.Color(4) = 0.015;
        hold on       
    end
    title(legends{i}, 'fontsize', 14);
    xlabel('Iteration', 'fontsize', 16);
    ylabel('log-value', 'fontsize', 16);
    ax = gca;
    ax.FontSize = 14;
    
    
end

end

