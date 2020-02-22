function [llh] = loglike_gaussian(logp, stdev, ICarray,measuredspidx, tspan, ...
	dataArray, genemodel, lognormvec)
% loglike_gaussian Gaussian log likelihood
%   Compute the log likelihood of the data given the model parameters and
%   the standard deviation. 

% Data is input as a array that is  N x M x nIC. N is number of time points in tspan, 
% M is the number of species species

mV = mean(dataArray(:,measuredspidx,:),1);
meanVals = mean(mV,3); %meanVals is a row vec.

wt = sum(meanVals)./meanVals; %hight mean = lower wt
relWt = wt/sum(wt);
CONC_temp = zeros(size(dataArray));


nIC = size(ICarray,1); % number of different initial conditions. 
% initial conditions are row vectors. Rows are different sets of initial conditions.

for i = 1:nIC
[~,CONC_temp(:,:,i)] = genemodel(logp, ICarray(i,:), tspan);
end

residuals = dataArray(:,measuredspidx,:) - CONC_temp(:,measuredspidx,:);
resi = bsxfun(@times, residuals, relWt);

res = resi(:);
llh = sum(lognormvec(res, stdev));



end

