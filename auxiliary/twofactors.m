function [num1, num2] = twofactors(n)
%twofactors decompose a number into approximately equal factors
%   if prime, just return 1 and number itself 
% Vipul singhal
pf = factor(n);
numpf = length(pf);
mid = ceil(numpf/2);
num1 = prod(pf(1:mid));
num2 = prod(pf(mid+1:end));
end

