function [out] = subtractCols(m, orig, sub, floorzero)
%subtractCols subtract certain columns of data from others. Useful for
%subtracting negative control data. 
% m = matrix of data, in #time points vs # columns
% orig, n element row vector of the indices of the original cols to be 
% subtracted form
% sub, n element row vector of columns indices for columns to stract from
% the original columns
%
% output = m(:,orig)-m(:,sub);


m(:,orig) = m(:,orig)-m(:,sub);
out = m;
if floorzero
        out(out<0)=0;
end



end

