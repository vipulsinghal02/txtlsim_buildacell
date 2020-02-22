function gd = dsg2014_d1
% dsg2014_d1 Dan / Zoltan's ACS paper data. 
% This file calls another file containing the raw data, and then organizes
% that into a groupedData format. 

[t, y, meta] = load_ACSDSG2014('RNAdeg');
Time = repmat(t*60, 9, 1);
RNA = reshape(y, numel(y), 1);
dosematrix = [[37.5, 75, 150, 200, 600, 700, 800, 900, 1000];
    zeros(length(t)-1, size(y, 2))]; 

Dose = reshape(dosematrix, numel(dosematrix), 1);
clear dosematrix
idmat = mtimes((1:9)', ones(1, length(t)));
ID = reshape(idmat', numel(idmat), 1);

dttable = table(ID, Time, RNA, Dose); 
gd = groupedData(dttable);

end

