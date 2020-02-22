function data_info = tierradataset(varargin)


% construct data info object:
    mn = {{'deGFP'}};
    dimlabels = ...
        {'time points', 'measured species',...
        'replicates', 'doses'};
    tv = [0:360:43200]';
    
    
    
if strcmp(varargin{1}, 'dataset11062018')
    dt = rawtierra2;
    eEC4_s70 = dt{2,1}(:, 28:39);
    eEC5_s70 = dt{2,2}(:, 28:39);
    eEC4_pTet = dt{2,1}(:, 16:27);
    eEC5_pTet = dt{2,2}(:, 16:27);
    eEC4_aTc = dt{2,1}(:, 1:15);
    eEC5_aTc = dt{2,2}(:, 1:15);
    eEC4_tetR = dt{2,1}(:,40:60); % tetR repression dataset. 
    eEC5_tetR = dt{2,2}(:,40:60); % tetR repression dataset. 
    
    
    % di element 7: the eEC4_tetR data
    dataInfo = 'Tierra Bio Data: eEC4_tetR';
    dv = [0.8, 0.4 0.2 0.1 0.05 0.025 0];
    nDoses=length(dv);
    currDataSet=eEC4_tetR;
    
    % set up the data for the different doses.
    da = zeros(121, 1, 3, nDoses);
    for i=1:3
        for j=1:nDoses
            da(:, 1, i, j) = currDataSet(:, (j-1)*3+i);
        end
    end
    
    di7 = struct('dataInfo', {dataInfo}, ...
        'timeVector', {tv}, ...
        'timeUnits', {'seconds'},...
        'dataArray', {da},...
        'measuredNames', {mn},...
        'dataUnits', {{'a.u.'}},...
        'dimensionLabels', {dimlabels}, ...
        'dosedNames', {{'tetR DNA'}},...
        'dosedVals', {dv}, ...
        'doseUnits', {'nM'});
    
    
    % di element 8: the eEC5_tetR data
    dataInfo = 'Tierra Bio Data: eEC5_tetR';
    dv = [0.8, 0.4 0.2 0.1 0.05 0.025 0];
    nDoses=length(dv);
    currDataSet=eEC5_tetR;
    
    % set up the data for the different doses.
    da = zeros(121, 1, 3, nDoses);
    for i=1:3
        for j=1:nDoses
            da(:, 1, i, j) = currDataSet(:, (j-1)*3+i);
        end
    end
    
    di8 = struct('dataInfo', {dataInfo}, ...
        'timeVector', {tv}, ...
        'timeUnits', {'seconds'},...
        'dataArray', {da},...
        'measuredNames', {mn},...
        'dataUnits', {{'a.u.'}},...
        'dimensionLabels', {dimlabels}, ...
        'dosedNames', {{'tetR'}},...
        'dosedVals', {dv}, ...
        'doseUnits', {'nM'});
    
else
    
    
    dt = rawtierra;
    eEC4_s70 = dt{2,1}(:, 28:39);
    eEC5_s70 = dt{2,2}(:, 28:39);
    eEC4_pTet = dt{2,1}(:, 16:27);
    eEC5_pTet = dt{2,2}(:, 16:27);
    eEC4_aTc = dt{2,1}(:, 1:15);
    eEC5_aTc = dt{2,2}(:, 1:15);
end

    
    % di element 1: the eEC4_s70 data
    dataInfo = 'Tierra Bio Data: eEC4_s70';
    dv = [1 2 4 8];
    nDoses=length(dv);
    currDataSet=eEC4_s70;
    
    % set up the data for the different doses.
    da = zeros(121, 1, 3, nDoses);
    for i=1:3
        for j=1:nDoses
            da(:, 1, i, j) = currDataSet(:, (j-1)*3+i);
        end
    end
    
    di1 = struct('dataInfo', {dataInfo}, ...
        'timeVector', {tv}, ...
        'timeUnits', {'seconds'},...
        'dataArray', {da},...
        'measuredNames', {mn},...
        'dataUnits', {{'a.u.'}},...
        'dimensionLabels', {dimlabels}, ...
        'dosedNames', {{'GFP DNA'}},...
        'dosedVals', {dv}, ...
        'doseUnits', {'nM'});
    
    % di element 2: the eEC5_s70 data
    dataInfo = 'Tierra Bio Data: eEC5_s70';
    dv = [1 2 4 8];
    nDoses=length(dv);
    currDataSet=eEC5_s70;
    
    % set up the data for the different doses.
    da = zeros(121, 1, 3, nDoses);
    for i=1:3
        for j=1:nDoses
            da(:, 1, i, j) = currDataSet(:, (j-1)*3+i);
        end
    end
    
    di2 = struct('dataInfo', {dataInfo}, ...
        'timeVector', {tv}, ...
        'timeUnits', {'seconds'},...
        'dataArray', {da},...
        'measuredNames', {mn},...
        'dataUnits', {{'a.u.'}},...
        'dimensionLabels', {dimlabels}, ...
        'dosedNames', {{'GFP DNA'}},...
        'dosedVals', {dv}, ...
        'doseUnits', {'nM'});
    
    % di element 3: the eEC4_pTet data
    dataInfo = 'Tierra Bio Data: eEC4_pTet';
    dv = [1 2 4 8];
    nDoses=length(dv);
    currDataSet=eEC4_pTet;
    
    % set up the data for the different doses.
    da = zeros(121, 1, 3, nDoses);
    for i=1:3
        for j=1:nDoses
            da(:, 1, i, j) = currDataSet(:, (j-1)*3+i);
        end
    end
    
    di3 = struct('dataInfo', {dataInfo}, ...
        'timeVector', {tv}, ...
        'timeUnits', {'seconds'},...
        'dataArray', {da},...
        'measuredNames', {mn},...
        'dataUnits', {{'a.u.'}},...
        'dimensionLabels', {dimlabels}, ...
        'dosedNames', {{'GFP DNA'}},...
        'dosedVals', {dv}, ...
        'doseUnits', {'nM'});
    
    
    % di element 4: the eEC5_pTet data
    dataInfo = 'Tierra Bio Data: eEC5_pTet';
    dv = [1 2 4 8];
    nDoses=length(dv);
    currDataSet=eEC5_pTet;
    
    % set up the data for the different doses.
    da = zeros(121, 1, 3, nDoses);
    for i=1:3
        for j=1:nDoses
            da(:, 1, i, j) = currDataSet(:, (j-1)*3+i);
        end
    end
    
    di4 = struct('dataInfo', {dataInfo}, ...
        'timeVector', {tv}, ...
        'timeUnits', {'seconds'},...
        'dataArray', {da},...
        'measuredNames', {mn},...
        'dataUnits', {{'a.u.'}},...
        'dimensionLabels', {dimlabels}, ...
        'dosedNames', {{'GFP DNA'}},...
        'dosedVals', {dv}, ...
        'doseUnits', {'nM'});
    
    
    % di element 5: the eEC4_aTc data
    dataInfo = 'Tierra Bio Data: eEC4_aTc';
    dv = [10000 1000 100 10 1];
    nDoses=length(dv);
    currDataSet=eEC4_aTc;
    
    % set up the data for the different doses.
    da = zeros(121, 1, 3, nDoses);
    for i=1:3
        for j=1:nDoses
            da(:, 1, i, j) = currDataSet(:, (j-1)*3+i);
        end
    end
    
    di5 = struct('dataInfo', {dataInfo}, ...
        'timeVector', {tv}, ...
        'timeUnits', {'seconds'},...
        'dataArray', {da},...
        'measuredNames', {mn},...
        'dataUnits', {{'a.u.'}},...
        'dimensionLabels', {dimlabels}, ...
        'dosedNames', {{'aTc'}},...
        'dosedVals', {dv}, ...
        'doseUnits', {'nM'});
    
    
    % di element 6: the eEC5_aTc data
    dataInfo = 'Tierra Bio Data: eEC4_aTc';
    dv = [10000 1000 100 10 1];
    nDoses=length(dv);
    currDataSet=eEC5_aTc;
    
    % set up the data for the different doses.
    da = zeros(121, 1, 3, nDoses);
    for i=1:3
        for j=1:nDoses
            da(:, 1, i, j) = currDataSet(:, (j-1)*3+i);
        end
    end
    
    di6 = struct('dataInfo', {dataInfo}, ...
        'timeVector', {tv}, ...
        'timeUnits', {'seconds'},...
        'dataArray', {da},...
        'measuredNames', {mn},...
        'dataUnits', {{'a.u.'}},...
        'dimensionLabels', {dimlabels}, ...
        'dosedNames', {{'aTc'}},...
        'dosedVals', {dv}, ...
        'doseUnits', {'nM'});
    
    
if strcmp(varargin{1}, 'dataset11062018')
    % added the tetR datasets
    data_info = [di1; di2; di3; di4; di5; di6; di7; di8];
    
else
    
    
    data_info = [di1; di2; di3; di4; di5; di6];
end
    
    
    
            
end
