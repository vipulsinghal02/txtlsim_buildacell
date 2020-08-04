function  [tstamp, specificprojdir, st] = project_init(varargin)
% this function is to be called at the beginning of the project script stored in the projects
% directory of the MCMC toolbox.
% The script calling this should be a project file (.m) that sets up the data, models,
% MCMC simulation, runs it, and optionally plots the results.
% The main function of this file is that every time a given project is run, a directory
% created (using a timestamp) where the results of that run are stored.
%
% This function can be run without any arguments, in which case it will check if a directory
% with the same name as the project that calls it exists in the projects folder. If it does,
% then all the simulation data is stored in timestamped subdirectories within this directory,
% and this directory is what gets checked first for the various things like models,
% experimental data, est_info structs, model_map structs, etc. If any of these things are
% not present in this directory, then the models or data directories in the main directory
% are checked. If those directories also do not contain the files required, an error is thrown.
% If the directory does not exist, then a directory with the same name as the project file
% gets created, and all the models/ data etc are looked for in the main directories. This
% directory will be where the simulation data will be stored.
%
% on the other hand, if a string argument is supplied, then the above applied, but instead of
% the project name, we use this string as the name.
%
% the output of this function is a string that specifies the timestamp for the purposes of
% identifying the directory created.



% get the name of the (project) file calling this function, and its path.
st = dbstack(1,'-completenames');
if ispc
    
    fp = st(1).file;
    slashes = regexp(fp, '\');
    projdir = fp(1:slashes(end)-1); % full path to the directory containing the
    % function\script that called this function, without the final slash ...\dirname
    projname = st(1).name; % name of the function \ script that called this function
    tstamp = datestr(now, 'yyyymmdd_HHMMSS');
    
    % get the optional input (the name of the directory where things will be stored)
    p = inputParser;
    addParameter(p, 'proj', projname, @ischar); %
    addParameter(p, 'saveStr', tstamp, @ischar); %
    p.parse(varargin{:});
    p=p.Results;
    specificprojdir = [projdir '\' p.proj];
    % display some output text
    disp(sprintf('############################################ \n'))
    disp(sprintf('File and directory info:\n'))
    disp(sprintf('Project name: \n ''%s'' \n', projname))
    disp(sprintf('Directory where the project file is stored: \n ''%s'' \n', projdir))
    disp(sprintf('Directory where data will be stored: \n ''%s'' \n', specificprojdir))
    disp(sprintf('Timestamp for this run (yyyymmdd_HHMMSS): \n ''%s'' \n', tstamp))
    
    % check if the project directory already exists, and create it if it doesnt.
    if exist( specificprojdir ,'dir')
        disp(sprintf(['Project directory already exists, using this to store data\n' ...
            ' (in a subdirectory named ''%s''). \n'], ['simdata_' p.saveStr]));
        
        addpath(specificprojdir);
        addpath(genpath(specificprojdir));
        rmpath(genpath(specificprojdir));
        mkdir([specificprojdir '\simdata_' p.saveStr]);
        addpath([specificprojdir '\simdata_' p.saveStr]);
    else
        mkdir(specificprojdir);
        mkdir([specificprojdir '\simdata_' p.saveStr]);
        addpath([specificprojdir '\simdata_' p.saveStr]);
        disp(sprintf(['Project directory does not exist. Creating it,' ...
            ' and using this to store data\n' ...
            ' (in a subdirectory named ''%s''). \n'], ['simdata_' p.saveStr]));
        
    end
    disp(sprintf('############################################ \n'))
    
    
else
    
    fp = st(1).file;
    slashes = regexp(fp, '/');
    projdir = fp(1:slashes(end)-1); % full path to the directory containing the
    % function/script that called this function, without the final slash .../dirname
    projname = st(1).name; % name of the function / script that called this function
    tstamp = datestr(now, 'yyyymmdd_HHMMSS');
    
    % get the optional input (the name of the directory where things will be stored)
    p = inputParser;
    addParameter(p, 'proj', projname, @ischar); %
    addParameter(p, 'saveStr', tstamp, @ischar); %
    p.parse(varargin{:});
    p=p.Results;
    specificprojdir = [projdir '/' p.proj];
    % display some output text
    disp(sprintf('############################################ \n'))
    disp(sprintf('File and directory info:\n'))
    disp(sprintf('Project name: \n ''%s'' \n', projname))
    disp(sprintf('Directory where the project file is stored: \n ''%s'' \n', projdir))
    disp(sprintf('Directory where data will be stored: \n ''%s'' \n', specificprojdir))
    disp(sprintf('Timestamp for this run (yyyymmdd_HHMMSS): \n ''%s'' \n', tstamp))
    
    % check if the project directory already exists, and create it if it doesnt.
    if exist( specificprojdir ,'dir')
        disp(sprintf(['Project directory already exists, using this to store data\n' ...
            ' (in a subdirectory named ''%s''). \n'], ['simdata_' p.saveStr]));
        
        addpath(specificprojdir);
        addpath(genpath(specificprojdir));
        rmpath(genpath(specificprojdir));
        mkdir([specificprojdir '/simdata_' p.saveStr]);
        addpath([specificprojdir '/simdata_' p.saveStr]);
    else
        mkdir(specificprojdir);
        mkdir([specificprojdir '/simdata_' p.saveStr]);
        addpath([specificprojdir '/simdata_' p.saveStr]);
        disp(sprintf(['Project directory does not exist. Creating it,' ...
            ' and using this to store data\n' ...
            ' (in a subdirectory named ''%s''). \n'], ['simdata_' p.saveStr]));
        
    end
    disp(sprintf('############################################ \n'))
    
    
end

end

