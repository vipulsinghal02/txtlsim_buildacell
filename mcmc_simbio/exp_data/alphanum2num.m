function [numix] = alphanum2num(startIx, endIx, cix)
%alphanum2num Convert the alphanumeric well ID to a numeric column number of data array   
%   A little background: when I copy the data from an excel file to a
%   MATLAB .m file, I skip the first two columns, which are the time and
%   temperature columns. After that starts the data. This function assumes
%   a single continuous block of data. 
% - startIx is the alphanumeric index of the top left well. 
% - endIx is the alphanumeric index of the bottom right well. 
% - cix is either a string of alphanumeric well ID that gets converted or it
%   is a cell array of strings, in which case the output numix is a matrix array of
%   the same dimensions as the cell array. 
%
% EXAMPLES
% alphanum2num('C8', 'H15', 'E11')
% ans =
%     20 
%
% alphanum2num('C8', 'H15', {'E11', 'E12', 'E13'; 'G13', 'G14', 'G15'})
% ans =
%     20    21    22
%     38    39    40
    

%check input data correctness, and parse out information from it
[salpha, snum, sa, sn] = parse_alphanum(startIx);
[ealpha, enum, ea, en] = parse_alphanum(endIx);

charformat = false;
cellformat = false;

if ischar(cix)
    charformat = true;
    [calpha, cnum, ca, cn] = parse_alphanum(cix);
elseif iscell(cix)
    cellformat = true;
    calpha = cell(size(cix));
    cnum = zeros(size(cix));
    ca = cell(size(cix));
    cn = cell(size(cix));
    for i = 1:size(cix, 1)
        for j = 1:size(cix, 2)
            [calpha{i,j} , cnum(i,j), ca{i,j}, cn{i, j}] = parse_alphanum(cix{i,j});
        end
    end
else
    error('The indices to convert must either be a single alphanumeric well ID or a cell array of alphanumeric well IDs')
end

% used the parsed information to compute the indices. Essentially we
% compute the index in a new coordinate system. 
letters = 'ABCDEFGHIJKLMNOP';
if charformat
    % COMPUTEIX FUNCTION
    numix = computeix(letters, calpha, cnum, salpha, ealpha, snum, enum);
elseif cellformat
    numix = zeros(size(cix));
    for i = 1:size(cix, 1)
        for j = 1:size(cix, 2)
            numix(i,j) = computeix(letters, calpha{i,j}, cnum(i,j), salpha,...
                ealpha, snum, enum);
        end
    end
end
end


function [alpha, num, a, n] = parse_alphanum(an)
% separate out the alpha and the numeric parts of the data
a = regexp(an, '[A-P]');
n = regexp(an, '[0-9]');
if length(an)<2 || length(an)>3 || isempty(a) || isempty(n)
    % this is not exactly the right check, but whatever. 
    error('Input must be in the form X##, where X is a letter from A to P, and ## is a number from 1 to 24')
   
end
alpha = an(a);
num = str2double(an(n));  
end

function ix = computeix(letters, calpha, cnum, salpha, ealpha, snum, enum)
    ocrIX = find(letters==calpha); %old coordinate row index of well to convert
    occIX = cnum; %old coordinate column index of well to convert
    
    ncsr = find(letters==salpha); % new coordinate system start row
    ncer = find(letters==ealpha); % new coordinate system end row. don't actually need this to transform coordinates
    ncsc = snum; %new coordinate start column
    ncec = enum; % new coordinate end column
    
    %new coordinate system width
    ncw = ncec - ncsc+1;
    ncr = ocrIX -ncsr+1; %row index of the input in the new coordinate system 
    ncc = occIX - ncsc + 1; %column index of the input in the new coordinate system
    
    % now compute the overall index of the input in the new coordinate
    % system assuming that we count row wise. 
    ix = (ncr-1)*ncw + ncc;
end
