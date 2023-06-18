function subject_list = BIDS_subjectfilter(toplevelBIDSdir,filter_string)

% apply a simple filter to a BIDS participants.tsv file
% and return ID of subject(s) that satisify filter
%
% input
%
% toplevelBIDSdir - top level BIDS dir (must contain participants.tsv)
% filter_string - filter specification (see below)
%
% output
%
% subject_list - list (cell array) of subjects satisifying filter
%                (can be passed to aas_processBIDS or other modules)
%
% filter string syntax currently supports three operations:
%
%   keyword==value      example: 'sex=F','age==25' 
%   keyword>value       example: 'age<25', 'age < 25'
%   kewword<value       example: 'GVTD>0.001'
%
% note both '=' or '==' are accepted, and the filter may include spaces
%
% if a vector is passed for "value" when using '=', it defines
% a range to test. Example: 'GVTD=[0.1 0.5]' tests 0.1 < GVTD < 0.5
%

subject_list = {};

% read participants data

if ~exist(fullfile(toplevelBIDSdir,'participants.tsv'),'file')
    error('Cannot find participants.tsv file. Aborting...');
end

% parse filter_string

filter_string = strrep(filter_string,'==','='); % change == to =

filter_tokens = split(filter_string,{'>','<','='});

if (length(filter_tokens) ~=2)
    error('Bad filter string (format: "keyword >|<|= value"). Aborting...');
end

filter_keyword = filter_tokens{1};
filter_value = filter_tokens{2};

temp = split(filter_string,{ filter_keyword,filter_value });
if (length(temp) ~=3)
    error('Bad filter string (format: "keyword >|<|= value"). Aborting...');
end
filter_op = temp{2};

% strip leading and trailing whitespace (if any)
filter_keyword = strtrim(filter_keyword);
filter_value = strtrim(filter_value);

if (~isempty(str2num(filter_value)))
    filter_value = str2num(filter_value);
end

pdata = tsvread(fullfile(toplevelBIDSdir,'participants.tsv'));

% sanity check: the filter string keyword must be one of the column headers

pdata_header = pdata(1,:);

column_index = find(ismember(pdata_header,filter_keyword));

if (isempty(column_index))
    error('Filter keyword is not a column in participants.tsv. Aborting...');
end
    
% loop over subjects, add those passing

subject_index = []; % accumulate subjects passing filter

for sindex = 2:size(pdata,1)
    
    this_subject_value = pdata(sindex,column_index);
    if filter_test(this_subject_value,filter_op,filter_value)
        subject_index = [ subject_index sindex ];
    end
    
end
    
% subject_index = subject_index + 1; % first row is header...
subject_list = pdata(subject_index,1);

end

function pass = filter_test(value1,op,value2)

   pass = false;
   
   % apply the test value1-op-value2 (e.g. 1>3 or 'F'=='M')
   % value1 is from tsv and is a string (may need conversion to numeric)
   % value2 is parsed from filter and assumed converted to correct type
   
   switch op
       
       case '='
           
           % equals might be numeric or string equality or range test
           
           if ischar(value2)
               if strcmp(value1,value2); pass = true; end
           elseif ~isscalar(value2) % range test
               if (value2(1) < str2double(value1) && str2double(value1) < value2(2)); pass = true; end
           else
               if (str2double(value1)==value2); pass = true; end
           end
              
           
       case '>'
           
           if (str2double(value1) > value2)
               pass = true;
           end
           
       case '<'
           if (str2double(value1) < value2)
               pass = true;
           end
           
       otherwise
           
           error('Unimplemented filter op. Aborting...');
           
   end
   
   
   
end



function LIST = tsvread(fname)
    filestr = fileread(fname);
    nLines = sum(double(filestr)==10);
    LIST = textscan(filestr,'%s','Delimiter','\t');
    LIST = LIST{1};
    % protect against extraneous whitespace in the tsv:
    LIST = LIST(~cellfun(@isempty,LIST));
    LIST = reshape(LIST,[],nLines)';
end


