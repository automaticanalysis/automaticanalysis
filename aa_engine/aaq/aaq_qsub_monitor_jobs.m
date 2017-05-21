function Nfinished = aaq_qsub_monitor_jobs(job_location, intervaltime, clearscreen, nloops)

if nargin == 0
    job_location = '/home/dp01/matlab/jobs/';
    intervaltime = 5;
    nloops = Inf;
    clearscreen = true;
end

if nargin < 2
   intervaltime = 5;
   nloops = Inf;
   clearscreen = true;
end

fclose('all');

status_complete = false;
mainloopind = 0;

while mainloopind < nloops
    C = {'' '' ''};
    mainloopind = mainloopind + 1;
    fclose('all');
    fid = 0;
    
    try
        T = splitstring(ls([job_location '*/*diary.txt']));
    catch %#ok<CTCH>
        disp('No Jobs Found!')
        pause(1) % allow time for ctrl+c to work
        continue
    end
    
    jobids = zeros(size(T));
    for ii = 1:length(T)
        [st, en] = regexp(T{ii}, 'Job[0-9]+');
        if ~isempty(st)
            jobids(ii) = str2double(T{ii}(st+3:en));
        end
    end
    [jobids, uind] = unique(jobids);
    T = T(uind);
    T(jobids == 0) = {};
    jobids(jobids == 0) = [];
    
    if length(fid) < max(jobids)
        fid(max(jobids)) = 0; %#ok<AGROW>
    end
    
    for ind = 1:length(jobids)
        ii = jobids(ind);
        jobpath = [fileparts(T{ind}) '/'];
        if exist([jobpath 'Task1.out.mat'], 'file')
            
            try
                load([jobpath 'Task1.out.mat'], 'finishtime','errormessage')
            catch %#ok<CTCH>
                % Best to restart the entire loop in case job folder was deleted
                disp('Trying again...')
                pause(1)
                continue
            end
            
            job_status = getlines([jobpath 'Task1.state.mat']);
            job_status = strtrim(job_status{end});
            C{ii,2} = job_status;
            if strcmp(job_status,'finished')
                status_complete(ii) = true;
            else
                status_complete(ii) = false;
            end
        end
        
        C{ii,1} = ['Job' num2str(ii)];
        
        if strcmp(C{ii,2},'running')
            if rem(ind, 5)
                tout = getlines([jobpath 'Task1.diary.txt']);
                C{ii,3} = strcat(tout{end});
            end
        elseif ~isempty(errormessage)
            C{ii,3} = errormessage;
            C{ii,2} = 'error';
            status_complete(ii) = 0;
        end
        
    end
    
    if clearscreen
        clc
    end
%     
%     status_complete(cellfun(@isempty, C(:,1), 'UniformOutput',true)) = 1;
%     if any(status_complete == 0)
%         disp('')
%         Cprint = fliplr(C(status_complete == 0,:)');
%         fprintf('%s \t %s \t %s \n',Cprint{:})
%         disp('')
%     else
%         disp('No jobs running')
%     end
    
    Nrunning = sum(strcmp(C(:,2),'running'));
    Nfinished = sum(strcmp(C(:,2),'finished'));
    Nerror = sum(strcmp(C(:,2),'error'));
    Npending = sum(strcmp(C(:,2),'pending'));
    disp('-----------------------------------------------')
    fprintf('Running: %d | Pending: %d | Finished: %d | Error: %d \n', Nrunning, Npending, Nfinished, Nerror)
    pause(intervaltime)
end

    function tout = getlines(fname)
        fid = fopen(fname, 'r');
        t = ''; tout = {''};
        while ischar(t)
            t = fgets(fid);
            if ischar(t)
                t = deblank(t);
                tout{end+1} = t; %#ok<AGROW>
            end
        end
        fclose(fid);
    end

    function rv = splitstring( str, varargin )
        %SPLITSTRING Split string into cell array
        %    ARRAY = SPLITSTRING( STR, DELIM, ALLOWEMPTYENTRIES ) splits the
        %    character string STR, using the delimiter DELIM (which must be a
        %    character array). ARRAY is a cell array containing the resulting
        %    strings. If DELIM is not specified, space delimiter is assumed (see
        %    ISSPACE documentation). ALLOWEMPTYENTRIES should be a logical single
        %    element, specifying weather empty elements should be included in the
        %    results. If not specified, the value of ALLOWEMPTYENTRIES is false.
        %
        %    Example:
        %         arr = splitstring( 'a,b,c,d', ',' )
        
        delim = '';
        AllowEmptyEntries = false;
        
        if numel(varargin) == 2
            delim = varargin{1};
            AllowEmptyEntries = varargin{2};
        elseif numel(varargin) == 1
            if islogical(varargin{1})
                AllowEmptyEntries = varargin{1};
            else
                delim = varargin{1};
            end
        end
        
        if isempty(delim)
            delim = ' ';
            ind2 = find( isspace( str ) );
        else
            ind2 = strfind( str, delim );
        end
        
        startpos = [1, ind2+length(delim)];
        endpos = [ind2-1, length(str)];
        
        rv = cell( 1, length(startpos) );
        for i=1:length(startpos)
            rv{i} = str(startpos(i):endpos(i));
        end
        
        if ~AllowEmptyEntries
            rv = rv( ~strcmp(rv,'') );
        end
    end

end

