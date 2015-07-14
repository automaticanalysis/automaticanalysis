function aa_close
% cleanup 

% restore path
global aacache;
if isstruct(aacache) && isfield(aacache,'bcp_path')
    path(aacache.bcp_path);
end

% kill running jobs
global taskqueue
if isstruct(taskqueue) && isfield(taskqueue,'scheduler')
    switch taskqueue.scheduler.Type
        case 'Torque'
            taskqueue.killall;
    end
end

clear global aa aacache aaparallel aaworker defaults taskqueue localtaskqueue;

end

