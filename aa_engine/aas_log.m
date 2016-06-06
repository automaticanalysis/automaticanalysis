% automatic analysis - log tool
% used to display messages
% if iserr is true, an error is generated
% function aas_log(aap,iserr,msg,style)

function aas_log(aap,iserr,msg,style)
global aa;
if ~isobject(aa) % aa is not running
    clear global aa;
    aa = aaClass('nopath','nogreet');
end

% figure out whether the caller is an engine (aaq_*)
ST = dbstack;
isEngine = numel(ST) > 1 && ~isempty(strfind(ST(2).name,'aaq'));

% aap defaults
if ~isstruct(aap) 
    aap.options.verbose = 2;
    aap.options.email = '';
	aap.gui_controls.usecolouroutput = false;
end
if ~isfield(aap.options,'verbose'), aap.options.verbose = 2; end

% don't attempt html tags if running outside of Matlab desktop
if nargin < 4 || ~aap.gui_controls.usecolouroutput
    style='text';
end;

if iserr % errors
    if aap.options.verbose > 0
        logitem('\n\n**** automatic analysis failed - see reason and line numbers below\n','red');
        
        % suppress e-mail from low level functions in case of cluster computing
        if ~isempty(aap.options.email) && (strcmp(aap.options.wheretoprocess,'localsingle') || isEngine)
            % In case the server is broken...
            try
                aas_finishedMail(aap.options.email, aap.acq_details.root, msg)
            catch
            end
        end
        logitem([msg '\n'],style);
        logitem('for help, see the ')
        logitem(sprintf('<a href="%s">aa wiki</a>\n',aa.aawiki),[0 0 1])
    end
    
    global aaworker
    if (isfield(aaworker,'errorflagname'))
        fid=fopen(aaworker.errorflagname,'w');
        fprintf(fid,'[%d.%02d.%02d %02d:%02d:%02d]  ',round(clock));
        fclose(fid);
    end
    if isfield(aap,'internal') && isfield(aap.internal,'pwd') && exist(aap.internal.pwd,'dir')
        cd(aap.internal.pwd)
    end
    
    if aap.options.verbose ~= -1, error(sprintf(['aa error:\n' msg '\n'])); end % undocumented, for devel only
else % warnings
    if aap.options.verbose == 2, logitem([msg '\n'],style); end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function logitem(msg,style)
global aaparallel;

if nargin < 2
    style='text';
end;

% output to screen
if ~isempty(strfind(msg,'parallel-manage-workers'))
    pmw_log.msg=msg;
    pmw_log.when=clock;
    if ~isfield(aaparallel,'manage_workers_log')
        aaparallel.manage_workers_log=pmw_log;
    else
        aaparallel.manage_workers_log(end+1)=pmw_log;
    end;
else
    % no colours when deployed... (or when on command line)
    if isdeployed || ~usejava('desktop') || ~exist('cprintf','file')
        fprintf(msg);
    else
        cprintf(style,msg);
    end;
end;
% and to worker file?
global aaworker
if isfield(aaworker,'logname')
    fid=fopen(aaworker.logname,'a');
    fprintf(fid,'[%d.%02d.%02d %02d:%02d:%02d]  ',round(clock));
    fprintf(fid,msg);
    fclose(fid);
    if isfield(aaworker,'master'), aas_propagateto(aaworker.master.hostname,aaworker.logname); end
end;
end
