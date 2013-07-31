% automatic analysis - log tool
% used to display messages
% if iserr is true, an error is generated
% function aas_log(aap,iserr,msg,style)

function aas_log(aap,iserr,msg,style)

try
    % don't attempt html tags if running outside of Matlab desktop
    if (~exist('style','var') || ~aap.gui_controls.usecolouroutput)
        style='text';
    end;
catch
    style='text';
end;

if (iserr)
    try
        logitem(aap,'\n\n**** automatic analysis failed - see reason and line numbers below\n','red');
        
        if ~isempty(aap.options.email)
            % In case the server is broken...
            try
                aas_finishedMail(aap.options.email, aap.acq_details.root, msg)
            catch
            end
        end
    catch
        error('internal aa error fatal -1');
    end;
end;

try
    verbose=aap.options.verbose || iserr;
catch
    verbose=1;
end;

if (verbose)
    try
        logitem(aap,[msg '\n'],style);
    catch
        error('internal aa error fatal -2');
    end;
end;
if (iserr)
    global aaworker
    if (isfield(aaworker,'errorflagname'))
        fid=fopen(aaworker.errorflagname,'w');
        fprintf(fid,'%d:%d:%d %d:%02d:%02d  ',round(clock));
        fclose(fid);
    end;
    try
        cd(aap.internal.pwd)
    catch
    end;
    disp('for help, see the <a href="http://www.cambridgeneuroimaging.com/aawiki">aa wiki</a>')
    error('aa error');
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function logitem(aap,msg,style)
global aaparallel;

if (~exist('style','var'))
    style='text';
end;

% output to screen
if (~isempty(strfind(msg,'parallel-manage-workers')))
    pmw_log.msg=msg;
    pmw_log.when=clock;
    if (~isfield(aaparallel,'manage_workers_log'))
        aaparallel.manage_workers_log=pmw_log;
    else
        aaparallel.manage_workers_log(end+1)=pmw_log;
    end;
else
    try
        % no colours when deployed... (or when on command line)
        if isdeployed || ~usejava('desktop')
            fprintf(msg);
        else
            cprintf(style,msg);
        end;
    catch
        fprintf(msg);
    end;
end;
% and to worker file?
global aaworker
if (isfield(aaworker,'logname'))
    fid=fopen(aaworker.logname,'a');
    fprintf(fid,'%d:%d:%d %d:%02d:%02d  ',round(clock));
    fprintf(fid,'%s',msg);
    fclose(fid);
    aas_propagateto(aaworker.master.hostname,aaworker.logname);
end;
% and to sdb
%   crypt=aaworker.aacc.crypt;



% No longer-this is now done in the python wrapper
% if (strcmp(aap.options.wheretoprocess,'aws'))
%     try
%         timenow=now;
%         randnum=round(rand()*1e6);
%         timestr=datestr(timenow,30);
%         itemname=sprintf('%s_%06d',timestr,randnum);
%         li=[];
%         try
%             li.workerid=num2str(aaparallel.processkey);
%             % This is time zone independent...
%             li.utctime=sprintf('%12.3f',utc_time());
%             li.analysisid=aap.directory_conventions.analysisid;
%             if (isnumeric(style))
%                 style=dec2hex(style*128);
%                 style=['#' style(:)'];
%             end;
%             li.style=style;
%             li.msg=crypt.tobase64(crypt.encrypt(int8(msg)));
%             li.msg(li.msg==10)=[];
%             li.source='matlab';
%         catch
%         end;
%     catch
%     end;
%     sdb_put_attributes(aap,aaworker.logqname,itemname,li);
%
% end;
end
