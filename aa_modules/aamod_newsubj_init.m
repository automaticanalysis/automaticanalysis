% AA module  - intialise a subject by making directory
% [aap,resp]=aamod_newsubj_init(aap,task,i)
% Initialises a subject
% Rhodri Cusack MRC CBU Cambridge Aug 2004
% djm, 170708: add directory for figures, diagnostic images etc
% @@@ DO WE NEED THIS IN AA4? @@@

function [aap,resp]=aamod_newsubj_init(aap,task,i)

if (ischar(aap))
    load(aap);
end;
if (ischar(i))
    i=str2num(i);
end;

resp='';

switch task
    case 'summary'
        subjpath=aas_getsubjpath(aap,i);
        resp=sprintf('Directory made for subject %s\n',subjpath);
        
    case 'report'
    case 'doit'
        if (strcmp(aap.directory_conventions.remotefilesystem,'none'))
%             % Make new subject directory
%             subjpath=aas_getsubjpath(aap,i);
%             if (~exist(subjpath,'dir'))
%                 aas_log(aap,0,sprintf('New directory for %s',subjpath));
%                 [s,w]=aas_shell(sprintf('mkdir %s',subjpath));
%                 if (s)
%                     aas_log(aap,1,w);
%                 end;
%             end;
            
%             %%% djm: create subdirectories here, so files can be added to them
%             %%% as necessary before starting main tasks.
%             % add subdirectory for figures, diagnostic images etc
%             if (~exist(fullfile(subjpath,'figures'),'dir'))
%                 [s,w]=aas_shell(sprintf('mkdir %s',fullfile(subjpath,'figures')));
%                 if (s); aas_log(aap,1,w); end;
%             end;
%             
%             % add subdirectory for event/epoch specifications etc
%             if (~exist(fullfile(subjpath,'events'),'dir'))
%                 [s,w]=aas_shell(sprintf('mkdir %s',fullfile(subjpath,'events')));
%                 if (s); aas_log(aap,1,w); end;
%             end
%             
%             % create/move to block directories
%             for b=aap.acq_details.selected_sessions
%                 blockdir=fullfile(subjpath,aap.acq_details.sessions(b).name);
%                 if ~exist(blockdir,'dir');
%                     [s,w]=aas_shell(sprintf('mkdir %s',blockdir));
%                     if (s); aas_log(aap,1,w); end;
%                 end
%                 cd(blockdir)
%                 
%                 % add subdirectory for figures, diagnostic images etc
%                 if (~exist(fullfile(blockdir,'figures'),'dir'))
%                     [s,w]=aas_shell(sprintf('mkdir %s',fullfile(blockdir,'figures')));
%                     if (s); aas_log(aap,1,w); end;
%                 end;
%                 
%                 % add subdirectory for event/epoch specifications etc
%                 if (~exist(fullfile(blockdir,'events'),'dir'))
%                     [s,w]=aas_shell(sprintf('mkdir %s',fullfile(blockdir,'events')));
%                     if (s); aas_log(aap,1,w); end;
%                 end;
%             end;
        end
        %%%
        

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;


