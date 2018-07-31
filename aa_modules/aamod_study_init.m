% AA module - makes study directory
% Rhodri Cusack MRC CBU Cambridge Aug 2004 
% Initialises a study

function [aap,resp]=aamod_study_init(aap,task)

if (ischar(aap))
    load(aap);
end;

resp='';

switch task
    case 'doit'
        if (strcmp(aap.directory_conventions.remotefilesystem,'none'))
             studydir=aap.acq_details.root;
             d=dir(studydir);
             if isempty(d) 
                 cmd=['mkdir ' studydir];
                 [s w]=aas_shell(cmd);
                 if (s~=0)
                     aas_log(aap,1,w);
                 end;    
             elseif (length(d)<2)
                 aas_log(aap,1,sprintf('There is a file with the same name as the desired study directory: %s',studydir));
             end;
% 
%             % djm: add subdirectory for event/epoch specifications etc
%              if (~exist(fullfile(studydir,'events'),'dir'))
%                  [s,w]=aas_shell(sprintf('mkdir %s',fullfile(studydir,'events')));
%                  if (s); aas_log(aap,1,w); end;
%              end
        end;
end;




