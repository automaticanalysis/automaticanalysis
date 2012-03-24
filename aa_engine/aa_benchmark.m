% Automatic analysis - display benchmark
% After some or all of an analysis is run, this command will dump the time
% taken to perform each of the stages
% Rhodri Cusack MRC CBU Cambridge 2004

function aa_benchmark(studyroot)

if (~exist('studyroot','var'))
    studyroot=pwd;
end;

% First, load AAP structure
aaploadfn=fullfile(studyroot,'aap_parameters');
if exist(aaploadfn,'file'); % djm: prompt if necessary
    aaploadfn=spm_select(1,'mat','Please select aap_parameters file:',pwd);
end
load(aaploadfn);

[s w]=system('hostname');
aas_log(aap,0,['BENCHMARK TIMINGS IN SECONDS ' datestr(now) ' MACHINE ' deblank(w)]);
aas_log(aap,0,['=============================================================']);

% THE MODULES IN AAP.TASKLIST.STAGES ARE RUN IF A CORRESPONDING DONE_ FLAG
% IS NOT FOUND. ONE IS CREATED AFTER SUCCESSFUL EXECUTION
% Now run stage-by-stage tasks
aap.internal.total=0;
aap.internal.stagesnotdone=0;
for k=1:length(aap.tasklist.stages)
    doneflag=fullfile(aap.acq_details.root,['done_' aap.tasklist.stages{k}]);

    % find out whether this module needs to be executed once per study, subject or session
    [aap,domain]=feval(aap.tasklist.stages{k},aap,'domain');


    switch (domain)
        case 'study'
            doneflag=fullfile(aap.acq_details.root,['done_' aap.tasklist.stages{k}]);
            aap=processflag(aap, doneflag);
        case 'subject'
            for i=1:length(aap.acq_details.subjects)
                doneflag=fullfile(aas_getsubjpath(aap,i),['done_' aap.tasklist.stages{k}]);
                aap=processflag(aap, doneflag);
            end;
        case 'session'

            for i=1:length(aap.acq_details.subjects)
                for j=aap.acq_details.selected_sessions
                    doneflag=fullfile(aas_getsesspath(aap,i,j),['done_' aap.tasklist.stages{k}]);
                    aap=processflag(aap, doneflag);
                end;
            end;
        otherwise
            aas_log(aap,1,sprintf('Unknown domain %s associated with stage %s',aap.tasklist.domain{k},aap.tasklist.stages{k}));
    end;
end;


fprintf('STAGES NOT COMPLETED: %d  TOTAL TIME SO FAR: %f\n',aap.internal.stagesnotdone,aap.internal.total);
return;

function [aap]=processflag(aap, doneflag)
fid=fopen(doneflag,'r');
if (fid==-1)
    aap.internal.stagesnotdone=aap.internal.stagesnotdone+1;
    fprintf('Stage %s not completed?\n',doneflag);
else
    tme=fscanf(fid,'%f');
    aap.internal.total=aap.internal.total+tme;
    fprintf('%f\ttaken to complete stage %s\n',tme,doneflag);
end;
return;