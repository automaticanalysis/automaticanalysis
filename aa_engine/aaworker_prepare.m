function aaworker_prepare(aa_workerid,hostname)
global aaworker; 

% Set up aaworker self description, including path where control files will
% be stored

aaworker=[];
aaworker.aap=[];
aaworker.parmpath=aaworker_getparmpath(aaworker.aap,aa_workerid);

% Log output
aaworker.errorflagname=fullfile(aaworker.parmpath,'errorflag.txt');
aaworker.logname=fullfile(aaworker.parmpath,'log.txt');
aaworker.diaryname=fullfile(aaworker.parmpath,'diary.txt');
diary(aaworker.diaryname)
timenow=now;

% Store parameters
aaworker.id=aa_workerid;
aaworker.master.hostname=char(hostname);

% Get rid of the splash screens
aaworker_setspmvisibility

% Send PID and hostname to aaworker file
[s w]=unix('hostname');
clear workerdetails
aaworker.hostname=w;
aaworker.pid=aas_getworkerpid(aaworker.aap,w,aa_workerid);
aas_log([],0,strcat('PIDs are',sprintf('\t%d',aaworker.pid)));

fn=fullfile(aaworker.parmpath,'aaworker.xml');
xml_write(fn,aaworker);
aas_propagateto(aaworker.master.hostname,fn);

% Now register as bored
fn=fullfile(aaworker.parmpath,'iambored.mat');
aas_log([],0,'Ready to go');
save(fn,'timenow');
aas_propagateto(aaworker.master.hostname,fn);
aaworker.lastexcitement=clock;

aaworker_pollsetup

