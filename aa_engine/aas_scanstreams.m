function aap=aas_garbagecollection(aap, modulestoscan)

if (~strcmp(aap.directory_conventions.remotefilesystem,'none'))
    aas_log(aap,true,'Remote file systems not currently supported by garbage collection');
end;

if (~exist('modulestoscan','var'))
    modulestoscan=1:length(aap.tasklist.main.module);
end;

for modind=modulestoscan
    inps={};
    outs={};
    % Get all the input and output stream files
    aap=aas_setcurrenttask(aap,modind);
    for streamname=aap.tasklist.currenttask.inputstreams.stream
        streamfn=sprintf('stream_%s_inputto_%s.txt',streamname{1},aap.tasklist.currenttask.name);
        inps=[inps findstreamfiles(aap,streamfn)];
    end;
    for streamname=aap.tasklist.currenttask.outputstreams.stream
        streamfn=sprintf('stream_%s_outputfrom_%s.txt',streamname{1},aap.tasklist.currenttask.name);
        outs=[outs findstreamfiles(aap,streamfn)];
    end;
    
    % Load up these stream files
    inpfn={};
    outfn={};
    for inpind=1:length(inps)
        inpfn=[inpfn loadstreamfile(aap,inps{inpind})];
    end;
    for outind=1:length(outs)
        outfn=[outfn loadstreamfile(aap,outs{outind})];
    end;
    
    % Find which inputs aren't an output
    % N^2 string comparisons - might be costly [seems fine, though]
    inpnotout=false(size(inpfn));
    for inpind=1:length(inpfn)
        inpnotout(inpind)=~any(strcmp(inpfn{inpind},outfn));         
    end;
    
    % Garbage collect
    if (~isempty(inpfn) && exist(aas_getstudypath(aap),'file'))
        garbagelog=fullfile(aas_getstudypath(aap),sprintf('garbage_collection.txt'));
        fid=fopen(garbagelog,'a');
        fprintf(fid,'Garbage collected: %s\n',datestr(now,31));
        fprintf(fid,'---following files deleted---\n');
        for fn=inpfn(inpnotout)
            fprintf(fid,'%s\n',fn{1});
        end;
        fprintf(fid,'---end---');
        fclose(fid);
    end;
end;


end


function [streampths]=findstreamfiles(aap,streamfn)
% Make a list of all the places the stream file might be
pthstocheck={aas_getstudypath(aap)};
for subjind=1:length(aap.acq_details.subjects)
    pthstocheck=[pthstocheck aas_getsubjpath(aap,subjind)];
    for sessind = 1:length(aap.acq_details.sessions)
        pthstocheck=[pthstocheck aas_getsesspath(aap,subjind,sessind)];        
    end;
end;
% Check to see which of these exist
streampths={};
for pthind=1:length(pthstocheck)
    if (exist(fullfile(pthstocheck{pthind},streamfn),'file'))
        streampths=[streampths fullfile(pthstocheck{pthind},streamfn)];
    end;
end;
end



function [imgfns]=loadstreamfile(aap,streampth)
fid=fopen(streampth,'r');
header=fgetl(fid);  
% If we don't see MD5 at top, something odd is wrong, so throw
% error
if (~strcmp(header(1:3),'MD5'))
    aas_log(aap,true,sprintf('Garbage collection found stream file that does not begin with MD5, so quitting. File is %s',streampth));
end;

% image files are relative to location of stream
[pth nme ext]=fileparts(streampth);

% now list all files
imgfns={};
while(~feof(fid))
    thisfn=fgetl(fid);
    if (~isempty(thisfn))
        imgfns=[imgfns fullfile(pth,thisfn)];
    end;
end;
end
