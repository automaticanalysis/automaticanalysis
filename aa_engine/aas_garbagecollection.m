% function aap=aas_garbagecollection(aap, actuallydelete, modulestoscan )
%  Tidy up unnecessary files - when output is being split by module, delete
%  any inputs to a module that are not also specified as an output.
%  Rhodri Cusack  www.cusacklab.org  March 2012
%  Tibor Auer MRC CBU Cambridge 2012-2013

function aap=aas_garbagecollection(aap, actuallydelete, modulestoscan, permanencethreshold)

fprintf('Garbage collection started...\n');

if ~isstruct(aap) % first parameter is studypath (for "standanole" use)
	if ~isempty(aap), cd(aap); end
end

% First, load AAP structure and make sure that the one loaded is already
% processed
if ~exist('aap','var') || ~isstruct(aap)
    if ~exist('aap_parameters.mat','file')
        error('aap structure not found');
    else
        load('aap_parameters');
    end
end
studypath = fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid);
if ~exist(studypath,'dir')
    aas_log(aap,true,sprintf('Study %s not found',studypath));
end
cd(studypath);
if ~isfield(aap,'internal')
    if ~exist('aap_parameters.mat','file')
        error('aap structure not found');
    else
        load('aap_parameters');
    end
end

if ~strcmp(aap.directory_conventions.remotefilesystem,'none')
    aas_log(aap,true,'Remote file systems not currently supported by garbage collection');
end

if ~strcmp(aap.directory_conventions.outputformat,'splitbymodule')
    aas_log(aap,true,sprintf('No garbage collection as aap.directory_conventions.outputformat is %s',aap.directory_conventions.outputformat));
end

if ~exist('modulestoscan','var') || isempty(modulestoscan)
    modulestoscan=1:length(aap.tasklist.main.module);
    %     exclude modelling modules, if any
    ind = cell_index(strfind({aap.tasklist.main.module.name}, 'level'),'1');
    if ind, modulestoscan=1:ind-1; end
end

if (~exist('permanencethreshold','var')) || isempty(permanencethreshold)
    permanencethreshold = 0;
end

if ~exist('actuallydelete','var')
    actuallydelete=false;
end;

% What things fall under the image category?
imgCategory = {'.IMA' '.dcm' '.nii' '.img' '.hdr'};

numdel=0;

garbagelog=fullfile(studypath,sprintf('garbage_collection.txt'));
fprintf('Logfile: %s\n',garbagelog);
fid=fopen(garbagelog,'w');
fprintf(fid,'Garbage collected: %s\n',datestr(now,31));

for modind=modulestoscan
    inps={};
    outs={};
    % Get all the input and output stream files
    aap=aas_setcurrenttask(aap,modind);
    
    fprintf('\nWorking on task %s\n', aap.tasklist.currenttask.name)
    
    for streamname=aap.tasklist.currenttask.inputstreams.stream
        if ~isstruct(streamname{1})
            streamfn=sprintf('stream_%s_inputto_%s.txt',streamname{1},aap.tasklist.currenttask.name);
        else
            streamfn=sprintf('stream_%s_inputto_%s.txt',streamname{1}.CONTENT,aap.tasklist.currenttask.name);
        end
        inps=[inps findstreamfiles(aap,streamfn)];
    end
    for streamname=aap.tasklist.currenttask.outputstreams.stream
        if ~isstruct(streamname{1})
            streamfn=sprintf('stream_%s_outputfrom_%s.txt',streamname{1},aap.tasklist.currenttask.name);
        else
            streamfn=sprintf('stream_%s_outputfrom_%s.txt',streamname{1}.CONTENT,aap.tasklist.currenttask.name);
        end        
        outs=[outs findstreamfiles(aap,streamfn)];
    end
    
    % Load up these stream files
    inpfn={};
    outfn={};
    for inpind=1:length(inps)
        inpfn=[inpfn loadstreamfile(aap,inps{inpind})];
    end
    for outind=1:length(outs)
        outfn=[outfn loadstreamfile(aap,outs{outind})];
    end
    
    % Find which inputs aren't an output
    % N^2 string comparisons - might be costly [seems fine, though]
    inpnotout=false(size(inpfn));
    for inpind=1:length(inpfn)
        inpnotout(inpind)=~any(strcmp(inpfn{inpind},outfn));
    end
    
    % Garbage collect outputs
    if ~isfield(aap.tasklist.currenttask.settings, 'permanenceofoutput')
        % To be on the safe side keep stuff without a permanenceofoutput flag!
        aap.tasklist.currenttask.settings.permanenceofoutput = Inf;
    end
    if abs(aap.tasklist.currenttask.settings.permanenceofoutput) <= permanencethreshold
        delete_outputs = 1;
    else
        delete_outputs = 0;
    end
    if aap.tasklist.currenttask.settings.permanenceofoutput < 0
        % If permanenceofoutput has negative value, delete only images
        % .nii .img .hdr files
        delete_onlyImg = 1;
    else
        % Otherwise delete also non-image outputs...
        delete_onlyImg = 0;
    end
    
    % Garbage collect inputs
    if ~isempty(inpfn)
        fprintf(fid,'Garbage collection in module: %s\n',aap.tasklist.currenttask.name);
        if actuallydelete
            fprintf(fid,'---following files deleted---\n');
        end
        for fn=inpfn(inpnotout)
            if exist(fn{1},'file')
                fprintf(fid,'%s\n',fn{1});
                if actuallydelete, delete(fn{1}); end
                numdel=numdel+1;
            end
        end
        fprintf(fid,'---end---\n\n');
    end
    
    % Garbage collect outputs
    if (~isempty(outfn) && exist(studypath,'dir') && delete_outputs)
        if actuallydelete
            fprintf(fid,'---following files deleted---\n');
        end
        for fn=outfn
            if (exist(fn{1},'file'))
                if actuallydelete
                    if exist(fn{1},'file')
                        fprintf(fid,'%s\n',fn{1});
                        if actuallydelete, delete(fn{1}); end
                        numdel=numdel+1;
                    end
                end
                numdel=numdel+1;
            end
        end
        if (actuallydelete)
            fprintf(fid,'---end---');
        end
    end
    
    
end

fprintf(fid,'---end---\n');
fclose(fid);

if (~actuallydelete)
    aas_log(aap,false,sprintf('Garbage collection would delete %d files',numdel));
else
    aas_log(aap,false,sprintf('Garbage collection deleted %d files',numdel));
end
end


function [streampths]=findstreamfiles(aap,streamfn)
% Make a list of all the places the stream file might be
pthstocheck={aas_getstudypath(aap)};
for subjind=1:length(aap.acq_details.subjects)
    pthstocheck=[pthstocheck aas_getsubjpath(aap,subjind)];
    for sessind = 1:length(aap.acq_details.sessions)
        pthstocheck=[pthstocheck aas_getsesspath(aap,subjind,sessind)];
    end
end
pthstocheck = unique(pthstocheck); % remove duplicates (e.g. structural "session")
% Check to see which of these exist
streampths={};
for pthind=1:length(pthstocheck)
    if (exist(fullfile(pthstocheck{pthind},streamfn),'file'))
        streampths=[streampths fullfile(pthstocheck{pthind},streamfn)];
    end
end
end


function [imgfns]=loadstreamfile(aap,streampth)
fid=fopen(streampth,'r');
header=fgetl(fid);
% If we don't see MD5 at top, something odd is wrong, so throw
% error
if (~strcmp(header(1:3),'MD5'))
    aas_log(aap,true,sprintf('Garbage collection found stream file that does not begin with MD5, so quitting. File is %s',streampth));
end

% image files are relative to location of stream
[pth nme ext]=fileparts(streampth);

% now list all files
imgfns={};
while(~feof(fid))
    thisfn=fgetl(fid);
    if (~isempty(thisfn))
        imgfns=[imgfns fullfile(pth,thisfn)];
    end
end
fclose(fid);
end