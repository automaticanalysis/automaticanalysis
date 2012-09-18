
%% CONNECT DATA PIPELINE BY IDENTIFYING SOURCE OF INPUT FROM OUTPUT STREAMS
%   Works back through dependencies to determine where inputs come from
function [aap]=aas_findinputstreamsources(aap)

% Make empty cell structures for input and output streams
aap.internal.inputstreamsources=cell(length(aap.tasklist.main.module),1);
aap.internal.outputstreamdestinations=cell(length(aap.tasklist.main.module),1);
for k1=1:length(aap.tasklist.main.module)
    aap.internal.inputstreamsources{k1}.stream=[];
    aap.internal.outputstreamdestinations{k1}.stream=[];
end;

% Now go through each module and find its input dependencies
%  then make bi-directional connections in inputstreamsources and
%  outputstreamdestinations
for k1=1:length(aap.tasklist.main.module)
    [stagepath stagename]=fileparts(aap.tasklist.main.module(k1).name);
    index=aap.tasklist.main.module(k1).index;
    
    % Find streams to be loaded remotely
    remotestream=aap.tasklist.main.module(k1).remotestream;
    
    if (isfield(aap.schema.tasksettings.(stagename),'inputstreams'))
        inputstreams=aap.schema.tasksettings.(stagename).inputstreams;
        streamlist=inputstreams.stream;
        if iscell(streamlist) && length(streamlist)>=1 && isstruct(streamlist{1}) && isfield(streamlist{1},'ATTRIBUTE') 
            streamlist=streamlist{1};
        end;
             
        for i=1:length(streamlist)
            if iscell(inputstreams.stream)  && ~isstruct(inputstreams.stream{1})
                inputstreamname=inputstreams.stream{i};
            else
                inputstreamname=inputstreams.stream{1}(i);
            end;
            ismodified=1;
            if isstruct(inputstreamname)
                if isfield(inputstreamname.ATTRIBUTE,'ismodified')
                    ismodified=inputstreamname.ATTRIBUTE.ismodified;
                end;
                inputstreamname=inputstreamname.CONTENT;
            end;
            findremote=[];
            if ~isempty(remotestream)
                findremote=find(strcmp(inputstreamname,{remotestream.stream}));
            end;
            if ~isempty(findremote)
                stream=[];
                stream.name=inputstreamname;
                stream.sourcenumber=-1;
                stream.sourcestagename=remotestream(findremote).stagetag;
                stream.sourcedomain=remotestream(findremote).sourcedomain;
                stream.depth=[];
                stream.host=remotestream(findremote).host;
                stream.aapfilename=remotestream(findremote).aapfilename;
                stream.ismodified=ismodified;
                if (isempty(aap.internal.inputstreamsources{k1}.stream))
                    aap.internal.inputstreamsources{k1}.stream=stream;
                else
                    aap.internal.inputstreamsources{k1}.stream(end+1)=stream;
                end;
                aas_log(aap,false,sprintf('Stage %s input %s comes from remote host %s stream %s',stagename,stream.name,stream.host,stream.sourcestagename));
            else
                
                [aap stagethatoutputs mindepth]=searchforoutput(aap,k1,inputstreamname,true,0,inf);
                if isempty(stagethatoutputs)
                    aas_log(aap,true,sprintf('Stage %s required input %s is not an output of any stage it is dependent on. You might need to add an aas_addinitialstream command or get the stream from a remote source.',stagename,inputstreamname));
                else
                    [sourcestagepath sourcestagename]=fileparts(aap.tasklist.main.module(stagethatoutputs).name);
                    sourceindex=aap.tasklist.main.module(stagethatoutputs).index;
                    aas_log(aap,false,sprintf('Stage %s input %s comes from %s which is %d dependencies prior',stagename,inputstreamname,sourcestagename,mindepth));
                    stream=[];
                    stream.name=inputstreamname;
                    stream.sourcenumber=stagethatoutputs;
                    stream.sourcestagename=sourcestagename;
                    stream.sourcedomain=[];
                    stream.depth=mindepth;
                    stream.host='';
                    stream.aapfilename='';
                    stream.ismodified=ismodified;
                    stream.sourcedomain=aap.schema.tasksettings.(sourcestagename)(sourceindex).ATTRIBUTE.domain;
                    if (isempty(aap.internal.inputstreamsources{k1}.stream))
                        aap.internal.inputstreamsources{k1}.stream=stream;
                    else
                        aap.internal.inputstreamsources{k1}.stream(end+1)=stream;
                    end;
                    
                    stream=[];
                    stream.name=inputstreamname;
                    stream.destnumber=k1;
                    stream.deststagename=stagename;
                    stream.depth=mindepth;
                    stream.destdomain=aap.schema.tasksettings.(stagename)(index).ATTRIBUTE.domain;
                    if (isempty(aap.internal.outputstreamdestinations{stagethatoutputs}.stream))
                        aap.internal.outputstreamdestinations{stagethatoutputs}.stream=stream;
                    else
                        aap.internal.outputstreamdestinations{stagethatoutputs}.stream(end+1)=stream;
                    end;
                end;
            end;
        end;
    end;
end;


% RECURSIVELY SEARCH DEPENDENCIES
%  to see which will have outputted each
%  input required for this stage
%  note that inputs can be affected by dependency map
function [aap,stagethatoutputs,mindepth]=searchforoutput(aap,currentstage,outputtype,notthislevelplease,depth,mindepth)

% is this branch ever going to do better than we already have?
if (depth>=mindepth)
    return;
end;

stagethatoutputs=[];

% Search the current level, see if it provides the required output
if (~notthislevelplease)
    depth=depth+1;
    [stagepath stagename]=fileparts(aap.tasklist.main.module(currentstage).name);
    stagetag=aas_getstagetag(aap,currentstage);
    index=aap.tasklist.main.module(currentstage).index;
    
    if (isfield(aap.schema.tasksettings.(stagename)(index),'outputstreams'))
        outputstreams=aap.schema.tasksettings.(stagename)(index).outputstreams;
        for i=1:length(outputstreams.stream)
            if (strcmp(outputtype,outputstreams.stream{i}) || strcmp(outputtype,[stagetag '.' outputstreams.stream{i}]))
                stagethatoutputs=currentstage;
                mindepth=depth;
            end;
        end;
    end;
end;

% If not found, search backwards further
if (isempty(stagethatoutputs))
    dependenton=aap.internal.dependenton{currentstage};
    for i=1:length(dependenton)
        [aap stagethatoutputs mindepth]=searchforoutput(aap,dependenton(i).stage,outputtype,false,depth,mindepth);
    end;
end;