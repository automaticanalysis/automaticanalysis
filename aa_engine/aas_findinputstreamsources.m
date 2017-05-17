
%% CONNECT DATA PIPELINE BY IDENTIFYING SOURCE OF INPUT FROM OUTPUT STREAMS
%   Works back through dependencies to determine where inputs come from
function [aap]=aas_findinputstreamsources(aap)

% save provenance
provfn = fullfile(aap.acq_details.root,aap.directory_conventions.analysisid,'aap_prov.dot');
pfid = fopen(provfn,'w');
fprintf(pfid,'digraph {\n\trankdir = LR;\n\tcharset="utf-8";\n\n');
dotR = {};

cmapfn = fullfile(aap.acq_details.root,aap.directory_conventions.analysisid,'aap_cmap.txt');
cfid = fopen(cmapfn,'w');

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
    
    if (isfield(aap.schema.tasksettings.(stagename)(index),'inputstreams'))
        inputstreams=aap.schema.tasksettings.(stagename)(index).inputstreams;
        
        
        % Random evil check?  Why is this here??
        %         if iscell(streamlist) && length(streamlist)>=1 && isstruct(streamlist{1}) && isfield(streamlist{1},'ATTRIBUTE')
        %             streamlist=streamlist{1};
        %         end;
        
        % RC 2014-04-11
        % Commenting out of the code above fixed one case but broke another
        % It is a nuisance, but the XML parser returns these stream elements of the
        % XML in different ways, depending on whether they have attributes or not.
        % [1] If they have no attributes, it returns a cell array of stream names (char)
        % [2] If the first has attributes, but the next doesn't, you get a cell array
        % containing a struct followed by a char
        % [3] If the first two have attributes, you get an cell array with one element,
        % which is a struct array.
        % The "evil" lines above fixed [3] but created a problem with [2]
        
        % Here's a refactored version of this code, which first unpacks
        % those structs
        streamlist={};
        for ind=1:length(inputstreams.stream)
            if isstruct(inputstreams.stream{ind}) && length(inputstreams.stream{ind})>1
                for structind=1:length(inputstreams.stream{ind})
                    streamlist{end+1}=inputstreams.stream{ind}(structind);
                end;
            else
                streamlist{end+1}=inputstreams.stream{ind};
            end;
        end;
        
        % Now we have one stream per cell
        for i=1:length(streamlist)
            inputstreamname=inputstreams.stream{i};
            ismodified=1; isessential=1;
            if isstruct(inputstreamname)
                if isfield(inputstreamname.ATTRIBUTE,'ismodified')
                    ismodified=inputstreamname.ATTRIBUTE.ismodified;
                end;
                if isfield(inputstreamname.ATTRIBUTE,'isessential')
                    isessential=inputstreamname.ATTRIBUTE.isessential;
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
                stream.sourcemodality=remotestream(findremote).sourcemodality;
                stream.depth=[];
                stream.host=remotestream(findremote).host;
                stream.aapfilename=remotestream(findremote).aapfilename;
                stream.ismodified=ismodified;
                stream.allowcache=remotestream(findremote).allowcache;
                stream.isessential=isessential;
                if (isempty(aap.internal.inputstreamsources{k1}.stream))
                    aap.internal.inputstreamsources{k1}.stream=stream;
                else
                    aap.internal.inputstreamsources{k1}.stream(end+1)=stream;
                end;
                aas_log(aap,false,sprintf('Stage %s input %s comes from remote host %s stream %s',stagename,stream.name,stream.host,stream.sourcestagename));
                % write provenance
                fprintf(pfid,'"R%s" -> "R%s_%05d" [ label="%s" ];',['Remote ' stream.host ': ' stream.sourcestagename],stagename,index,stream.name);
                dotR(end+1) = cellstr(['Remote ' stream.host ': ' stream.sourcestagename]);
                dotR(end+1) = cellstr(sprintf('%s_%05d',stagename,index));
				fprintf(cfid,'%s\t%s\t%s_%05d\n',['Remote ' stream.host ': ' stream.sourcestagename],stream.name,stagename,index);
            else
                
                [aap stagethatoutputs mindepth]=searchforoutput(aap,k1,inputstreamname,true,0,inf);
                if isempty(stagethatoutputs)
                    if isessential
                        aas_log(aap,true,sprintf('Stage %s required input %s is not an output of any stage it is dependent on. You might need to add an aas_addinitialstream command or get the stream from a remote source.',stagename,inputstreamname));
                    end
                else
                    [sourcestagepath sourcestagename]=fileparts(aap.tasklist.main.module(stagethatoutputs).name);
                    sourceindex=aap.tasklist.main.module(stagethatoutputs).index;
                    aas_log(aap,false,sprintf('Stage %s input %s comes from %s which is %d dependencies prior',stagename,inputstreamname,sourcestagename,mindepth));
                    % write provenance
                    fprintf(pfid,'"R%s_%05d" -> "R%s_%05d" [ label="%s" ];',sourcestagename,sourceindex,stagename,index,inputstreamname);
                    dotR(end+1) = cellstr(sprintf('%s_%05d',sourcestagename,sourceindex));
                    dotR(end+1) = cellstr(sprintf('%s_%05d',stagename,index));
					fprintf(cfid,'%s_%05d\t%s\t%s_%05d\n',sourcestagename,sourceindex,inputstreamname,stagename,index);
                    stream=[];
                    stream.name=inputstreamname;
                    stream.sourcenumber=stagethatoutputs;
                    stream.sourcestagename=sourcestagename;
                    stream.sourcedomain='';
                    stream.sourcemodality='';
                    stream.depth=mindepth;
                    stream.host='';
                    stream.aapfilename='';
                    stream.ismodified=ismodified;
                    stream.allowcache=false;
                    stream.isessential=isessential;
                    stream.sourcedomain=aap.schema.tasksettings.(sourcestagename)(sourceindex).ATTRIBUTE.domain;
                    stream.sourcemodality=aap.schema.tasksettings.(sourcestagename)(sourceindex).ATTRIBUTE.modality;
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
            if ~isempty(findremote) || ~isempty(stagethatoutputs)
                % change domain and modality if needed (due to input)
                currstage = aap.schema.tasksettings.(stagename)(index).ATTRIBUTE;
                currstream = aap.internal.inputstreamsources{k1}.stream(end);
                if (strcmp(currstage.domain,'*') && strcmp(currstage.modality,'MRI')) && ... % general-purpose module
                        (~isempty(strfind(currstream.sourcedomain,'session')) || ~strcmp(currstream.sourcemodality,'MRI')) % special
                    aap.schema.tasksettings.(stagename)(index).ATTRIBUTE.modality = currstream.sourcemodality;
                    aap.internal.aap_initial.schema.tasksettings.(stagename)(index).ATTRIBUTE.modality = currstream.sourcemodality;
                    aap.schema.tasksettings.(stagename)(index).ATTRIBUTE.domain = currstream.sourcedomain;
                    aap.internal.aap_initial.schema.tasksettings.(stagename)(index).ATTRIBUTE.domain = currstream.sourcedomain;
                end
            end
        end;
    end;
end;
fprintf(pfid,'\n\t// Resources\n');
dotR = unique(dotR,'stable');
for r = 1:numel(dotR)
    fprintf(pfid,'\t"R%s" [ label="%s", shape = ellipse, color = blue ]; \n',dotR{r},dotR{r});
end
fprintf(pfid,'\n\t// Anonymous nodes\n\n\t// Literals\n\n\tlabel="\\n\\nModel:\\n(Unknown)";\n}\n');
fclose(pfid);
% create provenance map
if ~unix('which dot')
    unix(sprintf('dot %s -Grankdir=TB -Tpng -o %s',provfn,strrep(provfn,'dot','png')));
end
fclose(cfid);

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
        outputstreams=aas_getstreams(aap,stagename,index,'output');
        
        for i=1:length(outputstreams)
            % get simple streamname
            streamsimple = '';
            dot_ind = find(outputstreams{i}=='.', 1, 'last');
            if ~isempty(dot_ind), streamsimple = outputstreams{i}(dot_ind+1:end); end
            
            if (strcmp(outputtype,outputstreams{i}) || ...
                    strcmp(outputtype,[stagetag '.' outputstreams{i}]) || ...
                    strcmp(outputtype,streamsimple))
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