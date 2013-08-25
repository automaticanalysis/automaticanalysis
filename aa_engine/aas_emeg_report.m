function myurl=aas_emeg_report(varargin)
% Generates seperate html tables, with links to images, for:
% - sensor timecourses
% - source timecourses
% - field topographies (from magnetometers)
% - TF analyses
% - spm group source analyses
% - various diagnostics
%
% Input can be aap structure or path to aap_parameters file
%
% Danny Mitchell 13/03/08; Inspired by aa_report by Rhodri Cusack MRC CBU Cambridge 2004

if isempty(varargin)
    % ask for aap_parameters file
    addpath /imaging/local/spm/spm5
    addpath /imaging/local/spm/spm5/cbu_updates/
    [aap ok]=spm_select(1,'any','Please select aap_parameters file or script to generate aap structure','',pwd,'^aa.*\.m.*$');
    if ~ok; return; end
    [pth nam ext]=fileparts(aap);
    if strcmp(ext,'.m')
        addpath(pth); run(aap); rmpath(pth);
    end
else aap=varargin{1};
end

if ~isstruct(aap);
    try load(aap);
    catch error('\nFailed to load aap structure.\n')
    end
end

if ~strcmp(aap.acq_details.root(end),filesep)
    aap.acq_details.root=[aap.acq_details.root filesep];
end

outdir=fullfile(aap.acq_details.root,'htmlsummary');
if ~exist(outdir,'dir');mkdir(outdir);end

fprintf('\nCollecting images');
% get all image files for each subject
images={};ipaths={};
for s=1:length(aap.acq_details.subjects)
    subdir=fullfile(aap.acq_details.root,aap.acq_details.subjects(s).megname);
    ipaths{s}=spm_select('FPList',fullfile(subdir,'figures'),'^.*\.(jpe?g|png)$');
    for b=aap.acq_details.selected_sessions
        fprintf('.')
        blockdir=fullfile(subdir,aap.acq_details.sessions(b).name,'figures');
        ipaths{s}=char(ipaths{s},spm_select('FPList',blockdir,'^.*\.(jpe?g|png)$'));
    end
    % collect anonymised image names (without extension)
    images=[images; regexprep(cellstr(ipaths{s}),{sprintf('%s',aap.acq_details.subjects(s).megname),'\....$'},{'#',''})];
    % remove empty paths
    ipaths{s}(strmatch('/ ',ipaths{s}),:)=[];
end
% add any group timecourse analyses
if exist(fullfile(aap.acq_details.root,'GroupAnalysis_Sensors','figures'),'dir')
    filt='^g\d+.*(\.jpe?g|\.png|\.avi)$';
    %filt=sprintf('^g%g.*\\.(jpe?g|png)$',length(aap.acq_details.subjects)); % only get averages across all subjects
    ipaths{end+1}=spm_select('FPList',fullfile(aap.acq_details.root,'GroupAnalysis_Sensors','figures'),filt);
    % modify paths to look like single subject's files
    % (without ext because topo avi and png have same filename and don't
    % want duplicate row)
    images=[images; regexprep(cellstr(ipaths{end}),'GroupAnalysis_.*/g\d+g?(.*N?ST?t?\d?_)(\D*)(\d)(.*_)(.*)\..*','#/$2/figures/$1$2$3$4$5')];
    % find group sizes
    temp=ipaths{end}';
    gs=unique(regexp(temp(:)','/g\d+','match'));
else gs=[];
end

% get unique names
images=unique(images);
if length(images{1})<4; images=images(2:end); end
if isempty(images); return; end

% sort into types
ind=~cellfun('isempty',regexp(images,'/Ai[^/]')); tables(1).images=images(ind);images(ind)=[];tables(1).type='Sensor_ICA';
ind=~cellfun('isempty',regexp(images,'/t\d_[^/]')); tables(4).images=images(ind);images(ind)=[];tables(4).type='Sensor_TF';
ind=~cellfun('isempty',regexp(images,'/[^/].*sensors_\d+')); tables(2).images=images(ind);images(ind)=[];tables(2).type='Sensor_Timecourses';
ind=~cellfun('isempty',regexp(images,'/[^/].*_topo')); tables(3).images=images(ind);images(ind)=[];tables(3).type='Sensor_Fields';
ind=~cellfun('isempty',regexp(images,'/[^/].*sources_\d+_\d+')); tables(5).images=images(ind);images(ind)=[];tables(5).type='Source_Timecourses';
tables(end+1).images=images; tables(end).type='Diagnostics';

% Get spms as seperate variables. There are lots of these, but only at the
% group level, so makes sense to put them in seperate table
figdir=fullfile(aap.acq_details.root,'GroupAnalysis','figures');
if exist(figdir,'dir')
    filt=sprintf('^.*\\.(jpe?g|png)$');
    sourcespms=cellstr(spm_select('FPList',figdir,filt));
    sourcestems=unique(spm_str_manip(regexprep(sourcespms,'_-?\d+-\d+ms.*',''),'E'));    
    sourceeffects=unique(regexprep(sourcespms,{'.*Effect','_\d+\..*'},{'',''}));
    sourcenums=regexprep(sourcespms,'.*_(\d+)\..*','$1');
    sourcemaxnum=max(str2num(char(unique(sourcenums)))); % vector
    sourcetfwins={};
    for ss=1:length(sourcestems)
        temp=~cellfun('isempty',regexp(sourcespms,sourcestems{ss}));
        sourcetfwins{ss}=unique(regexprep(sourcespms(temp),'.*_(-?\d+-\d+ms.+Hz_\d?).*','$1'));
        % if no tf window, then svd option is irrelevant
        for tf=1:length(sourcetfwins{ss})
           if ~isempty(regexp(sourcetfwins{ss}{tf},'_Hz'))
               sourcetfwins{ss}{tf}=regexprep(sourcetfwins{ss}{tf},'Hz_\d?','Hz');
           end
        end
        sourcetfwins{ss}=unique(sourcetfwins{ss});
    end
% might need to also split induced from evoked here
end

% Get sensor spms for another table
figdir=fullfile(aap.acq_details.root,'GroupAnalysis_Sensors','figures','spms');
if exist(figdir,'dir')
    filt=sprintf('^.*\\.(jpe?g|png)$');
    sensorspms=spm_select('FPList',figdir,filt);
    sensorspmgroups=unique(regexprep(cellstr(sensorspms),{'.*Effect','_\d+\..*'},{'',''}));
    nums=regexprep(cellstr(sensorspms),'.*_(\d+)\..*','$1');
    sensorspmmaxnum=max(str2num(char(unique(nums)))); % vector
    twins=unique(regexprep(cellstr(sensorspms),'.*/([^_]*)_.*','$1'));
end

fprintf('\nCreating html tables:\n')
style='style="font-family:times;font-size:11"';
% create links to other tables and generic headings
top='<p>';
for tab=1:length(tables)
    top=sprintf('%s<a href="table_%s.html" target="tableframe" %s>%s</a>&nbsp;&nbsp;',top, tables(tab).type, style, tables(tab).type);
end
top=[top '<a href="sourcespmtable.html" target="tableframe" ' style '>Source_SPMs</a></p>'];
top=[top '<a href="sensorspmtable.html" target="tableframe" ' style '>Sensor_SPMs</a></p>'];
top=[top '<table width="100%" border="1" cellspacing="0"><tr>'];

% add other rows depending on fig type
for tab=1:length(tables)
    % create first row of column headings
    html=[top '<td width="1%" scope="col">&nbsp;</td>'];
    if tab==length(tables); html=[html sprintf('<td width=1%%" scope="col" %s>Block</td>',style)]; end
    for s=1:length(aap.acq_details.subjects)
        html=[html '<td width="1%" scope="col" ' style '>s' num2str(s) '</td>']; % each subject
    end;
    for g=1:length(gs)
        html=[html '<td width="1%" scope="col" ' style '>' gs{g}(2:end) '</td>']; % group analyses
    end
    html=[html '</tr>'];
    % add other rows
    for i=1:length(tables(tab).images)
        imname=regexprep(tables(tab).images{i},'.*(#|/)',''); % remove path
        imname=regexprep(imname,'\..*',''); % remove extension
        html=[html '<tr><td width="1%" scope="row" ' style '>' imname '</td>']; % image name in 1st column
        block='&nbsp;';
        for b=aap.acq_details.selected_sessions
            if length(findstr(tables(tab).images{i},aap.acq_details.sessions(b).name))>0; block=aap.acq_details.sessions(b).name; end
        end
        if tab==length(tables); html=[html '<td ' style '>' block '</td>']; end; % block name in 2nd column
        for s=1:length(aap.acq_details.subjects)+length(gs)
            % does this image exist for this subject or group analysis?
            txt=['<td ' style '>--</td>'];
            try
                if s<=length(aap.acq_details.subjects);
                    ffilt=[block '.*' imname '\.'];
                    q=regexp(cellstr(ipaths{s}),ffilt,'once');
                    if all(cellfun('isempty',q)); q=regexp(cellstr(ipaths{s}),['.*' imname '\.'],'once'); end
                else
                    ffilt=['.*' gs{s-length(aap.acq_details.subjects)} imname '\.'];
                    q=regexp(cellstr(ipaths{end}),ffilt);
                end

                p=find(~cellfun('isempty',q));
                if length(p)>1; fprintf('warning: >1 file matching %s\n',ffilt); end

                if ~isempty(p)
                    % check for correct group analysis and format path for windows
                    if s>length(aap.acq_details.subjects)
                        if isempty(regexp(ipaths{end}(p(1),:),gs{s-length(aap.acq_details.subjects)},'ONCE'));
                            gocatch;
                        end
                        pth=strrep(deblank(strrep(ipaths{end}(p(1),:),aap.acq_details.root,'../')),'/','\');
                    else
                        pth=strrep(deblank(strrep(ipaths{s}(p(1),:),aap.acq_details.root,'../')),'/','\');
                    end
                    txt=['<td ' style '><a href="' pth '" target="viewframe">go</a></td>'];
                end
            catch
                % possibly no group analysis
            end
            html=[html txt];
        end
        html=[html '</tr>'];
        fprintf('.')
    end
    % finish off
    fprintf('saving\n')
    html=[html '</table>'];
    writehtml(html,fullfile(outdir,sprintf('table_%s.html',tables(tab).type)))
end % next table for different fig type

% create table for spm (and ROI?) figures
if exist('sourcetfwins','var')
    % create first row of inversion headings
    % html=[top '<td width="1%" scope="col" ' style '>Number:</td>'];
    html=[top '<td width="1%" scope="col" ' style '>Inversion:</td>'];
    for ss=1:length(sourcestems)
        html=sprintf('%s<td width="1%%" colspan="%g" %s>%s</td>',html,sourcemaxnum*length(sourcetfwins{ss}),style,sourcestems{ss});
    end;
    html=[html '</tr>'];
    % create 2nd row for time windowa
    html=[html '<td width="1%" scope="col" ' style '>Time windows (ms):</td>'];
    for ss=1:length(sourcestems)
        for tw=1:length(sourcetfwins{ss})
            html=sprintf('%s<td width="1%%" colspan="%g" %s>%s</td>',html,sourcemaxnum,style,regexprep(sourcetfwins{ss}{tw},'ms.*',''));
        end
    end;
    html=[html '</tr>'];
    % create 3rd row for frequency bands
    html=[html '<td width="1%" scope="col" ' style '>Frequencies (Hz):</td>'];
    for ss=1:length(sourcestems)
        for tw=1:length(sourcetfwins{ss})
            html=sprintf('%s<td width="1%%" colspan="%g" %s>%s</td>',html,sourcemaxnum,style,regexprep(sourcetfwins{ss}{tw},{'.*ms_','Hz_'},{'',' v'}));
        end
    end;
    html=[html '</tr>'];
    % create 4th row for effect subnums
    html=[html '<td width="1%" scope="col" ' style '>Number:</td>'];
    for ss=1:length(sourcestems)
        for tw=1:length(sourcetfwins{ss})
        for i=1:sourcemaxnum
            html=sprintf('%s<td width="1%%" scope="col" %s>%g</td>',html,style,i);
        end
        end
    end;
    html=[html '</tr>'];
    % create other rows for effects
    for g=1:length(sourceeffects)
        fprintf('.')
        html=[html '<tr><td width="1%" scope="row" ' style '>' regexprep(sourceeffects{g},{'.*/','Effect_','.\*ms_'},{'','',''}) '</td>']; % economised image name in 1st column
     for ss=1:length(sourcestems)
       for tw=1:length(sourcetfwins{ss})
            for i=1:sourcemaxnum
                % does this image exist for this inversion/event?
                txt=['<td ' style '>--</td>'];
                try
                    filt=sprintf('%s_%g.*',['.*' sourcestems{ss} '.*' sourcetfwins{ss}{tw} '.*' sourceeffects{g}(2:end)],i);
                    q=regexp(sourcespms,filt,'tokens');%q=regexp(cellstr(sensorspms),filt,'tokens');
                    ind=find(~cellfun('isempty',q));
                    if ~isempty(ind)
                        [junk shortest]=min(cellfun(@length,sourcespms(ind)));
                        % if multiple matches, assume the shortest is correct
                            % format path for windows
                            pth=strrep(deblank(strrep(sourcespms{ind(shortest)},aap.acq_details.root,'../')),'/','\');
                            txt=['<td ' style '><a href="' pth '" target="viewframe" ' style '>go</a></td>'];
                            %break
                    end
                catch
                    warning('aaMEG:Caught','wassup?');
                    keyboard
                end
                html=[html txt];
            end
        end
     end
    end
% finish off
fprintf('saving\n')
html=[html '</table>'];
writehtml(html,fullfile(outdir,'sourcespmtable.html'));
end

if exist('twins','var')
    % create table for sensor spm figures
    % create first row of time window headings
    % html=[top '<td width="1%" scope="col" ' style '>Number:</td>'];
    html=[top '<td width="1%" scope="col" ' style '>Time window (ms):</td>'];
    for tw=1:length(twins)
        html=sprintf('%s<td width="1%%" colspan="%g" %s>%s</td>',html,sensorspmmaxnum,style,twins{tw});
    end;
    html=[html '</tr>'];
    % create 2nd row
    html=[html '<td width="1%" scope="col" ' style '>Number:</td>'];
    for tw=1:length(twins)
        for i=1:sensorspmmaxnum
            html=sprintf('%s<td width="1%%" scope="col" %s>%g</td>',html,style,i);
        end
    end;
    html=[html '</tr>'];
    % create other rows
    for g=1:length(sensorspmgroups)
        fprintf('.')
        html=[html '<tr><td width="1%" scope="row" ' style '>' regexprep(sensorspmgroups{g},{'.*/','Effect_','.\*ms_'},{'','',''}) '</td>']; % economised image name in 1st column
        for tw=1:length(twins)
            for i=1:sensorspmmaxnum
                % does this image exist for this inversion/event?
                txt=['<td ' style '>--</td>'];
                try
                    %filt=sprintf('%s_%g.*',sensorspmgroups{g},i);
                    %filt=sprintf('%s_%g.*',strrep(sensorspmgroups{g},'.*ms',['.*' twins{tw} 'ms']),i);
                    filt=sprintf('%s_%g.*',['.*' twins{tw} '.*' sensorspmgroups{g}(2:end)],i);
                    q=regexp(cellstr(sensorspms),filt,'tokens');%q=regexp(cellstr(sensorspms),filt,'tokens');
                    ind=find(~cellfun('isempty',q));
                    if ~isempty(ind)
                        for p=ind(1) %find(~cellfun('isempty',q))
                            % format path for windows
                            pth=strrep(deblank(strrep(sensorspms(p,:),aap.acq_details.root,'../')),'/','\');
                            txt=['<td ' style '><a href="' pth '" target="viewframe" ' style '>go</a></td>'];
                            %break
                        end
                    end
                catch
                    warning('aaMEG:Caught','wassup?');
                    keyboard
                end
                html=[html txt];
            end
        end
    end
    % finish off
    fprintf('\nSaving final files...')
    html=[html '</table>'];
    writehtml(html,fullfile(outdir,'sensorspmtable.html'))
end

if ~exist(fullfile(outdir,'viewframe-unix.html'),'file');
    viewframe='<html><head><title>aaMEG viewframe</title></head><body></body></html>';
    writehtml(viewframe,fullfile(outdir,'viewframe.html'))
end

reportfile=fullfile(outdir,'Report-unix.html');
if ~exist(reportfile,'file');
    frameset='<html><head><title>aaMEG Report</title></head>';
    frameset=[frameset '<frameset rows="*" cols="50%,50%" framespacing="0" frameborder="yes">'];
    frameset=[frameset '<frame src="table_Sensor_Timecourses.html" name="tableframe" />'];
    frameset=[frameset '<frame src="viewframe.html" name="viewframe" scrolling="no" />'];
    frameset=[frameset '</frameset></html>'];
    writehtml(frameset,fullfile(outdir,'Report.html'))
end

fprintf('Done\n')
myurl=['file://' reportfile];
web(myurl,'-new');

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writehtml(html,fname)
% save for windows
fid=fopen(fname,'w');
fprintf(fid,'%s',html);
fclose(fid);
% save for unix
unixhtml=regexprep(html,{'\','\.html'},{'/','-unix.html'});
fid=fopen(regexprep(fname,'\.html','-unix.html'),'w');
fprintf(fid,'%s',unixhtml);
fclose(fid);
return