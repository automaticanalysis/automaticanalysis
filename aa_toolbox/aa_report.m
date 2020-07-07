% Automatic analysis - produces HTML summary of analysis
%
% For reporting a remote pipeline, copy the aap_parameters.mat file to a
% local directory and call aa_report from that directory without any argument
%
% Rhodri Cusack MRC CBU Cambridge 2004
% Tibor Auer MRC CBU Cambridge 2012-2016

function [aap]=aa_report(studyroot,stages)

fprintf('Fetching report started...\n');

% First, load AAP structure
if exist('studyroot','var')
    if isstruct(studyroot)
        studyroot=fullfile(aas_getstudypath(studyroot),studyroot.directory_conventions.analysisid);
    end
    cd(studyroot);
else
    studyroot=pwd;
end
if ~exist('aap_parameters.mat','file')
    error('ERROR: aap structure not found');
else
    dat = load('aap_parameters');
    aap = dat.aap;
end

aa_init(aap);

if ~exist('stages','var')
    stages={aap.tasklist.main.module.name};
end

if isfield(aap,'report'), aap = rmfield(aap,'report'); end
aap.report.numdependencies=0;
aap.report.dependency={};

aap.internal.total=0;
aap.internal.stagesnotdone=0;

% Provenance
aap.prov = aa_provenance(aap);

% Flags for special reports
% - collect module executables
mfstages = stages;
for m = 1:numel(mfstages)
    if ~exist(mfstages{m},'file')
        if ~isempty(aap.tasklist.main.module(m).aliasfor)
            mfstages{m} = aap.tasklist.main.module(m).aliasfor;
        else
            mfstages{m} = aap.schema.tasksettings.(mfstages{m}).ATTRIBUTE.mfile_alias;
        end
    end
end
has_motioncorrection = ~isempty(intersect(mfstages,{'aamod_realign' 'aamod_realignunwarp'}));
has_registration = ~isempty(intersect(mfstages,{'aamod_norm_write' 'aamod_norm_write_dartel'}));
has_contrast = ~isempty(intersect(mfstages,{'aamod_firstlevel_threshold'}));
has_meegepochs = ~isempty(intersect(mfstages,{'aamod_meeg_epochs'}));

% Main HTLMs
if (isfield(aap.directory_conventions,'reportname'))
    aap.report.html_main.fname=fullfile(studyroot,aap.directory_conventions.reportname);
else
    aap.report.html_main.fname=fullfile(studyroot,'report.htm');
end
aap.report.html_S00.fname = strrep(aap.report.html_main.fname,'.htm','_subjects.htm');
aap.report.subdir = fullfile(fileparts(aap.report.html_main.fname),'report_subjects');
if ~exist(aap.report.subdir,'dir'), mkdir(aap.report.subdir); end
if has_motioncorrection, aap.report.html_moco.fname = strrep(aap.report.html_main.fname,'.htm','_moco.htm'); end
if has_registration, aap.report.html_reg.fname = strrep(aap.report.html_main.fname,'.htm','_reg.htm'); end
if has_contrast
    aap.report.html_C00.fname = strrep(aap.report.html_main.fname,'.htm','_scon.htm');
    aap.report.condir = fullfile(fileparts(aap.report.html_main.fname),'report_scon');
end
if has_meegepochs, aap.report.html_er.fname = strrep(aap.report.html_main.fname,'.htm','_er.htm'); end

% Initialize HTLMs
copyfile(fullfile(aap.internal.aappath,'aa_toolbox','aa_styles.css'),fullfile(studyroot,'aa_styles.css'));
aap = aas_report_add(aap,[],'HEAD=aa Report');
aap = aas_report_add(aap,0,'HEAD=Subjects');
if ~exist(aap.report.subdir,'dir'), mkdir(aap.report.subdir); end
if has_motioncorrection, aap = aas_report_add(aap,'moco','HEAD=Motion correction summary'); end
if has_registration, aap = aas_report_add(aap,'reg','HEAD=Registration summary'); end
if has_contrast
    aap = aas_report_add(aap,'C00','HEAD=First level results');
    if ~exist(aap.report.condir,'dir'), mkdir(aap.report.condir); end
end
if has_meegepochs, aap = aas_report_add(aap,'er','HEAD=MEEG epoch summary'); end

aap.report.fbase = basename(aap.report.html_main.fname);

for k=1:numel(stages)
    fprintf('Fetching report for %s...\n',stages{k});
    % domain
    curr_aap = aas_setcurrenttask(aap,k);
    domain = curr_aap.tasklist.currenttask.domain;
    
    % mfile_alias
    if isfield(aap.tasklist.main.module(k),'aliasfor') && ~isempty(aap.tasklist.main.module(k).aliasfor)
        xml = xml_read([aap.tasklist.main.module(k).aliasfor '.xml']);
    else
        xml = xml_read([stages{k} '.xml']);
    end
    mfile_alias = stages{k};
    if ~exist(mfile_alias,'file'), mfile_alias = aap.tasklist.main.module(k).aliasfor; end
    if ~exist(mfile_alias,'file'), mfile_alias = xml.tasklist.currenttask.ATTRIBUTE.mfile_alias; end
        
    % Switch for stage
    all_stage = cell_index(stages, stages{k});
    istage = cell_index({aap.tasklist.main.module.name}, stages{k});
    istage = istage(all_stage==k);
    
    % add provenance
    if aap.prov.isvalid, aap.prov.addModule(istage); end
    
    % Skip stages of unknown domain - most like aamod_importfilesasstream
    if ~strcmp(domain,'[unknown]')
        stagerepname = stages{k};
        if ~isempty(aap.tasklist.main.module(k).extraparameters)
            stagerepname = [stagerepname aap.tasklist.main.module(k).extraparameters.aap.directory_conventions.analysisid_suffix];
        end
        
        % build dependency
        [dep, domaintree] = aas_dependencytree_allfromtrunk(curr_aap,domain);
        
        % Set inSession flag
        inSession = (numel(domaintree) >= 2) & ~isempty(strfind(domain,'session'));

        % run through
        for d = 1:numel(dep)
            indices = dep{d}{2};
            isdone = exist(aas_doneflag_getpath_bydomain(aap,domain,indices,k),'file');
            if numel(dep{d}{2}), subj = dep{d}{2}(1); else, subj = []; end % Subjects No
            if numel(dep{d}{2}) >= 2, sess = dep{d}{2}(2); end % Session/Occurrance No
            
            if inSession
                sessdomain = domain;
                if numel(domaintree) > 2, sessdomain = domaintree{2}; end
                [junk, iSess] = aas_getN_bydomain(curr_aap,sessdomain,subj);
                firstSess = iSess(1);
                lastSess = iSess(end);
            end
            
            if ~inSession || ((sess == firstSess)  && (d == find(cellfun(@(x) x{2}(1) == subj, dep),1,'first'))), aap = aas_report_add(aap,subj,...
                    ['<h2>Stage: ' stagerepname '</h2>']); 
            end
            
            % evaluate with handling sessions
            if inSession
                if (sess == firstSess) && (d == find(cellfun(@(x) x{2}(1) == subj, dep),1,'first')), aap = aas_report_add(aap,subj,'<table><tr>'); end % Open session
                aap = aas_report_add(aap,subj,'<td valign="top">');
                aap = aas_report_add(aap,subj,['<h3>Session: ' aap.acq_details.([domaintree{2} 's'])(dep{d}{2}(2)).name '</h3>']);
                if numel(dep{d}{2}) >= 3
                    descSubSession = strrep(domaintree{3},'_',' '); descSubSession(1) = upper(descSubSession(1));
                    aap = aas_report_add(aap,subj,sprintf('<h4>%s: %d</h4>',descSubSession,dep{d}{2}(3))); 
                end
            end
            if ~isdone
                aap = aas_report_add(aap,subj,'<h3>Not finished yet!</h3>');
            else
                % set local aap
                if ~isempty(subj)
                    run_aap = aas_setcurrenttask(aap,k,'subject',subj);
                else
                    run_aap = curr_aap;
                end
                
                % set fields
                run_aap.report = aap.report;
                run_aap.prov = aap.prov;
         
                run_aap = aa_feval_withindices(mfile_alias,run_aap,'report',indices);
                
                % save report
                aap.report = run_aap.report;
                aap.prov = run_aap.prov;
            end
            if inSession
                aap = aas_report_add(aap,subj,'</td>');
                if (sess == lastSess) && (d == find(cellfun(@(x) x{2}(1) == subj, dep),1,'last')), aap = aas_report_add(aap,subj,'</tr></table>'); end % Close session
            end
        end
    end
end

% Close files
aap = aas_report_add(aap,[],'EOF');
aap = aas_report_add(aap,0,'EOF');
fclose all;

% Provenance
aap.prov.serialise(studyroot);

% Show report
if ~isdeployed, web(['file://' aap.report.html_main.fname]); end
% Last, save AAP structure
save(fullfile(studyroot,'aap_parameters_reported.mat'), 'aap');

aa_close(aap);

end