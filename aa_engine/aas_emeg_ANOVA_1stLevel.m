function S=aas_emeg_ANOVA_1stLevel(aap,S,Overwrite)
% Work in progress. Should hopefully run 1st level contrasts of a
% 2- or 3-way repeated measures ANOVA with partitioned variance using a
% 2-level approach. See Rik's chapter.
%
% Currently somewhat specific to my experiment but potentially
% generalisable...
%
% S.voldir - directory (within each block) within which trialtype
%            directories are found
% S.fname - file name within trialtype directories
% S.diagnostic - flag for diagnostic plot
% S.factors - cell array of factor names (first factor being session type)
% S.levels  - vector of the number of S.levels per factor
% S.effects - cell array effects types
% e.g.   {...
%        'Com' 'Com'; ... % average of all conditions
%        'Dif' 'Com'; ... % main effect of factor 1
%        'Com' 'Dif'; ... % main effect of factor 2
%        'Dif' 'Dif'; ... % interaction of factors 1 and 2
%        'Sep' 'Dif'; ... % simple effect of factor 2 for each level of factor 1
%        };
%
% Replace a string with a vector for a specific contrast
% e.g. {'Sep' [-0.75 -0.25 0.25 0.75]} for a linear contrast across four
% levels of factor 2, at each level of factor 1.
%
% For mutliple regression, replace with a matrix:
% e.g. {'Sep' [-0.75 -0.25 0.25 0.75; ...
%              -1     0    0.5  0.5]} to get beta images for the multiple
%              regression with linear and plateauing predictors
%
% Following fields will be added to S, which will be saved as
% ANOVA_settings.mat in S.secondlevdir:
% .files - anonymised file names of component images
% .Knam - cell array of names of each 1st level effect
% .P - cell array of partitions of each Knam, for second level
% .Kfiles{x,e} - file names of first level contrasts e, for subjects x
% .secondlevdir - where this struct is saved and 2nd level analysis will go
% .datatype - 'erf', 'tf', or 'source'
% .nsub - number of subjects for whom 1st level contrasts were computed
%
% djm 04/08/08, 09/06/09

cd(aap.acq_details.root)
addpath /imaging/dm01/MoreTools

% load behavioural K estimates
warning off all
[Numeric,Txt]=xlsread('/imaging/dm01/MEG/17subs_0508/Kestimates.xls');
warning on all

% add marsbar paths
warning off all
root='/imaging/local/spm/marsbar/marsbar-0.40/';
[files dirs]=spm_select('List',root,'');
for d=1:size(dirs,1); addpath(strcat(root,dirs(d,:)),'-END'); end
warning on all

try diagnostic=S.diagnostic; catch diagnostic=true; end

fprintf('\nProcessing %s: %s: First level:',S.voldir,S.fname)

%% check data type (this might not be the best place or method!)
if ~isempty(regexp(S.voldir,'mt\d_','ONCE'))
    S.datatype='tf';
elseif ~isempty(regexp(S.fname,{'GS','IID','COH','LOR','ARD','MSP'},'ONCE'))
    S.datatype='source';
else S.datatype='erf';
end

%% get volume names for each cell of design
files=cell(prod(S.levels,2),1);
f=0;
for b=aap.acq_details.selected_sessions
    block=aap.acq_details.sessions(b).name;
    for dn=1:prod(S.levels(2:end))
        % assumes rotation of directories agrees with factor
        % specification! (Low to high frequency)
        f=f+1;
        files{f}=sprintf('%s/%s/trialtype%g/%s', ...
            block,regexprep(S.voldir,'BLOCK',block),dn,S.fname);
    end
end
clear b block f

%% set up 'initial component contrasts'
fprintf('\nDefining 1st level contrasts')
c=cell(length(S.factors),1);
d=c;
s=c;
for f=1:length(S.factors)
    c{f}=ones(S.levels(f),1); % common effects
    d{f}=-orth(diff(eye(S.levels(f)))'); % differential effects (orth is necessary to ensure the length of the contrast vectors are unity)
    s{f}=eye(S.levels(f))'; % seperate S.levels
end

%% set up contrasts for specified experimental effects
K=cell(size(S.effects,1),1); % contrasts for main effects and interactions etc.
Knam=K;
P=K;
t=cell(size(S.effects,1),length(S.factors));
p=t;
spm_figure('Getwin');
S.contrasts=S.effects; % user specified contrasts will be overwritten in latter
for e=1:size(S.effects,1)
    % terms to be kronned to get contrast matrix
    t(e,:)=c; % common effects by default
    usercon=find(cellfun(@isnumeric,S.effects)); % check for user-specified contrasts
    for uc=usercon % deal with any user-specified contrasts    
        t{e,uc}=orth(S.effects{e,uc}'); %(orth is necessary to ensure the length of the contrast vectors are unity)
        S.effects{e,uc}=sprintf('Con');
    end
    t(e,strmatch('Dif',S.effects(e,:)))=d(strmatch('Dif',S.effects(e,:)));
    t(e,strmatch('Sep',S.effects(e,:)))=s(strmatch('Sep',S.effects(e,:)));

    if length(S.factors)==3
        Knam{e}=sprintf('Effect_%s%s_%s%s_%s%s', ...
            S.effects{e,1},S.factors{1},S.effects{e,2},S.factors{2},S.effects{e,3},S.factors{3});
    elseif length(S.factors)==2
        Knam{e}=sprintf('Effect_%s%s_%s%s',S.effects{e,1},S.factors{1},S.effects{e,2},S.factors{2});
    else debugnow % currently assumes 2- or 3-way
    end
    % terms to be kronned to get partitions of 1st level design
    % matrix to be passed to 2nd level
    p(e,:)={1};
    for f=1:length(S.factors)
        if strcmp('Dif',S.effects(e,f)); p(e,f)={ones(1,size(d{f},2))}; end % F test to collapse over component contrasts
        if strcmp('Sep',S.effects(e,f)); p(e,f)=s(f); end
        if strncmp('Con',S.effects(e,f),3); p(e,f)={eye(size(t{e,f},2))}; end % Do seperate F/T tests for each predictor
    end
    K{e}=t{e,1};
    P{e}=p{e,1};
    for f=2:length(S.factors)
        K{e}=kron(K{e},t{e,f});
        P{e}=kron(P{e},p{e,f});
    end
    if diagnostic
        subplot(ceil(size(S.effects,1)),2,e*2-1)
        imagesc(K{e}',[-1 1])
        title(Knam{e})
        subplot(ceil(size(S.effects,1)),2,e*2)
        imagesc(P{e}',[-1 1])
        title('Partitions for 2nd level')
    end
end
drawnow
clear t e f p

%% load ROIs
ROIs=[];
if strcmpi(S.datatype,'source') && isfield(S,'ROIs') && iscell(S.ROIs) && ~isempty(S.ROIs);
    fprintf('\nLoading ROI definitions...')
    R=maroi(S.ROIs);
    dorois=true;
else dorois=false;
end

%% prepare directory for 2nd level contrasts
secondlevdir=fullfile(aap.acq_details.root,'GroupAnalysis',sprintf('FullFact_2ndLev_%s_%s',regexprep(S.voldir,{'-eeg','_BLOCK1'},{'',''}),S.fname));
try cd(secondlevdir);
catch mkdir(secondlevdir); cd(secondlevdir);
end

%% create 1st level contrasts for each subject
fprintf('\nCreating 1st level contrasts')
Kfiles={};

nsub=0;
for s=1:length(aap.acq_details.subjects)
    fprintf('\n%s: ',aap.acq_details.subjects(s).megname);
    subdir=fullfile(aap.acq_details.root,aap.acq_details.subjects(s).megname);
    fn=strcat(subdir,filesep,files);

    if ~all(cellfun(@exist,fn))
        fprintf(' - Files not found!'); continue
    end
    nsub=nsub+1;

    %% make (empty?) output directory
    firstlevdir=sprintf('%s/FullFact_1stLev_%s_%s/',subdir, ...
        regexprep(S.voldir,{'_BLOCK1'},{''}),S.fname);
    if ~exist(firstlevdir,'dir'); mkdir(firstlevdir);
        % else delete(fullfile(firstlevdir,'*'));
    end

    %% apply and save 1st level contrasts for main effects, contrasts and
    %% interactions of interest
    % Note that for orthgonalised predictors:
    % regress(Y,[ortho, ones(4,1)]) == [sum(Y.*ortho(:,1)); sum(Y.*ortho(:,2)); C]
    fprintf(' Applying contrasts')
    clear Y
    for e=1:length(K)
        for c=1:size(K{e},2)
            outfile=fullfile(firstlevdir,sprintf('%s_%g.img',Knam{e},c));
            if ~isempty(regexp(outfile,';-0')); debugnow; end %%%
            if c==1; Kfiles{nsub,e}=outfile;
            else Kfiles{nsub,e}=char(Kfiles{nsub,e},outfile);
            end
            if ~exist(outfile,'file') || Overwrite==1
                if ~exist('Y','var'); [Y Vout]=loadvolumes(fn); clear fn; end
                C=zeros(size(Y{1}));
                for y=1:length(Y)
                    C=C+Y{y}*K{e}(y,c);
                end
                % save
                Vout.descrip=regexprep(S.id,'.*_','');
                Vout.fname=outfile;
                try Vout=rmfield(Vout,{'private','pinfo'}); catch end
                spm_write_vol(Vout,C);
            end
        end
        fprintf('.')
    end

%     if dorois
%         % extract ROI data here for exporting to SPSS?
%         temp=fullfile(secondlevdir,aap.acq_details.subjects(s).megname);
%         for r=1:length(R)
%             marsy=get_marsy(R{r},char(fn),'mean');
%             roiData{r}=summary_data(marsy);
% 
%             roiTits{r}=textscan(sprintf('c%g,',1:length(roiData{r})),'%s','delimiter',',');
%             roiTits{r}=strcat(roiTits{r}{1},sprintf('r%g',r));
% 
%             % set 5th argument true to overwrite
%             aas_emeg_savestats(temp,roiTits{r},roiData{r},roiTits{r},true,sprintf('ROI_%g.mat',r));
%         end
%     end

    fprintf(' Done.');
end % next subject

S.Knam=Knam;
S.P=P;
S.Kfiles=Kfiles;
S.secondlevdir=secondlevdir;
S.nsub=nsub;

return

%%
function [Y Vout]=loadvolumes(fn)
% load volumes (only need to do this once)
fprintf('\n - Loading volumes')
V=spm_vol(fn);
Y=cell(length(V),1);
for v=1:length(V)
    Y{v}=single(spm_read_vols(V{v}));
    fprintf('.')
end
Vout=V{1};
fprintf('Done\n')
if ~iscell(Y); debugnow; end
return


