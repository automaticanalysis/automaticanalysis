function outputdirs=aas_emeg_ANOVA(aap,S,Overwrite)
% Work in progress. Should hopefully run 2- or 3-way repeated measures
% ANOVA with partitioned variance using a 2-level approach. Currently quite
% specific to my experiment but potentially generalisable...
% S.voldir - directory (within each block) within which trialtype
%            directories are found
% S.fname - file name with trialtype directories
% S.factors - cell array of factor names (first factor being session type)
% S.levels  - vector of the number of S.levels per factor
% S.effects - cell array effects types
% e.g.   {...
%        'Com' 'Com'; ... % average of all conditions
%        'Dif' 'Com'; ... % main effect of task
%        'Com' 'Dif'; ... % main effect of hand
%        'Dif' 'Dif'; ... % task x hand interaction
%        'Sep' 'Dif'; ... % simple effect of hand for each task
%        };
%
% djm 04/08/08
% See Rik's chapter

cd(aap.acq_details.root)
addpath /imaging/dm01/MoreTools
% load behavioural K estimates
warning off all
[Numeric,Txt]=xlsread('/imaging/dm01/MEG/17subs_0508/Kestimates.xls');
warning on all
% mean K-weighted contrast weights
meanX=[-0.476 -0.024 0.321 0.179];

% add marsbar paths
warning off all
    root='/imaging/local/spm/marsbar/marsbar-0.40/';
    [files dirs]=spm_select('List',root,'');
    for d=1:size(dirs,1); addpath(strcat(root,dirs(d,:)),'-END'); end
warning on all

fprintf('\nProcessing %s: %s: First level:',S.voldir,S.fname)

%% check data type
timefreq=false;
source=false;
erf=false;
if ~isempty(regexp(S.voldir,'mt\d_','ONCE'))
    timefreq=true;
elseif ~isempty(regexp(S.fname,{'GS','IID','COH','LOR','ARD','MSP'},'ONCE'))
    source=true;
else erf=true;
end

%% set up 'initial component contrasts'
fprintf('\nDefining 1st level contrasts')
c=cell(length(S.factors),1);
d=c;
s=c;
for f=1:length(S.factors)
    c{f}=ones(S.levels(f),1); % common effects
    d{f}=-orth(diff(eye(S.levels(f)))'); % differential effects
    s{f}=eye(S.levels(f))'; % seperate S.levels
end

%% set up contrasts for experimental effects
K=cell(size(S.effects,1),1); % contrasts for main effects and interactions etc.
Knam=K;
P=K;
t=cell(size(S.effects,1),length(S.factors));
p=t;
spm_figure('Getwin');
for e=1:size(S.effects,1)
    if length(S.factors)==3
        Knam{e}=sprintf('Effect_%s%s_%s%s_%s%s', ...
            S.effects{e,1},S.factors{1},S.effects{e,2},S.factors{2},S.effects{e,3},S.factors{3});
    elseif length(S.factors)==2
        Knam{e}=sprintf('Effect_%s%s_%s%s',S.effects{e,1},S.factors{1},S.effects{e,2},S.factors{2});
    else debugnow % currently assumes 2- or 3-way
    end
    % terms to be kronned to get design matrix
    t(e,:)=c;
    t(e,strmatch('Dif',S.effects(e,:)))=d(strmatch('Dif',S.effects(e,:)));
    t(e,strmatch('Sep',S.effects(e,:)))=s(strmatch('Sep',S.effects(e,:)));
    % terms to be kronned to get partitions of 1st level design
    % matrix to be passed to 2nd level
    p(e,:)={1};
    for f=1:length(S.factors)
        if strcmp('Dif',S.effects(e,f)); p(e,f)={ones(1,size(d{f},2))}; end
        if strcmp('Sep',S.effects(e,f)); p(e,f)=s(f); end
    end
    K{e}=t{e,1};
    P{e}=p{e,1};
    for f=2:length(S.factors)
        K{e}=kron(K{e},t{e,f});
        P{e}=kron(P{e},p{e,f});
    end
    %     %% diagnostics:
    %     subplot(ceil(size(S.effects,1)),2,e*2-1)
    %     imagesc(K{e}',[-1 1])
    %     title(Knam{e})
    %     subplot(ceil(size(S.effects,1)),2,e*2)
    %     imagesc(P{e}',[-1 1])
    %     %%%
    title('Partitions for 2nd level')
end
drawnow
clear t e f p

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

%% prepare directory for 2nd level contrasts
secondlevdir=fullfile(aap.acq_details.root,'GroupAnalysis',sprintf('FullFact_2ndLev_%s_%s',regexprep(S.voldir,{'-eeg','_BLOCK1'},{'',''}),S.fname));
try cd(secondlevdir);
catch mkdir(secondlevdir); cd(secondlevdir);
end

%% load binary FDR p<0.05 ROIs for SS effect for each task
% N='/imaging/dm01/MEG/17subs_0508/GroupAnalysis_Sensors/mace_pifST_vstm1_raw_SampleK_11/fdr.img';
ROIs=[];

%     if exist(secondlevdir,'dir')
%         % get ROIs from area of image
%         RN={};
%         RN{1}=fullfile(secondlevdir,'Effect_SepTask_ComCue_DifSetSize_1','0p05_FDR.img');
%         RN{2}=fullfile(secondlevdir,'Effect_SepTask_ComCue_DifSetSize_2','0p05_FDR.img');
%         ROIs=logical(spm_read_vols(spm_vol(char(RN))));
%         % or get ROIs as sensor voxels within this region
%         RN{1}=fullfile(secondlevdir,'Effect_SepTask_ComCue_DifSetSize_1','0p05_FDR_sensors.mat');
%         RN{2}=fullfile(secondlevdir,'Effect_SepTask_ComCue_DifSetSize_2','0p05_FDR_sensors.mat');
%         ROIs(:)=0;
%         for r=1:length(RN)
%             temp=zeros(size(ROIs(:,:,1,1)));
%             load(RN{r})
%             ind=sub2ind(size(temp),SensorVox(1,:),SensorVox(2,:));
%             temp(ind)=1;
%             ROIs(:,:,1,r)=logical(temp);
%         end
%         %% add flipped ROIs (to mirror lr although matrix has left at top)
%         ROIs=cat(4,ROIs,flipdim(ROIs,1));
%         RN=[RN, regexprep(RN,'\.','-flipped.')];
%         % imagesc(squeeze(sum(ROIs,4))); axis image
%         % ROIdata=zeros(length(aap.acq_details.subjects),8);
%         ROIdata=cell(length(aap.acq_details.subjects),prod(S.levels),length(RN));
%         % headings for columns suitable for spss
%         H={'Subject','vl1','vl2','vl4','vl6','vr1','vr2','vr4','vr6','el1','el2','el4','el6','er1','er2','er4','er6'};
%     end

% check if all been done before - may not be robust!
done=false;
[fils dirs]=spm_select('List',secondlevdir,'');
dirs=cellstr(dirs);
expected=sum(cell2mat(cellfun(@size,P,'UniformOutput', false)));
if length(dirs)==(2+expected(1)) && ~Overwrite
    for cdir=3:length(dirs)
        if ~exist(fullfile(secondlevdir,dirs{cdir},'0p05_none.img'),'file')
            break
        end
    end
    done=true;
end
%if done; continue; end % all 1st and 2nd level completed; haven't checked how robust this is!
clear done

if source && isfield(S,'ROIs') && iscell(S.ROIs) && ~isempty(S.ROIs);
    fprintf('\nLoading ROI definitions...')
    R=maroi(S.ROIs);
    dorois=true;
else dorois=false;
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

    if dorois
        % extract ROI data here for exporting to SPSS?
        %         marsy=get_marsy(R{:},char(fn),'mean');
        %         summary=summary_data(marsy);
        %         [a b]=ind2sub(size(summary),1:numel(summary));
        %         temp1=textscan(sprintf('c%g,',a),'%s','delimiter',',');
        %         temp2=textscan(sprintf('r%g,',b),'%s','delimiter',',');
        %         tits=strcat(temp1{1},temp2{1});
        %         vals=reshape(summary,1,[]);
        %         temp=regexprep(secondlevdir,'[^/]*N?ST?t\d_',''); % not ideal, but need shorter field names
        %         temp=fullfile(secondlevdir,aap.acq_details.subjects(s).megname);
        %         aas_emeg_savestats(temp,tits,vals,tits,false,'ROIdata.mat');
        temp=fullfile(secondlevdir,aap.acq_details.subjects(s).megname);
        for r=1:length(R)
            marsy=get_marsy(R{r},char(fn),'mean');
            summary=summary_data(marsy);
            %        [a b]=ind2sub(size(summary),1:numel(summary));

            summary(end+1)=mean([summary(1) summary(5)],2);
            summary(end+1)=mean([summary(2) summary(6)],2);
            summary(end+1)=mean([summary(3) summary(7)],2);
            summary(end+1)=mean([summary(4) summary(8)],2);
            summary(end+1)=mean([summary(1) summary(2) summary(3) summary(4)],2);
            summary(end+1)=mean([summary(5) summary(6) summary(7) summary(8)],2);

            summary(end+1)=mean([summary(9) summary(13)],2);
            summary(end+1)=mean([summary(10) summary(14)],2);
            summary(end+1)=mean([summary(11) summary(15)],2);
            summary(end+1)=mean([summary(12) summary(16)],2);
            summary(end+1)=mean([summary(9) summary(10) summary(11) summary(12)],2);
            summary(end+1)=mean([summary(13) summary(14) summary(15) summary(16)],2);

            tits=textscan(sprintf('c%g,',1:length(summary)),'%s','delimiter',',');
            tits=strcat(tits{1},sprintf('r%g',r));

            % true to overwrite
            aas_emeg_savestats(temp,tits,summary,tits,true,sprintf('ROI_%g.mat',r));
        end
    end

    %% make (empty?) output directory
    firstlevdir=sprintf('%s/FullFact_1stLev_%s_%s/',subdir, ...
        regexprep(S.voldir,{'_BLOCK1'},{''}),S.fname);
    if ~exist(firstlevdir,'dir'); mkdir(firstlevdir);
        % else delete(fullfile(firstlevdir,'*'));
    end

    %% apply and save 1st level contrasts for main effects and
    %% interactions of interest
    fprintf(' Applying contrasts')
    clear Y
    for e=1:length(K)
        for c=1:size(K{e},2)
            outfile=fullfile(firstlevdir,sprintf('%s_%g.img',Knam{e},c));
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
                Vout.fname=outfile;
                try Vout=rmfield(Vout,{'private','pinfo'}); catch end
                spm_write_vol(Vout,C);
            end
        end
        fprintf('.')
    end

    %% ROI summary from FDR-corrected SS effect, for each task
    if ~isempty(ROIs)
        for r=1:length(RN)
            ROIdata{s,1,r}=aap.acq_details.subjects(s).megname;
            for v=1:length(Y)
                ROIdata{s,v+1,r}=mean(Y{v}(ROIs(:,:,:,r)));
            end

            %% save ROI data to spreadsheet
            cell2csv(regexprep(RN{r},'/([^/]*)\....','_$1.csv'),[H;ROIdata(:,:,r)]);
        end
    end

    if ~isempty(findstr(S.voldir,'SampleK'))
        if ~isempty(regexp(firstlevdir,'eog|RMS','ONCE'))
            continue
        end
        %% individual K-weighted contrasts
        indiX=Numeric(strmatch(aap.acq_details.subjects(s).megname,Txt),:);
        indiX=indiX-mean(indiX);
        % meanX=[-0.476 -0.024 0.321 0.179]; moved outside loops

        conds={'VSTM_L', 'VSTM_R', 'ESTA_L', 'ESTA_R'};
        %named=false;
        for q=0:3 % each task and cue: VSTM_L VSTM_R ESTA_L ESTA_R
            %                 if ~named
            %                     Knam{end+1}=sprintf('K_%s',conds{q+1});
            %                     P{end+1}=1;
            %                     named=true;
            %                 end
            %                 Kfiles{s,e+q+1}=fullfile(firstlevdir,sprintf('%s.img',Knam{e+q+1}));

            outfile=fullfile(firstlevdir,sprintf('indiK_%s.img',conds{q+1}));
            if ~exist(outfile,'file') || Overwrite==1 ;
                if ~exist('Y','var'); [Y Vout]=loadvolumes(fn); clear fn; end
                indiKcon=zeros(size(Y{1}));
                for ss=1:4
                    indiKcon=indiKcon+indiX(ss)*Y{q*4+ss};
                end
                Vout.fname=outfile;
                try Vout=rmfield(Vout,{'private','pinfo'}); catch end
                spm_write_vol(Vout,indiKcon);
            end

            outfile=fullfile(firstlevdir,sprintf('meanK_%s.img',conds{q+1}));
            if ~exist(outfile,'file') || Overwrite==1 ;
                if ~exist('Y','var'); [Y Vout]=loadvolumes(fn); clear fn; end
                meanKcon=zeros(size(Y{1}));
                for ss=1:4
                    meanKcon=meanKcon+meanX(ss)*Y{q*4+ss};
                end
                Vout.fname=outfile;
                try Vout=rmfield(Vout,{'private','pinfo'}); catch end
                spm_write_vol(Vout,meanKcon);
            end
            fprintf('.')
            clear *Kcon
        end
        clear conds indiX
    end

    %         %% multiple regression of set-size data against individual K & SS at each voxel
    %         X=Numeric([1 strmatch(aap.acq_details.subjects(s).megname,Txt)],:);
    %         BetaK=zeros(size(Y{1}));
    %         BetaLinear=BetaK;
    %         conds={'VSTM_L', 'VSTM_R', 'ESTA_L', 'ESTA_R'};
    %         for q=0:3 % each task and cue: VSTM_L VSTM_R ESTA_L ESTA_R
    %             fprintf('.')
    %             if s==1
    %                 Knam{end+1}=sprintf('BetaK_%s',conds{q+1});
    %                 P{end+1}=1;
    %                 Knam{end+1}=sprintf('BetaLinear_%s',conds{q+1});
    %                 P{end+1}=1;
    %             end
    %             Kfiles{s,e+1+2*q}=fullfile(firstlevdir,sprintf('%s.img',Knam{e+1+2*q}));
    %             Kfiles{s,e+2+2*q}=fullfile(firstlevdir,sprintf('%s.img',Knam{e+2+2*q}));
    %             if ~exist(Kfiles{s,e+2+2*q},'file');
    %                 D=cat(3,Y{q*4+1},Y{q*4+2},Y{q*4+3},Y{q*4+4});
    %                 for x=1:size(D,1)
    %                     for y=1:size(D,2)
    %                         try
    %                             B=regress(squeeze(D(x,y,:)),[X' ones(4,1)]);
    %                             BetaLinear(x,y)=B(1);
    %                             BetaK(x,y)=B(2);
    %                         catch
    %                             BetaLinear(x,y)=NaN;
    %                             BetaK(x,y)=NaN;
    %                         end
    %                     end
    %                 end
    %                 Vout.fname=Kfiles{s,e+1+2*q};
    %                  try
    %                  Vout=rmfield(Vout,{'private','pinfo'});
    %                  catch end
    %                 spm_write_vol(Vout,BetaK);
    %                 Vout.fname=Kfiles{s,e+2+2*q};
    %                try Vout=rmfield(Vout,{'private','pinfo'}); catch end
    %                 spm_write_vol(Vout,BetaLinear);
    %             end
    %         end
    %     end

    %     %% save ROI data to spreadsheet
    %     for r=1:length(RN)
    %         cell2csv(regexprep(RN{r},'/([^/]*)\....','_$1.csv'),[H;ROIdata(:,:,r)]);
    %     end

    fprintf(' Done.');
end % next subject

if ~nsub; debugnow; end

%fprintf('Done\n(Not doing second level)\n'); return
fprintf('\nCompleted 1st Level. Starting 2nd level across %g subjects:',nsub)

%% for each event, get 1st level contrast files
outputdirs={};
for e=1:length(Knam) % effects of interest
    for kp=1:size(P{e},1); % each level of a 'separated' factor
        ind=repmat(P{e}(kp,:)',[size(Kfiles,1) 1]);
        files=char(Kfiles{:,e});
        files=files(ind==1,:);
        numcells=size(files,1)/size(Kfiles,1);
        condir=fullfile(secondlevdir,sprintf('%s_%g',Knam{e},kp));
        if ~exist(fullfile(condir,'ResMS.img'),'file') || Overwrite || exist(fullfile(condir,'ResI_0001.img'),'file')
            fprintf('\nCreating SPM for %s #%g',Knam{e},kp);
            clear factorial_design
            % create empty output directory & define 2nd level design
            factorial_design.dir{1}=condir;
            if ~exist(condir,'dir'); mkdir(condir);
            else delete(fullfile(condir,'*'));
            end
            factorial_design.dir{1}=condir;
            cd(factorial_design.dir{1})
            factorial_design.des.fd.fact.name='FirstLevelContrast';
            factorial_design.des.fd.fact.levels=numcells;
            factorial_design.des.fd.fact.dept=double(0);
            factorial_design.des.fd.fact.variance=double(1);
            factorial_design.des.fd.fact.gmsca=double(0);
            factorial_design.des.fd.fact.ancova=double(0);
            for c=1:numcells; % c=1:size(K{e},2)
                factorial_design.des.fd.icell(c).levels=c;
                factorial_design.des.fd.icell(c).scans=cellstr(files(mod(0:size(files,1)-1,numcells)+1==c,:));
            end
            factorial_design.masking.tm.tm_none=double([]);
            factorial_design.masking.im=double(1);
            %%% make explicit mask
            tv=spm_vol(factorial_design.des.fd.icell(c).scans{1});
            ty=spm_read_vols(tv);
            tv.fname=fullfile(condir,'explicit_mask.img');
            ts1=ty+flipdim(ty,1); % check for inverse symmetry in x
            ts2=ty-flipdim(ty,1); % check for symmetry in x
            if ~any(ts1(:)) || ~any(ts2(:))
                % data have symmetry so can mask out one side
                if timefreq; debugnow; end % shouldn't happen!
                ty(floor(end/2):end,:,:)=0;
            end
            ty=(~isnan(ty) & ty~=0); % mask any 0s or NaNs
            if source
                %By=spm_read_vols(spm_vol('/imaging/local/spm/spm5/apriori/brainmask.nii'));
                By=spm_read_vols(spm_vol('/imaging/local/spm/spm5/EEGtemplates/bsBinary_Template_Cortex_Mesh.img'));
                ty=ty.*By;
                brainvolume=true;
            else brainvolume=false;
            end
            if timefreq
                % ugly method!
                Ds=spm_select('FPList',fullfile(subdir,aap.acq_details.sessions(1).name), ...
                    regexprep(Kfiles{1,1},{'.*FullFact_1stLev_','_average-.*','(N?ST?t\d)'},{'.*','.*mat','$1.*'}));
                load(deblank(Ds(1,:)));
                ty=ty.*D.tf.mask';
                if size(ty,1)>300 %%% HACK FOR MY WORKING MEMORY DATA TO REDUCE SEARCH VOLUME
                    ty(1:89,:)=0; % mask out times prior to 200ms
                    ty(352:end,:)=0; % mask out times following 1250ms
                end
            end
            if size(ty,3)>300 %%% HACK FOR MY WORKING MEMORY DATA TO REDUCE SEARCH VOLUME
                ty(:,:,1:89)=0; % mask out times prior to 200ms
                ty(:,:,352:end)=0; % mask out times following 1250ms
            end
            try tv=rmfield(tv,{'private','pinfo'}); catch end
            spm_write_vol(tv,ty);
            %%%
            factorial_design.masking.em{1}=tv.fname;
            clear tv ty ts*
            factorial_design.globalc.g_omit=double([]);
            factorial_design.globalm.gmsca.gmsca_no=double([]);
            factorial_design.globalm.glonorm=double(1);
            factorial_design.dir{1}=condir;

            %% run job to build design
            clear jobs
            jobs{1}.stats{1}=struct('factorial_design',factorial_design);
            spm_jobman('run',jobs)
            clear factorial_design jobs files ind

            %% Estimate parameters
            spm_unlink(fullfile('.', 'mask.img')); % avoid overwrite dialog
            load SPM
            global defaults
            defaults.modality='EEG'; % otherwise non-sphericity correction might complain about no significant voxels?
            SPM = spm_spm(SPM);
            save('SPM.mat','SPM')
        else
            fprintf('\nLoading SPM for %s #%g',Knam{e},kp);
            load(fullfile(condir,'SPM.mat'));
        end

        %% initialise contrast fields
        indF=0; indPosT=0; indNegT=0;
        if ~isfield(SPM,'xCon') || Overwrite;
            SPM.xCon = [];
        else % check for exisitng contrasts
            for tc=1:length(SPM.xCon)
                if SPM.xCon(tc).STAT=='F' && strcmp(SPM.xCon(tc).name,Knam{e}); indF=tc; end
                if SPM.xCon(tc).STAT=='T' && strcmp(SPM.xCon(tc).name,Knam{e}); indPosT=tc; end
                if SPM.xCon(tc).STAT=='T' && strcmp(SPM.xCon(tc).name,['-' Knam{e}]); indNegT=tc; end
            end
        end
        
        % define contrasts
        cons2calc=[];
        if ~indF 
            if ~indF; indF=length(SPM.xCon)+1; end
            fprintf('\nDefining F contrast')
            F=eye(numcells);
            SPM.xCon(indF) = spm_FcUtil('Set',Knam{e},'F','c',F,SPM.xX.xKXs);
            cons2calc(end+1)=indF;
        else fprintf('\nFound F contrast')
        end
        if numcells==1
            if ~indPosT 
                if ~indPosT; indPosT=length(SPM.xCon)+1; end
                fprintf('\nDefining positive T contrast')
                SPM.xCon(indPosT) = spm_FcUtil('Set',Knam{e},'T','c',1,SPM.xX.xKXs);
                cons2calc(end+1)=indPosT;
            else fprintf('\nFound positive T contrast')
            end
            if ~indNegT
                if ~indNegT; indNegT=length(SPM.xCon)+1; end
                fprintf('\nDefining negative T contrast')
                SPM.xCon(indNegT) = spm_FcUtil('Set',['-' Knam{e}],'T','c',-1,SPM.xX.xKXs);
                cons2calc(end+1)=indNegT;
            else fprintf('\nFound negative T contrast')
            end
        end
        
        % calculate contrasts
        if ~isempty(cons2calc); 
            fprintf('\nCalculating contrasts:\n');
            SPM=spm_contrasts(SPM,cons2calc); 
        end
        
        outputdirs=[outputdirs; condir];
        
    end % next partition
end % next set of effects

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


