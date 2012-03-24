function aas_emeg_rmANOVA_partitioned(aap,sensorvoldirs,suffixes,Overwrite)
% Work in progress. Should hopefully run 2- or 3-way repeated measures 
% ANOVA with partitioned variance using a 2-level approach. Currently very 
% specific to my experiment but potentially generalisable...
% factors - cell array of factor names (first factor being session type)
% levels  - vector of the number of levels per factor
%
% See Rik's chapter
%
% djm 04/08/08

cd(aap.acq_details.root)
addpath /imaging/dm01/MoreTools

% load behavioural K estimates
warning off all
[Numeric,Txt]=xlsread('/imaging/dm01/MEG/17subs_0508/Kestimates.xls');
warning on all
% mean K-weighted contrast weights
meanX=[-0.476 -0.024 0.321 0.179];

%% get directories with sensor volumes -% now being passed in
%d=dir(fullfile(aap.acq_details.root cellstr(char(q{1,[q{4,:}]}))

%% sensor type suffixes also now being passed in

%% do it
for vd=1:length(sensorvoldirs)
    for suf=1:length(suffixes)
        suffix=suffixes{suf};
        if ~isempty(findstr(sensorvoldirs{vd},'Hand'))
            if ~isempty(regexp(sensorvoldirs{vd},'mt\d_','ONCE'))
                %fnames={['average' regexprep(suffix,'eeg','')]};
                fnames={['average' suffix]};
                timefreq=true; 
            else
                fnames={['avewin_-10-20ms' suffix],...
                ['lsaverage' suffix],...
                ['bsaverage' suffix]}; % sensor time volume
                timefreq=false;
            end
            factors={'Task','Hand'};
            levels=[2 2];
            E= {'Com' 'Com'; ... % average of all conditions
                'Dif' 'Com'; ... % main effect of task
                'Com' 'Dif'; ... % main effect of hand
                'Dif' 'Dif'; ... % task x hand interaction
                };
        else
            if ~isempty(regexp(sensorvoldirs{vd},'mt\d_','ONCE'))
                %fnames={['average' regexprep(suffix,'eeg','')]};
                fnames={['average' suffix]};
                timefreq=true; 
            else
                fnames={['avewin_150-190ms' suffix],...
                    ['avewin_550-1250ms' suffix],...
                    ['avewin_550-900ms' suffix],...
                    ['avewin_900-1250ms' suffix],...
                    ['lsaverage' suffix],...
                    ['bsaverage' suffix]}; % sensor time volume
                timefreq=false;
            end
            factors={'Task','Cue','SetSize'};
            levels=[2 2 4];
            E= {'Dif' 'Com' 'Com'; ... % main effect of task
                'Com' 'Dif' 'Com'; ... % main effect of cue
                'Com' 'Com' 'Dif'; ... % main effect of SS
                'Dif' 'Dif' 'Com'; ... % task x cue interaction
                'Dif' 'Com' 'Dif'; ... % task x SS interaction
                'Com' 'Dif' 'Dif'; ... % cue x SS interaction
                'Dif' 'Dif' 'Dif'; ... % task x cue x SS interaction
                'Com' 'Sep' 'Dif'; ... % main effect of SS for each cue
                'Sep' 'Com' 'Dif'; ... % main effect of SS for each task
                'Sep' 'Sep' 'Dif'; ... % main effect of SS for each task and cue
                'Sep' 'Dif' 'Com'; ... % main effect of cue for each task
                'Sep' 'Dif' 'Dif'; ... % cue x SS interaction for each task
                'Sep' 'Com' 'Com'; ... % average of all conditions for each task
                };
            %'Com' 'Com' 'Com'; ... % average of all conditions
            %'Com' 'Com' 'Sep'; ... % average effect at each SS
        end

        for tfn=1:length(fnames)
            fname=fnames{tfn};
            %% 1st level contrasts
            fprintf('\nProcessing %s: %s: First level:',sensorvoldirs{vd},fname)

            %% set up 'initial component contrasts'
            fprintf('\nDefining 1st level contrasts')
            c=cell(length(factors),1);
            d=c;
            s=c;
            for f=1:length(factors)
                c{f}=ones(levels(f),1); % common effects
                d{f}=-orth(diff(eye(levels(f)))'); % differential effects (A-B)
                s{f}=eye(levels(f))'; % seperate levels
            end

            %% set up contrasts for experimental effects
            K=cell(size(E,1),1); % contrasts for main effects and interactions etc.
            Knam=K;
            P=K;
            t=cell(size(E,1),length(factors));
            p=t;
            spm_figure('Getwin');
            for e=1:size(E,1)
                if length(factors)==3
                    Knam{e}=sprintf('Effect_%s%s_%s%s_%s%s',E{e,1},factors{1},E{e,2},factors{2},E{e,3},factors{3});
                elseif length(factors)==2
                    Knam{e}=sprintf('Effect_%s%s_%s%s',E{e,1},factors{1},E{e,2},factors{2});
                else debugnow % currently assumes 2- or 3-way
                end
                % terms to be kronned to get design matrix
                t(e,:)=c;
                t(e,strmatch('Dif',E(e,:)))=d(strmatch('Dif',E(e,:)));
                t(e,strmatch('Sep',E(e,:)))=s(strmatch('Sep',E(e,:)));
                % terms to be kronned to get partitions of 1st level design matrix to be passed to 2nd level
                %p(e,:)=cellfun(@transpose,c,'UniformOutput',false);
                p(e,:)={1};
                for f=1:length(factors)
                    if strcmp('Dif',E(e,f)); p(e,f)={ones(1,size(d{f},2))}; end
                    if strcmp('Sep',E(e,f)); p(e,f)=s(f); end
                end
                K{e}=t{e,1};
                P{e}=p{e,1};
                for f=2:length(factors)
                    K{e}=kron(K{e},t{e,f});
                    P{e}=kron(P{e},p{e,f});
                end
                %%% diagnostics:
%                 subplot(ceil(size(E,1)),2,e*2-1)
%                 imagesc(K{e}',[-1 1])
%                 title(Knam{e})
%                 subplot(ceil(size(E,1)),2,e*2)
%                 imagesc(P{e}',[-1 1])
                title('Partitions for 2nd level')
            end
            drawnow
            clear t e f p
            
            %% get volume names for each cell of design
            files=cell(prod(levels,2),1);
            f=0;
            for b=aap.acq_details.selected_sessions
                block=aap.acq_details.sessions(b).name;
                for dn=1:prod(levels(2:end))
                    % assumes rotation of directories agrees with factor
                    % specification!
                    f=f+1;
                    files{f}=sprintf('%s/%s/trialtype%g/%s.img', ...
                        block,regexprep(sensorvoldirs{vd},'BLOCK',block),dn, ...
                        regexprep(fname,'([^-])eeg$','$1'));
                end
            end
            clear b block f
            %% load binary FDR p<0.05 ROIs for SS effect for each task
            % N='/imaging/dm01/MEG/17subs_0508/GroupAnalysis_Sensors/mace_pifST_vstm1_raw_SampleK_11/fdr.img';
            ROIs=[];
            secondlevdir=fullfile(aap.acq_details.root,'GroupAnalysis_Sensors',sprintf('FullFact_2ndLev_%s_%s',regexprep(sensorvoldirs{vd},{'-eeg','_BLOCK1'},{'',''}),fname));

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
            %         ROIdata=cell(length(aap.acq_details.subjects),prod(levels),length(RN));
            %         % headings for columns suitable for spss
            %         H={'Subject','vl1','vl2','vl4','vl6','vr1','vr2','vr4','vr6','el1','el2','el4','el6','er1','er2','er4','er6'};
            %     end

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
            
            %% create 1st level contrasts for each subject
            fprintf('\nCreating 1st level contrasts')
            %Kfiles=cell(length(aap.acq_details.subjects),length(K));
            Kfiles={};

            nsub=0;
            for s=1:length(aap.acq_details.subjects)
                fprintf('\n%s',aap.acq_details.subjects(s).megname);
                subdir=fullfile(aap.acq_details.root,aap.acq_details.subjects(s).megname);
                fn=strcat(subdir,filesep,files);
                if ~all(cellfun(@exist,fn))
                    fprintf(' - Files not found!'); continue
                end
                nsub=nsub+1;
                
                %% make (empty?) output directory
                firstlevdir=sprintf('%s/FullFact_1stLev_%s_%s/',subdir, ...
                    regexprep(sensorvoldirs{vd},{'-eeg','_BLOCK1'},{'',''}),fname);
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
                
                if ~isempty(findstr(sensorvoldirs{vd},'SampleK'))
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

            %% prepare directory for 2nd level contrasts
            try cd(secondlevdir);
            catch mkdir(secondlevdir); cd(secondlevdir);
            end

            %% for each event, get 1st level contrast files
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
                        % make explicit mask
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
                        ty=(~isnan(ty) & ty~=0);
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
                        %disp(condir); spm_check_registration(tv.fname);
                        %beep; beep;beep;
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

                    %% add F-contrast
                    if ~isfield(SPM,'xCon') || ~isfield(SPM.xCon,'Vcon') || isempty(SPM.xCon.Vcon) || Overwrite
                        fprintf('\nComputing F contrast')
                        F=eye(numcells);
                        SPM.xCon = spm_FcUtil('Set',Knam{e},'F','c',F,SPM.xX.xKXs);
                        SPM=spm_contrasts(SPM);
                    else
                        fprintf('\nFound F contrast')
                    end
                    
                    if SPM.xY.VY(1).dim(3)>1; continue; end % for now just do stats and plot contours on 2D volumes

                    %% threshold at all levels and save volumes
                    if ~isempty(regexp(suffix,'eog|RMS','ONCE')) 
                        Meth={'none'};
                        getsensors=false;
                    else
                        Meth={'FWE','FDR','none'};
                        getsensors=true;
                    end

                    for m=1:length(Meth)
                        outname=fullfile(SPM.swd,sprintf('0p05_%s.img',Meth{m}));
                        if ~exist(outname,'file') || Overwrite
                            try
                            xSPM = struct('swd', fullfile(SPM.swd,'SPM.mat'), ...
                                'Ic',1,'Im',[],'pm',[],'Ex',[],'title','',...
                                'Mcp',Meth{m},'u',0.05,'k',0);
                            warning off all; % just figure warnings
                            [hReg,xSPM] = csl_getRes(xSPM);
                            warning on all;

                            %if ~isempty(regexp(suffix,'eog|RMS','ONCE'))
                            %else
                            % Set up output volume
                            Vo=struct('fname',outname,...
                                'dim',xSPM.DIM','mat',xSPM.M,'descrip','SPM-filtered',...
                                'dt',[spm_type('uint8') spm_platform('bigend')]);
                            % Reconstruct (filtered) image from XYZ & Z pointlist
                            Y      = zeros(xSPM.DIM(1:3)');
                            OFF    = xSPM.XYZ(1,:) + xSPM.DIM(1)*(xSPM.XYZ(2,:)-1 + xSPM.DIM(2)*(xSPM.XYZ(3,:)-1));
                            Y(OFF) = xSPM.Z.*(xSPM.Z > 0);
                            % Write the reconstructed volume
                            try Vo=rmfield(Vo,{'private','pinfo'}); catch end
                            Vo = spm_write_vol(Vo,Y);
                            %end

                            if ~isempty(regexp(suffix,'eog|RMS','ONCE'))
                                continue
                                %fprintf('\nCheck %s',mfilename); keyboard
                            end

                            if getsensors && ~timefreq

                                if strcmp(suffix,'-mags')||strcmp(suffix,'-grds')||strcmp(suffix,'-lats')||strcmp(suffix,'-longs')
                                    MEG=load('/imaging/local/spm/spm5/EEGtemplates/FIF306_setup.mat');
                                    %% get mags within 'significant' region and
                                    %% voxels associated with this sensor position
                                    Cpos=MEG.Cpos;
                                    Cnames=MEG.Cnames;
                                    locs=102;
                                    Sensors=find(diag(Y(round(Cpos(1,1:locs)*length(Y)),round(Cpos(2,1:locs)*length(Y))))>0);
                                    SensorVox=round(Cpos(:,Sensors)*length(Y));
                                    SensorNames=Cnames(Sensors);

                                    % convert to gradiometer names and indices
                                    % at this location if needed
                                    if ~strcmp(suffix,'-mags')
                                        % NEED TO WORK OUT HOW TO DO THIS!
                                        continue
                                    end
                                    if ~strcmp(suffix,'-lats')

                                    elseif ~strcmp(suffix,'-longs')

                                    elseif ~strcmp(suffix,'-grds')

                                    end

                                elseif strcmp(suffix,'-eeg')
                                    EEG=load('/imaging/local/spm/spm5/EEGtemplates/draft_eeg_montage_djm_mean.mat');
                                    Cpos=EEG.Cpos;
                                    Cnames=EEG.Cnames;
                                    locs=length(Cpos)-2;
                                    Sensors=find(diag(Y(round(Cpos(1,1:locs)*length(Y)),round(Cpos(2,1:locs)*length(Y))))>0);
                                    SensorVox=round(Cpos(:,Sensors)*length(Y));
                                    SensorNames=Cnames(Sensors);
                                else
                                    SensorNames=sprintf('Unknown for %s (see %s)',suffix,mfilename);
                                    fprintf('\nSensor names unknown for %s (see %s)\n',suffix,mfilename);
                                    keyboard;
                                    continue
                                end
                                clear EEG MEG

                                % plot diagnostic
                                clf; colormap jet; imagesc(Y); axis image; hold on
                                plot(Cpos(2,1:locs)*length(Y),Cpos(1,1:locs)*length(Y),'xk')
                                plot(SensorVox(2,:),SensorVox(1,:),'ow')
                                xlabel('Posterior <---> Anterior');
                                ylabel('Right <---> Left');
                                title(suffix);

                                % save
                                save(fullfile(SPM.swd,sprintf('0p05_%s_sensors.mat',Meth{m})), ...
                                    'Sensors','SensorVox','SensorNames');
                            end % find significant sensor locations
                            catch
                                % with tf analysis (FWE?), got error in
                                % csl_getSPM2 at 557:
                                % Q j dimensions non consistent
                            end
                        end % already checked significance?

                        %%
                    end % next threshold
                end % next partition
            end % next set of effects
        end % next file name (e.g. time window)
    end % next sensor type (suffix)
end % next sensorvoldir

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