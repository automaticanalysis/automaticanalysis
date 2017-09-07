function outputdirs=aas_emeg_ANOVA_2ndLevel(aap,S,Overwrite)
% Work in progress. Run 2nd level of 2- or 3-way repeated measures
% ANOVA with partitioned variance, using S from aas_emeg_ANOVA_1stLevel
% Currently somewhat specific to my experiment but potentially
% generalisable...
% S has fields...
% .files - anonymised file names of component images
% .Knam - cell array of names of each 1st level effect
% .P - cell array of partitions of each Knam, for second level
% .Kfiles - file names of first level contrasts
% .secondlevdir - where this struct is saved and 2nd level analysis will go
% .datatype - 'sensor, 'tf', or 'source'
% .nsub - number of subjects for whom first level contrasts were computed
% .PPM - whether to run PPMs
% djm 04/08/08, 09/06/09
% See Rik's chapter

cd(aap.acq_details.root)
addpath /imaging/dm01/MoreTools

% add marsbar paths
warning off all
root='/imaging/local/spm/marsbar/marsbar-0.40/';
[files dirs]=spm_select('List',root,'');
for d=1:size(dirs,1); addpath(strcat(root,dirs(d,:)),'-END'); end
warning on all

global defaults

fprintf('\nStarting 2nd level across %g subjects:',S.nsub)

try PPM=S.PPM; catch; PPM=false; end

%% for each event, get 1st level contrast files
outputdirs={};
for e=1:length(S.Knam) % effects of interest
    for kp=1:size(S.P{e},1); % each level of a 'separated' factor
        ind=repmat(S.P{e}(kp,:)',[size(S.Kfiles,1) 1]);
        files=char(S.Kfiles{:,e});
        files=files(ind==1,:);
        numcells=size(files,1)/size(S.Kfiles,1);
        condir=fullfile(S.secondlevdir,sprintf('%s_%g',S.Knam{e},kp));
        if ~isempty(regexp(condir,';-0')); debugnow; end %%%
        if ~exist(fullfile(condir,'ResMS.img'),'file') || Overwrite || exist(fullfile(condir,'ResI_0001.img'),'file')
            fprintf('\nCreating SPM for %s #%g',S.Knam{e},kp);
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
                ty(floor(end/2):end,:,:)=0;
            end
            ty=(~isnan(ty) & ty~=0); % mask any 0s or NaNs
            if strcmpi(S.datatype,'source')
                %By=spm_read_vols(spm_vol('/imaging/local/spm/spm5/apriori/brainmask.nii'));
                By=spm_read_vols(spm_vol('/imaging/local/spm/spm5/EEGtemplates/bsBinary_Template_Cortex_Mesh.img'));
                ty=ty.*By;
            elseif strcmpi(S.datatype,'tf')
                % ugly method!
                Ds=spm_select('FPList',fullfile(subdir,aap.acq_details.sessions(1).name), ...
                    regexprep(S.Kfiles{1,1},{'.*FullFact_1stLev_','_average-.*','(N?ST?t\d)'},{'.*','.*mat','$1.*'}));
                load(deblank(Ds(1,:)));
                ty=ty.*D.tf.mask';
                if size(ty,1)>300 %%% HACK FOR MY WORKING MEMORY DATA TO REDUCE SEARCH VOLUME
                    ty(1:89,:)=0; % mask out times prior to 200ms
                    ty(352:end,:)=0; % mask out times following 1250ms
                end
            elseif strcmpi(S.datatype,'erf')
                if size(ty,3)>300 %%% HACK FOR MY WORKING MEMORY DATA TO REDUCE SEARCH VOLUME
                    ty(:,:,1:89)=0; % mask out times prior to 200ms
                    ty(:,:,352:end)=0; % mask out times following 1250ms
                end
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
            defaults.modality='EEG'; % otherwise non-sphericity correction might complain about no significant voxels?
            SPM = spm_spm(SPM);
            save('SPM.mat','SPM')
        else
            fprintf('\nLoading SPM for %s #%g',S.Knam{e},kp);
            load(fullfile(condir,'SPM.mat'));
        end

        % add PPM
        if PPM
            cd(condir)
            if Overwrite || ~exist('Cbeta_0001.img','file') || ~isfield(SPM,'PPM')
                fprintf('\nCalculating PPM for %s...\n',pwd)
                SPM=spm_spm_Bayes(SPM);
            else fprintf('\nFound previously estimated PPM')
            end
        end

        %% initialise contrast fields
        indF=0; indPosT=0; indNegT=0; indF_Bayes=0; indPosT_Bayes=0; indNegT_Bayes=0;
        if ~isfield(SPM,'xCon') || Overwrite || 1==1
            SPM.xCon = [];
            SPM.PPM.xCon=[];
        else % check for exisitng contrasts
            for tc=1:length(SPM.xCon)
                if SPM.xCon(tc).STAT=='F' && strcmp(SPM.xCon(tc).name,[regexprep(S.Knam{e},'^(Effect_|)','Classical F test: ') sprintf('_%g',kp)]); indF=tc; end
                if SPM.xCon(tc).STAT=='T' && strcmp(SPM.xCon(tc).name,[regexprep(S.Knam{e},'^(Effect_|)','Classical +ve T: ') ,sprintf('_%g',kp)]); indPosT=tc; end
                if SPM.xCon(tc).STAT=='T' && strcmp(SPM.xCon(tc).name,[regexprep(S.Knam{e},'^(Effect_|)','Classical -ve T: '),sprintf('_%g',kp)]); indNegT=tc; end
                %if SPM.xCon(tc).STAT=='P' && strcmp(SPM.xCon(tc).name,[regexprep(S.Knam{e},'^(Effect_|)','Bayesian F test: '),sprintf('_%g',kp)]); indF_Bayes=tc; end
                if SPM.xCon(tc).STAT=='P' && strcmp(SPM.xCon(tc).name,[regexprep(S.Knam{e},'^(Effect_|)','Bayesian +ve T: '),sprintf('_%g',kp)]); indPosT_Bayes=tc; end
                if SPM.xCon(tc).STAT=='P' && strcmp(SPM.xCon(tc).name,[regexprep(S.Knam{e},'^(Effect_|)','Bayesian -ve T: '),sprintf('_%g',kp)]); indNegT_Bayes=tc; end
            end
        end

        % define contrasts
        cons2calc=[];
        if ~indF
            fprintf('\nDefining classical F test')
            F=eye(numcells);
            if isempty(SPM.xCon)
                SPM.xCon = spm_FcUtil('Set',[regexprep(S.Knam{e},'^(Effect_|)','Classical F test: '),sprintf('_%g',kp)],'F','c',F,SPM.xX.xKXs);
            else
                SPM.xCon(end+1) = spm_FcUtil('Set',[regexprep(S.Knam{e},'^(Effect_|)','Classical F test: '),sprintf('_%g',kp)],'F','c',F,SPM.xX.xKXs);
            end
            cons2calc(end+1)=length(SPM.xCon);
        else fprintf('\nFound classical F test')
        end
%         if PPM && ~indF_Bayes
%             debugnow
%         else fprintf('\nFound Bayesian F test')
%         end
        if numcells==1
            if ~indPosT
                fprintf('\nDefining classical positive T test')
                SPM.xCon(end+1) = spm_FcUtil('Set',[regexprep(S.Knam{e},'^(Effect_|)','Classical +ve T: '),sprintf('_%g',kp)],'T','c',1,SPM.xX.xKXs);
                cons2calc(end+1)=length(SPM.xCon);
            else fprintf('\nFound classical positive T test')
            end
            if ~indNegT
                fprintf('\nDefining classical negative T test')
                SPM.xCon(end+1) = spm_FcUtil('Set',[regexprep(S.Knam{e},'^(Effect_|)','Classical -ve T: '),sprintf('_%g',kp)],'T','c',-1,SPM.xX.xKXs);
                cons2calc(end+1)=length(SPM.xCon);
            else fprintf('\nFound classical negative T test')
            end
            if PPM && ~indPosT_Bayes
                fprintf('\nDefining Bayesian +ve T test')
                SPM.xCon(end+1) = spm_FcUtil('Set',[regexprep(S.Knam{e},'^(Effect_|)','Bayesian +ve T: '),sprintf('_%g',kp)],'P','c',1,SPM.xX.xKXs);
                SPM.xCon(end).eidf=0;
                SPM.PPM.xCon(length(SPM.xCon)).PSTAT='T';
                cons2calc(end+1)=length(SPM.xCon);
            else fprintf('\nFound Bayesian +ve T test')
            end
            if PPM && ~indNegT_Bayes
                fprintf('\nDefining Bayesian -ve T test')
                SPM.xCon(end+1) = spm_FcUtil('Set',[regexprep(S.Knam{e},'^(Effect_|)','Bayesian -ve T: '),sprintf('_%g',kp)],'P','c',-1,SPM.xX.xKXs);
                SPM.xCon(end).eidf=0;
                SPM.PPM.xCon(length(SPM.xCon)).PSTAT='T';
                cons2calc(end+1)=length(SPM.xCon);
            else fprintf('\nFound Bayesian -ve T test')
            end
        end

        % calculate contrasts
        if ~isempty(cons2calc);
            if PPM 
%                 for pc=1:length(SPM.PPM.xCon)
%                     SPM.PPM.xCon(pc).PSTAT='T'; % not sure why PPM changes according to whether all PPM contrasts are 'T'???
%                 end
                [SPM.PPM.xCon(:).PSTAT]=deal('T');
            end
            fprintf('\nCalculating contrasts:\n');
            SPM=spm_contrasts(SPM,cons2calc);
        end

        outputdirs=[outputdirs; condir];

    end % next partition
end % next set of effects

return
