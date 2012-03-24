function [aap resp]=aamod_emeg_groupsourceanalysis(varargin)
% Danny Mitchell 13/03/08

%% check task settings, subject, block etc
[aap subblock doit resp settings]=aa_emeg_checktasksettings(mfilename('fullpath'),varargin);  
if ~doit; return; end


%% get files to average
fprintf('\nFinding files from %g blocks of %g subjects', ...
    length(aap.acq_details.sessions),length(aap.acq_details.subjects))

%smoothsourcevolfilter='^sw_ace.*_\w+inv(1|2)w+_\d+\.nii$';
%smoothsourcevolfilter='^sw_mace.*inv(1|2).*_\d+\.nii$';
smoothsourcevolfilter='^.*ms.*Hz.*;.*\.nii$';

contrasts=[];cblocks=[];
for s=1:length(aap.acq_details.subjects)
    fprintf('.')
    for b=aap.acq_details.selected_sessions
        block=aap.acq_details.sessions(b).name;
        subdir=fullfile(aap.acq_details.root, ...
            aap.acq_details.subjects(s).megname,block,'sourcevols');

        % get all smoothed source contrasts
        temp=spm_select('List',subdir,smoothsourcevolfilter); % '^sw_.*Hz_\d*_\d\d\.nii$'
        contrasts=char(contrasts,temp);
        cblocks=char(cblocks,repmat(block,[size(temp,1) 1]));

    end % next block
end % next subject

try blank=find(contrasts(:,1)==' ');
catch
    aas_log(aap,true,sprintf('Found no data! Using filter: %s',smoothsourcevolfilter));
end
contrasts(blank,:)=[]; 
cblocks(blank,:)=[];

% get unique names
if isempty(contrasts); 
    aas_log(aap,true,sprintf('Found no data! Using filter: %s',smoothsourcevolfilter));
else
    [cnames I]=unique(contrasts,'rows');
    cblocks=cblocks(I,:);
end

%% fill (contrast x subject) matrix of full file paths
fprintf('\nFilling contrast x subject matrix');
cfiles={};
for s=1:length(aap.acq_details.subjects)
    fprintf('.')
    for c=1:size(cnames,1)
        %ocfiles(ismember(cnames,cellstr(contrasts{s})),s)=cellstr([sortrows(cblocks{s}) sortrows(contrasts{s})]);
        filename=deblank(fullfile(aap.acq_details.root, ...
            aap.acq_details.subjects(s).megname, cblocks(c,:), ...
            'sourcevols', cnames(c,:)));
        if exist(filename,'file'); cfiles(c,s)=cellstr(filename);
        else cfiles(c,s)={' '};
        end
    end
end
clear cblocks contrasts

%% confirm/create group analysis directory
cwd=pwd;
try cd(fullfile(aap.acq_details.root,'GroupAnalysis_Sources'));
catch
    mkdir(fullfile(aap.acq_details.root,'GroupAnalysis_Sources'));
    cd(fullfile(aap.acq_details.root,'GroupAnalysis_Sources'));
end
if ~exist(fullfile(aap.acq_details.root,'GroupAnalysis_Sources','figures'),'dir')
    mkdir(fullfile(aap.acq_details.root,'GroupAnalysis_Sources','figures'))
end
%% do it
% zuloc/zul: zscored unsigned localisation (of contrast) (just +ve tail)
% slocs: signed localisation of contrast (no longer used)
% scoul: signed contrast of unsigned localisation (both tails)
% zucosl: zscored rectified contrast of signed localisation (just +ve tail)

spm_defaults;
global defaults UFp
defaults.modality='FMRI';
UFp=0.001;

% try
%     load(efiles{1,1});
%     % vert seem to corespond for all subjects;
%     % Is (and data) sorted for each subject in source_plots prior to saving data,
%     % so data for all subjects should now be in same order as vert.
%     % Is=D.inv{1}.inverse.Is; % indices of vertices, different across subs.
%     vert=D.inv{1}.mesh.tess_mni.vert; % MNI coordinates of vertices (mm)
% catch
% end
sMRIfile = fullfile(spm('dir'),'templates','T2.nii');
Vin   = spm_vol(sMRIfile);

if strcmpi(settings.Stats,'auto')
    if length(aap.acq_details.subjects)>13; settings.Stats='SPM';
    else settings.Stats='SnPM';
    end
end

if strcmpi(settings.Stats,'SnPM')
    try addpath /imaging/local/spm/snpm5b/; snpm_defaults;
    catch
        aas_log(aap,true,'\n Failed to initialise SnPM software. \n')
    end
end

% Run RFX one-sample t-tests for each tf-contrast image
for c=1:size(cnames,1)
    clear SPM
    SPM.xY.P = char(cfiles{c,:}); % contrast files from each subject
    SPM.xY.P(all(SPM.xY.P==' ',2),:)=[]; % remove missing files
    if size(SPM.xY.P,1)<2;  % only continue if there are multiple subjects
        fprintf('\nFailed to find multiple subjects for contrast %s.\n', cnames(c,:))
        keyboard
        continue;
    elseif size(SPM.xY.P,1)~=length(aap.acq_details.subjects);
        fprintf('\nFound %g instead of %g subjects for contrast %s.\n', ...
            size(SPM.xY.P,1),length(aap.acq_details.subjects),cnames(c,:))
        disp(SPM.xY.P)
        %continue;
    end

    cname=deblank(cnames(c,:)); cname=cname(1:end-4);
    try cd(fullfile(aap.acq_details.root,'GroupAnalysis_Sources',cname));
    catch
        mkdir(fullfile(aap.acq_details.root,'GroupAnalysis_Sources',cname));
        cd(fullfile(aap.acq_details.root,'GroupAnalysis_Sources',cname));
    end

    if exist('SPM.mat','file')==0 || settings.Overwrite==1
        fprintf('\nRunning group stats on %s',cname)
        try SPM.xY.VY= spm_vol(SPM.xY.P);
        catch 
%             fprintf('\nFailed to load all of following files. Please debug.\n')
%             disp(SPM.xY.P); keyboard
        end
        SPM.nscan = size(SPM.xY.P,1); % number of subjects

        if strcmpi(settings.Stats,'SPM')
            % use standard parametric stats of SPM
            defaults.modality='FMRI'; % Some problems with the Results otherwise?
            % Assemble rest of SPM structure
            SPM.xX = struct('X',ones(size(SPM.xY.P,1),1),...
                'iH',1,'iC',zeros(1,0),'iB',zeros(1,0),'iG',zeros(1,0),...
                'name',{{'mean'}},'I',[(1:size(SPM.xY.P,1))' ones(size(SPM.xY.P,1),3)],...
                'sF',{{'obs'  ''  ''  ''}});
            SPM.xC = [];
            SPM.xGX = struct(...
                'iGXcalc',1,	'sGXcalc','omit',				'rg',[],...
                'iGMsca',9,	'sGMsca','<no grand Mean scaling>',...
                'GM',0,		'gSF',ones(size(SPM.xY.P,1),1),...
                'iGC',	12,	'sGC',	'(redundant: not doing AnCova)',	'gc',[],...
                'iGloNorm',9,	'sGloNorm','<no global normalisation>');
            SPM.xVi	= struct('iid',1,'V',speye(size(SPM.xY.P,1)));
            
            explicitmask='/imaging/dm01/MEG/aaMEG/GMmask.nii'; % perhaps better for PPMs?
            explicitmask='/imaging/local/spm/spm5/apriori/brainmask.nii';
            
            Mdes 	= struct(	'Analysis_threshold',	{'None (-Inf)'},...
                'Implicit_masking',	{'Yes: NaNs treated as missing'},...
                'Explicit_masking',	{sprintf('Yes: %s',explicitmask)});
            SPM.xM	= struct('T',-Inf,'TH',ones(size(SPM.xY.P,1)*2,1)*-Inf,...
                'I',1,'VM',spm_vol(explicitmask),'xs',Mdes);
            Pdes 	= {{'1 condition, +0 covariate, +0 block, +0 nuisance'; '1 total, having 1 degrees of freedom'; sprintf('leaving %g degrees of freedom from %g images',size(SPM.xY.P,1)-1,size(SPM.xY.P,1))}};
            SPM.xsDes = struct(	'Design',		{'One sample t-test'},...
                'Global_calculation',	{'omit'},...
                'Grand_mean_scaling',	{'<no grand Mean scaling>'},...
                'Global_normalisation',	{'<no global normalisation>'},...
                'Parameters',		Pdes);

            % Estimate parameters
            %===========================================================================
            spm_unlink(fullfile('.', 'mask.img')); % avoid overwrite dialog
            SPM = spm_spm(SPM);

            % setup contrasts (-ve tail 1st, so +ve stays in glass brain
            % when rendering)
            try
                SPM.xCon = spm_FcUtil('Set',['-' cname],'T','c',-1,SPM.xX.xKXs);
                SPM.xCon(end+1) = spm_FcUtil('Set',cname,'T','c',1,SPM.xX.xKXs);
                SPM=spm_contrasts(SPM);
%                 if SPM.xVol.DIM(3)==2
%                     % plot Vcon using example D structure
%                     % this is the mean localisation
%                     Dname=regexprep(regexprep(cfiles{c,2},'_timecourses.*','.mat'),'sw_','');
%                     source_plots(Dname,SPM.xCon(1).Vcon.fname);
%                 end
            catch
                % possibly no inmask voxels?
                fprintf('\nFailed to run contrast. Please debug.\n')
                keyboard
            end

            % transfer description from first image to spms
            for ind=1:length(SPM.xCon)
                SPM.xCon(ind).Vspm.descrip=SPM.xY.VY(1).descrip;
                SPM.xCon(ind).Vcon.descrip=SPM.xY.VY(1).descrip;
                save(fullfile(SPM.swd,'SPM.mat'),'SPM')
            end
            % render and print later...
        elseif strcmpi(settings.Stats,'SnPM')
            % use SnPM
            snpmP.fnames=SPM.xY.P;
            snpmP.vFWHM=10; % in mm for smoothing of local variance; 10 is recommended; 0 for real rather than pseudo-t
            try snpmP.iGloNorm=settings.GloNorm;
            catch snpmP.iGloNorm=1; % Global normalisation: 1=none; 2=proportional; 3=AnCova
            end
            snpmP.iGMsca=2; % no grand mean scaling
            snpmP.iTHRESH=1; % 1=no threshold masking; 2=proportional; 3=absolute
            % create design and config file in pwd
            snpm_ui_dm(snpmP);
            % perform computation
            snpm_cp(pwd);
            % transfer description from first image
        else
            error('Unknown value for field "Stats"')
        end
    end % already done?
end; % next contrast

% prepare to render and print significant results
if strcmpi(settings.MCC,'none')
    MCCs={'None','FDR'};
else MCCs={'None',settings.MCC};
end
Style=1; % Rendering Style: NaN=old, 1=normal, <1=brighter, 0.75="slightly", 0.25="lots"
Brain='/imaging/local/spm/spm5/rend/render_smooth_average.mat'; %'/imaging/local/spm/spm5/rend/render_single_subj.mat'
h = spm_figure('GetWin','Graphics');
set(h,'renderer','painters');
colormap gray
        
for c=1:size(cnames,1)
    cname=deblank(cnames(c,:)); cname=cname(1:end-4);
    if ~exist(fullfile(aap.acq_details.root,'GroupAnalysis_Sources','figures',[cname '.png']),'file')||settings.Overwrite==1
        if ~isempty(strfind(cname,'_timecourses_')); continue; end
        fprintf('\nRendering results for %s',cname);
        condir=strrep(deblank(cnames(c,:)),'.nii','');
        if strcmpi(settings.Stats,'SnPM')
            ppP.CWD=fullfile(aap.acq_details.root,'GroupAnalysis_Sources',condir);
            ppP.bNeg=2; % negative then positive
            ppP.MCC=settings.MCC;
            ppP.alpha=settings.thresh
            cd(ppP.CWD);
            [dat lab]=snpm_pp_dm(ppP);
        else
            ok2rend=[];
            clear dat
            for i=1:2 % -ve then +ve tails
                for mcctype=1:length(MCCs)
                xSPM = struct( ...
                    'swd', fullfile(aap.acq_details.root,'GroupAnalysis_Sources',condir,'SPM.mat'), ... % full path to SPM.mat file
                    'Ic', 1, ... % no of contrast (or contrasts for conjunction)
                    'Im', [],... % no of contrast to mask with. Empty for no masking
                    'pm', [],... % masking contrast uncorrected p
                    'Ex', [],... % whether masking is inclusive or exclusive
                    'title', '',... % if empty results in default contrast title
                    'Mcp', MCCs{mcctype},... % Mutiple comp method: FWE|FDR|none
                    'u',settings.thresh/2,... % threshold (corrected or uncorrected, as above)
                    'k', 0); % extent threshold
                try oSPM=load(xSPM.swd);
                catch
%                     fprintf('\nFailed to load 2nd level SPM.mat file. Please debug.\n')
%                     keyboard
                end
                xSPM.Ic=i; %Set the current contrast Index
                    warning off all; % just figure warnings
                    [hReg,xSPM] = csl_getRes(xSPM);
                    warning on all;
                    ok2rend=[ok2rend size(xSPM.XYZ,2)>1]; % rendering seems to fail for <2 voxels
%                     if xSPM.DIM(3)==2
%                         % remap dipoles to mni voxel coordinates
%                         % collapse across time & laterality, though this is retained in glassbrain view
%                         d=xSPM.XYZ(1,:); % significant dipole numbers
%                         p=vert(d,:)'; % mm coordinates of these dipoles. NB don't use Is.
%                         pv=inv(Vin.mat)*[p; ones(1,size(p,2))]; % their voxel coordinates
%                         xSPM.XYZ=pv(1:3,:);
%                         xSPM.M=Vin.mat;
%                         xSPM.dim=Vin.dim;
%                     end
                    if exist('dat','var')
                        dat(end+1) = struct( 'XYZ', xSPM.XYZ,'t', xSPM.Z','mat', xSPM.M,'dim', xSPM.DIM);
                    else dat(1) = struct( 'XYZ', xSPM.XYZ,'t', xSPM.Z','mat', xSPM.M,'dim', xSPM.DIM);
                    end
                end % next threshold
                end % next tail
            try lab=xSPM.Vspm.descrip; catch; lab='Description missing'; end
        end
        
        % do the rendering
        colors=[0 0 1; 0 1 0; 1 0 0; 0 1 0];
            if any(ok2rend) % do the rendering at all thresholds
                spm_render_dm(dat(logical(ok2rend)),Style,Brain,colors(logical(ok2rend),:));
                txt={sprintf('p<%g (each tail). Warm=+ve, cold=-ve; darker=uncorrected, lighter=FDR.',settings.thresh/2), ...
                    sprintf('Printed:...%s',datestr(now,0))};
                ann1=annotation('textbox',[0 .94 1 .06],'Color','r','String',txt,'edge','none');
            else
                    % probably no significant activation; try to show
                    % single subjects instead...
                    if ~isempty(strfind(cname,'_timecourses_')); continue; end
                    spm_check_registration(oSPM.SPM.xY.P);
                    % move to centre of ROI if one is specified
                    Dname=regexprep(deblank(oSPM.SPM.xY.P(1,:)), ...
                        {'/sourcevols','_\d+_\d+','sw_','_-?\d*--?\d*ms.*','\..*$'}, ...
                        {'','','','.mat','.mat'});
                    load(Dname);
                    if isfield(D,'ROIs')
                        if strcmpi('MNI',D.ROIs.space{1}); pos=D.ROIs.co(1,1:3);
                        else pos=tal2mni(D.ROIs.co(1,1:3));
                        end
                        spm_orthviews('reposition',pos);
                    end
                    ann1=annotation('textbox',[0 .94 1 .06],'Color','r','edge','none','String', ...
                        {sprintf('%s: No significant voxels (even uncorrected) at p<%g per tail.',lab,settings.thresh/2,settings.Stats), ...
                        sprintf('Showing single subjects. Crosshairs at %g, %g, %g.   Date printed: %s',round(pos(1)),round(pos(2)),round(pos(3)),datestr(now,0))});
            end

        % png -r82 slightly better than jpeg90; no difference of renderer on either format
        % for ps, painters is best, but little effect of resolution so keep default
        spm_print(fullfile(aap.acq_details.root,'GroupAnalysis_Sources','figures','SPMs.ps'));
        print('-dpng','-r82',fullfile(aap.acq_details.root,'GroupAnalysis_Sources','figures',[cname '.png']));
        try delete(ann1); catch end;
    end
end

return % ROI bit needs fixing for new file name format

%%
% Source level ROI analysis using Marsbar:
% Print each group of timewindows/ROIs to a separate page
% Within each group, plot different localisation methods and signages
% on different rows (as the scale of the signals is not comparable)
cgroups=unique(regexprep(cellstr(cnames),{'.*timecourses.*','(.*_[^_]+_[^_]+_[^_]+_[^_]+)_\d+.nii','T_[^_]+_'},{'1','$1','T_.*'}));
% one for each epoch, sensor type, contrast window and inversion

warning off all;
spm_figure('GetWin','Graphics'); clf;
% add paths
root='/imaging/local/spm/marsbar/marsbar-0.40/';
[files dirs]=spm_select('List',root,'');
for d=1:size(dirs,1); addpath(strcat(root,dirs(d,:)),'-END'); end
warning on all

for g=1:size(cgroups,1)
    clear roi
    % for each group, get all blocks and event contrasts
    filt=sprintf('%s_\\d*.nii',cgroups{g}); %sprintf('.*_%s_.*_%s_\\d*.nii',cgroups{g},inversions{i})
    good=regexp(cellstr(cnames),filt);
    timecourses=regexp(cellstr(cnames),'_timecourses_');
    subgroup=cnames(~cellfun('isempty',good) & cellfun('isempty',timecourses),:);
    outfile=fullfile(aap.acq_details.root,'GroupAnalysis_Sources','figures',[regexprep(cgroups{g},'\.\*','') '_ROIs.png']);
    if exist(outfile,'file') && settings.Overwrite==0
        continue
    end

    data=[]; clf
    clear scalefactor
    for c=1:size(subgroup,1)
        try
            if ~exist('roi','var') % get ROI centres and event/contrast name from D file
                t=1;
                while isempty(findstr(deblank(subgroup(c,:)),cfiles{t,1}));
                    t=t+1;
                end
                Dname=regexprep(cfiles{t,1},'sw_','');
                enum=regexp(Dname,'.*Hz_\w+_(\d+).*','tokens');
                try enum=str2double(enum{1}{1});
                catch
                    keyboard
                end
                block=regexp(Dname,'.*/(.*)/[^/]*/[^/]*','tokens');
                block=block{1}{1};
                Dname=regexprep(Dname,{'/sourcevols','_-?\d*--?\d*ms.*'},{'','.mat'});
                load('-mat',Dname);
                try ename=D.inv{1}.contrast{1}.names{enum}; % hmm, names only in 1st inversion?
                catch ename='missing';
                end

                if ~isfield(D,'ROIs')
                    warning off all
                    [Numeric,Txt]=xlsread(aap.acq_details.subjects(1).ContrastFile);
                    warning on all
                    D.ROIs.space=Txt(5,enum+2);
                    D.ROIs.co=str2num(char(Txt(6,strmatch(regexprep(ename,'\(.*\)',''),Txt(1,:)))));
                end

                if strcmpi(D.ROIs.space{1},'Tal'); roi=tal2mni(D.ROIs.co);
                else roi=D.ROIs.co;
                end
            end
            skip=false;
        catch
            fprintf('\nFailed to find ROI for %s!',subgroup(c,:))
            skip=true;
            continue
        end
        try
        cname=deblank(subgroup(c,:)); cname=cname(1:end-4);
        spm_name = fullfile(aap.acq_details.root,'GroupAnalysis_Sources',cname,'SPM.mat');
        load(spm_name)
        marsD  = mardo(SPM); % Make marsbar design object
        for r=1:size(roi,1)
            marsR = maroi_sphere(struct('centre',roi(r,:),'radius',10));  % create 10mm radius ROI
            marsY  = get_marsy(marsR, marsD, 'mean'); % Fetch data into marsbar data object
            summarydata=summary_data(marsY);
            % there might be missing subjects for some contrasts but not
            % others! Fill others with NaNs and use nanmean later...
            summarydata=[summarydata; nan(size(data,2)-length(summarydata),1)];
            
            %%% hack!! 
            if ~isempty(findstr(SPM.xCon(1).Vcon.descrip,'scoul')) && any(strcmp(cname(end-1:end),{'06','18','19','20','21'}))
               % scaled differently to other contrasts, so rescale to match max of 1st contrast
               if ~exist('scalefactor','var')
                   scalefactor=max(abs(summarydata))/max(max(abs(data(:,:,r))));
               end
               summarydata=summarydata/scalefactor;
            end
            %%%
            
            data(c,:,r)=summarydata;
        end
        catch
            fprintf('\nFailed to load 2nd level model for %s!',subgroup(c,:))
            skip=true;
            continue
        end
    end

    if skip; continue; end
    
    % get block names
    tmp=squeeze(struct2cell(aap.acq_details.sessions));
    tmp=[repmat('.*(',[length(tmp) 1]) char(tmp) repmat(').*',[length(tmp) 1])];
    blocks=regexprep(cellstr(subgroup),cellstr(tmp),'$1');

    % plot and save
    fid=fopen(strrep(outfile,'.png','.xls'),'w');
    for r=1:size(roi,1)
        ta=subplot(size(roi,1),1,r);
        b=bar(nanmean(data(:,:,r),2),0.8);box off;hold on
        set(b,'FaceColor',[1 1 1]); set(b,'LineWidth',2);
        plot(repmat((1:size(data,1))',1,size(data,2)),data(:,:,r),'.','color',[.7 .7 .7]);
        %set(ta,'xcolor',[1 1 1])
        set(ta,'color','none','xtickLabel',[])
        axis tight;
        if r==1
            title({sprintf('%s:',cgroups{g}), 'Data from 10mm radius ROIs, with 95%% CIs.'});
            % add rotated tick labels
            ylims=ylim;
            yoff=ylims(1)-(ylims(2)-ylims(1))/100;
            try labs=D.inv{1}.contrast{1}.names;
            catch labs=cellstr(repmat('missing',size(subgroup,1),1));
            end
            if isempty(labs);labs=cellstr(repmat('missing',size(subgroup,1),1));end
            for t=1:size(data,1)
                try text(t,yoff,labs{t},'rotation',90,'horizontalalignment','right');
                catch text(t,yoff,labs{t-length(labs)},'rotation',90,'horizontalalignment','right');
                end
            end
        end

        roilabel=sprintf('MNI: %.0f, %.0f, %.0f',roi(r,1),roi(r,2),roi(r,3));
        ylabel(roilabel);
        fprintf(fid,'\n\n%s:\n',roilabel);
        for b=1:length(blocks); fprintf(fid,'%s\t',blocks{b}); end
        fprintf(fid,'\n');
        for b=1:length(unique(blocks));
            for l=1:length(labs); fprintf(fid,'%s\t',labs{l}); end
        end
        fprintf(fid,'\n');
        for s=1:size(data,2);
            fprintf(fid,'%g\t',data(:,s,r));
            fprintf(fid,'\n');
        end

        % set(ta,'ytick',0);
        events=unique(regexprep(cellstr(subgroup),'.*_\w+_(\d+)(_U|_S|)\..*','$1'));
        set(ta,'xtick',1:c)
        % add error bars
        [h,p,ci]=ttest(data(:,:,r)');
        for b=1:length(h)
            if h(b)==1; plot([b b],ci(:,b),'r-','linewidth',2)
            else plot([b b],ci(:,b),'b-','linewidth',2)
            end
        end
        % separate blocks
        divs=length(events):length(events):size(subgroup,1);
        plot(repmat(divs,2,1)+0.5,ylim,'k-')
        drawnow
    end % next ROI on next row

    fclose(fid);

    % label blocks
    blocks=blocks(1:length(events):size(subgroup,1));
    pos=get(ta,'position');
    for b=1:length(blocks)
        annotation('textbox',[pos(1)+(b-1)*pos(3)/length(blocks) pos(2)-0.04 pos(3)/length(blocks) 0.02] ...
            ,'string',blocks{b},'horizontalalignment','center','edgecolor','none');
    end;
    % print
    spm_print(strrep(outfile,'.png','.ps')); % does something to avoid segmentation fault?
    print('-dpng', '-zbuffer','-r82', outfile); % low res, but filts in web browser
end

cd(cwd)

return