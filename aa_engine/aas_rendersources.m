function aas_rendersources(S)
%% render and print thresholded results
% Input structure S with fields:
% .condir - location of SPM.mat
% .Ic - contrast indices or 'F', 'T' or 'P'
% .thresh - threshold e.g. 0.025 for classical stats or 0.95 for Bayesian
% .MCCs - cell array of multiple comparison options e.g. {'None','FDR'}
% .brain - full path of mesh on which to render
% .style - NaN=old, 1=normal, <1=brighter, 0.75="slightly", 0.25="lots"
% .title - string; if empty results in default contrast title
% .overwrite - whether or not to overwrite
% .outfile - optional path of destination file
%
% Uses Ferath's modifications: csl_getRes
%
% djm 03/03/09, 10/06/09
%% fill in defaults
try S.Ic; catch S.Ic=1; end
try S.MCCs; catch; S.MCCs={'None','FDR'}; end
try S.brain; catch; S.brain='/imaging/local/spm/spm5/rend/render_single_subj.mat'; end
try S.style; catch; S.style=1; end
try S.overwrite; catch; S.overwrite=false; end
try S.title; catch; S.title='~~~'; end

%% check output directory
try outfile=S.outfile;
catch
    % if using aaMEG directory structure...
    outfile=fullfile(S.condir,'../../figures',[regexprep(S.condir,{'.nii','.img','.*2ndLev_','/'},{'','','','_'}) '.ps']);
end
[pth nam ext]=fileparts(outfile);
if ~exist(pth,'dir') % ... otherwise save to condir
    outfile=fullfile(condir,[nam ext]);
end
outfile=spm_select('Cpath',outfile);

%% find contrast indices
tails=1;
if ischar(S.Ic) 
    try oSPM=load(fullfile(S.condir,'SPM.mat'));
    catch; fprintf('\nFailed to load 2nd level SPM.mat file.\n'); debugnow
    end
    Ic=strmatch(S.Ic,char(oSPM.SPM.xCon.STAT));
    if isempty(Ic)
        fprintf('\nFound no contrasts of type "%s".\n',S.Ic); return;
    elseif length(Ic)>2; 
        fprintf('\nFound >2 contrasts of type "%s".\n',S.Ic); debugnow
    elseif length(Ic)==2 && strcmpi(S.Ic,'F')
        fprintf('\nFound >1 contrasts of type "F".\n'); debugnow
    elseif length(Ic)==2
        try
            if length([oSPM.SPM.xCon(Ic).c])==2 % put -ve tail first
                [c ord]=sort([oSPM.SPM.xCon(Ic).c]);
                Ic=Ic(ord);
            end
            if oSPM.SPM.xCon(Ic(1)).c~=-oSPM.SPM.xCon(Ic(2)).c; gocatch; end
            tails=2;
        catch fprintf('\n2 contrasts of type "%s" seem unbalanced.\n',S.Ic); debugnow
        end
    end
    switch S.Ic
        case 'T', 
            regexprep(outfile,'Effect','SPMT')
        case 'F', 
            regexprep(outfile,'Effect','SPMF')
        case 'P', 
            regexprep(outfile,'Effect','PPM')
    end
else Ic=S.Ic;
end    

%% check whether already done
if ~isempty(regexp(outfile,';-0')); debugnow; end %%%
if exist(outfile,'file') && ~S.overwrite
    fprintf('\nAlready rendered results to %s',outfile);
    return
end

%% get ready
fprintf('\nThresholding and rendering results to:\n %s',outfile);

h = spm_figure('GetWin','Graphics');
set(h,'renderer','painters');
colormap gray

if ~exist('oSMP','var')
    try oSPM=load(fullfile(S.condir,'SPM.mat'));
    catch; fprintf('\nFailed to load 2nd level SPM.mat file.\n'); debugnow
    end
end

warning off all; delete(outfile); warning on all

%% get thresholded results
ok2rend=[];
for ic=Ic(:)'
    if oSPM.SPM.xCon(ic).STAT=='P';
        try S.thresh; catch; S.thresh=0.99; end
        mccs={'none'};
    else
        try S.thresh; catch; S.thresh=0.05/tails; end
        mccs=S.MCCs;
    end
    for mcc=mccs
        xSPM = struct( ...
            'swd', S.condir, ... % directory containing SPM.mat file
            'Ic', ic, ... % no of contrast (or contrasts for conjunction)
            'Im', [],... % no of contrast to mask with. Empty for no masking
            'pm', [],... % masking contrast uncorrected p
            'Ex', [],... % whether masking is inclusive or exclusive
            'title', S.title,... % if empty results in default contrast title
            'thresDesc', mcc,... % Mutiple comp method: FWE|FDR|none
            'u',S.thresh,... % threshold (corrected or uncorrected, as above)
            'k', 0, ...  % extent threshold
            'auto',true); % to return Z for PPMs
        warning off all; % just figure warnings
        [hReg,xSPM] = aas_getRes(xSPM);
        warning on all;
        ok2rend=[ok2rend size(xSPM.XYZ,2)>1]; % rendering sometimes fails for <2 voxels ??

        if exist('dat','var')
            dat(end+1) = struct( 'XYZ', xSPM.XYZ,'t', xSPM.Z','mat', xSPM.M,'dim', xSPM.DIM);
        else dat(1) = struct( 'XYZ', xSPM.XYZ,'t', xSPM.Z','mat', xSPM.M,'dim', xSPM.DIM);
        end
    end
end

%% choose colours and labels for rendering
if tails==2 && length(mccs)==2 && oSPM.SPM.xCon(ic).STAT~='P';
    colors=[0 0 1; 0 1 0; 1 0 0; 0 1 0];
    txt={sprintf('t test: p<%g (each tail). Warm=+ve, cold=-ve; darker=%s, lighter=%s.', ...
        S.thresh,regexprep(mccs{1},'None','unc'),regexprep(mccs{2},'None','unc')), ...
        datestr(now,0)};
elseif tails==2 && oSPM.SPM.xCon(ic).STAT~='P';
    colors=[0 0 1; 1 0 0; 0 1 0];
    txt={sprintf('t test: p<%g (each tail; %s). Warm=+ve, cold=-ve',S.thresh,mcc),...
        datestr(now,0)};
elseif tails==2 && oSPM.SPM.xCon(ic).STAT=='P';
    colors=[0 0 1; 1 0 0; 0 1 0];
    txt={sprintf('PPM p>%g; effect size = %g; Warm=+ve, cold=-ve',...
        S.thresh, oSPM.SPM.xCon(ic).eidf),...
        datestr(now,0)};
else colors=[1 0 0; 0 0 1; 0 1 0];
    colnames={'Red','Blue','Green'}
    for cc=1:length(Ic)
        if oSPM.SPM.xCon(Ic(cc)).STAT=='P'
            statstr=sprintf('PPM (effect>%g), p>%g',oSPM.SPM.xCon(Ic(cc)).eidf,S.thresh);
        else
            statstr=sprintf('SPM %s, p<%g',oSPM.SPM.xCon(Ic(cc)).STAT,S.thresh);
        end
        txt{cc}=sprintf('%s: %s; %s',colnames{cc}, oSPM.SPM.xCon(Ic(cc)).name,statstr);
    end
    txt{end+1}=datestr(now,0);
end

%% render
if any(ok2rend) % do the rendering at all thresholds
    spm_render_dm(dat(logical(ok2rend)),S.style,S.brain,colors(logical(ok2rend),:));
    ann1=annotation('textbox',[0 .94 1 .06],'Color','r','String',txt,'edge','none');
else
    % ...ok, probably no significant activation;
    % show contrast values instead...
    try spm_check_registration(oSPM.SPM.xY.P);
        Txt='Showing 1st level contrasts from each subject.';
    catch % for factor with many levels there will be too many 1st level contrasts
        spm_check_registration(spm_select('List',S.condir,'(beta|ess).*img'))
        Txt='Showing beta and ess images.';
    end
    pos=[-40, -60, 40];
    spm_orthviews('reposition',pos);
    txt{end+1}=sprintf('No significant voxels. %s',Txt);
    txt{end+1}=sprintf('Crosshairs at %g, %g, %g.',round(pos(1)),round(pos(2)),round(pos(3)));
    ann1=annotation('textbox',[0 .94 1 .06],'Color','r','edge','none','String', txt);
end

%% print
% png -r82 slightly better than jpeg90; no difference of renderer on either format
% for ps, painters is best, but little effect of resolution so keep default
spm_print(outfile);
print('-dpng','-r82',strrep(outfile,'.ps','.png'));
try delete(ann1); catch end;

%% some old options
%     if mirrored % could fill in complement to make nicer figure?
%         V=spm_vol(xSPM.Vspm.fname);
%         Y=spm_read_vols(V);
%         Y=Y-flipdim(Y,1);
%         spm_write_vol(V,Y);
%     end
%
%     Could create hi-res version for nicer rendering in MRIcro?
%     template='/imaging/dm01/MyStructurals/MeanOf12_Segment&Normalise/wc1djm_MeanOf12.nii'; % Hi res!
%     Vtemplate=spm_vol(template);
%     try spm_check_orientations([Vtemplate; xSPM.Vspm.fname])
%     catch
%         fprintf('\nReslicing files')
%         reslice to template
%         clear flags
%         flags.mask=0;
%         flags.mean=0;
%         flags.which=1; % all but 1st
%         spm_reslice(char(template,xSPM.Vspm.fname),flags)
%     end

return