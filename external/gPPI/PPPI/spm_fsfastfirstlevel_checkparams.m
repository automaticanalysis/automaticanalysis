function [SPM,params]=spm_fsfastfirstlevel_checkparams(params)
%% This function simply checks the params file/structure to make sure the
%% inputs are correct. For a list of required params, see spm_firstlevel.
% Written by Donald McLaren (mclaren@nmr.mgh.harvard.edu)
% GRECC, Edith Norse Roers Memorial Veterans Hospital, Bedford, MA
% Department of Neurology, Massachusetts General Hospital and Harvard
%   Medical School
% 1/30/2011

%% Function begins here
if exist(params,'file')==2
    params=load(params);
elseif exist('SPM.mat','file')==2
    load SPM.mat
    if exist(fullfile(SPM.swd,'mask.img'),'file') == 2
        params.estimate=0;
        return;
    end
end
if isstruct(params)
else
    error('Params is not a file or a structure.')
end
while numel(fields(params))==1
    F=fieldnames(params);
    params=params.(F{1}); %Ignore coding error flag.
end

%% Set results directory
SPM.swd=[pwd filesep 'SPMana'];

%% Check if SPMmat has been created.
if exist(fullfile(SPM.swd,'SPM.mat'),'file') == 2 && exist(fullfile(SPM.swd,'mask.img'),'file') == 2
     params.estimate=0;
     return;
else
     params.estimate=1;
end

%% Check for nscan
try
    if isnumeric(params.runflac(1).flac.nruns)
        SPM.nscan=zeros(1,params.runflac(1).flac.nruns);
        for jj=1:params.runflac(1).flac.nruns
            SPM.nscan(jj) = size(params.runflac(jj).flac.X,1);
        end
    else
        invokecatchstatement
    end
catch
    error('Number of scans must be defined.')
end

%% Fill in SPM.xX.X
for jj = 1:numel(SPM.nscan)
    try
        tmp={};
        sessdir=fileparts(params.runflac(jj).flac.funcfspec);
        try
            fid = fopen([sessdir filesep params.runflac(jj).flac.parfile],'r');
        catch
            error(['Filename: ' sessdir filesep params.runflac(jj).flac.parfile ' does not exist. Please check filename. Did you move your data?']);
        end
        if fid==-1
           error(['Filename: ' sessdir filesep params.runflac(jj).flac.parfile ' is not readable. Please check filename.']);
        else
            fileline = 0;
            % read ascii-file line by line
            while fileline >= 0
                % read a line
                fileline = fgets(fid);
                if fileline < 0
                    break
                elseif isempty(strtrim(fileline))
                    fline(2)=0;
                    fileline=0;
                else
                    for ii=1:4
                        [fl fileline]=strtok(fileline);
                        fline(ii)=str2double(fl);
                    end
                    fileline=strtrim(fileline);
                end
                if fline(2)>0
                    kk=fline(2); %trial type index in FSFAST starts at 0.
                    try
                        SPM.Sess(jj).U(kk).name{1} = SPM.Sess(jj).U(kk).name{1};
                    catch
                        SPM.Sess(jj).U(kk).name{1}=fileline(~isspace(fileline));
                        %tmp{end+1}=fileline(~isspace(fileline));
                    end
                    if ~isfield(SPM.Sess(jj).U(kk),'ons'); SPM.Sess(jj).U(kk).ons=[];end
                    if ~isfield(SPM.Sess(jj).U(kk),'dur'); SPM.Sess(jj).U(kk).dur=[];end
                    SPM.Sess(jj).U(kk).ons(end+1) = fline(1);
                    SPM.Sess(jj).U(kk).dur(end+1) = fline(3); %% durations are modeled with a duration
                    %Parametric modulators are implemented as conditions in
                    %FSFAST, not as a modulator in SPM. See
                    %(https://surfer.nmr.mgh.harvard.edu/fswiki/FsFastParametricModulation).
                    %This script changes it to be like SPM (see below)
                    SPM.Sess(jj).U(kk).P.name = SPM.Sess(jj).U(kk).name{1};
                    if ~isfield(SPM.Sess(jj).U(kk),'P') || ~isfield(SPM.Sess(jj).U(kk).P,'P'); SPM.Sess(jj).U(kk).P.P=[];end
                    try
                        SPM.Sess(jj).U(kk).P.P(end+1,1) = fline(4);
                        SPM.Sess(jj).U(kk).P.h = 1;
                    catch
                        SPM.Sess(jj).U(kk).P.P(end+1,1) = 1;
                        SPM.Sess(jj).U(kk).P.h = 1;
                    end
                end
            end
        end
    catch
        error(['Timing parameter file: ' sessdir filesep params.runflac(jj).flac.parfile ' is not properly formatted.'])
    end
    for kk=1:numel(SPM.Sess(jj))     
        %Remove empty names
        if isempty(SPM.Sess(jj).U(kk).name{1})
           SPM.Sess(jj).U(kk)=[];
        end
    end
    for kk=1:numel(SPM.Sess(jj).U)
        %Parametric Modulators
        parmod=[];
        if all(SPM.Sess(jj).U(kk).P.P==1)
            SPM.Sess(jj).U(kk).P.P=[];
            SPM.Sess(jj).U(kk).P.h=0;
            SPM.Sess(jj).U(kk).P.name='none';
            SPM.Sess(jj).U(kk).P;
        else
            parmod(end+1)=kk;
        end
    end
    for kk=parmod
        for ii=1:numel(SPM.Sess(jj))
            if all(SPM.Sess(jj).U(ii).ons==SPM.Sess(jj).U(kk).ons) && ii~=kk
                SPM.Sess(jj).U(ii).P.P(:,end+1)=SPM.Sess(jj).U(kk).P.P;
                SPM.Sess(jj).U(kk).P.P=[];
                SPM.Sess(jj).U(ii).P.h=SPM.Sess(jj).U(kk).P.h;
                SPM.Sess(jj).U(kk).P.h=0;
                if strcmp(SPM.Sess(jj).U(kk).P.name{1},'none')
                    SPM.Sess(jj).U(kk).P.name=[];
                end
                SPM.Sess(jj).U(ii).P.name{end+1}=SPM.Sess(jj).U(kk).name{1};
                SPM.Sess(jj).U(kk).P.name{1}='none';
            end
        end
    end
end

for jj = 1:numel(SPM.nscan)
    SPM.xX.K(jj).HParam = Inf; % No highpass filter
end

for jj = 1:numel(SPM.nscan)
    try
        SPM.Sess(jj).C.C=params.runflac(jj).flac.X(:,params.runflac(jj).flac.indnuis);
        SPM.Sess(jj).C.C=SPM.Sess(jj).C.C(:,~all(SPM.Sess(jj).C.C==repmat(mean(SPM.Sess(jj).C.C),size(SPM.Sess(jj).C.C,1),1)));
        for kk=1:size(SPM.Sess(jj).C.C,2)
            SPM.Sess(jj).C.name{kk} = ['nuis' num2str(kk)];
        end
    catch
        SPM.Sess(jj).C.name = cell(0);
        SPM.Sess(jj).C.C = [];
    end
end

%%Define X using FSFAST
% SPM.xX.X=params.X;
% [tmpu,tmporder]=unique(tmp);
% SPM.xX.name=tmp(sort(tmporder));
% SPM.xX.name(numel(SPM.xX.name)+1:size(SPM.xX.X,2))=repmat({'nuis'},1,size(SPM.xX.X,2)-numel(SPM.xX.name));
% SPM.xX.iC=1:numel(tmpu);
% SPM.xX.iB=numel(SPM.xX.iC)+1:1:size(SPM.xX.X,2);
% SPM.xX.iG=[];

%% Check for TR
try
    if isnumeric(params.runflac(1).flac.TR)
        SPM.xY.RT=params.runflac(1).flac.TR;
    else
        invokecatchstatement
    end
catch
    error('TR is not specified.')
end

%% Obtaining proper file structure for input into SPM.xY fields
for jj=1:numel(SPM.nscan)
    params.P{jj}=params.runflac(jj).flac.mri.fspec;
end
try
    if numel(params.P)==numel(SPM.nscan)
        P=[];
        for ii=1:numel(SPM.nscan)
            for kk=1:SPM.nscan(ii) 
                q = [params.P{ii} ',' num2str(kk)];
                P = strvcat(P,q);
            end
            SPM.Sess(ii).row=[sum(SPM.nscan(1:ii-1))+1:1:sum(SPM.nscan(1:ii))];
        end
        SPM.xY.P=P;
        if size(SPM.xY.P,1)~=sum(SPM.nscan)
            invokecatch
        end
    else
        invokecatch
    end
catch
    P=[];
    for ii=1:length(SPM.nscan)
        q = spm_select(SPM.nscan(ii),'image',['Select images for session ' num2str(ii)],{},pwd,'.*',['1:' num2str(SPM.nscan(ii))]);
        P = strvcat(P,q);
        SPM.Sess(ii).row=[sum(SPM.nscan(1:ii-1))+1:1:sum(SPM.nscan(1:ii))];
    end
    SPM.xY.P=P;
end

%% Optional Settings
spm('defaults','FMRI')
global defaults
try 
    Err=0;
    params.olddefs.stats.fmri.fmri_t=spm_get_defaults('stats.fmri.fmri_t');
    if isnumeric(params.MicrotimeRes)
        defaults.stats.fmri.t=params.MicrotimeRes;
    else
        Err=1;
        invokecatchstatement
    end
catch
    if ~Err
        defaults.stats.fmri.t=spm_get_defaults('stats.fmri.fmri_t');
    else
        error('MicrotimeRes not specified correctly.')
    end
end
try
    params.olddefs.stats.fmri.fmri_t0=spm_get_defaults('stats.fmri.fmri_t0');
    if isnumeric(MicrotimeOnset)
        defaults.stats.fmri.t0=params.MicrotimeOnset;
    else
        Err=1;
        invokecatchstatement
    end
catch
    if ~Err
        defaults.stats.fmri.t0=spm_get_defaults('stats.fmri.fmri_t0');
    else
        error('MicrotimeOnset not specified correctly.')
    end   
end

%SET UP BF OPTIONS
SPM.xBF.UNITS = 'secs';
SPM.xBF.dt = params.runflac(1).flac.TR/defaults.stats.fmri.t;
SPM.xBF.T  = defaults.stats.fmri.t;
SPM.xBF.T0 = defaults.stats.fmri.t0;
if params.runflac(1).flac.ana.gammafit==1
    SPM.xBF.name = 'FSFAST Gamma Function';
    %t=0:params.runflac(1).flac.TR:params.runflac(1).flac.ana.timewindow; %
    %Does not match FSFAST, probably because of the u variable in
    %SPM.Sess.U.
    t=0:SPM.xBF.dt:params.runflac(1).flac.ana.timewindow;
    delta=params.runflac(1).flac.ana.gamdelay;
    tau=params.runflac(1).flac.ana.gamtau;
    alpha=params.runflac(1).flac.ana.gamexp;
    if numel(delta)~=1 && numel(tau)~=1 && numel(alpha)~=1; error('Only 1 value of delta. alpha, and tau is allowed.'); end
    h = ( ( ((t - delta)./tau).^alpha) .* exp(-((t - delta)./tau)) );
    h(t<delta) = zeros(1,sum(t<delta));
    SPM.xBF.bf = h/((alpha.^alpha)*exp(-alpha));
    [n m]=size(SPM.xBF.bf);
    if m>n
        SPM.xBF.bf=SPM.xBF.bf';
    end
    SPM.xBF.dt=params.runflac(1).flac.TR/defaults.stats.fmri.t;
    SPM.xBF.length=params.runflac(1).flac.ana.timewindow;
    SPM.xBF.order=1;
    
elseif params.runflac(1).flac.ana.spmhrffit==1
    if params.runflac(1).flac.ana.nspmhrfderiv==0
        SPM.xBF.name='hrf';
    elseif params.runflac(1).flac.ana.nspmhrfderiv==1
        SPM.xBF.name='hrf (with time derivative)';
    elseif params.runflac(1).flac.ana.nspmhrfderiv==2
        SPM.xBF.name='hrf (with time and dispersion derivatives)';
    else
        error('derivatives are not set correctly.')
    end
    SPM.xBF = spm_get_bf(SPM.xBF);
else
    error('firfit is not a valid option at this time.')
end
SPM.xBF.Volterra = 1;  %% This is not modeling volterra interactions

%Restore defaults, these were changed by spm_firstlevel_checkstruct
spm_get_defaults('stats.fmri.fmri_t',params.olddefs.stats.fmri.fmri_t); % Restore old timing
spm_get_defaults('stats.fmri.fmri_t0',params.olddefs.stats.fmri.fmri_t0); % parameters
SPM.xGX.iGXcalc = 'None';
SPM.xVi.form = 'none';
SPM.xGX.sGXcalc = 'mean voxel value';
SPM.xGX.sGMsca =  'session specific';

%Save SPM
try
    save([SPM.swd filesep 'SPM.mat'],'SPM')
catch
    mkdir(SPM.swd)
    save([SPM.swd filesep 'SPM.mat'],'SPM')
end

%Create and estimate the model
%save([SPM.swd filesep 'SPM.mat'],'SPM')
%spm_fmri_spm_ui(SPM);
%spm_spm(SPM);
%load([SPM.swd filesep 'SPM.mat']);
%Setup beta, mask, residual files....
%SPM.Vbeta
%SPM.Vbeta=spm_vol([params.runflac(1).flac.sess filesep 'bold' filesep params.runflac(1).flac.name filesep 'beta.nii']); %local change
%SPM.VResMS
%SPM.VResMS=spm_vol([params.runflac(1).flac.sess filesep 'bold' filesep params.runflac(1).flac.name filesep 'rvar.nii']); %local change
%resss=spm_read_vols(SPM.VResMS);
%resms=resss/params.DOF;
%SPM.VResMS.fname=[SPM.swd filesep 'ResMS.nii'];
%spm_write_vol(SPM.VResMS,resms);
%Mask file
%SPM.xM.T(:) = -Inf;  %% disable threshold masking
%try
%   SPM.xM.VM = spm_vol(params.runflac(1).flac.maskfspec);
%catch
%   try
%       SPM.xM.VM = spm_vol([params.runflac(1).flac.maskfspec '.nii']);
%   catch
%        error('No mask file.');
%  end
%end
%SPM.VM=SPM.xM.VM;
%save([SPM.swd filesep 'SPM.mat'],'SPM')
end