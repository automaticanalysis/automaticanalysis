function [SPM,params]=spm_firstlevel_checkparams(params)
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
end
if isstruct(params)
else
    error('Params is not a file or a structure.')
end
while numel(fields(params))==1
    F=fieldnames(params);
    params=params.(F{1}); %Ignore coding error flag.
end

%% Check for nscan
try
    if isnumeric(params.nscan)
        SPM.nscan = params.nscan;
    else
        invokecatchstatement
    end
catch
    error('Number of scans must be defined.')
end

%% Check for TR
try
    if isnumeric(params.TR)
        SPM.xY.RT=params.TR;
    else
        invokecatchstatement
    end
catch
    error('TR is not specified.')
end

%% Obtaining proper file structure for input into SPM.xY fields
try
    params.P
    if numel(params.P)==1 && (numel(params.P)~=numel(params.nscan) || exist(params.P,'file')~=2)
        tmp = dir(params.P);
        path=fileparts(params.P);
        if isempty(path)
            for ii=1:numel(tmp)
                params.P{ii,1} = [pwd filesep tmp(ii).name];
            end
        else
           for ii=1:numel(tmp)
                params.P{ii,1} = [pwd filesep path filesep tmp(ii).name];
            end
        end
    end
    if numel(params.P)==sum(params.nscan)
        if isempty(strfind(params.P{1},','))
           for ii=1:numel(params.P)
               q = [params.P{ii} ',' num2str(ii)];
               P = strvcat(P,q);
           end            
        end
    elseif numel(params.P)==numel(params.nscan)
        P=[];
        for ii=1:numel(params.nscan)
            for kk=1:SPM.nscan(ii) 
                q = [params.P{ii} ',' num2str(kk)];
                P = strvcat(P,q);
            end
        end
        SPM.xY.P=P;
        if numel(SPM.xY.P)~=sum(params.nscan)
            invokecatch
        end
    else
        invokecatch
    end
catch
    P=[];
    for ii=1:length(params.nscan)
        q = spm_select(params.nscan(ii),'image',['Select images for session ' num2str(ii)],{},pwd,'.*',['1:' num2str(params.nscan(ii))]);
        P = strvcat(P,q);
    end
    SPM.xY.P=P;
end

%% Set results directory
try
    params.dir
catch
    filepath=fileparts(char(SPM.xY.P(1,:)));
    params.dir=[filepath filesep 'results'];
    clear tmp filepath
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
try
    if strcmp('secs',params.Units) || strcmp('scans',params.Units)
        SPM.xBF.UNITS = params.Units;
    else
        Err=1;
        invokecatchstatement
    end
catch
    if ~Err
        SPM.xBF.UNITS = 'secs'; % default setting
    else
       error('Units not specified correctly.')
    end  
end
SPM.xBF.dt = params.TR/defaults.stats.fmri.t;
SPM.xBF.T     = defaults.stats.fmri.t;
SPM.xBF.T0    = defaults.stats.fmri.t0;
try
    if strcmp('hrf',params.basis) || ...
               strcmp('hrf (with time derivative)',params.basis) || ...
               strcmp('hrf (with time and dispersion derivatives)',params.basis) || ...
               strcmp('Fourier set',params.basis) || ...
               strcmp('Fourier set (Hanning)',params.basis) || ...
               strcmp('Gamma functions',params.basis) || ...
               strcmp('Finite Impulse Response',params.basis)
        SPM.xBF.name = params.basis;
    else
        Err=1;
        invokecatchstatement
    end
catch
    if ~Err
        SPM.xBF.name = 'hrf'; % default
    else
        error('Basis function is not specified correctly.')
    end 
end
SPM.xBF = spm_get_bf(SPM.xBF);
try
    if params.Volterra==1 || params.Volterra==2
        SPM.xBF.Volterra = params.Volterra;
    else
        Err=1;
        invokecatchstatement
    end
catch
    if ~Err
        SPM.xBF.Volterra = 1;  %% This is not modeling volterra interactions
    else
        error('Volterra was not set to 1 or 2')
    end
end
%Restore defaults, these were changed by spm_firstlevel_checkstruct
spm_get_defaults('stats.fmri.fmri_t',params.olddefs.stats.fmri.fmri_t); % Restore old timing
spm_get_defaults('stats.fmri.fmri_t0',params.olddefs.stats.fmri.fmri_t0); % parameters
try
    if strcmp('Scaling',params.iGXcalc) || strcmp('Scaling',params.iGXcalc)
        SPM.xGX.iGXcalc = params.iGXcalc;
    else
        
    end
catch
    if ~Err=1
        SPM.xGX.iGXcalc = 'None';
    else
        error()
    end
end
try
    SPM.xVi.form = params.xVi;
catch
    SPM.xVi.form = 'AR(1)';
end

try
    if isumeric(params.HP)
    else
       Err=1;
       invokecatchstatement
    end
catch
    if ~Err
        params.HP=128; % Default setting
    else
        error()
    end
end

params.mask

try 
    params.onsets(jj)
    
catch
    if Err
        error('onsets, names, or durations not specified correctly.')
    end
end
try
    params.regressors(jj).files;
catch
end

try
        params.regressors(jj).names;
catch
end
