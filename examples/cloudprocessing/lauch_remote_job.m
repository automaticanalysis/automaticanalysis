aa_ver4_devel

fn=dir('/mridata/cbu/*MR10025');

import camneuro.aacc.*;
if (~exist('cc','var') || isempty(cc))
    cc=Aacc();
end;
if (~cc.isSessionAlive())
    cc.login();
end;

cc.bucketname='imaging-analysis';
cc.rawdicombucketname='imaging-rawdata';

% DEFINE STUDY SPECIFIC PARAMETERS
aap.options.aa_minver=4.0; % will only work on aa version 4.0 or above
% This is an example of a basic analysis
aap=aarecipe('aap_parameters_defaults.xml','aap_tasklist_typical_fmri.xml');

aap=aas_localconfig(aap);

% Directory for analysed data
aap.directory_conventions.analysisid='analysis_asn_077';

aap.tasksettings.aamod_slicetiming.sliceorder=[1:32];

% The name of your scanner
aap.directory_conventions.rawdatadir='Machine_MRC-CBU_MRC35119_TrioTim';


% upload whatever has not be uploaded, or whatever has changed
for subjind=1:length(fn)
    visit_fn=dir(fullfile('/mridata/cbu',fn(subjind).name,'2010*'));
    sess_fn=dir(fullfile('/mridata/cbu',fn(subjind).name,visit_fn(1).name,'*MEPI*'));
    allsess=[];
    for sessind=1:length(sess_fn)
        allsess(sessind)=str2num(sess_fn(sessind).name(8:10));
    end;
    
    % Add subject & sessions
    if (~isempty(allsess))
        aap=aas_addsubject(aap,sprintf('Patient_%s',fullfile(fn(subjind).name,visit_fn(1).name)),{allsess});
    end;
    % SET ANY OTHER PARAMETERS YOU WOULD LIKE TO BE DIFFERENT FROM THE DEFAULTS
    
    
end;

aap=aas_addsession(aap,'movie');

save('/imaging/rhodri/camneuro2/asn_test/aap_movie.mat','aap');
cc.uploadAap('/imaging/rhodri/camneuro2/asn_test/aap_movie.mat',aap.directory_conventions.analysisid)
% DO PROCESSING
%        aa_doprocessing_camneuro(aap,'rhodri','rhodricusacklori','camneuroasneu','activestatenetworks')
