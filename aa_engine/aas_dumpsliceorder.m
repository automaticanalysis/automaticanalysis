% Function to check slice order
% Loads up aap_parameters.mat and goes through all sessions for all
% subjects checking dicom_headers.mat
% field CSASeriesInfo, item MrPhoenixProtocol, ASCCONV block,
% sSliceArray.ucMode parameter
% Rhodri Cusack MRC CBU Cambridge May 2008

function aas_dumpsliceorder(fn)

if (isstruct(fn))
    aap=fn;
else
    if (~exist('fn','var'))
        fn='.';
    end;

    load(fullfile(fn,'aap_parameters.mat'));
end;

for i=1:length(aap.acq_details.subjects)
    fprintf('Subject %s\n',aap.acq_details.subjects(i).mriname);
    for j=1:length(aap.acq_details.sessions)
        sesspth=aas_getsesspath(aap,i,j);
        fprintf('Session %s ',aap.acq_details.sessions(i).name);
        aas_dumpsliceorder_fromdicom(fullfile(sesspth,'dicom_headers.mat'));

    end;
end;

end


    