% Automatic analysis - check for missing scans
% A problem with data transfer once caused a scan to be missing in the raw
% DICOM data. You may use this function to scan all of the raw data
% directories in a study for missing scans.
% Rhodri Cusack MRC CBU Cambridge 2004

function [aap]=aa_checkformissingscans(studyroot)

if (~exist('studyroot','var'))
    studyroot=pwd;
end;


% First, load AAP structure
aaploadfn=fullfile(studyroot,'aap_parameters');
load(aaploadfn);

wasmissing=false;
for i=1:length(aap.acq_details.subjects)
    fprintf('Checking %s\n',aap.acq_details.subjects{i});
    for j=aap.acq_details.selected_sessions
        fles=aas_getimages(aap,i,j,'f');
        for k=1:size(fles,1)
            pos=findstr('-',fles(k,:));
            ind=str2num(fles(k,(pos(end-1)+1):(pos(end)-1)));
            if (k~=1)
                if (ind~=(oldind+1))
                    fprintf('Missing in session %s scan %d\n',aap.acq_details.sessions{j},ind);
                    wasmissing=true;
                end;
            end;
            oldind=ind;
        end;
    end;
end;
        

if (~wasmissing)
    fprintf('No scans missing\n');
end;