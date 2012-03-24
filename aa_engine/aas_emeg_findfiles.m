% Automatic Analysis function to retrieve MEG data
% could also consider using aas_getimages or aas_findfiles...

function [files]=aas_emeg_findfiles(aap,filt,subblock)

if length(subblock)==1
    subdir=aas_getsubjpath(aap,subblock{1});
elseif length(subblock)==2
    filt=regexprep(filt,'BLOCK',aap.acq_details.sessions(subblock{2}).name);
    subdir=aas_getsesspath(aap,subblock{1},subblock{2});
end

files=cellstr(spm_select('FPList',subdir,filt));
if exist(files{1},'file')~=2; files(1)=[]; end % could be'/'
    
% ignore any files storing HPI points
files=files(cellfun('isempty',(regexp(files,'HPI.mat'))));

return

