% Prepare the diagnostic folder/name, etc.

function subjname = aas_prepare_diagnostic(aap,subj)

if nargin < 2
    try subj = aap.subj; catch; end;
end

% Save graphical output to common diagnostics directory
if ~exist(fullfile(aap.acq_details.root, 'diagnostics'), 'dir')
    mkdir(fullfile(aap.acq_details.root, 'diagnostics'))
end
if exist('subj', 'var')
    subjname = aap.acq_details.subjects(subj).subjname;
else
    subjname = [];
end

%try figure(spm_figure('FindWin', 'Graphics')); catch; figure(1); end;
%set(gcf,'PaperPositionMode','auto')