function [aap,resp] = aamod_binaryfromlabels(aap, task, subjind)
% AAMOD_BINARYFROMLABELS Create binary GM image from a labeled brain.
%
% [aap,resp] = AAMOD_BINARYFROMLABELS(aap, task, subjind)
%
% input stream:     labeled_brain
% output stream:   labeled_(descriptor) [probably "grey"]
%
% The label numbers probably integers) that are included in the mask are
% specified in settings.labels (i.e.
% aap.tasklist.currenttask.settings.labels).

resp='';

switch task
    case 'domain'
        resp='subject';
    case 'description'
        resp='Create binary mask from a labeled brain'
    case 'doit'

        
        labels = aap.tasklist.currenttask.settings.labels;     % what numbers get combined
        outlabel = aap.tasklist.currenttask.settings.outputlabel; % label (e.g., 'grey')
        
        labelImg = aas_getfiles_bystream(aap, subjind, 'labeled_brain');
        
        [pth, nm, ext] = fileparts(labelImg);
        
        V = spm_vol(labelImg);       
        Vout = V;                    % create copy V for output
        Vout.descrip = 'Binarized using aamod_binaryfromlabels based on labeled brain';
        
        % Read in labeled brain
        [Y, xyz] = spm_read_vols(V);
        
        % Binarize data based on labels: 1 if it is in list, 0 otherwise
        Y = ismember(Y, labels);
        
        % write out
        Vout.fname = fullfile(pth, [nm '_' outlabel ext]);
        spm_write_vol(Vout, Y);


        % describe output
        aap = aas_desc_outputs(aap, subjind, sprintf('labeled_%s', outlabel), Vout.fname);

end