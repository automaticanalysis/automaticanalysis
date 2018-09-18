function [aap,resp] = aamod_binarizeimage(aap, task, subjind)
% AAMOD_BINARIZEIMAGE Create binary image from a continuous one
%
% [aap,resp] = AAMOD_BINARIZEIMAGE(aap, task, subjind)
%
% input stream:    (continuous image, e.g. native_grey)
% output stream:   (binarized image, e.g. binarized_native_grey)
%
% The threshold for binarizing is in:
%
% aap.tasklist.currenttask.settings.thresh

resp='';

switch task
    case 'domain'
        resp='subject';
    case 'description'
        resp='Create binary image from a continuous valued image';
    case 'doit'
        
        thresh = aap.tasklist.currenttask.settings.thresh;
        %inputlabel = aap.tasklist.currenttask.settings.inputlabel;
        inputlabel = aap.tasklist.currenttask.inputstreams(1).stream{1};
        outputlabel = sprintf('%s_binarized', inputlabel);
        
        
        inputImg = aas_getfiles_bystream(aap, subjind, inputlabel);
        
        [pth, nm, ext] = fileparts(inputImg);
        
        V = spm_vol(inputImg);
        Vout = V;                    % create copy V for output
        Vout.descrip = sprintf('Binarized using aamod_binarizeimage using threshold of %g', thresh);
        
        % Read in labeled brain
        [Y, xyz] = spm_read_vols(V);
        
        % Binarize data based on labels: 1 if it is in list, 0 otherwise
        Y = Y > thresh;
        
        % write out
        Vout.fname = fullfile(pth, [nm '_' outputlabel ext]);
        spm_write_vol(Vout, Y);
        
        
        % describe output
        aap = aas_desc_outputs(aap, subjind, outputlabel, Vout.fname);
        
end