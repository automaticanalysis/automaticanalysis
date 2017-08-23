% This is a tempate for a module code processing an MRI session

function [aap,resp]=aamod_ASL(aap,task,subj,sess)
resp='';

switch task
    case 'report'
        
    case 'doit'
        
        % obtain input filename(s)
        hdrs = load(aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,[subj sess],'ASL_dicom_header')); hdrs = hdrs.dcmhdr;
        imgs = cellstr(aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,[subj sess],'ASL'));

        
        % do stuff
        subseries = unique(spm_file(imgs,'path'))';
        for s = 1:numel(subseries)

            if ~isempty(strfind(hdrs{s}.ImageType,'ORIGINAL')), subsess = 'raw';
            elseif ~isempty(strfind(hdrs{s}.ImageType,'DERIVED'))
                if ~isempty(strfind(hdrs{s}.ImageType,'RELCBF')), subsess = 'CBF'; 
                else subsess = 'PWI'; end
            end
            
            inputfnames = char(imgs(cell_index(imgs,subseries{s})));
            
            switch subsess
                case 'CBF'
                    sc = aas_getsetting(aap,'scaleCBF',s);
                    V = spm_vol(inputfnames);
                    dat = spm_read_vols(V);
                    dat = dat  * sc;
                    outputfnames = spm_file(inputfnames,'prefix','sc_');
                    V = V(1);
                    V.dt = spm_type('float32');
                    nifti_write(outputfnames,dat,'scCBF',V);
                otherwise
                    outputfnames = inputfnames;
            end            
            
            % Describe outputs
            aap=aas_desc_outputs(aap,aap.tasklist.currenttask.domain,[subj sess],['ASL-' subsess],outputfnames);
        end
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end