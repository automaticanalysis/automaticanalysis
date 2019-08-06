% AA module - Converts fieldmaps maps to NIFTI format
% Rhodri Cusack MRC CBU Cambridge Nov 2005

function [aap,resp]=aamod_convert_fieldmaps(aap,task,subj,sess)

resp='';

switch task
    case 'report'
    case 'doit'
        domain = aap.tasklist.currenttask.domain;
        
        [aap, convertedfns, dcmhdr]=aas_convertseries_fromstream(aap,domain,[subj,sess],'dicom_fieldmap');

        sesspath=aas_getsesspath(aap,subj,sess);
        
        % Restructure outputs
        if iscell(convertedfns{1}), convertedfns = vertcat(convertedfns{:}); end % multiseries
        outstream = char(convertedfns{:});
        if iscell(dcmhdr{1}), dcmhdr = horzcat(dcmhdr{:}); end % multiple series
        
        V = spm_vol(outstream);
        toMask = false;
        if numel(V) ~= numel(dcmhdr) 
            if numel(dcmhdr) == 1 && isfield(dcmhdr{1}, 'PerFrameFunctionalGroupsSequence') % Philips new stacked (precalc fieldmap)
                H = cell2mat(dcmhdr{1}.PerFrameFunctionalGroupsSequence);
                pos = cell2mat([H.FrameContentSequence]);
                pos = reshape([pos.DimensionIndexValues],4,[])';
                H = H(pos(:,2) == 1); % first slices
                ftype = cell2mat([H.MRImageFrameTypeSequence]);
                Vfm = V(cell_index({ftype.FrameType},'FIELD_MAP'));
                Vmag = V(cell_index({ftype.FrameType},'\M')); Vmag = Vmag(1);

                Y = spm_read_vols(Vmag);
                Vmag.fname = fullfile(sesspath,aap.directory_conventions.fieldmapsdirname,'mag.nii'); 
                Vmag.n = [1 1];
                spm_write_vol(Vmag,Y);
                aas_runfslcommand(aap,sprintf('bet2 %s %s -m -n -f 0.7',Vmag.fname,spm_file(Vmag.fname,'ext','')))
                Ymag = spm_read_vols(spm_vol(spm_file(Vmag.fname,'suffix','_mask')));
                toMask = true;

            else
                aas_log(aap,true,sprintf('ERROR: Unhandled mismatch between volumes (n=%d) and header (n=%d).',numel(V),numel(dcmhdr)));
            end
        end
                
        outstream = cell(0,0);
        for v = 1:numel(Vfm)
            Y = spm_read_vols(Vfm(v));
            
            if toMask, Y = Y.*Ymag; end
            
            Vfm(v).fname = fullfile(sesspath,aap.directory_conventions.fieldmapsdirname,sprintf('serie%02d.nii',v)); 
            Vfm(v).n = [1 1];
            spm_write_vol(Vfm(v),Y);
            outstream(v) = cellstr(Vfm(v).fname);
        end        
        aap=aas_desc_outputs(aap,domain,[subj,sess],'fieldmap',outstream);

        dcmhdrfn=fullfile(sesspath,'fieldmap_dicom_header.mat');
        save(dcmhdrfn,'dcmhdr');
        aap=aas_desc_outputs(aap,domain,[subj,sess],'fieldmap_dicom_header',dcmhdrfn);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,true,sprintf('Unknown task %s',task));
end