% AA module - fieldmap from NIFTI

function [aap,resp]=aamod_fieldmapfromnifti(aap,task,subj)

resp='';

switch task
    case 'report'
    case 'doit'
        forder = {...
            'rawmag'...
            'rawmag'...
            'rawphase'...
            };
       
        sesspth=fullfile(aas_getsubjpath(aap,subj),aap.directory_conventions.fieldmapsdirname);
        aas_makedir(aap,sesspth);
               
        if ~iscell(aap.acq_details.subjects(subj).fieldmaps)
            aas_log(aap,true,'Was exepcting list of filenames in cell array');
        end;
        niftistruct = aap.acq_details.subjects(subj).fieldmaps{1};
        niftifile = niftistruct.fname; % checks first only
        hdrfile = niftistruct.hdr;
        if ~exist(niftifile{1},'file')
            niftisearchpth=aas_findvol(aap,'');
            if ~isempty(niftisearchpth)
                niftifile = spm_file(niftifile,'path',niftisearchpth);
                hdrfile = spm_file(hdrfile,'path',niftisearchpth);
            end
        end
        
        % images
        for f = 1:numel(niftifile)
            comp = strcmp(spm_file(niftifile{1},'Ext'),'gz');
            if comp
                gunzip(niftifile{f});
                niftifile{f} = niftifile{f}(1:end-3);
            end
            
            aas_makedir(aap,fullfile(sesspth,forder{f}));
            V=spm_vol(niftifile{f});
            Y=spm_read_vols(V(1));
            V.fname = spm_file(niftifile{f},'path',fullfile(sesspth,forder{f}));
            spm_write_vol(V,Y);
            fn{f} = V.fname;
            
            if comp, delete(niftifile{f}); end
        end
        aap=aas_desc_outputs(aap,subj,'fieldmap',char(fn));
        
        % header
        hdr = loadjson(hdrfile);
        dcmhdr{1}.EchoTime = hdr.EchoTime;
        dcmhdrfn=fullfile(aas_getsubjpath(aap,subj),'fieldmap_dicom_header.mat');
        save(dcmhdrfn,'dcmhdr');
        aap=aas_desc_outputs(aap,subj,'fieldmap_dicom_header',dcmhdrfn);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end