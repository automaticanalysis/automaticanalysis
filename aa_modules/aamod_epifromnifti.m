% AA module - realignment
% Realignment using SPM5
% i=subject num
% Rhodri Cusack MRC CBU 2004-6 based on original by Matthew Brett
% Modified by Rik Henson 2006-8 to accept reslice "which" option
% 	(plus more defaults can be passed)

function [aap,resp]=aamod_epifromnifti(aap,task,i)

resp='';

switch task
    case 'report'
    case 'doit'
        if (not(iscell(aap.acq_details.subjects(i).seriesnumbers)))
            aas_log(aap,true,'Was exepcting list of filenames in cell array instead of series numbers, check aas_addsubject command in user script');
        end;
        rejectimages=aap.tasklist.currenttask.settings.rejectimages;        
        if ischar(rejectimages)
            rejectimages= str2num(rejectimages);
        end;
        for j=1:length(aap.acq_details.subjects(i).seriesnumbers)
            allfn={};
            if (iscell(aap.acq_details.subjects(i).seriesnumbers{j}))
                % lots of files 3D
                for fn=aap.acq_details.subjects(i).seriesnumbers{j}
                    allfn=[allfn fulfile(aap.directory_conventions.rawdatadir,fn)];
                end;
            else
                %Only one file, assume 4D
                V=spm_vol(fullfile(aap.directory_conventions.rawdatadir,aap.acq_details.subjects(i).seriesnumbers{j}));
                sesspth=aas_getsesspath(aap,i,j);
                aas_makedir(aap,sesspth);
                [pth fle ext]=fileparts(aap.acq_details.subjects(i).seriesnumbers{j});
                for fileind=1:length(V)
                    Y=spm_read_vols(V(fileind));
                    fn=fullfile(sesspth,[sprintf('%s_%04d',fle,fileind) ext]);
                    
                    V(fileind).fname=fn;
                    V(fileind).n=[1 1];
                    spm_write_vol(V(fileind),Y);
                    if (not(any(rejectimages==fileind)))
                        allfn=[allfn fn];
                    end;
                end;
            end;
            % Write out the files, now likely in 3d
            aap=aas_desc_outputs(aap,i,j,'epi',allfn);
            
        end;
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;














