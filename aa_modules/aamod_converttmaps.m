% AA module - Converts real time t maps to NIFTI format
% For all real time t-maps identified, convert them and put them in directory specified by
% aap.directory_conventions.tmapsdirname
% Rhodri Cusack MRC CBU Cambridge Nov 2005

function [aap,resp]=aamod_converttmaps(aap,task,i)

resp='';

switch task
    case 'domain'
        resp='subject';   % this module needs to be run once per subject
        
    case 'description'
        resp='Convert real time t statistic outputs';
        
    case 'summary'
        if (length(aap.acq_details.subjects(i).tmaps)==0)
            resp=sprintf('No t maps for subject %s\n',aap.acq_details.subjects(i).mriname);
        else
            resp=sprintf('Converted t maps for subject %s to %s\n', aap.acq_details.subjects(i).mriname,aap.directory_conventions.centralstore_structurals);
        end;
	
    case 'inandout'
  	aap.tasklist.inputs={};
	anoutput.processingstream='realtime';
	anoutput.prefix='f';
	aap.tasklist.outputs={anoutput}; 

    case 'report'
    case 'doit'
        subjpath=aas_getsubjpath(aap,i);
        tmapsdir=fullfile(subjpath,aap.directory_conventions.tmapsdirname);
        aas_makedir(aap,tmapsdir);
        
        for j=1:length(aap.acq_details.subjects(i).tmaps)
            seriesnum=aap.acq_details.subjects(i).tmaps(j);
            outputpath=fullfile(aas_getsubjpath(aap,i),aap.directory_conventions.tmapsdirname);
            [aap]=aas_convertseries(aap,i,seriesnum,outputpath);
            [pth fle ext]=fileparts(outputpath);
            [pth fle ext]=fileparts(pth);
            mapfn=dir(fullfile(outputpath,sprintf('f%s-%04d*nii',fle,seriesnum)));    
            if (length(mapfn)>0)
                Vtmap=spm_vol(fullfile(outputpath,mapfn.name));
                Ytmap=spm_read_vols(Vtmap);
                Ytmap=(Ytmap-2048)/256; %Rescale!
                spm_write_vol(Vtmap,Ytmap);
            end;
        end;
        
                
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
