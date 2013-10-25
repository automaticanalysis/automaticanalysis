% AA module - run melodic tensor ICA from
%

function [aap,resp]=aamod_tensor_ica(aap,task)

resp='';

switch task
    case 'report'
        
    case 'doit'
        
        
        % Get inputs
        allfn = '';
        allstructfn='';
        
        % Get TR from first subject
        fn = aas_getfiles_bystream(aap,1,1,'epi_dicom_header');
        D = load(fn);
        TR = D.DICOMHEADERS{1}.RepetitionTime;
        
        for i=1:length(aap.acq_details.subjects)
% Set up EPI for fsf file
            % set feat_files(1) "/home/vjohn/imaging/vivek_first_preprocesswithspm/aamod_3dto4d_00001/CBU110254_MR10033_CC310008/movie/warfMR10033_CC310008-0010-00_4d"
            % Run on first selected sessions
            subjfn=aas_getfiles_bystream(aap, i, aap.acq_details.selected_sessions(1), 'epi');            
            allfn=[allfn sprintf('set feat_files(%d) "%s"\n',i,subjfn)];

% Set up structrual for fsf file
% set highres_files(1) "/home/vjohn/imaging/vivek_first_preprocesswithspm/aamod_bet_modified_00001/CBU110254_MR10033_CC310008/structurals/bet_mmsMR10033_CC310008-0002-00001-000192-01"

            structfn=aas_getfiles_bystream(aap,i,'structural');
            allstructfn=[allstructfn sprintf('set highres_files(%d) "%s"\n',i,structfn)];
            
            % Check TR same in all other subjects
            fn = aas_getfiles_bystream(aap,i,1,'epi_dicom_header');
            D = load(fn);                       
            if TR~=D.DICOMHEADERS{1}.RepetitionTime;
                aas_log(aap,true,'found different TR in subjects');
            end
        end
        
        % Number of volumes
        [s w]=aas_runfslcommand(aap,sprintf('fslinfo %s',subjfn));
        if (s)
            aas_log(aap,true,sprintf('Error getting number of volumes with fslinfo %s',w));
        else
            pos=strfind(w,[10 'dim4']);
            pos2=strfind(w((pos(1)+5):end),10);
            totalvolumes=w((pos+5):(pos+pos2(1)+4));
            totalvolumes=str2double(totalvolumes);
        end;
        
        % Now run this fsf
        studypth=aas_getstudypath(aap);
        fsffn=fullfile(studypth,'melodic.fsf');
        
        parms=[];
        parms.totalvolumes=totalvolumes;
        parms.TR=TR/1000;  % ms to seconds
        parms.feat_file_list=allfn;
        parms.numsubjects=length(aap.acq_details.subjects);
        parms.highres_files=allstructfn;
        parms.outputdirectory=[fullfile(studypth,sprintf('melodicoutput_%s', datestr(now,30))) '/'];
        
        aap=aas_template(aap,'aamod_tensorica_template.fsf', parms, fsffn);
        
        [s w]=aas_runfslcommand(aap,sprintf('feat %s',fsffn));
        if (s)
           aas_log(aap,true,sprintf('Error from FEAT while trying to run melodic %s',w));
        end;
end;
