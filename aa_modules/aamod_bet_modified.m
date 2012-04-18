% AA module
% Runs BET on structural


function [aap,resp]=aamod_bet(aap,task,i)

resp='';

switch task        
    case 'summary'
        subjpath=aas_getsubjpath(i);
        resp=sprintf('Align %s\n',subjpath);
        
    case 'report'
        
    case 'doit'
        
        warning off
        
        tic
        
        % Let us use the native space...
        Sfn = aas_getfiles_bystream(aap,i,'structural');
        
        % Which file is considered, as determined by the structural parameter!
        if size(Simg,1) > 1
            Simg = deblank(Simg(aap.tasklist.currenttask.settings.structural, :));
            fprintf('\tWARNING: Several structurals found, considering: %s\n', Simg)
        end
        
        [pth nme ext]=fileparts(Sfn);
        
        outStruct=fullfile(pth,['bet_' nme ext]);
        % Run BET [-R Using robust setting to avoid neck!]   
        fprintf('Initial BET pass (recursive) to find optimal centre of gravity\n')
        [junk, w]=aas_runfslcommand(aap, ...
            sprintf('bet %s %s -f %f -v -R',Sfn,outStruct, ...
            aap.tasklist.currenttask.settings.bet_f_parameter));
        
        % This outputs last centre of gravity from recursive command...
        indxS = strfind(w, 'c-of-g');
        indxS = indxS(end) + 7;
        indxE = strfind(w(indxS:end), 'mm');
        indxE = indxE(1) - 3;
        COG = w(indxS:indxS+indxE);
        
        fprintf('Second BET pass extracting also brain masks\n')
        % Run BET [-A Now let's get the brain masks and meshes!!]   
        [junk, w]=aas_runfslcommand(aap, ...
            sprintf('bet %s %s -f %f -c %s -v -A',Sfn,outStruct, ...
            aap.tasklist.currenttask.settings.bet_f_parameter, COG)...
            );
        
        
        % Find the mask images
        D = dir(fullfile(pth, 'bet*mask*'));
        outMask = '';
        for d = 1:length(D)
            outMask = strvcat(outMask, fullfile(pth, D(d).name));
        end

        % Find the meshes
        D = dir(fullfile(pth, 'bet*mesh*'));
        outMesh = '';
        for d = 1:length(D)
            outMesh = strvcat(outMesh, fullfile(pth, D(d).name));
        end
       
        % DESCRIBE OUTPUTS!
        % Structural image after BETting
        aap=aas_desc_outputs(aap,i,'structural',outStruct);
        aap=aas_desc_outputs(aap,i,'BETmask',outMask);
        aap=aas_desc_outputs(aap,i,'BETmesh',outMesh);
        
        time_elapsed
end
