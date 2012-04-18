% AA module
% Runs MVPA classification using the Donders Machine Learning Toolbox
% by Dr Marcel van Gerven
% Cognitive Artificial Intelligence, Room B.02.04
% Donders Institute for Brain, Cognition and Behaviour
% Tel: +31 (0) 24 3615606
% m.vangerven@donders.ru.nl
%
% This is WORK IN PROGRESS and thus UNFINISHED!

function [aap,resp]=aamod_decodeDMLT(aap,task,p)

resp='';

switch task
    case 'domain'
        resp='subject';  % this module needs to be run once per subject
        
    case 'description'
        resp='SPM5 align';
        
    case 'summary'
        subjpath=aas_getsubjpath(p);
        resp=sprintf('Align %s\n',subjpath);
        
    case 'report'
        
    case 'doit'
                
        
        %% PLAN
        % 1) Get the raw data (or betas - usually)
        % 2) Extract the trials x voxels data
        % 3) Get the column of conditions from the aas_addevent settings
        % 4) Choose which type of analysis to do
        % 5) Train/classify (MAT file)
        % 6) Generate weight maps
        % 7) Describe the outputs
        
        % X = trials x voxels
        % Y = trails x 1 [column vector of conditions]
        
        %% IN PROGRESS
        
        % Let us use the native space...
        EPIfn = aas_getfiles_bystream(aap,p,s,'epi');
        
        
        
        %% DESCRIBE OUTPUTS!
        
        % DMLT outputs
        aap=aas_desc_outputs(aap,p,'DMLT',outDMLT);
        
end