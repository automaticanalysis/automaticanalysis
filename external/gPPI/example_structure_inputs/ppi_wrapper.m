function ppi_wrapper(regionnumber,firstsubject,lastsubject)

%User input required (directories and subjects)
addpath('PPPIdirectory')
addpath('spm8directory')
Subjects={'subject1' 'subject2'};

%User input required (region files)
regionfile={'region1.nii'... 
'region2.nii'};

%User input required (region names)
region={'region1'... 
'region2'};

%User input required (master template)
load('ppi_master_template.mat')

P.VOI=regionfile{regionnumber};
P.Region=region{regionnumber};

%User input required (change directory to where the input structure should
%be saved)
save(['directory' region{regionnumber} '.mat'],'P');

for i=firstsubject:lastsubject;
    try
        %User input required directory of SPM.mat files
        Directory=['subjectdirectory'];
        cd(Directory)
        
        %User input required (directory same as line 23 above)
        load(['directory' region{regionnumber} '.mat']);
        
        P.subject=Subjects{i};
        P.directory=Directory;
        
        %User input required (change analysis to be more specific)
        save([Subjects{i} '_analysis_' region{regionnumber} '.mat'],'P');
        PPPI([Subjects{i} '_analysis_' region{regionnumber} '.mat']);
    catch
        disp(['Failed: ' Subjects{i}])
    end
end
end
