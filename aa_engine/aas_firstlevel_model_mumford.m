function aas_firstlevel_model_mumford(aap, anadir, coreSPM, files, allfiles, ...
    model, modelC, eventNumber, sessNumber, numReg, nDest, ...
    movementRegs, compartmentRegs, physiologicalRegs, spikeRegs, GLMDNregs)


Tmodel = cell(size(model));


% Make temporary directory inside folder
Tanadir = (fullfile(anadir, sprintf('temp_%03d', numReg)));


try rmdir(Tanadir); catch; end
mkdir(Tanadir)


cols_nuisance=[];
cols_interest=[];
sessnuminspm=0;
currcol=1;
rSPM = coreSPM;


for sess = aap.acq_details.selected_sessions
    sessnuminspm=sessnuminspm+1;
    
    rSPM.nscan(sessnuminspm) = size(files{sess},1);
    rSPM.xX.K(sessnuminspm).HParam = aap.tasklist.currenttask.settings.highpassfilter;
    
    % * Get noise regressor numbers
    noiseRegs = eventNumber(sessNumber == sess);
    if ismember(numReg, noiseRegs)
        noiseRegs(numReg) = [];
    end
    
    % No need to check model{sess} again, did already before...
    Tmodel{sess} = model{sess};
    Tmodel{sess}.event = [];
    % Event names, rather simple...
    Tmodel{sess}.event(1).name = 'Noise';
    % I don't see a way to do parameteric stuff here
    Tmodel{sess}.event(1).parametric = [];
    % Easy to set up the Regressor event
    try
        Tmodel{sess}.event(2).ons = model{sess}.event(numReg).ons;
        Tmodel{sess}.event(2).dur = model{sess}.event(numReg).dur;
        Tmodel{sess}.event(2).name = 'Reg';
        Tmodel{sess}.event(2).parametric = [];
    catch
        % If we can't model that event... don't model it
    end
    % Trickier for the Noise event
    Tons = [];
    Tdur = [];
    for r = noiseRegs
        % Make sure we don't include the modelled regressor in
        % the noise regressor
        Tons = [Tons; model{sess}.event(r).ons(:)];
        Tdur = [Tdur; model{sess}.event(r).dur(:)];
    end
    
    % Sort the onsets, and apply same reordering to durations
    [Tons, ind]=sort(Tons);
    if (length(Tdur)>1)
        Tdur=Tdur(ind);
    end
    Tmodel{sess}.event(1).ons = Tons;
    Tmodel{sess}.event(1).dur = Tdur;
    
    rSPM.Sess(sessnuminspm).C.C = [];
    rSPM.Sess(sessnuminspm).C.name = {};
    
    [rSPM, cols_interest, cols_nuisance, currcol] = aas_firstlevel_model_define(aap, sess, sessnuminspm, rSPM, Tmodel, modelC, ...
        cols_interest, cols_nuisance, currcol, ...
        movementRegs, compartmentRegs, physiologicalRegs, spikeRegs, GLMDNregs);
    
end
cd(Tanadir)


% DEBUG
%{
if numReg == 1
    subplot(3,1,1); plot(SPMdes.xX.X(1:100,1:sessRegs(end))); title('Normal model')
    subplot(3,1,2); plot(rSPMdes.xX.X(1:100,1:2)); title('New model')
    subplot(3,1,3); plot([SPMdes.xX.X(1:100,1) sum(SPMdes.xX.X(1:100,2:sessRegs(end)),2)]); title('Old model, summed')
end
%}
rSPM.swd = Tanadir;
rSPM.xY.P = allfiles;
rSPMdes = spm_fmri_spm_ui(rSPM);


% now check real covariates and nuisance variables are
% specified correctly
rSPMdes.xX.iG=cols_nuisance;
rSPMdes.xX.iC=cols_interest;


% Turn off masking if requested
if ~aap.tasklist.currenttask.settings.firstlevelmasking
    rSPMdes.xM.I=0;
    rSPMdes.xM.TH=-inf(size(rSPMdes.xM.TH));
end


spm_unlink(fullfile('.', 'mask.img')); % avoid overwrite dialog
rSPMest = spm_spm(rSPMdes);


% After this, move the betas to correct location
% Find the regressors we wish to move...
% They contain the string Reg...
nOrig = find(~cellfun('isempty', strfind(rSPMest.xX.name, 'Reg')));


% Now move the actual files
for f = 1:length(nOrig)
    unix(['mv ' fullfile(Tanadir, sprintf('beta_%04d.img', nOrig(f))) ...
        ' ' fullfile(anadir, sprintf('beta_%04d.img', nDest(f)))]);
    unix(['mv ' fullfile(Tanadir, sprintf('beta_%04d.hdr', nOrig(f))) ...
        ' ' fullfile(anadir, sprintf('beta_%04d.hdr', nDest(f)))]);
end

cd(anadir);

% Delete the Tanadir
unix(['rm -rf ' Tanadir]);