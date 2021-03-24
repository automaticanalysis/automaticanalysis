function [aap, resp] = aamod_meg_convert(aap,task,subj,sess)

resp='';

switch task
    case 'report'
        
    case 'doit'
        %% Initialise
        [junk, SPMtool] = aas_cache_get(aap,'spm');
        SPMtool.addCollection('meeg');

        sessdir = aas_getsesspath(aap,subj,sess);
        inputstreams = aas_getstreams(aap,'input');
        infname = aas_getfiles_bystream(aap,'meeg_session',[subj sess],inputstreams{2});        
        outfname = fullfile(sessdir,['spm12_' basename(infname)]); % specifying output filestem
        
        chan=load(aas_getfiles_bystream(aap,'channellabels'));
        
        S=[];
        S.dataset       = infname;
        S.outfile       = outfname;
        S.save          = 0;
        S.reviewtrials  = 0;
        S.channels      = chan.label;  
        S.continuous    = 1;
        S.checkboundary = 0;
        
        D = spm_eeg_convert(S);
        
        % Relevant to triggers, and no triggers at rest!
        tc = find(strcmp(D.chanlabels,'STI101'));
        if ~isempty(tc)
            D(tc,:) = D(tc,:)/1e6;  % 1e6 new scaling introduced by SPM12!
            D = units(D,tc,'V');
            D.save;
        end
        
        %% Rename any channels
        rename_chans = aap.tasklist.currenttask.settings.ChannelRename;
        if size(rename_chans,1) > 0
            for e = 1:size(rename_chans,1)
                ci = indchannel(D,rename_chans{e,1});
                if isempty(ci)
                    aas_log(aap,false,sprintf('WARNING: No channel %s found',rename_chans{e,1}));
                    continue;
                end
                D = chanlabels(D, ci, rename_chans{e,2});
                if ~isempty(rename_chans{e,3})
                    D = chantype(D, ci, rename_chans{e,3});
                end
            end
            D.save
        end
        
        SPMtool.rmCollection('meeg');

        %% Outputs
        aap=aas_desc_outputs(aap,subj,sess,'meeg',char([outfname '.dat'],[outfname '.mat']));        
    case 'checkrequirements'
        if ~aas_cache_get(aap,'spm'), aas_log(aap,true,'SPM is not found'); end
end
