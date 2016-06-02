function [aap, resp] = aamod_meg_grandmean(aap,task)

resp='';

switch task
    case 'doit'
        fname = {};
        switch aap.internal.inputstreamsources{aap.tasklist.currenttask.modulenumber}.stream.sourcedomain
            case 'subject' % merged data
                for subj = 1:numel(aap.acq_details.subjects)
                    megdata = aas_getfiles_bystream(aap,'subject',subj,'meg');
                    fname{end+1} = spm_file(megdata(1,:),'ext','mat');
                end
            case 'meg_session' % non-merged data
                for subj = 1:numel(aap.acq_details.subjects)
                    [junk, megsess] = aas_getN_bydomain(aap,'meg_session',subj);
                    for sess = megsess
                        megdata = aas_getfiles_bystream(aap,'meg_session',[subj sess],'meg');
                        fname{end+1} = spm_file(megdata(1,:),'ext','mat');
                    end
                end
        end

        S.D = char(fname);
        S.weighted = aas_getsetting(aap,'weighted');
        S.outfile = fullfile(aas_getstudypath(aap), aas_getsetting(aap,'outfile'));
        spm_eeg_grandmean(S);
        
        %% Outputs
        aap=aas_desc_outputs(aap,'study',[],'meg_grandmean',...
            char(spm_file(S.outfile,'ext','mat'),spm_file(S.outfile,'ext','dat')));
end