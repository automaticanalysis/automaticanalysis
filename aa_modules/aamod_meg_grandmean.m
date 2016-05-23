function [aap, resp] = aamod_meg_grandmean(aap,task)

resp='';

switch task
    case 'doit'
        fname = {};
        switch aap.internal.inputstreamsources{aap.tasklist.currenttask.modulenumber}.stream.sourcedomain
            case 'subject' % merged data
                for subj = 1:numel(aap.acq_details.subjects)
                    fname{end+1} = aas_getfiles_bystream(aap,'subject',subj,'meg');
                end
            case 'meg_session' % non-merged data
                for subj = 1:numel(aap.acq_details.subjects)
                    [junk, megsess] = aas_getN_bydomain(aap,'meg_session',subj);
                    for sess = megsess
                        fname{end+1} = aas_getfiles_bystream(aap,'session',[subj sess],'meg');
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