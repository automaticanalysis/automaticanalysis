% AA module - smoothing
% [aap,resp]=aamod_smooth(aap,task,subj,sess)
% Gaussian smoothing images, important for Gaussian Field Theory
% Kernel size determined by aap.spm_analysis.FWHM
% Rhodri Cusack MRC CBU Cambridge Jan 2006- Aug 2007
% Now once per session for improved performance when running in parallel

function [aap,resp]= aamod_ggmfit(aap,task,subj,sess)
resp='';

switch task
    case 'doit'
        subjname = aas_prepare_diagnostic(aap, subj);
        
        % Is session specified in task header?
        if (isfield(aap.tasklist.currenttask.settings,'session'))
            sess = aap.tasklist.currenttask.settings.session;
        end
        
        streams=aap.tasklist.currenttask.inputstreams.stream;
        
        for streamind=1:length(streams)
            
            % Images to ggmfit
            if (exist('sess','var'))
                P = aas_getfiles_bystream(aap,subj,sess,streams{streamind});
            else
                P = aas_getfiles_bystream(aap,subj,streams{streamind});
            end
            
            % now ggmfit
            for p = 1:size(P,1)
                V = spm_vol(P(p,:));
                Y = spm_read_vols(V);
                
                % Only non-zero non-nan values of image...
                M = and(isfinite(Y), Y~=0);
                
                % First we try the ggm model (Gaussian/Gamma)
                ggmmix = ggmfit(Y(M)', 3, 'ggm');
                
                % If this does not work, we try with 2 mixtures
                if ~isfinite(ggmmix.mus(1)) || ggmmix.mus(1) == 0 ...
                        || ~isfinite(ggmmix.sig(1)) || ggmmix.sig(1) == 0
                    aas_log(aap,0,'Error in ggm, mu and/or sigma are NaN, trying ggm with 2 mixtures...')
                    ggmmix = ggmfit(Y(M)', 2, 'ggm');
                    
                    if isnan(ggmmix.mus(1)) || isnan(ggmmix.sig(1))
                        aas_log(aap,1,'Error in ggm, mu and/or sigma are NaN')
                    end
                end
                
                % Standardise image: subtract mean, divide by stdev
                Y(M) = (Y(M) - ggmmix.mus(1)) ./ ggmmix.sig(1);
                Y(~M) = NaN;
                
                % Put back into image...
                spm_write_vol(V,Y);
            end
            
            % DEBUG (uncomment if you have trouble)            
            h = img2hist(P, [], 'Contrast distributions');
            legend('off')
            print('-depsc2', fullfile(aap.acq_details.root, 'diagnostics', ...
                [mfilename '_' subjname '.eps']));
            try close(h); catch; end
            
            % Describe outputs
            if (exist('sess','var'))
                aap=aas_desc_outputs(aap,subj,sess,streams{streamind},P);
            else
                aap=aas_desc_outputs(aap,subj,streams{streamind},P);
            end
            
        end;
        % All done
        spm_progress_bar('Clear');
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;



