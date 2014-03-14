% aa module - calculates spectrum in time of different tissue compartments
% Rhodri Cusack BMI Western Nov 2013

function [aap,resp]=aamod_tissue_wavelets_summarize(aap,task)
resp='';

splitbyage=false;

switch task
    case 'report'
        
    case 'doit'
        
        %% Load up and summarize wavelets
        allages=[];
        C={};
        nsubj=length(aap.acq_details.subjects);
        for sessind=aap.acq_details.selected_sessions
            for subjind=1:nsubj
                tw= load(aas_getfiles_bystream(aap,'session',[subjind,sessind],'tiss_wavelets'));
                tw=tw.tiss_wavelet;
                for tissind=1:length(tw)
                    C{subjind}{tissind}={tw{tissind}.ca tw{tissind}.cd{:}};
                    if subjind==1 && sessind==1
                        lev=tw{1}.level;
                        minC=inf(lev+1,1);
                        maxC=-inf(lev+1,1);
                        ntiss=length(tw);
                    end;
                    for levind=1:lev+1
                        minC(levind)=min(minC(levind),prctile(C{subjind}{tissind}{levind}(:),0.1));
                        maxC(levind)=max(maxC(levind),prctile(C{subjind}{tissind}{levind}(:),99.9));
                    end;
                end;
            end;
            nbins=1000;
            hists=zeros(nsubj,ntiss,lev+1,nbins);
            binX=[];
            for levind=1:lev+1
                
                binX(levind,:)=linspace(minC(levind),maxC(levind),nbins);
                for subjind=1:nsubj
                    for tissind=1:ntiss
                        hists(subjind,tissind,levind,:)=hist(C{subjind}{tissind}{levind}(:),binX(levind,:));
                    end;
                end;
            end;
            figure(9+sessind);
            ind=1;
            for tissind=1:ntiss
                for levind=1:lev+1
                    subplot(ntiss,lev+1,ind);
                    plot(binX(levind,:),squeeze(hists(:,tissind,levind,:)));
                    ind=ind+1;
                    
                    
                end;
            end;
            
         % Describe outputs
        ps_fn=fullfile(aas_getstudypath(aap),'tiss_wavelet_summarize.mat');
        save(ps_fn,'hists','binX');
        aap=aas_desc_outputs(aap,'study',[],'tiss_wavelet_summarize',ps_fn);
            
        end;
       
        
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;



