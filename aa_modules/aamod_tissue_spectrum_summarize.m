% aa module - calculates spectrum in time of different tissue compartments
% Rhodri Cusack BMI Western Nov 2013

function [aap,resp]=aamod_tissue_spectrum_summarize(aap,task)
resp='';

splitbyage=false;

switch task
    case 'report'
        
    case 'doit'
        
        %% Load up and summarize
        powerspect=[];
        allages=[];
        for sessind=aap.acq_details.selected_sessions
            for subjind=1:length(aap.acq_details.subjects)
                thispowerspect= load(aas_getfiles_bystream(aap,'session',[subjind,sessind],'tiss_powerspect'));
                deltaf= thispowerspect.freq(2)- thispowerspect.freq(1);
                powerspect(subjind,sessind,:,:)=thispowerspect.tiss_powerspect/deltaf; %/sqrt(thispowerspect.tiss_powerspect(1));
                if splitbyage
                    dcm=load(aas_getfiles_bystream(aap,'session',[subjind, sessind],'epi_header'));
                    age=dcm.DICOMHEADERS{1}.PatientsAge;
                    allages(subjind)=str2num(age(1:end-1));
                end;
            end;
            
            freq=thispowerspect.freq;
            period=1./freq;
            
            polyfitS=[];
            polyfitP=[];
            polyfitS_ind=[];
            polyfitP_ind=[];
            % Plot mean
            for tissind=1:size(powerspect,3)
                figure(10+sessind);
                subplot(2,size(powerspect,3),tissind);
                loglog( freq(2:end),squeeze(powerspect(:,sessind,tissind,2:end))');
                title('Split by subject');
                xlabel('Frequency (Hz)')
                subplot(2,size(powerspect,3),tissind+size(powerspect,3));
                loglog( freq(2:end),squeeze(mean(powerspect(:,sessind,tissind,2:end),1))');
                xlabel('Frequency (Hz)')
                title('Mean across subjects');
                
                % Period on x-axis
                figure(20+sessind);
                % Split by subj
                subplot(2,size(powerspect,3),tissind);
                plot( period(2:end),squeeze(powerspect(:,sessind,tissind,2:end))');
                title('Split by subject');
                xlabel('Period (s)')
                for subjind=1:size(powerspect,1)
                    [P S]=polyfit(period(2:end),squeeze(powerspect(subjind,sessind,tissind,2:end))',1);
                    polyfitP_ind(subjind,tissind,:)=P;
                    if subjind==1 && tissind==1
                        polyfitS_ind=S;
                    else
                        polyfitS_ind(subjind,tissind)=S;
                    end;
                end;
                
                % Mean across subj
                subplot(2,size(powerspect,3),tissind+size(powerspect,3));
                plot( period(2:end),squeeze(mean(powerspect(:,sessind,tissind,2:end),1))');
                xlabel('Period (s)')
                title('Mean across subjects');
                
                [P S]=polyfit(period(2:end),squeeze(mean(powerspect(:,sessind,tissind,2:end),1))',1);
                polyfitP(tissind,:)=P;
                if tissind==1
                    polyfitS=S;
                else
                    polyfitS(tissind)=S;
                end;
            end;
            if splitbyage
                % Now by decile
                meanbydecile=[];
                for decind=1:7
                    decile=8+10*decind;
                    for tissind=1:size(powerspect,3)
                        subjmask=(allages>=decile) & (allages<(decile+10));
                        meanbydecile(tissind,decind,:)=mean(powerspect(subjmask,sessind,tissind,:),1);
                    end;
                end;
                
                
                %% Graph deciles
                figure(20);
                for tissind=1:size(meanbydecile,1)
                    subplot(size(meanbydecile,1),1,tissind);
                    plot([1:7],squeeze(meanbydecile(tissind,:,2:end)));
                end;
            end;
        end;
        % Describe outputs
        ps_fn=fullfile(aas_getstudypath(aap),'tiss_powerspect_summarize.mat');
        save(ps_fn,'powerspect','polyfitP','polyfitS','polyfitP_ind','polyfitS_ind','period','freq');
        aap=aas_desc_outputs(aap,'study',[],'tiss_powerspect_summarize',ps_fn);
        
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;



