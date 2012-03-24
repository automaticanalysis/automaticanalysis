function [aap resp]=aamod_emeg_rfx(varargin)
% Danny Mitchell 13/03/08;
% Parallelised 11/08
% Adapted for source volumes 04/09

%% check task settings, subject, block etc
[aap subblock doit resp settings]=aa_emeg_checktasksettings(mfilename('fullpath'),varargin);
if ~doit; return; end

%% Parallelise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gfile=fullfile(aap.acq_details.root,...
    sprintf('%s_%g.parallel.mat',mfilename,(aap.tasklist.currenttask.index)));
if isempty(varargin) || strcmp(varargin{2},'parallelise')
    % prepare globals and looping variable
    %% get files to average
    fprintf('\nFinding files from %g blocks of %g subjects', ...
        length(aap.acq_details.sessions),length(aap.acq_details.subjects))

    fnames={};warnings=[];
    for s=1:length(aap.acq_details.subjects)
        fprintf('.')
        for b=aap.acq_details.selected_sessions
            bdir=aas_getsesspath(aap,s,b);
            ds=dir(bdir);
            di=find(~cellfun('isempty',regexp({ds.name},settings.DirFilter)) & [ds.isdir]);
            for d=1:length(di)
                % check trial numbers agree with level specification
                ts=dir(fullfile(bdir,ds(di(d)).name));
                ti=find(~cellfun('isempty',regexp({ts.name},'trialtype')) & [ts.isdir]);
                if length(ti)==prod(settings.levels(2:end));
                    warnings(end+1)=0;
                else warnings(end+1)=length(ti);
                end
                % find files
                tempfile=spm_select('List',fullfile(bdir,ds(di(d)).name,'trialtype1'),settings.FileFilter);
                tempdir=regexprep(ds(di(d)).name,aap.acq_details.sessions(b).name,'BLOCK');
                fnames=[fnames; strcat(tempdir,'/TRIAL/',cellstr(tempfile))];
            end
        end % next block
    end % next subject

    if all(warnings)
        % could check whether block has been forgotten and automatically fix?
        % (length(aap.acq_details.sessions)~=settings.levels(1))
        error('Specified levels are not always consistent with the number of trialtype directories.')
    elseif any(warnings)
        error('Specified levels are not always consistent with the number of trialtype directories.')
    elseif length(settings.factors)~=length(settings.levels)
        error('Number of factor labels differs from length of levels vector')
    end

    fnames=unique(fnames);
    voldirs=regexprep(fnames,'/.*','');
    files=regexprep(fnames,'.*/','');

    loopvar={};
    for v=1:length(fnames)
        for e=1:size(settings.effects,1)
            loopvar{end+1}=struct('effects',{settings.effects(e,:)});
            
            % reformat any user-specified contrasts to string for id
            temp=settings.effects;
            usercon=find(cellfun(@isnumeric,settings.effects(e,:)));
            for uc=usercon
                temp{e,uc}=regexprep(mat2str(temp{e,uc}),{'\.',' '},{'p','.'});
            end           
            loopvar{end}.id=[regexprep(fnames{v},{'/','.nii'},{'.',''}) '_' cell2mat(temp(e,:))];
            
            loopvar{end}.voldir=voldirs{v};
            loopvar{end}.fname=files{v};
            loopvar{end}.factors=settings.factors;
            loopvar{end}.levels=settings.levels;
            if isfield(settings,'ROIs') && ~isempty(settings.ROIs);
                loopvar{end}.ROIs=settings.ROIs;
            end
        end
    end

    %% confirm/create group analysis directory
    try cd(fullfile(aap.acq_details.root,'GroupAnalysis'));
    catch
        mkdir(fullfile(aap.acq_details.root,'GroupAnalysis'));
        cd(fullfile(aap.acq_details.root,'GroupAnalysis'));
    end
    if ~exist('figures','dir'); mkdir('figures'); end

    %% save variables for parallelisation
    save(gfile,'loopvar');

    if ~isempty(varargin); return; end
end

%% do it %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try load(gfile); loopvar;
catch aas_log(aap,1,sprintf('\n***Failed to prepare parallelisation structure!\n%s',resp));
end;

if isempty(subblock); subblock=1:length(loopvar);
else subblock=subblock{1};
end

for e=subblock
    
    % run repeated measures ANOVA
    S=aas_emeg_ANOVA_1stLevel(aap,loopvar{e},settings.Overwrite);
    if settings.PPM; S.PPM=true; end
    outputdirs=aas_emeg_ANOVA_2ndLevel(aap,S,settings.Overwrite);
    
    if settings.rendersources && strcmpi(S.datatype,'source') % Render source level SPM/PPM
        for d=1:length(outputdirs)        
            % See aas_emeg_rmANOVApartitioned for
            % old implementation of sensor level significance checking
            R.condir=outputdirs{d};
            R.Ic='P';
            R.brain='/imaging/local/spm/spm5/rend/render_no_cereb.mat';
            R.style=1;
            R.overwrite=true;
            R.title=S.fname;
            aas_rendersources(R);
        end % next effect        
    end   
    
end % next set of effects

% aas_emeg_plot_sensor_spms(aap,settings.Overwrite);

return