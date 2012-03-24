function [aap resp]=aamod_emeg_importfif(varargin)
% Imports EMEG data from fif file to spm/matlab format
% Calls spm_eeg_rdata_FIF_dm.m
% If data is continuous, default is to load it all (then epoch with spm_epoch);
% If data is averaged, each set will be saved to a different file
%
% Danny Mitchell 03/03/08

%% check task settings, subject, block etc
[aap subblock doit resp settings]=aa_emeg_checktasksettings(mfilename('fullpath'),varargin);
if ~doit; return; end

%% find files and decide whether to run task;
files=aas_emeg_findfiles(aap,settings.InputFilter,subblock);
if isempty(files); aas_log(aap,1,sprintf('\nFound no data! (Input filter is %s)\n',settings.InputFilter)); end

%% run task for each file
for f=1:length(files);

%% prepare parameters
    S=settings;
    S.Fdata=files{f};
    try S.HPIfile=regexprep(files{f},'raw.*','raw_HPI.mat');catch end;
        
%% import averaged data from fif file (OLD!! - NOT SURE IF IT STILL WORKS)
    if ~isempty(findstr(S.InputFilter,'_av'))
       % averaged data
       load(strrep(files{f},'.fif','_nave')); % load time windows & number of averages from mat file
       N=loadfif(files{f},'sets');
       % import all as one file with time window as large as necessary 
       S.twin(1)=min(windows(:,1));
       S.twin(2)=max(windows(:,2));
       S.conds=(1:N);
       D=spm_eeg_rdata_FIF(S); % import data
       D.events.repl = nave; % add number of averages from mat file
       save(fullfile(D.path, D.fname), 'D');
    else
%% import continuous data from fif file
        [pth name]=fileparts(files{f});
        outfile=fullfile(aas_getsesspath(aap,subblock{1},subblock{2}),[name '.mat']);
        
        [na,ki]=channames(files{f});
        if sum(ki==2)>2;
            % also check eeg has been imported
            checkoutfiles=char(outfile,strrep(outfile,'.mat','-eeg.mat'));
        else
            checkoutfiles=outfile;
        end
        if ~all(cellfun(@exist,cellstr(checkoutfiles))) || S.Overwrite==1;
            fprintf('\nImporting continuous data...')
            [temp, D.Radc] = rawchannels(files{f},S.trig_chan);
            S.Pout=outfile;
            S.twin = ([1 length(temp)]-1)/D.Radc*1000; % import it all (in ms!)
            [D B H]=spm_eeg_rdata_FIF_dm(S);
        else
            try
                load(outfile);
            catch
                fprintf('\n !! Please debug %s !!\n',mfilename); keyboard
            end
        end
    end % continuous or averaged
    
%% Add # harmonic components and sphere details from MaxFilter log if it exists
    try 
        D.MaxFilter='';
        %logfile=regexprep(files{f},'\.fif$','_MaxFilterLog.txt');
        logfiles=spm_select('FPlist',D.path,[regexprep(D.fname,'\.mat$','') '_MaxFilterLog\d?.txt']);
        logfile=deblank(logfiles(end,:));
        fid = fopen(logfile,'rt');
        tline='';warnings={};collectwarnings=false;
        while 1==1
            tline = fgetl(fid);
            if ~ischar(tline); break; end;
            if ~isempty(regexp(tline,'Using .* harmonic components','once'))
                nin=regexp(tline,'Using (\d*)/\d* inside','tokens');
                D.MaxFilter.HarmonicComponents=str2double(nin{1}{1});
                D.MaxFilter.HarmonicComponentsTxt=tline;
            elseif ~isempty(regexp(tline,'Inside origin.* head frame','once'))
                ori=regexp(tline,'origin:..(.*) mm.*head frame','tokens');
                D.MaxFilter.InsideOriginInHeadFrame=str2num(ori{1}{1}); % vector
            elseif ~isempty(regexp(tline,'min,max,ave distance to origin','once'))
                radii=regexp(tline,'distance to origin (\d*), (\d*)','tokens');
                D.MaxFilter.MinDistanceToOrigin=str2double(radii{1}{1});
                D.MaxFilter.MaxDistanceToOrigin=str2double(radii{1}{2});
            elseif ~isempty(regexp(tline,'auto bad channels','once'))
                abc=regexp(tline,'channels \((.*)\)','tokens');
                D.MaxFilter.AutoBadChans=str2num(abc{1}{1}); % vector    
            elseif ~isempty(regexp(tline,'user bad channels','once'))
                ubc=regexp(tline,'channels \((.*)\)','tokens');
                D.MaxFilter.UserBadChans=str2num(ubc{1}{1}); % vector 
            elseif ~isempty(regexp(tline,'Retrying with force','once'))
                collectwarnings=false;
            elseif ~isempty(regexp(tline,'Warning:','once')) || collectwarnings
                collectwarnings=true;
                warnings=[warnings;tline];
            end
            
            if all(isfield(D.MaxFilter,{'HarmonicComponents','InsideOriginInHeadFrame','MinDistanceToOrigin','MaxDistanceToOrigin'}));
                save(fullfile(D.path,D.fname),'D');
                break            
            end
        end
        fclose(fid);
    catch
        aas_log(aap,0,sprintf('\nFailed to extract all details from MaxFilter log: %s\n',logfile))
    end;
    
%% add spheres to diagnostic and print
    F=spm_figure('findwin');
    if ~isempty(F) && exist('H','var');
        set(F,'renderer','opengl')
        try
            ori=D.MaxFilter.InsideOriginInHeadFrame;
            oritext=sprintf('SSS expansion origin = %g, %g, %g, in head space',ori(1),ori(2),ori(3));
        catch oritext='SSS expansion origin unknown/unused';
        end
        try
            rmax=D.MaxFilter.MaxDistanceToOrigin;
            rmin=D.MaxFilter.MinDistanceToOrigin;
            radiustext=sprintf('Inner radius = %g mm; outer radius = %g mm',rmin, rmax);
        catch radiustext ='';       
        end
        try
            HarmonicComponentsTxt=D.MaxFilter.HarmonicComponentsTxt;
        catch HarmonicComponentsTxt='';      
        end
        try
            badchantext=sprintf('User-marked bad channels: %s \nAuto-marked bad channels: %s', ...
                ubc{1}{1}, abc{1}{1});
        catch badchantext='';
        end
        
        if exist('ori','var') && exist('rmin','var')
            [x y z]=sphere(30);
            for a=1:3
                axes(H(a));
                plot3(ori(1),ori(2),ori(3),'k*','markersize',10,'linewidth',2)
                inner=surf(x*rmin+ori(1),y*rmin+ori(2),z*rmin+ori(3), ...
                    'facecolor','k','facealpha',0.1,'edgecolor','none', ...
                    'facelighting','phong','specularexponent',0.4,'specularstrength',2);
                outer=surf(x*rmax+ori(1),y*rmax+ori(2),z*rmax+ori(3), ...
                    'facecolor','k','facealpha',0.1,'edgecolor','none', ...
                    'facelighting','phong','specularexponent',0.3,'specularstrength',2);
                L=camlight('headlight');
                axis tight
            end
        end
        anntop=annotation(F,'textbox',[0 0.88 1 .12],'String', ...
            {'Geometry of sensors, scalp & MaxFilter origin,', ...
            sprintf('imported into %s',fullfile(D.path,D.fname)), ...
            oritext, radiustext, HarmonicComponentsTxt, badchantext...
            sprintf('Date:...%s',datestr(now,0))}, ...
            'edge','none','fontsize',10); 
        
        if ~isempty(warnings)
            annbot=annotation(F,'textbox',[0 0 1 .15],'String',warnings);
        end
        
        print('-dpng', '-r82',fullfile(D.path,'figures',sprintf('%s_SSS_sensors_HPI.png',name))); % better for web browser
        print('-dpng', '-r150',fullfile(D.path,'figures',sprintf('%s_SSS_sensors_HPI.hires',name)));
        delete(F);
    end
    
    fprintf('.');
end

return
