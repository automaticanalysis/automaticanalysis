function [P, error]=PPPI_checkstruct(inputstructure,structfile)
% This function check that all inputs are correctly specified AND that the
% values are valid.
%
%   See PPPI for required and optional fields for inputstructure
%
% License:
%   Copyright (c) 2011, Donald G. McLaren and Aaron Schultz
%   All rights reserved.
%
%    Redistribution, with or without modification, is permitted provided that the following conditions are met:
%    1. Redistributions must reproduce the above copyright
%        notice, this list of conditions and the following disclaimer in the
%        documentation and/or other materials provided with the distribution.
%    2. All advertising materials mentioning features or use of this software must display the following acknowledgement:
%        This product includes software developed by the Harvard Aging Brain Project.
%    3. Neither the Harvard Aging Brain Project nor the
%        names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
%    4. You are not permitted under this Licence to use these files
%        commercially. Use for which any financial return is received shall be defined as commercial use, and includes (1) integration of all 	
%        or part of the source code or the Software into a product for sale or license by or on behalf of Licensee to third parties or (2) use 	
%        of the Software or any derivative of it for research with the final aim of developing software products for sale or license to a third 	
%        party or (3) use of the Software or any derivative of it for research with the final aim of developing non-software products for sale 
%        or license to a third party, or (4) use of the Software to provide any service to an external organisation for which payment is received.
%
%   THIS SOFTWARE IS PROVIDED BY DONALD G. MCLAREN (mclaren@nmr.mgh.harvard.edu) AND AARON SCHULTZ (aschultz@nmr.mgh.harvard.edu)
%   ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
%   FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
%   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
%   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR 
%   TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
%   USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
%   DAMAGE.
%   
%
%   Created on 9/27/2011 by Donald G. McLaren
%   (mclaren@nmr.mgh.harvard.edu)
%   GRECC, Bedford VAMC
%   Department of Neurology, Massachusetts General Hospital and Havard
%       Medical School

%% Begin program
if nargin<2 && ~isstruct(inputstructure) && ~exist(inputstructure,'file')
    disp('Wrong inputs. File name must be specified as the 1st argument.')
    return;
else
    if isstruct(inputstructure)
    elseif exist(inputstructure,'file')
        try 
            structfile;
            inputstructure=load(inputstructure);
        catch
            [path file ext]=fileparts(inputstructure);
            inputstructure=load(inputstructure);
            if isempty(path)
                path=pwd;
            end
            try 
                structfile=[path filesep inputstructure.subject '_' file ext];
            catch
                structfile=[path filesep 'invalidPPI_' file ext];
            end
        end
    elseif nargin==2
        try
            inputstructure=load(structfile);
            [path file ext]=fileparts(structfile);
            if isempty(path)
                path=pwd;
            end
            try 
                structfile=[path filesep inputstructure.subject '_' file ext];
            catch
                structfile=[path filesep 'invalidPPI_' file ext];
            end
        catch
            disp('Wrong inputs. Structfile could not be loaded.')
            return;
        end
    else
        disp('Wrong inputs. File name or parameter structure must be specified as the 1st argument.')
        return;
    end
end

while numel(fields(inputstructure))==1
    F=fieldnames(inputstructure);
    inputstructure=inputstructure.(F{1}); %Ignore coding error flag.
end

error={};
if isstruct(inputstructure)
    %Check subject
    try
        P.subject=inputstructure.subject;
    catch
        error{end+1}='Program will exit. Subject is not specified.';
    end
    if exist('P','var') && isfield(P,'subject') && ~ischar(P.subject);
            error{end+1}='Program will exit. Subject must be a string.';
    end
    
    %Check SPM directory
    try
        P.directory=inputstructure.directory;
        if isempty(P.directory)
           invokecatchstatement
        end
        if exist([P.directory filesep 'SPM.mat'],'file')~=2
            error{end+1}='Program will exit. Directory is incorrect or does not contain an SPM.mat file.';
            SPM.swd=pwd;
            SPMexist=0;
        else
            SPMexist=1;
            load([P.directory filesep 'SPM.mat']);
        end
    catch
        SPMexist=0;
        error{end+1}='Program will exit. Directory is not specified.';
    end

    %Check VOI
    try
        P.VOI=inputstructure.VOI;
        if isstruct(P.VOI)
            if ~isfield(P.VOI,'VOI')
                error{end+1}='Program will exit. VOI not specified.';
                invokecatchstatement
            elseif exist(P.VOI.VOI,'file')~=2
                error{end+1}='Program will exit. VOI file does not exist.';
            end
            if ~any([~isempty(strfind(P.VOI.VOI,'.img')) ~isempty(strfind(P.VOI.VOI,'.nii'))])
                error{end+1}='Program will exit. VOI file is not a img or nii file.';
            end
            if ~isfield(P.VOI,'masks')
                P.VOI.masks=[];   
            end
            for jj=1:numel(P.VOI.masks)
                if isempty(fileparts(P.VOI.masks{jj}))
                   P.VOI.masks{jj}=[SPM.swd filesep P.VOI.masks{jj}]; 
                end
                if exist(P.VOI.masks{jj},'file')~=2
                    error{end+1}=['Program will exit. Mask statistic file ' P.VOI.masks{jj} ' does not exist (in VOI.masks field).'];
                end
                if ~any([~isempty(strfind(P.VOI.masks{jj},'.img')) ~isempty(strfind(P.VOI.masks{jj},'.nii'))])
                    error{end+1}='Program will exit. Mask statistic file is not a .img or .nii file.';
                end
            end
            if ~isfield(P.VOI,'exact') || P.VOI.exact~=1
                P.VOI.exact=0;
            else
                if isempty(which('peak_nii'))
                    error{end+1}='Program will exit. peak_nii.m not found. Please add it to the path.';
                end
            end
            if ~isfield(P.VOI,'thresh')
                P.VOI.thresh=[];
            elseif numel(P.VOI.thresh)~=numel(P.VOI.masks)
                error{end+1}='Program will exit. Thresholds do not match number of masks.';
            elseif iscell(P.VOI.thresh)
               error{end+1}='Program will exit. Thresholds cannot be in a cell array.';
            elseif ~all(P.VOI.thresh>=0)
                error{end+1}='Program will exit. Thresholds are not all greater than 0.';
            end
        else
            if exist(P.VOI,'file')~=2
                error{end+1}='Program will exit. VOI file does not exist.';
            end
            if ~any([~isempty(strfind(P.VOI,'.img')) ~isempty(strfind(P.VOI,'.nii')) ~isempty(strfind(P.VOI,'.mat'))])
                error{end+1}='Program will exit. VOI file does not a .img, .nii, or .mat file.';
            end
        end
    catch
        error{end+1}='Program will exit. VOI is not specified.';
    end
  
    
    %Check Region name
    try
        P.Region=inputstructure.Region;
        if ~ischar(P.Region)
            try 
                [path file]=fileparts(P.VOI);
                P.Region=file;
                clear path file
            catch
                error{end+1}='Program will exit. Region name must be a string.'; 
            end
        end
    catch
        try 
           [path file]=fileparts(P.VOI);
           P.Region=file;
           clear path file
        catch
            error{end+1}='Program will exit. Region name is not specified.';
        end
    end
   
    %Check Estimate, 0 means do not estimate PPI model
    try
        P.Estimate=inputstructure.Estimate;
        if isnumeric(P.Estimate)
            if P.Estimate~=0 && P.Estimate~=1 && P.Estimate~=2
                P.Estimate=0;
            end
        else
            error{end+1}='Program will exit. Estimate was not numeric. You must specify a number.';
        end
    catch
        P.Estimate=0;
    end

    
    if P.Estimate~=2
        %Check if contrast adjustment is present. Defaults to omnibus.
        if SPMexist==1
            try
                P.contrast=inputstructure.contrast;
                if (~iscellstr(P.contrast) && ~isnumeric(P.contrast)) || isempty(P.contrast)
                    P.contrast={'Omnibus F-test for PPI Analyses'};
                end
            catch
                P.contrast={'Omnibus F-test for PPI Analyses'};
            end
        end
        
        %Check extraction method, defaults to eigenvariate
        try
            P.extract=inputstructure.extract;
            if ischar(P.extract)
                if isempty(strfind(P.extract,'eig')) && isempty(strfind(P.extract,'mean'));
                    P.extract='eig';
                end
            else
                error{end+1}='Program will exit. Extract is not a string. You must specify P.extract as a string.';
            end
        catch
            P.extract='eig';
        end

        %Check Tasks
        try
            P.Tasks=inputstructure.Tasks;
            if ~iscellstr(P.Tasks)
                error{end+1}='Program will exit. Tasks are not defined in an array.';
            end
        catch
            P.Tasks={};
        end

        %Check Weights
        try
            P.Weights=inputstructure.Weights;
            if ~isnumeric(P.Weights)
                P.Weights=[];
            end
        catch
            P.Weights=[];
        end
        
        %Check for location to store seed region mask, only used with .mat
        %files
        try
            P.maskdir=inputstructure.maskdir;
        catch
            P.maskdir=P.directory;
        end
        

        %Check Equal ROIs
        try
            P.equalroi=inputstructure.equalroi;
            if P.equalroi~=0 && P.equalroi~=1
                P.equalroi=1;
            end
        catch
            P.equalroi=1;
        end

        if P.equalroi==0
            %Check for restricting ROI with first level mask -- use of
            %FLmask must be explicit
            try
                P.FLmask=inputstructure.FLmask;
                if P.FLmask~=0 && P.FLmask~=1
                    P.FLmask=0;
                end
            catch
                P.FLmask=0;
            end
        else
            try 
                P.FLmask=inputstructure.FLmask; 
                if P.FLmask~=0 && P.FLmask~=1
                    P.FLmask=0;
                end
            catch
                P.FLmask=0;
            end
        end
    end
    
    %Check for VOI2
    try
       P.VOI2=inputstructure.VOI2;
       if isstruct(P.VOI2)
           if ~isfield(P.VOI2,'VOI')
               error{end+1}='Program will exit. VOI2 file not specified.';
               invokecatchstatement
           elseif exist(P.VOI2.VOI,'file')~=2
               error{end+1}='Program will exit. VOI2 file does not exist.';
           end
           if ~any([~isempty(strfind(P.VOI2.VOI,'.img')) ~isempty(strfind(P.VOI2.VOI,'.nii'))])
               error{end+1}='Program will exit. VOI2 file is not a img or nii file.';
           end
           if ~isfield(P.VOI2,'masks')
               P.VOI2.masks=[];
           end
           for jj=1:numel(P.VOI2.masks)
               if isempty(fileparts(P.VOI2.masks{jj}))
                   P.VOI2.masks{jj}=[SPM.swd filesep P.VOI2.masks{jj}];
               end
               if exist(P.VOI2.masks{jj},'file')~=2
                   error{end+1}='Program will exit. Mask statistic file does not exist.';
               end
               if ~any([~isempty(strfind(P.VOI2.masks{jj},'.img')) ~isempty(strfind(P.VOI2.masks{jj},'.nii'))])
                   error{end+1}='Program will exit. Mask statistic file is not a .img or .nii file.';
               end
           end
           if ~isfield(P.VOI2,'exact') || P.VOI2.exact~=1
               P.VOI.exact=0;
           else
               if isempty(which('peak_nii'))
                   error{end+1}='Program will exit. peak_nii.m not found. Please add it to the path.';
               end
           end
           if ~isfield(P.VOI2,'thresh')
               P.VOI.thresh=[];
           elseif numel(P.VOI2.thresh)~=numel(P.VOI2.masks)
               error{end+1}='Program will exit. Thresholds do not match number of masks.';
           elseif iscell(P.VOI2.thresh)
               error{end+1}='Program will exit. Thresholds cannot be in a cell array.';
           elseif ~all(P.VOI2.thresh>=0)
                error{end+1}='Program will exit. Thresholds are not all greater than 0.';
           end
       else
            if exist(P.VOI2,'file')~=2
                error{end+1}='Program will exit. VOI2 file does not exist.';
            end
            if ~any([~isempty(strfind(P.VOI2,'.img')) ~isempty(strfind(P.VOI2,'.nii')) ~isempty(strfind(P.VOI2,'.mat'))])
                error{end+1}='Program will exit. VOI2 file does not a .img, .nii, or .mat file.';
            end
       end
    catch
       P.VOI2={};
    end
    
    %Check analysis
    try
        P.analysis=inputstructure.analysis;
        if ischar(P.analysis)
            if isempty(strfind(P.analysis,'psy')) && isempty(strfind(P.analysis,'phy')) && isempty(strfind(P.analysis,'psyphy'))
                error{end+1}='Program will exit. Analysis is not specified correctly, must be psy, phy, or psyphy.';
            end
        else
            error{end+1}='Program will exit. Analysis is not specified correctly, must be a string.';
        end
    catch
        error{end+1}='Program will exit. Analysis is not specified.'; 
    end
    
    %Check PPI method. Defaults to conditional.
    try
        P.method=inputstructure.method;
        if ischar(P.method)
            if isempty(strfind(P.method,'cond')) && isempty(strfind(P.method,'trad'))
                error{end+1}=('Program will exit. Method not correct.');
            end
        else
            error{end+1}=('Program will exit. Method must be a string.');
        end
    catch
        error{end+1}=('Program will exit. Method is not specified.');
    end

    %Check Compute Contrasts, 0 means do not compute contrasts
    if P.Estimate~=0
        try
            P.CompContrasts=inputstructure.CompContrasts;
            if isnumeric(P.CompContrasts)
                if mod(P.CompContrasts,1) || P.CompContrasts<0 || P.CompContrasts>4
                    P.CompContrasts=0;
                end
            else
                error{end+1}=('Program will exit. CompContrasts must be a number.');
            end
        catch
            P.CompContrasts=0;
        end
        
        %Check Contrasts
        if P.CompContrasts==1 || P.CompContrasts==3
            %Check Weighted, Defaults to treating all runs equally (e.g. 0)
            try
                P.Weighted=inputstructure.Contrasts.Weighted;
                if ~isnumeric(P.Weighted)
                    P.Weighted=0;
                end
            catch
                try 
                    P.Weighted;
                    if ~isnumeric(P.Weighted)
                        P.Weighted=0;
                    end
                catch
                    P.Weighted=0;
                end
            end
  
            %Check Contrasts
            try
                P.Contrasts=inputstructure.Contrasts;
                if ~isstruct(P.Contrasts) && iscell(P.Contrasts)
                    A = struct('left',[],'right',[],'name',[],'Weighted',[],'STAT',[],'c',[],'Contrail',[],'Prefix',[],'MinEvents',[],'MinEventsPer',[]);
                    A(1:numel(P.Contrasts))=A; clear t;
                    for ii=1:numel(P.Contrasts)
                        for jj=1:8
                            try
                                A(ii).left=P.Contrasts{ii}.left;
                            catch
                            end
                            try
                                A(ii).right=P.Contrasts{ii}.right;
                            catch
                            end
                            try
                                A(ii).name=P.Contrasts{ii}.name;
                            catch
                            end
                            try
                                A(ii).Weighted=P.Contrasts{ii}.Weighted;
                            catch
                            end
                            try
                                A(ii).STAT=P.Contrasts{ii}.STAT;
                            catch
                            end
                            try
                                A(ii).c=P.Contrasts{ii}.c;
                            catch
                            end
                            try
                                A(ii).Contrail=P.Contrasts{ii}.Contrail;
                            catch
                            end
                            try
                                A(ii).Prefix=P.Contrasts{ii}.Prefix;
                            catch
                            end
                            try
                                A(ii).MinEvents=P.Contrasts{ii}.MinEvents;
                            catch
                            end
                            try
                                A(ii).MinEventsPer=P.Contrasts{ii}.MinEventsPer;
                            catch
                            end
                        end
                    end
                    P.Contrasts=A;
                    for ii=1:numel(P.Contrasts)
                        if ~isfield(P.Contrasts(ii),'left') || isempty(P.Contrasts(ii).left) || ~iscell(P.Contrasts(ii).left)
                            error{end+1}=['P.Contrasts(' num2str(ii) ').left is not formatted correctly.'];
                        end
                        if ~isfield(P.Contrasts(ii),'right') || isempty(P.Contrasts(ii).right) || ~iscell(P.Contrasts(ii).right)
                            error{end+1}=['P.Contrasts(' num2str(ii) ').right is not formatted correctly.'];
                        end
                        if ~isfield(P.Contrasts(ii),'STAT') || (~strcmp(P.Contrasts(ii).STAT,'T') && ~strcmp(P.Contrasts(ii).STAT,'F'))
                            error{end+1}=['P.Contrasts(' num2str(ii) ').STAT is not formatted correctly.'];
                        end
                    end
                elseif isstruct(P.Contrasts) 
                    for ii=1:numel(P.Contrasts)
                        if (~isfield(P.Contrasts(ii),'left') || isempty(P.Contrasts(ii).left) || ~iscell(P.Contrasts(ii).left)) && ~(isnumeric(P.Contrasts(ii).left) && isempty(P.Contrasts(ii).right))
                            error{end+1}=['P.Contrasts(' num2str(ii) ').left is not formatted correctly.'];
                        end
                        if (~isfield(P.Contrasts(ii),'right') || isempty(P.Contrasts(ii).right) || ~iscell(P.Contrasts(ii).right)) && ~(isnumeric(P.Contrasts(ii).left) && isempty(P.Contrasts(ii).right))
                            error{end+1}=['P.Contrasts(' num2str(ii) ').right is not formatted correctly.'];
                        end
                        if ~isfield(P.Contrasts(ii),'STAT') || (~strcmp(P.Contrasts(ii).STAT,'T') && ~strcmp(P.Contrasts(ii).STAT,'F'))
                            error{end+1}=['P.Contrasts(' num2str(ii) ').STAT is not formatted correctly.'];
                        end
                    end
                else
                    invokecatchstatement
                end
                if numel(P.Contrasts)==1 && ((~isfield(P.Contrasts,'left') || isempty(P.Contrasts.left) || ~iscell(P.Contrasts.left)) || (~isfield(P.Contrasts,'right') || isempty(P.Contrasts.right) || ~iscell(P.Contrasts.right)));
                    if isfield(P.Contrasts,'left') && isnumeric(P.Contrasts.left) && isfield(P.Contrasts,'right') && isempty(P.Contrasts.right)
                    else
                        try
                            P.Contrasts=defContrasts([P.directory filesep 'SPM.mat'],P.Weighted,1,P.Contrasts.left);
                        catch
                            try
                                P.Contrasts=defContrasts([P.directory filesep 'SPM.mat'],P.Weighted,1,P.Contrasts.right);
                            catch
                                invokecatchstatement
                            end
                        end
                    end
                end
                for j=1:numel(P.Contrasts)
                    if ~isfield(P.Contrasts(j),'name') || isempty(P.Contrasts(j).name)
                        tmp=[P.Contrasts(j).left 'minus' P.Contrasts(j).right];
                        tmp=sprintf('%s_',tmp{:});
                        tmp=tmp(1:end-1);
                        P.Contrasts(j).name=tmp;
                        P.Contrasts(j).Weighted=P.Weighted;
                        clear tmp
                    end
                end
            catch
                try
                    if ~isempty(P.Tasks)
                        if isnumeric(P.Tasks{1})
                            tmpTasks=P.Tasks{2:end};
                        else
                            tmpTasks=P.Tasks;
                        end
                        P.Contrasts=defContrasts([P.directory filesep 'SPM.mat'],P.Weighted,1,tmpTasks);
                    else
                        P.Contrasts=defContrasts([P.directory filesep 'SPM.mat'],P.Weighted);
                    end
                    for j=1:numel(P.Contrasts)
                        tmp=[P.Contrasts(j).left 'minus' P.Contrasts(j).right];
                        tmp=sprintf('%s_',tmp{:});
                        tmp=tmp(1:end-1);
                        P.Contrasts(j).name=tmp;
                        P.Contrasts(j).Weighted=P.Weighted;
                        clear tmp
                    end
                catch
                    keyboard
                    if isnumeric(P.Contrasts(1).left)
                        disp(['Possible Error: [' num2str(P.Contrasts(1).left) ']'])
                    else
                        disp(['Possible Error: ' P.Contrasts(1).left])
                        error{end+1}=('Program will exit. Contrasts could not be defined automatically.');
                    end
                end
            end
        end
    end
    if isstruct(P.VOI)
        if ~isfield(inputstructure,'VOImin') || ~isnumeric(inputstructure.VOImin) || isempty(inputstructure.VOImin) || inputstructure.VOImin<0 
            P.VOImin=0;
        else
            P.VOImin=inputstructure.VOImin;
        end
        if P.VOI.exact && P.VOImin==0 % Fix this and P.VOI in the case that the field is missing...
            error{end+1}='P.VOImin cannot be zero if P.VOI.exact or P.VOI2.exact was set to 1';
        end
    end
    if isstruct(P.VOI2)
        if ~isfield(inputstructure,'VOImin') || ~isnumeric(inputstructure.VOImin) || isempty(inputstructure.VOImin) || inputstructure.VOImin<0 
            P.VOImin=0;
        else
            P.VOImin=inputstructure.VOImin;
        end
        if P.VOI2.exact && P.VOImin==0 % Fix this and P.VOI in the case that the field is missing...
            error{end+1}='P.VOImin cannot be zero if P.VOI.exact or P.VOI2.exact was set to 1';
        end
    end
    try
        save(structfile,'P')
    catch
        try
            save([pwd filesep inputstructure.subject '_PPIstructure.mat'],'P');
        catch
            save([pwd filesep 'invalidPPIstructure.mat'],'P');
        end
    end
    if isempty(error)
        try
            [P, tmperror]=PPPIinputsvalid(P);
            if isempty(tmperror)
                P.correct=1;
                error={};
            else
                error{end+1}=tmperror;
                error{end+1}='PPPI inputs were specified correctly, but failed error checking for dependencies. See errorval{1}.';
            end
        catch
            error{end+1}='PPPI inputs were specified correctly, but failed error checking for dependencies. This is likely a bug in the program.';
        end
    end
else
    try
        error{end+1}=[inputstructure 'does not contain a structure variable.']; 
        error{end+1}='PPPI_checkstruct.m will now exit. input structure file did not have a structure variable.';
        P={};
    catch
        P={};
        error{end+1}='PPPI_checkstruct.m will now exit. input structure was not a file with a structure.';
    end
end
if ~isempty(error)
    for ii=1:numel(error)
        disp(['ERROR ' num2str(ii) ': ' error{ii}])
    end
    return;
end


    
