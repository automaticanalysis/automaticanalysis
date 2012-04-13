function [P,error]=PPPIinputsvalid(P)
% Checks for file dependencies
% Do not run in isolation
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
%   Last modified on 11/6/2010 by Donald G. McLaren
%   (mclaren@nmr.mgh.harvard.edu)
%   Wisconsin Alzheimer's Disease Research Center - Imaging Core, Univ. of
%   Wisconsin - Madison
%   Neuroscience Training Program and Department of Medicine, Univ. of
%   Wisconsin - Madison
%   GRECC, William S. Middleton Memorial Veteren's Hospital, Madison, WI
%   GRECC, Bedford VAMC
%   Department of Neurology, Massachusetts General Hospital and Havard
%       Medical School

error=[];
% New inputs....
[region1,region2]=strtok(P.Region);
region2=region2(2:end);
cd(P.directory);
load SPM.mat

%% Build a log file
tmp=clock; time{1}=num2str(tmp(1)); time{2}=num2str(tmp(2));time{3}=num2str(tmp(3)); clear tmp
if exist([P.subject '_PPPI_' time{2} '_' time{3} '_' time{1} '.log'],'file')
    logs=dir([P.subject '_PPPI_' time{2} '_' time{3} '_' time{1} '*.log']);
    R=length(logs)+1;
    logfile=[P.subject '_PPPI_' time{2} '_' time{3} '_' time{1} '_' num2str(R) '.log'];
else
    logfile=[P.subject '_PPPI_' time{2} '_' time{3} '_' time{1} '.log'];
end
disp(['Log File: ' logfile])
diary(logfile)
diary on
disp(' ')
disp('PPPI Version: 7.6-22-2011')
disp(' ')
disp(' ')
disp('Parameters used:')
disp(['     Processing subject: ' P.subject])
try 
    disp(['     VOI file          : ' P.VOI])
catch
    disp(['     VOI file          : ' P.VOI.VOI])
end
if ~isempty(P.VOI2)
    disp(['     VOI2 file         : ' P.VOI2])
    if isempty(region2)
        error='Program will exit. VOI2 has no region name associated with it.';  %%%%%%%%%%%Error Check
        return;
    end
    if strcmp(P.analysis,'phy')
        disp(['     Output filenames will follow the following format:' ' ' P.subject '_' region1 '_and_' region2 '_session#_PPI_regressors.txt'])
    elseif strcmp(P.analysis,'psyphy')
        disp(['     Output filenames will follow the following format:' ' ' P.subject '_' region1 '_and_' region2 '_session#_' P.method '_PPPI_regressors.txt'])
    else
        disp('Program will exit. Only one VOI can be specified for Psychophysiological Interactions')
        return;
    end
else
    if strcmp(P.analysis,'psy')
        disp(['     Output file will be:' ' ' P.subject '_' region1 '_session#_' P.method '_PPI_regressors.txt'])
    else
        disp('Program will exit. VOI2 was not specified and is required for selected analysis.')
        return;
    end
end
if isfield(P,'contrast')
    if isnumeric(P.contrast) && P.contrast==0
        disp('     Contrast          : No adjustment')
    elseif isnumeric(P.contrast) && P.contrast<=length(SPM.xCon)
        disp(['     Contrast          : ' SPM.xCon(P.contrast).name])
        P.contrast=P.contrast;
    elseif isnumeric(P.contrast) && P.contrast>length(SPM.xCon)
        error=['Program will exit. Your contrast: ' num2str(P.contrast) ' does not exist.'];
        return
    else
        if iscellstr(P.contrast)
            match=0;
            for i=1:length(SPM.xCon)
                if strcmp(P.contrast,SPM.xCon(i).name)
                    disp(['     Contrast          : ' SPM.xCon(i).name])
                    P.contrast=i;
                    match=1;
                    break;
                end
            end
            if ~match
                for i=1:length(SPM.xCon)
                    if strcmpi('Omnibus F-test for PPI Analyses',SPM.xCon(i).name)
                        disp(['     Contrast          : ' SPM.xCon(i).name])
                        P.contrast=i;
                        match=1;
                        break;
                    end
                end
            end
            if ~match
                try
                    P.contrast=defContrasts(SPM,0,-1);
                    xCon = spm_FcUtil('Set',P.contrast.name,P.contrast.STAT,'c',P.contrast.c,SPM.xX.xKXs);
                    init=length(SPM.xCon);
                    if init~=0
                        SPM.xCon(init+1) = xCon;
                    elseif init==0
                        SPM.xCon = xCon;
                    else
                        triggercatchstatement
                    end
                    SPM = spm_contrasts(SPM,init+1);
                    P.contrast=init+1;
                    disp(['     Contrast          : ' SPM.xCon(P.contrast).name])
                catch
                    error='Program will exit. Contrast cannot be defined.';
                    return;
                end
            end
        else
            error='Program will exit. This should not be possible.';
            return;
        end
    end
    
    if strcmp(P.analysis,'psy')
        disp('     Analysis          : Psychophysiological Interactions');
    elseif strcmp(P.analysis,'phy')
        disp('     Analysis          : Physiophysiological Interactions');
    else
        disp('     Analysis          : Psychophysiophysiological Interactions');
    end
    
    if strcmp(P.extract,'eig')
        disp('     Extraction        : eigenvariate')
    else
        disp('     Extraction        : mean')
    end
    
    if length(P.Tasks)>str2double(max(SPM.xsDes.Trials_per_session))
        error='Program will now exit. You specified too many tasks.';
        return
    elseif isempty(P.Tasks)
        disp('     Tasks             : ALL')
    else
        tmp=['     Tasks             : ' P.Tasks];
        tmp=sprintf('%s_',tmp{:});
        tmp=tmp(1:end-1);
        disp(tmp)
        clear tmp
    end
    if strcmp(P.method,'trad') && ~strcmp(P.analysis,'phy')
        if isempty(P.Weights)
            error='Program will exit. Weight vector not specified.';
            return;
        else
            if length(P.Weights)~=length(P.Tasks) && length(P.Weights)~=max(str2num(SPM.xsDes.Trials_per_session))
                error=(['Program will exit. Length of Weight vector is not equal to number of tasks'...
                    'OR the number of tasks varies from session to session, which is not allowed' ...
                    '  with the traditional approach.']);
                return;
            end
        end
        disp('     Method            : Traditional')
        disp(['     Weight vector     : ' num2str(P.Weights)])
    elseif ~isempty(P.Tasks)
        if str2double(P.Tasks{1})~=0 && str2double(P.Tasks{1})~=1
            error=('Program will exit. Tasks must begin with a 1 or 0 to specify if all tasks are needed or not.');
            return;
        end
        disp('     Method            : Condition Specific')
    else
        disp('     Method            : Condition Specific')
    end
    
    if isempty(P.Tasks)
        tasks={};
        for ss=1:numel(SPM.Sess)
            for tt=1:numel(SPM.Sess(ss).U)
                tasks=[tasks SPM.Sess(ss).U(tt).name];
            end
        end
        [Tasks]=unique(tasks);
        if strcmp(P.method,'trad')
            P.Tasks=Tasks;
        else
            P.Tasks=['0' Tasks];
        end
    end
end
diary off