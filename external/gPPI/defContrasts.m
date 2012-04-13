function Contrasts=defContrasts(SPM,Weighted,opt,Tasks)
% This function builds T contrasts structure or the omibus F-contrast.
% opt==1, compute omnibus for PPI adjustment
% opt==2, ignore task number restriction
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

if nargin<4
    Tasks={};
end
if nargin<3
    opt=0;
end
if nargin<2
    Weighted=0;
end
if ~isstruct(SPM)
    if exist(SPM,'file')
        load(SPM);
    else
        error('SPM.mat does not exist.');
    end
end

%% Estimate Omnibus F-contrast for PPI Adjustment
% Handles temporal and dispersion derivatives
if opt==-1
    tasks={};
    dur=[];
    trials=[];
    for ss=1:numel(SPM.Sess)
        for tt=1:numel(SPM.Sess(ss).U)
            tasks=[tasks SPM.Sess(ss).U(tt).name];
            dur(end+1)=max(SPM.Sess(ss).U(tt).dur)<Weighted;
            trials(end+1)=length(SPM.Sess(ss).U(tt).ons);
        end
    end
    [Tasks]=unique(tasks);

    for ii=1:length(Tasks)*3
        Contrasts.left{ii,1}=Tasks{ceil(ii/3)};
        Contrasts.right{ii,1}={'none'};    
    end
    Contrasts.STAT='F';
    Contrasts.name='Omnibus F-test for PPI Analyses';
    Contrasts.c=createVec(Contrasts.left,Contrasts.right,SPM,Weighted,' ',repmat({'bf(1)';'bf(2)';'bf(3)'},length(Tasks),1))';
    return;
end

%% Define Combinations
if ~iscellstr(Tasks) || isempty(Tasks)
    Tasks={};
    for ii=1:numel(SPM.Sess)
        for jj=1:numel(SPM.Sess(ii).U)
            for kk=1:numel(SPM.Sess(ii).U(jj).name)
                Tasks(end+1)=SPM.Sess(ii).U(jj).name(kk);
            end
        end
    end
else 
    if size(Tasks,1)~=1 && size(Tasks,2)~=1
        disp('Tasks formatted incorrectly');
        return;
    end
end
Tasks=unique(Tasks);
numTasks=numel(Tasks);
if opt==0
    if numTasks>4
        Contrasts.left={'More than four tasks detected. Please specify the tasks.'};
        return;
    elseif numTasks==0
        Contrasts.left={'No Tasks Detected. Something went wrong.'};
        return;
    end
elseif numTasks==0
    Contrasts.left={'No Tasks Detected. Something went wrong.'};
    return;
end

% Generate all Possible Combinations unless opt=1 OR contrast vectors are defined
c=mat2cell(repmat([-1 0 1],1,numTasks),1,3*ones(1,numTasks)); % Each task can either be positive, negative, or 0 weighted
n = length(c);
[c{:}]=ndgrid(c{:});
c=cat(n+1,c{:});
c = reshape(c,[],n);

ind=[];
for i=1:size(c,1)
    Contrasts(i).STAT='T';
    Contrasts(i).left=Tasks(c(i,:)>0);
    if isempty(Contrasts(i).left)
        Contrasts(i).left={'none'};
    end
    Contrasts(i).right=Tasks(c(i,:)<0);
    if isempty(Contrasts(i).right)
        Contrasts(i).right={'none'};
    end
    tmp=[Contrasts(i).left 'minus' Contrasts(i).right];
    tmp=sprintf('%s_',tmp{:});
    tmp=tmp(1:end-1);
    Contrasts(i).name=tmp;
    if strcmp(Contrasts(i).left,Contrasts(i).right)
        ind(end+1)=i;
    end
end
Contrasts(ind)=[];
end
