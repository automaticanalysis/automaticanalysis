function RunSPM8ana(F,opt,pars)
%%% This is a wrapper around the SPM8 batch tools for running second level
%%% analyses.  The input is the output from CreateDesign.m (plus a Scans
%%% field).  Just a simple way of running second level SPM8 analyses
%%%
%%% Inputs:  
%%%
%%% F = Design Structure from CreateDesign.m. There is also an
%%% additional optinon to specify an F.con structure with an spm8 job
%%% mananger style contrast specificaiton field that will generate the
%%% specified contrasts
%%%
%%% opt = an option flag to open the spm batch confuration GUI with the
%%% specified design and scans fille in.  1 = Yes. 0 = No
%%%
%%% pars = spm second level settings.  If you don't want to use the default
%%% pars.? settings below you can make you own pars data structure and pass
%%% it in.
%%%
%%% Written by Aaron Schultz - aschultz@martinos.org
%%% Copyright (C) 2011,  Aaron P. Schultz
%%%
%%% Supported in part by the NIH funded Harvard Aging Brain Study (P01AG036694) and NIH R01-AG027435 
%%%
%%% This program is free software: you can redistribute it and/or modify
%%% it under the terms of the GNU General Public License as published by
%%% the Free Software Foundation, either version 3 of the License, or
%%% any later version.
%%% 
%%% This program is distributed in the hope that it will be useful,
%%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%% GNU General Public License for more details.

in = F.IN;
if ~exist('opt');
    opt = 0;
end
if ~exist('pars');
    pars.dept = 0;  %% Assume Dependence?       0=no or 1=yes
    pars.var = 1;   %% Assume Unequal Variance?   0=no or 1=yes
    pars.tm.tm_none = 1;    %% No Threshold Masking
%     pars.tm.tma.athresh= 100;   %% Absolute theshold masking
%     pars.tm.tmr.rthresh= 0.80;  %% Relative theshold masking
    pars.im = 0;    %% Implicit Masking   1 = yes, 0 = no
    pars.em = {''}; %% Explicit Masks     list mask files e.g. {'/autofs/space/plato_002/users/APS_MATLAB/spm8/templates/epi.nii,1'};
    
    pars.g_omit = 1;    %% Global Calculation:   1 is the default setting
    pars.gmsca.gmsca_no = 1;  %% No Grand Mean Scaling.
%     pars.gmsca.gmsca_yes.gmscv = 50;
    pars.glonorm = 1;   %% Global Normalization.  1 = No. 2 = Proportional 3 = ANCOVA
end

try
    job.dir{1} = in.swd;
catch
    job.dir{1} = pwd;
end
%% Flexible Factorial
cc = 0;
if ~isempty(in.Within);
    
    if ~isempty(in.Between);
        for ii = 1:length(in.Between)
            cc = cc+1;
            job.des.fblock.fac(cc).name = in.FactorLabs{cc};
            job.des.fblock.fac(cc).dept = 1-in.Independent(cc);
            job.des.fblock.fac(cc).variance = 1-in.EqualVar(cc);
            job.des.fblock.fac(cc).gmsca = 0;
            job.des.fblock.fac(cc).ancova = 0;
        end
    end

    for ii = 1:length(in.Within)
        cc = cc+1;
        job.des.fblock.fac(cc).name = in.FactorLabs{cc};
        job.des.fblock.fac(cc).dept = 1-in.Independent(cc);
        job.des.fblock.fac(cc).variance = 1-in.EqualVar(cc);
        job.des.fblock.fac(cc).gmsca = 0;
        job.des.fblock.fac(cc).ancova = 0;
    end
    
    cc = cc+1;
    job.des.fblock.fac(cc).name = 'Subjects';
    job.des.fblock.fac(cc).dept = 0;
    job.des.fblock.fac(cc).variance = 0;
    job.des.fblock.fac(cc).gmsca = 0;
    job.des.fblock.fac(cc).ancova = 0;
    
else 
    if ~isempty(in.Between);
        for ii = 1:length(in.Between)
            cc = cc+1;
            job.des.fblock.fac(cc).name = in.FactorLabs{cc};
            job.des.fblock.fac(cc).dept = 1-in.Independent(cc);
            job.des.fblock.fac(cc).variance = 1-in.EqualVar(cc);
            job.des.fblock.fac(cc).gmsca = 0;
            job.des.fblock.fac(cc).ancova = 0;
        end
    end
end

if ~isempty(in.Within);
% %     I = [ones(size(F.FM,1),1) F.FM];
% %     I(:,end+1:4) = 1;
%     I = [ones(size(F.FM,1),1) F.FM(:,2:end) F.FM(:,1)];
%     I(:,end+1:4) = 1;
    I = [F.FM F.FM(:,1)];
    I(:,end+1:4) = 1;
else
    I = F.FM;
    I(:,end+1:4) = 1;
end

job.des.fblock.fsuball.specall.scans = F.Scans;
job.des.fblock.fsuball.specall.imatrix = I;

if isempty(in.Within)
    job.des.fblock.maininters = [];
    for ii = 1:length(job.des.fblock.fac)
        job.des.fblock.maininters{end+1}.fmain.fnum = ii;
    end
    
    if numel(in.Interactions)>0
        for ii = 1:length(in.Interactions)
            job.des.fblock.maininters{end+1}.inter.fnums = in.Interactions{ii};
        end
    end
else
    job.des.fblock.maininters = [];
    for ii = 1:length(job.des.fblock.fac)-1
        job.des.fblock.maininters{end+1}.fmain.fnum = ii;
    end
    
    if numel(in.Interactions)>0
        for ii = 1:length(in.Interactions)
            job.des.fblock.maininters{end+1}.inter.fnums = in.Interactions{ii};
        end
    end
    
    %%do subjects stuff
    if numel(in.Within)==1
        job.des.fblock.maininters{end+1}.fmain.fnum = length(job.des.fblock.fac);
    elseif numel(in.Within)==2
        job.des.fblock.maininters{end+1}.fmain.fnum = length(job.des.fblock.fac);
%         for ii = numel(in.Between)+1:numel(in.Within);
%             job.des.fblock.maininters{end+1}.inter.fnums = [ii length(job.des.fblock.fac)];
%         end
    elseif numel(in.Within)>2
        error('This script is not designed to handle more than two repeated factors');
    end
    
end

if isfield(F, 'cov')
    job.cov = F.cov;
else
    job.cov = struct([]);
    % job.cov.c=[];
    % job.cov.cname=[];
    % job.cov.iCFI=[];
    % job.cov.iCC=[];
end

job.masking.tm = pars.tm;
job.masking.im = pars.im;
job.masking.em = pars.em;

job.globalc.g_omit = pars.g_omit;

job.globalm.gmsca = pars.gmsca;
job.globalm.glonorm = pars.glonorm;


job2.spmmat{1} = [job.dir{1} filesep 'SPM.mat'];
job2.method.Classical = 1;

if isfield(F,'con')
    job3.spmmat = {[job.dir{1} filesep 'SPM.mat']};
    for ii = 1:length(F.con.name)
        if F.con.type{ii} == 't'
            job3.consess{ii}.tcon.name = F.con.name{ii};
            job3.consess{ii}.tcon.convec = F.con.vect{ii};
            job3.consess{ii}.tcon.sessrep = 'none';
        end
        if in.con.type{ii} == 'f'
            job3.consess{ii}.fcon.name = F.con.name{ii};
            job3.consess{ii}.fcon.convec = {F.con.vect{ii}};
            job3.consess{ii}.fcon.sessrep = 'none';
        end
    end
    job3.delete = 1;
else
    job3 = [];
end



if opt == 1;
    cfg_ui;
    clear batch;
    batch{1}.spm.stats.factorial_design = job;
    batch{2}.spm.stats.fmri_est = job2;
    if isfield(in,'con')
        batch{3}.spm.stats.con = job3;
    end
    spm_jobman('interactive', batch); 
    return
end

save jobs.mat job job2 job3;

spm_run_factorial_design(job);
spm_run_fmri_est(job2);
if isfield(in,'con')
    spm_run_con(job3);
end


