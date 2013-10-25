function [B,D] = do_mumford(Y,TR,motion,onsets,durations)
%
% estimate betas for each trial independently using the Mumford approach
%
% Requires:
% 
% Y 
%
% fMRI data as an N x K matrix where N is the number of TRs and K the number of voxels
%
% TR
% 
% Repetition time
%
% motion
%
% An NxK (typically K=6) matrix of motion regressors (the rp_*.txt file in SPM)
%
% onsets
% 
% An Mx1 vector of onsets for all stimuli; slice-timing is assumed to
% be taken into account. E.g. if slice 15 is the reference slice for 30 slices then 
% onsets = oldonsets - 15*TR/30
%
% durations
%
% An Mx1 vector of durations for each of the stimuli
%
% output
%
% M x K Beta estimates for each of the trials and each of the voxels
% D example of the design matrix constructed for the final trial
%
% Test:
% 
% Y = randn(200,10);
% TR = 2;
% onsets = (0:5:345)';
% durations = ones(70,1);
% motion = 1e-1*randn(200,6);
% [B,D] = do_mumford(Y,TR,motion,onsets,durations)

 res = 0.1; % time resolution (just choose a small number << TR)
    
 N = size(Y,1); % number of volumes
 M = numel(onsets); % number of trials

 Y = zscore(Y); % create zero mean and unit SD data
 
 % estimate betas
 B = [];
 for i=1:M
      
   fprintf('processing trial %d of %d\n',i,M);
   
   % set time resolution and index for TR in that vector
   t = 0:res:(N*TR);
   ttr = closest(0:TR:((N-1)*TR),t);
   
   % generate stimulus matrix (stick functions)
   S = zeros(numel(t),2);
   
   % create first regressor
   tbeg = find(t >= onsets(i),1);
   tend = find(t > onsets(i) + durations(i),1);
   S(tbeg:tend,1) = 1;
     
   % create second regressor
   for k=[1:(i-1) (i+1):M]
     tbeg = find(t >= onsets(k),1);
     tend = find(t > onsets(k) + durations(k),1);
     S(tbeg:tend,2) = 1;
   end
   
   % create default HRF with res second resolution
   xBF.dt = res;
   xBF.name = 'hrf';
   %xBF.length = 8;
   %xBF.order = 4;
   hrf = spm_get_bf(xBF);
   hrf = hrf.bf;
   nbasis = size(hrf,2);
   
   % create regressors by convolving with HRF
   X = zeros(N,size(S,2)*nbasis);
   for j=1:size(S,2)
     for k=1:nbasis
       xc = conv(S(:,j),hrf(:,k));
       xc = xc(1:numel(t));
       X(:,(j-1)*nbasis+k) = xc(ttr);
     end
   end
   
   % downsample stimulus matrix
   S = S(ttr,:);

   % create DCT basis functions
   dur = max(durations); % high-pass filter cutoff in seconds (max(128,4*max block length))
   if isempty(dur), dur = 0; end
   HPFc = max(128,4*dur); % set high-pass filter cutoff
   HPFk = fix(min(2*(N*TR)/HPFc,N/2)); % determine high-pass filter order
   % note: SPM first takes out these contributions and only then
   % estimates the remainder of the design matrix wrt the residual
   % This is different from our approach; no idea what works best
   C = spm_dctmtx(N,HPFk);
   
   % create final design matrix;
   D = [X motion C];
   
   % least squares estimate of the model
   if size(D,2) > size(D,1)
     fprintf('p > n; using pseudo-inverse to compute least-squares solution\n');
     Betas = pinv(D) * Y;
   else
     Betas = D\Y;
   end
   
   B = [B; Betas(1,:)];
   
 end
 
end