function [Out,ICs] = detect_ICA_artefacts(S)

% Version 2.0 of ICA artifact detection for SPM8/12
%                  Rik.Henson@mrc-cbu.cam.ac.uk, March 2013, with thanks to Nitin Williams and Jason Taylor
%
% Updated 23-4-15 to allow IC's to be re-created from previously calculated IC weights (eg if want to use function to re-threshold)
% Updated 6-2-15 to gather outputs into single "Out" structure (and keep separate ICs identified by each reference)
%
% Needs EEGLAB on Matlab path to perform ICA
% Needs SPM/Fieldtrip on Matlab path to perform any filtering (though could be replaced by Matlab filtfilt)
%
% Main purpose is to return a FieldTrip/SPM TraMat to be applied later by,
% for example, spm_eeg_montage, to project out of data topographic patterns 
% related to artifacts, as defined by correlations between ICs and
% user-provided reference signals (eg recorded EOG/ECG). 
% 
% ICA is applied and those ICs that are artifacts are defined 1) temporally 
% by correlating with references, using boot-strapping testing if 
% necessary, and/or 2) spatially, by correlating with user-supplied
% topographies for classic artifactual sources, eg ocular or cardiac.
% Note that temporal correlation on its own can remove true brain signals
% that happen to be picked up by EOG/ECG electrodes, which is where spatial 
% correlation may help disambiguate. Conversely, some topographies (eg
% cardiac) seem to vary across individuals (others like blinks do not), 
% which is why temporal correlation helps. We have some standard Neuromag
% Vectorview Magnetometer and Planar Gradiometer topographies for blinks, 
% horizontal eye-movements and cardiac created by manual inespection and 
% averaging across people after transforming their data to Device space 
% using MaxFilter, which can be obtained from authors on request (assuming 
% you can also transform your data the same way). Relative (Z-scored) 
% thresholds are also helpful, assuming a given artifact does not split 
% into more than one IC, eg owing to movement half-way through a recording 
% - which is where MaxFilter's movement compensation also helps (Taylor & 
% Henson, 2014, Biomag conference).
%
% Fields of Out:
%   allrem   = cell array of IC numbers believed to be artifacts across all references 
%   bothrem  = cell array of IC numbers believed to be artifacts for each reference (spatial+temporal)
%   temprem  = T x N cell array of IC numbers believed to be artifacts for each of N temporal references, where:
%              T = 1 corresponds to single IC with max correlation
%              T = 2 corresponds to ICs with correlations surviving absolute p-value
%              T = 3 corresponds to ICs with correlations surviving relative Z-score across correlations
%              T = 4 (optional on Nperm below) corresponds to ICs with correlations surviving boot-strapped, phase-permuted p-value
%   spatrem  = T x N cell array of IC numbers believed to be artifacts for each of N spatial references, like for temrem above (except no bootstrap option)
%   weights  = full ICA weight matrix
%   TraMat   = channel trajectory matrix to project artifacts from data
%   temcor   = matrix of Pearson temporal correlations between each IC and 
%              each reference channel
%   spacor   = matrix of Pearson spatial correlations between each IC topo  
%              and each reference topography
%   compvars = vector of IC loadings 
%   varexpl  = vector of variance explained by each IC (simple function of compvars)
%
% Additional (large) output if requested:
%   all_remove = ICs that survive strictest definition of all thresholds for both temporal and spatial definitions
%   ICs        = ICs themselves! (may be big) (all optional outputs)
%
% Inputs (fields of argument structure "S"):
%   d          = Channel (Sensor) X Time matrix of EEG/MEG data to be corrected
%   refs.tem   = N cell array of 1 X Time vectors for N reference channels (eg EOG)
%   refs.spa   = N cell array of 1 X Channel vector for N reference topographies 
%                (order *must match* refs.tem; with empty cells if not available!)
%   PCA_dim    = number of PCs for data reduction priot to ICA (eg 60)
%   VarThr     = percentage of variance required to be an important artifactual IC (eg 1/PCA_dim, or could be 0) (DEFAULT possible)
%   FiltPars   = 1x2 or 1x3 matrix for low- or band- pass filtering before calculating correlation ([]=none=default), where third element is Sampling Rate
%   TemAbsPval = p-value threshold for absolute temporal correlation to define artifact (eg 0.01)
%   TemRelPval = additional Z-value for relative threshold of absolute temporal correlation to define artifact (eg 3)
%   SpaAbsPval = p-value threshold for absolute spatial correlation to define artifact (eg 0.01)
%   SpaRelPval = additional Z-value for relative threshold of absolute spatial correlation to define artifact (eg 2)
%   PermPval   = p-value for boot-strapped absolute temporal correlation to define artifact (eg 0.05)
%   Nperm      = number of permutations to create null distribution of Pearson (eg 1000)
%   thresholding = artifact thresholding (one or more of 1:4) (see T above)
%   remove     = artifacts to be removed ('all'|'both'|'temp'|'spat') (see fields *rem of Out above)
%
%   ICs        = (optional) if already have ICs and just want to re-check correlations with references
%   weights    = (optional) IC weights if want this function to recreate ICs (and also necessary if passed own ICs as above)
%   compvars   = (optional) IC weights if want this function to recreate ICs (and also necessary if passed own ICs as above)
%
% Could be extended to bootstrap spatial correlations too, though would
% need to do 2D FT to allow for spatial correlation in topographies...
%
% Could be extended to correlate with template power spectra (eg muscle artifact)

%% Inputs
try refs.tem = S.refs.tem;    catch, refs.tem = []; end
try refs.spa = S.refs.spa;    catch, refs.spa = []; end

if isempty(refs.tem) 
    if isempty(refs.spa)
        error('No temporal or spatial references for artifacts provided')
    else
        refs.tem = cell(1,length(refs.spa));
    end
% else
%     if isempty(refs.spa)
%         refs.spa = cell(1,length(refs.tem));
%     elseif length(refs.tem) ~= length(refs.spa)
%         error('Number of temporal references needs to match number of spatial references')
%     end
end

Nrefs(1) = length(refs.tem);
Nrefs(2) = length(refs.spa);

try TemAbsPval = S.TemAbsPval;  catch, TemAbsPval = 1; end
try SpaAbsPval = S.SpaAbsPval;  catch, SpaAbsPval = 1; end
try TemRelZval = S.TemRelZval;  catch, TemRelZval = 0; end
try SpaRelZval = S.SpaRelZval;  catch, SpaRelZval = 0; end

if isfield(S,'ICs')
    try
        weights = S.weights;
        compvars = S.compvars;
    catch
        error('If passing pre-computed ICs, then must also pass weights and compvars')
    end
end

try d = S.d; catch, error('Must provide Channel x Time MEG data matrix'); end

if isfield(S,'weights')
    try
        compvars = S.compvars;
    catch
        error('If passing pre-computed IC weights, then must also pass compvars')
    end
    try
        S.ICs = S.weights * d;
    catch
        error('Pre-specified IC weights do not match data size')
    end
end

try PCA_dim = S.PCA_dim;    catch, PCA_dim = round(0.75*size(d,1)); end
try VarThr = S.VarThr;      catch, VarThr = 0; end % could default to 100/PCA_dim?
try FiltPars = S.FiltPars;  catch, FiltPars = []; end

try Rseed = S.Rseed;        catch, Rseed = []; end  % If want reproducibility

% For permutation testing of temporal correlations (to turn off, pass S.Nperm = 0)
try PermPval = S.PermPval;  catch, PermPval = .05/(2*PCA_dim);                   end % Kind of 2-tailed Bonferonni!
try Nperm = S.Nperm;        catch, Nperm = round((4*sqrt(CorPval)/CorPval)^2);   end % http://www.epibiostat.ucsf.edu/biostat/sen/statgen/permutation.html#_how_many_permutations_do_we_need


%% Initialize
Out = [];
Nsamp = size(d,2);
if rem(Nsamp,2)==0  % Must be odd number of samples for fft below
    d     = d(:,2:Nsamp);
    for r=1:length(refs.tem)
        refs.tem{r}  = refs.tem{r}(:,2:Nsamp);
    end
    Nsamp = Nsamp-1;
end


%% Run ICA
try
    ICs = S.ICs;
    ICs = ICs(:,1:Nsamp);
catch
    try
        [weights,sphere,compvars,bias,signs,lrates,ICs] = rik_runica(d,'pca',PCA_dim,'extended',1,'maxsteps',800,'rseed',Rseed); % Just local copy where rand seed can be passed
    catch
        if ~isempty(Rseed); warning('Random seed requested but no facility with standard runica?'); end
        [weights,sphere,compvars,bias,signs,lrates,ICs] = runica(d,'pca',PCA_dim,'extended',1,'maxsteps',800);
    end
end
Out.weights = weights;

%% Filtering (if any) (and transposition for speed)
if length(FiltPars) == 3
    fprintf('Bandpass filtering from %d to %d Hz (warning - can fail)\n',FiltPars(1), FiltPars(2));
    for r=1:Nrefs(1)
        refs.tem{r} = ft_preproc_bandpassfilter(refs.tem{r}, FiltPars(3), FiltPars(1:2),  [], 'but','twopass','reduce');
    end
    ICs  = ft_preproc_bandpassfilter(ICs,  FiltPars(3), FiltPars(1:2),  [], 'but','twopass','reduce')';
elseif length(FiltPars) == 2
    fprintf('Lowpass filtering to %d Hz\n',FiltPars(1));
    for r=1:Nrefs(1)
        refs.tem{r} = ft_preproc_lowpassfilter(refs.tem{r}, FiltPars(2), FiltPars(1),  5, 'but','twopass','reduce');
    end
    ICs  = ft_preproc_lowpassfilter(ICs,  FiltPars(2), FiltPars(1),  5, 'but','twopass','reduce')';
else
    ICs  = ICs'; % Faster if pre-transpose once (check tic;toc)
end
% figure; for i=1:Nrefs; plot(refs.tem{r}'); title(i); pause; end
% figure; for i=1:PCA_dim; plot(ICs(:,i)); title(i); pause; end
% figure; for i=1:PCA_dim; plot(ICs(30000:40000,i)); title(i); pause; end

iweights  = pinv(Out.weights);
Out.temcor = zeros(Nrefs(1),PCA_dim); tempval = zeros(Nrefs(1),PCA_dim);
Out.spacor = zeros(Nrefs(2),PCA_dim); spapval = zeros(Nrefs(2),PCA_dim);
%reltemcor = zeros(Nrefs(1),PCA_dim); relspacor = zeros(Nrefs(2),PCA_dim);
temremove  = cell(4,Nrefs(1)); sparemove  = cell(3,Nrefs(2));

for r = 1:Nrefs(1)
    
    %% Check temporal correlation with any reference channels
    for k = 1:PCA_dim
        [Out.temcor(r,k),tempval(r,k)] = corr(refs.tem{r}',ICs(:,k));
    end
    
    [~,temremove{1,r}] = max(abs(Out.temcor(r,:)));
    
    temremove{2,r} = find(tempval(r,:) < TemAbsPval);
    
    temremove{3,r} = find(abs(zscore(Out.temcor(r,:))) > TemRelZval);
    
    if Nperm > 0
        permcor = zeros(1,PCA_dim);
        maxcor  = zeros(Nperm,1);
        
        ff = fft(refs.tem{r}',Nsamp);
        mf = abs(ff);
        wf = angle(ff);
        hf = floor((length(ff)-1)/2);
        rf = mf;
        
        for l = 1:Nperm % could parfor...
            rf(2:hf+1)=mf(2:hf+1).*exp((0+1i)*wf(randperm(hf)));    % randomising phases (preserve mean, ie rf(1))
            rf((hf+2):length(ff))=conj(rf((hf+1):-1:2));            % taking complex conjugate
            btdata = ifft(rf,Nsamp);                                % Inverse Fourier transform
            
            for k = 1:PCA_dim
                permcor(k) = corr(btdata,ICs(:,k));
            end
            maxcor(l) = max(abs(permcor));
            fprintf('.');
        end
        fprintf('\n')
        %         figure,hist(maxcor)
        
        temremove{4,r} = find(abs(Out.temcor(r,:)) > prctile(maxcor,100*(1-PermPval)));
    end
end

for r = 1:Nrefs(2)        
    %% Check spatial correlation with any reference channels
        for k = 1:PCA_dim
            [Out.spacor(r,k),spapval(r,k)] = corr(refs.spa{r}',iweights(:,k));
        end
        
        [~,sparemove{1,r}] = max(abs(Out.spacor(r,:)));
        
        sparemove{2,r} = find(spapval(r,:) < SpaAbsPval);

        sparemove{3,r} = find(abs(zscore(Out.spacor(r,:))) > SpaRelZval);
end

%% Variance Thresholding
Out.compvars = compvars;
Out.varexpl = 100*compvars/sum(compvars);
varenough   = find(Out.varexpl > VarThr);

Out.temrem = temremove;
Out.sparem = sparemove;

thresholding = S.thresholding;
if ~Nperm, thresholding(thresholding == 4) = []; end % remove thresholding if not permutation
for r = 1:Nrefs(1)
    remove = 1:PCA_dim;
    for t = thresholding
        remove = intersect(remove,temremove{t,r});
    end
    remove = intersect(remove,varenough);  % plus sufficient Variance Explained (if required)
    Out.bothrem{1,r} = remove;
end

thresholding(thresholding == 4) = []; % there is no option for permutation in spatial
for r = 1:Nrefs(2)
    remove = 1:PCA_dim;
    for t = thresholding
        remove = intersect(remove,sparemove{t,r});
    end    
    remove = intersect(remove,varenough);  % plus sufficient Variance Explained (if required)
    Out.bothrem{2,r} = remove;
end

Out.allrem = 1:PCA_dim;
if Nrefs(1) > 0
    Out.allrem = intersect(Out.allrem,cat(2,Out.bothrem{1,:}));
end
if Nrefs(2) > 0
    Out.allrem = intersect(Out.allrem,cat(2,Out.bothrem{2,:}));
end
 
toremove = Out.([S.remove 'rem']);

if ~isempty(toremove)
    finalics  = setdiff(1:PCA_dim,toremove);
    Out.TraMat    = iweights(:,finalics) * Out.weights(finalics,:);
else
    Out.TraMat    = eye(size(d,1));
end

return

%% Some code for manual inspection, assuming D structure loaded from SPM file
% f=figure('color',[1 1 1],'deleteFcn',@dFcn);
% in.ParentAxes = axes('parent',f);
% in.f = f;
% in.type ='MEGMAG'; 
% %in.type ='MEGPLANAR';
% chanind = D.indchantype(in.type);
% 
% d = D(chanind,:);
% w = Out.weights; iw = pinv(w);
% d = w*d;
% 
% %toinspect = remove; % !!! If just want to look at sig correlations with EOG/ECG
% toinspect = 1:size(d,1);       % !!! If want to look at all
% 
% for i=1:length(toinspect)
%     ii = toinspect(i);
%     fprintf('Component %d (%d)\n',ii,i)
%     figure(3); plot(d(ii,:)); title(ii); zoom on
% %    figure(4); plot([0:4:12000],d(ii,80000:83000)); title(ii); % zoom in on 12sec section
%     [dum1,pow,dum2]=pow_spec(d(ii,:)',1/250,1,0,5); axis([0 80 min(pow) max(pow)]);
%     spm_eeg_plotScalpData(iw(:,ii),D.coor2D(chanind),D.chanlabels(chanind),in);
%     pause
% end

                        