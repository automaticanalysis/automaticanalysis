function [remove,weights,TraMat,temcor,spacor,varexpl,ICs] = detect_ICA_artefacts(S);

% Version 2.0 of ICA artifact detection for SPM8/12
%
% Main purpose is to return a FieldTrip/SPM TraMat to be applied later by,
% for example, spm_eeg_montage, to project out of data topographic patterns 
% related to artifacts. 
% 
% ICA is applied and those ICs that are artifacts are defined 1) temporally 
% by correlating with recorded EOG or ECG, using boot-strapping testing if 
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
% Outputs:
%   remove   = cell array of IC numbers believed to be artifacts
%   weights  = full ICA weight matrix
%   TraMat   = channel trajectory matrix to project artifacts from data
%   temcor   = matrix of Pearson temporal correlations between each IC and 
%              each reference channel
%   spacor   = matrix of Pearson spatial correlations between each IC topo  
%              and each reference topography
%   varexpl  = vector of variance explained by each IC
%   ICs      = ICs themselves! (may be big) (all optional outputs)
%
% Inputs (fields of argument structure "S"):
%   d          = Channel (Sensor) X Time matrix of EEG/MEG data to be corrected
%   refs.tem    = N cell array of 1 X Time vectors for N reference channels (eg EOG)
%   refs.spa    = N cell array of 1 X Channel vector for N reference topographies 
%                (order *must match* refs.tem; with empty cells if not available!)
%   PCA_dim    = number of PCs for data reduction priot to ICA (eg 60)
%   VarThr     = percentage of variance required to be an important artifactual IC (eg 1/PCA_dim, or could be 0) (DEFAULT possible)
%   FiltPars   = 1x2 or 1x3 matrix for low- or band- pass filtering before calculating correlation ([]=none=default), where third element is Sampling Rate
%   TemAbsPval = p-value threshold for absolute temporal correlation to define artifact (eg 0.01)
%   TemRelPval = Z-value for relative threshold of absolute temporal correlation to define artifact (eg 3)
%   SpaAbsPval = p-value threshold for absolute spatial correlation to define artifact (eg 0.01)
%   SpaRelPval = Z-value for relative threshold of absolute spatial correlation to define artifact (eg 2)
%   PermPval   = p-value for boot-strapped absolute temporal correlation to define artifact (e.g, 0.05)
%   Nperm      = number of permutations to create null distribution of Pearson (eg 1000)
%
%   ICs        = (optional) if already have ICs and just want to re-check correlations with references
%   weights    = IC weights; only necessary if passed own ICs as above
%   compvars   = IC variances; only necessary if passed own ICs as above
%
% Needs EEGLAB on Matlab path to perform ICA
% Needs SPM/Fieldtrip on Matlab path to perform any filtering (though could be replaced by Matlab filtfilt)
%
% Could be extended to bootstrap spatial correlations too, though would
% need to do 2D FT to allow for spatial correlation in topographies...
%
% Could be extended to correlate with template power spectra (eg muscle artifact)
%
% Rik.Henson@mrc-cbu.cam.ac.uk, March 2013, with thanks to Nitin Williams and Jason Taylor

%% Inputs
try TemAbsPval = S.TemAbsPval;  catch TemAbsPval = 1; end
try SpaAbsPval = S.SpaAbsPval;  catch SpaAbsPval = 1; end
try TemRelZval = S.TemRelZval;  catch TemRelZval = 0; end
try SpaRelZval = S.SpaRelZval;  catch SpaRelZval = 0; end
  
if TemAbsPval == 1 & SpaAbsPval == 1 & TemRelZval == 0 & SpaRelZval == 0 
    error('Must have at least one temporal or spatial absolute or relative threshold for artifacts')
end

try refs.tem = S.refs.tem;    catch refs.tem = []; end
try refs.spa = S.refs.spa;    catch refs.spa = []; end

if ~isempty(refs.tem) & TemAbsPval == 1 & TemRelZval == 0
    error('No thresholds for temporal reference channels')
end
if ~isempty(refs.spa) & SpaAbsPval == 1 & SpaRelZval == 0
    error('No thresholds for spatial reference topographies')
end

if isempty(refs.tem) 
    if isempty(refs.spa)
        error('No temporal or spatial references for artifacts provided')
    else
        refs.tem = cell(1,length(refs.spa));
    end
else
    if isempty(refs.spa)
        refs.spa = cell(1,length(refs.tem));
    elseif length(refs.tem) ~= length(refs.spa)
        error('Number of temporal references needs to match number of spatial references')
    end
end
Nrefs = length(refs.tem);

if isfield(S,'ICs')
    try
        weights = S.weights;
        compvars = S.compvars;
    catch
        error('If passing pre-computed ICs, then must also pass weights and compvars')
    end
end

try d = S.d;                catch error('Must provide Channel x Time MEG data matrix'); end
try PCA_dim = S.PCA_dim;    catch PCA_dim = round(0.75*size(d,1)); end
try VarThr = S.VarThr;      catch VarThr = 0; end % could default to 100/PCA_dim?
try FiltPars = S.FiltPars;  catch FiltPars = []; end

try Rseed = S.Rseed;        catch Rseed = []; end  % If want reproducibility

% For permutation testing of temporal correlations (to turn off, pass S.Nperm = 0)
try PermPval = S.PermPval;  catch PermPval = .05/(2*PCA_dim);                   end % Kind of 2-tailed Bonferonni!
try Nperm = S.Nperm;        catch Nperm = round((4*sqrt(CorPval)/CorPval)^2);   end % http://www.epibiostat.ucsf.edu/biostat/sen/statgen/permutation.html#_how_many_permutations_do_we_need


%% Initialize
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
    if ~isempty(Rseed); warning('Random seed requested'); end
    [weights,sphere,compvars,bias,signs,lrates,ICs] = runica(d,'pca',PCA_dim,'extended',1,'maxsteps',800);
end


%% Filtering (if any) (and transposition for speed)
if length(FiltPars) == 3
    fprintf('Bandpass filtering from %d to %d Hz (warning - can fail)\n',FiltPars(1), FiltPars(2));
    for r=1:Nrefs
        refs.tem{r} = ft_preproc_bandpassfilter(refs.tem{r}, FiltPars(3), FiltPars(1:2),  [], 'but','twopass','reduce');
    end
    ICs  = ft_preproc_bandpassfilter(ICs,  FiltPars(3), FiltPars(1:2),  [], 'but','twopass','reduce')';
elseif length(FiltPars) == 2
    fprintf('Lowpass filtering to %d Hz\n',FiltPars(1));
    for r=1:Nrefs
        refs.tem{r} = ft_preproc_lowpassfilter(refs.tem{r}, FiltPars(2), FiltPars(1),  5, 'but','twopass','reduce');
    end
    ICs  = ft_preproc_lowpassfilter(ICs,  FiltPars(2), FiltPars(1),  5, 'but','twopass','reduce')';
else
    ICs  = ICs'; % Faster if pre-transpose once (check tic;toc)
end
% figure; for i=1:Nrefs; plot(refs.tem{r}'); title(i); pause; end
% figure; for i=1:PCA_dim; plot(ICs(:,i)); title(i); pause; end
% figure; for i=1:PCA_dim; plot(ICs(30000:40000,i)); title(i); pause; end

iweights  = pinv(weights);
remove = {};
temcor = zeros(Nrefs,PCA_dim); tempval = zeros(Nrefs,PCA_dim);
spacor = zeros(Nrefs,PCA_dim); spapval = zeros(Nrefs,PCA_dim);
reltemcor = zeros(Nrefs,PCA_dim); relspacor = zeros(Nrefs,PCA_dim);
temremove  = cell(1,Nrefs); sparemove  = cell(1,Nrefs);
for r = 1:Nrefs
    %parfor r = 1:Nrefs        % If want to parallelise (could parallelise Nperm loop below instead)
    
    %% Check temporal correlation with any reference channels
    if ~isempty(refs.tem{r})
        for k = 1:PCA_dim
            [temcor(r,k),tempval(r,k)] = corr(refs.tem{r}',ICs(:,k));
        end
        reltemcor(r,:) = zscore(temcor(r,:));
        
        if Nperm > 0
            permcor = zeros(1,PCA_dim);
            maxcor  = zeros(Nperm,1);
            
            ff = fft(refs.tem{r}',Nsamp);
            mf = abs(ff);
            wf = angle(ff);
            hf = floor((length(ff)-1)/2);
            rf = mf;
            
            for l = 1:Nperm
                btdata = zeros(size(ff));
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
            
            temremove{r} = find(abs(temcor(r,:)) > prctile(maxcor,100*(1-TemPval)));
            
        else
            temremove{r} = find(tempval(r,:) < TemAbsPval);
        end
        
        temremove{r} = [0 intersect(temremove{r},find(abs(reltemcor(r,:)) > TemRelZval))];  % 0 for no temp comps found, rather than none asked for (ie, passed empty)
    end
    
    %% Check spatial correlation with any reference channels
    if ~isempty(refs.spa{r})
        for k = 1:PCA_dim
            [spacor(r,k),spapval(r,k)] = corr(refs.spa{r}',iweights(:,k));
        end
        relspacor(r,:) = zscore(spacor(r,:));
        
        sparemove{r} = find(spapval(r,:) < SpaAbsPval);
        sparemove{r} = intersect(sparemove{r},find(abs(relspacor(r,:)) > SpaRelZval));
        
        if ~isempty(temremove{r}) % where a 0 entry is important for case where none found, rather than none asked for
            remove{r} = intersect(temremove{r},sparemove{r});
        else
            remove{r} = sparemove{r};
        end
    else
        if ~isempty(temremove{r}) 
            temremove{r}(find(temremove{r}==0))=[]; % get rid of any 0s from above
            remove{r} = temremove{r}; 
        end
    end
end
%figure,hist(temcor')
%figure,hist(reltemcor')
%figure,hist(spacor')
%figure,hist(log10(spapval+eps'))
%figure,hist(relspacor')

% checking/making each element of remove{} a row vector

rmat=[];
for ridx=1:length(remove),
   
    if iscolumn(remove{ridx}),
    
    rmat=[rmat,remove{ridx}']; 
    
    else
        
    rmat=[rmat,remove{ridx}];
    
    end
    
end

remove = unique(rmat);  % find unique ICs

%% Variance Thresholding
varexpl   = 100*compvars/sum(compvars);
varenough = find(varexpl > VarThr);
remove    = intersect(remove,varenough);  % plus sufficient Variance Explained (if required)

% figure; for i=remove; plot(ICs(30000:40000,i)); title(i); pause; end

if ~isempty(remove)
    finalics  = setdiff([1:PCA_dim],remove);
    TraMat    = iweights(:,finalics) * weights(finalics,:);
else
    TraMat    = eye(size(d,1));
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
% w = weights; iw = pinv(w);
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

                        