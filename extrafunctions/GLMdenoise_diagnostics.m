% Inputs from GLMdenoise
% resultsGD = resultsGD from denoised regressors
% dataGD = dataOrig after GLMdenoising
% resultsOrig = resultsGD without PC denoising
% dataOrig = original dataOrig

function [h, figName] = GLMdenoise_diagnostics(resultsGD, dataGD, resultsOrig, dataOrig)
names = {'Original', 'Denoised'};

% Check dimensionality of the dataOrig:
dims = length(resultsGD.inputs.datasize{1}) - 1;

h = [];
figName = {};

%% Inspect outputs

% The outputs of GLMdenoisedata are contained in
% the variables 'resultsGD' and 'dataGD'.
% Here we do some basic inspections of the outputs.

% Select a voxel to inspect.  This is done by finding voxels that have
% cross-validated R^2 values between 0% and 5% under the initial model (no PCs),
% and then selecting the voxel that shows the largest improvement when
% using the final model.

if dims == 3;
    ix = find(resultsGD.pcR2(:,:,:,1) > 0 & resultsGD.pcR2(:,:,:,1) < 5);
    improvement = resultsGD.pcR2(:,:,:,1+resultsGD.pcnum) - resultsGD.pcR2(:,:,:,1);
    [mm,ii] = max(improvement(ix));
    ix = ix(ii);
    [xx,yy,zz] = ind2sub(resultsGD.inputs.datasize{1}(1:3),ix);
elseif dims == 1;
    ix = find(resultsGD.pcR2(:,1) > 0 & resultsGD.pcR2(:,1) < 5);
    improvement = resultsGD.pcR2(:,1+resultsGD.pcnum) - resultsGD.pcR2(:,1);
    [mm,ii] = max(improvement(ix));
    ix = ix(ii);
end

% Inspect amplitude estimates for the example voxel.  The first row shows resultsGD
% obtained from the GLMdenoisedata call that does not involve global noise regressors;
% the second row shows resultsGD obtained from the regular GLMdenoisedata call.
% The left panels show amplitude estimates from individual bootstraps, whereas
% the right panels show the median and 68% interval across bootstraps (this
% corresponds to the final model estimate and the estimated error on the
% estimate, respectively).
h(end + 1) = figure;
figName{end + 1} = 'Regressors';
set(gcf,'Units','points','Position',[100 100 700 500]);
for p=1:2
    if p==1
        if dims == 3
            ampboots = squeeze(resultsOrig.models{2}(xx,yy,zz,:,:));  % conditions x boots
            amp = flatten(resultsOrig.modelmd{2}(xx,yy,zz,:));        % 1 x conditions
            ampse = flatten(resultsOrig.modelse{2}(xx,yy,zz,:));      % 1 x conditions
        elseif dims == 1
            ampboots = squeeze(resultsOrig.models{2}(ix,:,:));  % conditions x boots
            amp = flatten(resultsOrig.modelmd{2}(ix,:));        % 1 x conditions
            ampse = flatten(resultsOrig.modelse{2}(ix,:));      % 1 x conditions
        end
    else
        if dims == 3
            ampboots = squeeze(resultsGD.models{2}(xx,yy,zz,:,:));     % conditions x boots
            amp = flatten(resultsGD.modelmd{2}(xx,yy,zz,:));           % 1 x conditions
            ampse = flatten(resultsGD.modelse{2}(xx,yy,zz,:));         % 1 x conditions
        elseif dims == 1
            ampboots = squeeze(resultsGD.models{2}(ix,:,:));     % conditions x boots
            amp = flatten(resultsGD.modelmd{2}(ix,:));           % 1 x conditions
            ampse = flatten(resultsGD.modelse{2}(ix,:));         % 1 x conditions
        end
    end
    n = length(amp);
    subplot(2,2,(p-1)*2 + 1); hold on;
    plot(ampboots);
    straightline(0,'h','k-');
    xlabel('Condition number');
    ylabel('BOLD signal (% change)');
    title(sprintf('%s: Amplitude estim. (ind. bootstraps)', names{p}));
    subplot(2,2,(p-1)*2 + 2); hold on;
    bar(1:length(amp),amp,1);
    errorbar2(1:length(amp),amp,ampse,'v','r-');
    if p==1
        ax = axis; axis([0 n+1 ax(3:4)]); ax = axis;
    end
    xlabel('Condition number');
    ylabel('BOLD signal (% change)');
    title(sprintf('%s: Amplitude estim. (median, 68% interval)', names{p}));
end
for p=1:2
    subplot(2,2,(p-1)*2 + 1); axis(ax);
    subplot(2,2,(p-1)*2 + 2); axis(ax);
end
%%

% Compare SNR before and after the use of global noise regressors.
% To focus on voxels that related to the experiment, we select voxels with
% cross-validated R^2 values that are greater than 0% under the initial model.
% To ensure that the SNR values reflect only changes in the noise level,
% we ignore the SNR computed in each individual GLMdenoisedata call and
% re-compute SNR, holding the numerator (the signal) constant across the
% two calls.  (Note: GLMdenoisedata automatically writes out a figure that
% shows a before-and-after SNR comparison (SNRcomparebeforeandafter.png).
% What is shown in this script does the comparison manually just for sake
% of example.)

h(end + 1) = figure;
figName{end + 1} = 'SNRdiff';

if length(resultsGD.inputs.datasize{1}) == 4;
    ok = resultsGD.pcR2(:,:,:,1) > 0;
else
    ok = resultsGD.pcR2(:,1) > 0;
end
signal = mean([resultsGD.signal(ok) resultsOrig.signal(ok)],2);
snrOrig = signal ./ resultsOrig.noise(ok);
snrGD = signal ./ resultsGD.noise(ok);

hold on;
scatter(snrOrig, snrGD,'r.');
ax = axis;
mx = max(ax(3:4));
axissquarify;
axis([0 mx 0 mx]);
xlabel('SNR before');
ylabel('SNR after');

[p,ht,stats] = signrank(snrGD, snrOrig);
title(sprintf('Signed rank: Z: %0.2f, p: %0.6f, med-diff: %0.2f', stats.zval, p, median(snrGD-snrOrig)))
%%

% An alternative to using the GLM estimates provided by GLMdenoisedata is
% to use 'dataGD', which contains the original time-series dataOrig but
% with the component of the dataOrig that is estimated to be due to the global
% noise regressors subtracted off.  Here we inspect the denoised dataOrig for
% the same example voxel examined earlier.  Note that the dataOrig components
% that are present in 'dataGD' can be customized (see opt.denoisespec
% in GLMdenoisedata.m).  The default parameters leave in the estimated baseline
% signal drift, which explains the drift in the plotted time-series.
h(end + 1) = figure;
figName{end + 1} = 'Timecourse';
hold on;

set(gcf,'Units','points','Position',[100 100 700 250]);
if dims == 3
    data1 = flatten(dataOrig{1}(xx,yy,zz,:));
    data2 = flatten(dataGD{1}(xx,yy,zz,:));
elseif dims == 1
    data1 = flatten(dataOrig{1}(ix,:));
    data2 = flatten(dataGD{1}(ix,:));
end
n = length(data1);
h1 = plot(data1,'ro-');
h2 = plot(data2,'bo-');
ax = axis; axis([0 n+1 ax(3:4)]);
legend([h1 h2],{'Original' 'Denoised'});
xlabel('Time point');
ylabel('MR signal');
%%