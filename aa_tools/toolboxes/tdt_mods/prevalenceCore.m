function [results, params] = prevalenceCore(a, P2, alpha)

% permutation-based prevalence inference using the minimum statistic, core
% implementation of the method proposed by Allefeld, Goergen and Haynes (2016)
%
% [results, params] = prevalenceCore(a, P2 = 1e6, alpha = 0.05)
%
% a:            three-dimensional array of test statistic values
%               (voxels x subjects x first-level permutations)
%               a(:, :, 1) must contain actual values
% P2:           number of second-level permutations to generate
% alpha:        significance level
% results:      per-voxel analysis results
%   .puGN         uncorrected p-values for global null hypothesis         (Eq. 24)
%   .pcGN         corrected p-values for global null hypothesis           (Eq. 26)
%   .puMN         uncorrected p-values for majority null hypothesis       (Eq. 19)
%   .pcMN         corrected p-values for majority null hypothesis         (Eq. 21)
%   .gamma0c      corrected prevalence lower bounds                       (Eq. 23)
%   .aTypical     median values of test statistic where pcMN <= alpha     (Fig. 4b)
% params:        analysis parameters and properties
%   .V            number of voxels
%   .N            number of subjects
%   .P1           number of first-level permutations
%   .P2           number of second-level permutations actually generated
%   .alpha        significance level
%   .pcMNMin      smallest possible corrected p-value for majority H0
%   .gamma0cMax   largest possible corrected prevalence lower bound       (Eq. 27)
%
% The 'majority null hypothesis' referenced here is a special case of the
% prevalence null hypothesis (Eq. 17), where the critical value is gamma0 =
% 0.5. It describes the situation where there is no effect in the majority
% of subjects in the population. Rejecting it allows to infer that there is
% an effect in the majority of subjects in the population. aTypical is only
% defined where the (spatially extended) majority null hypothesis can be
% rejected. Compare Fig. 4b and 'What does it mean for an effect to be
% "typical" in the population?' in the Discussion of Allefeld, Goergen and
% Haynes (2016).
%
% The function opens a figure window that shows results based on the
% permutations generated so far and is regularly updated. If this window is
% closed, the computation is stopped (not aborted) and results are returned
% based on the permutations so far.
%
% See also prevalence.
%
%
% This file is part of v1.0.1 of prevalence-permutation, see
% https://github.com/allefeld/prevalence-permutation/releases
%
% Copyright (C) 2016 Carsten Allefeld
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.

% default parameter values
if (nargin < 2) || isempty(P2)
    P2 = 1e6;
end
if (nargin < 3) || isempty(alpha)
    alpha = 0.05;
end

% init
[V, N, P1] = size(a); % V: voxels, N: subjects, P1: first-level permutations
fprintf('generating %d of %d second-level permutations\n', P2, P1 ^ N)
if P2 > P1 ^ N
    error('Monte Carlo implementation is inadequate!')
    % implement full enumeration of permutations? (Issue #1)
end
fprintf('the computation can be stopped at any time by closing the output window\n\n')

% prepare plot window
fh = figure('Name', 'permutation-based prevalence inference using the minimum statistic');
text(0.5, 0.5, {'please wait for results', '', ...
    'close window to stop computation at any time'}, 'HorizontalAlignment', 'center')
axis off
drawnow

% generate second-level permutations
uRank = zeros(V, 1);
cRank = zeros(V, 1);
nPermsReport = 1;
tic
for j = 1 : P2
    % select first-level permutations
    if j == 1
        % select neutral permutations
        sp = ones(N, 1);
    else
        % select permutations randomly (Monte Carlo)
        sp = randi(P1, N, 1);
    end
    % indices of permutation values for each subject
    ind = sub2ind([N, P1], (1 : N)', sp);
    
    % test statistic: minimum across subjects
    m = min(a(:, ind), [], 2);
    % store result of neutral permutation (actual value) for each voxel
    if j == 1
        m1 = m;
    end
    
    % compare actual value with permutation value for each voxel separately,
    % determines uncorrected p-values for global null hypothesis (see below)
    uRank = uRank + (m >= m1);          % part of Eq. 24
    % compare actual value at each voxel with maximum across voxels,
    % determines corrected p-values for global null hypothesis (see below)
    cRank = cRank + (max(m) >= m1);     % Eq. 25 & part of Eq. 26
    
    % calibrate reporting interval to be between 2 and 5 seconds
    if (nPermsReport == 1) && (toc >= 2)
        nPermsReport = 10 ^ floor(log10(j)) * [1 2 5 10];
        nPermsReport = min(nPermsReport(nPermsReport >= j));
        nPermsReport = max(nPermsReport, 2);
    end
    
    % at intervals
    if ((nPermsReport > 1) && (mod(j, nPermsReport) == 0)) || (j == P2)
        drawnow
        stop = ~ishandle(fh);
        
        % compute results,
        % based on permutations performed so far: j plays the role of P2
        % uncorrected p-values for global null hypothesis
        puGN = uRank / j;                               % part of Eq. 24
        % corrected p-values for global null hypothesis
        pcGN = cRank / j;                               % part of Eq. 26
        % significant voxels for global null hypothesis
        sigGN = (pcGN <= alpha);
        % * Step 5a: compute p-values for given prevalence bound
        % (here specifically gamma0 = 0.5, i.e the majority null hypothesis)
        % uncorrected p-values for majority null hypothesis
        puMN = ((1 - 0.5) * puGN .^ (1/N) + 0.5) .^ N;  % Eq. 19
        % corrected p-values for majority null hypothesis
        pcMN = pcGN + (1 - pcGN) .* puMN;               % Eq. 21
        % significant voxels for majority null hypothesis
        sigMN = (pcMN <= alpha);
        % lower bound on corrected p-values for majority null hypothesis
        puMNMin = ((1 - 0.5) * 1/j .^ (1/N) + 0.5) .^ N;
        pcMNMin = 1/j + (1 - 1/j) .* puMNMin;
        % * Step 5b: compute lower bounds for prevalence
        % lower bounds for prevalence
        alphac = (alpha - pcGN) ./ (1 - pcGN);          % Eq. 22
        gamma0c = (alphac .^ (1/N) - puGN .^ (1/N)) ./ (1 - puGN .^ (1/N)); % Eq. 23
        gamma0c(alphac < puGN) = nan;                   % undefined
        % upper bound for lower bounds
        alphacMax = (alpha - 1/j) / (1 - 1/j);          % Eq. 27
        gamma0cMax = (alphacMax .^ (1/N) - 1/j .^ (1/N)) ./ (1 - 1/j .^ (1/N)); % Eq. 27
        
        % print summary
        fprintf('  %d of %d permutations\n', j, P2)
        fprintf('    minimal rank\n')
        fprintf('      uncorrected:  %d,  reached at %d voxels\n', ...
            min(uRank), sum(uRank == min(uRank)))
        fprintf('      corrected:    %d,  reached at %d voxels\n', ...
            min(cRank), sum(cRank == min(cRank)))
        fprintf('    minimal p-value for global null hypothesis\n')
        fprintf('      uncorrected:  %g\n', ...
            min(puGN))
        fprintf('      corrected:    %g\n', ...
            min(pcGN))
        fprintf('    minimal p-value for majority null hypothesis\n')
        fprintf('      uncorrected:  %g\n', ...
            min(puMN))
        fprintf('      corrected:    %g\n', ...
            min(pcMN))
        fprintf('    number of voxels (of %d) at which (corrected)\n', V)
        fprintf('      global null hypothesis is rejected:    %d\n', ...
            sum(sigGN))
        fprintf('      majority null hypothesis is rejected:  %d\n', ...
            sum(sigMN))
        fprintf('      prevalence bound is defined:           %d\n', ...
            sum(~isnan(gamma0c)))
        fprintf('    largest corrected prevalence bound:  %g\n', max(gamma0c))
        fprintf('\n')
        
        % graphical display
        if stop
            fh = figure('Name', 'permutation-based prevalence inference using the minimum statistic');
        else
            % make figure current without getting in the way
            set(groot, 'CurrentFigure', fh)
            clf
        end
        % prevalence bounds
        if V > 0
            subplot(2, 1, 1)
            plot([0.5, V + 0.5], 0.5 * [1 1], 'k')
            hold all
            plot([0.5, V + 0.5], gamma0cMax * [1 1], 'r')
            plot(gamma0c, 'b.')
            axis([0.5, V + 0.5, 0, 1])
            title('prevalence lower bounds')
            ylabel('\gamma_0')
            if V > 200
                set(gca, 'XTick', [])
            else
                set(gca, 'XTick', 1 : V)
                if V > 20
                    set(gca, 'XTickLabel', [])
                end
            end
            if j < P2
                text(0.5, 1, sprintf(' %.0f %%', j / P2 * 100), ...
                    'Color', [0.8 0.8 0.8], 'HorizontalAlignment', 'left', ...
                    'VerticalAlignment', 'top')
                text(V + 0.5, 1, sprintf('%.1f / %.1f min ', ...
                    toc / 60, toc / 60 * P2 / j), ...
                    'Color', [0.8 0.8 0.8], 'HorizontalAlignment', 'right', ...
                    'VerticalAlignment', 'top')
            end
            % majority null hypothesis p-values
            subplot(2, 1, 2)
            semilogy([0.5, V + 0.5], alpha * [1 1], 'k')
            hold all
            plot([0.5, V + 0.5], pcMNMin * [1 1], 'r')
            plot(pcMN, '.b')
            xlim([0.5, V + 0.5])
            title('p-values majority null hypothesis')
            xlabel('voxels')
            ylabel('p^*_N(m | \gamma \leq 0.5)')
            if V > 200
                set(gca, 'XTick', [])
            else
                set(gca, 'XTick', 1 : V)
                if V > 20
                    set(gca, 'XTickLabel', [])
                end
            end
            drawnow
        end
        
        if stop
            s = warning('query', 'backtrace');
            warning off backtrace
            warning('computation has been stopped')
            warning(s)
            P2 = j;
            break
        end
    end
    
end

% where majority null hypothesis can be rejected, typical value of test statistic
aTypical = nan(V, 1);
aTypical(sigMN) = median(a(sigMN, :, 1), 2);

% collect return values
params = struct;
params.V = V;
params.N = N;
params.P1 = P1;
params.P2 = P2;
params.alpha = alpha;
params.pcMNMin = pcMNMin;
params.gamma0cMax = gamma0cMax;
results = struct;
results.puGN = puGN;
results.pcGN = pcGN;
results.puMN = puMN;
results.pcMN = pcMN;
results.gamma0c = gamma0c;
results.aTypical = aTypical;
