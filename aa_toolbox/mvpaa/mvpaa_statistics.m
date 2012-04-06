% MVPAa_STATISTICS Sort how to do MVPA statistics depending on data
% C = correlation data
% EP.triangulation = how are the block*subblock comparisons structured

function [Stats meanSimil] = mvpaa_statistics(aap, Simil)

% Rename settings to keep easier track...
EP = aap.tasklist.currenttask.settings;

%C => [EP.blocks*EP.sessions, EP.blocks*EP.sessions, EP.conditions, EP.conditions]

ind=0;
if EP.triangulation == 1
    uSimil = zeros(EP.conditions, ...
        EP.conditions, ...
        ((EP.blocks ...
        * EP.sessions)^2 ...
        - EP.blocks ...
        * EP.sessions) / 2);
    % Triangulate, such that comparisons across EP.blocks occur only once
    for i = 1:(EP.blocks*EP.sessions)
        for j = (i+1):(EP.blocks*EP.sessions)
            ind = ind + 1;
            uSimil(:,:,ind) = Simil(i,j,:,:);
        end
    end
elseif EP.triangulation == 0
    uSimil = zeros(EP.conditions, ...
        EP.conditions, ...
        (EP.blocks ...
        * EP.sessions)^2 ...
        - EP.blocks ...
        * EP.sessions);
    % Don't triangulate (useful for sparse data... and unusual designs)
    for i = 1:(EP.blocks*EP.sessions)
        for j = 1:(EP.blocks*EP.sessions) % Only one triangle of matrix carries information...
            if i == j
                % Ignore within block-subblock comparisons
                continue
            end
            ind = ind + 1;
            uSimil(:,:,ind) = Simil(i,j,:,:);
        end
    end
elseif EP.triangulation == 2
    % One single big matrix!
    % Already sorted out in mvpaa_correlation script.
    uSimil = zeros(EP.conditions ...
        * EP.blocks ...
        * EP.sessions, ...
        EP.conditions ...
        * EP.blocks ...
        * EP.sessions, 1);
    uSimil(:,:,1) = Simil;    
elseif EP.triangulation == 3
    % Triangulation but only consider within block comparisons...
    uSimil = zeros(EP.conditions, ...
        EP.conditions, ...
        EP.sessions ...
        * ((EP.blocks)^2 ...
        - EP.blocks) / 2);
    % Triangulate, such that comparisons across EP.blocks within each block occur only once...    
    for b = EP.sessions
        for i = 1:(EP.blocks)
            for j = (i+1):(EP.blocks)
                ind = ind + 1;
                uSimil(:,:,ind) = Simil((b-1)*EP.blocks + i, (b-1)*EP.blocks + j,:,:);
            end
        end
    end
elseif EP.triangulation == 4
    % This averages over all subblock/block comparisons (excluding the
    % equivalent ones, but does not exclue the leading diagonal, instead
    % naning it when i and j are equal. Finally, it nanmeans the matrix.
    % Your contrasts should probably only consider the tril or triu of the
    % contrast matrix (@@@ TO BE IMPLEMENTED @@@);
    uSimil = zeros(EP.conditions, ...
        EP.conditions, ...
        (EP.blocks * ...
        EP.sessions)^2);
    for i = 1:(EP.blocks*EP.sessions)
        for j = i:(EP.blocks*EP.sessions)
            ind = ind + 1;
            if i == j
                % Ignore within block-subblock comparisons
                tmp = squeeze(Simil(i,j,:,:));
                tmp(logical(eye(size(tmp)))) = NaN;
                uSimil(:,:,ind) = tmp;
            else
                uSimil(:,:,ind) = Simil(i,j,:,:);
            end
        end
    end
    uSimil = nanmean(uSimil,3);
end

% Useful for ROI analysis visualization...
meanSimil = mean(uSimil,3);

Stats = mvpaa_stats1st(aap, uSimil);