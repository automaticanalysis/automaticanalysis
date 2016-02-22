function [moveRegs Rnames] = aas_movPars(aap, subjI, moveMat, volterraMovement)
% Select what type of movement parameters to use, and create them!
%
% [moveRegs, Rnames] = aas_movPars(aap, p, moveMat, volterraMovement)
% moveMat selects which orders (squared, cubed, etc.)...
% and derivatives (gradient, gradient(gradient), etc.) we use!
% Good default: moves, moves.^2, gradient(moves) and gradient(moves.^2)
% but no "spin history", which is (fairly highly) correlated with grad.
%
% A longer, explicit example...:
% [moves    grad(moves)    grad(grad(moves))    ...  spinhist(moves)
%  moves^2  grad(moves^2)  grad(grad(moves^2))  ...  spinhist(moves^2)
%  moves^3  grad(moves^3)  grad(grad(moves^3))  ...  spinhist(moves^3)
%  ...      ...            ...                  ...  ...]
%
% volterraMovement includes a Volterra expansion of the movemenet
% parameters in addition to other factors. If not specified, not done.


if nargin<4 || isempty(volterraMovement)
    volterraMovement = 0;
end

if nargin<3 || isempty(moveMat)
    moveMat = [1 1 0; ...
        1 1 0];    
end


% Order and derivative start at 0
maxO = size(moveMat, 1) - 1; % Maximal order you want
maxD = size(moveMat,2) - 2; % Maximal derivative you want
%NOTE: last derivative [maxD+2] will be the spin history...

% We can now have missing sessions per subject, so we're going to use only
% the sessions that are common to this subject and selected_sessions
[numSess, sessInds] = aas_getN_bydomain(aap, 'session', subjI);
subjSessionI = intersect(sessInds, aap.acq_details.selected_sessions);
numSess = numel(subjSessionI);

% Get the movement parameters for each session separately...
moveRegs = cell(numSess,1);

for s = 1:numSess
        
    % Linear movement parameters text file...
    Mfn = aas_getfiles_bystream(aap, subjI, subjSessionI(s), 'realignment_parameter');
    for f = 1:size(Mfn,1)
        fn = deblank(Mfn(f,:));
        if strcmp(fn(end-2:end), 'txt')
            Mfn = fn;
            break
        end
    end
    moves = load(Mfn);
    
    % If more than 6 columns, just use first 6 (e.g. for HCP data)
    if size(moves,2) > 6
        moves = moves(:,1:6);
    end
    
    movesN = cell(size(moveMat));
    
    for o = 1:(maxO+1)
        movesN{o,1} = moves.^o;
        % Movement derivatives
        for d = 2:(maxD+1)
            [junk, movesN{o,d}] = gradient(movesN{o,d-1});
        end
        % Spin history is difference between current and old...
        % It is also awfully similar to the gradient...
        movesN{o, size(moveMat,2)} = [zeros(1,6); movesN{o,1}(2:end, :) - movesN{o,1}(1:end-1,:)];
    end
    
    for o = 1:(maxO+1)
        for d = 1:(maxD+2)
            if moveMat(o,d)
                moveRegs{s} = [moveRegs{s}, movesN{o,d}];
            end
        end
    end
    
    
    % Volterra expansion, if requested
    if volterraMovement
        M = spm_detrend(moves);
        U=[];
        for c=1:6
            U(c).u = M(:,c);
            U(c).name{1}='c';
        end
        aM = spm_Volterra(U,[1 0 0; 1 -1 0; 0 1 -1]', 2);
        
        moveRegs{s} = [moveRegs{s} aM];                
    end
    
    % Show an image of correlated timecourses... DEBUG
    %corrTCs(moveRegs, Rnames)   
end



% Basic names
bnames = {'x' 'y' 'z' 'r' 'p' 'j'};

Rnames = {};

% Create final names for movement regressors
for o = 1:(maxO+1)
    for d = 1:(maxD+2)
        if moveMat(o,d)
            Rnames = [Rnames bnames{:}];
            % Put order in names
            if o > 1
                for n = 0:5
                    Rnames{end-n} = [Rnames{end-n} '^' num2str(o)];
                end
            end
            if d > 1 && d < size(moveMat,2)
                for n = 0:5
                    for t = 2:d
                        Rnames{end-n} = ['g(' Rnames{end-n} ')'];
                    end
                end
            elseif d == size(moveMat,2)
                for n = 0:5
                    Rnames{end-n} = ['sh(' Rnames{end-n} ')'];
                end
            end
        end
    end
end

if volterraMovement
    for v = 1:size(aM, 2)
        Rnames{end+1} = sprintf('Volterra movement %d', v);
    end
end % adding names for Volterra


