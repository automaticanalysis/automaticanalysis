function [moveRegs Rnames] = aas_movPars(aap,p, moveMat)
% Select what type of movement parameters to use, and create them!

if isempty(moveMat)
    % moveMat selects which orders (squared, cubed, etc.)...
    % and derivatives (gradient, gradient(gradient), etc.) we use!
    % Good default: moves, moves.^2, gradient(moves) and gradient(moves.^2)
    % but no "spin history", which is (fairly highly) correlated with grad.
    moveMat = [1 1 0; ...
        1 1 0];
    
    % A longer, explicit example...:
    % [moves    grad(moves)    grad(grad(moves))    ...  spinhist(moves)
    %  moves^2  grad(moves^2)  grad(grad(moves^2))  ...  spinhist(moves^2)
    %  moves^3  grad(moves^3)  grad(grad(moves^3))  ...  spinhist(moves^3)
    %  ...      ...            ...                  ...  ...]
end

% Basic names
bnames = {'x' 'y' 'z' 'r' 'p' 'j'};

% Order and derivative start at 0
maxO = size(moveMat, 1) - 1; % Maximal order you want
maxD = size(moveMat,2) - 2; % Maximal derivative you want
%NOTE: last derivative [maxD+2] will be the spin history...

Rnames = {};

% Create final names
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

% Get the movement parameters for each session separately...
moveRegs = cell(size(aap.acq_details.sessions));

for s = aap.acq_details.selected_sessions
        
    % Linear movement parameters text file...
    Mfn = aas_getfiles_bystream(aap,p,s,'realignment_parameter');
    for f = 1:size(Mfn,1)
        fn = deblank(Mfn(f,:));
        if strcmp(fn(end-2:end), 'txt')
            Mfn = fn;
            break
        end
    end
    moves = load(Mfn);
    
    movesN = cell(size(moveMat));
    
    for o = 1:(maxO+1)
        movesN{o,1} = moves.^o;
        % Movement derivatives
        for d = 2:(maxD+1)
            [~, movesN{o,d}] = gradient(movesN{o,d-1});
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
    
    % Show an image of correlated timecourses... DEBUG
    %corrTCs(moveRegs, Rnames)   
end