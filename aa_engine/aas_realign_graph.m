function aas_realign_graph(mv)
% mv can be an array (strvcat or cell) of file names of movement parameters...
% ... or a cell array of matrices containing the parameters...
if ~iscell(mv)
    mv = strvcat2cell(mv);
end

movePars = [];
for s = 1:length(mv)
    if isstr(mv{s})
        % Typically, movePars will be a text file...
        tmp = load(deblank(mv{s}));
        movePars = [movePars; tmp];
    else
        movePars = mv{s};
    end
end

% Movements
movePars=movePars-repmat(movePars(1,:),[size(movePars,1) 1]);

movePars(:,4:6)=movePars(:,4:6)*180/pi; % convert to degrees!

mvmean=mean(movePars);
mvmax=max(movePars);
mvstd=std(movePars);

try close(2); catch; end
figure(2)
scrsz = get(0,'ScreenSize');
set(2, 'Position',[1 1 scrsz(3) scrsz(4)])

DmovePars = movePars(:,1:3);
RmovePars = movePars(:,4:6);

subplot(2,1,1)
plot(DmovePars)
xlim([0, size(DmovePars,1)])
ylim([-5 5])
title('Displacement (mm) [x: blue; y: green; z:red]')

subplot(2,1,2)
plot(RmovePars)
xlim([0, size(RmovePars,1)])
ylim([-3 3])
title('Rotation (deg) [r: blue; p: green; j:red]')
    
    
       