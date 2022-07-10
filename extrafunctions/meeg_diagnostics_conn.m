function fig = meeg_diagnostics_conn(data,diag,varargin)

if nargin >= 3, figtitle = varargin{1}; else, figtitle = 'Sample'; end
if nargin >= 4, savepath = varargin{2}; else, savepath = ''; end

%% config
if isfield(diag,'snapshotfwoi'), snapshotwoi = diag.snapshotfwoi; fieldofinterest = 'band';
elseif isfield(diag,'snapshotfwoiphase') && ~isempty(diag.snapshotfwoiphase), snapshotwoi = diag.snapshotfwoiphase; fieldofinterest = 'bandlow';
elseif isfield(diag,'snapshotfwoiamplitude') && ~isempty(diag.snapshotfwoiamplitude), snapshotwoi = diag.snapshotfwoiamplitude; fieldofinterest = 'bandhigh';
else, aas_log([],true,'no valid snapshot specification found');
end
if iscellstr(snapshotwoi)
    if numel(snapshotwoi) < size(data{1}.stat.stat,2)
        [~,indTWOI] = intersect(data{1}.(fieldofinterest),snapshotwoi);
    elseif numel(snapshotwoi) > size(data{1}.stat.stat,2)
        aas_log([],true,'snapshot(f)woi has a larger range than the output of statistics');
    end
end

[~, FT] = aas_cache_get([],'fieldtrip');
FT.load;
FT.addExternal('brewermap');
cmapbase = flipud(brewermap(128,'RdBu'));

[~, BNV] = aas_cache_get([],'bnv');
BNV.load;

inputfiles.surf = fullfile(BNV.toolPath,'Data','SurfTemplate','BrainMesh_Ch2_smoothed.nv');
inputfiles.map = '';

groupStat = data{1};

atlas = table(groupStat.elec.elecpos,ones(numel(groupStat.elec.label),2),strrep(groupStat.elec.label,' ','.'));
nROI = size(atlas,1);
fnNode = [savepath '_net.node'];
writetable(atlas,fnNode,'FileType','text','WriteVariableNames',false,'Delimiter',' ');
inputfiles.node = fnNode;

if isfield(groupStat,'freq'), freq = groupStat.freq;
elseif isfield(groupStat,'freqlow'), freq = groupStat.freqlow;
elseif isfield(groupStat,'freqhigh'), freq = groupStat.freqhigh;
else, aas_log([],true,'no fieldname found for frequency'); 
end

labelcmb = arrayfun(@(l) strjoin(groupStat.labelcmb(l,:),'-'),1:size(groupStat.labelcmb,1),'UniformOutput',false);
labeltickstep = diff(find(strcmp(groupStat.labelcmb(:,2),groupStat.labelcmb{1,2}))); labeltickstep = labeltickstep(1);
labeltick = 1:labeltickstep:size(groupStat.labelcmb,1);

%% Get data
mask = groupStat.stat.mask & ~isnan(groupStat.stat.stat);
stat = groupStat.stat.stat .* mask;
if ~any(stat,'all')
    delete(fnNode);
    BNV.unload;
    return; 
end
[cmap, cmaprange] = cmap_adjust(cmapbase,stat(mask));
if diff(cmaprange) == 0, cmaprange = sort(cmaprange .* [0.9 1.1]); end % add range in case of uniform values
if all(stat(stat~=0)>0), sfx = 'p';
elseif all(stat(stat~=0)<0), sfx = 'n';
else, sfx = 'pn';
end

%% Matrix
% image
h = figure;
set(h,'Position',[0,0,1920,1080]);
set(h,'PaperPositionMode','auto');
ax = axes(h);
imAlpha=ones(size(stat));
imAlpha(~stat)=0;
imagesc(stat,'AlphaData',imAlpha);
colormap(cmap); caxis(cmaprange); colorbar;
set(ax, 'Color', [0.75 0.75 0.75]);
for a = ax
    set(a,'YTick',labeltick,'YGrid','on'); set(a,'YTickLabel',{});
    set(a,'XTick',1:size(stat,2)); 
    if iscellstr(snapshotwoi)
        set(a,'XTickLabel',snapshotwoi);
    else
        set(a,'XTickLabel',round(freq(get(a,'XTick'))));
    end
end
set(ax(1),'YTickLabel',labelcmb(labeltick));
set(h,'Name',figtitle);
if ~isempty(savepath)
    print(h,'-noui',[savepath '_matrix.jpg'],'-djpeg','-r300');
    close(h);
end

% txt
fnTxt = [savepath '_matrix.txt'];
idx = find(mask(:));
if ~isempty(idx)
    indx = cell(1,2);
    [indx{:}]=ind2sub(size(mask),idx);
    [indx{1},sind] = sort(indx{1});
    indx{2} = indx{2}(sind);
    
    fid = fopen(fnTxt,'w');
    prevlc = 0;
    for l=1:numel(idx)
        if l == 1
            fprintf(fid,'%s -> %s: %1.2f', ...
                groupStat.labelcmb{indx{1}(l),1}, groupStat.labelcmb{indx{1}(l),2}, freq(indx{2}(l)));
        elseif indx{1}(l) ~= prevlc
            fprintf(fid,'\n%s -> %s: %1.2f', ...
                groupStat.labelcmb{indx{1}(l),1}, groupStat.labelcmb{indx{1}(l),2}, freq(indx{2}(l)));
        else
            fprintf(fid,', %1.2f', freq(indx{2}(l)));
        end
        prevlc = indx{1}(l);
    end
    fclose(fid);
end

%% Network (for each FOI)
if numel(sfx) > 1, aas_log([],true,'mixed stat - NYI'); end

for f = 1:size(snapshotwoi,1)
    mat = zeros(nROI,nROI);
    for roi = unique(groupStat.labelcmb(:,1),'stable')'
        roiind = strcmp(atlas.Var3,strrep(roi{1},' ','.'));
        lcoi = strcmp(groupStat.labelcmb(:,1),roi{1});
        if iscellstr(snapshotwoi)
            boiind = [f f];
            fnEdge = sprintf('%s_net_%s_%c.edge',savepath,snapshotwoi{f},sfx);
        else
            boiind = [find(freq >= snapshotwoi(f,1),1,'first') find(freq <= snapshotwoi(f,2),1,'last')];
            fnEdge = sprintf('%s_net_%d-%d_%c.edge',savepath,snapshotwoi(f,:),sfx);
        end
        meas = mean(stat(lcoi,boiind(1):boiind(2)),2);
        if nROI > numel(meas)
            mat(:,roiind) = [meas(1:find(roiind)-1); 0; meas(find(roiind):end)]; % auto-connectivity
        else
            mat(:,roiind) = meas;
        end
        % mat(mat(:,roiind)<cmaprange(1)/2,roiind) = 0; % at least half of the band
    end
    if ~any(mat,'all'), continue; end
    
    dlmwrite(fnEdge,mat,'\t');    
    inputfiles.edge = fnEdge;
    fig = BrainNet(inputfiles,BNV.BNVSettings);
    set(fig,'Name',figtitle);
    if ~isempty(savepath)
        print(fig,'-noui',spm_file(fnEdge,'ext','jpg'),'-djpeg','-r300');
        close(fig);
    end
end
if isempty(savepath)
    delete(fnNode);
    delete(fnEdge);
end
BNV.unload;
end

function [cmap, cmaprange] = cmap_adjust(cmapbase,stat)
minval = prctile(stat,1,'all');
maxval = prctile(stat,99,'all');
% colormaps
if (minval < 0) && (maxval > 0)
    r = maxval/-minval;
    if r > 1, cmap = [cmapcold(cmapbase,round(64/r)); cmaphot(cmapbase,64)];
    else, cmap = [cmapcold(cmapbase,64); cmaphot(cmapbase,round(r*64))];
    end
elseif minval < 0, cmap = cmapcold(cmapbase,64);
else
    cmap = cmaphot(cmapbase,64);
end
cmaprange = [minval, maxval];
end

function cmap = cmaphot(cmapbase,n)
cmap = cmapbase(65:(end-(64-n)),:);
end

function cmap = cmapcold(cmapbase,n)
cmap = cmapbase((64-n+1):64,:);
end