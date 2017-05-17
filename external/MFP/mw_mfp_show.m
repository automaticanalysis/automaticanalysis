function mw_mfp_show(p)

% Init
if ~nargin, p = pwd; end
fn=dir(fullfile(p,'rp_*.txt'));
pr = load(fullfile(p,fn(1).name));
load(fullfile(p,'mw_motion.mat'));
fn=dir(fullfile(p,'mw_mfp_*.txt'));
curr_tc = load(fullfile(p,fn(1).name));
shifted = logical(str2double(fn(1).name(strfind(fn(1).name,'sh')+2)));
squared = logical(str2double(fn(1).name(strfind(fn(1).name,'sq')+2)));
mysteps = 40;
pres = 150;

% get handle & average voxel size
fn=dir(fullfile(p,'mw_motmask.*'));
V = spm_vol(fullfile(p,fn(1).name));
vs = (abs(V.mat(1,1))+abs(V.mat(2,2))+abs(V.mat(3,3))) / 3;

% create graphics
gcf1 = figure;
set(0,'Units','pixels');
scnsize = get(0,'ScreenSize');
nrows = 5; % do_mfp, do_fancy
pos1 = [round(scnsize(3)/40), round(scnsize(4)/40), round(scnsize(3)/10*6), round(scnsize(4)*nrows/5.75)];
set(gcf1, 'Position', pos1, 'Color', 'w', 'Renderer', 'zbuffer');
set(gcf1, 'Name', ['Motion parameters from ' p], 'Visible', 'on');

% optimize scaling
mot = abs(pr(:,1:3));
deg = abs((pr(:,4:6)*180/pi));
scaleme = ([mot deg td' sts']);
scaleme = [ceil(max(abs(scaleme(:)))) * -1 ceil(max(abs(scaleme(:))))];


% first, the usual, well-known suspects: translations...
subplot(nrows,3,[1 2 3],'replace');
plot(pr(:,1:3));
grid on;
xlim([0 size(pr,1)]);
ylim(scaleme);
title('Realignment parameters: shifts [top, in mm] and rotations [middle, in dg]');


% ... and rotations
subplot(nrows,3,[4 5 6],'replace');
plot(pr(:,4:6)*180/pi);
grid on;
xlim([0 size(pr,1)]);
ylim(scaleme);
title(['(data from ' p ')'], 'interpreter', 'none');


% ... and now, total replacement at d_avg: total & scan-to-scan
sp3 = subplot(nrows,3,[7 8 9],'replace');
plot(td);
hold on;
bar(sts);
grid on;
xlim([0 size(pr,1)]);
ylim(scaleme);
title(['Total displacement: overall (line) and scan-to-scan (bars) at d_a_v_g (here: ' sprintf('%02.1f', davg) ' mm)']);


% use the old rule of thumb to show "one voxel size"
if max(scaleme) > vs
    
    props=get(sp3);
    lh1 = line(props.XLim,[-vs -vs],'color','r', 'linestyle', '--');
    lh2 = line(props.XLim,[vs   vs],'color','r', 'linestyle', '--');
    drawnow;
end;

% annotate
xm = find(mot == max(mot(:)));
xd = find(deg == max(deg(:)));
[junk, xt] = find(td == max(td(:)));
[junk, xs] = find(sts == max(sts(:)));
axes('Position',[0.005,0.005,0.1,0.1],'Visible','off');
text(1,0.1,['MaxVal@scan: T:' sprintf('%01.1f', max(mot(:))) ' mm @ scan ' num2str(xm) '; R: ' sprintf('%01.1f', max(deg(:))) ' � @ scan ' num2str(xd) '; TD: '  sprintf('%01.1f', max(td(:))) ' mm @ scan ' num2str(xt) '; STS: '  sprintf('%01.1f', max(sts(:))) ' mm @ scan ' num2str(xs) '.'],'FontSize',8,'interpreter','none');

% ... motion fingerprint...
set(0,'CurrentFigure',gcf1);
subplot(nrows,3,[10 11 12],'replace');
plot(curr_tc);
hold on;
grid on;
xlim([0 size(pr,1)]);
ylim([floor(min(curr_tc(:))) ceil(max(curr_tc(:)))]);
if shifted == 0, withshift  = 'no';  else  withshift  = 'incl.';  end;
if squared == 0, withsquare = 'no';  else  withsquare = 'incl.';  end;
temp = oris_id(inds,:);
show = temp(1,:); for j = 2:size(temp,1), show = [show '/' temp(j,:)]; end
title(['Motion fingerprint (from ' show '), ' withshift ' shifted, ' withsquare ' squared versions']);

% add fancy plot?
% ... seems like we have a go...
set(0,'CurrentFigure',gcf1);


% compute combined motion in all directions
temp = zeros(size(pr,1),3);
for ii = 1:size(pr,1);
    temp(ii,1) = pr(ii,1) + (pr(ii,4) .* davg);
    temp(ii,2) = pr(ii,2) + (pr(ii,5) .* davg);
    temp(ii,3) = pr(ii,3) + (pr(ii,6) .* davg);
end;
maxd = max(abs(temp(:)));
if maxd < vs,  mybound = vs;  else mybound = ceil(maxd);  end;


% generate colormap, remove blue values (lower 50%)
cm = colormap(jet(mysteps+round(mysteps/2)));
cm(1:round(mysteps/2),:) = [];


% scale results to colormap, set values above vs to red and avoid <1
thr = vs;
store = round(td.*mysteps/thr);
store(store >= mysteps) = mysteps;
store(store < 1) = 1;


% ... and show
if nrows == 4
    
    subplot(nrows,3,10,'replace');
    
elseif nrows == 5
    
    subplot(nrows,3,13,'replace');
end;
gca1 = gca;
set(gca1,'XTick',[-round(vs) 0 round(vs)],'YTick',[-round(vs) 0 round(vs)],'DataAspectRatio',[1 1 1]);
xlim(gca1,[-round(mybound) round(mybound)]);
ylim(gca1,[-round(mybound) round(mybound)]);
box(gca1,'on');
grid(gca1,'on');
hold(gca1,'all');
for ii = 1:size(temp,1)
    
    scatter(temp(ii,2),temp(ii,3),[],cm(store(ii),:),'filled');
    line([0 temp(ii,2)],[0 temp(ii,3)], 'LineStyle',':','Color',[0.8 0.8 0.8]);
end;
title('TD in Y & Z...');


% ... and show
if nrows == 4
    
    subplot(nrows,3,11,'replace');
    
elseif nrows == 5
    
    subplot(nrows,3,14,'replace');
end;
gca1 = gca;
set(gca1,'XTick',[-round(vs) 0 round(vs)],'YTick',[-round(vs) 0 round(vs)],'DataAspectRatio',[1 1 1]);
xlim(gca1,[-round(mybound) round(mybound)]);
ylim(gca1,[-round(mybound) round(mybound)]);
box(gca1,'on');
grid(gca1,'on');
hold(gca1,'all');
for ii = 1:size(temp,1)
    
    scatter(temp(ii,1),temp(ii,3),[],cm(store(ii),:),'filled');
    line([0 temp(ii,1)],[0 temp(ii,3)], 'LineStyle',':','Color',[0.8 0.8 0.8]);
end;
title('... X & Z ...');


% ... and show
if nrows == 4
    
    subplot(nrows,3,12,'replace');
    
elseif nrows == 5
    
    subplot(nrows,3,15,'replace');
end;
gca1 = gca;
set(gca1,'XTick',[-round(vs) 0 round(vs)],'YTick',[-round(vs) 0 round(vs)],'DataAspectRatio',[1 1 1]);
xlim(gca1,[-round(mybound) round(mybound)]);
ylim(gca1,[-round(mybound) round(mybound)]);
box(gca1,'on');
grid(gca1,'on');
hold(gca1,'all');
for ii = 1:size(temp,1)
    
    scatter(temp(ii,1),temp(ii,2),[],cm(store(ii),:),'filled');
    line([0 temp(ii,1)],[0 temp(ii,2)], 'LineStyle',':','Color',[0.8 0.8 0.8]);
end;
title('... and X & Y.');

% print
motname = [p filesep 'mw_motion.jpg'];
osu = get(gcf1,'Units'); opu = get(gcf1,'PaperUnits'); opp = get(gcf1,'PaperPosition');
set(gcf1,'Units','pixels'); scrpos = get(gcf1,'Position'); newpos = scrpos/100;
set(gcf1,'PaperUnits','inches','PaperPosition',newpos);
print(gcf1, '-djpeg', '-noui', ['-r' num2str(pres)], motname);
close(gcf1);