function SliceView(IN)

% IN.IM = {which('defaultUnderlay.nii') '02_DefaultMode.nii'};
% IN.H = {[] []};
% IN.TH = {[-inf inf] [100  inf]};
% IN.LIMS = {[20 120] [90 inf]};
% IN.TRANS = {1 1};
% IN.CM = {'gray' 'hot'};
% IN.Coords = -52:6:84;
% IN.opt = [1 2 3];

if ~isfield(IN,'opt') || isempty(IN.opt)
    IN.opt = 1:3;
end


coords = IN.Coords;
n = numel(coords);

MH = spm_vol(which('defaultUnderlay.nii'));

for ii = 1:numel(IN.IM)
    if ischar(IN.IM{ii})
        h = spm_vol(IN.IM{ii});
        [N newTarg] = SliceAndDice3(h,MH,h,h,[3 NaN],[]);
        h.mat = newTarg;
        IN.H{ii} = h;
        IN.IM{ii} = N;
    end
    
    i1 = find(IN.TH{ii}==inf);
    i2 = find(IN.TH{ii}==-inf);
    
    IN.TH{ii}(i1) = max(IN.IM{ii}(:));
    IN.TH{ii}(i2) = min(IN.IM{ii}(:));
    
    i1 = find(IN.LIMS{ii}==inf);
    i2 = find(IN.LIMS{ii}==-inf);
    
    IN.LIMS{ii}(i1) = max(IN.IM{ii}(:));
    IN.LIMS{ii}(i2) = min(IN.IM{ii}(:));
    
    sss = size(IN.IM{ii});
    pos = [];
    t1 = [1 1 1 1; sss(1) 1 1 1]*IN.H{ii}.mat';
    pos{1} = t1(1,1):(t1(2,1)-t1(1,1))/(sss(1)-1):t1(2,1);
    t1 = [1 1 1 1; 1 sss(2) 1 1]*IN.H{ii}.mat';
    pos{2} = t1(1,2):(t1(2,2)-t1(1,2))/(sss(2)-1):t1(2,2);
    t1 = [1 1 1 1; 1 1 sss(3) 1]*IN.H{ii}.mat';
    pos{3} = t1(1,3):(t1(2,3)-t1(1,3))/(sss(3)-1):t1(2,3);
    IN.POS{ii} =  pos;
end

opt = IN.opt(:)';

nn1 = ceil(sqrt(n));
nn2 = ceil(n/nn1);
inc1 = 1/nn1;
inc2 = 1/nn2;

for zz = opt
    nf = figure;
    set(nf,'color','k');
    ss = get(0,'ScreenSize');
    sss = size(IN.IM{1});

    switch zz
        case 1
            set(nf,'Position',[50 50 (sss(1)/sss(2))*(ss(4)-150) ss(4)-150],'Visible','off');
        case 2
            set(nf,'Position',[50 50 floor((sss(2)/sss(3))*(ss(4)-150)) ss(4)-150],'Visible','off');
        case 3
            set(nf,'Position',[50 50 (sss(1)/sss(3))*(ss(4)-150) ss(4)-150],'Visible','off');
        otherwise
            warning('opt must be 1 (Axial), 2 (Sagittal), or 3 (Coronal)');
            return
    end
    
    txt = [];
    c = 0;
    for ii = coords
        c = c+1;
        for jj = 1:numel(IN.IM)
            
            tx(c) = subplot(nn2,nn1,c);
            
            switch zz
                case 1
                    j = [0 0 coords(c) 1]/(IN.H{jj}.mat');
                    j = round(j(3));
                    dat  = flip(rot90(squeeze(IN.IM{jj}(:,:,j)),1),1);
                case 2
                    j = [coords(c) 0 0 1]/(IN.H{jj}.mat');
                    j = round(j(1));
                    dat = flip(flip(rot90(squeeze(IN.IM{jj}(j,:,:)),1),1),2);
                case 3
                    j = [0 coords(c) 0 1]/(IN.H{jj}.mat');
                    j = round(j(2));
                    dat = flip(rot90(squeeze(IN.IM{jj}(:,j,:)),1),1);
                otherwise
                    %close(nf);
                    warning('opt must be 1 (Axial), 2 (Sagittal), or 3 (Coronal)');
                    return
            end
            
            
            
            
            
            if numel(IN.TH{jj})==2
                ind = find(dat<IN.TH{jj}(1) | dat>IN.TH{jj}(2));
            else
                ind = find(dat<IN.TH{jj}(1) | dat>IN.TH{jj}(4) | (dat>IN.TH{jj}(2) & dat<IN.TH{jj}(3)));
            end
            dat(ind) = NaN;
            
            [cols cm cc] = cmap(dat, IN.LIMS{jj}, IN.CM{jj});
            col = reshape(cols,[size(dat) 3]);
            
            
            switch zz
                case 1
                    xx = image(IN.POS{jj}{1}, IN.POS{jj}{2}, col); hold on;
                case 2
                    xx = image(IN.POS{jj}{2}, IN.POS{jj}{3}, col); hold on;
                case 3
                    xx = image(IN.POS{jj}{1}, IN.POS{jj}{3}, col); hold on;
                otherwise
                    warning('opt must be 1 (Axial), 2 (Sagittal), or 3 (Coronal)');
                    return
            end
            
            set(xx,'AlphaData', double(~isnan(dat))*IN.TRANS{jj});
            set(gca,'XTick',[],'YTick',[]);
            set(gca, 'YDir','Normal','XDir','Normal');
            
            axis off;
            axis tight;
            axis equal;
            
        end
        ax = axis;
        cc = [ax(1)+(diff(ax(1:2))*0) ax(4)-(diff(ax(3:4))*.04)];
        txt(c) = text(cc(1),cc(2),1,num2str(ii),'color','w','fontsize',18,'fontweight','bold','fontname','times new roman');
    end

    
    c = 0;
    for nn = nn2:-1:1
        for mm = 1:nn1
            c = c+1;
            pos = [inc1*(mm-1) inc2*(nn)-inc2 inc1 inc2];
            if c>numel(tx)
                set(nf,'Visible','on');
                continue;
            end
            set(tx(c),'Units', 'Normalized');
            set(tx(c),'Position',pos);
        end
    end
    set(nf,'Visible','on');
end
