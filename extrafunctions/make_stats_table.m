function make_stats_table(varargin)
%
% display and optionally save to jpeg an SPM thresholded stats table
%
% Usage:
%
%  make_stats_table(SPM, [fname, Ic, u, thresDesc, k, Im]);
%
% INPUT
%
%	SPM			- SPM struct containing thresholded results
%
% optional
%
%	fname		- jpeg name or path (omit or pass in [] to skip save)
%	Ic			- contrast # to examine
%	u			- threshold
%	thresDesc	- thresholding ('none' or 'FWE')
%	k			- cluster extent
%	Im			- masking option (currently ignored)
%
% if an optional parameter is missing, the following defaults are used:
%
%	no save (fname = [])
%	contrast #1 (Ic = 1)
%	significance = 0.001 (u = 0.001)
%	uncorrected (thresDesc = 'none')
%	voxel extent (k) = 0
%	no mask (Im = [])
%		
% sanity check

if (nargin < 1)
	fprintf('Usage: make_stats_table(SPM, [fname, Ic, u, thresDesc, k, Im])\n');
	return;
end

SPM = varargin{1};

% sanity check -- make sure the passed SPM has results to display 

if (~isfield(SPM,'xCon') || isempty(SPM.xCon))
	fprintf('%s: SPM struct does not contain contrasts. Aborting...\n', mfilename);
	return;
end

% defaults

fname = [];
SPM.Ic = 1;
SPM.u = 0.001;
SPM.thresDesc = 'none';
% 	SPM.u = 0.05;
% 	SPM.thresDesc = 'FWE';
SPM.k = 0;
SPM.Im = [];


if (nargin > 1)
	fname = varargin{2};
end

if (nargin > 2)
	SPM.Ic = varargin{3};
end

if (nargin > 3)
	SPM.u = varargin{4};
end

if (nargin > 4)
	SPM.thresDesc = varargin{5};
end

if (nargin > 5)
	SPM.k = varargin{6};
end

% Im is a ruse: currently we can only handle "no mask"

% 	if (nargin > 6)
%  		SPM.Im = varargin{7};
% 	end


% spm_list_display_noUI can only handle 'none' and 'FWE'

if (~strcmp(SPM.thresDesc,'none') && ~strcmp(SPM.thresDesc,'FWE'))
	fprintf('%s: Can only handle ''none'' and ''FWE'' thresholding. Aborting...\n', mfilename);
	return;
end

% extract info and make table

[~,xSPM] = spm_getSPM(SPM);
TabDat = spm_list('Table',xSPM);
hreg = spm_figure('CreateSatWin');
spm_list_display_noUI(TabDat,hreg);

% save?

if (~isempty(fname))

	set(hreg,'Renderer','opengl'); % I think this is the default

	% workaround for font rescaling weirdness

	set(findall(hreg,'Type','text'),'FontUnits','normalized');

	% tweak paper size to account for landscape layout

	set(hreg,'PaperUnits','inches','PaperPosition',[0 0 5 4]);

	% force jpg suffix otherwise there can be weirdness
	% also, bump up resolution to 200 bc numbers

	[p,n,~] = fileparts(fname);
	fname = fullfile(p,[n '.jpg']);

	print(hreg, '-djpeg', '-r200', fname);

	close(hreg);

end

end


function spm_list_display_noUI(TabDat,hReg)
%
% this is essentially spm_list('Display',...) w/ the parts that expect
% the SPM interactive windows to be up and the parts that don't
% play nice with save-to-jpeg gutted. See make_stats_table for
% usage
%

    %-Setup Graphics panel
    %----------------------------------------------------------------------
    Fgraph = spm_figure('FindWin','Satellite');
    if ~isempty(Fgraph)
        spm_figure('Focus',Fgraph);
        ht = 0.85; bot = 0.14;
    else
        Fgraph = spm_figure('GetWin','Graphics');
        ht = 0.4; bot = 0.1;
    end
    spm_results_ui('Clear',Fgraph)
    FS     = spm('FontSizes');           %-Scaled font sizes
    PF     = spm_platform('fonts');      %-Font names (for this platform)
    
    %-Table axes & Title
    %----------------------------------------------------------------------
    hAx   = axes('Parent',Fgraph,...
                 'Position',[0.025 bot 0.9 ht],...
                 'DefaultTextFontSize',FS(8),...
                 'DefaultTextInterpreter','Tex',...
                 'DefaultTextVerticalAlignment','Baseline',...
                 'Tag','SPMList',...
                 'Units','points',...
                 'Visible','off');

    AxPos = get(hAx,'Position'); set(hAx,'YLim',[0,AxPos(4)])
    dy    = FS(9);
    y     = floor(AxPos(4)) - dy;

% this is not playing well with jpeg save
%
%     text(0,y,['Statistics:  \it\fontsize{',num2str(FS(9)),'}',TabDat.tit],...
%               'FontSize',FS(11),'FontWeight','Bold');   y = y - dy/2;

	line([0 1],[y y],'LineWidth',3,'Color','r'),        y = y - 9*dy/8;
    
    %-Display table header
    %----------------------------------------------------------------------
    set(hAx,'DefaultTextFontName',PF.helvetica,'DefaultTextFontSize',FS(8))

    Hs = []; Hc = []; Hp = [];
    h  = text(0.01,y, [TabDat.hdr{1,1} '-level'],'FontSize',FS(9)); Hs = [Hs,h];
    h  = line([0,0.11],[1,1]*(y-dy/4),'LineWidth',0.5,'Color','r'); Hs = [Hs,h];
    h  = text(0.02,y-9*dy/8,    TabDat.hdr{3,1});              Hs = [Hs,h];
    h  = text(0.08,y-9*dy/8,    TabDat.hdr{3,2});              Hs = [Hs,h];
    
    h = text(0.22,y, [TabDat.hdr{1,3} '-level'],'FontSize',FS(9));    Hc = [Hc,h];
    h = line([0.14,0.44],[1,1]*(y-dy/4),'LineWidth',0.5,'Color','r'); Hc = [Hc,h];
    h  = text(0.15,y-9*dy/8,    TabDat.hdr{3,3});              Hc = [Hc,h];
    h  = text(0.24,y-9*dy/8,    TabDat.hdr{3,4});              Hc = [Hc,h];
    h  = text(0.34,y-9*dy/8,    TabDat.hdr{3,5});              Hc = [Hc,h];
    h  = text(0.39,y-9*dy/8,    TabDat.hdr{3,6});              Hc = [Hc,h];
    
    h = text(0.64,y, [TabDat.hdr{1,7} '-level'],'FontSize',FS(9));    Hp = [Hp,h];
    h = line([0.48,0.88],[1,1]*(y-dy/4),'LineWidth',0.5,'Color','r'); Hp = [Hp,h];
    h  = text(0.49,y-9*dy/8,    TabDat.hdr{3,7});              Hp = [Hp,h];
    h  = text(0.58,y-9*dy/8,    TabDat.hdr{3,8});              Hp = [Hp,h];
    h  = text(0.67,y-9*dy/8,    TabDat.hdr{3,9});              Hp = [Hp,h];
    h  = text(0.75,y-9*dy/8,    TabDat.hdr{3,10});             Hp = [Hp,h];
    h  = text(0.82,y-9*dy/8,    TabDat.hdr{3,11});             Hp = [Hp,h];
    
    text(0.92,y - dy/2,TabDat.hdr{3,12},'Fontsize',FS(8));

    %-Move to next vertical position marker
    %----------------------------------------------------------------------
    y     = y - 7*dy/4;
    line([0 1],[y y],'LineWidth',1,'Color','r')
    y     = y - 5*dy/4;
    y0    = y;

    %-Table filtering note
    %----------------------------------------------------------------------
    text(0.5,4,TabDat.str,'HorizontalAlignment','Center',...
        'FontName',PF.helvetica,'FontSize',FS(8),'FontAngle','Italic')

    %-Footnote with SPM parameters (if classical inference)
    %----------------------------------------------------------------------
    line([0 1],[0.01 0.01],'LineWidth',1,'Color','r')
    if ~isempty(TabDat.ftr)
        set(gca,'DefaultTextFontName',PF.helvetica,...
            'DefaultTextInterpreter','None','DefaultTextFontSize',FS(8))
        
        fx = repmat([0 0.5],ceil(size(TabDat.ftr,1)/2),1);
        fy = repmat((1:ceil(size(TabDat.ftr,1)/2))',1,2);
        for i=1:size(TabDat.ftr,1)
            text(fx(i),-fy(i)*dy,sprintf(TabDat.ftr{i,1},TabDat.ftr{i,2}),...
                'UserData',TabDat.ftr{i,2},...
                'ButtonDownFcn','get(gcbo,''UserData'')');
        end
    end
    
    %-Characterize excursion set in terms of maxima
    % (sorted on Z values and grouped by regions)
    %======================================================================
    if isempty(TabDat.dat)
        text(0.5,y-6*dy,'no suprathreshold clusters',...
            'HorizontalAlignment','Center',...
            'FontAngle','Italic','FontWeight','Bold',...
            'FontSize',FS(16),'Color',[1,1,1]*.5);
        return
    end
    
    %-Table proper
    %======================================================================

    %-Column Locations
    %----------------------------------------------------------------------
    tCol = [ 0.01      0.08 ...                                %-Set
             0.15      0.24      0.33      0.39 ...            %-Cluster
             0.49      0.58      0.65      0.74      0.83 ...  %-Peak
             0.92];                                            %-XYZ
    
    %-Pagination variables
    %----------------------------------------------------------------------
    hPage = [];
    set(gca,'DefaultTextFontName',PF.courier,'DefaultTextFontSize',FS(7));

    %-Set-level p values {c} - do not display if reporting a single cluster
    %----------------------------------------------------------------------
    if isempty(TabDat.dat{1,1}) % Pc
        set(Hs,'Visible','off');
    end
    
    if TabDat.dat{1,2} > 1 % c
        h     = text(tCol(1),y,sprintf(TabDat.fmt{1},TabDat.dat{1,1}),...
                    'FontWeight','Bold', 'UserData',TabDat.dat{1,1},...
                    'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];
        h     = text(tCol(2),y,sprintf(TabDat.fmt{2},TabDat.dat{1,2}),...
                    'FontWeight','Bold', 'UserData',TabDat.dat{1,2},...
                    'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];
    else
        set(Hs,'Visible','off');
    end
    
    %-Cluster and local maxima p-values & statistics
    %----------------------------------------------------------------------
    HlistXYZ   = [];
    HlistClust = [];
    for i=1:size(TabDat.dat,1)
        
        %-Paginate if necessary
        %------------------------------------------------------------------
        if y < dy
            h = text(0.5,-5*dy,...
                sprintf('Page %d',spm_figure('#page',Fgraph)),...
                        'FontName',PF.helvetica,'FontAngle','Italic',...
                        'FontSize',FS(8));
            spm_figure('NewPage',[hPage,h])
            hPage = [];
            y     = y0;
        end
        
        %-Print cluster and maximum peak-level p values
        %------------------------------------------------------------------
        if  ~isempty(TabDat.dat{i,5}), fw = 'Bold'; else fw = 'Normal'; end
        
        for k=3:11
            h = text(tCol(k),y,sprintf(TabDat.fmt{k},TabDat.dat{i,k}),...
                     'FontWeight',fw,...
                     'UserData',TabDat.dat{i,k},...
                     'ButtonDownFcn','get(gcbo,''UserData'')');
            hPage = [hPage, h];
            if k == 5
                HlistClust = [HlistClust, h];
                set(h,'UserData',struct('k',TabDat.dat{i,k},'XYZmm',TabDat.dat{i,12}));
                set(h,'ButtonDownFcn','getfield(get(gcbo,''UserData''),''k'')');
            end
        end
        
        % Specifically changed so it properly finds hMIPax
        %------------------------------------------------------------------
        tXYZmm = TabDat.dat{i,12};
        BDFcn  = [...
            'spm_mip_ui(''SetCoords'',get(gcbo,''UserData''),',...
                'findobj(''tag'',''hMIPax''));'];
        BDFcn = 'spm_XYZreg(''SetCoords'',get(gcbo,''UserData''),hReg,1);';
        h = text(tCol(12),y,sprintf(TabDat.fmt{12},tXYZmm),...
            'FontWeight',fw,...
            'Tag','ListXYZ',...
            'ButtonDownFcn',BDFcn,...
            'Interruptible','off',...
            'BusyAction','Cancel',...
            'UserData',tXYZmm);

        HlistXYZ = [HlistXYZ, h];
        hPage  = [hPage, h];
        y      = y - dy;
    end
    
    %-Number and register last page (if paginated)
    %----------------------------------------------------------------------
    if spm_figure('#page',Fgraph)>1
        h = text(0.5,-5*dy,sprintf('Page %d/%d',spm_figure('#page',Fgraph)*[1,1]),...
            'FontName',PF.helvetica,'FontSize',FS(8),'FontAngle','Italic');
        spm_figure('NewPage',[hPage,h])
	end 
end

