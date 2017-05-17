function varargout = cg_tfce_list(varargin)
% Display an analysis of SPM{.}
% FORMAT TabDat = cg_tfce_list('List',SPM,hReg,[Num,Dis,Str])
% Summary list of local maxima for entire volume of interest
% FORMAT TabDat = cg_tfce_list('ListCluster',SPM,hReg,[Num,Dis,Str])
% List of local maxima for a single suprathreshold cluster
%
% SPM    - structure containing SPM, distribution & filtering details
%        - required fields are:
% .swd   - SPM working directory - directory containing current SPM.mat
% .Z     - minimum of n Statistics {filtered on u and k}
% .n     - number of conjoint tests
% .STAT  - distribution {Z, T, X or F}
% .df    - degrees of freedom [df{interest}, df{residual}]
% .u     - height threshold
% .k     - extent threshold {voxels}
% .XYZ   - location of voxels {voxel coords}
% .S     - search Volume {voxels}
% .R     - search Volume {resels}
% .FWHM  - smoothness {voxels}
% .M     - voxels - > mm matrix
% .VOX   - voxel dimensions {mm}
% .DIM   - image dimensions {voxels}
% .units - space units
% .VRpv  - filehandle - Resels per voxel
% .Ps    - uncorrected P values in searched volume (for voxel FDR)
% .Pp    - uncorrected P values of peaks (for peak FDR)
% .Pc    - uncorrected P values of cluster extents (for cluster FDR)
% .uc    - 0.05 critical thresholds for FWEp, FDRp, FWEc, FDRc
% .thresDesc - description of height threshold (string)
%
% (see spm_getSPM for further details of xSPM structures)
%
% hReg   - Handle of results section XYZ registry (see spm_results_ui.m)
%
% Num    - number of maxima per cluster
% Dis    - distance among clusters (mm)
% Str    - header string
%
% TabDat - Structure containing table data
%        - fields are
% .tit   - table Title (string)
% .hdr   - table header (2x12 cell array)
% .fmt   - fprintf format strings for table data (1x12 cell array)
% .str   - table filtering note (string)
% .ftr   - table footnote information (5x2 cell array)
% .dat   - table data (Nx12 cell array)
%
%                           ----------------
%
% FORMAT cg_tfce_list('TxtList',TabDat,c)
% Prints a tab-delimited text version of the table
% TabDat - Structure containing table data (format as above)
% c      - Column of table data to start text table at
%          (E.g. c=3 doesn't print set-level results contained in columns 1 & 2)
%                           ----------------
%
% FORMAT cg_tfce_list('SetCoords',xyz,hAx,hC)
% Highlighting of table co-ordinates (used by results section registry)
% xyz    - 3-vector of new co-ordinate
% hAx    - table axis (the registry object for tables)
% hReg   - Handle of caller (not used)
%__________________________________________________________________________
%
% cg_tfce_list characterizes SPMs (thresholded at u and k) in terms of
% excursion sets (a collection of face, edge and vertex connected
% subsets or clusters).  The corrected significance of the results are
% based on set, cluster and voxel-level inferences using distributional
% approximations from the Theory of Gaussian Fields.  These
% distributions assume that the SPM is a reasonable lattice
% approximation of a continuous random field with known component field
% smoothness.
%
% The p values are based on the probability of obtaining c, or more,
% clusters of k, or more, resels above u, in the volume S analysed =
% P(u,k,c).  For specified thresholds u, k, the set-level inference is
% based on the observed number of clusters C, = P(u,k,C).  For each
% cluster of size K the cluster-level inference is based on P(u,K,1)
% and for each voxel (or selected maxima) of height U, in that cluster,
% the voxel-level inference is based on P(U,0,1).  All three levels of
% inference are supported with a tabular presentation of the p values
% and the underlying statistic:
%
% Set-level     - c    = number of suprathreshold clusters
%               - P    = prob(c or more clusters in the search volume)
%
% Cluster-level - k    = number of voxels in this cluster
%               - Pc   = prob(k or more voxels in the search volume)
%               - Pu   = prob(k or more voxels in a cluster)
%               - Qc   = lowest FDR bound for which this cluster would be
%                        declared positive
%
% Peak-level    - T/F  = Statistic upon which the SPM is based
%               - Ze   = The equivalent Z score - prob(Z > Ze) = prob(t > T)
%               - Pc   = prob(Ze or higher in the search volume)
%               - Qp   = lowest FDR bound for which this peak would be
%                        declared positive
%               - Pu   = prob(Ze or higher at that voxel)
%
% Voxel-level   - Qu   = Expd(Prop of false positives among voxels >= Ze)
%
% x,y,z (mm)    - Coordinates of the voxel
%
% The table is grouped by regions and sorted on the Ze-variate of the
% primary maxima.  Ze-variates (based on the uncorrected p value) are the
% Z score equivalent of the statistic. Volumes are expressed in voxels.
%
% Clicking on values in the table returns the value to the Matlab
% workspace. In addition, clicking on the co-ordinates jumps the
% results section cursor to that location. The table has a context menu
% (obtained by right-clicking in the background of the table),
% providing options to print the current table as a text table, or to
% extract the table data to the Matlab workspace.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston & Andrew Holmes
% $Id: cg_tfce_list.m 85 2015-12-01 10:55:55Z gaser $

% always voxel-wise FDR
topoFDR = false;

%==========================================================================
switch lower(varargin{1}), case 'list'                            %-List
%==========================================================================
% FORMAT TabDat = cg_tfce_list('list',SPM,hReg)

    %-Tolerance for p-value underflow, when computing equivalent Z's
    %----------------------------------------------------------------------
    tol = eps*10;

    %-Parse arguments and set maxima number and separation
    %----------------------------------------------------------------------
    if nargin < 2,  error('insufficient arguments'),     end
    if nargin < 3,  hReg = []; else  hReg = varargin{3}; end
    xSPM = varargin{2};
    if isempty(xSPM), varargout = {{}}; return;          end

    %-Get current location (to highlight selected voxel in table)
    %----------------------------------------------------------------------
    xyzmm     = spm_results_ui('GetCoords');

    %-Extract data from xSPM
    %----------------------------------------------------------------------
    S      = varargin{2}.S;
    VOX    = varargin{2}.VOX;
    DIM    = varargin{2}.DIM;
    n      = varargin{2}.n;
    XYZ    = varargin{2}.XYZ;
    XYZmm  = varargin{2}.XYZmm;
    STAT   = varargin{2}.STAT;
    df     = varargin{2}.df;
    u      = varargin{2}.u;
    M      = varargin{2}.M;
    k      = varargin{2}.k;
    n_perm = varargin{2}.n_perm;
    try, QPs = varargin{2}.Ps; end
    try, QPp = varargin{2}.Pp; end
    try, QPc = varargin{2}.Pc; end
    try
        thresDesc = sprintf('{%s}', varargin{2}.thresDesc);
    catch
        thresDesc = '';
    end
    
    if STAT~='P'
        R     = full(varargin{2}.R);
        FWHM  = full(varargin{2}.FWHM);
    end
    try
        units = varargin{2}.units;
    catch
        units = {'mm' 'mm' 'mm'};
    end
    units{1}  = [units{1} ' '];
    units{2}  = [units{2} ' '];

    DIM   = DIM > 1;                  % non-empty dimensions
    if strcmp(spm('ver'),'SPM12')
        if spm_mesh_detect(xSPM.Vspm)
            DIM   = true(1,3);        % non-empty dimensions
        end
    end
    VOX       = VOX(DIM);             % scaling

    if STAT~='P'
        R     = full(xSPM.R);         % Resel counts
        FWHM  = FWHM(DIM);            % Full width at max/2
        FWmm  = FWHM.*VOX;            % FWHM {units}
        V2R   = 1/prod(FWHM);         % voxels to resels
        k     = k*V2R;                % extent threshold in resels
        R     = R(1:find(R~=0,1,'last')); % eliminate null resel counts
        try, QPs = sort(QPs(:)); end  % Needed for voxel   FDR
        try, QPp = sort(QPp(:)); end  % Needed for peak    FDR
        try, QPc = sort(QPc(:)); end  % Needed for cluster FDR
    end

    %-get number and separation for maxima to be reported
    %----------------------------------------------------------------------
    if length(varargin) > 3
        Num    = varargin{4};         % number of maxima per cluster
        Dis    = varargin{5};         % distance among clusters (mm)
    else
        Num    = 3;
        Dis    = 8;
    end
    if length(varargin) > 5
        Title  = varargin{6};
    else
        Title  = 'nonparametric p-values adjusted for search volume';
    end

    %-Setup graphics panel
    %----------------------------------------------------------------------
    spm('Pointer','Watch')
    Fgraph = spm_figure('FindWin','Satellite');
    if ~isempty(Fgraph)
        figure(Fgraph);
        ht = 0.85; bot = 0.14;
    else
        Fgraph = spm_figure('GetWin','Graphics');
        ht = 0.4; bot = 0.1;
    end
    spm_results_ui('Clear',Fgraph)
    FS    = spm('FontSizes');           %-Scaled font sizes
    PF    = spm_platform('fonts');      %-Font names (for this platform)


    %-Table header & footer
    %======================================================================

    %-Table axes & Title
    %----------------------------------------------------------------------

    if STAT == 'P'
        Title = 'Posterior Probabilities';
    end

    hAx   = axes('Position',[0.025 bot 0.9 ht],...
                'DefaultTextFontSize',FS(8),...
                'DefaultTextInterpreter','Tex',...
                'DefaultTextVerticalAlignment','Baseline',...
                'Tag','SPMList',...
                'Units','points',...
                'Visible','off');

    AxPos = get(hAx,'Position'); set(hAx,'YLim',[0,AxPos(4)])
    dy    = FS(9);
    y     = floor(AxPos(4)) - dy;

    text(0,y,['Statistics:  \it\fontsize{',num2str(FS(9)),'}',Title],...
              'FontSize',FS(11),'FontWeight','Bold');   y = y - dy/2;
    line([0 1],[y y],'LineWidth',3,'Color','b'),        y = y - 9*dy/8;

    %-Construct table header
    %----------------------------------------------------------------------
    set(gca,'DefaultTextFontName',PF.helvetica,'DefaultTextFontSize',FS(8))

    Hc = [];
    Hp = [];
    text(0.22,y,        'cluster-level','FontSize',FS(9));
    line([0.14,0.44],[1,1]*(y-dy/4),'LineWidth',0.5,'Color','b');
    h  = text(0.34,y-9*dy/8,    '\itk\rm_E');
    
    if STAT=='TFCE'
        text(0.55,y,        'combined peak-cluster-level','FontSize',FS(9));
    else
        text(0.55,y,        'peak-level','FontSize',FS(9));
    end
    line([0.48,0.88],[1,1]*(y-dy/4),'LineWidth',0.5,'Color','b');
    h  = text(0.49,y-9*dy/8,    '\itp\rm_{FWE-corr}');     Hp = [Hp,h];
    h  = text(0.58,y-9*dy/8,        '\itq\rm_{FDR-corr}'); Hp = [Hp,h];
    h  = text(0.82,y-9*dy/8,    '\itp\rm_{uncorr}');       Hp = [Hp,h];
    h  = text(0.67,y-9*dy/8,     sprintf('\\it%c',STAT));
%    h  = text(0.75,y-9*dy/8,    '(\itZ\rm_\equiv)');
    
    text(0.92,y - dy/2,[units{:}],'Fontsize',FS(8));


    %-Headers for text table...
    %-----------------------------------------------------------------------
    TabDat.tit = Title;
    TabDat.hdr = {...
        'set',      'p';...
        'set',      'c';...
        'cluster',  'p(FWE-cor)';...
        'cluster',  'p(FDR-cor)';...
        'cluster',  'equivk';...
        'cluster',  'p(unc)';...
        'peak',     'p(FWE-cor)';...
        'peak',     'p(FDR-cor)';...
        'peak',      STAT;...
        'peak',     'equivZ';...
        'peak',     'p(unc)';...
        '',         'x,y,z {mm}'}';...

    %-Coordinate Precisions
    %----------------------------------------------------------------------
    if isempty(XYZmm) % empty results
        xyzfmt = '%3.0f %3.0f %3.0f';
        voxfmt = repmat('%0.1f ',1,nnz(DIM));
    elseif ~any(strcmp(units{3},{'mm',''})) % 2D data
        xyzfmt = '%3.0f %3.0f %3.0f';
        voxfmt = repmat('%0.1f ',1,nnz(DIM));
    else % 3D data, work out best precision based on voxel sizes and FOV
        xyzsgn = min(XYZmm(DIM,:),[],2) < 0;
        xyzexp = max(floor(log10(max(abs(XYZmm(DIM,:)),[],2)))+(max(abs(XYZmm(DIM,:)),[],2) >= 1),0);
        voxexp = floor(log10(abs(VOX')))+(abs(VOX') >= 1);
        xyzdec = max(-voxexp,0);
        voxdec = max(-voxexp,1);
        xyzwdt = xyzsgn+xyzexp+(xyzdec>0)+xyzdec;
        voxwdt = max(voxexp,0)+(voxdec>0)+voxdec;
        tmpfmt = cell(size(xyzwdt));
        for i = 1:numel(xyzwdt)
            tmpfmt{i} = sprintf('%%%d.%df ', xyzwdt(i), xyzdec(i));
        end
        xyzfmt = [tmpfmt{:}];
        tmpfmt = cell(size(voxwdt));
        for i = 1:numel(voxwdt)
            tmpfmt{i} = sprintf('%%%d.%df ', voxwdt(i), voxdec(i));
        end
        voxfmt = [tmpfmt{:}];
    end
    TabDat.fmt = {  '%-0.3f','%g',...                          %-Set
        '%0.3f', '%0.3f','%0.0f', '%0.3f',...                  %-Cluster
        '%0.3f', '%0.3f', '%6.2f', '%5.2f', '%0.3f',...        %-Peak
        '%3.0f %3.0f %3.0f'};                                  %-XYZ

    %-Column Locations
    %----------------------------------------------------------------------
    tCol = [ 0.01      0.08 ...                                %-Set
             0.15      0.24      0.33      0.39 ...            %-Cluster
             0.49      0.58      0.65      0.74      0.83 ...  %-Peak
             0.92];                                            %-XYZ

    %-Move to next vertical position marker
    %----------------------------------------------------------------------
    y     = y - 7*dy/4;
    line([0 1],[y y],'LineWidth',1,'Color','b')
    y     = y - 5*dy/4;
    y0    = y;


    %-Table filtering note
    %----------------------------------------------------------------------
    if isinf(Num)
        TabDat.str = sprintf('table shows all local maxima > %.1fmm apart',Dis);
    else
        TabDat.str = sprintf(['table shows %d local maxima ',...
            'more than %.1fmm apart'],Num,Dis);
    end
    text(0.5,4,TabDat.str,'HorizontalAlignment','Center','FontName',PF.helvetica,...
        'FontSize',FS(8),'FontAngle','Italic')


    %-Volume, resels and smoothness (if classical inference)
    %----------------------------------------------------------------------
    line([0 1],[0 0],'LineWidth',1,'Color','b')
        
    %-Volume, resels and smoothness (if classical inference)
    %----------------------------------------------------------------------
    line([0 1],[0 0],'LineWidth',1,'Color','b')
    if STAT ~= 'P'
        %-Footnote with SPM parameters
        %------------------------------------------------------------------
        set(gca,'DefaultTextFontName',PF.helvetica,...
            'DefaultTextInterpreter','None','DefaultTextFontSize',FS(8))
        vx = 'voxels';
        if strcmp(spm('ver'),'SPM12')
            if spm_mesh_detect(xSPM.Vspm)
                vx = 'vertices'; 
            end
        end
        TabDat.ftr    = cell(5,2);
        TabDat.ftr{1} = ...
            sprintf('Degrees of freedom = [%0.1f, %0.1f]',df);

        TabDat.ftr{2,1} = ...
            ['FWHM = ' voxfmt units{:} '; ' voxfmt '{' vx '}'];
        TabDat.ftr{2,2} = [FWmm FWHM];
        TabDat.ftr{3,1} = ...
            ['Volume: %0.0f = %0.0f ' vx ' = %0.1f resels'];
        TabDat.ftr{3,2} = [S*prod(VOX),S,R(end)];
        TabDat.ftr{4,1} = ...
            ['Voxel size: ' voxfmt units{:} '; (resel = %0.2f ' vx ')'];
        TabDat.ftr{4,2} = [VOX,prod(FWHM)];

        if strcmp(spm('ver'),'SPM12')
            if spm_mesh_detect(xSPM.Vspm)
                TabDat.ftr{2,1} = ...
                    ['FWHM = ' voxfmt '{' vx '}'];
                TabDat.ftr{2,2} = FWHM;
                TabDat.ftr{3,1} = ['Volume: %0.0f ' vx ' = %0.1f resels'];
                TabDat.ftr{3,2} = [S,R(end)];
                TabDat.ftr{4,1} = ['(resel = %0.2f ' vx ')'];
                TabDat.ftr{4,2} = prod(FWHM);
            end
        end

        TabDat.ftr{5} = ...
            sprintf('Permutations = %d',n_perm);

        set(gca,'DefaultTextFontName',PF.helvetica,...
            'DefaultTextInterpreter','None','DefaultTextFontSize',FS(8));
        fx = repmat([0 0.5],ceil(size(TabDat.ftr,1)/2),1);
        fy = repmat((1:ceil(size(TabDat.ftr,1)/2))',1,2);
        for i=1:size(TabDat.ftr,1)
            text(fx(i),-fy(i)*dy,sprintf(TabDat.ftr{i,1},TabDat.ftr{i,2}),...
                'UserData',TabDat.ftr{i,2},...
                'ButtonDownFcn','get(gcbo,''UserData'')',...
                'FontName',PF.helvetica,'FontSize',FS(8));
        end
    else
        TabDat.ftr = {};
    end


    %-Characterize excursion set in terms of maxima
    % (sorted on Z values and grouped by regions)
    %======================================================================
    if isempty(varargin{2}.Z)
        text(0.5,y-6*dy,'no suprathreshold clusters',...
            'HorizontalAlignment','Center',...
            'FontAngle','Italic','FontWeight','Bold',...
            'FontSize',FS(16),'Color',[1,1,1]*.5);
        TabDat.dat = cell(0,12);
        varargout  = {TabDat};
        spm('Pointer','Arrow')
        return
    end

    %-Workaround in spm_max for conjunctions with negative thresholds
    %----------------------------------------------------------------------
    minz          = abs(min(min(varargin{2}.Z)));
    Z       = 1 + minz + varargin{2}.Z;
    if strcmp(spm('ver'),'SPM12')
        if ~spm_mesh_detect(xSPM.Vspm)
            [N,Z,XYZ,A]  = spm_max(Z,XYZ);
        else
            [N,Z,XYZ,A]  = spm_mesh_max(Z,varargin{2}.XYZ,xSPM.G);
        end
    else 
        [N,Z,XYZ,A]  = spm_max(Z,XYZ);
    end
    
    Z             = Z - minz - 1;

    % find corresponding p-values for Z
    if strcmp(spm('ver'),'SPM12')
        Qu  = spm_data_read(varargin{2}.VQu,'xyz',XYZ);
        Pz  = spm_data_read(varargin{2}.VPz,'xyz',XYZ);
        Pu  = spm_data_read(varargin{2}.VPu,'xyz',XYZ);
    else
        Qu  = spm_get_data(varargin{2}.VQu,XYZ);
        Pz  = spm_get_data(varargin{2}.VPz,XYZ);
        Pu  = spm_get_data(varargin{2}.VPu,XYZ);
    end
    Qu(find(Qu<0)) = 0;
    Pz(find(Pz<0)) = 0;
    Pu(find(Pu<0)) = 0;
    
    % convert from -log10
    Qu = 10.^-Qu;
    Pu = 10.^-Pu;
    Pz = 10.^-Pz;

    %-Convert maxima locations from voxels to mm
    %----------------------------------------------------------------------
    XYZmm = M(1:3,:)*[XYZ; ones(1,size(XYZ,2))];
    if strcmp(spm('ver'),'SPM12')
        if spm_mesh_detect(xSPM.Vspm)
            XYZmm = xSPM.G.vertices(XYZ(1,:),:)';
        end
    end



    %-Table proper (& note all data in cell array)
    %======================================================================

    %-Pagination variables
    %----------------------------------------------------------------------
    hPage = [];
    set(gca,'DefaultTextFontName',PF.courier,'DefaultTextFontSize',FS(7))

    TabLin     = 1;                 %-Table data line


    %-Local maxima p-values & statistics
    %----------------------------------------------------------------------
    HlistXYZ = [];
    while numel(find(isfinite(Z)))

        % Paginate if necessary
        %------------------------------------------------------------------
        if y < min(Num + 1,3)*dy

            % added Fgraph term to paginate on Satellite window
            %--------------------------------------------------------------
            h     = text(0.5,-6*dy,...
                sprintf('Page %d',spm_figure('#page',Fgraph)),...
                'FontName',PF.helvetica,'FontAngle','Italic',...
                'FontSize',FS(8));

            spm_figure('NewPage',[hPage,h])
            hPage = [];
            y     = y0;
        end

        %-Find largest remaining local maximum
        %------------------------------------------------------------------
        [U,i]   = max(Z);           % largest maxima
        j       = find(A == A(i));  % maxima in cluster


        %-Compute cluster {k} and peak-level {u} p values for this cluster
        %------------------------------------------------------------------
        if STAT ~= 'P'
            Pk      = [];
            Pn      = [];
            Qc      = [];
            Qp      = [];

            if Pz < tol                               % Equivalent Z-variate
                Ze  = Inf;                            % (underflow => can't compute)
            else
                Ze  = spm_invNcdf(1 - Pz);
            end
        else
            
            Pz      = [];
            Pu      = [];
            Qu      = [];
            Pk      = [];
            Pn      = [];
            Qc      = [];
            Qp      = [];
            Ze      = spm_invNcdf(U);
        end


        %-Print cluster and maximum peak-level p values {Z}
        %------------------------------------------------------------------
        h     = text(tCol(3),y,sprintf(TabDat.fmt{3},Pk),'FontWeight','Bold',...
            'UserData',Pk,'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];
        h     = text(tCol(4),y,sprintf(TabDat.fmt{4},Qc),'FontWeight','Bold',...
            'UserData',Qc,'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];
        h     = text(tCol(5),y,sprintf(TabDat.fmt{5},N(i)),'FontWeight','Bold',...
            'UserData',N(i),'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];
        h     = text(tCol(6),y,sprintf(TabDat.fmt{6},Pn),'FontWeight','Bold',...
            'UserData',Pn,'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];
        h     = text(tCol(7),y,sprintf(TabDat.fmt{7},Pu(i)),'FontWeight','Bold',...
            'UserData',Pu(i),'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];
        if topoFDR
            h     = text(tCol(8),y,sprintf(TabDat.fmt{8},Qp),'FontWeight','Bold',...
            'UserData',Qp,'ButtonDownFcn','get(gcbo,''UserData'')');
        else
            h     = text(tCol(8),y,sprintf(TabDat.fmt{8},Qu(i)),'FontWeight','Bold',...
            'UserData',Qu(i),'ButtonDownFcn','get(gcbo,''UserData'')');
        end
        hPage = [hPage, h];
        h     = text(tCol(9),y,sprintf(TabDat.fmt{9},U),'FontWeight','Bold',...
            'UserData',U,'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];
        if 0
        h     = text(tCol(10),y,sprintf(TabDat.fmt{10},Ze(i)),'FontWeight','Bold',...
            'UserData',Ze,'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];
        end
        h     = ...
            text(tCol(11),y,sprintf(TabDat.fmt{11},Pz(i)),'FontWeight','Bold',...
            'UserData',Pz(i),'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];

        % Specifically changed so it properly finds hMIPax
        %------------------------------------------------------------------
        tXYZmm = XYZmm(DIM,i);
        h     = text(tCol(12),y,sprintf(TabDat.fmt{12},tXYZmm),...
            'FontWeight','Bold',...
            'Tag','ListXYZ',...
            'ButtonDownFcn',[...
            'hMIPax = findobj(''tag'',''hMIPax'');',...
            'spm_mip_ui(''SetCoords'',get(gcbo,''UserData''),hMIPax);'],...
            'Interruptible','off','BusyAction','Cancel',...
            'UserData',XYZmm(:,i));

        HlistXYZ = [HlistXYZ, h];
        if spm_XYZreg('Edist',xyzmm,XYZmm(:,i))<tol && ~isempty(hReg)
            set(h,'Color','r')
        end
        hPage  = [hPage, h];

        y      = y - dy;

        if topoFDR
            [TabDat.dat{TabLin,3:12}] = deal(Pk,Qc,N(i),Pn,Pu(i),Qp,U,Ze,Pz(i),XYZmm(:,i));
        else
            [TabDat.dat{TabLin,3:12}] = deal(Pk,Qc,N(i),Pn,Pu(i),Qu,U,Ze,Pz(i),XYZmm(:,i));
        end
        TabLin = TabLin + 1;

        %-Print Num secondary maxima (> Dis mm apart)
        %------------------------------------------------------------------
        [l q] = sort(-Z(j));                % sort on Z value
        D     = i;
        for i = 1:length(q)
            d = j(q(i));
            if min(sqrt(sum((XYZmm(:,D)-XYZmm(:,d)*ones(1,size(D,2))).^2)))>Dis;

                if length(D) < Num

                    % Paginate if necessary
                    %------------------------------------------------------
                    if y < dy
                        h = text(0.5,-6*dy,sprintf('Page %d',...
                            spm_figure('#page',Fgraph)),...
                            'FontName',PF.helvetica,...
                            'FontAngle','Italic',...
                            'FontSize',FS(8));

                        spm_figure('NewPage',[hPage,h])
                        hPage = [];
                        y     = y0;
                    end

                    % voxel-level p values {Z}
                    %------------------------------------------------------
                    if STAT ~= 'P'
                        Qp    = [];
                        if 0 
                            Ze    = spm_invNcdf(Z(d));
                        end
                    else
                        Pz    = [];
                        Pu    = [];
                        Qu    = [];
                        Qp    = [];
                        Ze    = spm_invNcdf(Z(d));
                    end

                    h     = text(tCol(7),y,sprintf(TabDat.fmt{7},Pu(d)),...
                        'UserData',Pu(d),...
                        'ButtonDownFcn','get(gcbo,''UserData'')');
                    hPage = [hPage, h];

                    if topoFDR
                        h     = text(tCol(8),y,sprintf(TabDat.fmt{8},Qp),...
                        'UserData',Qp,...
                        'ButtonDownFcn','get(gcbo,''UserData'')');
                    else
                        h     = text(tCol(8),y,sprintf(TabDat.fmt{8},Qu(d)),...
                        'UserData',Qu(d),...
                        'ButtonDownFcn','get(gcbo,''UserData'')');
                    end
                    hPage = [hPage, h];
                    h     = text(tCol(9),y,sprintf(TabDat.fmt{9},Z(d)),...
                        'UserData',Z(d),...
                        'ButtonDownFcn','get(gcbo,''UserData'')');
                    hPage = [hPage, h];
                    if 0
                        h     = text(tCol(10),y,sprintf(TabDat.fmt{10},Ze),...
                            'UserData',Ze,...
                            'ButtonDownFcn','get(gcbo,''UserData'')');
                    end
                    hPage = [hPage, h];
                    h     = text(tCol(11),y,sprintf(TabDat.fmt{11},Pz(d)),...
                        'UserData',Pz(d),...
                        'ButtonDownFcn','get(gcbo,''UserData'')');
                    hPage = [hPage, h];

                    % specifically modified line to use hMIPax
                    %------------------------------------------------------
                    tXYZmm = XYZmm(DIM,d);
                    h     = text(tCol(12),y,...
                        sprintf(TabDat.fmt{12},tXYZmm),...
                        'Tag','ListXYZ',...
                        'ButtonDownFcn',[...
                        'hMIPax = findobj(''tag'',''hMIPax'');',...
                        'spm_mip_ui(''SetCoords'',',...
                        'get(gcbo,''UserData''),hMIPax);'],...
                        'Interruptible','off','BusyAction','Cancel',...
                        'UserData',XYZmm(:,d));

                    HlistXYZ = [HlistXYZ, h];
                    if spm_XYZreg('Edist',xyzmm,XYZmm(:,d))<tol && ...
                            ~isempty(hReg)
                        set(h,'Color','r')
                    end
                    hPage = [hPage, h];
                    D     = [D d];
                    y     = y - dy;
                    if topoFDR
                        [TabDat.dat{TabLin,7:12}] = ...
                            deal(Pu(d),Qp,Z(d),Ze,Pz(d),XYZmm(:,d));
                    else
                        [TabDat.dat{TabLin,7:12}] = ...
                            deal(Pu(d),Qu,Z(d),Ze,Pz(d),XYZmm(:,d));
                    end
                    TabLin = TabLin+1;
                end
            end
        end
        Z(j) = NaN;     % Set local maxima to NaN
    end             % end region


    %-Number and register last page (if paginated)
    %-Changed to use Fgraph for numbering
    %----------------------------------------------------------------------
    if spm_figure('#page',Fgraph)>1
        h = text(0.5,-6*dy,sprintf('Page %d/%d',spm_figure('#page',Fgraph)*[1,1]),...
            'FontName',PF.helvetica,'FontSize',FS(8),'FontAngle','Italic');
        spm_figure('NewPage',[hPage,h])
    end

    %-End: Store TabDat in UserData of axes & reset pointer
    %======================================================================
    h = uicontextmenu('Tag','TabDat',...
        'UserData',TabDat);
    set(gca,'UIContextMenu',h,...
        'Visible','on',...
        'XColor','w','YColor','w')
    uimenu(h,'Label','Print text table',...
        'Tag','TD_TxtTab',...
        'CallBack',...
        'spm_list(''txtlist'',get(get(gcbo,''Parent''),''UserData''),3)',...
        'Interruptible','off','BusyAction','Cancel');
    uimenu(h,'Separator','off','Label','Extract table data structure',...
        'Tag','TD_Xdat',...
        'CallBack','TabDat=get(get(gcbo,''Parent''),''UserData'')',...
        'Interruptible','off','BusyAction','Cancel');
    if ispc
        uimenu(h,'Separator','off','Label','Export to Excel',...
        'Tag','TD_Xdat',...
        'CallBack',@export2excel,...
        'Interruptible','off','BusyAction','Cancel');
    end
    uimenu(h,'Separator','on','Label','Help',...
        'Tag','TD_Xdat',...
        'CallBack','spm_help(''spm_list'')',...
        'Interruptible','off','BusyAction','Cancel');

    %-Setup registry
    %----------------------------------------------------------------------
    set(hAx,'UserData',struct('hReg',hReg,'HlistXYZ',HlistXYZ))
    spm_XYZreg('Add2Reg',hReg,hAx,'cg_tfce_list');

    %-Return TabDat structure & reset pointer
    %----------------------------------------------------------------------
    varargout = {TabDat};
    spm('Pointer','Arrow')


    %======================================================================
    case 'listcluster'                      %-List for current cluster only
    %======================================================================
    % FORMAT TabDat = cg_tfce_list('listcluster',SPM,hReg)

        spm('Pointer','Watch')

        %-Parse arguments
        %------------------------------------------------------------------
        if nargin < 2,  error('insufficient arguments'),     end
        if nargin < 3,  hReg = []; else hReg = varargin{3}; end
        SPM    = varargin{2};

        %-get number and separation for maxima to be reported
        %------------------------------------------------------------------
        if length(varargin) > 3

            Num    = varargin{4};       % number of maxima per cluster
            Dis    = varargin{5};       % distance among clusters (mm)
        else
            Num    = 32;
            Dis    = 4;
        end


        %-if there are suprathreshold voxels, filter out all but current cluster
        %------------------------------------------------------------------
        if ~isempty(SPM.Z)

            %-Jump to voxel nearest current location
            %--------------------------------------------------------------
            [xyzmm,i] = spm_XYZreg('NearestXYZ',...
                spm_results_ui('GetCoords'),SPM.XYZmm);
            spm_results_ui('SetCoords',SPM.XYZmm(:,i));

            %-Find selected cluster
            %--------------------------------------------------------------
            A         = spm_clusters(SPM.XYZ);
            j         = find(A == A(i));
            SPM.Z     = SPM.Z(j);
            SPM.Qu    = SPM.Qu(j);
            SPM.Pu    = SPM.Pu(j);
            SPM.Pz    = SPM.Pz(j);
            SPM.XYZ   = SPM.XYZ(:,j);
            SPM.XYZmm = SPM.XYZmm(:,j);
            if isfield(SPM,'Rd'), SPM.Rd = SPM.Rd(:,j); end
        end

        %-Call 'list' functionality to produce table
        %------------------------------------------------------------------
        varargout = {cg_tfce_list('list',SPM,hReg,Num,Dis)};


    %======================================================================
    case 'txtlist'                                 %-Print ASCII text table
    %======================================================================
    % FORMAT cg_tfce_list('TxtList',TabDat,c)

        if nargin<2, error('Insufficient arguments'), end
        if nargin<3, c=1; else c=varargin{3}; end
        TabDat = varargin{2};

        %-Table Title
        %------------------------------------------------------------------
        fprintf('\n\nSTATISTICS: %s\n',TabDat.tit)
        fprintf('%c',repmat('=',1,80)), fprintf('\n')

        %-Table header
        %------------------------------------------------------------------
        fprintf('%s\t',TabDat.hdr{1,c:end-1}), fprintf('%s\n',TabDat.hdr{1,end})
        fprintf('%s\t',TabDat.hdr{2,c:end-1}), fprintf('%s\n',TabDat.hdr{2,end})
        fprintf('%c',repmat('-',1,80)), fprintf('\n')

        %-Table data
        %------------------------------------------------------------------
        for i = 1:size(TabDat.dat,1)
            for j=c:size(TabDat.dat,2)
                fprintf(TabDat.fmt{j},TabDat.dat{i,j})
                fprintf('\t')
            end
            fprintf('\n')
        end
        for i=1:max(1,12-size(TabDat.dat,1)), fprintf('\n'), end
        fprintf('%s\n',TabDat.str)
        fprintf('%c',repmat('-',1,80)), fprintf('\n')

        %-Table footer
        %------------------------------------------------------------------
        for i=1:size(TabDat.ftr,1)
            fprintf([TabDat.ftr{i,1} '\n'],TabDat.ftr{i,2});
        end
        fprintf('%c',repmat('=',1,80)), fprintf('\n\n')



        %==================================================================
    case 'setcoords'                                   %-Co-ordinate change
        %==================================================================
        % FORMAT cg_tfce_list('SetCoords',xyz,hAx,hReg)
        if nargin<3, error('Insufficient arguments'), end
        hAx      = varargin{3};
        xyz      = varargin{2};
        UD       = get(hAx,'UserData');
        HlistXYZ = UD.HlistXYZ(ishandle(UD.HlistXYZ));

        %-Set all co-ord strings to black
        %------------------------------------------------------------------
        set(HlistXYZ,'Color','k')

        %-If co-ord matches a string, highlight it in red
        %------------------------------------------------------------------
        XYZ      = get(HlistXYZ,'UserData');
        if iscell(XYZ), XYZ = cat(2,XYZ{:}); end
        [null,i,d] = spm_XYZreg('NearestXYZ',xyz,XYZ);
        if d<eps
            set(HlistXYZ(i),'Color','r')
        end

        %==================================================================
    otherwise                                       %-Unknown action string
        %==================================================================
        error('Unknown action string')
end
%==========================================================================

%==========================================================================
function export2excel(obj,evd,h)
TabDat     = get(get(obj,'Parent'),'UserData');
d          = [TabDat.hdr;TabDat.dat];
xyz        = d(3:end,end);
xyz        = num2cell([xyz{:}]');
d(:,end+1) = d(:,end);
d(:,end+1) = d(:,end);
d(3:end,end-2:end) = xyz;
tmpfile    = [tempname '.xls'];
xlswrite(tmpfile, d);
winopen(tmpfile);
