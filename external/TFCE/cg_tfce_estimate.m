function cg_tfce_estimate(job)

global Y initialized

initialized = 0;

% use debug for displaying permuted design matrix
debug = 1;

% allow to test permutations with analyzing data
test_mode = 0;

% define stepsize for tfce
n_steps_tfce = 100;
    
% define histogram bins
% use value > 1000 to reliable estimate p<0.001 levels
n_hist_bins = 1100;
    
% colors and alpha levels
col   = [0 0 1; 0 0.5 0; 1 0 0];
alpha = [0.05   0.01     0.001];
    
% give same results each time
rand('state',0);

% tolerance for comparing real numbers
tol = 1e-4;
    
% check spm version
if strcmp(spm('ver'),'SPM12')
  spm12 = 1;
else
  spm12 = 0;
end

load(job.spmmat{1});
cwd = fileparts(job.spmmat{1});

%-Check that model has been estimated
if ~isfield(SPM, 'xVol')
    str = { 'This model has not been estimated.';...
              'Would you like to estimate it now?'};
    if spm_input(str,1,'bd','yes|no',[1,0],1)
        cd(cwd)
        SPM = spm_spm(SPM);
    else
        return
    end
end
    
Ic0 = job.conspec.contrasts;
try
  xCon0 = SPM.xCon(Ic0(1));
catch
  [Ic0,xCon] = spm_conman(SPM,'T',Inf,...
        '  Select contrast(s)...',' ',1);
  SPM.xCon = xCon;
end

% check that no temporal filter was used
if isstruct(SPM.xX.K)
    fprintf('ERROR: No first level analysis with temporal correlations allowed.\n');
    return
end
    
% whitening matrix
if isfield(SPM.xX,'W')
    W = SPM.xX.W;
    if any(W~=1)
        fprintf('Whitening of the data is not yet supported.\n');
        W = speye(size(SPM.xX.X,1));
    end
else
    W = speye(size(SPM.xX.X,1));
end
    
% get some parameters from SPM
xX     = SPM.xX;
TH     = SPM.xM.TH;
VY     = SPM.xY.VY;
n_data = size(xX.X,1);

% sometimes xX.iB and xX.iH are not correct and cannot be used to reliably recognize the design
xX = correct_xX(xX);

% find exchangeability block labels for longitudinal designs (paired t-test, flexible factorial)
repeated_anova = ~isempty(xX.iB);
if repeated_anova
    [rw,cl] = find(xX.I == length(xX.iB)); % find column which codes subject factor (length(xX.iB) -> n_subj)
    exch_block_labels = xX.I(:,cl(1));      % column from above contains the subject factor

    % check that labels are defined for each block
    for i=1:n_data
        groupListed(exch_block_labels(i)) = true;
    end
    for i = 1: max(exch_block_labels)
        if ~groupListed(i) 
            fprintf('Error: block %d must be assigned to at least one design row in the blocks file.\n',i);
            return
        end
    end
    fprintf('Please note that permutation is only done within subjects for repeated Anova.\n',i);
else
    exch_block_labels = ones(1,n_data);
end

% go through all contrasts
for con = 1:length(Ic0)

    Ic = Ic0(con);
    xCon = SPM.xCon(Ic);
    
    if ~strcmp(xCon.STAT,'T')
      fprintf('ERROR: Only T-contrasts allowed.\n');
      return
    end
    
    n_perm = job.conspec.n_perm(1);
    if numel(job.conspec.n_perm) > 1
      n_perm_break = job.conspec.n_perm(2);
    end
    
    vFWHM = job.conspec.vFWHM;
    
    if length(Ic) > 1
      fprintf('ERROR: No conjunction allowed.\n');
      return
    end
    
    % get contrast and name
    c = xCon.c;
    ind_con = find(c~=0)';
    c_name  = deblank(xCon.name);

    % find exchangeability blocks using contrasts without zero values
    [exch_blocks, I, J]   = unique(c(ind_con),'first');
    n_exch_blocks = length(exch_blocks);
    
    % check for exchangeability blocks and design matrix
    if n_exch_blocks == 1
        n_cond = length(find(xX.iH==ind_con)); % check whether the contrast is defined at columns for condition effects
    else
        n_cond = 0;
        n_data_cond = [];
        for k=1:length(xX.iH)
          n_data_cond = [n_data_cond sum(xX.X(:,xX.iH(k)))];
        end
        for j=1:n_exch_blocks
          c_exch_blocks = find(c==exch_blocks(j));
          for k=1:length(c_exch_blocks)
            n_cond = n_cond + length(find(xX.iH==c_exch_blocks(k)));
          end
        end       
    end

    use_half_permutations = 0;
    % check if sample size is equal for both conditions
    if n_cond == 2
        if sum(n_data_cond(find(c==exch_blocks(1)))) == sum(n_data_cond(find(c==exch_blocks(2))))
            use_half_permutations = 1;
            fprintf('Equal sample sizes: half of permutations are used.\n');
        end
    end

    % repated Anova does not allow to use only half of the permutions
    if repeated_anova
      use_half_permutations = 0;
    end
    
    ind_exch_blocks = cell(n_exch_blocks,1);
    for j=1:n_exch_blocks
      ind_exch_blocks{j} = find(c==exch_blocks(j));
    end
        
    % check design
    switch n_cond
    case 0 % correlation
      if n_exch_blocks == 2
        warning('Interaction design between two regressors found. This should work, but is not yet fully tested and I am unsure whether the # of max. permutations is correctly estimated.')
      else
        if repeated_anova
          fprintf('Repeated Anova with contrast for covariate found\n');
        else
          fprintf('Multiple regression design found\n');
        end
      end
      label = 1:n_data;
    case 1 % one-sample t-test
      fprintf('One sample t-test found\n');
      % use exchangeability blocks for labels
      label = zeros(1,n_data);
      for j=1:n_exch_blocks
        for k=1:length(ind_exch_blocks{j})
          label(find(xX.X(:,ind_exch_blocks{j}(k))==1)) = j;
        end
      end
    otherwise  % Anova with at least 2 groups
      if repeated_anova
        fprintf('Repeated Anova found\n');
      else
        fprintf('Anova found\n');
      end
      % use exchangeability blocks for labels
      label = zeros(1,n_data);
      for j=1:n_exch_blocks
        for k=1:length(ind_exch_blocks{j})
          label(find(xX.X(:,ind_exch_blocks{j}(k))==1)) = j;
        end
      end
    end
    
    fprintf('\n')
    
    % get index for label values > 0
    ind_label = find(label > 0);
    
    n_data_with_contrast = length(find(label > 0));
    
    % estimate # of permutations
    % Anova/correlation: n_perm = (n1+n2+...+nk)!/(n1!*n2!*...*nk!)
    if n_cond ~=1  % Anova/correlation
      n_perm_full = factorial(n_data_with_contrast);
      single_subject = 0;
      
      for i=1:n_cond
        % check whether only a single subject is in one group
        if length(find(label == i)) == 1
          single_subject = 1;
        end
        n_perm_full = n_perm_full/factorial(length(find(label == i)));
      end
      
      if isnan(n_perm_full)
        % correct number of permutations for large samples when factorial is not working
        if (n_cond == 2) & (single_subject == 1)
          n_perm_full = n_data_with_contrast;
        else
          n_perm_full = realmax;
        end
      end
      
      % find where data are defined for that contrast
      if ~isempty(find(xX.iH == ind_con(1)))
        % first checking whether contrasts are defined for iH
        ind_data_defined = find(any(xX.X(:,xX.iH(ind_con)),2));
      else
        % if not check iC
        if length(ind_con) == 1 % for one contrast covariate column use all subjects
          ind_data_defined = 1:size(xX.X,1);
        else % check wehere contrast is defined
          ind_data_defined = find(any(xX.X(:,ind_con),2));
        end
      end
      
      % and restrict exchangeability block labels to those rows
      exch_block_labels_new = exch_block_labels(ind_data_defined);

      % Repated Anova: n_perm = n_cond1!*n_cond2!*...*n_condk!
      % for a full model where each condition is defined fro all subjects the easier
      % estimation is: n_perm = (n_cond!)^n_subj
      % check that no regression analysis inside repeated anova is used
      if repeated_anova & n_cond~=0
        n_subj = max(exch_block_labels_new);
        n_perm_full = 1;
        for k=1:n_subj
          n_cond_subj = length(find(exch_block_labels_new == k));
          n_perm_full = n_perm_full*factorial(n_cond_subj);
        end
      else
        n_perm_full = round(n_perm_full);
      end
      
    else  % one-sample t-test: n_perm = 2^n
      n_perm_full = 2^n_data_with_contrast;
      exch_block_labels_new = exch_block_labels;
    end

    if debug
      fprintf('# full permutations: %d\n',n_perm_full);
      if use_half_permutations
        fprintf('Equal sample size found: Use half permutations.\n');
      end
      fprintf('exchangeability blocks: ');
      for j=1:n_exch_blocks
        fprintf('%d ',ind_exch_blocks{j});
      end
      fprintf('\n');
      fprintf('# of conditions: %d\n',n_cond);
    end
    
    n_perm = min([n_perm n_perm_full]);
        
    % check for mask image that should exist for any analysis
    if exist(fullfile(cwd, 'mask.img'))
      file_ext = '.img';
    elseif exist(fullfile(cwd, 'mask.nii'))
      file_ext = '.nii';
    elseif exist(fullfile(cwd, 'mask.gii'))
      file_ext = '.gii';
    else
      fprintf('ERROR: No mask file found.\n');
      return
    end

    % check for meshes
    if spm12
      if spm_mesh_detect(VY)
        mesh_detected = 1;
        if vFWHM > 0
          warning('Variance smoothing not yet prepared!');
        end
      else
        mesh_detected = 0;
      end
    else
      mesh_detected = 0;
    end
    
    if ~test_mode
      % get mask file
      if isempty(job.mask)
        maskname = fullfile(cwd,['mask' file_ext]);
      else
        maskname = job.mask{1};
      end
    
      % load mask
      try
        if spm12
          Vmask = spm_data_hdr_read(maskname);
        else
          Vmask = spm_vol(maskname);
        end
      catch
        if mesh_detected
          maskname = spm_select(1,'mesh','select mask image');
        else
          maskname = spm_select(1,'image','select mask image');
        end
        if spm12
          Vmask = spm_data_hdr_read(maskname);
        else
          Vmask = spm_vol(maskname);
        end
      end
    
      % if first image was not found you have to select all files again
      if ~exist(VY(1).fname);
        n = size(SPM.xY.VY,1);
        if mesh_detected
          P = spm_select(size(SPM.xY.VY,1),'mesh','select images');
        else
          P = spm_select(size(SPM.xY.VY,1),'image','select images');
        end
        if spm12
          VY = spm_data_hdr_read(P);
        else
          VY = spm_vol(P);
        end
        SPM.xY.VY = VY;
      
        % update SPM
        if size(SPM.xY.VY,1)==n
          save(job.spmmat{1},'SPM');
        else
          fprintf('ERROR: Number of files is not correct\n');
          return
        end
      end
        
      % check whether mask images fits to the data
      if mesh_detected, dim_index = 1; else dim_index=1:3; end
      if sum(sum((Vmask.mat-VY(1).mat).^2)) > 1e-6 || any(Vmask.dim(dim_index) ~= VY(1).dim(dim_index))
        error('Mask must have the same dimensions and orientation as the data.');
      end

      % compute unpermuted t-map
      t0 = calc_glm(VY,xX.X,xCon.c,Vmask,vFWHM,TH,W,job.openmp);
    
      % get deltaT for unpermuted map
      deltaT = max(abs(t0(:)))/n_steps_tfce;
    
      % calculate tfce of unpermuted t-map
      if mesh_detected
        tfce0 = tfce_mesh(SPM.xVol.G.faces, t0, deltaT);
      else
        if job.openmp
          tfce0 = tfceMex(t0, deltaT, job.tbss, 1);
        else
          tfce0 = tfceMex_noopenmp(t0, deltaT, job.tbss);
        end
      end
    
      % get largest tfce
      tfce0_max = max(tfce0(:));
      t0_max    = max(t0(:));
    
      % get vector for histogram bins
      tfce_bins = linspace(0, max(abs(tfce0(:))), n_hist_bins);
      t_bins    = linspace(0, max(abs(t0(:))), n_hist_bins);
    
      % prepare countings
      t_hist       = zeros(1, n_hist_bins);
      tfce_hist    = zeros(1, n_hist_bins);
      t_max        = [];
      t_max_th     = [];
      t_th         = [];
      tfce_max     = [];
      tfce_max_th  = [];
      tfce_th      = [];
    
    end % test_mode
    
    label_matrix = [];
    
    % general initialization
    try % use try commands to allow batch mode without graphical output
      Fgraph = spm_figure('GetWin','Graphics');
      spm_figure('Clear',Fgraph);
      figure(Fgraph)
    
      h = axes('position',[0.45 0.95 0.1 0.05],'Units','normalized','Parent',...
        Fgraph,'Visible','off');
        
      text(0.5,0.6,c_name,...
        'FontSize',spm('FontSize',10),...
        'FontWeight','Bold',...
        'HorizontalAlignment','Center',...
        'VerticalAlignment','middle')
    
      text(0.5,0.25,spm_str_manip(SPM.swd,'a50'),...
        'FontSize',spm('FontSize',8),...
        'HorizontalAlignment','Center',...
        'VerticalAlignment','middle')
    end
    
    % check that label has correct dimension
    sz = size(label);
    if sz(1)>sz(2)
      label = label';
    end
    
    stopStatus = false;
    try, spm_progress_bar('Init',n_perm,'Calculating','Permutations'); end
    cg_progress('Init',n_perm,'Calculating','Permutations')
    
    % update interval for progress bar
    progress_step = max([1 round(n_perm/1000)]);
    
    i = 1;
    
    while i<=n_perm
    
      % randomize subject vector
      if i==1 % first permutation is always unpermuted model
        if n_cond == 1 % one-sample t-test
          rand_label = ones(1,n_data_with_contrast);
        else % correlation or Anova
          rand_order = 1:n_data_with_contrast;
          rand_label = label(ind_label(rand_order));
        end
      else
        % init permutation and
        % check that each permutation is used only once
        if n_cond == 1 % one-sample t-test
          rand_label = sign(randn(1,n_data_with_contrast));
          while any(ismember(label_matrix,rand_label,'rows'))
            rand_label = sign(randn(1,n_data_with_contrast));
          end
        else % correlation or Anova
          
          % permute inside exchangeability blocks only
          rand_order = ones(1,n_data_with_contrast);
          rand_label = ones(1,n_data_with_contrast);
          for k = 1:max(exch_block_labels_new)
            ind_block   = find(exch_block_labels_new == k);
            n_per_block = length(ind_block);
            rand_order(ind_block) = ind_block(randperm(n_per_block));
            rand_label(ind_block) = label(ind_label(rand_order(ind_block)));
          end

          % check whether this permutation was already used
          while any(ismember(label_matrix,rand_label,'rows'))
            for k = 1:max(exch_block_labels_new)
              ind_block   = find(exch_block_labels_new == k);
              n_per_block = length(ind_block);
              rand_order(ind_block) = ind_block(randperm(n_per_block));
              rand_label(ind_block) = label(ind_label(rand_order(ind_block)));
            end
          end
        end    
      end   
      
      % add Stop button after 20 iterations
      try % use try commands to allow batch mode without graphical output
        if i==21
          hStopButton = uicontrol(Fgraph,...
            'position',[10 10 70 20],...
            'style','toggle',...
            'string','Stop',...
            'backgroundcolor',[1 .5 .5]); % light-red
        end
      
        if i>=21
          stopStatus = get(hStopButton,'value');
        end
      
        % check Stop status
      if (stopStatus == true)
          fprintf('Stopped after %d iterations.\n',i);
          break; % stop the permutation loop
        end
      end
        
      % change design matrix according to permutation order
      % only permute columns, where contrast is defined
      Xperm = xX.X;
      Xperm_debug = xX.X;
      
      if debug
        % scale covariates and nuisance variables to a range 0.8..1
        % to properly display these variables with indicated colors
        if ~isempty(xX.iC)
          val = Xperm_debug(:,xX.iC);
          mn = repmat(min(val),length(val),1); mx = repmat(max(val),length(val),1);
          val = 0.8 + 0.2*(val-mn)./(mx-mn);
          Xperm_debug(:,xX.iC) = val;
        end
        if ~isempty(xX.iG)
          val = Xperm_debug(:,xX.iG);
          mn = repmat(min(val),length(val),1); mx = repmat(max(val),length(val),1);
          val = 0.8 + 0.2*(val-mn)./(mx-mn);
          Xperm_debug(:,xX.iG) = val;
        end
        
        % use different colors for indicated columns
        Xperm_debug(:,xX.iH) = 16*Xperm_debug(:,xX.iH);
        Xperm_debug(:,xX.iC) = 24*Xperm_debug(:,xX.iC);
        Xperm_debug(:,xX.iB) = 32*Xperm_debug(:,xX.iB);
        Xperm_debug(:,xX.iG) = 48*Xperm_debug(:,xX.iG);
      end
      
      if n_cond==1 % one-sample t-test
        % only change sign in the design matrix
        for j=1:n_data_with_contrast
          Xperm(ind_label(j),ind_con) = rand_label(j)*Xperm(ind_label(j),ind_con);
          Xperm_debug(ind_label(j),ind_con) = Xperm(ind_label(j),ind_con);
          
          if rand_label(j) > 0
            Xperm_debug(ind_label(j),ind_con) = 60*rand_label(j)*Xperm_debug(ind_label(j),ind_con);
          else
            Xperm_debug(ind_label(j),ind_con) = 56*rand_label(j)*Xperm_debug(ind_label(j),ind_con);
          end
        end
      else % correlation or Anova
        Xperm(ind_label,ind_con) = Xperm(ind_label(rand_order),ind_con);
        Xperm_debug(ind_label,ind_con) = Xperm(ind_label,ind_con);
        
        % scale exchangeability blocks also to values 0.8..1
        val = Xperm_debug(:,ind_con);
        ind0 = find(val==0);
        mn = repmat(min(val),length(val),1); mx = repmat(max(val),length(val),1);
        val = 0.8 + 0.2*(val-mn)./(mx-mn);
        
        % rescue zero entries
        val(ind0) = 0;
        
        Xperm_debug(:,ind_con) = 60*val;
      end
      
      % display permuted design matrix
      try
        if debug
          subplot(2,2,3);
          image(Xperm_debug); axis off
          title('Permuted design matrix','FontWeight','bold');
          drawnow
        
          % use different colormap for permuted design matrix
          cmap = jet(64);
        
          % zero values should be always black
          cmap(1,:) = [0 0 0];
          colormap(cmap)
        
          subplot(2,2,4); axis off
        
          % color-coded legend
          y = 1.0;
          text(-0.2,y, 'Columns of design matrix: ', 'Color',cmap(1, :),'FontWeight','Bold','FontSize',10); y = y - 0.10;
          text(-0.2,y,['Exchangeability blocks: ' num2str_short(sort(cell2mat(ind_exch_blocks))')], 'Color',cmap(60,:),'FontWeight','Bold','FontSize',10); y = y - 0.05;
          if ~isempty(xX.iH)
            text(-0.2,y, ['iH - Indicator variables: ' num2str_short(xX.iH)], 'Color',cmap(16,:),'FontWeight','Bold','FontSize',10);
            y = y - 0.05; 
          end
          if ~isempty(xX.iC)
            text(-0.2,y, ['iC - Covariates: ' num2str_short(xX.iC)], 'Color',cmap(24,:),'FontWeight','Bold','FontSize',10);
            y = y - 0.05;
          end
          if ~isempty(xX.iB)
            text(-0.2,y, ['iB - Block effects: ' num2str_short(xX.iB)], 'Color',cmap(32,:),'FontWeight','Bold','FontSize',10);
            y = y - 0.05;
          end
          if ~isempty(xX.iG)
            text(-0.2,y, ['iG - Nuisance variables: ' num2str_short(xX.iG)], 'Color',cmap(48,:),'FontWeight','Bold','FontSize',10);
            y = y - 0.05;
          end
        end
      end
      
      if ~test_mode
        % calculate permuted t-map
        if i==1
          t = t0;
          tfce = tfce0;
        else
          % -----------------------------------------------------
          % -----------------------------------------------------
          % What about W for whitening the data??? Should this also permuted???
          % -----------------------------------------------------
          % -----------------------------------------------------
        
          t = calc_glm(VY,Xperm,xCon.c,Vmask,vFWHM,TH,W,job.openmp);
        
          % compute tfce
          if mesh_detected
            tfce = tfce_mesh(SPM.xVol.G.faces, t, deltaT);
          else
            if job.openmp
              tfce = tfceMex(t, deltaT, job.tbss);
            else
              tfce = tfceMex_noopenmp(t0, deltaT, job.tbss);
            end
          end
        end
      end % test_mode
      
      % update label_matrix and order_matrix for checking of unique permutations
      if use_half_permutations
        label_matrix = [label_matrix; rand_label; 3-rand_label];
    
        if ~test_mode
          % maximum statistic
          tfce_max = [tfce_max max(tfce(:)) -min(tfce(:))];
          t_max    = [t_max max(t(:)) -min(t(:))];
        end
      else
        label_matrix = [label_matrix; rand_label];

        if ~test_mode
          % maximum statistic
          tfce_max = [tfce_max max(tfce(:))];
          t_max    = [t_max max(t(:))];
        end
      end
        
      if ~test_mode
        % cummulate histogram
        tfce_gt0 = tfce(find(tfce>0));
        if ~isempty(tfce_gt0)
          tfce_hist = tfce_hist + hist(tfce_gt0, tfce_bins);
        end
        t_gt0 = t(find(t>0));
        if ~isempty(t_gt0)
          t_hist = t_hist + hist(t_gt0, t_bins);
        end
      
        if use_half_permutations
          tfce_lt0 = tfce(find(tfce<0));
          if ~isempty(tfce_lt0)
            tfce_hist = tfce_hist + hist(-tfce_lt0, tfce_bins);
          end
          t_lt0 = t(find(t<0));
          if ~isempty(t_lt0)
            t_hist = t_hist + hist(-t_lt0, t_bins);
          end
        end
      
        % use cummulated sum to find threshold
        tfce_max = sort(tfce_max);
        t_max    = sort(t_max);
    
        % find corrected thresholds
        ind_max  = ceil((1-alpha).*length(t_max));
        t_max_th = [t_max_th; t_max(ind_max);];
        if use_half_permutations
          t_max_th = [t_max_th; t_max(ind_max);];
        end
        
        ind_max     = ceil((1-alpha).*length(tfce_max));
        tfce_max_th = [tfce_max_th; tfce_max(ind_max);];
        if use_half_permutations
          tfce_max_th = [tfce_max_th; tfce_max(ind_max);];
        end
    
        % plot thresholds and histograms
        try
          figure(Fgraph)
          h1 = axes('position',[0 0 1 0.95],'Parent',Fgraph,'Visible','off');
          plot_distribution(tfce_max, tfce_max_th, 'tfce', alpha, col, 1, tfce0_max);
          if ~debug
            plot_distribution(t_max, t_max_th, 't-value', alpha, col, 2, t0_max);
          end
          drawnow
        end
      
        if numel(job.conspec.n_perm) > 1
          if i > n_perm_break
            if isempty(find(tfce0_max > tfce_max_th(50:end,1)))
              fprintf('No FWE-corrected suprathreshold value after %d permutations found\n', n_perm_break);
              i = n_perm;
            end
          end  
        end
      end % test_mode

      
      if use_half_permutations
        if ~rem(i,progress_step) || ~rem(i+1,progress_step)
          cg_progress('Set',i,Fgraph)
          try, spm_progress_bar('Set',i); end
        end
      else
        if ~rem(i,progress_step)
          cg_progress('Set',i,Fgraph)
          try, spm_progress_bar('Set',i); end
        end
      end
    
      if use_half_permutations  
        i = i + 2;
      else
        i = i + 1;
      end
    
    end
    
    cg_progress('Clear',Fgraph)
    try
      spm_progress_bar('Clear');
      delete(hStopButton);
      spm_print;
    end
    
    if ~test_mode
      % get correct number of permutations in case that process was stopped
      n_perm = length(tfce_max);

      %---------------------------------------------------------------
      % corrected threshold based on permutation distribution
      %---------------------------------------------------------------
    
      % allow thresholds depending on # of permutations
      n_alpha = 3;
      if n_perm < 1000, n_alpha = 2; end
      if n_perm <  100, n_alpha = 1; end
    
      % pepare output files
      Vt = VY(1);
      Vt.dt(1)    = 16;
      Vt.pinfo(1) = 1;
        
      mask_0 = find(t0 == 0);
      mask_P = find(t0 > 0);
      mask_N = find(t0 < 0);

      %---------------------------------------------------------------
      % save unpermuted TFCE map
      %---------------------------------------------------------------
      name = sprintf('TFCE_%04d',Ic);
      Vt.fname = fullfile(cwd,[name file_ext]);
      if vFWHM > 0
        Vt.descrip = sprintf('TFCE FWHM=%.1fmm, Contrast %04d',vFWHM,Ic);
      else
        Vt.descrip = sprintf('TFCE Contrast %04d',Ic);
      end
      if spm12
        Vt = spm_data_hdr_write(Vt);
        spm_data_write(Vt,tfce0);
      else
        spm_write_vol(Vt,tfce0);
      end
    
      %---------------------------------------------------------------
      % save unpermuted T map
      %---------------------------------------------------------------
      name = sprintf('T_%04d',Ic);
      Vt.fname = fullfile(cwd,[name file_ext]);
      if vFWHM > 0
        Vt.descrip = sprintf('T FWHM=%.1fmm, Contrast %04d',vFWHM,Ic);
      else
        Vt.descrip = sprintf('T Contrast %04d',Ic);
      end
      if spm12
        Vt = spm_data_hdr_write(Vt);
        spm_data_write(Vt,t0);
      else
        spm_write_vol(Vt,t0);
      end
    
      % save ascii file with number of permutations
      fid = fopen(fullfile(cwd,[name '.txt']),'w');
      fprintf(fid,'%d\n',n_perm);
      fclose(fid);
    
      %---------------------------------------------------------------
      % save uncorrected p-values for TFCE
      %---------------------------------------------------------------
      fprintf('Save uncorrected p-values.\n');

      name = sprintf('TFCE_log_p_%04d',Ic);
      Vt.fname = fullfile(cwd,[name file_ext]);
      if vFWHM > 0
        Vt.descrip = sprintf('TFCE FWHM=%.1fmm, Contrast %04d',vFWHM,Ic);
      else
        Vt.descrip = sprintf('TFCE Contrast %04d',Ic);
      end
    
      nPtfce = zeros(size(tfce0));
    
      % estimate p-values
      tfce_cumsum = cumsum(tfce_hist);
      for j=n_hist_bins:-1:1
        tmp = min(find(tfce_cumsum>=ceil(j/n_hist_bins*sum(tfce_hist))));
        nPtfce((tfce0  > tfce_bins(tmp)) & (nPtfce==0)) =  j/n_hist_bins;
        nPtfce((-tfce0 > tfce_bins(tmp)) & (nPtfce==0)) = -j/n_hist_bins;
      end
    
      nPtfce(mask_P) = 1 - nPtfce(mask_P);
      nPtfce(mask_N) = 1 + nPtfce(mask_N);
      nPtfce(mask_0) = NaN;

      nPtfcelog10 = zeros(size(tfce0));
      nPtfcelog10(mask_P) = -log10(nPtfce(mask_P));
      nPtfcelog10(mask_N) =  log10(nPtfce(mask_N));
      nPtfcelog10(mask_0) = NaN;

      if spm12
        Vt = spm_data_hdr_write(Vt);
        spm_data_write(Vt,nPtfcelog10);
      else
        spm_write_vol(Vt,nPtfcelog10);
      end
    
      %---------------------------------------------------------------
      % save uncorrected p-values for T
      %---------------------------------------------------------------
      name = sprintf('T_log_p_%04d',Ic);
      Vt.fname = fullfile(cwd,[name file_ext]);
      if vFWHM > 0
        Vt.descrip = sprintf('T FWHM=%.1fmm, Contrast %04d',vFWHM,Ic);
      else
        Vt.descrip = sprintf('T Contrast %04d',Ic);
      end
    
      nPt = zeros(size(t0));
      nPtlog10 = zeros(size(t0));
    
      % estimate p-values
      t_cumsum = cumsum(t_hist);
      for j=n_hist_bins:-1:1
        tmp = min(find(t_cumsum>=ceil(j/n_hist_bins*sum(t_hist))));
        nPt((t0  > t_bins(tmp)) & (nPt==0)) =  j/n_hist_bins;
        nPt((-t0 > t_bins(tmp)) & (nPt==0)) = -j/n_hist_bins;
      end
     
      nPt(mask_P) = 1 - nPt(mask_P);
      nPt(mask_N) = 1 + nPt(mask_N);
      nPt(mask_0) = NaN;

      nPtlog10(mask_P) = -log10(nPt(mask_P));
      nPtlog10(mask_N) =  log10(nPt(mask_N));
      nPtlog10(mask_0) = NaN;

      if spm12
        Vt = spm_data_hdr_write(Vt);
        spm_data_write(Vt,nPtlog10);
      else
        spm_write_vol(Vt,nPtlog10);
      end
    
      %---------------------------------------------------------------
      % save corrected p-values for TFCE
      %---------------------------------------------------------------
      fprintf('Save corrected p-values.\n');

      name = sprintf('TFCE_log_pFWE_%04d',Ic);
      Vt.fname = fullfile(cwd,[name file_ext]);
      if vFWHM > 0
        Vt.descrip = sprintf('TFCE FWHM=%.1fmm, Contrast %04d',vFWHM,Ic);
      else
        Vt.descrip = sprintf('TFCE Contrast %04d',Ic);
      end

      corrP = zeros(size(t));
    
      for t2 = tfce_max
        %-FWE-corrected p is proportion of randomisation greater or
        % equal to statistic.
        %-Use a > b -tol rather than a >= b to avoid comparing
        % two reals for equality.
        corrP = corrP + (t2 > tfce0 - tol);
        corrP = corrP - (t2 > -tfce0 - tol);
      end
      corrP = corrP / n_perm;  

      corrP(mask_P) = 1 + corrP(mask_P);
      corrP(mask_N) = 1 - corrP(mask_N);
      corrP(mask_0) = NaN;

      corrPlog10 = zeros(size(tfce0));
      corrPlog10(mask_P) = -log10(corrP(mask_P));
      corrPlog10(mask_N) =  log10(corrP(mask_N));
      corrPlog10(mask_0) = NaN;

      if spm12
        Vt = spm_data_hdr_write(Vt);
        spm_data_write(Vt,corrPlog10);
      else
        spm_write_vol(Vt,corrPlog10);
      end
    
      %---------------------------------------------------------------
      % save corrected p-values for T
      %---------------------------------------------------------------
      name = sprintf('T_log_pFWE_%04d',Ic);
      Vt.fname = fullfile(cwd,[name file_ext]);
      if vFWHM > 0
        Vt.descrip = sprintf('T FWHM=%.1fmm, Contrast %04d',vFWHM,Ic);
      else
        Vt.descrip = sprintf('T Contrast %04d',Ic);
      end

      corrP = zeros(size(t));
    
      for t2 = t_max
        %-FWE-corrected p is proportion of randomisation greater or
        % equal to statistic.
        %-Use a > b -tol rather than a >= b to avoid comparing
        % two reals for equality.
        corrP = corrP + (t2 > t0 - tol);
        corrP = corrP - (t2 > -t0 - tol);
      end
      corrP = corrP / n_perm;  

      corrP(mask_P) = 1 + corrP(mask_P);
      corrP(mask_N) = 1 - corrP(mask_N);
      corrP(mask_0) = NaN;

      corrPlog10 = zeros(size(tfce0));
      corrPlog10(mask_P) = -log10(corrP(mask_P));
      corrPlog10(mask_N) =  log10(corrP(mask_N));
      corrPlog10(mask_0) = NaN;
    
      if spm12
        Vt = spm_data_hdr_write(Vt);
        spm_data_write(Vt,corrPlog10);
      else
        spm_write_vol(Vt,corrPlog10);
      end
    
      %---------------------------------------------------------------
      % save corrected FDR-values for TFCE
      %---------------------------------------------------------------
      fprintf('Save corrected FDR-values.\n');
    
      name = sprintf('TFCE_log_pFDR_%04d',Ic);
      Vt.fname = fullfile(cwd,[name file_ext]);
      if vFWHM > 0
        Vt.descrip = sprintf('TFCE FWHM=%.1fmm, Contrast %04d',vFWHM,Ic);
      else
        Vt.descrip = sprintf('TFCE Contrast %04d',Ic);
      end

      if ~isempty(mask_P)
        [snP_pos,I_pos] = sort(nPtfce(mask_P));
        corrPfdr_pos = snpm_P_FDR([],[],'P',[],snP_pos);
        corrPfdr_pos(I_pos) = corrPfdr_pos;
      end
    
    
      if ~isempty(mask_N)
        [snP_neg,I_neg] = sort(nPtfce(mask_N));
        corrPfdr_neg = snpm_P_FDR([],[],'P',[],snP_neg);
        corrPfdr_neg(I_neg) = corrPfdr_neg;
      end
    
      corrPfdr = NaN(size(t));
      if ~isempty(mask_P)
        corrPfdr(mask_P) = corrPfdr_pos;
      end
      if ~isempty(mask_N)
        corrPfdr(mask_N) = corrPfdr_neg;
      end
    
      corrPfdrlog10 = zeros(size(tfce0));
      corrPfdrlog10(mask_P) = -log10(corrPfdr(mask_P));
      corrPfdrlog10(mask_N) =  log10(corrPfdr(mask_N));
      corrPfdrlog10(mask_0) = NaN;

      if spm12
        Vt = spm_data_hdr_write(Vt);
        spm_data_write(Vt,corrPfdrlog10);
      else
        spm_write_vol(Vt,corrPfdrlog10);
      end
    
      %---------------------------------------------------------------
      % save corrected FDR-values for T
      %---------------------------------------------------------------
      name = sprintf('T_log_pFDR_%04d',Ic);
      Vt.fname = fullfile(cwd,[name file_ext]);
      if vFWHM > 0
        Vt.descrip = sprintf('T FWHM=%.1fmm, Contrast %04d',vFWHM,Ic);
      else
        Vt.descrip = sprintf('T Contrast %04d',Ic);
      end

      if ~isempty(mask_P)
        [snP_pos,I_pos] = sort(nPt(mask_P));
        corrPfdr_pos = snpm_P_FDR([],[],'P',[],snP_pos);
        corrPfdr_pos(I_pos) = corrPfdr_pos;
      end
    
      if ~isempty(mask_N)
        [snP_neg,I_neg] = sort(nPt(mask_N));
        corrPfdr_neg = snpm_P_FDR([],[],'P',[],snP_neg);
        corrPfdr_neg(I_neg) = corrPfdr_neg;
      end
    
      corrPfdr = NaN(size(t));
      if ~isempty(mask_P)
        corrPfdr(mask_P) = corrPfdr_pos;
      end
      if ~isempty(mask_N)
        corrPfdr(mask_N) = corrPfdr_neg;
      end  
    
      corrPfdrlog10 = zeros(size(tfce0));
      corrPfdrlog10(mask_P) = -log10(corrPfdr(mask_P));
      corrPfdrlog10(mask_N) =  log10(corrPfdr(mask_N));
      corrPfdrlog10(mask_0) = NaN;

      if spm12
        Vt = spm_data_hdr_write(Vt);
        spm_data_write(Vt,corrPfdrlog10);
      else
        spm_write_vol(Vt,corrPfdrlog10);
      end
    end % test_mode
  end

colormap(gray)

%---------------------------------------------------------------

function plot_distribution(val_max,val_th,name,alpha,col,order,val0_max)

corr = 1;

n = length(val_th);
sz_val_max = length(val_max);

% allow other thresholds depending on # of permutations
n_alpha = 3;
if sz_val_max < 1000, n_alpha = 2; end
if sz_val_max <  100, n_alpha = 1; end

% with 20 values we have the lowest possible alpha of 0.05
if sz_val_max >= 20
  alpha  = alpha(1:n_alpha);
  val_th = val_th(:,1:n_alpha);

  [hmax, xmax] = hist(val_max, 20);
    
  subplot(2,2,(2*order)-1)
  
  h = bar(xmax,hmax);
  set(h,'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5]);

  avg_h = mean(hmax);
  max_h = max(hmax);
  lim_x = xlim;

  % plot maximum observed value for unpermuted model
  hl = line([val0_max val0_max], [0 max_h]);
  set(hl,'Color','m','LineWidth',2);
  text(0.95*lim_x(2),0.95*max_h,'Max. observed value ',...
    'Color','m','HorizontalAlignment','Right','FontSize',8)
    
  % plot thresholds
  for j=1:n_alpha
    hl = line([val_th(n,j) val_th(n,j)], [0 max_h]);
    set(hl,'Color',col(j,:),'LineStyle','--');
    text(0.95*lim_x(2),(0.5+0.1*j)*max_h,['p<' num2str(alpha(j))],...
      'Color',col(j,:),'HorizontalAlignment','Right','FontSize',8)
  end
  
  ylabel('Frequency');
  xlabel(['Max ' name]);
  if corr
    title(['Distribution of maximum ' name],'FontWeight','bold');
  else
    title(['Distribution of ' name],'FontWeight','bold');
  end
  
  subplot(2,2,2*order)
  
  val_min = min(min(val_th(1:n,:)));
  val_max = max(max(val_th(1:n,:)));
  if val_max/val_min > 10
    semilogy(1:n,val_th(1:n,:))
    yl = log10(ylim);
    ylim(10.^[floor(yl(1)) ceil(yl(2))])
  else
    plot(1:n,val_th(1:n,:))
  end
  if corr
    title(['Corr. threshold of ' name],'FontWeight','bold')
  else
    title(['Uncorr. threshold of ' name],'FontWeight','bold')
  end
  ylabel('Threshold')
  xlabel('Permutations')   
end


%---------------------------------------------------------------
function t = calc_glm(V,X,c,Vmask,vFWHM,TH,W,openmp)
% compute t-statistic using GLM
%
% V      - memory mapped files
% X      - design structure
% c      - contrast
% Vmask  - memory mapped mask image
% vFWHM  - filter width for variance smoothing
% TH     - threshold vector
% W      - whitening matrix
% openmp - use multi-threading

if strcmp(spm('ver'),'SPM12')
  if spm_mesh_detect(V)
    mesh_detected = 1;
  else
    mesh_detected = 0;
  end
else
  mesh_detected = 0;
end

n_data = size(X,1);
n_beta = size(X,2);

X = W*X;
pKX  = pinv(X);
Bcov = pinv(X'*X);
trRV = n_data - rank(X);

vx  = sqrt(sum(V(1).mat(1:3,1:3).^2));

% try mex-file, which is faster
if mesh_detected
   [Beta, ResSS] = calc_beta_mesh(V,Vmask,X,pKX,TH,W);
   dim = V(1).dim(1:2);
else
  dim = V(1).dim(1:3);
  if openmp
    [Beta, ResSS] = cg_glm_get_Beta_ResSS(V,Vmask,X,pKX,TH,W);
    ResSS = reshape(ResSS,dim);
  else
    [Beta, ResSS] = cg_glm_get_Beta_ResSS_noopenmp(V,Vmask,X,pKX,TH,W);
    ResSS = reshape(ResSS,dim);
  end
end

Beta = reshape(Beta,[prod(dim) n_beta]);
ResMS = ResSS/trRV;

con = Beta*c;
clear Beta
con = reshape(con,dim);

% Code for variance smoothing was shameless taken from snpm3 code...
% Blurred mask is used to truncate kernal to brain; if not
% used variance at edges would be underestimated due to
% convolution with zero activity out side the brain.

if vFWHM > 0
  mask = spm_read_vols(Vmask);
  Q = find(mask > 0);
  SmResMS   = zeros(dim);
  SmMask    = zeros(dim);
  TmpVol    = zeros(dim);
  TmpVol(Q) = ones(size(Q));
  spm_smooth(TmpVol,SmMask,vFWHM./vx);
  TmpVol(Q) = ResMS(Q);
  spm_smooth(TmpVol,SmResMS,vFWHM./vx);
  ResMS(Q)  = SmResMS(Q)./SmMask(Q);
end

t = con./(eps+sqrt(ResMS*(c'*Bcov*c)));

%---------------------------------------------------------------
% 
function [beta, ResSS] = calc_beta_mesh(VY,Vm,X,pKX,TH,W);

global Y initialized

n = size(VY,1);
m = size(pKX,1);
beta = zeros([Vm.dim(1:2) m]);
ResSS = zeros(Vm.dim(1:2));

mask = spm_data_read(Vm);
ind = find(mask>0);
 
if ~isempty(ind)
  pXY = zeros([prod(Vm.dim(1:2)) m],'single');
  res = zeros(Vm.dim(1:2));
  
  % load data only the first time
  if ~initialized
    Y = zeros([length(ind) n],'single');

    for i=1:n
      tmp = spm_data_read(VY(i));
      Y(:,i) = single(full(W(i,i)))*tmp(ind);
    end
    initialized = 1;
  end
  
  pXY(ind,:) = Y*single(pKX'); 
  
  res0 = pXY(ind,:)*single(X') - Y; %-Residuals
  res(ind) = double(sum(res0.^2,2));
  ResSS(:,:) = res;   %-Residual SSQ

  pXY = reshape(pXY,[Vm.dim(1:2) m]);
  beta(:,:,:) = double(pXY);
end
 
%---------------------------------------------------------------

function xX = correct_xX(xX)

% vector of covariates and nuisance variables
iCG = [xX.iC xX.iG];
iHB = [xX.iH xX.iB];

% set columns with covariates and nuisance variables to zero
X = xX.X;
X(:,iCG) = 0;

ncol = size(X,2);

% calculate sum of columns
% The idea behind this is that for each factor the sum of all of its columns should be "1".
Xsum = zeros(size(X));
for i=1:ncol
  % only sum up columns without covariates and nuisance variables
  if isempty(find(iCG==i))
    Xsum(:,i) = sum(X(:,1:i),2);
  end
end

% find columns where all entries are constant except zeros entries
% that indicate columns with covariates and nuisance variables
ind = find(any(diff(Xsum))==0 & sum(Xsum)>0);

% no more than 2 factors expected
if length(ind) > 2
  error('Weird design was found that cannot be analyzed correctly.');
end

% correction is only necessary if 2 factors (iH/iB) were found
if length(ind) > 1
  iF = cell(length(ind),1);

  j = 1;
  % skip columns with covariates and nuisance variables
  while find(iCG==j),  j = j + 1; end

  for i=j:length(ind)
    iF{i} = [j:ind(i)];
  
    j = ind(i)+1;
    % skip columns with covariates and nuisance variables
    while find(iCG==j), j = j + 1; end
  end
  
  % not sure whether this will always work but usually iB (subject effects) should be longer than iH (time effects)
  if length(iF{1}) > length(iF{2})
    xX.iB = iF{1};
    xX.iH = iF{2};
  else
    xX.iB = iF{2};
    xX.iH = iF{1};
  end
end

%---------------------------------------------------------------

function str = num2str_short(num)
% get shorther strings for continuous numbers with length > 4

if length(num) > 4
  % check whether vector consist of continuous numbers
  if all(diff(num)==1)
    str = [num2str(num(1)) ':' num2str(num(end))];
  else
    str = num2str(num);
  end
else
  str = num2str(num);
end


%---------------------------------------------------------------
% not used anymore
function [beta, ResSS] = calc_beta(VY,Vm,X,pKX,TH,W);

n = size(VY,1);
m = size(pKX,1);
beta = zeros([Vm.dim(1:3) m]);
ResSS = zeros(Vm.dim(1:3));

for j=1:Vm.dim(3),

 M  = spm_matrix([0 0 j]);
 mask = spm_slice_vol(Vm,M,Vm.dim(1:2),[1 0]);
 ind = find(mask>0);
 
 if ~isempty(ind)
  Y = zeros([length(ind) n],'single');
  pXY = zeros([prod(Vm.dim(1:2)) m],'single');
  res = zeros(Vm.dim(1:2));

  % Load slice j from all images
  for i=1:n
   tmp = spm_slice_vol(VY(i),M,Vm.dim(1:2),[1 0]);
   Y(:,i) = single(W(i,i)*tmp(ind));
  end
  
  pXY(ind,:) = Y*single(pKX'); 
  
  res0 = pXY(ind,:)*single(X') - Y; %-Residuals
  res(ind) = double(sum(res0.^2,2));
  ResSS(:,:,j) = res;   %-Residual SSQ

  pXY = reshape(pXY,[Vm.dim(1:2) m]);
  beta(:,:,j,:) = double(pXY);
 end
end
