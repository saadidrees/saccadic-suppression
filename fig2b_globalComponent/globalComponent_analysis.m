% Fig. 2b, Suppl Fig. 5c, 5d, 7e

% loads data from all periphery saccade experiments. The main data used is the peak response of flashes at different delays after saccade
% and for different conditions. in general, data is organized as [delays,masks,sc,drugs,cells].
% The data file also contains info about each cell's RF and how much of it was covered by the mask.
%% initiate variables
clear all; clc;
load data_peripherySaccade.mat


sc = [25,50,300];
masks = [600,800,1000];
conds_mask = {'full-field','center-600µm','surround-600µm','surround-800µm','surround-1000µm'};
flash_onset = [17,100,250,500,2000];
drugs = {'none','gabazine','gabazine+ptx'};
numSaccs = 40;
dataLen_pre = 400;
dataLen_post = 3600;
saccOnset = 2000;
idx_data = saccOnset-dataLen_pre+1:saccOnset+dataLen_post;

baseFlash_idx = length(flash_onset);      % take the last flash i.e. at 1867ms as the baseline 

%% Calculate modulation indices

% Modulation Index - refined peaks - average across 39 saccade-flash sequences
    mi_avg_flash_ON = (peaks_ref_avg_flash_ON - peaks_ref_avg_baseFlash_ON)./(peaks_ref_avg_flash_ON + peaks_ref_avg_baseFlash_ON);     % [delays,masks,sc,drugs,cells]
    mi_avg_flash_OFF = (peaks_ref_avg_flash_OFF - peaks_ref_avg_baseFlash_OFF)./(peaks_ref_avg_flash_OFF + peaks_ref_avg_baseFlash_OFF);
    
% Modulation Index - refined peaks - for each individual saccade-flash sequence
    mi_allSacc_flash_ON = (peaks_ref_allSacc_flash_ON - peaks_ref_allSacc_baseFlash_ON)./(peaks_ref_allSacc_flash_ON + peaks_ref_allSacc_baseFlash_ON);
    mi_allSacc_flash_OFF = (peaks_ref_allSacc_flash_OFF - peaks_ref_allSacc_baseFlash_OFF)./(peaks_ref_allSacc_flash_OFF + peaks_ref_allSacc_baseFlash_OFF);

        
%% RF stuff
% 1. Calculate which RFs are dirty so to exclude them from further analysis
% 2. Calculate the Masking Factor defined as multiple of sigma of the two dimensional Gaussian fit for which
% the ellipse just touched the mask boundary (Fig. S7d). Cells with RF centers outside the mask are defined to have negative factors
% Cells located close to the mask's edge have masking factors between -1 to +1
% these factors were calculated for each cell in a seperate analysis and are being pooled here

rf_dia_sig = 1;
rf_dia_ON = round(2*rf_dia_sig*nanmax(rf_fit_sigma_ON,[],1)*um2pix);   % rf dia in µm
rf_dia_OFF = round(2*rf_dia_sig*nanmax(rf_fit_sigma_OFF,[],1)*um2pix);   % rf dia in µm

rgb_rf_fit_center_ON = repmat(rf_fit_center_ON,[1,1,size(mask_corn_ON,1)]); rgb_rf_fit_center_ON = permute(rgb_rf_fit_center_ON,[3,1,2]);
idx_outside_ON = squeeze(any([(rgb_rf_fit_center_ON < mask_corn_ON(:,1:2,:)),(rgb_rf_fit_center_ON > mask_corn_ON(:,1:2,:)+mask_corn_ON(:,3:4,:))],2));
assert(any(idx_outside_ON(5,:)==1 & (idx_outside_ON(4,:)==0 | idx_outside_ON(3,:)==0),2),'something wrong in idx_outside_ON')

rgb_rf_fit_center_OFF = repmat(rf_fit_center_OFF,[1,1,size(mask_corn_OFF,1)]); rgb_rf_fit_center_OFF = permute(rgb_rf_fit_center_OFF,[3,1,2]);
idx_outside_OFF = squeeze(any([(rgb_rf_fit_center_OFF < mask_corn_OFF(:,1:2,:)),(rgb_rf_fit_center_OFF > mask_corn_OFF(:,1:2,:)+mask_corn_OFF(:,3:4,:))],2));
assert(any(idx_outside_OFF(5,:)==1 & (idx_outside_OFF(4,:)==0 | idx_outside_OFF(3,:)==0),2),'something wrong in idx_outside_OFF')


rf_dia_ratio_cutOff = 3;
rf_dia_cutOff = 90;    %um

% rf dirty
    idx_rf_dirty_ON = rf_dia_ON < rf_dia_cutOff;        % discar RF less than rf_dia_cutOff µm as it is probably noise
    idx_rf_dirty_ON = idx_rf_dirty_ON | repmat(((nanmax(rf_fit_sigma_ON,[],1)./nanmin(rf_fit_sigma_ON,[],1))>rf_dia_ratio_cutOff),size(idx_rf_dirty_ON,1),1);       % ellipticity of RF
    
    idx_rf_dirty_OFF = rf_dia_OFF < rf_dia_cutOff;        % discar RF less than rf_dia_cutOff µm as it is probably noise
    idx_rf_dirty_OFF = idx_rf_dirty_OFF | repmat(((nanmax(rf_fit_sigma_OFF,[],1)./nanmin(rf_fit_sigma_OFF,[],1))>rf_dia_ratio_cutOff),size(idx_rf_dirty_OFF,1),1);       % ellipticity of RF



% rf outside
    rgb = repmat(~idx_rf_dirty_ON,[size(idx_outside_ON,1),1]) & idx_outside_ON;         % RF outside mask
    idx_rf_outside_ON = rgb;

    rgb = repmat(~idx_rf_dirty_OFF,[size(idx_outside_OFF,1),1]) & idx_outside_OFF;         % RF outside mask
    idx_rf_outside_OFF = rgb;
    
% rf fac combine
    rf_fac_sig_ON = rf_fac_sig_withinMask_ON;
    rf_fac_sig_ON(idx_rf_outside_ON) = -1*rf_fac_sig_outsideMask_ON(idx_rf_outside_ON);
    rf_fac_sig_OFF = rf_fac_sig_withinMask_OFF;
    rf_fac_sig_OFF(idx_rf_outside_OFF) = -1*rf_fac_sig_outsideMask_OFF(idx_rf_outside_OFF);
    
 

%% Perform statistical tests on the refined peaks to check whether the 39 peaks are really modulated or just show random variations
 
 % A sign test is performed here for each cell at each time point and sc to
 % test whether the modulation index is significantly below zero or above
 % it. The sign test is performed on the selected saccs of a single cell.
 % In case of testing for suppression, a left tail test is performed. In
 % case of enhancement, right tail. Power of each test is also computed to
 % see whether the effect we see with the given sample size is actually
 % powerful or not. The 'stats' variable from signtest always gives the number of
 % saccs where mi was above zero. In case of testing for suppression, we
 % consider the left tail as we want the probability for having fewer +ve
 % mi saccs (higher -ve mi saccs)(upper limit is the number of +ve mi saccs). In case of testing for enhancement, we
 % consider right tail as we want the probability of having at least the
 % +ve num or more saccs (num of +ve mi saccs is now the lower limit).
 % http://www.real-statistics.com/binomial-and-related-distributions/statistical-power-binomial-distribution/
    
    p_avg_flash_ON = nan(length(flash_onset)-1,length(conds_mask),length(sc),length(drugs),length(uname_ON));
    pow_p_avg_flash_ON = nan(length(flash_onset)-1,length(conds_mask),length(sc),length(drugs),length(uname_ON));
    sampNum_p_avg_flash_ON = nan(length(flash_onset)-1,length(conds_mask),length(sc),length(drugs),length(uname_ON));
    num_sacc_ON = [];
    
    
    p_avg_flash_OFF = nan(length(flash_onset)-1,length(conds_mask),length(sc),length(drugs),length(uname_OFF));
    pow_p_avg_flash_OFF = nan(length(flash_onset)-1,length(conds_mask),length(sc),length(drugs),length(uname_OFF));
    sampNum_p_avg_flash_off = nan(length(flash_onset)-1,length(conds_mask),length(sc),length(drugs),length(uname_OFF));
    num_sacc_OFF = [];
    
    idx_onCells = false(length(flash_onset)-1,length(conds_mask),length(sc),length(drugs),length(uname_ON));
    idx_offCells = false(length(flash_onset)-1,length(conds_mask),length(sc),length(drugs),length(uname_OFF));
    
    for m = 1:size(mi_allSacc_flash_OFF,6)      % cells
        for l = 1:size(mi_allSacc_flash_OFF,5)      % drugs
            for k = 1:size(mi_allSacc_flash_OFF,4)      % sc
                for j = 1:size(mi_allSacc_flash_OFF,3)  % conds_mask
                    for i = 1:size(mi_allSacc_flash_OFF,1)  % time point
                        rgb = squeeze(mi_allSacc_flash_OFF(i,:,j,k,l,m));

                        idx_toRemove = isnan(rgb) | rgb==inf | rgb==-inf;
                        rgb(idx_toRemove) = [];
                        num_sacc_OFF(i,j,k,l,m) = length(rgb);
                         if ~isempty(rgb)
                            idx_offCells(i,j,k,l,m) = true;
                            if mi_avg_flash_OFF(i,j,k,l,m) < 0
                                [p_avg_flash_OFF(i,j,k,l,m),~,stats] = signtest(rgb,0,'tail','left');
                                pow_p_avg_flash_OFF(i,j,k,l,m) = round(100*sampsizepwr('p',0.50,(stats.sign/num_sacc_OFF(i,j,k,l,m))+eps,[],num_sacc_OFF(i,j,k,l,m),'tail','left'));        
                            else 
                                [p_avg_flash_OFF(i,j,k,l,m),~,stats] = signtest(rgb,0,'tail','right');
                                pow_p_avg_flash_OFF(i,j,k,l,m) = round(100*sampsizepwr('p',0.50,(stats.sign/num_sacc_OFF(i,j,k,l,m))-eps,[],num_sacc_OFF(i,j,k,l,m),'tail','right'));   % 1-binocdf(critcL,num_sacc_OFF(i,j,k),stats.sign/num_sacc_OFF(i,j,k)) where critclL = 25;

                            end

                         end 
                    end
                end
            end
        end
    end
    

    for m = 1:size(mi_allSacc_flash_ON,6)      % cells
        for l = 1:size(mi_allSacc_flash_ON,5)      % drugs
            for k = 1:size(mi_allSacc_flash_ON,4)   % sc
                for j = 1:size(mi_allSacc_flash_ON,3)   % conds_mask
                    for i = 1:size(mi_allSacc_flash_ON,1)   % delays
                        rgb = squeeze(mi_allSacc_flash_ON(i,:,j,k,l,m));
                        idx_toRemove = isnan(rgb) | rgb==inf | rgb==-inf;
                        rgb(idx_toRemove) = [];
                        num_sacc_ON(i,j,k,l,m) = length(rgb);

                        if ~isempty(rgb)
                            idx_onCells(i,j,k,l,m) = true;
                            if mi_avg_flash_ON(i,j,k,l,m) < 0
                                [p_avg_flash_ON(i,j,k,l,m),~,stats] = signtest(rgb,0,'tail','left');
                                pow_p_avg_flash_ON(i,j,k,l,m) = round(100*sampsizepwr('p',0.50,(stats.sign/num_sacc_ON(i,j,k,l,m))+eps,[],num_sacc_ON(i,j,k,l,m),'tail','left'));
                            else 
                                [p_avg_flash_ON(i,j,k,l,m),~,stats] = signtest(rgb,0,'tail','right');
                                pow_p_avg_flash_ON(i,j,k,l,m) = round(100*sampsizepwr('p',0.50,(stats.sign/num_sacc_ON(i,j,k,l,m))-eps,[],num_sacc_ON(i,j,k,l,m),'tail','right'));
                            end

                        end

                    end
                end
            end
        end
    end
  
%% Fig 2b - Line plots

    % use ; to seperate conditions across rows
    rf_fac_sig_cutOff_lowBound = [2];%2;
    rf_fac_sig_cutOff_upBound = [inf]; %inf;
    timePoints = 1:4;
    mask_select = [5];      % conds_mask = {'full-field','center-600µm','surround-600µm','surround-800µm','surround-1000µm'};
    drugs_select = [1;2];     % 1=None,2=gabazine;3=gabazine+ptx
    sc_select = [3];        % sc = [25,50,300];
    sc_text = cellstr(num2str(sc'));
    same_cells = 0;
    ON_OFF_select = 'ALL'; % 'ON' 'OFF' 'ALL'

    plots_idx = 1;
    lim_x = [-100,2100];
    bins = [-1:0.1:1];
    
    col_on = 'r';
    col_off = 'b';
    LINE_WIDTH = 2;
    LINE_TYPES = {'-','--',':','-.'};
    
    
    h_main=figure; suptitle(['ON & OFF Cells - MF: ', num2str(rf_fac_sig_cutOff_lowBound),' to ',num2str(rf_fac_sig_cutOff_upBound)])    
    h_m = [];
    
    s = 1; m = 1; d = 1; t = 0;

    idx_valid_ON = squeeze(all(all(all(~isnan(mi_avg_flash_ON(:,mask_select,sc_select,drugs_select,:)),2),3),4));
    idx_valid_OFF = squeeze(all(all(all(~isnan(mi_avg_flash_OFF(:,mask_select,sc_select,drugs_select,:)),2),3),4));

    p_pop_ON = [];
    p_pop_OFF = [];
    pop_across_data_ON = {}
    pop_across_data_OFF = {}

    counter = 0;
    for c = 1:size(rf_fac_sig_cutOff_lowBound,1)
        for s = 1:size(sc_select,1)
            for m = 1:size(mask_select,1)
                for d = 1:size(drugs_select,1)
                    
                    points_x_ON = [];     
                    for j = 1:size(mi_avg_flash_ON,5)
                        points_x_ON(j,:) = squeeze(mi_avg_flash_ON(:,mask_select(m,1),sc_select(s),drugs_select(d),j));
                    end

                    points_x_OFF = [];     
                    for j = 1:size(mi_avg_flash_OFF,5)
                        points_x_OFF(j,:) = squeeze(mi_avg_flash_OFF(:,mask_select(m,1),sc_select(s),drugs_select(d),j));
                    end


                    idx_respDrugs_ON = all(respConds_ON(drugs_select(d),:),1); idx_respDrugs_ON(isnan(idx_respDrugs_ON)) = 0;
                    idx_respDrugs_ON = logical(idx_respDrugs_ON);

                    idx_respDrugs_OFF = all(respConds_OFF(drugs_select(d),:),1); idx_respDrugs_OFF(isnan(idx_respDrugs_OFF)) = 0;
                    idx_respDrugs_OFF = logical(idx_respDrugs_OFF);

                    if mask_select(m)>2
                        
                        rgb = isnan(rf_fac_sig_ON);       % undefined sigmas
                        rgb = rgb | ((rf_fac_sig_ON < rf_fac_sig_cutOff_lowBound(c)) | (rf_fac_sig_ON > rf_fac_sig_cutOff_upBound(c)));        % smaller sigma than cut off
                        rgb = rgb | repmat(idx_rf_dirty_ON,[size(idx_outside_ON,1),1]); % dirty RFs
                        idx_rf_masked_ON = ~rgb;
                        idx_toKeep_ON = idx_respDrugs_ON & idx_rf_masked_ON(mask_select(m),:) & ~isnan(points_x_ON');

                        rgb = isnan(rf_fac_sig_OFF);       % undefined sigmas
                        rgb = rgb | ((rf_fac_sig_OFF < rf_fac_sig_cutOff_lowBound(c)) | (rf_fac_sig_OFF > rf_fac_sig_cutOff_upBound(c)));        % smaller sigma than cut off
                        rgb = rgb | repmat(idx_rf_dirty_OFF,[size(idx_outside_OFF,1),1]); % dirty RFs
                        idx_rf_masked_OFF = ~rgb;
                        idx_toKeep_OFF = idx_respDrugs_OFF & idx_rf_masked_OFF(mask_select(m),:) & ~isnan(points_x_OFF');
                        
                        
                        
                    else
                        idx_toKeep_ON = repmat(idx_respDrugs_ON,size(points_x_ON,2),1) & ~isnan(points_x_ON');
                        idx_toKeep_OFF = repmat(idx_respDrugs_OFF,size(points_x_OFF,2),1) & ~isnan(points_x_OFF');

                    end
                    
                    if same_cells == 1
                        sameCells_ON = ~isnan(mi_avg_flash_ON(:,mask_select,sc_select,drugs_select,:));   % 1 where the selected masks have responses.
                        sameCells_ON = squeeze(all(all(all(sameCells_ON,2),3),4));
                        
                        sameCells_OFF = ~isnan(mi_avg_flash_OFF(:,mask_select,sc_select,drugs_select,:));   % 1 where the selected masks have responses.
                        sameCells_OFF = squeeze(all(all(all(sameCells_OFF,2),3),4));
                        
                    else
                        sameCells_ON = idx_toKeep_ON;
                        sameCells_OFF = idx_toKeep_OFF;
                    end
                    
                    idx_toKeep_ON = idx_toKeep_ON & sameCells_ON;
                    idx_toKeep_OFF = idx_toKeep_OFF & sameCells_OFF;
                    
                    plot_x_ON  = points_x_ON;
                    plot_x_ON(~idx_toKeep_ON') = nan;
                    std_plot_x_ON = nanstd(plot_x_ON,[],1);
                    sem_plot_x_ON = std_plot_x_ON./sqrt(nansum(~isnan(plot_x_ON),1));
                    rgb_num_ON = nansum(~isnan(plot_x_ON),1);
                    rgb_num_ON = nanmax(rgb_num_ON);
                    mean_plot_x_ON = nanmean(plot_x_ON,1);
                    
                    plot_x_OFF  = points_x_OFF;
                    plot_x_OFF(~idx_toKeep_OFF') = nan;
                    std_plot_x_OFF = nanstd(plot_x_OFF,[],1);
                    sem_plot_x_OFF = std_plot_x_OFF./sqrt(nansum(~isnan(plot_x_OFF),1));
                    rgb_num_OFF = nansum(~isnan(plot_x_OFF),1);
                    rgb_num_OFF = nanmax(rgb_num_OFF);
                    mean_plot_x_OFF = nanmean(plot_x_OFF,1);

                   

                    h_m(1) = subplot(1,1,plots_idx(1));
                    hold on
                    counter = counter+1;

                        switch ON_OFF_select
                            case 'ON'
                                errorbar(flash_onset,[mean_plot_x_ON,0],[sem_plot_x_ON,0],'Color',col_on,'LineWidth',LINE_WIDTH,'LineStyle',LINE_TYPES{counter},'DisplayName',['ON(',num2str(rgb_num_ON),')|',sc_text{sc_select(s)},'µm|',conds_mask{mask_select(m)},'|',drugs{drugs_select(d)}]);
                            case 'OFF'
                                errorbar(flash_onset,[mean_plot_x_OFF,0],[sem_plot_x_OFF,0],'Color',col_off,'LineWidth',LINE_WIDTH,'LineStyle',LINE_TYPES{counter},'DisplayName',['OFF(',num2str(rgb_num_OFF),')|',sc_text{sc_select(s)},'µm|',conds_mask{mask_select(m)},'|',drugs{drugs_select(d)}]);
                            otherwise
                                errorbar(flash_onset,[mean_plot_x_ON,0],[sem_plot_x_ON,0],'Color',col_on,'LineWidth',LINE_WIDTH,'LineStyle',LINE_TYPES{counter},'DisplayName',['ON(',num2str(rgb_num_ON),')|',sc_text{sc_select(s)},'µm|',conds_mask{mask_select(m)},'|',drugs{drugs_select(d)}]);
                                errorbar(flash_onset,[mean_plot_x_OFF,0],[sem_plot_x_OFF,0],'Color',col_off,'LineWidth',LINE_WIDTH,'LineStyle',LINE_TYPES{counter},'DisplayName',['OFF(',num2str(rgb_num_OFF),')|',sc_text{sc_select(s)},'µm|',conds_mask{mask_select(m)},'|',drugs{drugs_select(d)}]);
                        end

                    
                    legend('-DynamicLegend');
                    
                    ylim([-0.6 0.3])
                    ylabel('Modulation index')
                    xlabel('Time from saccade onset (ms)')
                    xlim(lim_x)
                    
                    for i = 1:size(plot_x_ON,2)
                        temp = plot_x_ON(~isnan(plot_x_ON(:,i)),i);
                        [p_pop_ON(i,m,s,d)] = signrank(temp,0);
                        pop_across_data_ON{i,m,s,d} = temp;
                        temp = plot_x_OFF(~isnan(plot_x_OFF(:,i)),i);
                        [p_pop_OFF(i,m,s,d)] = signrank(temp,0);
                        pop_across_data_OFF{i,m,s,d} = temp;
                    end
                end
            end
        end
    end
    plot([-500,2500],[0,0],'--k')
    
    % Compute statistics across two conditions. Note make sure only 2 conditions are selected in the above plot else results will not be right
    p_pop_ON = squeeze(p_pop_ON);   % [squeeze the variable which has dimensions [flash_delays selected,masks selected,sc selected,drugs selected]
    p_pop_OFF = squeeze(p_pop_OFF);
           
    pop_across_data_ON = squeeze(pop_across_data_ON);   % squeeze the variable which has dimensions flash_delays selected,masks selected,sc selected,drugs selected]
    pop_across_data_OFF = squeeze(pop_across_data_OFF);
    
    p_across_ON = cell2mat(cellfun(@(x,y) ranksum(x,y,'tail','both'),pop_across_data_ON(:,1),pop_across_data_ON(:,2),'UniformOutput',0));       % [flash delay,mask,sc,drugs] perform statistics for the two conditions across ON RGCs
    p_across_OFF = cell2mat(cellfun(@(x,y) ranksum(x,y,'tail','both'),pop_across_data_OFF(:,1),pop_across_data_OFF(:,2),'UniformOutput',0));    % same as above for OFF RGCs
    p_across = [flash_onset(timePoints);p_across_ON';p_across_OFF'] % just concatenates values in a convenient form
    
   
%% Suppl Fig. 5c, 5d, 7e -  Scatter Plots
    % flash times are wrt to saccade offset. Add 100 ms for paper version which is wrt to onset
    % use [x,y] for x and y-axis var and ; to seperate conditions across rows. 
    
    for_fig = 'S5c'; % Select one: ['S5c','S5d','S7e']
    
    
    rfDist_measure = 'sig';     % {'sig','dist'}
    timePoints = 1:4;
    sc_select = [3,3];
    ON_OFF_select = 'ALL';   % {'ON','OFF','ALL'}
    

    switch for_fig
        case 'S5c'
            mask_select = [1,5];            
            drugs_select = [1,1];
            rf_fac_sig_cutOff_lowBound = [2];  
            rf_fac_sig_cutOff_upBound = [inf];  
    
        case 'S5d'
            mask_select = [5,5];             
            drugs_select = [1,2];
            rf_fac_sig_cutOff_lowBound = [2];  
            rf_fac_sig_cutOff_upBound = [inf]; 
            
        case 'S7e'
            mask_select = [6,6];             
            drugs_select = [1,1];
            rf_fac_sig_cutOff_lowBound = [-3];  
            rf_fac_sig_cutOff_upBound = [5]; 

    end
    
    
    sc_text = cellstr(num2str(sc'));

    plots_idx = [1:(length(timePoints)*size(mask_select,1)*size(sc_select,1)*size(drugs_select,1)*size(rf_fac_sig_cutOff_lowBound,1))]';
    plots_idx = reshape(plots_idx',length(timePoints),size(mask_select,1),size(sc_select,1),size(drugs_select,1),size(rf_fac_sig_cutOff_lowBound,1));        % [time,mask_conds,sc]
    
    lim_x = [-1.2,1.2];
    lim_y = [-1.2,1.2];
    bins = [-1:0.1:1];
    baseVal = -1.2;
    
    col_on = [0.6,0.6,0.6];
    col_off = [0.2,0.2,0.2];
    
    assert(size(rf_fac_sig_cutOff_lowBound,2) == 1 & size(rf_fac_sig_cutOff_upBound,2),'allowed only as rows')
    assert(size(rf_fac_sig_cutOff_lowBound,1) == size(rf_fac_sig_cutOff_upBound,1),'low and up bound should be in pairs')

    h_main=figure; suptitle(['Scatter Plots - masked cells'])    
    h_m = [];
    
    s = 1; m = 1; d = 1; t = 1;
 
    for c = 1:size(rf_fac_sig_cutOff_lowBound,1)
        for d = 1:size(drugs_select,1)
            for s = 1:size(sc_select,1)
                for m = 1:size(mask_select,1)
                    for t = timePoints;


                        if mask_select(m,1) > 5
                            points_x_ON = rf_fac_sig_ON(5,:)';
                            points_x_OFF = rf_fac_sig_OFF(5,:)';

                            points_y_ON = squeeze(mi_avg_flash_ON(t,5,sc_select(s,2),drugs_select(d,2),:));
                            points_y_OFF = squeeze(mi_avg_flash_OFF(t,5,sc_select(s,2),drugs_select(d,2),:));                       

                        else
                            points_x_ON = squeeze(mi_avg_flash_ON(t,mask_select(m,1),sc_select(s,1),drugs_select(d,1),:));
                            points_y_ON = squeeze(mi_avg_flash_ON(t,mask_select(m,2),sc_select(s,2),drugs_select(d,2),:));

                            points_x_OFF = squeeze(mi_avg_flash_OFF(t,mask_select(m,1),sc_select(s,1),drugs_select(d,1),:));
                            points_y_OFF = squeeze(mi_avg_flash_OFF(t,mask_select(m,2),sc_select(s,2),drugs_select(d,2),:));
                        end


                        idx_respDrugs_ON = respConds_ON(drugs_select(d,:),:); idx_respDrugs_ON(isnan(idx_respDrugs_ON)) = 0;
                        idx_respDrugs_ON = all(logical(idx_respDrugs_ON),1);

                        idx_respDrugs_OFF = respConds_OFF(drugs_select(d,:),:); idx_respDrugs_OFF(isnan(idx_respDrugs_OFF)) = 0;
                        idx_respDrugs_OFF = all(logical(idx_respDrugs_OFF),1);

                        maxMask = max(mask_select(m,:),[],2);
                        if mask_select(m,1)>2 || mask_select(m,2)>2
                            rgb = isnan(rf_fac_sig_ON); 
                            rgb = rgb | (rf_fac_sig_ON < rf_fac_sig_cutOff_lowBound(c)) | (rf_fac_sig_ON > rf_fac_sig_cutOff_upBound(c));        % smaller sigma than cut off
                            rgb = rgb | repmat(idx_rf_dirty_ON,[size(idx_outside_ON,1),1]); % dirty RFs
                            idx_rf_masked_ON = ~rgb;

                            rgb = isnan(rf_fac_sig_OFF); 
                            rgb = rgb | (rf_fac_sig_OFF < rf_fac_sig_cutOff_lowBound(c)) | (rf_fac_sig_OFF > rf_fac_sig_cutOff_upBound(c));        % smaller sigma than cut off
                            rgb = rgb | repmat(idx_rf_dirty_OFF,[size(idx_outside_OFF,1),1]); % dirty RFs
                            idx_rf_masked_OFF = ~rgb;

                            if maxMask < 6
                                idx_toKeep_ON = idx_respDrugs_ON & idx_rf_masked_ON(maxMask,:) & ~isnan(points_x_ON') & ~isnan(points_y_ON');
                                idx_toKeep_OFF = idx_respDrugs_OFF & idx_rf_masked_OFF(maxMask,:) & ~isnan(points_x_OFF') & ~isnan(points_y_OFF');
                            else
                                idx_toKeep_ON = idx_respDrugs_ON & idx_rf_masked_ON(5,:) & ~isnan(points_x_ON') & ~isnan(points_y_ON');
                                idx_toKeep_OFF = idx_respDrugs_OFF & idx_rf_masked_OFF(5,:) & ~isnan(points_x_OFF') & ~isnan(points_y_OFF');
                            end

                        else
                            idx_toKeep_ON = idx_respDrugs_ON & ~isnan(points_x_ON') & ~isnan(points_y_ON');
                            idx_toKeep_OFF = idx_respDrugs_OFF & ~isnan(points_x_OFF') & ~isnan(points_y_OFF');
                        end

                        plot_x_ON  = points_x_ON(idx_toKeep_ON);
                        plot_y_ON  = points_y_ON(idx_toKeep_ON);

                        plot_x_OFF  = points_x_OFF(idx_toKeep_OFF);
                        plot_y_OFF  = points_y_OFF(idx_toKeep_OFF);
                        
                        switch ON_OFF_select
                            case 'ON'
                                plot_x_OFF = [];
                                plot_y_OFF = [];
                                
                            case 'OFF'
                                plot_x_ON = [];
                                plot_y_ON = [];
                        end


                        plots_uname_ON = uname_ON(idx_toKeep_ON');
                        plots_uname_OFF = uname_OFF(idx_toKeep_OFF');

                        rf_dia_ON_temp = rf_dia_ON(idx_toKeep_ON')*1.5;
                        rf_dia_OFF_temp = rf_dia_OFF(idx_toKeep_OFF')*1.5;

                        rf_bins = [100,200,300,500,inf];
                        [~,~,c1] = histcounts(rf_dia_ON_temp,rf_bins);
                        z1 = rf_bins(c1);
                        [~,~,c2] = histcounts(rf_dia_OFF_temp,rf_bins);
                        z2 = rf_bins(c2);

                        colormap_rf = [1,0,0;0,1,0;0,0,1;1,0,1];%;0,0,0;1,1,1];
                        
                        h_m(t,m,s,d,c) = subplot(size(rf_fac_sig_cutOff_lowBound,1)*size(mask_select,1)*size(sc_select,1)*size(drugs_select,1),length(timePoints),plots_idx(t,m,s,d,c));
                        hold on
                        scatter(plot_x_ON,plot_y_ON,'MarkerEdgeColor',col_on,'MarkerFaceColor',col_on) 
                        scatter(plot_x_OFF,plot_y_OFF,'MarkerEdgeColor',col_off,'MarkerFaceColor',col_off)



                        if maxMask < 6

                            [f_off_all,x] = hist(plot_y_OFF,bins);
                            h = barh(x,(f_off_all/trapz(f_off_all))+baseVal);
                            set(h,'FaceColor',col_off+.05,'EdgeColor','k');
                            h.BaseValue = baseVal;


                            [f_on_all,x] = hist(plot_y_ON,bins);
                            h = barh(x,(f_on_all/trapz(f_on_all))+baseVal);
                            set(h,'FaceColor',col_on+.05,'EdgeColor','k');
                            h.BaseValue = baseVal;
                            set(h,'FaceAlpha',0.6)

                            [f_off_all,x] = hist(plot_x_OFF,bins);
                            h = bar(x,(f_off_all/trapz(f_off_all))+baseVal);
                            set(h,'FaceColor',col_off+.05,'EdgeColor','k');
                            h.BaseValue = baseVal;

                            [f_on_all,x] = hist(plot_x_ON,bins);
                            h = bar(x,(f_on_all/trapz(f_on_all))+baseVal);
                            set(h,'FaceColor',col_on+.05,'EdgeColor','k');
                            h.BaseValue = baseVal;
                            set(h,'FaceAlpha',0.6)

                            plot([-1,1],[0,0],'--','color',[0.8,0.8,0.8])
                            plot([0,0],[-1,1],'--','color',[0.8,0.8,0.8])
                            plot([-1,1],[-1,1],'k')
                            xlim(lim_x)
                            ylim(lim_y)

                        else
                            plot([min([plot_x_ON;plot_x_OFF]),max([plot_x_ON;plot_x_OFF])],[0,0],'--','color','r')
                            switch rfDist_measure
                                case 'dist'
                                    plot([500,500],[-1,1],'--','color','b')
    %                                 xlim([0 4]);%max(plot_x_ON)+.2]);
                                case 'sig'
                                    plot([0,0],[-1,1],'--','color','b')
    %                                 rgb1 = repmat([1:10],2,1);
    %                                 rgb2 = repmat([-1,1]',1,10);
    %                                 plot(rgb1,rgb2,'--','color','y')
                                    xlim([rf_fac_sig_cutOff_lowBound-.5,rf_fac_sig_cutOff_upBound+.5]) %([-5 10]);%max(plot_x_ON)+.2]);
                            end
    %                         xlim([0 4]);%max(plot_x_ON)+.2]);
                            ylim(lim_y)                       
                        end

                        legend({['ON ',num2str(length(plot_x_ON))],['OFF: ',num2str(length(plot_x_OFF))]})

                        axis square
                        title([num2str(flash_onset(t)),'ms | masking factor: [',num2str(rf_fac_sig_cutOff_lowBound(c)),',',num2str(rf_fac_sig_cutOff_upBound(c)),']'])

                    end


                    if maxMask < 6
                        xlabel(['sc: ',sc_text{sc_select(s,1)},' µm | ', conds_mask{mask_select(m,1)},' | drugs: ',drugs{drugs_select(d,1)}])
                        ylabel(['sc: ',sc_text{sc_select(s,2)},' µm | ', conds_mask{mask_select(m,2)},' | drugs: ',drugs{drugs_select(d,2)}])
                    else
                        xlabel(['rf metric: ',rfDist_measure])
                        ylabel(['sc: ',sc_text{sc_select(s,2)},' µm | ',drugs{drugs_select(d,2)}])
                    end
                    
   
                end
            end
        end
    end


%% Fig. S5b
% plot RFs using mask_facs - 

    % use ; to seperate conditions across rows
    timePoints = 1:4;
    mask_select = [5];      % conds_mask = {'full-field','center-600µm','surround-600µm','surround-800µm','surround-1000µm'};
    drugs_select = [1];     % drugs = {'none','gabazine','gabazine+ptx'};
    sc_select = [3];        % sc = [25,50,300];
    facs = [2];             % masking factors
    rf_fac_sig_cutOff_lowBound = [2]; %[-inf;-2;-1;0;1;2;3];
    rf_fac_sig_cutOff_upBound = [inf];

    sc_text = cellstr(num2str(sc'));

   
    plots_idx = [1,2;3,4];
    
    lim_x = [-1.2,1.2];
    bins = [-1:0.1:1];
    
    col_on = [0.6,0.6,0.6];
    col_off = [0.2,0.2,0.2];
    
    assert(size(rf_fac_sig_cutOff_lowBound,2) == 1 & size(rf_fac_sig_cutOff_upBound,2),'allowed only as rows')
    assert(size(rf_fac_sig_cutOff_lowBound,1) == size(rf_fac_sig_cutOff_upBound,1),'low and up bound should be in pairs')
 
%     h_main=figure; suptitle(['ON & OFF Cells'])    
%     h_m = [];
%     
    num_ONCells = [];
    num_ONCells_sig = [];
    
    s = 1; m = 1; d = 1; t = 1; c = 1;
    cat_idx_toKeep_ON = [];
    cat_idx_toKeep_OFF = [];


    for c = 1:size(rf_fac_sig_cutOff_lowBound,1)
        for s = 1:size(sc_select,1)
            for m = 1:size(mask_select,1)
                for d = 1:size(drugs_select,1)
                    for t = timePoints

                        points_x_ON = [];     
                        for j = 1:size(mi_avg_flash_ON,5)
                            points_x_ON(j,:) = squeeze(mi_avg_flash_ON(t,mask_select(m,1),sc_select(s),drugs_select(d),j));
                        end

                        points_x_OFF = [];     
                        for j = 1:size(mi_avg_flash_OFF,5)
                            points_x_OFF(j,:) = squeeze(mi_avg_flash_OFF(t,mask_select(m,1),sc_select(s),drugs_select(d),j));
                        end

                        idx_respDrugs_ON = respConds_ON(drugs_select(d),:); idx_respDrugs_ON(isnan(idx_respDrugs_ON)) = 0;
                        idx_respDrugs_ON = logical(idx_respDrugs_ON);

                        idx_respDrugs_OFF = respConds_OFF(drugs_select(d),:); idx_respDrugs_OFF(isnan(idx_respDrugs_OFF)) = 0;
                        idx_respDrugs_OFF = logical(idx_respDrugs_OFF);

                        if mask_select(m)>2
                            rgb = isnan(rf_fac_sig_ON);       % undefined sigmas
                            rgb = rgb | ((rf_fac_sig_ON < rf_fac_sig_cutOff_lowBound(c)) | (rf_fac_sig_ON > rf_fac_sig_cutOff_upBound(c)));        % smaller sigma than cut off
                            rgb = rgb | repmat(idx_rf_dirty_ON,[size(idx_outside_ON,1),1]); % dirty RFs
                            idx_rf_masked_ON = ~rgb;
                            idx_toKeep_ON = idx_respDrugs_ON & idx_rf_masked_ON(mask_select(m),:) & ~isnan(points_x_ON');

                            rgb = isnan(rf_fac_sig_OFF);       % undefined sigmas
                            rgb = rgb | ((rf_fac_sig_OFF < rf_fac_sig_cutOff_lowBound(c)) | (rf_fac_sig_OFF > rf_fac_sig_cutOff_upBound(c)));        % smaller sigma than cut off
                            rgb = rgb | repmat(idx_rf_dirty_OFF,[size(idx_outside_OFF,1),1]); % dirty RFs
                            idx_rf_masked_OFF = ~rgb;
                            idx_toKeep_OFF = idx_respDrugs_OFF & idx_rf_masked_OFF(mask_select(m),:) & ~isnan(points_x_OFF');

                         else
                            idx_toKeep_ON = idx_respDrugs_ON & ~isnan(points_x_ON');
                            idx_toKeep_OFF = idx_respDrugs_OFF & ~isnan(points_x_OFF');

                        end
                        
                        cat_idx_toKeep_ON(:,t,m,s,d,c) = idx_toKeep_ON;
                        cat_idx_toKeep_OFF(:,t,m,s,d,c) = idx_toKeep_OFF;
                    end
                end
            end
        end
    end
    
        
    mf_idx_toKeep_ON = squeeze(all(all(all(all(cat_idx_toKeep_ON,2),3),4),5));
    mf_idx_toKeep_OFF = squeeze(all(all(all(all(cat_idx_toKeep_OFF,2),3),4),5));
    
%     assert(sum(mf_idx_toKeep_ON(:,1) & mf_idx_toKeep_ON(:,2)) == 0,'overlapping ON cells');
%     assert(sum(mf_idx_toKeep_OFF(:,1) & mf_idx_toKeep_OFF(:,2)) == 0,'overlapping OFF cells');

    t = -pi:0.01:pi;


    idx_firstNonNan = find(~isnan(rf_fac_sig_ON(mask_select,:)),1,'first');

    hypoth_win = [0,0,4000,2000];
    start_corn = [1500,500];

    figure;hold on

   
   for c = 1:size(rf_fac_sig_cutOff_lowBound,1)
 
        mf_idx_ON = find(mf_idx_toKeep_ON(:,c));
        mf_idx_OFF = find(mf_idx_toKeep_OFF(:,c));
        subplot(2,2,plots_idx(c,1));hold on
        rectangle('position',hypoth_win);
        rectangle('position',[start_corn,mask_corn_ON(5,3:4,idx_firstNonNan)*um2pix]);
        
        for u = 1:length(mf_idx_ON)
            ellipse_coord = {};
            for i = 1:length(facs)
                x = facs(i)*rf_fit_sigma_ON(1,mf_idx_ON(u)).*cos(t')*um2pix;
                y = facs(i)*rf_fit_sigma_ON(2,mf_idx_ON(u)).*sin(t')*um2pix;
                transformMatrix = [cosd(rf_fit_angle_ON(mf_idx_ON(u))),sind(rf_fit_angle_ON(mf_idx_ON(u)));-sind(rf_fit_angle_ON(mf_idx_ON(u))),cosd(rf_fit_angle_ON(mf_idx_ON(u)))]; % elipse roation

                rf_rotated = [x,y]*transformMatrix;
                center_relative = rf_fit_center_ON(:,mf_idx_ON(u))'*um2pix - mask_corn_ON(5,1:2,mf_idx_ON(u))*um2pix;
                ellipse_coord{i}= bsxfun(@plus,rf_rotated,center_relative+start_corn);

        %                 rectangle('position',[start_corn,mask_corn_ON(5,3:4,u)]);
            end
            h_l = plot(ellipse_coord{1}(:,1),ellipse_coord{1}(:,2),'linewidth',1);
        %         plot(ellipse_coord{2}(:,1),ellipse_coord{2}(:,2),'--','linewidth',1,'Color',h_l.Color)
            set(gca,'Ydir','reverse')
            

        end
        xlim([-100 4100])
        ylim([-100 2300])
        title([num2str(u),' ON Cells | low: ',num2str(rf_fac_sig_cutOff_lowBound(c)),' | up: ',num2str(rf_fac_sig_cutOff_upBound(c))])
        
        
        subplot(2,2,plots_idx(c,2));hold on
        rectangle('position',hypoth_win);
        rectangle('position',[start_corn,mask_corn_ON(5,3:4,idx_firstNonNan)*um2pix]);

        for u = 1:length(mf_idx_OFF)
            ellipse_coord = {};
            for i = 1:length(facs)
                x = facs(i)*rf_fit_sigma_OFF(1,mf_idx_OFF(u)).*cos(t')*um2pix;
                y = facs(i)*rf_fit_sigma_OFF(2,mf_idx_OFF(u)).*sin(t')*um2pix;
                transformMatrix = [cosd(rf_fit_angle_OFF(mf_idx_OFF(u))),sind(rf_fit_angle_OFF(mf_idx_OFF(u)));-sind(rf_fit_angle_OFF(mf_idx_OFF(u))),cosd(rf_fit_angle_OFF(mf_idx_OFF(u)))]; % elipse roation

                rf_rotated = [x,y]*transformMatrix;
                center_relative = rf_fit_center_OFF(:,mf_idx_OFF(u))'*um2pix - mask_corn_OFF(5,1:2,mf_idx_OFF(u))*um2pix;
                ellipse_coord{i}= bsxfun(@plus,rf_rotated,center_relative+start_corn);

        %                 rectangle('position',[start_corn,mask_corn_ON(5,3:4,u)]);
            end
            h_l = plot(ellipse_coord{1}(:,1),ellipse_coord{1}(:,2),'linewidth',1);
%             text(ellipse_coord{1}(1,1),ellipse_coord{1}(1,2),num2str(u))
        %         plot(ellipse_coord{2}(:,1),ellipse_coord{2}(:,2),'--','linewidth',1,'Color',h_l.Color)
            set(gca,'Ydir','reverse')

        end
        xlim([-100 4100])
        ylim([-100 2300])
        title([num2str(u),' OFF Cells | low: ',num2str(rf_fac_sig_cutOff_lowBound(c)),' | up: ',num2str(rf_fac_sig_cutOff_upBound(c))])
        
   end

























