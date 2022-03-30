
%% Fig. 4 and Suppl Fig 8
% load data from luminance step experiments. The peak variables have peak responses of flashes after subtracting the luminance step response.
% These responses are arranged as [delays,contr_type,drugs,units]

%%
clear all; clc; set(0,'defaultTextInterpreter','none');
load data_luminanceSteps_assocIndex.mat
dataLen_pre = 100;
dataLen_post = 2500;
saccOnset = 2000;
idx_data = saccOnset-dataLen_pre+1:saccOnset+dataLen_post;
drugs = {'none','gabazine+ptx+str'};

bins_contr_pos = [0:0.025:0.5];        % low,medium, high
N_bins_contr = length(bins_contr_pos)-1;

flash_onset = [17,33,50,100,250,500,1000,2000];
baseFlash_idx = length(flash_onset);

nds = [5];
contrast_polarity = {'Dark','Bright'};

%% neg and positive contrasts pooled

peaks_ref_allCombContr_flash_ON = squeeze(cat(2,peaks_ref_allSepContr_flash_ON(:,:,1,:,:),peaks_ref_allSepContr_flash_ON(:,:,2,:,:)));
peaks_ref_allCombContr_flash_OFF = squeeze(cat(2,peaks_ref_allSepContr_flash_OFF(:,:,1,:,:),peaks_ref_allSepContr_flash_OFF(:,:,2,:,:)));

peaks_ref_allCombContr_baseFlash_ON = squeeze(cat(2,peaks_ref_allSepContr_baseFlash_ON(:,:,1,:,:),peaks_ref_allSepContr_baseFlash_ON(:,:,2,:,:)));
peaks_ref_allCombContr_baseFlash_OFF = squeeze(cat(2,peaks_ref_allSepContr_baseFlash_OFF(:,:,1,:,:),peaks_ref_allSepContr_baseFlash_OFF(:,:,2,:,:)));

area_allCombContr_sacc_ON = squeeze(cat(2,area_allSepContr_sacc_ON(:,:,1,:,:),area_allSepContr_sacc_ON(:,:,2,:,:)));
area_allCombContr_sacc_OFF = squeeze(cat(2,area_allSepContr_sacc_OFF(:,:,1,:,:),area_allSepContr_sacc_OFF(:,:,2,:,:)));

%% Modulation Index - (a-b)/(a+b)
    % seperate negative and positive contrasts
    mi_avgSepContr_ON = (peaks_ref_avgSepContr_flash_ON - peaks_ref_avgSepContr_baseFlash_ON)./(peaks_ref_avgSepContr_flash_ON + peaks_ref_avgSepContr_baseFlash_ON);       % [delays,contr_type,drugs,units]
    mi_avgSepContr_OFF = (peaks_ref_avgSepContr_flash_OFF - peaks_ref_avgSepContr_baseFlash_OFF)./(peaks_ref_avgSepContr_flash_OFF + peaks_ref_avgSepContr_baseFlash_OFF);  % [delays,contr_type,drugs,units]
    
    mi_allSepContr_ON = (peaks_ref_allSepContr_flash_ON - peaks_ref_allSepContr_baseFlash_ON)./(peaks_ref_allSepContr_flash_ON + peaks_ref_allSepContr_baseFlash_ON);
    mi_allSepContr_OFF = (peaks_ref_allSepContr_flash_OFF - peaks_ref_allSepContr_baseFlash_OFF)./(peaks_ref_allSepContr_flash_OFF + peaks_ref_allSepContr_baseFlash_OFF);
        
    % pooled negative and poisitive contrasts
    mi_avgCombContr_ON = (peaks_ref_avgCombContr_flash_ON - peaks_ref_avgCombContr_baseFlash_ON)./(peaks_ref_avgCombContr_flash_ON + peaks_ref_avgCombContr_baseFlash_ON);
    mi_avgCombContr_OFF = (peaks_ref_avgCombContr_flash_OFF - peaks_ref_avgCombContr_baseFlash_OFF)./(peaks_ref_avgCombContr_flash_OFF + peaks_ref_avgCombContr_baseFlash_OFF);
    
    mi_allCombContr_ON = (peaks_ref_allCombContr_flash_ON - peaks_ref_allCombContr_baseFlash_ON)./(peaks_ref_allCombContr_flash_ON + peaks_ref_allCombContr_baseFlash_ON);
    mi_allCombContr_OFF = (peaks_ref_allCombContr_flash_OFF - peaks_ref_allCombContr_baseFlash_OFF)./(peaks_ref_allCombContr_flash_OFF + peaks_ref_allCombContr_baseFlash_OFF);
    
%% Association index for association/correlation between modulation index and saccade response
 
 % ON CElls   
    assoc_R_sepContr_ON = nan(size(mi_avgSepContr_ON));
    assoc_p_sepContr_ON = nan(size(mi_avgSepContr_ON));
    
    assoc_R_combContr_ON = nan(size(mi_avgCombContr_ON));
    assoc_p_combContr_ON = nan(size(mi_avgCombContr_ON));
    
    for l = 1:size(area_allSepContr_sacc_ON,5)       % units
        for k = 1:size(area_allSepContr_sacc_ON,4)       % drugs
            for i = 1:length(flash_onset)-1
                for j = 1:size(area_allSepContr_sacc_ON,3)       % contr_type
                    a = squeeze(area_allSepContr_sacc_ON(i,:,j,k,l));a = a(:);      % sacc area
                    b = squeeze(mi_allSepContr_ON(i,:,j,k,l)); b = b(:);  % flashes
                    idx = (isnan(a) | isnan(b));
                    a(idx) = []; b(idx) = [];

                    if ~isempty(a) || ~isempty(b)
                        [assoc_R_sepContr_ON(i,j,k,l),assoc_p_sepContr_ON(i,j,k,l)] = corr(a,b,'type','Spearman');
                    end
                end
                    a = squeeze(area_allCombContr_sacc_ON(i,:,k,l));a = a(:);      % sacc area
                    b = squeeze(mi_allCombContr_ON(i,:,k,l)); b = b(:);  % flashes
                    idx = (isnan(a) | isnan(b));
                    a(idx) = []; b(idx) = [];

                    if ~isempty(a) || ~isempty(b)
                        [assoc_R_combContr_ON(i,k,l),assoc_p_sepContr_ON(i,k,l)] = corr(a,b,'type','Spearman');
                    end
                
            end
        end
    end
    
 % OFF CElls   
    assoc_R_sepContr_OFF = nan(size(mi_avgSepContr_OFF));
    assoc_p_sepContr_OFF = nan(size(mi_avgSepContr_OFF));
    
    assoc_R_combContr_OFF = nan(size(mi_avgCombContr_OFF));
    assoc_p_combContr_OFF = nan(size(mi_avgCombContr_OFF));
    
    for l = 1:size(area_allSepContr_sacc_OFF,5)       % units
        for k = 1:size(area_allSepContr_sacc_OFF,4)       % drugs
            for i = 1:length(flash_onset)-1
                for j = 1:size(area_allSepContr_sacc_OFF,3)       % contr_type
                    a = squeeze(area_allSepContr_sacc_OFF(i,:,j,k,l));a = a(:);      % sacc area
                    b = squeeze(mi_allSepContr_OFF(i,:,j,k,l)); b = b(:);  % flashes
                    idx = (isnan(a) | isnan(b));
                    a(idx) = []; b(idx) = [];

                    if ~isempty(a) || ~isempty(b)
                        [assoc_R_sepContr_OFF(i,j,k,l),assoc_p_sepContr_OFF(i,j,k,l)] = corr(a,b,'type','Spearman');
                    end
                end
                    a = squeeze(area_allCombContr_sacc_OFF(i,:,k,l));a = a(:);      % sacc area
                    b = squeeze(mi_allCombContr_OFF(i,:,k,l)); b = b(:);  % flashes
                    idx = (isnan(a) | isnan(b));
                    a(idx) = []; b(idx) = [];

                    if ~isempty(a) || ~isempty(b)
                        [assoc_R_combContr_OFF(i,k,l),assoc_p_sepContr_OFF(i,k,l)] = corr(a,b,'type','Spearman');
                    end
                
            end
        end
    end

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
    
 % all steps
    p_avgSepContr_ON = nan(length(flash_onset)-1,length(contrast_polarity),length(drugs),length(uname_ON));
    pow_p_avgSepContr_ON = nan(length(flash_onset)-1,length(contrast_polarity),length(drugs),length(uname_ON));
    sampNum_p_avgSepContr_ON = nan(length(flash_onset)-1,length(contrast_polarity),length(drugs),length(uname_ON));
    num_sacc_avgSepContr_ON = [];
    

    p_avgSepContr_OFF = nan(length(flash_onset)-1,length(contrast_polarity),length(drugs),length(uname_OFF));
    pow_p_avgSepContr_OFF = nan(length(flash_onset)-1,length(contrast_polarity),length(drugs),length(uname_OFF));
    sampNum_p_avgSepContr_OFF = nan(length(flash_onset)-1,length(contrast_polarity),length(drugs),length(uname_OFF));
    num_sacc_avgSepContr_OFF = [];
    
    idx_cells_avgSepContr_ON = false(length(flash_onset)-1,length(contrast_polarity),length(drugs),length(uname_ON));
    idx_cells_avgSepContr_OFF = false(length(flash_onset)-1,length(contrast_polarity),length(drugs),length(uname_OFF));
 
% combined contrast polarity
    p_avgCombContr_ON = nan(length(flash_onset)-1,length(drugs),length(uname_ON));
    pow_p_avgCombContr_ON = nan(length(flash_onset)-1,length(drugs),length(uname_ON));
    sampNum_p_avgCombContr_ON = nan(length(flash_onset)-1,length(drugs),length(uname_ON));
    num_sacc_avgCombContr_ON = [];
    

    p_avgCombContr_OFF = nan(length(flash_onset)-1,length(drugs),length(uname_OFF));
    pow_p_avgCombContr_OFF = nan(length(flash_onset)-1,length(drugs),length(uname_OFF));
    sampNum_p_avgCombContr_OFF = nan(length(flash_onset)-1,length(drugs),length(uname_OFF));
    num_sacc_avgCombContr_OFF = [];
    
    idx_cells_avgCombContr_ON = false(length(flash_onset)-1,length(drugs),length(uname_ON));
    idx_cells_avgCombContr_OFF = false(length(flash_onset)-1,length(drugs),length(uname_OFF));



    for l = 1:size(mi_allSepContr_ON,5)      % units   
        for k = 1:size(mi_allSepContr_ON,4)      % drugs
            for i = 1:size(mi_allSepContr_ON,1)  % delays
                for j = 1:size(mi_allSepContr_ON,3)  % contrast type

                    rgb = squeeze(mi_allSepContr_ON(i,:,j,k,l));
                    idx_toRemove = isnan(rgb) | rgb==inf | rgb==-inf;
                    rgb(idx_toRemove) = [];
                    num_sacc_avgSepContr_ON(i,j,k,l) = length(rgb);

                    if ~isempty(rgb)
                        idx_cells_avgSepContr_ON(i,j,k,l) = true;
                        if mi_avgSepContr_ON(i,j,k,l) < 0
                            [p_avgSepContr_ON(i,j,k,l),~,stats] = signtest(rgb,0,'tail','left');
                            pow_p_avgSepContr_ON(i,j,k,l) = round(100*sampsizepwr('p',0.50,(stats.sign/num_sacc_avgSepContr_ON(i,j,k,l))+eps,[],num_sacc_avgSepContr_ON(i,j,k,l),'tail','left'));
                        else 
                            [p_avgSepContr_ON(i,j,k,l),~,stats] = signtest(rgb,0,'tail','right');
                            pow_p_avgSepContr_ON(i,j,k,l) = round(100*sampsizepwr('p',0.50,(stats.sign/num_sacc_avgSepContr_ON(i,j,k,l))-eps,[],num_sacc_avgSepContr_ON(i,j,k,l),'tail','right'));
                        end
                    end

                end
                
                rgb = squeeze(mi_allCombContr_ON(i,:,k,l));
                idx_toRemove = isnan(rgb) | rgb==inf | rgb==-inf;
                rgb(idx_toRemove) = [];
                num_sacc_avgCombContr_ON(i,k,l) = length(rgb);

                if ~isempty(rgb)
                    idx_cells_avgCombContr_ON(i,k,l) = true;
                    if mi_avgCombContr_ON(i,k,l) < 0
                        [p_avgCombContr_ON(i,k,l),~,stats] = signtest(rgb,0,'tail','left');
                        pow_p_avgCombContr_ON(i,k,l) = round(100*sampsizepwr('p',0.50,(stats.sign/num_sacc_avgCombContr_ON(i,k,l))+eps,[],num_sacc_avgCombContr_ON(i,k,l),'tail','left'));
                    else 
                        [p_avgCombContr_ON(i,k,l),~,stats] = signtest(rgb,0,'tail','right');
                        pow_p_avgCombContr_ON(i,k,l) = round(100*sampsizepwr('p',0.50,(stats.sign/num_sacc_avgCombContr_ON(i,k,l))-eps,[],num_sacc_avgCombContr_ON(i,k,l),'tail','right'));
                    end
                end
            end
        end
    end
    
    for l = 1:size(mi_allSepContr_OFF,5)      % units   
        for k = 1:size(mi_allSepContr_OFF,4)      % drugs
            for i = 1:size(mi_allSepContr_OFF,1)  % delays
                for j = 1:size(mi_allSepContr_OFF,3)  % contrast type

                    rgb = squeeze(mi_allSepContr_OFF(i,:,j,k,l));
                    idx_toRemove = isnan(rgb) | rgb==inf | rgb==-inf;
                    rgb(idx_toRemove) = [];
                    num_sacc_avgSepContr_OFF(i,j,k,l) = length(rgb);

                    if ~isempty(rgb)
                        idx_cells_avgSepContr_OFF(i,j,k,l) = true;
                        if mi_avgSepContr_OFF(i,j,k,l) < 0
                            [p_avgSepContr_OFF(i,j,k,l),~,stats] = signtest(rgb,0,'tail','left');
                            pow_p_avgSepContr_OFF(i,j,k,l) = round(100*sampsizepwr('p',0.50,(stats.sign/num_sacc_avgSepContr_OFF(i,j,k,l))+eps,[],num_sacc_avgSepContr_OFF(i,j,k,l),'tail','left'));
                        else 
                            [p_avgSepContr_OFF(i,j,k,l),~,stats] = signtest(rgb,0,'tail','right');
                            pow_p_avgSepContr_OFF(i,j,k,l) = round(100*sampsizepwr('p',0.50,(stats.sign/num_sacc_avgSepContr_OFF(i,j,k,l))-eps,[],num_sacc_avgSepContr_OFF(i,j,k,l),'tail','right'));
                        end
                    end

                end
                
                rgb = squeeze(mi_allCombContr_OFF(i,:,k,l));
                idx_toRemove = isnan(rgb) | rgb==inf | rgb==-inf;
                rgb(idx_toRemove) = [];
                num_sacc_avgCombContr_OFF(i,k,l) = length(rgb);

                if ~isempty(rgb)
                    idx_cells_avgCombContr_OFF(i,k,l) = true;
                    if mi_avgCombContr_OFF(i,k,l) < 0
                        [p_avgCombContr_OFF(i,k,l),~,stats] = signtest(rgb,0,'tail','left');
                        pow_p_avgCombContr_OFF(i,k,l) = round(100*sampsizepwr('p',0.50,(stats.sign/num_sacc_avgCombContr_OFF(i,k,l))+eps,[],num_sacc_avgCombContr_OFF(i,k,l),'tail','left'));
                    else 
                        [p_avgCombContr_OFF(i,k,l),~,stats] = signtest(rgb,0,'tail','right');
                        pow_p_avgCombContr_OFF(i,k,l) = round(100*sampsizepwr('p',0.50,(stats.sign/num_sacc_avgCombContr_OFF(i,k,l))-eps,[],num_sacc_avgCombContr_OFF(i,k,l),'tail','right'));
                    end
                end
            end
        end
    end    

%% Fig. 4

select_fig = '4c'

switch select_fig
    case '4a'
        contrastPolarity_select = [5];  % {'neg'    'pos'    'preferred'    'anti-preferred'    'combined'}
        drugs_select = [1];   % drugs = {'none','gabazine+ptx+str'};
        
    case '4b'
        contrastPolarity_select = [2;1];  % {'neg'    'pos'    'preferred'    'anti-preferred'    'combined'}
        drugs_select = [1];   % drugs = {'none','gabazine+ptx+str'};
        
    case '4c'
        contrastPolarity_select = [2;1];  % {'neg'    'pos'    'preferred'    'anti-preferred'    'combined'}
        drugs_select = [1;2];   % drugs = {'none','gabazine+ptx+str'};
end
    timePoints = 1:7;
    
    num_onCells = nan(length(flash_onset));
    num_offCells = nan(length(flash_onset));
    
    exps_select = [1:length(dates)];        % select which experiments data to use
%     contrastPolarity_select = [2];  % {'neg'    'pos'    'preferred'    'anti-preferred'    'combined'}
%     drugs_select = [1;2];   % drugs = {'none','gabazine+ptx+str'};
    
    contrast_conds = {'neg','pos','preferred','anti-preferred','combined'};
    idx_contrast_conds = [1,1;2,2;2,1;1,2]; % [off,on] rows are the contrast conditions [neg,pos,preferred,antipreferred]
    
    plots_idx = 1;
    
    lim_x = [-200,2100];
    bins = [-1:0.1:1];
    
    col_on = 'r';%[0.6,0.6,0.6];
    col_off = 'b';%[0.2,0.2,0.2];
    LINE_WIDTH = 2;
    LINE_TYPES = {'-','--',':','-.'};

    
    s = 1; d = 1; t = 1; cc = 1;
    counter = 0;
    
    
    idx_cells_expsSelect_ON = ismember(unitExpIds_ON,exps_select);
    idx_cells_expsSelect_OFF = ismember(unitExpIds_OFF,exps_select);
    
    idx_valid_ON = squeeze(all(~isnan(mi_avgCombContr_ON(:,drugs_select,:)),2));
    idx_valid_OFF = squeeze(all(~isnan(mi_avgCombContr_OFF(:,drugs_select,:)),2));

    h_main=figure; suptitle(['Contrast steps'])    
    h_m = [];
    f_on_all = [];
    f_off_all = [];
    p_pop_ON = [];
    p_pop_OFF = [];

    for cc = 1:size(contrastPolarity_select,1)
        for d = 1:size(drugs_select,1)

            switch contrastPolarity_select(cc)
                case 5
                    points_x_ON = squeeze(mi_avgCombContr_ON(:,drugs_select(d),:));
                    points_x_OFF = squeeze(mi_avgCombContr_OFF(:,drugs_select(d),:));
                    
                otherwise
                    points_x_ON = squeeze(mi_avgSepContr_ON(:,idx_contrast_conds(contrastPolarity_select(cc,1),2),drugs_select(d),:));
                    points_x_OFF = squeeze(mi_avgSepContr_OFF(:,idx_contrast_conds(contrastPolarity_select(cc,1),1),drugs_select(d),:));
            end

            idx_respDrugs_ON = all(respConds_ON(drugs_select,:),1); idx_respDrugs_ON(isnan(idx_respDrugs_ON)) = 0;
            idx_respDrugs_ON = logical(idx_respDrugs_ON);

            idx_respDrugs_OFF = all(respConds_OFF(drugs_select,:),1); idx_respDrugs_OFF(isnan(idx_respDrugs_OFF)) = 0;
            idx_respDrugs_OFF = logical(idx_respDrugs_OFF);

            idx_toKeep_ON = repmat(idx_respDrugs_ON,size(points_x_ON,1),1) & ~isnan(points_x_ON) & repmat(idx_cells_expsSelect_ON',size(points_x_ON,1),1);
            idx_toKeep_OFF = repmat(idx_respDrugs_OFF,size(points_x_OFF,1),1) & ~isnan(points_x_OFF) & repmat(idx_cells_expsSelect_OFF',size(points_x_OFF,1),1);
            
            idx_toKeep_ON = idx_toKeep_ON & idx_valid_ON;
            idx_toKeep_OFF = idx_toKeep_OFF & idx_valid_OFF;

            plot_x_ON  = points_x_ON;
            plot_x_ON(~idx_toKeep_ON) = nan;
            std_plot_x_ON = nanstd(plot_x_ON,[],2);
            sem_plot_x_ON = std_plot_x_ON./sqrt(nansum(~isnan(plot_x_ON),2));
            rgb_num_ON = nansum(~isnan(plot_x_ON),2);
            mean_plot_x_ON = nanmedian(plot_x_ON,2);

            plot_x_OFF  = points_x_OFF;
            plot_x_OFF(~idx_toKeep_OFF) = nan;
            std_plot_x_OFF = nanstd(plot_x_OFF,[],2);
            sem_plot_x_OFF = std_plot_x_OFF./sqrt(nansum(~isnan(plot_x_OFF),2));
            rgb_num_OFF = nansum(~isnan(plot_x_OFF),2);
            mean_plot_x_OFF = nanmedian(plot_x_OFF,2);

            

            h_m(1) = subplot(1,1,plots_idx(1));
            hold on
            counter = counter+1;

            rgb = unique([flash_onset(timePoints),flash_onset(end)]);
            errorbar(rgb,[mean_plot_x_ON;0],[sem_plot_x_ON;0],'Color',col_on,'LineWidth',LINE_WIDTH,'LineStyle',LINE_TYPES{counter},'DisplayName',['ON(',num2str(max(rgb_num_ON)),') |',contrast_conds{contrastPolarity_select(cc)},' |',drugs{drugs_select(d)}]);
            errorbar(rgb,[mean_plot_x_OFF;0],[sem_plot_x_OFF;0],'Color',col_off,'LineWidth',LINE_WIDTH,'LineStyle',LINE_TYPES{counter},'DisplayName',['OFF(',num2str(max(rgb_num_OFF)),') |',contrast_conds{contrastPolarity_select(cc)},' |',drugs{drugs_select(d)}]);

            legend('-DynamicLegend');


            ylim([-1 0.4])
            ylabel('Modulation index')
            xlabel('Time from saccade offset (ms)')
            xlim(lim_x)
            
            for i = 1:size(plot_x_ON,1)
                temp = plot_x_ON(i,~isnan(plot_x_ON(i,:)));
                [p_pop_ON(i,cc,d)] = signrank(temp,0);
                pop_across_data_ON{i,cc,d} = temp;
                temp = plot_x_OFF(i,~isnan(plot_x_OFF(i,:)));
                [p_pop_OFF(i,cc,d)] = signrank(temp,0);
                pop_across_data_OFF{i,cc,d} = temp;
            end

        end
    end
    plot([-500,2500],[0,0],'--k')
    
    %
%     p_pop_ON = squeeze(p_pop_ON);
%     p_pop_OFF = squeeze(p_pop_OFF);
%     
%     a = [flash_onset(timePoints);p_pop_ON']
%     b = [flash_onset(timePoints);p_pop_OFF']
%     
%     pop_across_data_ON = squeeze(pop_across_data_ON);
%     pop_across_data_OFF = squeeze(pop_across_data_OFF);
%     
%     p_across_ON = cell2mat(cellfun(@(x,y) ranksum(x,y,'tail','both'),pop_across_data_ON(:,1),pop_across_data_ON(:,2),'UniformOutput',0));
%     p_across_OFF = cell2mat(cellfun(@(x,y) ranksum(x,y,'tail','both'),pop_across_data_OFF(:,1),pop_across_data_OFF(:,2),'UniformOutput',0));
%     p_across = [flash_onset(timePoints);p_across_ON';p_across_OFF']
%     
%     p_across_ON_OFF = cell2mat(cellfun(@(x,y) ranksum(x,y,'tail','both'),pop_across_data_ON(:,1),pop_across_data_OFF(:,1),'UniformOutput',0));
%     p_across_ON_OFF = [flash_onset(timePoints);p_across_ON_OFF']

%% Supp Fig. 8 Histograms

    plot_sig = 1;
    
    timePoints = 1:7;
    
    
    exps_select = [1:length(dates)];        % select which experiments data to use
    contrastPolarity_select = [1;2;5];  % {'neg'    'pos'    'preferred'    'anti-preferred'    'combined'}
    drugs_select = [1];   % drugs = {'none','gabazine+ptx+str'};
    
    contrast_conds = {'neg','pos','preferred','anti-preferred','combined'};
    idx_contrast_conds = [1,1;2,2;2,1;1,2]; % [off,on] rows are the contrast conditions [neg,pos,preferred,antipreferred]
    
    plots_idx = [1:(length(timePoints)*size(drugs_select,1)*size(contrastPolarity_select,1))]';
    plots_idx = reshape(plots_idx',length(timePoints),size(drugs_select,1),size(contrastPolarity_select,1));
%     plots_idx = reshape(plots_idx',1,(length(timePoints)*length(contrast_conds)));

    lim_x = [-1.2,1.2];
    lim_y = [0,0.5];
    bins = [-1:0.1:1];
    
    d = 1; cc = 1; t = 1; i = 1;
    
    idx_cells_expsSelect_ON = ismember(unitExpIds_ON,exps_select);
    idx_cells_expsSelect_OFF = ismember(unitExpIds_OFF,exps_select);
    
    h_main=figure; suptitle(['Binned contrasts - ',contrast_conds{contrastPolarity_select(cc)}])    
    h_m = [];
    f_on_all = [];
    f_off_all = [];
    alpha_p  = 0.05;
    pow_thresh = 80;
    num_onCells = nan(1,length(flash_onset));
    num_offCells = nan(1,length(flash_onset));
    num_onCells_sig = nan(1,length(flash_onset));
    num_offCells_sig = nan(1,length(flash_onset));


    for cc = 1:size(contrastPolarity_select,1)
        for d = 1:size(drugs_select,1)
            for t = timePoints;


                rgb_on = [];
                rgb_off = [];
                switch contrast_conds{contrastPolarity_select(cc)}
                    case 'combined'
                        rgb_on = squeeze(mi_avgCombContr_ON(t,drugs_select(d),:));
                        rgb_off = squeeze(mi_avgCombContr_OFF(t,drugs_select(d),:));
                        
                        sig_cells_ON = logical(squeeze(p_avgCombContr_ON(t,drugs_select(d),:)<alpha_p)) & logical(squeeze(pow_p_avgCombContr_ON(t,drugs_select(d),:)>=pow_thresh));
                        sig_cells_OFF = logical(squeeze(p_avgCombContr_OFF(t,drugs_select(d),:)<alpha_p)) & logical(squeeze(pow_p_avgCombContr_OFF(t,drugs_select(d),:)>=pow_thresh));


                    otherwise
                        rgb_on = squeeze(mi_avgSepContr_ON(t,idx_contrast_conds(contrastPolarity_select(cc,1),2),drugs_select(d),:));
                        rgb_off = squeeze(mi_avgSepContr_OFF(t,idx_contrast_conds(contrastPolarity_select(cc,1),1),drugs_select(d),:));
                        
                        sig_cells_ON = logical(squeeze(p_avgSepContr_ON(t,idx_contrast_conds(contrastPolarity_select(cc,1),2),drugs_select(d),:)<alpha_p)) & logical(squeeze(pow_p_avgSepContr_ON(t,idx_contrast_conds(contrastPolarity_select(cc,1),2),drugs_select(d),:)>=pow_thresh));
                        sig_cells_OFF = logical(squeeze(p_avgSepContr_OFF(t,idx_contrast_conds(contrastPolarity_select(cc,1),2),drugs_select(d),:)<alpha_p)) & logical(squeeze(pow_p_avgSepContr_OFF(t,idx_contrast_conds(contrastPolarity_select(cc,1),2),drugs_select(d),:)>=pow_thresh));
                end

                idx_respDrugs_ON = respConds_ON(drugs_select(d),:); idx_respDrugs_ON(isnan(idx_respDrugs_ON)) = 0;
                idx_respDrugs_ON = logical(idx_respDrugs_ON);

                idx_respDrugs_OFF = respConds_OFF(drugs_select(d),:); idx_respDrugs_OFF(isnan(idx_respDrugs_OFF)) = 0;
                idx_respDrugs_OFF = logical(idx_respDrugs_OFF);

                idx_toKeep_ON = idx_respDrugs_ON & ~isnan(rgb_on') & idx_cells_expsSelect_ON';
                idx_toKeep_OFF = idx_respDrugs_OFF & ~isnan(rgb_off') & idx_cells_expsSelect_OFF';

                plot_x_ON  = rgb_on(idx_toKeep_ON);
                plot_x_OFF  = rgb_off(idx_toKeep_OFF);
                
                plot_x_ON_sig = rgb_on(idx_toKeep_ON & sig_cells_ON');
                plot_x_OFF_sig = rgb_off(idx_toKeep_OFF & sig_cells_OFF');

                h_m(t,d,cc) = subplot(size(contrastPolarity_select,1)*size(drugs_select,1),length(timePoints),plots_idx(t,d,cc)); hold on

                [f_on_all,x] = hist(plot_x_ON,bins);
                h = bar(x,f_on_all/trapz(f_on_all),'histc');
%                 h = histogram(plot_x_ON,bins,'Normalization','probability');
                set(h,'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
                
                if plot_sig == 1
                    [f_on_sig,x] = hist(plot_x_ON_sig,bins);
                    h = bar(x,f_on_sig/trapz(f_on_all),'histc');
                    set(h,'FaceColor',[0 0.6 0],'EdgeColor',[0 0.6 0]);
                end
                

                [f_off_all,x] = hist(plot_x_OFF,bins);
                h = bar(x,f_off_all/trapz(f_off_all),'histc');
%                 h = histogram(plot_x_OFF,bins,'Normalization','probability');
                set(h,'FaceColor',[0.1,0.1,0.1],'EdgeColor',[0.1 0.1 0.1]);
%                 h.FaceAlpha = 0.6;
                if plot_sig == 1
                    [f_off_sig,x] = hist(plot_x_OFF_sig,bins);
                    h = bar(x,f_off_sig/trapz(f_off_all),'histc');
                     set(h,'FaceColor',[0.6 0 0],'EdgeColor',[0.6 0 0]);
                end



                plot([0,0],[0,1],'--k')

                hold off
                xlim(lim_x)
                ylim(lim_y)

                num_onCells(cc,t,d) = sum(~isnan(plot_x_ON));
                num_offCells(cc,t,d) = sum(~isnan(plot_x_OFF));
                
                num_onCells_sig(cc,t,d) = sum(~isnan(plot_x_ON_sig));
                num_offCells_sig(cc,t,d) = sum(~isnan(plot_x_OFF_sig));
                

                if plot_sig == 1
                    h_l = legend({['ON: ',num2str(sum(~isnan(plot_x_ON)))],...
                                  ['ON\_s: ',num2str(round(sum(~isnan(plot_x_ON_sig))/sum(~isnan(plot_x_ON))*100)),'%'],...
                                  ['OFF: ',num2str(sum(~isnan(plot_x_OFF)))],...
                                  ['OFF\_s: ',num2str(round(sum(~isnan(plot_x_OFF_sig))/sum(~isnan(plot_x_OFF))*100)),'%']});
                else
                    h_l = legend({['ON: ',num2str(sum(~isnan(plot_x_ON)))],...
                                  ['OFF: ',num2str(sum(~isnan(plot_x_OFF)))]});
                end
                
                xlabel([drugs{drugs_select(d)}])

            end
        end
    end
    
    

        % 
    for t = timePoints
        title(h_m(t,1,1),[num2str(flash_onset(t)),'ms'])
    end

    for cc = 1:size(contrastPolarity_select,1)
        for d = 1:size(drugs_select,1)
            ylabel(h_m(1,d,cc),['contrast: ',contrast_conds(contrastPolarity_select(cc))])
            xlabel(h_m(1,d,cc),['drugs: ',drugs{drugs_select(d)}])
        end
    end
    
 %% Supp Fig. 10 (a=bottom row; b = top row)
    timePoints = 1:7;
    exps_select = 1:length(dates);        % select which experiments data to use
    contrastStep_select = [1;2];        % contrast_conds = {'neg','pos','preferred','anti-preferred','combined'};
    drugs_select = [1];   % drugs = {'none','gabazine+ptx+str'};
    
    idx_contrast_conds = [1,1;2,2;2,1;1,2]; % [on,off] rows are the conditions
    contrast_conds = {'neg','pos','preferred','anti-preferred','combined'};
    
    plots_idx = [1:(length(timePoints)*size(contrastStep_select,1)*size(drugs_select,1))]';
    plots_idx = reshape(plots_idx',length(timePoints),size(drugs_select,1),size(contrastStep_select,1));

    lim_x = [-1.2,1.2];
    lim_y = [-1.2,1.2];
    bins = [-1:0.1:1];
    baseVal = -1.2;
    p_sig = 0.05;
    pow_thresh = 80;
    
    col_on = [0.6,0.6,0.6];
    col_off = [0.2,0.2,0.2];
    
    d = 1; cc = 1; t = 1; i = 1;
        
    h_main=figure; suptitle(['Association index based on all steps'])    
    h_m = [];
    f_on_all = [];
    f_off_all = [];
    
    num_onCells = nan(length(flash_onset));
    num_offCells = nan(length(flash_onset));
    
    idx_cells_expsSelect_ON = ismember(unitExpIds_ON,exps_select);
    idx_cells_expsSelect_OFF = ismember(unitExpIds_OFF,exps_select);

        
    for cc = 1:size(contrastStep_select,1)

            for d = 1:size(drugs_select,1)
                    for t = timePoints;


                        points_x_ON = [];
                        points_y_ON = [];
                        points_x_OFF = [];
                        points_y_OFF = [];
                        
                        switch contrast_conds{contrastStep_select(cc)}
                            case 'combined'
%                                 points_x_ON = squeeze(mi_allCombContr_ON(t,contrast_Nbins(i,1),drugs_select(d,1),:));
%                                 points_y_ON = squeeze(mi_allCombContr_ON(t,contrast_Nbins(i,2),drugs_select(d,2),:));
%                                 
%                                 points_x_OFF = squeeze(mi_allCombContr_OFF(t,contrast_Nbins(i,1),drugs_select(d,1),:));
%                                 points_y_OFF = squeeze(mi_allCombContr_OFF(t,contrast_Nbins(i,2),drugs_select(d,2),:));
                                
                            otherwise
                                points_x_ON = squeeze(mi_avgSepContr_ON(t,idx_contrast_conds(contrastStep_select(cc),2),drugs_select(d),:));
                                points_y_ON = squeeze(assoc_R_sepContr_ON(t,idx_contrast_conds(contrastStep_select(cc),2),drugs_select(d),:));
                                
                                points_x_OFF = squeeze(mi_avgSepContr_OFF(t,idx_contrast_conds(contrastStep_select(cc),1),drugs_select(d),:));
                                points_y_OFF = squeeze(assoc_R_sepContr_OFF(t,idx_contrast_conds(contrastStep_select(cc),1),drugs_select(d),:));

                        end

        %                 rgb_on = rgb_on(idx_cells_lowTrials_ON);    
        %                 rgb_off = rgb_off(idx_cells_lowTrials_OFF);
                        idx_respDrugs_ON = respConds_ON(drugs_select(d),:); idx_respDrugs_ON(isnan(idx_respDrugs_ON)) = 0;
                        idx_respDrugs_ON = logical(idx_respDrugs_ON);

                        idx_respDrugs_OFF = respConds_OFF(drugs_select(d),:); idx_respDrugs_OFF(isnan(idx_respDrugs_OFF)) = 0;
                        idx_respDrugs_OFF = logical(idx_respDrugs_OFF);
                        
                        idx_sig_supp_ON = (squeeze(p_avgSepContr_ON(t,idx_contrast_conds(contrastStep_select(cc),2),drugs_select(d),:)) < p_sig) & (squeeze(pow_p_avgSepContr_ON(t,idx_contrast_conds(contrastStep_select(cc),2),drugs_select(d),:)) > pow_thresh);
                        idx_sig_supp_OFF = (squeeze(p_avgSepContr_OFF(t,idx_contrast_conds(contrastStep_select(cc),2),drugs_select(d),:)) < p_sig) & (squeeze(pow_p_avgSepContr_OFF(t,idx_contrast_conds(contrastStep_select(cc),2),drugs_select(d),:)) > pow_thresh);
                        
                        idx_sig_assoc_ON = squeeze(assoc_p_sepContr_ON(t,idx_contrast_conds(contrastStep_select(cc),2),drugs_select(d),:)) < p_sig;
                        idx_sig_assoc_OFF = squeeze(assoc_p_sepContr_OFF(t,idx_contrast_conds(contrastStep_select(cc),2),drugs_select(d),:)) < p_sig;
                        
                        
                        idx_toKeep_ON = idx_respDrugs_ON & ~isnan(points_x_ON') & ~isnan(points_y_ON') & idx_cells_expsSelect_ON';% & idx_sig_supp_ON' & idx_sig_assoc_ON';
                        idx_toKeep_OFF = idx_respDrugs_OFF & ~isnan(points_x_OFF') & ~isnan(points_y_OFF') & idx_cells_expsSelect_OFF';% & idx_sig_supp_OFF' & idx_sig_assoc_OFF';
                        
                        plot_x_ON  = points_x_ON(idx_toKeep_ON);
                        plot_y_ON  = points_y_ON(idx_toKeep_ON);
                        plot_x_OFF  = points_x_OFF(idx_toKeep_OFF);
                        plot_y_OFF  = points_y_OFF(idx_toKeep_OFF);
                        
%                         plot_x_OFF = [];        % to only plot ON
%                         plot_y_OFF = [];
% 
%                         plot_x_ON = [];        % to only plot OFF
%                         plot_y_ON = [];

                        h_m(t,d,cc) = subplot(size(contrastStep_select,1)*size(drugs_select,1),length(timePoints),plots_idx(t,d,cc));
                        hold on

                        scatter(plot_x_ON,plot_y_ON,'MarkerEdgeColor','g','MarkerFaceColor',col_on+.05) 
                        scatter(plot_x_OFF,plot_y_OFF,'MarkerEdgeColor','g','MarkerFaceColor',col_off+.05)
%                         
%                         scatter(points_x_ON(idx_toKeep_ON & idx_sig_supp_ON'),points_y_ON(idx_toKeep_ON & idx_sig_supp_ON'),'MarkerEdgeColor',[1,0,0],'MarkerFaceColor',col_on+.05) 
%                         scatter(points_x_OFF(idx_toKeep_OFF & idx_sig_supp_OFF'),points_y_OFF(idx_toKeep_OFF & idx_sig_supp_OFF'),'MarkerEdgeColor',[0.9,0,0],'MarkerFaceColor',col_off+.05) 
                        
                        scatter(points_x_ON(idx_toKeep_ON & idx_sig_supp_ON' & idx_sig_assoc_ON'),points_y_ON(idx_toKeep_ON & idx_sig_supp_ON' & idx_sig_assoc_ON'),'MarkerEdgeColor',[0,0,1],'MarkerFaceColor',col_on+.05) 
                        scatter(points_x_OFF(idx_toKeep_OFF & idx_sig_supp_OFF' & idx_sig_assoc_OFF'),points_y_OFF(idx_toKeep_OFF & idx_sig_supp_OFF' & idx_sig_assoc_OFF'),'MarkerEdgeColor',[0,0,0.9],'MarkerFaceColor',col_off+.05) 

                        [f_off_all,x] = hist(plot_y_OFF,bins);
                        h = barh(x,(f_off_all/trapz(f_off_all))+baseVal);
                        set(h,'FaceColor',col_off,'EdgeColor','k');
                        h.BaseValue = baseVal;

                        [f_on_all,x] = hist(plot_y_ON,bins);
                        h = barh(x,(f_on_all/trapz(f_on_all))+baseVal);
                        set(h,'FaceColor',col_on,'EdgeColor','k');
                        h.BaseValue = baseVal;
                        set(h,'FaceAlpha',0.6)

                        [f_off_all,x] = hist(plot_x_OFF,bins);
                        h = bar(x,(f_off_all/trapz(f_off_all))+baseVal);
                        set(h,'FaceColor',col_off,'EdgeColor','k');
                        h.BaseValue = baseVal;

                        [f_on_all,x] = hist(plot_x_ON,bins);
                        h = bar(x,(f_on_all/trapz(f_on_all))+baseVal);
                        set(h,'FaceColor',col_on,'EdgeColor','k');
                        h.BaseValue = baseVal;
                        set(h,'FaceAlpha',0.6)

                        plot([-1,1],[0,0],'--','color',[0.8,0.8,0.8])
                        plot([0,0],[-1,1],'--','color',[0.8,0.8,0.8])
                        plot([-1,1],[-1,1],'k')
                        xlim(lim_x)
                        ylim(lim_y)

                        axis square

                        hold off
                        num_onCells(t,d,cc) = sum(~isnan(plot_x_ON));
                        num_offCells(t,d,cc) = sum(~isnan(plot_x_OFF));

                        h_l = legend({['ON: ',num2str(sum(~isnan(plot_x_ON)))],['OFF: ',num2str(sum(~isnan(plot_x_OFF)))]});

                    end
                    
                    xlabel(['MI | step: ',contrast_conds{contrastStep_select(cc)},' | ',drugs{drugs_select(d)}])
                    ylabel(['AI | step: ',contrast_conds{contrastStep_select(cc)},' | ',drugs{drugs_select(d)}])
            end
            
    end

    for t = timePoints
        title(h_m(t,1,1),[num2str(flash_onset(t)),'ms'])
    end
        
    

