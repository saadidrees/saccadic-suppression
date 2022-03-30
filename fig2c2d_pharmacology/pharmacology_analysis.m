
%% Figures 2c, 2d and Supp Fig. 6
% loads data from all full-field experiments where pharmacological agents were also used. 
% The main data used is the peak response of flashes at different delays after saccade
% and for different conditions. in general, data is organized as [delays,sc,drugs,cells].

%%
clear all; clc; set(0,'defaultTextInterpreter','none');

load('data_pharmacology.mat')

sc = [25,150];
flash_onset = [-265,17,33,50,100,250,500,1000,1990];        % wrt saccade offset
drugs = {'none','gabazine','ptx','str','gabazine+ptx','gabazine+ptx+str'};
numSaccs = 40;
dataLen_pre = 400;
dataLen_post = 3600;
saccOnset = 2000;
idx_data = saccOnset-dataLen_pre+1:saccOnset+dataLen_post;

baseFlash_idx = length(flash_onset);      % take the 9th flash i.e. at 1867ms as the baseline flash


%% Calculate modulation indices

% Modulation Index - refined peaks - avg
    mi_avg_flash_ON = (peaks_ref_avg_flash_ON - peaks_ref_avg_baseFlash_ON)./(peaks_ref_avg_flash_ON + peaks_ref_avg_baseFlash_ON);     % [delays,sc,drugs,cells]
    mi_avg_flash_OFF = (peaks_ref_avg_flash_OFF - peaks_ref_avg_baseFlash_OFF)./(peaks_ref_avg_flash_OFF + peaks_ref_avg_baseFlash_OFF);
    
% Modulation Index - refined peaks - all saccs
    mi_allSacc_flash_ON = (peaks_ref_allSacc_flash_ON - peaks_ref_allSacc_baseFlash_ON)./(peaks_ref_allSacc_flash_ON + peaks_ref_allSacc_baseFlash_ON);
    mi_allSacc_flash_OFF = (peaks_ref_allSacc_flash_OFF - peaks_ref_allSacc_baseFlash_OFF)./(peaks_ref_allSacc_flash_OFF + peaks_ref_allSacc_baseFlash_OFF);

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
    
    p_avg_flash_ON = nan(length(flash_onset)-1,length(sc),length(drugs),length(uname_ON));
    pow_p_avg_flash_ON = nan(length(flash_onset)-1,length(sc),length(drugs),length(uname_ON));
    sampNum_p_avg_flash_ON = nan(length(flash_onset)-1,length(sc),length(drugs),length(uname_ON));
    num_sacc_ON = [];
    
    
    p_avg_flash_OFF = nan(length(flash_onset)-1,length(sc),length(drugs),length(uname_OFF));
    pow_p_avg_flash_OFF = nan(length(flash_onset)-1,length(sc),length(drugs),length(uname_OFF));
    sampNum_p_avg_flash_off = nan(length(flash_onset)-1,length(sc),length(drugs),length(uname_OFF));
    num_sacc_OFF = [];
    
    idx_onCells = false(length(flash_onset)-1,length(sc),length(drugs),length(uname_ON));
    idx_offCells = false(length(flash_onset)-1,length(sc),length(drugs),length(uname_OFF));
    
    for l = 1:size(mi_allSacc_flash_OFF,5)      % cells
        for k = 1:size(mi_allSacc_flash_OFF,4)      % drugs
            for j = 1:size(mi_allSacc_flash_OFF,3)      % sc
                for i = 1:size(mi_allSacc_flash_OFF,1)  % time point
                    rgb = squeeze(mi_allSacc_flash_OFF(i,:,j,k,l));

                    idx_toRemove = isnan(rgb) | rgb==inf | rgb==-inf;
                    rgb(idx_toRemove) = [];
                    num_sacc_OFF(i,j,k,l) = length(rgb);
                     if ~isempty(rgb)
                        idx_offCells(i,j,k,l) = true;
                        if mi_avg_flash_OFF(i,j,k,l) < 0
                            [p_avg_flash_OFF(i,j,k,l),~,stats] = signtest(rgb,0,'tail','left');
                            pow_p_avg_flash_OFF(i,j,k,l) = round(100*sampsizepwr('p',0.50,(stats.sign/num_sacc_OFF(i,j,k,l))+eps,[],num_sacc_OFF(i,j,k,l),'tail','left'));        
                        else 
                            [p_avg_flash_OFF(i,j,k,l),~,stats] = signtest(rgb,0,'tail','right');
                            pow_p_avg_flash_OFF(i,j,k,l) = round(100*sampsizepwr('p',0.50,(stats.sign/num_sacc_OFF(i,j,k,l))-eps,[],num_sacc_OFF(i,j,k,l),'tail','right'));   % 1-binocdf(critcL,num_sacc_OFF(i,j,k),stats.sign/num_sacc_OFF(i,j,k)) where critclL = 25;

                        end

                     end 
                end
            end
        end
    end
    

    for l = 1:size(mi_allSacc_flash_ON,5)      % cells
        for k = 1:size(mi_allSacc_flash_ON,4)      % drugs
            for j = 1:size(mi_allSacc_flash_ON,3)   % sc
                for i = 1:size(mi_allSacc_flash_ON,1)   % delays
                    rgb = squeeze(mi_allSacc_flash_ON(i,:,j,k,l));
                    idx_toRemove = isnan(rgb) | rgb==inf | rgb==-inf;
                    rgb(idx_toRemove) = [];
                    num_sacc_ON(i,j,k,l) = length(rgb);

                    if ~isempty(rgb)
                        idx_onCells(i,j,k,l) = true;
                        if mi_avg_flash_ON(i,j,k,l) < 0
                            [p_avg_flash_ON(i,j,k,l),~,stats] = signtest(rgb,0,'tail','left');
                            pow_p_avg_flash_ON(i,j,k,l) = round(100*sampsizepwr('p',0.50,(stats.sign/num_sacc_ON(i,j,k,l))+eps,[],num_sacc_ON(i,j,k,l),'tail','left'));
                        else 
                            [p_avg_flash_ON(i,j,k,l),~,stats] = signtest(rgb,0,'tail','right');
                            pow_p_avg_flash_ON(i,j,k,l) = round(100*sampsizepwr('p',0.50,(stats.sign/num_sacc_ON(i,j,k,l))-eps,[],num_sacc_ON(i,j,k,l),'tail','right'));
                        end

    %                     sampNum_p_avg_flash_on(i,j,k) = sampsizepwr('p',0.50,(stats.sign/num_sacc_on(i,j,k))+eps,0.8);
    % 
                    end

                end
            end
        end
    end
  

%% Population significance modulation for specific condition
time_select = [2,4:8];   % flash_onset = [-265,17,33,50,100,250,500,1000,1990];        % wrt saccade offset
drugs_select = [1,6];       % drugs = {'none','gabazine','ptx','str','gabazine+ptx','gabazine+ptx+str'};
sc_select = [2];        % [25,150]
same_cells = 1;

if same_cells == 1
    rgb = respConds_ON; rgb(isnan(rgb)) = 0;
    idx_sameCells_ON = all(rgb(drugs_select,:));
    rgb = respConds_OFF; rgb(isnan(rgb)) = 0;
    idx_sameCells_OFF = all(rgb(drugs_select,:));
else
    idx_sameCells_ON = logical(ones(1,size(mi_avg_flash_ON,4)));
    idx_sameCells_OFF = logical(ones(1,size(mi_avg_flash_OFF,4)));
end

p_pop_ON = nan(length(drugs_select),length(time_select));
p_pop_OFF = nan(length(drugs_select),length(time_select));


for i = 1:length(time_select)
    for j = 1:length(drugs_select)
        try
            [p_pop_ON(j,i)] = signrank(squeeze(mi_avg_flash_ON(time_select(i),sc_select,drugs_select(j),idx_sameCells_ON)),0,'tail','both');
            [p_pop_OFF(j,i)] = signrank(squeeze(mi_avg_flash_OFF(time_select(i),sc_select,drugs_select(j),idx_sameCells_OFF)),0,'tail','both');
        end
    end
end

p_sig_supp_ON = [flash_onset(time_select);p_pop_ON];
p_sig_supp_OFF = [flash_onset(time_select);p_pop_OFF];


%% population significance modulation for specific condition with and without drugs 
time_select = [2,4:8];   % flash_onset = [-265,17,33,50,100,250,500,1000,1990];        % wrt saccade offset
drugs_select = [1,6];       % drugs = {'none','gabazine','ptx','str','gabazine+ptx','gabazine+ptx+str'};
sc_select = [2];        % [25,150]

rgb = respConds_ON; rgb(isnan(rgb)) = 0;
idx_sameCells_ON = all(rgb(drugs_select,:));
rgb = respConds_OFF; rgb(isnan(rgb)) = 0;
idx_sameCells_OFF = all(rgb(drugs_select,:));

p_pop_ON = nan(1,length(time_select));
p_pop_OFF = nan(1,length(time_select));


for i = 1:length(time_select)
    try
        [p_pop_ON(i)] = ranksum(squeeze(mi_avg_flash_ON(time_select(i),sc_select,drugs_select(1),idx_sameCells_ON)),squeeze(mi_avg_flash_ON(time_select(i),sc_select,drugs_select(2),idx_sameCells_ON)),'tail','both');
        [p_pop_OFF(i)] = ranksum(squeeze(mi_avg_flash_OFF(time_select(i),sc_select,drugs_select(1),idx_sameCells_OFF)),squeeze(mi_avg_flash_OFF(time_select(i),sc_select,drugs_select(2),idx_sameCells_OFF)),'tail','both');
    end
end

p_sig_acrossDrugs = [flash_onset(time_select);p_pop_ON;p_pop_OFF];

%% Fig. 2c and Fig. 2d - Line plots of modulation index

    select_fig = '2d'; %{'2c','2d'} select one
    
    switch select_fig
        case '2c'
            drugs_select = [1;5];       % drugs = {'none','gabazine','ptx','str','gabazine+ptx','gabazine+ptx+str'};
        case '2d'
            drugs_select = [1;6];       % drugs = {'none','gabazine','ptx','str','gabazine+ptx','gabazine+ptx+str'};
    end

    exps_select = [1:length(dates)];
    time_select = [2,4:8];   % flash_onset = [-265,17,33,50,100,250,500,1000,1990];        % wrt saccade offset
    sc_select = [2];
    same_cells = 1;

    idx_cells_expsSelect_ON = ismember(unitExpIds_ON,exps_select);
    idx_cells_expsSelect_OFF = ismember(unitExpIds_OFF,exps_select);
    sc_text = cellstr(num2str(sc'));

   
    plots_idx = 1;
    
    lim_x = [-200,2100];
    bins = [-1:0.1:1];
    
    col_on = 'r';%[0.6,0.6,0.6];
    col_off = 'b';%[0.2,0.2,0.2];
    LINE_WIDTH = 2;
    LINE_TYPES = {'-','--',':','-.'};

    h_main=figure; suptitle(['Pharmacology'])    
    h_m = [];
    
    s = 1; d = 1; t = 1;
    counter = 0;
    
    for s = 1:size(sc_select,1)
        for d = 1:size(drugs_select,1)

                points_x_ON = [];     
                for j = 1:size(mi_avg_flash_ON,4)
                    points_x_ON(j,:) = squeeze(mi_avg_flash_ON(time_select,sc_select(s),drugs_select(d),j));
                end

                points_x_OFF = [];     
                for j = 1:size(mi_avg_flash_OFF,4)
                    points_x_OFF(j,:) = squeeze(mi_avg_flash_OFF(time_select,sc_select(s),drugs_select(d),j));
                end


                idx_respDrugs_ON = respConds_ON(drugs_select(d),:); idx_respDrugs_ON(isnan(idx_respDrugs_ON)) = 0;
                idx_respDrugs_ON = logical(idx_respDrugs_ON);

                idx_respDrugs_OFF = respConds_OFF(drugs_select(d),:); idx_respDrugs_OFF(isnan(idx_respDrugs_OFF)) = 0;
                idx_respDrugs_OFF = logical(idx_respDrugs_OFF);
                
                rgb = respConds_ON; rgb(isnan(rgb)) = 0;
                idx_sameCells_ON = all(rgb(drugs_select,:));
                rgb = respConds_OFF; rgb(isnan(rgb)) = 0;
                idx_sameCells_OFF = all(rgb(drugs_select,:));

                idx_toKeep_ON = repmat(idx_respDrugs_ON,size(points_x_ON,2),1) & ~isnan(points_x_ON') & repmat(idx_cells_expsSelect_ON',size(points_x_ON,2),1) & repmat(idx_sameCells_ON,size(points_x_ON,2),1);
                idx_toKeep_OFF = repmat(idx_respDrugs_OFF,size(points_x_OFF,2),1) & ~isnan(points_x_OFF') & repmat(idx_cells_expsSelect_OFF',size(points_x_OFF,2),1) & repmat(idx_sameCells_OFF,size(points_x_OFF,2),1);

                plot_x_ON  = points_x_ON;
                plot_x_ON(~idx_toKeep_ON') = nan;
                std_plot_x_ON = nanstd(plot_x_ON,[],1);
                sem_plot_x_ON = std_plot_x_ON./sqrt(nansum(~isnan(plot_x_ON),1));
                rgb_num_ON = nansum(~isnan(plot_x_ON),1);
                mean_plot_x_ON = nanmean(plot_x_ON,1);

                plot_x_OFF  = points_x_OFF;
                plot_x_OFF(~idx_toKeep_OFF') = nan;
                std_plot_x_OFF = nanstd(plot_x_OFF,[],1);
                sem_plot_x_OFF = std_plot_x_OFF./sqrt(nansum(~isnan(plot_x_OFF),1));
                rgb_num_OFF = nansum(~isnan(plot_x_OFF),1);
                mean_plot_x_OFF = nanmean(plot_x_OFF,1);


                h_m(1) = subplot(1,1,plots_idx(1));
                hold on
                counter = counter+1;

                rgb = unique([flash_onset(time_select),flash_onset(end)]);
                errorbar(rgb,[mean_plot_x_ON,0],[sem_plot_x_ON,0],'Color',col_on,'LineWidth',LINE_WIDTH,'LineStyle',LINE_TYPES{counter},'DisplayName',['ON (',num2str(nanmax(rgb_num_ON)),') |',sc_text{sc_select(s)},'µm |',drugs{drugs_select(d)}]);
                errorbar(rgb,[mean_plot_x_OFF,0],[sem_plot_x_OFF,0],'Color',col_off,'LineWidth',LINE_WIDTH,'LineStyle',LINE_TYPES{counter},'DisplayName',['OFF (',num2str(nanmax(rgb_num_OFF)),') |',sc_text{sc_select(s)},'µm |',drugs{drugs_select(d)}]);

                legend('-DynamicLegend');




                ylim([-0.7 0.3])
                ylabel('Modulation index')
                xlabel('Time from saccade offset (ms)')
                xlim(lim_x)
                title(['ON = ',num2str(nanmax(rgb_num_ON)),' | OFF = ',num2str(nanmax(rgb_num_OFF))])
                

        end
    end
    plot([-500,2500],[0,0],'--k')
    

%% Supp Fig. S6 a(top row) b(bottom row) - Scatter Plots

    time_select = [2,4:8];   % flash_onset = [-265,17,33,50,100,250,500,1000,1990];        % wrt saccade offset
    drugs_select = [1,5;1,6];       % drugs = {'none','gabazine','ptx','str','gabazine+ptx','gabazine+ptx+str'};
    sc_select = [2,2];
    exps_select = 1:length(dates);
    
    sc_text = cellstr(num2str(sc'));

    plots_idx = [1:(length(time_select)*size(sc_select,1)*size(drugs_select,1))]';
    plots_idx = reshape(plots_idx',length(time_select),size(sc_select,1),size(drugs_select,1));        % [time,mask_conds,sc]
    
    lim_x = [-1.2,1.2];
    lim_y = [-1.2,1.2];
    bins = [-1:0.1:1];
    baseVal = -1.2;
    
    col_on = [0.6,0.6,0.6];
    col_off = [0.2,0.2,0.2];
    
    idx_cells_expsSelect_ON = ismember(unitExpIds_ON,exps_select);
    idx_cells_expsSelect_OFF = ismember(unitExpIds_OFF,exps_select);

    h_main=figure; suptitle(['Scatter Plots - pharmacology'])    
    h_m = [];
    
    s = 1; m = 1; d = 1; t = 1;

    for d = 1:size(drugs_select,1)
        for s = 1:size(sc_select,1)
            for t = 1:size(time_select,2)


                points_x_ON = squeeze(mi_avg_flash_ON(time_select(t),sc_select(s,1),drugs_select(d,1),:));
                points_y_ON = squeeze(mi_avg_flash_ON(time_select(t),sc_select(s,2),drugs_select(d,2),:));

                points_x_OFF = squeeze(mi_avg_flash_OFF(time_select(t),sc_select(s,1),drugs_select(d,1),:));
                points_y_OFF = squeeze(mi_avg_flash_OFF(time_select(t),sc_select(s,2),drugs_select(d,2),:));


                h_m(t,s,d) = subplot(size(sc_select,1)*size(drugs_select,1),length(time_select),plots_idx(t,s,d));
                hold on

                idx_respDrugs_ON = respConds_ON(drugs_select(d,:),:); idx_respDrugs_ON(isnan(idx_respDrugs_ON)) = 0;
                idx_respDrugs_ON = all(logical(idx_respDrugs_ON),1);

                idx_respDrugs_OFF = respConds_OFF(drugs_select(d,:),:); idx_respDrugs_OFF(isnan(idx_respDrugs_OFF)) = 0;
                idx_respDrugs_OFF = all(logical(idx_respDrugs_OFF));

                idx_toKeep_ON = idx_respDrugs_ON & ~isnan(points_x_ON') & ~isnan(points_y_ON') & idx_cells_expsSelect_ON';
                idx_toKeep_OFF = idx_respDrugs_OFF & ~isnan(points_x_OFF') & ~isnan(points_y_OFF') & idx_cells_expsSelect_OFF';

                plot_x_ON  = points_x_ON(idx_toKeep_ON);
                plot_y_ON  = points_y_ON(idx_toKeep_ON);

                plot_x_OFF  = points_x_OFF(idx_toKeep_OFF);
                plot_y_OFF  = points_y_OFF(idx_toKeep_OFF);

                plots_uname_ON = uname_ON(idx_toKeep_ON');
                plots_uname_OFF = uname_OFF(idx_toKeep_OFF');

                scatter(plot_x_ON,plot_y_ON,'MarkerEdgeColor',col_on,'MarkerFaceColor',col_on) 
                scatter(plot_x_OFF,plot_y_OFF,'MarkerEdgeColor',col_off,'MarkerFaceColor',col_off)

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
                legend({['ON ',num2str(length(plot_x_ON))],['OFF: ',num2str(length(plot_x_OFF))]})
            end
            

            xlabel(['sc: ',sc_text{sc_select(s,1)},' µm | ', drugs{drugs_select(d,1)}])
            ylabel(['sc: ',sc_text{sc_select(s,2)},' µm | ',drugs{drugs_select(d,2)}])
  
        end
    end
    
    for i = 1:size(time_select,2)
        title(h_m(i,1,1,1),[num2str(flash_onset(time_select(i))),'ms'])
    end
    


