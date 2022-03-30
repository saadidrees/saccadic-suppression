
%% Fig. 9a, c, d
%%
clear all; clc; set(0,'defaultTextInterpreter','none');

load data_nonpreferredFlashes.mat

delays = [-233,-217,-184,-167,-150,-133,-117,-50,17,50,100,250,500,1000,2000];
flash_onset = delays+100;
baseFlash_idx = length(flash_onset);
sc = [25,50,150,300];
flash_polarity = {'dark flash','bright flash'};
numSaccs = 40;
dataLen_pre = 400;
dataLen_post = 3600;
saccOnset = 2000;
idx_data = saccOnset-dataLen_pre+1:saccOnset+dataLen_post;



%% Calculate modulation indices

% Modulation Index - refined peaks - avg - TEMPLATE MATCHING
    mi_tm_avg_flash_ON = (sf_tm_avg_flash_ON - 1)./(sf_tm_avg_flash_ON + 1);     % [delays,sc,drugs,cells]
    mi_tm_avg_flash_OFF = (sf_tm_avg_flash_OFF - 1)./(sf_tm_avg_flash_OFF + 1);     % [delays,sc,drugs,cells]
    
    mi_tm_avg_flash_ON(mi_tm_avg_flash_ON<-1 | mi_tm_avg_flash_ON>1) = nan;
    mi_tm_avg_flash_OFF(mi_tm_avg_flash_OFF<-1 | mi_tm_avg_flash_OFF>1) = nan;
    
    faulty_cells_OFF = 19;
    mi_tm_avg_flash_OFF(:,:,:,faulty_cells_OFF) = nan;


% Modulation Index - refined peaks - all saccs
    mi_allSacc_flash_ON = (peaks_ref_allSacc_flash_ON - peaks_ref_allSacc_baseFlash_ON)./(peaks_ref_allSacc_flash_ON + peaks_ref_allSacc_baseFlash_ON);
    mi_allSacc_flash_OFF = (peaks_ref_allSacc_flash_OFF - peaks_ref_allSacc_baseFlash_OFF)./(peaks_ref_allSacc_flash_OFF + peaks_ref_allSacc_baseFlash_OFF);

%% Peak locations

% relative to flash onset
idx_start_peakExtraction = dataLen_pre;

loc_relToFlash_avg_flash_ON = bsxfun(@minus,loc_avg_flash_ON,flash_onset(1:end-1)')-idx_start_peakExtraction;
loc_relToFlash_avg_baseFlash_ON = bsxfun(@minus,loc_avg_baseFlash_ON,flash_onset(1:end-1)')-idx_start_peakExtraction;

loc_relToFlash_avg_flash_OFF = bsxfun(@minus,loc_avg_flash_OFF,flash_onset(1:end-1)')-idx_start_peakExtraction;
loc_relToFlash_avg_baseFlash_OFF = bsxfun(@minus,loc_avg_baseFlash_OFF,flash_onset(1:end-1)')-idx_start_peakExtraction;

% for template matching modulation index
loc_relToFlash_tm_avg_flash_ON = bsxfun(@minus,loc_tm_avg_flash_ON,(flash_onset(1:end-1)'+dataLen_pre));
loc_relToFlash_tm_avg_flash_OFF = bsxfun(@minus,loc_tm_avg_flash_OFF,(flash_onset(1:end-1)'+dataLen_pre));

latency_relToFlash_tm_avg_baseFlash_ON = loc_tm_avg_baseFlash_ON;
latency_relToFlash_tm_avg_baseFlash_OFF = loc_tm_avg_baseFlash_OFF;


lat_range_OFF = [];
for i = 1:size(loc_relToFlash_tm_avg_flash_OFF,4)
    rgb = loc_relToFlash_tm_avg_flash_OFF(:,:,3,i);
    lat_range_OFF(:,i) = range(rgb);
end
    
    
% -- RELATIVE TO SACCADE PEAK --%
idx_start_sacc = dataLen_pre;

loc_relToSacc_avg_flash_ON = loc_avg_flash_ON - loc_avg_sacc_ON;
loc_relToSacc_avg_flash_OFF = loc_avg_flash_OFF - loc_avg_sacc_OFF;


    % another method
loc_relToSacc_tm_avg_flash_ON = loc_tm_avg_flash_ON - loc_avg_sacc_ON;      % both these variables are from start of the matrix
loc_relToSacc_tm_avg_flash_OFF = loc_tm_avg_flash_OFF - loc_avg_sacc_OFF;


loc_relToSaccOnset_tm_avg_flash_ON = loc_tm_avg_flash_ON - dataLen_pre;
loc_relToSaccOnset_tm_avg_flash_OFF = loc_tm_avg_flash_OFF - dataLen_pre;

  
%% Latency distributions
    saccOffset = 2100-saccOnset+dataLen_pre;       % delays are wrt saccade offset
    saccDur = 100;
    latency_peak_ON = bsxfun(@minus,loc_avg_flash_ON,[saccOffset + flash_onset(1:end-1)]');
    latency_peak_ON(latency_peak_ON==1) = nan;
    timeFromSacc_ON = loc_avg_flash_ON - (saccOffset-saccDur);
    
    latency_peak_OFF = bsxfun(@minus,loc_avg_flash_OFF,[saccOffset + flash_onset(1:end-1)]');
    latency_peak_OFF(latency_peak_ON==1) = nan;
    timeFromSacc_OFF = loc_avg_flash_OFF - (saccOffset-saccDur);
    
 

%% Fig. 9a - Line Plots
    timePoints = [2:5,8:14];
    sc_select = [3];        % sc = [25,50,150,300];
    sc_text = cellstr(num2str(sc'));
    same_cells = 1;
    flashPolarity_select = [1;2];
    rgc_select = 'ALL'; % 'ON' 'OFF' 'ALL'
    mi_select = 'tm'; % {orig,pv,tm} 

       
    plots_idx = 1;
    lim_x = [-200,flash_onset(timePoints(end))+100];
    bins = [-1:0.1:1];
    
    col_on = [0.6,0.6,0.6];
    col_off = [0.2,0.2,0.2];
    col_type = [0.2,0.2,0.2;0.6,0.6,0.6];
    LINE_WIDTH = 2;
    LINE_TYPES = {'-','--',':','-.'};
    
    mi_ON = mi_tm_avg_flash_ON;
    mi_OFF = mi_tm_avg_flash_OFF;
    text_title = 'template matching';
            
    h_main=figure; suptitle(['Modulation index: ',text_title])    
    h_m = [];
    
    s = 1; t = 0; f = 2;


    p_pop_ON = []; %nan(size(mi_avg_flash_ON,1),size(mi_avg_flash_ON,2),size(mi_avg_flash_ON,3),size(mi_avg_flash_ON,4));
    p_pop_OFF = []; %nan(size(mi_avg_flash_OFF,1),size(mi_avg_flash_OFF,2),size(mi_avg_flash_OFF,3),size(mi_avg_flash_OFF,4));

    counter = 0;
    for s = 1:size(sc_select,1)
        for f = 1:size(flashPolarity_select,1)

            points_x_ON = [];     
            for j = 1:size(mi_ON,4)
                points_x_ON(j,:) = squeeze(mi_ON(:,flashPolarity_select(f,1),sc_select(s),j));
            end

            points_x_OFF = [];     
            for j = 1:size(mi_OFF,4)
                points_x_OFF(j,:) = squeeze(mi_OFF(:,flashPolarity_select(f,1),sc_select(s),j));
            end


                idx_toKeep_ON = ~isnan(points_x_ON');
                idx_toKeep_OFF = ~isnan(points_x_OFF');


            if same_cells == 1
                sameCells_ON = ~isnan(mi_ON(:,flashPolarity_select,sc_select,:));   % 1 where the selected masks have responses.
                sameCells_ON = squeeze(all(all(sameCells_ON,2),3));

                sameCells_OFF = ~isnan(mi_OFF(:,flashPolarity_select,sc_select,:));   % 1 where the selected masks have responses.
                sameCells_OFF = squeeze(all(all(sameCells_OFF,2),3));


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
            mean_plot_x_ON = nanmedian(plot_x_ON,1);

            plot_x_OFF  = points_x_OFF;
            plot_x_OFF(~idx_toKeep_OFF') = nan;
            std_plot_x_OFF = nanstd(plot_x_OFF,[],1);
            sem_plot_x_OFF = std_plot_x_OFF./sqrt(nansum(~isnan(plot_x_OFF),1));
            rgb_num_OFF = nansum(~isnan(plot_x_OFF),1);
            rgb_num_OFF = floor(nanmax(rgb_num_OFF));
            mean_plot_x_OFF = nanmedian(plot_x_OFF,1); 

            h_m(1) = subplot(1,1,plots_idx(1));
            hold on
            counter = counter+1;

                switch rgc_select
                    case 'ON'
                        errorbar(flash_onset([timePoints,length(flash_onset)]),[mean_plot_x_ON(timePoints),0],[sem_plot_x_ON(timePoints),0],'Color',col_type(f,:),'LineWidth',LINE_WIDTH,'DisplayName',['ON RGC(',num2str(rgb_num_ON),')|',flash_polarity{flashPolarity_select(f)},'|',sc_text{sc_select(s)},'??m']);

                    case 'OFF'
                        errorbar(flash_onset([timePoints,length(flash_onset)]),[mean_plot_x_OFF(timePoints),0],[sem_plot_x_OFF(timePoints),0],'Color',col_type(f,:),'LineWidth',LINE_WIDTH,'DisplayName',['OFF RGC(',num2str(rgb_num_OFF),')|',flash_polarity{flashPolarity_select(f)},'|',sc_text{sc_select(s)},'??m']);

                    case 'All'
                        if flashPolarity_select(f,1)==2
                            errorbar(flash_onset([timePoints]),[mean_plot_x_OFF(timePoints)],[sem_plot_x_OFF(timePoints)],'Color',col_type(1,:),'LineWidth',LINE_WIDTH,'DisplayName',['OFF RGC(',num2str(rgb_num_OFF),')|',flash_polarity{flashPolarity_select(f)},'|',sc_text{sc_select(s)},'??m']);
                        elseif flashPolarity_select(f,1)==1
                            errorbar(flash_onset([timePoints]),[mean_plot_x_ON(timePoints)],[sem_plot_x_ON(timePoints)],'Color',col_type(2,:),'LineWidth',LINE_WIDTH,'DisplayName',['ON RGC(',num2str(rgb_num_ON),')|',flash_polarity{flashPolarity_select(f)},'|',sc_text{sc_select(s)},'??m']);
                        end
                            

                    otherwise
                        errorbar(flash_onset([timePoints,length(flash_onset)]),[mean_plot_x_ON(timePoints),0],[sem_plot_x_ON(timePoints),0],'Color',col_on,'LineWidth',LINE_WIDTH,'LineStyle',LINE_TYPES{counter},'DisplayName',['ON RGC(',num2str(rgb_num_ON),')|',flash_polarity{flashPolarity_select(f)},'|',sc_text{sc_select(s)},'??m']);
                        errorbar(flash_onset([timePoints,length(flash_onset)]),[mean_plot_x_OFF(timePoints),0],[sem_plot_x_OFF(timePoints),0],'Color',col_off,'LineWidth',LINE_WIDTH,'LineStyle',LINE_TYPES{counter},'DisplayName',['OFF RGC(',num2str(rgb_num_OFF),')|',flash_polarity{flashPolarity_select(f)},'|',sc_text{sc_select(s)},'??m']);
                end
                plot([0,0],[-1,1],'--k')


            legend('-DynamicLegend');

            ylim([-.8 0.3])
            ylabel('Modulation index')
            xlabel('Time from saccade onset (ms)')
            xlim(lim_x)

        end
    end
    plot([-500,2500],[0,0],'--k')
    
    
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

%% Fig. 9c Histograms - flash response latency
    rgc_select = 'ON'; % select 'ON' for ON RGCs plot or 'OFF' for OFF RGCs plot

    timePoints = [2:5,8:14];
    sc_select = [3];        % sc = [25,50,150,300];
    sc_text = cellstr(num2str(sc'));
    same_cells = 1;
    flashPolarity_select = [1;2];
    mi_select = 'tm'; % {orig,pv,tm} 

       
    plots_idx = [1:size(flashPolarity_select,1)]';
    lim_x = [-200,flash_onset(timePoints(end))+100];
    bins = [40:20:600];
    lim_x = [bins(1)-10,bins(end)+10];
    
    col_on = [0.6,0.6,0.6];
    col_off = [0.2,0.2,0.2];
    col_flashPolarity = [0.2,0.2,0.2;0.6,0.6,0.6];
    LINE_WIDTH = 2;
    LINE_TYPES = {'-','--',':','-.'};
    
    mi_ON = mi_tm_avg_flash_ON;
    mi_OFF = mi_tm_avg_flash_OFF;
    text_title = 'template matching';
            
    
    h_main=figure; suptitle(['Flash response latency: ',text_title])    
    h_m = [];
    
    s = 1; t = 0; f = 2;


    p_pop_ON = []; %nan(size(mi_avg_flash_ON,1),size(mi_avg_flash_ON,2),size(mi_avg_flash_ON,3),size(mi_avg_flash_ON,4));
    p_pop_OFF = []; %nan(size(mi_avg_flash_OFF,1),size(mi_avg_flash_OFF,2),size(mi_avg_flash_OFF,3),size(mi_avg_flash_OFF,4));

    counter = 0;
    for s = 1:size(sc_select,1)
        for f = 1:size(flashPolarity_select,1)

            if flashPolarity_select(f)==3
                flashPolarity_select_ON = 1;
                flashPolarity_select_OFF = 2;
            else
                flashPolarity_select_ON = flashPolarity_select(f);
                flashPolarity_select_OFF = flashPolarity_select(f);
            end

            points_x_ON = [];     
            for j = 1:size(mi_ON,4)
                points_x_ON(j,:) = squeeze(mi_ON(:,flashPolarity_select(f,1),sc_select(s),j));
            end

            points_x_OFF = [];     
            for j = 1:size(mi_OFF,4)
                points_x_OFF(j,:) = squeeze(mi_OFF(:,flashPolarity_select(f,1),sc_select(s),j));
            end


                idx_toKeep_ON = ~isnan(points_x_ON');
                idx_toKeep_OFF = ~isnan(points_x_OFF');


            if same_cells == 1
                sameCells_ON = ~isnan(mi_ON(:,flashPolarity_select,sc_select,:));   % 1 where the selected masks have responses.
                sameCells_ON = squeeze(all(all(sameCells_ON,2),3));

                sameCells_OFF = ~isnan(mi_OFF(:,flashPolarity_select,sc_select,:));   % 1 where the selected masks have responses.
                sameCells_OFF = squeeze(all(all(sameCells_OFF,2),3));


            else
                sameCells_ON = idx_toKeep_ON;
                sameCells_OFF = idx_toKeep_OFF;

            end

            idx_toKeep_ON = idx_toKeep_ON & sameCells_ON;
            idx_toKeep_OFF = idx_toKeep_OFF & sameCells_OFF;

            idx_latency_ON = any(idx_toKeep_ON,1);
            idx_latency_OFF = any(idx_toKeep_OFF,1);

            plot_x_ON  = squeeze(latency_relToFlash_tm_avg_baseFlash_ON(end,flashPolarity_select(f,1),sc_select(s),:)); %points_x_ON;
            plot_x_ON(~idx_latency_ON') = nan;
            std_plot_x_ON = nanstd(plot_x_ON,[],1);
            sem_plot_x_ON = std_plot_x_ON./sqrt(nansum(~isnan(plot_x_ON),1));
            rgb_num_ON = nansum(~isnan(plot_x_ON),1);
            rgb_num_ON = nanmax(rgb_num_ON);
            mean_plot_x_ON = nanmean(plot_x_ON,1);

            plot_x_OFF  = squeeze(latency_relToFlash_tm_avg_baseFlash_OFF(end,flashPolarity_select(f,1),sc_select(s),:)); %points_x_OFF;
            plot_x_OFF(~idx_latency_OFF') = nan;
            std_plot_x_OFF = nanstd(plot_x_OFF,[],1);
            sem_plot_x_OFF = std_plot_x_OFF./sqrt(nansum(~isnan(plot_x_OFF),1));
            rgb_num_OFF = nansum(~isnan(plot_x_OFF),1);
            rgb_num_OFF = floor(nanmax(rgb_num_OFF));
            mean_plot_x_OFF = nanmean(plot_x_OFF,1);


            h_m(f) = subplot(1,size(flashPolarity_select,1),plots_idx(f));
            hold on
            counter = counter+1;
                switch rgc_select
                    case 'ON'
                        histogram(plot_x_ON,bins,'EdgeColor',col_on,'FaceColor',col_flashPolarity(flashPolarity_select_ON,:),'DisplayName',['ON RGC | ',flash_polarity{flashPolarity_select_ON},' | N = ',num2str(rgb_num_ON)]) 
                    case 'OFF'
                        histogram(plot_x_OFF,bins,'EdgeColor',col_off,'FaceColor',col_flashPolarity(flashPolarity_select_ON,:),'DisplayName',['OFF RGC | ',flash_polarity{flashPolarity_select_OFF},' | N = ',num2str(rgb_num_OFF)]) 
                end


            legend('-DynamicLegend');

            ylabel('Num of cells')
            xlabel('Flash response latency (ms)')
            xlim(lim_x)

        end
    end
    plot([-500,2500],[0,0],'--k')
    
    
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

%% Fig. 9d - bottom row is top row in paper figure
    rgc_select = 'ON'; % select here 'ON' or 'OFF'


    timePoints =[2:5]% [2:5]%,8:14];
    sc_select = [3];        % sc = [25,50,150,300];
    sc_text = cellstr(num2str(sc'));
    same_cells = 1;
    flashPolarity_select = [1;2]; % [1=OFF; 2=ON; 3=anti=preferred]
    select_metric = 'rel_sacc_onset'; % ['rel_flash','rel_sacc_peak','rel_sacc_onset']
    mi_select = 'tm'; % {'orig','pv','tm'}
    transpose_figure = 0;
    
    idx_cell_mark = 0%162;

 
    plots_idx = [1:(length(timePoints)*size(flashPolarity_select,1))]';
    
    if transpose_figure
        plots_idx = reshape(plots_idx',size(flashPolarity_select,1),length(timePoints));        % [time,mask_conds,sc]
    else
        plots_idx = reshape(plots_idx',length(timePoints),size(flashPolarity_select,1));        % [time,mask_conds,sc]
    end
    
    lim_y = [-1.1,0.5];
    bins = [-1:0.1:1];
    baseVal = -1.2;
    
    col_on = [0.6,0.6,0.6];
    col_off = [0.2,0.2,0.2];
    col_flashPolarity = [0.2,0.2,0.2;0.6,0.6,0.6];
    LINE_WIDTH = 2;
    LINE_TYPES = {'-','--',':','-.'};
    
    
    h_main=figure; suptitle([rgc_select,' RGC'])    
    h_m = [];
    
    s = 1; t = 0; f = 2;
    
            

    mi_ON = mi_tm_avg_flash_ON;
    mi_OFF = mi_tm_avg_flash_OFF;
    text_title = 'templation matching';
            

            
    switch select_metric
        case 'rel_flash'
            if strcmp('orig',mi_select)
                peak_loc_ON = loc_relToFlash_avg_flash_ON;
                peak_loc_OFF = loc_relToFlash_avg_flash_OFF;
            else
                peak_loc_ON = loc_relToFlash_tm_avg_flash_ON;
                peak_loc_OFF = loc_relToFlash_tm_avg_flash_OFF;
            end

            lim_x = [-200,600];
            x_line = [0,0];
            x_label = 'Time of flash peak from flash onset';

            
        case 'rel_sacc_peak'
            if strcmp('orig',mi_select)
                peak_loc_ON = loc_relToSacc_avg_flash_ON;
                peak_loc_OFF = loc_relToSacc_avg_flash_OFF;
            else
                peak_loc_ON = loc_relToSacc_tm_avg_flash_ON;
                peak_loc_OFF = loc_relToSacc_tm_avg_flash_OFF;
            end
            lim_x = [-700,1600];
            x_line = [1,1];
            x_label = 'Time of flash peak from saccade peak';
            
        case 'rel_sacc_onset'
            if strcmp('orig',mi_select)
                peak_loc_ON = loc_relToSaccOnset_avg_flash_ON;
                peak_loc_OFF = loc_relToSaccOnset_avg_flash_OFF;
            else
                peak_loc_ON = loc_relToSaccOnset_tm_avg_flash_ON;
                peak_loc_OFF = loc_relToSaccOnset_tm_avg_flash_OFF;
            end
            lim_x = [-700,1600];
            x_line = [1,1];
            x_label = 'Time of flash peak from saccade onset';


    end
            
    counter = 0;
    for s = 1:size(sc_select,1)
        for f = 1:size(flashPolarity_select,1)
            for t=1:length(timePoints)
            
                if flashPolarity_select(f)==3
                    flashPolarity_select_ON = 1;
                    flashPolarity_select_OFF = 2;
                else
                    flashPolarity_select_ON = flashPolarity_select(f);
                    flashPolarity_select_OFF = flashPolarity_select(f);
                end

                points_x_ON = [];     
                points_y_ON = [];
                for j = 1:size(mi_ON,4)
                    points_x_ON(j,:) = squeeze(peak_loc_ON(timePoints(t),flashPolarity_select_ON,sc_select(s),j));
                    points_y_ON(j,:) = squeeze(mi_ON(timePoints(t),flashPolarity_select_ON,sc_select(s),j));
                end

                points_x_OFF = [];     
                points_y_OFF = []; 
                for j = 1:size(mi_OFF,4)
                    points_x_OFF(j,:) = squeeze(peak_loc_OFF(timePoints(t),flashPolarity_select_OFF,sc_select(s),j));
                    points_y_OFF(j,:) = squeeze(mi_OFF(timePoints(t),flashPolarity_select_OFF,sc_select(s),j));
                end


                idx_toKeep_ON = ~isnan(points_y_ON');
                idx_toKeep_OFF = ~isnan(points_y_OFF');


                if same_cells == 1
                    sameCells_ON = ~isnan(mi_ON(:,flashPolarity_select,sc_select,:));   % 1 where the selected masks have responses.
                    sameCells_ON = squeeze(all(all(sameCells_ON,2),3));

                    sameCells_OFF = ~isnan(mi_OFF(:,flashPolarity_select,sc_select,:));   % 1 where the selected masks have responses.
                    sameCells_OFF = squeeze(all(all(sameCells_OFF,2),3));

                else
                    sameCells_ON = idx_toKeep_ON;
                    sameCells_OFF = idx_toKeep_OFF;

                end

                idx_toKeep_ON = idx_toKeep_ON & sameCells_ON(timePoints(t),:);
                idx_toKeep_OFF = idx_toKeep_OFF & sameCells_OFF(timePoints(t),:);

                plot_x_ON  = points_x_ON;
                plot_x_ON(~idx_toKeep_ON') = nan;
                plot_y_ON = points_y_ON;
                plot_y_ON(~idx_toKeep_ON') = nan;
                rgb_num_ON = nansum(~isnan(plot_x_ON),1);
                rgb_num_ON = nanmax(rgb_num_ON);
                mean_plot_x_ON = nanmedian(plot_x_ON,1);

                plot_x_OFF  = points_x_OFF;
                plot_x_OFF(~idx_toKeep_OFF') = nan;
                plot_y_OFF  = points_y_OFF;
                plot_y_OFF(~idx_toKeep_OFF') = nan;
                rgb_num_OFF = nansum(~isnan(plot_x_OFF),1);
                rgb_num_OFF = nanmax(rgb_num_OFF);
                mean_plot_x_OFF = nanmedian(plot_x_OFF,1);
                
                if transpose_figure == 1
                    h_m(t,f) = subplot(length(timePoints),size(flashPolarity_select,1),plots_idx(f,t));
                else
                    h_m(t,f) = subplot(size(flashPolarity_select,1),length(timePoints),plots_idx(t,f));
                end
                hold on

                switch rgc_select
                    case 'ON'
                        scatter(plot_x_ON,plot_y_ON,'MarkerEdgeColor',col_on,'MarkerFaceColor',col_flashPolarity(flashPolarity_select_ON,:),'DisplayName',[flash_polarity{flashPolarity_select_ON},' | N = ',num2str(rgb_num_ON)]) 
                    case 'OFF'
                        scatter(plot_x_OFF,plot_y_OFF,'MarkerEdgeColor',col_off,'MarkerFaceColor',col_flashPolarity(flashPolarity_select_ON,:),'DisplayName',[flash_polarity{flashPolarity_select_OFF},' | N = ',num2str(rgb_num_OFF)]) 

                    otherwise
                        scatter(plot_x_ON,plot_y_ON,'MarkerEdgeColor',col_on,'MarkerFaceColor',col_on) 
                        scatter(plot_x_OFF,plot_y_OFF,'MarkerEdgeColor',col_off,'MarkerFaceColor',col_off)
                end
                
                if idx_cell_mark>0
                    scatter(plot_x_OFF(idx_cell_mark),plot_y_OFF(idx_cell_mark),'r')
                end

                plot(x_line,[-1,1],'--','color',[0.3,0.3,0.3])
                plot([-200,2500],[0,0],'--k')
                xlim(lim_x)
                ylim(lim_y)
                axis square
                
                legend('-DynamicLegend','location','best')
            end

        end
    end
    
    
    
    ylabel(h_m(1,1),'Modulation Index')
    xlabel(h_m(2,1),x_label)
    
    ylabel(h_m(1,2),'Modulation Index')
    xlabel(h_m(2,2),x_label)

    
    for i = 1:size(timePoints,2)
        title(h_m(i,1),[num2str(flash_onset(timePoints(i))),'ms'])
    end

    
