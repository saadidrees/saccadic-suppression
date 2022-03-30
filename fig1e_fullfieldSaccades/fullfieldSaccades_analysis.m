% Fig 1e and Suppl Fig 3

% Loads full-field data from all full-field saccade experiments, globalcomponent experiments and the local component experiments.
% The data in these files include modualtion index for each cell at each flash and spatial scale. Also contains the corresponding peak firing rates

%%

clear; clc

matFiles = dir(['*_standardParadigm_*.mat']);   % datafiles containing the peak firing rates for different conditions

flash_onset = [-265,-233,-217,-184,-167,-150,-133,-117,-50,17,33,50,100,250,500,1000];      % times at which flashes were presented wrt to saccade offset [add 100 to get wrt to saccade onset]
baseFlash_onset = 2000;     % time of control flash or the baseline flash
sc = [25,50,150,300];       % spatial scale vector


% initialize variables to store data in from all the data files
    % All the variables below will have 3 dimensions as follows:
    % [flash_onset,sc,cells]

mi_avg_flash_ON = []; 
num_onCells = zeros(length(sc),length(flash_onset));
p_avg_flash_ON = [];
pow_p_avg_flash_ON = [];
peaks_ref_avg_flash_ON = [];
peaks_ref_avg_baseFlash_ON = [];
spikeC_flash_ON = []; 
spikeC_baseFlash_ON = [];
latency_peak_ON = [];
timeFromSacc_ON = []; 
assoc_R_ON = [];
assoc_p_ON = [];

mi_avg_flash_OFF = [];
num_offCells = zeros(length(sc),length(flash_onset));
p_avg_flash_OFF = [];
pow_p_avg_flash_OFF = [];
peaks_ref_avg_flash_OFF = [];
peaks_ref_avg_baseFlash_OFF = [];
spikeC_flash_OFF = [];
spikeC_baseFlash_OFF = []; 
latency_peak_OFF = [];
timeFromSacc_OFF = [];
assoc_R_OFF = [];
assoc_p_OFF = [];

unameON_all = {};
unameOFF_all = {};
experiments_name = {};


for m = 1:size(matFiles,1)
% Load Experiment Files
    sprintf('loading %s\n',matFiles(m).name)
    rgb = load(fullfile(pwd,matFiles(m).name));     % most of the variables used from this data file have 3 dimensions [flash_onset,sc,cells]
    
    [~,sc_idx] = ismember(rgb.sc,sc);
    if rgb.flash_onset(end) > 1800 
        rgb.flash_onset = rgb.flash_onset(1:end-1);
    end
    [~,delay_idx] = ismember(rgb.flash_onset,flash_onset);
    
    % ON Cells
        temp = nan(length(flash_onset),length(sc),size(rgb.mi_avg_flash_ON,3));
        temp(delay_idx,sc_idx,:) = rgb.mi_avg_flash_ON;
        mi_avg_flash_ON = cat(3,mi_avg_flash_ON,temp);

        temp = nan(length(sc),length(flash_onset));
        temp(sc_idx,delay_idx) = rgb.num_onCells;
        num_onCells = nansum(cat(3,num_onCells,temp),3);

        temp = nan(length(flash_onset),length(sc),size(rgb.p_avg_flash_ON,3));
        temp(delay_idx,sc_idx,:) = rgb.p_avg_flash_ON;
        p_avg_flash_ON = cat(3,p_avg_flash_ON,temp);

        temp = nan(length(flash_onset),length(sc),size(rgb.pow_p_avg_flash_ON,3));
        temp(delay_idx,sc_idx,:) = rgb.pow_p_avg_flash_ON;
        pow_p_avg_flash_ON = cat(3,pow_p_avg_flash_ON,temp);

        temp = nan(length(flash_onset),length(sc),size(rgb.peaks_refined_avg_flash_ON,3));
        temp(delay_idx,sc_idx,:) = rgb.peaks_refined_avg_flash_ON;
        peaks_ref_avg_flash_ON = cat(3,peaks_ref_avg_flash_ON,temp);

        temp = nan(length(flash_onset),length(sc),size(rgb.peaks_refined_avg_baseFlash_ON,3));
        temp(delay_idx,sc_idx,:) = rgb.peaks_refined_avg_baseFlash_ON;
        peaks_ref_avg_baseFlash_ON = cat(3,peaks_ref_avg_baseFlash_ON,temp);

        temp = nan(length(flash_onset),length(sc),size(rgb.spikeC_flash_ON,3));
        temp(delay_idx,sc_idx,:) = rgb.spikeC_flash_ON;
        spikeC_flash_ON = cat(3,spikeC_flash_ON,temp);

        temp = nan(length(flash_onset),length(sc),size(rgb.spikeC_baseFlash_ON,3));
        temp(delay_idx,sc_idx,:) = rgb.spikeC_baseFlash_ON;
        spikeC_baseFlash_ON = cat(3,spikeC_baseFlash_ON,temp);

        temp = nan(length(flash_onset),length(sc),size(rgb.latency_peak_ON,3));
        temp(delay_idx,sc_idx,:) = rgb.latency_peak_ON;
        latency_peak_ON = cat(3,latency_peak_ON,temp);

        temp = nan(length(flash_onset),length(sc),size(rgb.timeFromSacc_ON,3));
        temp(delay_idx,sc_idx,:) = rgb.timeFromSacc_ON;
        timeFromSacc_ON = cat(3,timeFromSacc_ON,temp);
        
        temp = nan(length(flash_onset),length(sc),size(rgb.assoc_R_ON,3));
        temp(delay_idx,sc_idx,:) = rgb.assoc_R_ON;
        assoc_R_ON = cat(3,assoc_R_ON,temp);
        
        temp = nan(length(flash_onset),length(sc),size(rgb.assoc_p_ON,3));
        temp(delay_idx,sc_idx,:) = rgb.assoc_p_ON;
        assoc_p_ON = cat(3,assoc_p_ON,temp);

        
    % OFF Cells
        temp = nan(length(flash_onset),length(sc),size(rgb.mi_avg_flash_OFF,3));
        temp(delay_idx,sc_idx,:) = rgb.mi_avg_flash_OFF;
        mi_avg_flash_OFF = cat(3,mi_avg_flash_OFF,temp);

        temp = nan(length(sc),length(flash_onset));
        temp(sc_idx,delay_idx) = rgb.num_offCells;
        num_offCells = nansum(cat(3,num_offCells,temp),3);

        temp = nan(length(flash_onset),length(sc),size(rgb.p_avg_flash_OFF,3));
        temp(delay_idx,sc_idx,:) = rgb.p_avg_flash_OFF;
        p_avg_flash_OFF = cat(3,p_avg_flash_OFF,temp);

        temp = nan(length(flash_onset),length(sc),size(rgb.pow_p_avg_flash_OFF,3));
        temp(delay_idx,sc_idx,:) = rgb.pow_p_avg_flash_OFF;
        pow_p_avg_flash_OFF = cat(3,pow_p_avg_flash_OFF,temp);

        temp = nan(length(flash_onset),length(sc),size(rgb.peaks_refined_avg_flash_OFF,3));
        temp(delay_idx,sc_idx,:) = rgb.peaks_refined_avg_flash_OFF;
        peaks_ref_avg_flash_OFF = cat(3,peaks_ref_avg_flash_OFF,temp);

        temp = nan(length(flash_onset),length(sc),size(rgb.peaks_refined_avg_baseFlash_OFF,3));
        temp(delay_idx,sc_idx,:) = rgb.peaks_refined_avg_baseFlash_OFF;
        peaks_ref_avg_baseFlash_OFF = cat(3,peaks_ref_avg_baseFlash_OFF,temp);

        temp = nan(length(flash_onset),length(sc),size(rgb.spikeC_flash_OFF,3));
        temp(delay_idx,sc_idx,:) = rgb.spikeC_flash_OFF;
        spikeC_flash_OFF = cat(3,spikeC_flash_OFF,temp);

        temp = nan(length(flash_onset),length(sc),size(rgb.spikeC_baseFlash_OFF,3));
        temp(delay_idx,sc_idx,:) = rgb.spikeC_baseFlash_OFF;
        spikeC_baseFlash_OFF = cat(3,spikeC_baseFlash_OFF,temp);

        temp = nan(length(flash_onset),length(sc),size(rgb.latency_peak_OFF,3));
        temp(delay_idx,sc_idx,:) = rgb.latency_peak_OFF;
        latency_peak_OFF = cat(3,latency_peak_OFF,temp);

        temp = nan(length(flash_onset),length(sc),size(rgb.timeFromSacc_OFF,3));
        temp(delay_idx,sc_idx,:) = rgb.timeFromSacc_OFF;
        timeFromSacc_OFF = cat(3,timeFromSacc_OFF,temp);
        
        temp = nan(length(flash_onset),length(sc),size(rgb.assoc_R_OFF,3));
        temp(delay_idx,sc_idx,:) = rgb.assoc_R_OFF;
        assoc_R_OFF = cat(3,assoc_R_OFF,temp);
        
        temp = nan(length(flash_onset),length(sc),size(rgb.assoc_p_OFF,3));
        temp(delay_idx,sc_idx,:) = rgb.assoc_p_OFF;
        assoc_p_OFF = cat(3,assoc_p_OFF,temp);

        
        unameON_all = cat(1,unameON_all,rgb.unameON_all);
        unameOFF_all = cat(1,unameOFF_all,rgb.unameOFF_all);
        experiments_name = cat(1,experiments_name,rgb.experiments);
end
flash_onset = flash_onset+100;  % now flash onset is wrt saccade onset

% delays to remove because of low N numbers
idx_delaysToRemove = [7,8,11];    % 33

flash_onset(idx_delaysToRemove) = [];

mi_avg_flash_ON(idx_delaysToRemove,:,:) = [];
p_avg_flash_ON(idx_delaysToRemove,:,:) = [];
pow_p_avg_flash_ON(idx_delaysToRemove,:,:) = [];
peaks_ref_avg_flash_ON(idx_delaysToRemove,:,:) = [];
peaks_ref_avg_baseFlash_ON(idx_delaysToRemove,:,:) = [];
spikeC_flash_ON(idx_delaysToRemove,:,:) = [];
spikeC_baseFlash_ON(idx_delaysToRemove,:,:) = []; 
latency_peak_ON(idx_delaysToRemove,:,:) = [];
timeFromSacc_ON(idx_delaysToRemove,:,:) = [];
assoc_R_ON(idx_delaysToRemove,:,:) = [];
assoc_p_ON(idx_delaysToRemove,:,:) = [];
num_onCells(:,idx_delaysToRemove) = [];

mi_avg_flash_OFF(idx_delaysToRemove,:,:) = [];
p_avg_flash_OFF(idx_delaysToRemove,:,:) = [];
pow_p_avg_flash_OFF(idx_delaysToRemove,:,:) = [];
peaks_ref_avg_flash_OFF(idx_delaysToRemove,:,:) = [];
peaks_ref_avg_baseFlash_OFF(idx_delaysToRemove,:,:) = [];
spikeC_flash_OFF(idx_delaysToRemove,:,:) = [];
spikeC_baseFlash_OFF(idx_delaysToRemove,:,:) = []; 
latency_peak_OFF(idx_delaysToRemove,:,:) = [];
timeFromSacc_OFF(idx_delaysToRemove,:,:) = [];
assoc_R_OFF(idx_delaysToRemove,:,:) = [];
assoc_p_OFF(idx_delaysToRemove,:,:) = [];
num_offCells(:,idx_delaysToRemove) = [];

% Mean across all cells
mean_mi_avg_flash_ON = nanmean(mi_avg_flash_ON,3);
std_mi_avg_flash_ON = nanstd(mi_avg_flash_ON,[],3);
sem_mi_avg_flash_ON = std_mi_avg_flash_ON./sqrt(num_onCells');
sig_ON = p_avg_flash_ON<0.05 & pow_p_avg_flash_ON > 80;

mean_mi_avg_flash_OFF = nanmean(mi_avg_flash_OFF,3);
std_mi_avg_flash_OFF = nanstd(mi_avg_flash_OFF,[],3);
sem_mi_avg_flash_OFF = std_mi_avg_flash_OFF./sqrt(num_offCells');
sig_OFF = p_avg_flash_OFF < 0.05 & pow_p_avg_flash_OFF > 80;


% Pool ON and OFF RGCs
mi_avg_flash_combined = cat(3,mi_avg_flash_ON,mi_avg_flash_OFF);
numCells_combined = num_onCells + num_offCells;
mean_mi_avg_flash_combined = nanmean(mi_avg_flash_combined,3);
std_mi_avg_flash_combined = nanstd(mi_avg_flash_combined,[],3);
sem_mi_avg_flash_combined = std_mi_avg_flash_combined./sqrt(numCells_combined');
median_mi_avg_flash_combined = nanmedian(mi_avg_flash_combined,3);
peaks_flash_combined = cat(3,peaks_ref_avg_flash_ON,peaks_ref_avg_flash_OFF);
peaks_base_combined = cat(3,peaks_ref_avg_baseFlash_ON,peaks_ref_avg_baseFlash_OFF);
spikeC_flash_combined = cat(3,spikeC_flash_ON,spikeC_flash_OFF);
spikeC_baseFlash_combined = cat(3,spikeC_baseFlash_ON,spikeC_baseFlash_OFF);
percSpikeC_combined = 100*(spikeC_flash_combined./spikeC_baseFlash_combined);
latency_peak_combined = cat(3,latency_peak_ON,latency_peak_OFF);
timeFromSacc_combined = cat(3,timeFromSacc_ON,timeFromSacc_OFF);
sig_combined = cat(3,sig_ON,sig_OFF);
uname_all = cat(1,unameON_all,unameOFF_all);
ci_mi_combined = cat(3,mean_mi_avg_flash_combined-(1.95*sem_mi_avg_flash_combined),mean_mi_avg_flash_combined+(1.95*sem_mi_avg_flash_combined));

p_avg_flash_combined = cat(3,p_avg_flash_ON,p_avg_flash_OFF);
pow_p_avg_flash_combined = cat(3,pow_p_avg_flash_ON,pow_p_avg_flash_OFF);

assoc_R_combined = cat(3,assoc_R_ON,assoc_R_OFF);
assoc_p_combined = cat(3,assoc_p_ON,assoc_p_OFF);


timePoints_all = 1:length(flash_onset);
sig_mi_ON = [];
sig_mi_OFF = [];
nonsig_mi_ON = [];
nonsig_mi_OFF = [];

for i = 1:length(timePoints_all)
    for j = 1:length(sc)
        rgb = squeeze(mi_avg_flash_ON(i,j,:));
        sig_mi_ON(i,j) = nanmean(rgb(sig_ON(i,j,:)));
        nonsig_mi_ON(i,j) = nanmean(rgb(~sig_ON(i,j,:)));
                        
        rgb = squeeze(mi_avg_flash_OFF(i,j,:));
        sig_mi_OFF(i,j) = nanmean(rgb(sig_OFF(i,j,:)));
        nonsig_mi_OFF(i,j) = nanmean(rgb(~sig_OFF(i,j,:)));
    end
end


% significance testing - signrank test - combined

p_sig_sc = [];
select_sc = [1,3];  % select the index of the spatial scales to compare. For eg, [1,3] will test whether the population modulation index differs across the spatial scales, for each time point.

for i = 1:length(flash_onset)
    p_sig_sc(i) = signrank(squeeze(mi_avg_flash_combined(i,select_sc(1),:)),squeeze(mi_avg_flash_combined(i,select_sc(2),:)));
end
p_sig_sc



% significance testing across ON and OFF RGCs
p_sig_ON_OFF = nan(length(sc),length(flash_onset));
for j = 1:length(sc)
    for i = 1:length(flash_onset)
        try
            p_sig_ON_OFF(j,i) = ranksum(squeeze(mi_avg_flash_ON(i,j,:)),squeeze(mi_avg_flash_OFF(i,j,:)));
        end
    end
end
p_sig_ON_OFF = [flash_onset;p_sig_ON_OFF];

% Significance testing for suppression - ON and OFF

p_sig_supp_ON = [];
p_sig_supp_OFF = [];

p_sig_supp_combined = [];
h = [];
effect_size = [];

for i = 1:size(mi_avg_flash_ON,1)
    for j = 1:size(mi_avg_flash_ON,2)
        rgb = squeeze(mi_avg_flash_ON(i,j,:));
        rgb(isnan(rgb)) = [];
        if ~isempty(rgb)
           [p_sig_supp_ON(j,i),~,stats] = signrank(rgb);
        else
            p_sig_supp_ON(j,i) = nan;
        end

        rgb = squeeze(mi_avg_flash_OFF(i,j,:));
        rgb(isnan(rgb)) = [];
        if ~isempty(rgb)
           [p_sig_supp_OFF(j,i),~,stats] = signrank(rgb);
        else
            p_sig_supp_OFF(j,i) = nan;
        end


        rgb = squeeze(mi_avg_flash_combined(i,j,:));
        rgb(isnan(rgb)) = [];
        if ~isempty(rgb)
            [p_sig_supp_combined(j,i)] = signrank(rgb,0);
%             [~,p_sig_supp_combined(j,i)] = ttest(rgb,0);
        else
            p_sig_supp_combined(j,i) = nan;
        end
    end
end


%% Fig 1e and Suppl Fig. 3
% LinePlots - Figure 1e
timePoints_pre = 1:5;   % select the index for pre-sacc flashes to plot
timePoints_post = [6:13]; % select the index for post-sacc flashes to plot flash_onset = [-265,-233,-217,-184,-167,-150,-50,17,50,100,250,500,1000] % from saccade onset
sc_select = [1,2,3,4];  % select the index of spatial scales to plot
sc_text = {'25µm','50µm','150µm','300µm'};

figure;

for s = sc_select
    subplot(2,1,1)
    hold on
    
    h_l = errorbar(repmat(flash_onset(timePoints_post)',1,length(s)),mean_mi_avg_flash_ON(timePoints_post,s),sem_mi_avg_flash_ON(timePoints_post,s),'DisplayName',sc_text{s});
%     legend(sc_text(s))
    errorbar(repmat(flash_onset(timePoints_pre)',1,length(s)),mean_mi_avg_flash_ON(timePoints_pre,s),sem_mi_avg_flash_ON(timePoints_pre,s),'LineStyle',':','Color',h_l.Color,'DisplayName',sc_text{s})
    ylim([-0.6 0.3])
    title(['ON Cells - N = ',num2str(size(mi_avg_flash_ON,3))])
    ylabel('Modulation index')
    xlabel('Time from saccade onset (ms)')
    legend('-DynamicLegend')
    subplot(2,1,2)
    hold on
    h_l = errorbar(repmat(flash_onset(timePoints_post)',1,length(s)),mean_mi_avg_flash_OFF(timePoints_post,s),sem_mi_avg_flash_OFF(timePoints_post,s),'DisplayName',sc_text{s});
    legend('-DynamicLegend')
    errorbar(repmat(flash_onset(timePoints_pre)',1,length(s)),mean_mi_avg_flash_OFF(timePoints_pre,s),sem_mi_avg_flash_OFF(timePoints_pre,s),'LineStyle',':','Color',h_l.Color,'DisplayName',sc_text{s})
    ylim([-0.6 0.3])
    title(['OFF Cells - N = ',num2str(size(mi_avg_flash_OFF,3))])
    ylabel('Modulation index')
    xlabel('Time from saccade onset (ms)')
    
end

% Histograms - Supplemental Figure 3

    timePoints = [1:5,6:13];  % flash_onset = [-265,-233,-217,-184,-167,-150,-133,-117,-50,17,33,50,100,250,500,1000];
    sc_select = [1;2;3;4];        % sc = [25,50,300];
    sc_text = cellstr(num2str(sc'));

   
    plots_idx = [1:(length(timePoints)*size(sc_select,1))]';
    plots_idx = reshape(plots_idx',length(timePoints),size(sc_select,1));        % [time,mask_conds,sc]
    
    lim_x = [-1.2,1.2];
    bins = [-1:0.1:1];
    
    col_on = [0.6,0.6,0.6];
    col_off = [0.2,0.2,0.2];
    
    
    h_main=figure; suptitle(['ON & OFF Cells - standard paradigm'])    
    h_m = [];
    
    s = 1; m = 1; d = 1; t = 1;


    
    for s = 1:size(sc_select,1)
                for t = 1:length(timePoints)
                    
                    points_x_ON = [];     
                    for j = 1:size(mi_avg_flash_ON,3)
                        points_x_ON(j,:) = squeeze(mi_avg_flash_ON(timePoints(t),sc_select(s),j));
                    end

                    points_x_OFF = [];     
                    for j = 1:size(mi_avg_flash_OFF,3)
                        points_x_OFF(j,:) = squeeze(mi_avg_flash_OFF(timePoints(t),sc_select(s),j));
                    end


                    h_m(t,s) = subplot(size(sc_select,1),length(timePoints),plots_idx(t,s));
                    hold on


                    plot_x_ON  = points_x_ON(~isnan(points_x_ON));
                    plot_x_OFF  = points_x_OFF(~isnan(points_x_OFF));
                    

                    [f_off_all,x] = hist(plot_x_OFF,bins);
                    h = bar(x,(f_off_all/trapz(f_off_all)));
                    set(h,'FaceColor',col_off,'EdgeColor','k');

                    [f_on_all,x] = hist(plot_x_ON,bins);
                    h = bar(x,(f_on_all/trapz(f_on_all)));
                    set(h,'FaceColor',col_on,'EdgeColor','k');
                    set(h,'FaceAlpha',0.6)


                    plot([0,0],[0,.4],'--k')

                    xlim(lim_x)
                    ylim([0,.5])

                    axis square
                    legend({['OFF ',num2str(length(plot_x_OFF))],['ON: ',num2str(length(plot_x_ON))]})

                end

              ylabel(['sc: ',sc_text{sc_select(s)},' Âµm'])

    end
    for t = 1:length(timePoints)
        title(h_m(t,1),[num2str(flash_onset(timePoints(t))),'ms'])
    end
