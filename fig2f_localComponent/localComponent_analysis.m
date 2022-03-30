
%% Fig 2f and Suppl Fig. 7c
% loads data from all checkerboard mask experiments. The main data used is the peak response of flashes at different delays after saccade
% and for different conditions. in general, data is organized as [delays,maskRegions,maskSize,sc,drugs,cells].
% Visual stimulation regions are three in total. Set of all squares assigned label '3'. Then one set of alternate checkers '1' and the remaining checkers as region '2'
% Depending on the experiment, data could be collected at 5 different settings:
% 1. sacc+flash in all regions - regSacc_3-regFlash_3
% 2. sacc in region 1 flash in region 2 - regSacc_1-regFlash_2
% 3. sacc in 2 and flash in 1 - regSacc_2-regFlash_1
% 4. sacc in 1 and flash in 1 - regSacc_1-regFlash_1  % this condition is not used
% 5. sacc in 2 and flash in 2 - regSacc_2-regFlash_2  % this condition is not used

%% Load data 
clear all; clc; set(0,'defaultTextInterpreter','none');

load('data_checkerboardMask.mat')

nds = [5];
sc = [25,150];
maskSize = [100,150];  % size of checkerboard in µm
conds_mask = {'regSacc_3-regFlash_3','regSacc_1-regFlash_2','regSacc_2-regFlash_1','regSacc_1-regFlash_1','regSacc_2-regFlash_2'};
regions = {'flash+sacc in all','no sacc in center','sacc in center'};
flash_onset = [17,50,100,250,500,2000];
drugs = {'none','gabazine+ptx'};
numSaccs = 40;

dataLen_pre = 400;
dataLen_post = 3600;
saccOnset = 2000;
idx_data = saccOnset-dataLen_pre+1:saccOnset+dataLen_post;

baseFlash_idx = length(flash_onset);      % take the 9th flash i.e. at 1867ms as the baseline flash
medThresh = [0.6];

electrode_coordinates = load('MEA60_elec_positions.mat','electrode_coordinates');
electrode_coordinates = (electrode_coordinates.electrode_coordinates)/3.75;

%% Calculate modulation indices

% Modulation Index - refined peaks - avg
    mi_avg_flash_ON = (peaks_ref_avg_flash_ON - peaks_ref_avg_baseFlash_ON)./(peaks_ref_avg_flash_ON + peaks_ref_avg_baseFlash_ON);     % [delays,maskRegions,maskSize,sc,drugs,cells]
    mi_avg_flash_OFF = (peaks_ref_avg_flash_OFF - peaks_ref_avg_baseFlash_OFF)./(peaks_ref_avg_flash_OFF + peaks_ref_avg_baseFlash_OFF);
    
% Modulation Index - refined peaks - all saccs
    mi_allSacc_flash_ON = (peaks_ref_allSacc_flash_ON - peaks_ref_allSacc_baseFlash_ON)./(peaks_ref_allSacc_flash_ON + peaks_ref_allSacc_baseFlash_ON);
    mi_allSacc_flash_OFF = (peaks_ref_allSacc_flash_OFF - peaks_ref_allSacc_baseFlash_OFF)./(peaks_ref_allSacc_flash_OFF + peaks_ref_allSacc_baseFlash_OFF);

    
% Normalize saccade areas
    areaNorm_avg_sacc_ON = bsxfun(@rdivide,area_avg_sacc_ON,area_avg_sacc_ON(:,1,:,:,:,:));     % [delays,maskRegions,maskSize,,sc,drugs,cells]
    areaNorm_avg_sacc_OFF = bsxfun(@rdivide,area_avg_sacc_OFF,area_avg_sacc_OFF(:,1,:,:,:,:));
    
%% Area stuff - to see how much of cell's rf is in a specific stimulation region
rf_dia_sig = 1;
rf_dia_ON = round(2*rf_dia_sig*nanmax(rf_fit_sigma_ON,[],1).*pix2um_ON);   % rf dia in µm
rf_dia_OFF = round(2*rf_dia_sig*nanmax(rf_fit_sigma_OFF,[],1).*pix2um_OFF);   % rf dia in µm


regIdx_sacc_ON = 3*ones(size(regIdx_flash_ON));regIdx_sacc_ON(regIdx_flash_ON==3) = 2;
regIdx_sacc_OFF = 3*ones(size(regIdx_flash_OFF));regIdx_sacc_OFF(regIdx_flash_OFF==3) = 2;

regArea_sacc_ON = []; regArea_flash_ON = []; regArea_total_ON = [];
regArea_sacc_OFF = []; regArea_flash_OFF = []; regArea_total_OFF = [];

for i = 1:length(uname_ON)
    regArea_sacc_ON(:,i) = nansum(nansum(rf_plot_ON(:,:,i).*double(mask_regions_ON{i}{regIdx_flash_ON(i)})/255))*100;       % white regions in the mask is where saccades are being shown..!Black region is always for flash
    regArea_flash_ON(:,i) = nansum(nansum(rf_plot_ON(:,:,i).*double(mask_regions_ON{i}{regIdx_sacc_ON(i)})/255))*100;       % white regions in the mask is where saccades are being shown..!Black region is always for flash
    regArea_total_ON(:,i) = nansum(nansum(rf_plot_ON(:,:,i).*double(mask_regions_ON{i}{1})/255))*100;
end
for i = 1:length(uname_OFF)
    regArea_sacc_OFF(:,i) = nansum(nansum(rf_plot_OFF(:,:,i).*double(mask_regions_OFF{i}{regIdx_flash_OFF(i)})/255))*100;       % white regions in the mask is where saccades are being shown..!Black region is always for flash
    regArea_flash_OFF(:,i) = nansum(nansum(rf_plot_OFF(:,:,i).*double(mask_regions_OFF{i}{regIdx_sacc_OFF(i)})/255))*100;       % white regions in the mask is where saccades are being shown..!Black region is always for flash
    regArea_total_OFF(:,i) = nansum(nansum(rf_plot_OFF(:,:,i).*double(mask_regions_OFF{i}{1})/255))*100;
end

%% RF location relative to mask

dist_rfFromSacc_ON = [];
parfor j = 1:length(uname_ON)   
    rgb_dist = bsxfun(@minus,electrode_coordinates,rf_fit_center_ON(:,j)');
    temp = sqrt(sum(rgb_dist.^2,2));
    [min_dist,min_dist_idx] = min(temp);
    dist_rfFromSacc_ON(j,:) = rgb_dist(min_dist_idx,:);
end

dist_rfFromSacc_OFF = [];
parfor j = 1:length(uname_OFF)   
    rgb_dist = bsxfun(@minus,electrode_coordinates,rf_fit_center_OFF(:,j)');
    temp = sqrt(sum(rgb_dist.^2,2));
    [min_dist,min_dist_idx] = min(temp);
    dist_rfFromSacc_OFF(j,:) = rgb_dist(min_dist_idx,:);
end
    
    
%% Masking factors
%     j = 1
%     h1 = figure;
%   
%     combinedImage = imfuse(rf_plot_ON(:,:,i),uint8(mask_regions_ON{i}{regIdx_flash_ON(i)}),'blend');
%     imagesc(combinedImage);colormap(gray);
% 
t = -pi:0.01:pi;
fac_sig_inc = .01;
fac_sig_all = [0.1:fac_sig_inc:2];

maskingFac_ON = [];
maskingFac_OFF = [];
    
parfor j = 1:length(uname_ON)
    r_x = fac_sig_all*rf_fit_sigma_ON(1,j);
    r_y = fac_sig_all*rf_fit_sigma_ON(2,j);  
    transformMatrix = [cosd(rf_fit_ellipseAngle_ON(j)),sind(rf_fit_ellipseAngle_ON(j));-sind(rf_fit_ellipseAngle_ON(j)),cosd(rf_fit_ellipseAngle_ON(j))]; % elipse roation

    ellipse_coord = {};

    for i = 1:length(fac_sig_all)
        x = r_x(i).*cos(t');
        y = r_y(i).*sin(t');
        rf_rotated = [x,y]*transformMatrix;
        ellipse_coord{i} = bsxfun(@plus,rf_rotated,rf_fit_center_ON(:,j)');

%             plot(ellipse_coord{i}(:,1),ellipse_coord{i}(:,2),'linewidth',1)
%             set(gca,'Ydir','reverse')
    end
    
    
    img_mask_withinMask = mask_regions_ON{j}{regIdx_flash_ON(j)};
    img_segments = {};
    img_segments{1} = poly2mask(ellipse_coord{1}(:,1),ellipse_coord{1}(:,2),size(img_mask_withinMask,1),size(img_mask_withinMask,2));

    for i = 2:length(fac_sig_all)
        img_segments{i} = poly2mask(ellipse_coord{i}(:,1),ellipse_coord{i}(:,2),size(img_mask_withinMask,1),size(img_mask_withinMask,2)) - poly2mask(ellipse_coord{i-1}(:,1),ellipse_coord{i-1}(:,2),size(img_mask_withinMask,1),size(img_mask_withinMask,2));
    end
    
    perc_sacc_all_withinMask = [];
    for i = 1:length(fac_sig_all)
        numPixs_intersect_withinMask = sum(sum(img_mask_withinMask&img_segments{i}));
        perc_sacc_all_withinMask(i) = 100*(numPixs_intersect_withinMask/sum(sum(img_segments{i})));
    end
%     figure;plot(perc_sacc_all_withinMask,'o')
%     figure;plot(diff(perc_sacc_all_withinMask),'o')
    
    [~,max_diff] = max(diff(perc_sacc_all_withinMask));
    rgb = sub2ind(size(perc_sacc_all_withinMask),max_diff);
    perc_sacc_withinMask = perc_sacc_all_withinMask(rgb);

    fac_sig_withinMask = (max_diff*fac_sig_inc);
    maskingFac_ON(j) = fac_sig_withinMask;
end

parfor j = 1:length(uname_OFF)
    r_x = fac_sig_all*rf_fit_sigma_OFF(1,j);
    r_y = fac_sig_all*rf_fit_sigma_OFF(2,j);  
    transformMatrix = [cosd(rf_fit_ellipseAngle_OFF(j)),sind(rf_fit_ellipseAngle_OFF(j));-sind(rf_fit_ellipseAngle_OFF(j)),cosd(rf_fit_ellipseAngle_OFF(j))]; % elipse roation

    ellipse_coord = {};

    for i = 1:length(fac_sig_all)
        x = r_x(i).*cos(t');
        y = r_y(i).*sin(t');
        rf_rotated = [x,y]*transformMatrix;
        ellipse_coord{i} = bsxfun(@plus,rf_rotated,rf_fit_center_OFF(:,j)');

%             plot(ellipse_coord{i}(:,1),ellipse_coord{i}(:,2),'linewidth',1)
%             set(gca,'Ydir','reverse')
    end
    
    
    img_mask_withinMask = mask_regions_OFF{j}{regIdx_flash_OFF(j)};
    img_segments = {};
    img_segments{1} = poly2mask(ellipse_coord{1}(:,1),ellipse_coord{1}(:,2),size(img_mask_withinMask,1),size(img_mask_withinMask,2));

    for i = 2:length(fac_sig_all)
        img_segments{i} = poly2mask(ellipse_coord{i}(:,1),ellipse_coord{i}(:,2),size(img_mask_withinMask,1),size(img_mask_withinMask,2)) - poly2mask(ellipse_coord{i-1}(:,1),ellipse_coord{i-1}(:,2),size(img_mask_withinMask,1),size(img_mask_withinMask,2));
    end
    
    perc_sacc_all_withinMask = [];
    for i = 1:length(fac_sig_all)
        numPixs_intersect_withinMask = sum(sum(img_mask_withinMask&img_segments{i}));
        perc_sacc_all_withinMask(i) = 100*(numPixs_intersect_withinMask/sum(sum(img_segments{i})));
    end
%     figure;plot(perc_sacc_all_withinMask,'o')
%     figure;plot(diff(perc_sacc_all_withinMask),'o')
    
    [~,max_diff] = max(diff(perc_sacc_all_withinMask));
    rgb = sub2ind(size(perc_sacc_all_withinMask),max_diff);
    perc_sacc_withinMask = perc_sacc_all_withinMask(rgb);

    fac_sig_withinMask = (max_diff*fac_sig_inc);
    maskingFac_OFF(j) = fac_sig_withinMask;
end

%% Calculate region idx

    regIdx_sacc_flash_ON = 5*ones(size(regIdx_flash_ON));regIdx_sacc_flash_ON(regIdx_flash_ON==3) = 4;
    regIdx_sacc_flash_OFF = 5*ones(size(regIdx_flash_OFF));regIdx_sacc_flash_OFF(regIdx_flash_OFF==3) = 4;
    
    regs_cat_ON = [ones(size(regIdx_flash_ON));regIdx_flash_ON;regIdx_sacc_flash_ON];
    regs_cat_OFF = [ones(size(regIdx_flash_OFF));regIdx_flash_OFF;regIdx_sacc_flash_OFF];
    
    
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
    
    p_avg_flash_ON = nan(length(flash_onset)-1,length(conds_mask),length(maskSize),length(sc),length(drugs),length(uname_ON));
    pow_p_avg_flash_ON = nan(length(flash_onset)-1,length(conds_mask),length(maskSize),length(sc),length(drugs),length(uname_ON));
    sampNum_p_avg_flash_ON = nan(length(flash_onset)-1,length(conds_mask),length(maskSize),length(sc),length(drugs),length(uname_ON));
    num_sacc_ON = [];
    
    
    p_avg_flash_OFF = nan(length(flash_onset)-1,length(conds_mask),length(maskSize),length(sc),length(drugs),length(uname_OFF));
    pow_p_avg_flash_OFF = nan(length(flash_onset)-1,length(conds_mask),length(maskSize),length(sc),length(drugs),length(uname_OFF));
    sampNum_p_avg_flash_off = nan(length(flash_onset)-1,length(conds_mask),length(maskSize),length(sc),length(drugs),length(uname_OFF));
    num_sacc_OFF = [];
    
    idx_onCells = false(length(flash_onset)-1,length(conds_mask),length(maskSize),length(sc),length(drugs),length(uname_ON));
    idx_offCells = false(length(flash_onset)-1,length(conds_mask),length(maskSize),length(sc),length(drugs),length(uname_OFF));
    
    for n = 1:size(mi_allSacc_flash_OFF,7)      % cells
        for m = 1:size(mi_allSacc_flash_OFF,6)      % drugs
            for l = 1:size(mi_allSacc_flash_OFF,5)      % sc
                for k = 1:size(mi_allSacc_flash_OFF,4)  % mask size
                    for j = 1:size(mi_allSacc_flash_OFF,3)  % conds_mask
                        for i = 1:size(mi_allSacc_flash_OFF,1)  % time point
                            rgb = squeeze(mi_allSacc_flash_OFF(i,:,j,k,l,m,n));

                            idx_toRemove = isnan(rgb) | rgb==inf | rgb==-inf;
                            rgb(idx_toRemove) = [];
                            num_sacc_OFF(i,j,k,l,m,n) = length(rgb);
                             if ~isempty(rgb)
                                idx_offCells(i,j,k,l,m,n) = true;
                                if mi_avg_flash_OFF(i,j,k,l,m,n) < 0
                                    [p_avg_flash_OFF(i,j,k,l,m,n),~,stats] = signtest(rgb,0,'tail','left');
                                    pow_p_avg_flash_OFF(i,j,k,l,m,n) = round(100*sampsizepwr('p',0.50,(stats.sign/num_sacc_OFF(i,j,k,l,m,n))+eps,[],num_sacc_OFF(i,j,k,l,m,n),'tail','left'));        
                                else 
                                    [p_avg_flash_OFF(i,j,k,l,m,n),~,stats] = signtest(rgb,0,'tail','right');
                                    pow_p_avg_flash_OFF(i,j,k,l,m,n) = round(100*sampsizepwr('p',0.50,(stats.sign/num_sacc_OFF(i,j,k,l,m,n))-eps,[],num_sacc_OFF(i,j,k,l,m,n),'tail','right'));   % 1-binocdf(critcL,num_sacc_OFF(i,j,k),stats.sign/num_sacc_OFF(i,j,k)) where critclL = 25;

                                end

                             end 
                        end
                    end
                end
            end
        end
    end
    

    p_avg_flash_ON = nan(length(flash_onset)-1,length(conds_mask),length(maskSize),length(sc),length(drugs),length(uname_ON));
    pow_p_avg_flash_ON = nan(length(flash_onset)-1,length(conds_mask),length(maskSize),length(sc),length(drugs),length(uname_ON));
    sampNum_p_avg_flash_ON = nan(length(flash_onset)-1,length(conds_mask),length(maskSize),length(sc),length(drugs),length(uname_ON));
    num_sacc_ON = [];
    
    
    p_avg_flash_ON = nan(length(flash_onset)-1,length(conds_mask),length(maskSize),length(sc),length(drugs),length(uname_ON));
    pow_p_avg_flash_ON = nan(length(flash_onset)-1,length(conds_mask),length(maskSize),length(sc),length(drugs),length(uname_ON));
    sampNum_p_avg_flash_off = nan(length(flash_onset)-1,length(conds_mask),length(maskSize),length(sc),length(drugs),length(uname_ON));
    num_sacc_ON = [];
    
    idx_onCells = false(length(flash_onset)-1,length(conds_mask),length(maskSize),length(sc),length(drugs),length(uname_ON));
    idx_offCells = false(length(flash_onset)-1,length(conds_mask),length(maskSize),length(sc),length(drugs),length(uname_ON));
    
    for n = 1:size(mi_allSacc_flash_ON,7)      % cells
        for m = 1:size(mi_allSacc_flash_ON,6)      % drugs
            for l = 1:size(mi_allSacc_flash_ON,5)      % sc
                for k = 1:size(mi_allSacc_flash_ON,4)  % mask size
                    for j = 1:size(mi_allSacc_flash_ON,3)  % conds_mask
                        for i = 1:size(mi_allSacc_flash_ON,1)  % time point
                            rgb = squeeze(mi_allSacc_flash_ON(i,:,j,k,l,m,n));

                            idx_toRemove = isnan(rgb) | rgb==inf | rgb==-inf;
                            rgb(idx_toRemove) = [];
                            num_sacc_ON(i,j,k,l,m,n) = length(rgb);
                             if ~isempty(rgb)
                                idx_offCells(i,j,k,l,m,n) = true;
                                if mi_avg_flash_ON(i,j,k,l,m,n) < 0
                                    [p_avg_flash_ON(i,j,k,l,m,n),~,stats] = signtest(rgb,0,'tail','left');
                                    pow_p_avg_flash_ON(i,j,k,l,m,n) = round(100*sampsizepwr('p',0.50,(stats.sign/num_sacc_ON(i,j,k,l,m,n))+eps,[],num_sacc_ON(i,j,k,l,m,n),'tail','left'));        
                                else 
                                    [p_avg_flash_ON(i,j,k,l,m,n),~,stats] = signtest(rgb,0,'tail','right');
                                    pow_p_avg_flash_ON(i,j,k,l,m,n) = round(100*sampsizepwr('p',0.50,(stats.sign/num_sacc_ON(i,j,k,l,m,n))-eps,[],num_sacc_ON(i,j,k,l,m,n),'tail','right'));   % 1-binocdf(critcL,num_sacc_ON(i,j,k),stats.sign/num_sacc_ON(i,j,k)) where critclL = 25;

                                end

                             end 
                        end
                    end
                end
            end
        end
    end  

%% Fig. 2f - Line Plots

% use ; to seperate conditions across rows
cutOff_saccArea = 15;   % in percentage
timePoints = 1:length(flash_onset)-1;
reg_select = [1;2];     % regions = {'flash+sacc in all','no sacc in center','sacc in center'};
maskSize_select = [1];     % [100,150]
drugs_select = [1];
sc_select = [2];  % [25,150]
ON_OFF_select = 'ALL';    % 'ON' 'OFF' 'ALL'

sc_text = cellstr(num2str(sc'));

plots_idx = 1;
lim_x = [-100,2100];
bins = [-1:0.1:1];

col_on = 'r';%[0.6,0.6,0.6];
col_off = 'b';%[0.2,0.2,0.2];
LINE_WIDTH = 2;
LINE_TYPES = {'-','--',':','-.'};


h_main=figure; suptitle(['ON & OFF Cells - sacc area cutoff: ',num2str(cutOff_saccArea),'%'])    
h_m = [];

s = 1; m = 1; d = 1; t = 0;

idx_valid_ON = squeeze(all(all(all(all(~isnan(mi_avg_flash_ON(:,reg_select,maskSize_select,sc_select,drugs_select,:)),2),3),4),5));
idx_valid_OFF = squeeze(all(all(all(all(~isnan(mi_avg_flash_OFF(:,reg_select,maskSize_select,sc_select,drugs_select,:)),2),3),4),5));

p_pop_ON = []; %nan(size(mi_avg_flash_ON,1),size(mi_avg_flash_ON,2),size(mi_avg_flash_ON,3),size(mi_avg_flash_ON,4));
p_pop_OFF = []; %nan(size(mi_avg_flash_OFF,1),size(mi_avg_flash_OFF,2),size(mi_avg_flash_OFF,3),size(mi_avg_flash_OFF,4));
pop_across_data_ON = {};
pop_across_data_OFF = {};

counter = 0;
for d = 1:size(drugs_select,1)
    for s = 1:size(sc_select,1)
        for z = 1:size(maskSize_select)
            for m = 1:size(reg_select,1)

                points_x_ON = [];     
                for j = 1:size(mi_avg_flash_ON,6)
                    points_x_ON(j,:) = squeeze(mi_avg_flash_ON(:,regs_cat_ON(reg_select(m),j),maskSize_select(z),sc_select(s),drugs_select(d),j));
                end

                points_x_OFF = [];     
                for j = 1:size(mi_avg_flash_OFF,6)
                    points_x_OFF(j,:) = squeeze(mi_avg_flash_OFF(:,regs_cat_OFF(reg_select(m),j),maskSize_select(z),sc_select(s),drugs_select(d),j));
                end


                idx_respDrugs_ON = all(respConds_ON(drugs_select,:),1); idx_respDrugs_ON(isnan(idx_respDrugs_ON)) = 0;
                idx_respDrugs_ON = logical(idx_respDrugs_ON);

                idx_respDrugs_OFF = all(respConds_OFF(drugs_select,:),1); idx_respDrugs_OFF(isnan(idx_respDrugs_OFF)) = 0;
                idx_respDrugs_OFF = logical(idx_respDrugs_OFF);

                if any(ismember(reg_select,2))   % condition where there should be no saccade in RF center
                    idx_toKeep_ON = repmat(idx_respDrugs_ON,size(points_x_ON,2),1) & ~isnan(points_x_ON') & repmat(regArea_sacc_ON,size(points_x_ON,2),1) <= cutOff_saccArea;
                    idx_toKeep_OFF = repmat(idx_respDrugs_OFF,size(points_x_OFF,2),1) & ~isnan(points_x_OFF') & repmat(regArea_sacc_OFF,size(points_x_OFF,2),1) <= cutOff_saccArea;
                else
                    idx_toKeep_ON = repmat(idx_respDrugs_ON,size(points_x_ON,2),1) & ~isnan(points_x_ON');
                    idx_toKeep_OFF = repmat(idx_respDrugs_OFF,size(points_x_OFF,2),1) & ~isnan(points_x_OFF');
                end


                idx_toKeep_ON = idx_toKeep_ON & idx_valid_ON;
                idx_toKeep_OFF = idx_toKeep_OFF & idx_valid_OFF;

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
                        errorbar(flash_onset,[mean_plot_x_ON,0],[sem_plot_x_ON,0],'Color',col_on,'LineWidth',LINE_WIDTH,'LineStyle',LINE_TYPES{counter},'DisplayName',['ON (',num2str(rgb_num_ON),')|',sc_text{sc_select(s)},'µm|',regions{reg_select(m)},'|',drugs{drugs_select(d)}]);
                    case 'OFF'
                        errorbar(flash_onset,[mean_plot_x_OFF,0],[sem_plot_x_OFF,0],'Color',col_off,'LineWidth',LINE_WIDTH,'LineStyle',LINE_TYPES{counter},'DisplayName',['OFF (',num2str(rgb_num_OFF),')|',sc_text{sc_select(s)},'µm|',regions{reg_select(m)},'|',drugs{drugs_select(d)}]);
                    otherwise
                        errorbar(flash_onset,[mean_plot_x_ON,0],[sem_plot_x_ON,0],'Color',col_on,'LineWidth',LINE_WIDTH,'LineStyle',LINE_TYPES{counter},'DisplayName',['ON (',num2str(rgb_num_ON),')|',sc_text{sc_select(s)},'µm|',regions{reg_select(m)},'|',drugs{drugs_select(d)}]);
                        errorbar(flash_onset,[mean_plot_x_OFF,0],[sem_plot_x_OFF,0],'Color',col_off,'LineWidth',LINE_WIDTH,'LineStyle',LINE_TYPES{counter},'DisplayName',['OFF (',num2str(rgb_num_OFF),')|',sc_text{sc_select(s)},'µm|',regions{reg_select(m)},'|',drugs{drugs_select(d)}]);
                end

                legend('-DynamicLegend');

                ylim([-0.6 0.3])
                ylabel('Modulation index')
                xlabel('Time from saccade onset (ms)')
                xlim(lim_x)
                title(['ON = ',num2str(nanmax(rgb_num_ON)),' | OFF = ',num2str(nanmax(rgb_num_OFF))])

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

% population significance
p_pop_ON = squeeze(p_pop_ON);
p_pop_OFF = squeeze(p_pop_OFF);

a = [flash_onset(timePoints);p_pop_ON']
b = [flash_onset(timePoints);p_pop_OFF']

pop_across_data_ON = squeeze(pop_across_data_ON);
pop_across_data_OFF = squeeze(pop_across_data_OFF);

p_across_ON = cell2mat(cellfun(@(x,y) ranksum(x,y,'tail','both'),pop_across_data_ON(:,1),pop_across_data_ON(:,2),'UniformOutput',0));
p_across_OFF = cell2mat(cellfun(@(x,y) ranksum(x,y,'tail','both'),pop_across_data_OFF(:,1),pop_across_data_OFF(:,2),'UniformOutput',0));
p_across = [flash_onset(timePoints);p_across_ON';p_across_OFF']


%% Supp Fig. 7c - Scatter Plots
    % colums for x and y axis of scatter | rows for another row of scatter plot.
    plotMaskingFacs = 0;
    cutOff_saccArea = 15;   % in percentage
    timePoints = 1:length(flash_onset)-1;
    reg_select = [1,2];     % regions = {'flash+sacc in all','no sacc in center','sacc in center'};
    maskSize_select = [1,1];     % [100,150]
    drugs_select = [1,1];
    sc_select = [2,2];  % [25,150]
    
    sc_text = cellstr(num2str(sc'));

    plots_idx = [1:(length(timePoints)*size(reg_select,1)*size(maskSize_select,1)*size(sc_select,1)*size(drugs_select,1))]';
    plots_idx = reshape(plots_idx',length(timePoints),size(reg_select,1),size(maskSize_select,1),size(sc_select,1),size(drugs_select,1));        % [time,mask_conds,sc]
    
    lim_x = [-1.2,1.2];
    lim_y = [-1.2,1.2];
    bins = [-1:0.1:1];
    baseVal = -1.2;
    
    col_on = [0.6,0.6,0.6];
    col_off = [0.2,0.2,0.2];

    h_main=figure; suptitle(['Scatter Plots - sacc area cutoff: ',num2str(cutOff_saccArea),'%'])    
    h_m = [];
    
    s = 1; m = 1; d = 1; t = 1;

    for d = 1:size(drugs_select,1)
        for s = 1:size(sc_select,1)
            for z = 1:size(maskSize_select)
                for m = 1:size(reg_select,1)
                    for t = timePoints;

                        points_x_ON = []; points_y_ON = [];
                        for j = 1:size(mi_avg_flash_ON,6)
                            points_x_ON(j,:) = squeeze(mi_avg_flash_ON(t,regs_cat_ON(reg_select(m,1),j),maskSize_select(z,1),sc_select(s,1),drugs_select(d,1),j));
                            points_y_ON(j,:) = squeeze(mi_avg_flash_ON(t,regs_cat_ON(reg_select(m,2),j),maskSize_select(z,2),sc_select(s,2),drugs_select(d,2),j));
                        end

                        points_x_OFF = []; points_y_OFF = [];
                        for j = 1:size(mi_avg_flash_OFF,6)
                            points_x_OFF(j,:) = squeeze(mi_avg_flash_OFF(t,regs_cat_OFF(reg_select(m,1),j),maskSize_select(z,1),sc_select(s,1),drugs_select(d,1),j));
                            points_y_OFF(j,:) = squeeze(mi_avg_flash_OFF(t,regs_cat_OFF(reg_select(m,2),j),maskSize_select(z,2),sc_select(s,2),drugs_select(d,2),j));
                        end

                        h_m(t,m,z,s,d) = subplot(size(reg_select,1)*size(reg_select,1)*size(sc_select,1)*size(drugs_select,1),length(timePoints),plots_idx(t,m,z,s,d));
                        hold on

                        idx_respDrugs_ON = respConds_ON(drugs_select(d),:); idx_respDrugs_ON(isnan(idx_respDrugs_ON)) = 0;
                        idx_respDrugs_ON = logical(idx_respDrugs_ON);

                        idx_respDrugs_OFF = respConds_OFF(drugs_select(d),:); idx_respDrugs_OFF(isnan(idx_respDrugs_OFF)) = 0;
                        idx_respDrugs_OFF = logical(idx_respDrugs_OFF);

                        if reg_select(m,1) || reg_select(m,2) == 2   % condition where there should be no saccade in RF center
                            idx_toKeep_ON = idx_respDrugs_ON & ~isnan(points_x_ON') & ~isnan(points_y_ON') & regArea_sacc_ON <= cutOff_saccArea;
                            idx_toKeep_OFF = idx_respDrugs_OFF & ~isnan(points_x_OFF') & ~isnan(points_y_OFF') & regArea_sacc_OFF <= cutOff_saccArea;
                        else
                            idx_toKeep_ON = idx_respDrugs_ON & ~isnan(points_x_ON') & ~isnan(points_y_ON');
                            idx_toKeep_OFF = idx_respDrugs_OFF & ~isnan(points_x_OFF') & ~isnan(points_y_OFF');
                        end

                        plot_x_ON  = points_x_ON(idx_toKeep_ON);
                        plot_y_ON  = points_y_ON(idx_toKeep_ON);

                        plot_x_OFF  = points_x_OFF(idx_toKeep_OFF);
                        plot_y_OFF  = points_y_OFF(idx_toKeep_OFF);

                        plots_uname_ON = uname_ON(idx_toKeep_ON');
                        plots_uname_OFF = uname_OFF(idx_toKeep_OFF');
                        
                        if plotMaskingFacs == 1
                            assert(reg_select(m,1)==reg_select(m,2) & maskSize_select(z,1)==maskSize_select(z,2) & sc_select(s,1)==sc_select(s,2) & drugs_select(d,1)==drugs_select(d,2),'combinations dont match')
                            plot_x_ON = maskingFac_ON(idx_toKeep_ON);
                            plot_x_OFF = maskingFac_OFF(idx_toKeep_OFF);
                            
                            baseVal = 0;
                        end

                        scatter(plot_x_ON,plot_y_ON,'MarkerEdgeColor',col_on+.01,'MarkerFaceColor',col_on+.01) 
                        scatter(plot_x_OFF,plot_y_OFF,'MarkerEdgeColor',col_off+.01,'MarkerFaceColor',col_off+.01)


                        [f_off_all,x] = hist(plot_y_OFF,bins);
                        h = barh(x,(f_off_all/trapz(f_off_all))+baseVal);
                        set(h,'FaceColor',col_off,'EdgeColor','k');
                        h.BaseValue = baseVal;

                        [f_on_all,x] = hist(plot_y_ON,bins);
                        h = barh(x,(f_on_all/trapz(f_on_all))+baseVal);
                        set(h,'FaceColor',col_on,'EdgeColor','k');
                        h.BaseValue = baseVal;
                        set(h,'FaceAlpha',0.6)

                        if plotMaskingFacs ~= 1
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
                        else
                            max_x = max([plot_x_ON,plot_x_OFF])+.1;
                            plot([0,max_x],[0,0],'--','color',[0.8,0.8,0.8])
                            xlim([0,max_x])
                            ylim(lim_y)
                        end

                        axis square

                    end
                    legend({['ON ',num2str(length(plot_x_ON))],['OFF: ',num2str(length(plot_x_OFF))]})

                    if plotMaskingFacs ~= 1
                        xlabel(['sc: ',sc_text{sc_select(s,1)},' µm | ', regions{reg_select(m,1)},' | ',drugs{drugs_select(d,1)}])
                        ylabel(['sc: ',sc_text{sc_select(s,2)},' µm | ', regions{reg_select(m,2)},' | ',drugs{drugs_select(d,2)}])
                    else
                        xlabel(['masking factor'])
                        ylabel(['sc: ',sc_text{sc_select(s,2)},' µm | ', regions{reg_select(m,2)},' | ',drugs{drugs_select(d,2)}])
                    end

                end
            end
        end
    end
    
    for i = timePoints
        title(h_m(i,1,1,1),[num2str(flash_onset(i)),'ms'])
    end
    
    
%% plot RFs using mask_facs

    % use ; to seperate conditions across rows
    cutOff_saccArea = 15;   % in percentage
    timePoints = 1:length(flash_onset)-1;
    reg_select = [1;2]; % regions = {'flash+sacc in all','no sacc in center','sacc in center'};
    maskSize_select = [1];  % [100,150]
    drugs_select = [1]; % [no drugs]
    sc_select = [2];    % [25,150];
    sc_text = cellstr(num2str(sc'));

   
    plots_idx = [1,2;3,4];
    
    col_on = [0.6,0.6,0.6];
    col_off = [0.2,0.2,0.2];
    
    s = 1; m = 1; d = 1; t = 1;
    cat_idx_toKeep_ON = [];
    cat_idx_toKeep_OFF = [];

    for s = 1:size(sc_select,1)
        for z = 1:size(maskSize_select,1)
            for m = 1:size(reg_select,1)
                for d = 1:size(drugs_select,1)
                    for t = timePoints

                        points_x_ON = [];     
                        for j = 1:size(mi_avg_flash_ON,6)
                            points_x_ON(j,:) = squeeze(mi_avg_flash_ON(t,regs_cat_ON(reg_select(m,1),j),maskSize_select(z),sc_select(s),drugs_select(d),j));
                        end

                        points_x_OFF = [];     
                        for j = 1:size(mi_avg_flash_OFF,6)
                            points_x_OFF(j,:) = squeeze(mi_avg_flash_OFF(t,regs_cat_OFF(reg_select(m,1),j),maskSize_select(z),sc_select(s),drugs_select(d),j));
                        end



                        idx_respDrugs_ON = respConds_ON(drugs_select(d),:); idx_respDrugs_ON(isnan(idx_respDrugs_ON)) = 0;
                        idx_respDrugs_ON = logical(idx_respDrugs_ON);

                        idx_respDrugs_OFF = respConds_OFF(drugs_select(d),:); idx_respDrugs_OFF(isnan(idx_respDrugs_OFF)) = 0;
                        idx_respDrugs_OFF = logical(idx_respDrugs_OFF);

                        if m == 2   % condition where there should be no saccade in RF center
                            idx_toKeep_ON = idx_respDrugs_ON & ~isnan(points_x_ON') & regArea_sacc_ON <= cutOff_saccArea;
                            idx_toKeep_OFF = idx_respDrugs_OFF & ~isnan(points_x_OFF') & regArea_sacc_OFF <= cutOff_saccArea;
                        else
                            idx_toKeep_ON = idx_respDrugs_ON & ~isnan(points_x_ON');
                            idx_toKeep_OFF = idx_respDrugs_OFF & ~isnan(points_x_OFF');
                        end
                        
                        cat_idx_toKeep_ON(:,t,m,s,d) = idx_toKeep_ON;
                        cat_idx_toKeep_OFF(:,t,m,s,d) = idx_toKeep_OFF;

                    end

                end
            end
        end
    end
    
  % ON RGCs
    pix2um = 1;%pix2um_ON(1);

    hypoth_win = [0,0,800,600];
    start_corn = [400,300];
    
    img_mask_withinMask = mask_regions_ON{1}{regIdx_sacc_ON(1)};
   
    mf_idx_toKeep_ON = squeeze(all(all(all(all(cat_idx_toKeep_ON,2),3),4),5));
    mf_idx_toKeep_ON(40) = false;
    mf_idx_ON = find(mf_idx_toKeep_ON);
    
    t = -pi:0.01:pi;
    facs = [1];
    figure;hold on
    imagesc(img_mask_withinMask);colormap(gray);hold on
    rectangle('position',hypoth_win);


    for u = 1:length(mf_idx_ON)
        rgb_sigma_x = rf_fit_sigma_ON(1,mf_idx_ON(u));
        rgb_sigma_y = rf_fit_sigma_ON(2,mf_idx_ON(u));
        sig_rat = rgb_sigma_y/rgb_sigma_x;
        
        if sig_rat > 2
            rgb_sigma_y = 2*rgb_sigma_x;
        end
        
        x = facs*rgb_sigma_x.*cos(t')*pix2um;
        y = facs*rgb_sigma_y.*sin(t')*pix2um;
        transformMatrix = [cosd(rf_fit_ellipseAngle_ON(mf_idx_ON(u))),sind(rf_fit_ellipseAngle_ON(mf_idx_ON(u)));-sind(rf_fit_ellipseAngle_ON(mf_idx_ON(u))),cosd(rf_fit_ellipseAngle_ON(mf_idx_ON(u)))]; % elipse roation

        rf_rotated = [x,y]*transformMatrix;
        center_relative = dist_rfFromSacc_ON(mf_idx_ON(u),:)'*pix2um;
        ellipse_coord= bsxfun(@plus,rf_rotated,center_relative'+start_corn);

    h_l = plot(ellipse_coord(:,1),ellipse_coord(:,2),'linewidth',1);
    set(gca,'Ydir','reverse')

    end
    plot(start_corn(1),start_corn(2),'wx')
    xlim([-10 810])
    ylim([-10 610])
    title([num2str(u),' ON Cells'])
        
    % OFF RGCS
    mf_idx_toKeep_OFF = squeeze(all(all(all(all(cat_idx_toKeep_OFF,2),3),4),5));
    mf_idx_OFF = find(mf_idx_toKeep_OFF);
    
    figure;hold on
    imagesc(img_mask_withinMask);colormap(gray);hold on
    rectangle('position',hypoth_win);


    for u = 1:length(mf_idx_OFF)
        
        rgb_sigma_x = rf_fit_sigma_OFF(1,mf_idx_OFF(u));
        rgb_sigma_y = rf_fit_sigma_OFF(2,mf_idx_OFF(u));
        sig_rat = rgb_sigma_y/rgb_sigma_x;
        
        if sig_rat > 2
            rgb_sigma_y = 2*rgb_sigma_x;
        end
        
        x = facs*rgb_sigma_x.*cos(t')*pix2um;
        y = facs*rgb_sigma_y.*sin(t')*pix2um;
        transformMatrix = [cosd(rf_fit_ellipseAngle_OFF(mf_idx_OFF(u))),sind(rf_fit_ellipseAngle_OFF(mf_idx_OFF(u)));-sind(rf_fit_ellipseAngle_OFF(mf_idx_OFF(u))),cosd(rf_fit_ellipseAngle_OFF(mf_idx_OFF(u)))]; % elipse roation

        rf_rotated = [x,y]*transformMatrix;
        center_relative = dist_rfFromSacc_OFF(mf_idx_OFF(u),:)'*pix2um;
        ellipse_coord= bsxfun(@plus,rf_rotated,center_relative'+start_corn);

    h_l = plot(ellipse_coord(:,1),ellipse_coord(:,2),'linewidth',1);
    set(gca,'Ydir','reverse')

    end
    plot(start_corn(1),start_corn(2),'wx')
    xlim([-10 810])
    ylim([-10 610])
    title([num2str(u),' OFF Cells'])
        
                    
                    
                    
                    
                    
                    
                    



