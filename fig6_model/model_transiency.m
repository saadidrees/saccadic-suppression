%% Fig 6g and 6h
% loads the data file containing model RGC modulation indices for flashes presented after luminance step for different
% model parameters: for 0.5<tau_d<15 --> in paper this is normalized. 0.5 tau_f is 1 and 15 is taken as 0.
% 0<outputThreshold<1. This is not a normalized scale.


%%

clear; clc;

load data_model_transiency.mat

%% Fig 6g-h - scatter plot transiency

select_fig = '6h'   %{'6g','6h'}

switch select_fig
    case '6g' % keep nonlinearity fixed at 0.1 i.e. index 2
        select_tr_thresh = find(outputThresholdFast==0.1);
        select_tr_filt = [1:length(tr_filt_bilobe)];
        select_tr_filt = fliplr(select_tr_filt);
        
    case '6h' % as a function of nonlinearity. Keep filter transiency fixed at the sustained setting i.e. last index of tr_filt_bilobe
        select_tr_thresh = 1:length(outputThresholdFast);
        select_tr_filt = length(tr_filt_bilobe);
end

lims_y = [-1.1,0.5];
timePoints_select = [17,33,50]; % the flash time points to plot
[~,~,idx_tPoints] = intersect(timePoints_select,flash_onset);
idx_plots = 1:length(timePoints_select)*2;
idx_plots = reshape(idx_plots,length(timePoints_select),2)';
idx_plots = flipud(idx_plots);

h1 = figure; suptitle(select_fig)
    
    for c = 1:length(contr_arr) % loop over luminance step contrast. Top row is positive contrast
        counter = 0;
        for i = idx_tPoints'
            counter = counter+1;
            if length(select_tr_filt)>1 & length(select_tr_thresh)==1
                x_val = tr_filt_bilobe(select_tr_filt);
%                 x_val = fliplr(x_val);
                x_lim = [0, x_val(end)+1];
            elseif length(select_tr_thresh)>1 & length(select_tr_filt)==1
                x_val = outputThresholdFast(select_tr_thresh);
                x_lim = [0, x_val(end)+0.1];
            else
                error('select only one x variable')
            end
            subplot(2,length(timePoints_select),idx_plots(c,counter));hold on; title(['contrast step: ',num2str(contr_arr(c)),' | ',num2str(flash_onset(i)),' ms'])
                plot(x_val,squeeze(mi_flash_ON(i,c,select_tr_filt,select_tr_thresh)),'-co');
                plot(x_val,squeeze(mi_flash_OFF(i,c,select_tr_filt,select_tr_thresh)),'-ro');
%                 plot([0,x_val(end)],[0,0],'k')
%                 xlim(x_lim)
                ylim(lims_y)
                legend({'ON RGCs','OFF RGCs'})
                xlabel('Filter param')
                ylabel('Modulation index')
                if select_fig == '6g'
                    set(gca,'XDir','reverse')
                end
                axis square
        end
    end
        
%     info = 'mi_flash_ON dimensions: [flash_onset,contr_arr,tr_filt_bilobe,outputThresholdFast]';
%     save data_model_transiency.mat mi_flash_ON mi_flash_OFF outputThresholdFast tr_filt_bilobe flash_onset contr_arr contr_flash flash_type info

