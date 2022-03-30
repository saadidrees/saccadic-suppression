%% Fig 6d
% loads the data file containing model RGC modulation index for flashes presented after luminance steps
% modulation index is based on response to flash presented 2000 ms after the luminance step

%%
clear; clc;

load data_model_modulationIndex.mat

figure; suptitle('Model RGCs modulation index as a function of flash time')
select_tau = 1%find(tr_filt_bilobe == 50.64);
select_thresh = find(outputThresholdFast==0.1);

for c = 1:length(contr_arr)
    subplot(2,1,c);hold on

        plot(flash_onset(1:end-1),squeeze(mi_flash_OFF(1:end-1,c,select_tau,select_thresh)),'-x','Color','b','DisplayName',['OFF | ',num2str(contr_arr(c))]);        % OFF
        plot(flash_onset(1:end-1),squeeze(mi_flash_ON(1:end-1,c,select_tau,select_thresh)),'-x','Color','r','DisplayName',['ON | ',num2str(contr_arr(c))]);        % ON
%         plot([preSacc_dur_all{c},preSacc_dur_all{c}],[-1,1],'k')
        plot([0,flash_onset(end-1)],[0,0],'k')
        xlabel('Time of flash from luminance step (ms)')
        ylabel('Modulation index')
        title(['Contrast step: ',num2str(contr_arr(c))])
        legend('-DynamicLegend');

end

