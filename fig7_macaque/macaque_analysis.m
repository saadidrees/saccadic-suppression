clear all; clc; close all; format long g;set(0,'defaultTextInterpreter','none');

load data_macaque.mat

colors_plots = [1,0,1;0,0,1;0,1,0;1,0,0;0.2,0.5,0.8;0,1,1;0.8,0.8,0;0.7,0.2,0.2];
% delays = [17,50,100,250,500,2000] % flash onset

%% population stats

for i = 1:size(mi_neg_ON,1)
    p_neg_ON(i) = signrank(mi_neg_ON(i,:),0);
    p_pos_ON(i) = signrank(mi_pos_ON(i,:),0);
    p_neg_OFF(i) = signrank(mi_neg_OFF(i,:),0);
    p_pos_OFF(i) = signrank(mi_pos_OFF(i,:),0);
end
    
%% plotting
    
    contr_type = {'Negative','Positive'};

    figure;
    for i = 1:2
        subplot(2,1,i); hold on
        if i ==1
            plot(delays,mi_neg_ON,'ro')
        elseif i==2
            plot(delays,mi_pos_ON,'ro')
        end

        errorbar(delays,mi_avg_ON(:,i),mi_std_ON(:,i),'-ro')
        errorbar(delays,mi_avg_OFF(:,i),mi_std_OFF(:,i),'-bo')

        plot([0,2000],[0,0],'k');
        legend({['ON RGCs N = ',num2str(num_ON(i))],['OFF RGCs N = ',num2str(num_OFF(i))]})
        xlabel('Time of flash after step (ms)')
        ylabel('Modulation index')
        title([contr_type{i},'-contrast luminance step'])
        ylim([-1 1])
    end
        
        
        


