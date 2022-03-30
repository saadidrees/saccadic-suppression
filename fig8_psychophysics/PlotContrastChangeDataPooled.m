%% Fig. 8

load('PercentCorrectDataContrastChange')
times=[-24,-12,36,72,108];
lumdif=[-98   -84   -70   -56   -42   -28   -14 14    28    42    56    70    84    98];
abslumchange=abs(lumdif*0.56);
color = {[0.4940 0.1840 0.5560],[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880]};
greys = (linspace(60,200,7))'./256;
grey = repmat(greys,1,3);

    
    figure
    orderCS=[1 2 4 3 ];
    for kkk=[1,4] % the conditions in the paper
        
        subplot(2,2,orderCS(kkk))
        if kkk<3
            a1=area([-40 0],[1 1]);
            hold on
            a2=area([0 120],[50 50]);
            a2.FaceAlpha=0.2;
            a1.FaceAlpha=0.2;
            a2.FaceColor=[.4 .4 .4];
            a1.FaceColor=[.9 .9 .9];
            
        else
            a1=area([-40 0],[1 1]);
            hold on
            a2=area([0 120],[50 50]);
            a2.FaceAlpha=0.2;
            a1.FaceAlpha=0.2;
            a1.FaceColor=[.4 .4 .4];
            a2.FaceColor=[.9 .9 .9];
            
        end
             if kkk==1
                title('Contrast Steps between -0.56 and -0.3 Michelson Contrast','FontSize',12 )
             elseif kkk==2
                title('Contrast Steps between -0.3 and -0.053 Michelson Contrast','FontSize',12 )
             elseif kkk==4
                title('Contrast Steps between 0.053 and 0.3 Michelson Contrast','FontSize',12 )
                xlabel('Time [ms]')
             elseif kkk==3
                title('Contrast Steps between 0.3 and 0.56 Michelson Contrast','FontSize',12 )
        
             end
        errorbar(times(1:2),cell2mat(percentCorrectCS(kkk,1:2,1)),...
        cell2mat(percentCorrectCS(kkk,1:2,1))-cell2mat(cilowCS(kkk,1:2,1)),...
        cell2mat(ciupCS(kkk,1:2,1))-cell2mat(percentCorrectCS(kkk,1:2,1)),'-xk','LineWidth',2)
        hold on
        darkflash=plot(times(3:5),cell2mat(percentCorrectCS(kkk,3:5,1)),'-xk')
        shadedplot(times(3:5),cell2mat(cilowCS(kkk,3:5,1)),cell2mat(ciupCS(kkk,3:5,1)),'k','k')
        hold on
        alpha(0.3)
        errorbar(times(1:2),cell2mat(percentCorrectCS(kkk,1:2,2)),...
            cell2mat(percentCorrectCS(kkk,1:2,2))-cell2mat(cilowCS(kkk,1:2,2)),...
            cell2mat(ciupCS(kkk,1:2,2))-cell2mat(percentCorrectCS(kkk,1:2,2)),'-xb','LineWidth',2)
        hold on
        brightflash=plot(times(3:5),cell2mat(percentCorrectCS(kkk,3:5,2)),'-xb')
        shadedplot(times(3:5),cell2mat(cilowCS(kkk,3:5,2)),cell2mat(ciupCS(kkk,3:5,2)),'b','b')
        alpha(0.3)
        
        
            if kkk==1 || kkk==4
                ylabel('Percent Correct Trials')
            elseif kkk==3 || kkk==4
                xlabel('Time [ms]')
            end
        
        ylim([0 1])
        xlim([-30 115])
            if kkk==3
            legend([brightflash, darkflash],'Stim Increase','Stim Decrease')
            end
    end