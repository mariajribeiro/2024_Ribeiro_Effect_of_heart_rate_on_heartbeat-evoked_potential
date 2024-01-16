% plot and compare RT interval

load RT_interval

 plot_all_data_2conds_2grps(RT_interval{1, 1}(:, 2), RT_interval{1, 2}(:, 2), RT_interval{2, 1}(:, 2), RT_interval{2, 2}(:, 2),...
     'RT interval (ms)', '', {'Young', 'Older'})
 
 %% save data in excel file for spss analyses
 
 % create table
% group, RT_interv_simpleRT, RT_interv_gng
% IBI_avg_std{group, 2}(old, t, :) = [mean(IBI), std(IBI)];
RT_interval_table = table([RT_interval{1, 1}(:, 1); RT_interval{2, 1}(:, 1)], [ones(length(RT_interval{1, 1}), 1);ones(length(RT_interval{2, 1}), 1)*2 ], ...
    [RT_interval{1, 1}(:, 2); RT_interval{2, 1}(:, 2)], [RT_interval{1, 2}(:, 2); RT_interval{2, 2}(:, 2)],...
    'VariableNames',{'participant' 'group', 'RT_interv_simpleRT', 'RT_interv_gng'}); 

filename = 'RT_interval.xlsx';
writetable(RT_interval_table,filename) 
 

%% plot data 2 conditions

function plot_all_data_2conds_2grps(grp1_cond1, grp1_cond2, grp2_cond1, grp2_cond2,...
    y_label_text, title_text, x_label_text)
 
figure;
    for y=1:length(grp1_cond2)
        plot([1+rand*0.2-0.1, 2+rand*0.2-0.1], [grp1_cond1(y) grp1_cond2(y)] ,'-o',  'color', [.5 .5 .5], ...
        'MarkerEdgeColor',[.5 .5 .5],'MarkerSize',12, 'LineWidth', 1.5);
        hold on;
    end
    
    % plot data for group 1 cond1
    yMean1=nanmean(grp1_cond1);
    y_se1 = std(grp1_cond1,'omitnan')/sqrt(length(grp1_cond1));
    
        % plot data for group 1 cond2
    yMean2=nanmean(grp1_cond2);
    y_se2 = std(grp1_cond2,'omitnan')/sqrt(length(grp1_cond2));
    
    errorbar([1, 2], [yMean1, yMean2],[y_se1, y_se2], '-s', 'color', 'k','MarkerSize',12, 'MarkerFaceColor', 'k',  'LineWidth', 3);
    
    % group 2
    for y=1:length(grp2_cond2)
        plot([4+rand*0.2-0.1, 5+rand*0.2-0.1], [grp2_cond1(y) grp2_cond2(y)] ,'-o',  'color', [.5 .5 .5], ...
        'MarkerEdgeColor',[.5 .5 .5],'MarkerSize',12, 'LineWidth', 1.5);
        hold on;
    end
    
    % plot data for group 2 cond1
    yMean1=nanmean(grp2_cond1);
    y_se1 = std(grp2_cond1,'omitnan')/sqrt(length(grp2_cond1));
    
        % plot data for group 2 cond2
    yMean2=nanmean(grp2_cond2);
    y_se2 = std(grp2_cond2,'omitnan')/sqrt(length(grp2_cond2));
    
    errorbar([4, 5], [yMean1, yMean2],[y_se1, y_se2], '-s', 'color', 'k','MarkerSize',12, 'MarkerFaceColor', 'k',  'LineWidth', 3);
    
    % axes('XColor','none');
    hold off;
    box off;
    axis([0 6 200 340]);
    ax = gca;
    ax.LineWidth = 2.5; 
    ax.FontSize = 24;
    ax.FontName = 'Arial';
    ax.Color = 'none';
%     ax.XTickLabel= [];
    ax.XAxis.FontSize = 26;
    xticks([1.5 4.5])
    ax.XTickLabel= x_label_text;

    ylabel(y_label_text, 'FontSize', 30, 'FontWeight','normal')
    title(title_text, 'FontSize', 32, 'FontWeight','normal')

%     x0=10;
%     y0=10;
%     width=400;
%     height=400;
%     set(gcf,'position',[x0,y0,width,height])
end


