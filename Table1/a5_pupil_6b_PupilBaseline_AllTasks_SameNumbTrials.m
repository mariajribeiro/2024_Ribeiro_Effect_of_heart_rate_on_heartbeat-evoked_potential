% effect of time on task - within-subject correlation analysis
% Maria Ribeiro June 2016 adapted october 2016
% changed 29 July2017 and again 12 May 2020

clear; close all
% younger=[4 9 10 13 15 16 25 26 28 31 33 34 36 42 44 45 46 50 51 53 54 56 59 62 66 68 72 74 76 78 80 81 82 84 85];
% older=[7 8	11	12	14	17	19	20	21	22	23	32	35	37	41	43	47 48 49 52 55	57	58	60	61	63	64  65	67	69	70	71	73	75  77  79 83   86];

% include only participants that were included in analyses of T-locked
% heartbeat-evoked potentials (HEP)
older = [11 12 14 20 21 22 23 32 37 38 41 43 47 48 49 52 55 57 58 63 64 65 67 69 7 70 71 75 8 83 86];
young = [13 15 16 25 26 28 31 33 34 36 4 42 44 45 46 50 51 53 54 56 59 6 62 66 68 72 74 76 78 80 82 84 85 9];
% exclude participants with no T peaks
older = setdiff(older, [21, 50, 55]);
young = setdiff(young, [21, 50, 55]);

% exclude participant with no pupil data
young = setdiff(young, 6);
older = setdiff(older, 38);

participants = {young older};
pupil_baseline = cell(2, 1);

for group = 1:2
    for p=1:length(participants{group})
        directory=strcat('G:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\AB', num2str(participants{group}(p)), '\ET\');

        participant=strcat('AB', num2str(participants{group}(p)));
        % only correct trials with errors and trials after errors excluded
        %load pupil baseline data created in pupil_analysis_5_pupil_baseline_RejectedEpochs.m 
        load([directory, 'PupilBaselinePerTaskRun']); % 5 task runs - passive task, 2 runs simple RT and 2 runs gng
        
        
        % pupil baseline participant number, passive, simple RT and gng - same number of
        % trials; + simple RT and gng - all trials
        pupil_simpleRT = [PupilBaselinePerTaskRun{2} PupilBaselinePerTaskRun{3}];
        pupil_gng = [PupilBaselinePerTaskRun{4} PupilBaselinePerTaskRun{5}];
        pupil_baseline{group} = [pupil_baseline{group}; participants{group}(p), group,...
            mean(PupilBaselinePerTaskRun{1}), ...
            mean(pupil_simpleRT(1:length(PupilBaselinePerTaskRun{1}))), mean(pupil_gng(1:length(PupilBaselinePerTaskRun{1}))), ...
            mean(pupil_simpleRT), mean(pupil_gng)];
   
    end
end

%% plot pupil baseline data
% plot_all_data_points(data_grp1_task1, data_grp1_task2, data_grp1_task3, data_grp2_task1, data_grp2_task2, data_grp2_task3, y_label_text)

% plot pupil baseline mean values
plot_all_data_points(pupil_baseline{1}(:, 3), pupil_baseline{1}(:, 4), pupil_baseline{1}(:, 5),...
    pupil_baseline{2}(:, 3), pupil_baseline{2}(:, 4), pupil_baseline{2}(:, 5),'Pupil baseline (mm)')


%% plot data from simple RT and gng - all trials
% plot_2conditions(grp1_cond1, grp1_cond2, grp2_cond1, grp2_cond2, y_label_text, x_label, title_text)
plot_2conditions(pupil_baseline{1}(:, 6), pupil_baseline{1}(:, 7), pupil_baseline{2}(:, 6), pupil_baseline{2}(:, 7),...
    'Pupil baseline (mm)', {'Simple RT' 'Go/no-go' 'Simple RT' 'Go/no-go'}, '')


%% create excel file for SPSS analysis

pupil_baseline4table = [pupil_baseline{1}; pupil_baseline{2}];


T = array2table(pupil_baseline4table, ...
    'VariableNames',{'part_id', 'group', 'mean_passive', 'mean_simpleRT', 'mean_gng','mean_simpleRT_alltrials', 'mean_gng_alltrials',});

filename = 'pupil_baseline_alltasks_samenumbtrials.xlsx';
writetable(T, filename,'Sheet',1,'Range','A1')


%% function to plot all data points 3 tasks 2 groups

function plot_all_data_points(data_grp1_task1, data_grp1_task2, data_grp1_task3, data_grp2_task1, data_grp2_task2, data_grp2_task3, y_label_text)

% young group
% define colour of data points
%     for i=1:length(data_grp1_task1)
%         cmap(i, :) = [floor(255/length(data_grp1_task1)*i)/255 floor(255/length(data_grp1_task1)*i)/255 floor(255/length(data_grp1_task1)*i)/255];
%     end

    % plot data for young group - simple RT and go/nogo task

        figure;
    for y=1:length(data_grp1_task1)
        plot([1 2 3]+rand*0.2-0.1, [data_grp1_task1(y) data_grp1_task2(y) data_grp1_task3(y)] ,'-o', 'color', [.8 .8 .8], ...
            'MarkerFaceColor',[.8 .8 .8], 'MarkerEdgeColor','k','MarkerSize',8, 'LineWidth', 1);
        hold on;
    end

    %plot the median line
    yMean=nanmean([data_grp1_task1 data_grp1_task2 data_grp1_task3], 1);
    plot([1 2 3],[yMean(1) yMean(2) yMean(3)] ,'Color','k','LineWidth',1.5);
    plot([1-0.3 1+0.3],[yMean(1) yMean(1)] ,'Color','k','LineWidth',5);
    plot([2-0.3 2+0.3],[yMean(2) yMean(2)] ,'Color','k','LineWidth',5);
    plot([3-0.3 3+0.3],[yMean(3) yMean(3)] ,'Color','k','LineWidth',5);
    

    % older group
%     for i=1:length(data_grp2_task1)
%         cmap(i, :) = [floor(255/length(data_grp2_task1)*i)/255 0 0];
%     end
    % plot data for older group - simple RT and go/nogo task

    
%     for y=1:length(data_grp2_task1)
%         plot([4 5 6]+rand*0.2-0.1, [data_grp2_task1(y) data_grp2_task2(y) data_grp2_task3(y)] ,'-o', 'color', cmap(y,:), ...
%             'MarkerFaceColor',cmap(y,:), 'MarkerEdgeColor','k','MarkerSize',8, 'LineWidth', 1);
%         hold on;  
%     end
    for y=1:length(data_grp2_task1)
        plot([4 5 6]+rand*0.2-0.1, [data_grp2_task1(y) data_grp2_task2(y) data_grp2_task3(y)] ,'-o', 'color', [1 .8 .8], ...
            'MarkerFaceColor',[1 .8 .8], 'MarkerEdgeColor','k','MarkerSize',8, 'LineWidth', 1);
        hold on;  
    end

     %plot the median line
    yMean=nanmean([data_grp2_task1 data_grp2_task2 data_grp2_task3], 1);
    plot([4 5 6],[yMean(1) yMean(2) yMean(3)] ,'Color','k','LineWidth',1.5);
    plot([4-0.3 4+0.3],[yMean(1) yMean(1)] ,'Color','k','LineWidth',5);
    plot([5-0.3 5+0.3],[yMean(2) yMean(2)] ,'Color','k','LineWidth',5);
    plot([6-0.3 6+0.3],[yMean(3) yMean(3)] ,'Color','k','LineWidth',5);
    

    % axes('XColor','none');
    hold off;
    axis([0 7 -inf inf]);
    ax = gca;
    c = ax.Color;
    ax.FontSize = 18;
    ax.FontName = 'Arial';
    ax.Color = 'none';
    ax.XTickLabel=[1 2 3 1 2 3];
    xticks([1 2 3 4 5 6])
    xticklabels({'1','2','3','1','2','3'})
    ylabel(y_label_text, 'FontSize', 24, 'FontWeight','normal')
    
%     x0=10;
%     y0=10;
%     width=400;
%     height=400;
%     set(gcf,'position',[x0,y0,width,height])
end



%% plot_2conditions(grp1_cond1, grp1_cond2, grp2_cond1, grp2_cond2, y_label_text, x_label, title_text)


function plot_2conditions(grp1_cond1, grp1_cond2, grp2_cond1, grp2_cond2, y_label_text, x_label, title_text)
    figure;
    index_nan1 = find(isnan(grp1_cond1));
    index_nan2 = find(isnan(grp1_cond2));

    grp1_cond1(unique([index_nan1, index_nan2])) = [];
    grp1_cond2(unique([index_nan1, index_nan2])) = [];

    % plot data for condition 1 group 1
    yMean1=nanmean(grp1_cond1);
    y_se = std(grp1_cond1,'omitnan')/sqrt(length(grp1_cond1));
    
        %plot the mean+-SEM box
        %   RECTANGLE('Position',pos) creates a rectangle in 2-D coordinates.
%   Specify pos as a four-element vector of the form [x y w h] in data
%   units. The x and y elements determine the location and the w and h
%   elements determine the size. The function plots into the current axes
%   without clearing existing content from the axes.
    box off
    rectangle('Position',[1-0.3,yMean1-y_se, 0.6, 2*y_se ],'FaceColor',[.7 .7 .7],'EdgeColor', [.7 .7 .7],'LineWidth',0.1);
    hold on
    %plot the mean line    
    plot([1-0.3 1+0.3],[yMean1 yMean1] ,'Color','k','LineWidth',5);
    
    % condition 2 group 1
    yMean2=nanmean(grp1_cond2);
    y_se = std(grp1_cond2,'omitnan')/sqrt(length(grp1_cond2));
    rectangle('Position',[2-0.3,yMean2-y_se, 0.6, 2*y_se ],'FaceColor',[.7 .7 .7],'EdgeColor', [.7 .7 .7],'LineWidth',0.1)
    %plot the mean line
    plot([2-0.3 2+0.3],[yMean2 yMean2] ,'Color','k','LineWidth',5);
    
    for y=1:length(grp1_cond1)
        plot([1+rand*0.2-0.1, 2+rand*0.2-0.1], [grp1_cond1(y), grp1_cond2(y)],'-o', 'color', [.8 .8 .8], ...
            'MarkerEdgeColor','k','MarkerSize',8, 'LineWidth', 1);
        hold on;
    end
    
    % plot data for group 2
    index_nan1 = find(isnan(grp2_cond1));
    index_nan2 = find(isnan(grp2_cond2));

    grp2_cond1(unique([index_nan1, index_nan2])) = [];
    grp2_cond2(unique([index_nan1, index_nan2])) = [];

    % plot data for condition 1
    yMean1=nanmean(grp2_cond1);
    y_se = std(grp2_cond1,'omitnan')/sqrt(length(grp2_cond1));
    
        %plot the mean+-SEM box
        %   RECTANGLE('Position',pos) creates a rectangle in 2-D coordinates.
%   Specify pos as a four-element vector of the form [x y w h] in data
%   units. The x and y elements determine the location and the w and h
%   elements determine the size. The function plots into the current axes
%   without clearing existing content from the axes.
    box off
    rectangle('Position',[4-0.3,yMean1-y_se, 0.6, 2*y_se ],'FaceColor',[.7 .7 .7],'EdgeColor', [.7 .7 .7],'LineWidth',0.1);
    hold on
    %plot the mean line    
    plot([4-0.3 4+0.3],[yMean1 yMean1] ,'Color','k','LineWidth',5);
    
    % condition 2 - group 2
    yMean2=nanmean(grp2_cond2);
    y_se = std(grp2_cond2,'omitnan')/sqrt(length(grp2_cond2));
    rectangle('Position',[5-0.3,yMean2-y_se, 0.6, 2*y_se ],'FaceColor',[.7 .7 .7],'EdgeColor', [.7 .7 .7],'LineWidth',0.1)
    %plot the mean line
    plot([5-0.3 5+0.3],[yMean2 yMean2] ,'Color','k','LineWidth',5);
    
    for y=1:length(grp2_cond1)
        plot([4+rand*0.2-0.1, 5+rand*0.2-0.1], [grp2_cond1(y), grp2_cond2(y)],'-o', 'color', [.8 .8 .8], ...
            'MarkerEdgeColor', 'k','MarkerSize',8, 'LineWidth', 1);
        hold on;
    end

    % axes('XColor','none');
    hold off;
    axis([0 6 -inf inf]);
    xticks([1 2 4 5])
    
    ax = gca;
    ax.LineWidth = 2.5; 
    ax.FontSize = 24;
    ax.FontName = 'Arial';
    ax.Color = 'none';
    ax.XTickLabel= x_label;
    xtickangle(45)
    ax.XAxis.FontSize = 18;

%     ax.XTickLabel= [];

    ylabel(y_label_text, 'FontSize', 26, 'FontWeight','normal')
    title(title_text, 'FontSize', 28, 'FontWeight','normal')

%     x0=10;
%     y0=10;
%     width=400;
%     height=400;
%     set(gcf,'position',[x0,y0,width,height])
end