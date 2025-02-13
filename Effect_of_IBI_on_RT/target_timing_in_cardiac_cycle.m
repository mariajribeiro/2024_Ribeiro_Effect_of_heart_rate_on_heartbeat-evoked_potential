% calculate HEP 2s interval before cue and 2 s interval after cue
% exclude period of auditory ERP
clear; close all;
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

older = [11 12 14 20 21 22 23 32 37 38 41 47 48 49 52 55 57 58 63 64 65 67 69 7 70 71 75 8 83 86];
young = [13 15 16 25 26 28 31 33 34 36 4 42 44 45 46 50 51 53 54 56 59 6 62 66 68 72 74 76 78 80 82 84 85 9];
% all files
files_folder = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\eeg_derivatives_05Hzfilt_cardiacICs';
full_files_folder = [files_folder filesep 'no_removed_segments'];
folder_content = dir(full_files_folder);

% subjects included
subjects = {};
for f = 1:length(folder_content)
    if contains(folder_content(f).name, 'eeg_full_Rpeak_Tpeak')
        subjects = [subjects, folder_content(f).name(26:29)];
    end
end
subjects = unique(subjects);

tasks = {'simpleRT', 'gonogo'};

RT_vs_cardiac_phase = cell(2,2);
RT_vs_IBI = cell(2,2);
yng = 0; old = 0;
for s = 1:length(subjects)
    if strcmp(subjects{s}(4), '_')
        sbj = 3;
    else
        sbj = 3:4;
    end
    
    subj_number = str2num(subjects{s}(sbj));
    if ismember(subj_number, young)
        group = 1; yng = yng + 1;
    else
        group = 2; old = old + 1;
    end
    
    RT_vs_cardiac_phase{group, 1} = [RT_vs_cardiac_phase{group, 1}; subj_number];
    RT_vs_IBI{group, 1} = [RT_vs_IBI{group, 1}; subj_number];
       
    for t = 1:length(tasks)
        % clear eeglab
        STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];  

        % load data
%         filename = {};
        for f = 1:length(folder_content)
            if contains(folder_content(f).name, subjects(s)) && contains(folder_content(f).name, tasks(t))
                filename = folder_content(f).name;
                EEG = pop_loadset('filename', filename, 'filepath', full_files_folder );
                EEG = pop_select( EEG, 'nochannel',{'EKG'});
                [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            end
        end

        % append data from both runs
        if length(ALLEEG) == 2
            EEG = pop_mergeset( ALLEEG, [1  2], 0);
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off'); 
        end
        
        %% R peaks before and after target
        EEG = pop_epoch( EEG, {  '2'  }, [-2        1], 'newname', 'epochs b4cue', 'epochinfo', 'yes');
        % calculate heart rate (or IBI) before cue
        IBI_b4target_tmp = []; Rpeaks = []; reaction_time = []; cardiac_phase = [];
        for ep = 1:length(EEG.epoch)
%             if ep == 49
%                 a = 1;
%             end
            % calculate IBI before target
            Rpeaks = find(strcmp(EEG.epoch(ep).eventtype, 'Rpeak'));
            R_latency = cell2mat(EEG.epoch(ep).eventlatency(Rpeaks));
            Rpeaks_b4target = find(R_latency < 0);
            if isempty(Rpeaks_b4target) || length(Rpeaks_b4target) == 1
                IBI_b4target_tmp(ep) = NaN;
            elseif length(Rpeaks_b4target) > 1
                IBI_b4target_tmp(ep) = R_latency(end) - R_latency(end-1);
            else
                stop = 1 
            end
            
            % calculate timing of target in relation to cardiac cycle
            Rpeaks_aftertarget = find(R_latency > 0);
            if ~isempty(Rpeaks_b4target) && ~isempty(Rpeaks_aftertarget)
                RR = R_latency(Rpeaks_aftertarget(1))- R_latency(Rpeaks_b4target(end));
                R2Target = abs(R_latency(Rpeaks_b4target(end)));
                cardiac_phase(ep) = R2Target/RR*2*pi;
            else
                cardiac_phase(ep) = NaN;
            end

            % calculate reaction time
            response = find(strcmp(EEG.epoch(ep).eventtype, '5'));
            if ~isempty(response) && length(response) == 1 && cell2mat(EEG.epoch(ep).eventlatency(response)) > 0
                reaction_time(ep) = cell2mat(EEG.epoch(ep).eventlatency(response));
            else
                reaction_time(ep) = NaN;
            end
        end
                
        
        % regression between reaction time and sin and cos of target timing in cardiac
        % cycle
        
        % [b,bint,r,rint,stats] = regress(y,X)
        y = reaction_time';
        x = [ones(length(reaction_time), 1), [sin(cardiac_phase)]', [cos(cardiac_phase)]'];
        [b,bint,r,rint,stats] = regress(y, x);
        
        % if s == 33
        %    a = 0; 
        % end
        % correlation between reaction time and IBI b4 target
        % find outliers in IBI
        zscore_IBI = (IBI_b4target_tmp - mean(IBI_b4target_tmp, 'omitnan'))/std(IBI_b4target_tmp, 'omitnan');
        IBI_outliers = find(abs(zscore_IBI) > 4);
        IBI_b4target_tmp(IBI_outliers) = NaN;
        
        reaction_time_nan = find(isnan(reaction_time));
        IBI_nan = find(isnan(IBI_b4target_tmp));
        reaction_time_without_nan = reaction_time;
        reaction_time_without_nan([reaction_time_nan, IBI_nan]) = [];
        
        IBI_b4target_tmp([reaction_time_nan, IBI_nan]) = [];
        
        [R,P] = corrcoef(reaction_time_without_nan, IBI_b4target_tmp);
        
        % figure; plot(reaction_time_without_nan, IBI_b4target_tmp, 'o')
        % regression including reaction time = dependent variable and
        % target timing and IBI and target timing x IBI as independent
        % variables
        
        % svae out put of regression and correlation
        if group == 1
            RT_vs_cardiac_phase{group, 2}(yng, t, :) = b;
            RT_vs_IBI{group, 2}(yng, t) = R(1, 2);
        else
            RT_vs_cardiac_phase{group, 2}(old, t, :) = b;
            RT_vs_IBI{group, 2}(old, t) = R(1, 2);
        end
    end       
end

save RT_vs_cardiac_phase RT_vs_cardiac_phase
save RT_vs_IBI RT_vs_IBI

%% plot figures
load RT_vs_IBI; load RT_vs_cardiac_phase;
% tasks = {'simpleRT', 'gonogo'};
% plot_2conditions(grp1_cond1, grp1_cond2, grp2_cond1, grp2_cond2, y_label_text, x_label, title_text)
plot_2conditions(RT_vs_IBI{1, 2}(:, 1), RT_vs_IBI{1, 2}(:, 2),...
    RT_vs_IBI{2, 2}(:, 1), RT_vs_IBI{2, 2}(:, 2), 'Correlation coefficient r', {'Young' 'Older'},...
    'RT vs IBI b4 target')
tasks = {'simpleRT', 'gonogo'};    
for t = 1:2
    plot_2conditions(RT_vs_cardiac_phase{1, 2}(:, t, 2), RT_vs_cardiac_phase{1, 2}(:, t, 3),...
        RT_vs_cardiac_phase{2, 2}(:, t, 2), RT_vs_cardiac_phase{2, 2}(:, t, 3),'Regression coefficients', {'Young' 'Older'},...
        tasks{t})
end


%% test coefficients dofference from zero

for t = 1:2
    for grp = 2
        [h,p,ci,stats] = ttest(RT_vs_cardiac_phase{grp, 2}(:, t, 2))
        [h,p,ci,stats] = ttest(RT_vs_cardiac_phase{grp, 2}(:, t, 3))
    end
end

%% mean std
for t = 1:2
    for grp = 2
        mean(RT_vs_cardiac_phase{grp, 2}(:, t, 2))
        std(RT_vs_cardiac_phase{grp, 2}(:, t, 2))

        mean(RT_vs_cardiac_phase{grp, 2}(:, t, 3))
        std(RT_vs_cardiac_phase{grp, 2}(:, t, 3))
    end
end

%%
for t = 1:2
    [h,p,ci,stats] = ttest([RT_vs_cardiac_phase{1, 2}(:, t, 2); RT_vs_cardiac_phase{1, 2}(:, t, 2)])

end
%% save data in excel file to perform stats in SPSS
% load RT_vs_IBI; load RT_vs_cardiac_phase; % CARRY ON FOR SPSS STATS|«!!!!
% % create table
% % participant, group, IBI_b4_passive,IBI_after_passive, IBI_b4_simpleRT,IBI_after_simpleRT, IBI_b4_gng,IBI_after_gng, 
% % IBI_avg_std{group, 2}(old, t, :) = [mean(IBI), std(IBI)];
% IBI_table = table([IBI_b4cue{1, 1};IBI_b4cue{2, 1}],...
%     [ones(length(IBI_b4cue{1, 1}), 1); ones(length(IBI_b4cue{2, 1}), 1)*2], ...
%     [IBI_b4cue{1, 2}(:, 1); IBI_b4cue{2, 2}(:, 1)], ...
%     [IBI_aftercue{1, 2}(:, 1); IBI_aftercue{2, 2}(:, 1)], ...
%     [IBI_b4cue{1, 2}(:, 2); IBI_b4cue{2, 2}(:, 2)], ...
%     [IBI_aftercue{1, 2}(:, 2); IBI_aftercue{2, 2}(:, 2)], ...
%     [IBI_b4cue{1, 2}(:, 3); IBI_b4cue{2, 2}(:, 3)], ...
%     [IBI_aftercue{1, 2}(:, 3); IBI_aftercue{2, 2}(:, 3)], ...
%     'VariableNames',{'participant', 'group', ' IBI_b4_passive' , 'IBI_after_passive',...
%     'IBI_b4_simpleRT', 'IBI_after_simpleRT', 'IBI_b4_gng', 'IBI_after_gng'});
% 
% filename = 'IBI_b4_after_cue.xlsx';
% writetable(IBI_table,filename)

%% functions

function plot_2conditions(grp1_cond1, grp1_cond2, grp2_cond1, grp2_cond2, y_label_text, x_label, title_text)
    figure;
    % line at zero
    plot([0 6], [0 0], '-k')
    hold on
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
            'MarkerEdgeColor','k','MarkerSize',6, 'LineWidth', .5);
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
            'MarkerEdgeColor','k','MarkerSize',6, 'LineWidth', .5);
        hold on;
    end
    
    
    % axes('XColor','none');
    hold off;
    axis([0 6 -inf inf]);
    xticks([1.5 4.5])
    
    ax = gca;
    ax.XTickLabel= x_label;
    ax.LineWidth = 2.5; 
    ax.FontSize = 24;
    ax.FontName = 'Arial';
    ax.Color = 'none';
%     ax.XTickLabel= [];
    ax.XAxis.FontSize = 24;
    ylabel(y_label_text, 'FontSize', 26, 'FontWeight','normal')
    title(title_text, 'FontSize', 28, 'FontWeight','normal')

%     x0=10;
%     y0=10;
%     width=400;
%     height=400;
%     set(gcf,'position',[x0,y0,width,height])
end