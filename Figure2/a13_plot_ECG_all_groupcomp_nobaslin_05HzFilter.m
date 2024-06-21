% plot ECG - all time course
% all participants together and group differences
clear
load ECG4plot_avg;
load('D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\expected_chanlocs_both.mat');

% average across tasks
for grp =1:2
    for prt = 1:size(ECG4plot_avg{grp, 1}, 1)
        clear tmp
        tmp(:, 1) = squeeze(ECG4plot_avg{grp, 1}(prt, :));
        tmp(:, 2) = squeeze(ECG4plot_avg{grp, 2}(prt, :));
        ECG{grp}(prt, :) = mean(tmp, 2); 
    end
end

%% plot ECG both groups together, avg both tasks

ECG_all = cat(1, ECG{1}, ECG{2});
ECG_all_mean = squeeze(mean(ECG_all, 1));
ECG_all_se = squeeze(std(ECG_all, [], 1)/sqrt(size(ECG_all, 1)));

% level2_onesample_dir = 'D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\ECG_grp_task_effect\level2_onesample_ECG_all_participants';

x_axis=-.199:0.002:0.4; % in seconds
figure;
plot(x_axis, ECG_all_mean, 'color', 'k',  'LineWidth',1.5) %plot(light blue)
hold on
jbfill(x_axis, [ECG_all_mean + ECG_all_se], [ECG_all_mean - ECG_all_se], 'k', 'k', 1, 0.2)
hold on
plot(x_axis, zeros(length(x_axis), 1), '--k')
%     hold on
% 
%     %     % horizontal line marking where young group is different from zero
%     %     % find mask for channel chosen   
%     load([level2_onesample_dir, filesep, 'mask.mat'])
%     hold on
%     % plot horizontal line where the groups differ - mask starts on .1
%     % after R peak = .250 ms = 125 points
%     plot(x_axis(mask(chn, :) > 0) + .250, ones(length(find(mask(chn, :) > 0)), 1) * -1.15, 's', 'MarkerSize', 4.5, 'color', [.5 .5 .5], ...
%         'MarkerFaceColor', [.5 .5 .5])

hold off
box off
ax = gca;
%     axis([-inf inf -1.5 2.25])
ax.LineWidth = 2.5;
ax.FontSize = 28;
ax.FontName = 'Arial';
xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
title('T-locked ECG', 'FontSize', 32, 'FontWeight','normal')


%%
young_all_mean = squeeze(mean(ECG{1}, 1));
young_all_se = squeeze(std(ECG{1}, [], 1))/sqrt(size(ECG{1}, 1));

older_all_mean = squeeze(mean(ECG{2}, 1));
older_all_se = squeeze(std(ECG{2}, [], 1))/sqrt(size(ECG{2}, 1));

%% comparison across groups - including time windows with group differences
% level2_RM_ANOVA_dir = [pwd filesep 'level2_RM_ANOVA_ECG_grp_task_effect'];
% load([level2_RM_ANOVA_dir, filesep, 'mask_group_effect.mat'])
   
x_axis=-.199:0.002:0.4; % in seconds
figure;
plot(x_axis, young_all_mean, 'color', 'k',  'LineWidth',1.5) %plot(light blue)
hold on
jbfill(x_axis, [young_all_mean + young_all_se], [young_all_mean - young_all_se], 'k', 'k', 1, 0.2)
hold on
plot(x_axis, older_all_mean, '--r',  'LineWidth',1.5) 
hold on
jbfill(x_axis, [older_all_mean + older_all_se], [older_all_mean - older_all_se], 'r', 'r', 1, 0.2)
hold on
plot(x_axis, zeros(length(x_axis), 1), '--k')
hold on

%     % horizontal line marking where young group is different from zero
%     % find mask for channel chosen   
%     level2_onesample_dir = [pwd, '\level2_onesample_ECG'];
%     load([level2_onesample_dir, filesep, 'young', filesep, 'mask.mat'])
%     plot(x_axis(mask(chn, :) == 1), ones(length(find(mask(chn, :) == 1)), 1) * -1.25, 'k', 'LineWidth', 4.5)
%     hold on
%     load([level2_onesample_dir, filesep, 'older', filesep, 'mask.mat'])
%     plot(x_axis(mask(chn, :) > 0), ones(length(find(mask(chn, :) > 0)), 1) * -1.35, 's', 'MarkerSize', 4.5, 'color', 'r', ...
%         'MarkerFaceColor', 'r')
%     hold on
%     % plot horizontal line where the groups differ - mask starts on .1
%     % after R peak = .250 ms = 125 points
%     plot(x_axis(mask(chn, :) > 0) + .250, ones(length(find(mask(chn, :) > 0)), 1) * -1.15, 's', 'MarkerSize', 4.5, 'color', [.5 .5 .5], ...
%         'MarkerFaceColor', [.5 .5 .5])

hold off
box off
ax = gca;
%     axis([-inf inf -1.5 2.25])
ax.LineWidth = 2.5;
ax.FontSize = 28;
ax.FontName = 'Arial';
xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
title('T-locked ECG', 'FontSize', 32, 'FontWeight','normal')
        